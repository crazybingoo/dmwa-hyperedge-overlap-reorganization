import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import sys,itertools,json,time,argparse
from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx
from numba import njit
from scipy.signal import butter,sosfiltfilt
from concurrent.futures import ProcessPoolExecutor,as_completed
from kernels import evaluate
H=Path(__file__).resolve().parents[1]/'outputs'/'simulation'
(H/'data').mkdir(parents=True,exist_ok=True)
from primary_pipeline import phase_scores,geometry,threshold_curve
MET=['S','D_G','D_H','D_PR','f_simple','f_recip','f_dense']
FS=1024;N=18;LEVELS=np.array([0.,2.,4.,6.,8.,12.]);MIX=np.array([0.,.15,.30,.45,.60,.75])
pair,tri,inc,tp=geometry(N)
masks=np.array([sum(1<<int(j) for j in np.flatnonzero(x)) for x in inc],dtype=np.uint64)

@njit(cache=True)
def kuramoto(freq,theta,A,noise,fs,levels,burn=20,segment=24):
 dt=1/fs;n=len(freq);out=np.zeros((n,len(noise)-burn*fs));cur=theta.copy()
 for t in range(len(noise)):
  lev=0 if t<burn*fs else min((t-burn*fs)//(segment*fs),len(levels)-1)
  nxt=np.empty(n)
  for i in range(n):
   interaction=0.
   for j in range(n):
    if A[i,j]:interaction+=np.sin(cur[j]-cur[i])
   nxt[i]=cur[i]+dt*(2*np.pi*freq[i]+levels[lev]*interaction/6)+np.sqrt(.6*dt)*noise[t,i]
  cur=nxt%(2*np.pi)
  if t>=burn*fs:
   for i in range(n):out[i,t-burn*fs]=np.cos(cur[i])
 return out

def diagnostics(edges,m2):
 vals,extra,q,qs,qx,d=evaluate(edges,m2)
 m=len(edges);m3=m-m2;deg2=np.zeros(N);deg3=np.zeros(N);nest=0
 for i,e in enumerate(edges):
  for j in range(N):
   if int(e)&(1<<j):
    if i<m2:deg2[j]+=1
    else:deg3[j]+=1
 for a in edges[:m2]:nest+=sum((a&b)==a for b in edges[m2:])
 cv=lambda x:float(x.std()/x.mean()) if x.mean()>0 else np.nan
 return dict(zip(MET,vals))|{'m2':m2,'m3':m3,'density2':m2/len(pair),'density3':m3/len(tri),'cv2':cv(deg2),'cv3':cv(deg3),'overlap_cv':cv(d),'nesting':nest/(3*m3) if m3 else np.nan}

def signal_run(model,seed,pilot=False):
 prefix='pilot' if pilot else 'main';out=H/'data'/f'{prefix}_{model}_{seed}.csv'
 if out.exists():return str(out)
 rng=np.random.default_rng(seed);T=6*24*FS
 if model=='kuramoto':
  G=nx.random_regular_graph(6,N,seed=seed);A=nx.to_numpy_array(G)
  freq=rng.normal(8,.5,N);theta=rng.uniform(0,2*np.pi,N)
  noise=rng.normal(size=(T+20*FS,N));x=kuramoto(freq,theta,A,noise,FS,LEVELS)
  x+=.2*rng.normal(size=x.shape);levels=LEVELS
 elif model=='common_input':
  # Stationary Gaussian carriers, unit variance per full realization.
  sos=butter(4,[4,30],fs=FS,btype='bandpass',output='sos')
  carriers=sosfiltfilt(sos,rng.normal(size=(N+1,T+4*FS)),axis=1)[:,2*FS:-2*FS]
  carriers/=carriers.std(axis=1,keepdims=True);x=np.zeros((N,T));levels=MIX
  for l,a in enumerate(MIX):
   sl=slice(l*24*FS,(l+1)*24*FS)
   x[:,sl]=np.sqrt(1-a)*carriers[:N,sl]+np.sqrt(a)*carriers[N,sl]
  x+=.2*rng.normal(size=x.shape)
 else:raise ValueError(model)
 p=[];c=[];ix=[]
 for l in range(6):
  for w in range(16):
   start=(24*l+6+w)*FS
   pp,cc=phase_scores(x[:,start:start+3*FS],pair,tri);p.append(pp);c.append(cc);ix.append((l,w))
 p=np.array(p);c=np.array(c)
 th2,_,f2=threshold_curve(p);th3,_,f3=threshold_curve(c)
 rows=[]
 for j,(l,w) in enumerate(ix):
  for branch in ['adaptive','fixed_counts']:
   if branch=='adaptive':ids2=np.flatnonzero(p[j]>=th2);ids3=np.flatnonzero(c[j]>=th3)
   else:ids2=np.argsort(-p[j],kind='stable')[:round(.2*len(pair))];ids3=np.argsort(-c[j],kind='stable')[:round(.05*len(tri))]
   edges=np.r_[masks[ids2],masks[len(pair)+ids3]]
   rows.append({'model':model,'seed':seed,'branch':branch,'level':l,'strength':levels[l],'window':w,'theta2':th2,'theta3':th3,'mean_plv':p[j].mean(),'mean_coc':c[j].mean(),**diagnostics(edges,len(ids2))})
 pd.DataFrame(rows).to_csv(out,index=False)
 np.savez_compressed(H/'data'/f'{prefix}_{model}_{seed}_scores.npz',plv=p,coc=c,theta2=th2,theta3=th3)
 if seed in [101,1001]:np.savez_compressed(H/'data'/f'{prefix}_{model}_{seed}_signal.npz',x=x,fs=FS)
 return str(out)

def tri_sequence(degrees,rng):
 for attempt in range(2000):
  left=degrees.copy();edges=[];seen=set();failed=False
  while left.sum():
   if (left>0).sum()<3:failed=True;break
   ok=False
   for trial in range(60):
    nodes=np.sort(rng.choice(N,size=3,replace=False,p=left/left.sum()));key=tuple(nodes)
    if key not in seen:ok=True;break
   if not ok:failed=True;break
   edges.append(nodes);seen.add(key);left[nodes]-=1
  if not failed:return np.array(edges,dtype=np.int64)
 raise RuntimeError('Could not generate simple degree sequence')

@njit(cache=True)
def mix_triples(t,A,beta,proposals,seed):
 np.random.seed(seed);m=len(t);out=t.copy();bits=np.zeros(m,np.uint64)
 for i in range(m):
  for j in range(3):bits[i]|=np.uint64(1)<<np.uint64(out[i,j])
 accepted=0
 for step in range(proposals):
  i=np.random.randint(m);j=np.random.randint(m)
  if i==j:continue
  a=np.random.randint(3);b=np.random.randint(3);u=out[i,a];v=out[j,b]
  if u==v or v in out[i] or u in out[j]:continue
  bi=bits[i]^(np.uint64(1)<<np.uint64(u))^(np.uint64(1)<<np.uint64(v))
  bj=bits[j]^(np.uint64(1)<<np.uint64(v))^(np.uint64(1)<<np.uint64(u))
  duplicate=False
  for z in range(m):
   if z!=i and z!=j and (bits[z]==bi or bits[z]==bj):duplicate=True;break
  if duplicate:continue
  delta=0
  for z in range(3):
   if z!=a:delta+=A[v,out[i,z]]-A[u,out[i,z]]
   if z!=b:delta+=A[u,out[j,z]]-A[v,out[j,z]]
  if delta>=0 or np.random.random()<np.exp(beta*delta):
   out[i,a]=v;out[j,b]=u;bits[i]=bi;bits[j]=bj;accepted+=1
 return out,accepted

def graph_run(seed,profile,pilot=False):
 prefix='pilot' if pilot else 'main';out=H/'data'/f'{prefix}_graph_{profile}_{seed}.csv'
 if out.exists():return str(out)
 rng=np.random.default_rng(seed);base=np.tile([1,1,2],6) if profile=='heterogeneous' else np.ones(N,int)
 rows=[]
 for level in range(1,6):
  dg=base*level;G=nx.havel_hakimi_graph(dg.tolist())
  nx.double_edge_swap(G,nswap=30*G.number_of_edges(),max_tries=2000*G.number_of_edges(),seed=int(rng.integers(2**30)))
  A=nx.to_numpy_array(G,dtype=np.int64);ps=np.array(sorted(tuple(sorted(e)) for e in G.edges()))
  pm=np.array([sum(1<<int(v) for v in e) for e in ps],dtype=np.uint64)
  ts=tri_sequence(3*dg,rng)
  for beta in [0.,1.5,3.]:
   ts,acc=mix_triples(ts,A,beta,100*len(ts),int(rng.integers(2**30)))
   tm=np.array([sum(1<<int(v) for v in e) for e in ts],dtype=np.uint64)
   assert len(np.unique(tm))==len(tm)
   assert np.array_equal(np.bincount(ts.ravel(),minlength=N),3*dg)
   assert np.array_equal(np.array([G.degree(i) for i in range(N)]),dg)
   rows.append({'model':'graph','seed':seed,'profile':profile,'level':level,'beta':beta,'accepted_swaps':acc,**diagnostics(np.r_[pm,tm],len(pm))})
 pd.DataFrame(rows).to_csv(out,index=False);return str(out)

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--pilot',action='store_true');args=parser.parse_args()
 seeds=range(101,105) if args.pilot else range(1001,1041)
 jobs=[]
 with ProcessPoolExecutor(max_workers=6) as pool:
  for model in ['kuramoto','common_input']:
   for seed in seeds:jobs.append(pool.submit(signal_run,model,seed,args.pilot))
  for profile in ['heterogeneous','homogeneous']:
   for seed in (range(101,105) if args.pilot else range(1001,1061)):jobs.append(pool.submit(graph_run,seed,profile,args.pilot))
  for f in as_completed(jobs):print('DONE',f.result(),flush=True)
if __name__=='__main__':
 sys.stdout.reconfigure(encoding='utf8');main()
