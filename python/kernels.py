import os
for x in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[x]='1'
import sys
from pathlib import Path
import numpy as np
from numba import njit
from bit_census import popcount

@njit(cache=True)
def evaluate(edges,m2,tau=.85):
    m=len(edges);d=np.zeros(m,np.int64);overlap=np.zeros((m,m),np.int8)
    for i in range(m):
        for j in range(i+1,m):
            c=int(popcount(np.uint64(edges[i]&edges[j])))
            if c:overlap[i,j]=c;d[i]+=1;d[j]+=1
    nb=(m+63)//64;U=np.zeros((m,nb),np.uint64);V=np.zeros((m,nb),np.uint64)
    state=np.zeros((m,m),np.int8);du=np.zeros(m,np.int64);dv=np.zeros(m,np.int64)
    net=np.zeros(m);total=np.zeros(m);ns=np.zeros(m);nx=np.zeros(m)
    low=0.;high=0.
    for i in range(m):
        ki=2 if i<m2 else 3
        for j in range(i+1,m):
            c=overlap[i,j]
            if c==0:continue
            kj=2 if j<m2 else 3
            if ki==kj:
                a=c/ki*d[j]/(d[i]+d[j]);b=c/kj*d[i]/(d[i]+d[j])
            else:a=c/ki;b=c/kj
            v=a-b;t=a+b;net[i]+=v;net[j]-=v;total[i]+=t;total[j]+=t
            if ki==kj:ns[i]+=v;ns[j]-=v
            else:nx[i]+=v;nx[j]-=v
            low+=min(a,b);high+=max(a,b)
            U[i,j//64]|=np.uint64(1)<<np.uint64(j%64);U[j,i//64]|=np.uint64(1)<<np.uint64(i%64)
            du[i]+=1;du[j]+=1
            if min(a,b)/max(a,b)+1e-12>=tau:state[i,j]=2
            else:
                state[i,j]=1;dv[i]+=1;dv[j]+=1
                V[i,j//64]|=np.uint64(1)<<np.uint64(j%64);V[j,i//64]|=np.uint64(1)<<np.uint64(i%64)
    out=np.full(7,np.nan)
    q=net/(total+np.finfo(np.float64).eps);qs=ns/(total+np.finfo(np.float64).eps);qx=nx/(total+np.finfo(np.float64).eps)
    if high>0:out[0]=low/high
    if m>1:
        s=np.maximum(q,0);t=np.maximum(-q,0)
        if s.sum()>0 and t.sum()>0:
            vals=np.zeros((2,3))
            for v in range(2):
                x=np.sort(s if v==0 else t);tot=x.sum();wg=0.;h=0.;sq=0.
                for k in range(m):
                    p=x[k]/tot;wg+=(k+1)*x[k];sq+=p*p
                    if p>0:h-=p*np.log(p)
                vals[v,0]=2*wg/(m*tot)-(m+1)/m;vals[v,1]=h/np.log(m);vals[v,2]=1/(sq*m)
            out[1]=vals[0,0]-vals[1,0];out[2]=vals[1,1]-vals[0,1];out[3]=vals[1,2]-vals[0,2]
    # Exact induced connected triples. U=all overlaps; V=nonreciprocal dyads.
    # Wedges_U-2*T_U counts connected triples; RVV removes reciprocal closings.
    tu3=0;tv3=0;rvv=0
    for i in range(m):
        for j in range(i+1,m):
            st=state[i,j]
            if not st:continue
            cu=0;cv=0
            for b in range(nb):
                cu+=int(popcount(U[i,b]&U[j,b]));cv+=int(popcount(V[i,b]&V[j,b]))
            tu3+=cu
            if st==1:tv3+=cv
            else:rvv+=cv
    tu=tu3//3;tv=tv3//3;wu=np.sum(du*(du-1))//2;wv=np.sum(dv*(dv-1))//2
    connected=wu-2*tu
    if connected:
        simple=wv-3*tv-rvv;nonmutual=wv-2*tv-rvv
        out[4]=simple/connected;out[5]=(connected-nonmutual)/connected;out[6]=(tu-tv)/connected
    # Per-order equal-hyperedge averages and four mutually exclusive states.
    extra=np.full(16,np.nan)
    for layer in range(2):
        lo=0 if layer==0 else m2;hi=m2 if layer==0 else m
        if hi<=lo:continue
        extra[layer*8:layer*8+8]=0.
        for i in range(lo,hi):
            b=layer*8;extra[b]+=q[i];extra[b+1]+=qs[i];extra[b+2]+=qx[i];extra[b+3]+=d[i]
            if d[i]==0:extra[b+7]+=1
            elif q[i]>1e-12:extra[b+4]+=1
            elif q[i]<-1e-12:extra[b+5]+=1
            else:extra[b+6]+=1
        extra[layer*8:layer*8+8]/=(hi-lo)
    return out,extra,q,qs,qx,d
