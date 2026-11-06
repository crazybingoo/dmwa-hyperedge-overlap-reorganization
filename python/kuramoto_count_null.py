import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import numpy as np,pandas as pd
from concurrent.futures import ProcessPoolExecutor,as_completed
from simulate import H,pair,tri,masks,MET
from kernels import evaluate

def run(seed,B=40):
 out=H/'data'/f'kuramoto_null_{seed}.npz'
 if out.exists():return str(out)
 d=pd.read_csv(H/'data'/f'main_kuramoto_{seed}.csv').query("branch=='adaptive'").sort_values(['level','window'])
 rng=np.random.default_rng(700000+seed);res=np.full((B,6,16,7),np.nan)
 for _,r in d.iterrows():
  a=int(r.m2);b=int(r.m3);l=int(r.level);w=int(r.window)
  for z in range(B):
   ids2=rng.choice(len(pair),a,replace=False);ids3=rng.choice(len(tri),b,replace=False)
   edges=np.r_[masks[ids2],masks[len(pair)+ids3]]
   res[z,l,w]=evaluate(edges,a)[0]
 np.savez_compressed(out,values=res,seed=seed,B=B)
 return str(out)
if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=6) as p:
  for f in as_completed([p.submit(run,s) for s in range(1001,1041)]):print('DONE',f.result(),flush=True)
