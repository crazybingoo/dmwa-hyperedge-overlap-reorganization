"""Synthetic tests: no clinical inputs or identifiers are used."""
from pathlib import Path
import sys,itertools
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'python'))
from primary_pipeline import geometry,phase_scores,threshold_curve
from dmwa import build_metrics,LUT,motif_census
from bit_census import bit_census
from kernels import evaluate
from hypergraph_controls import sample_count,degree,swaps,simulate
from state_space import standardize_feature_array,exact_signflip_pvalue,bh_adjust
from simulate import kuramoto

rng=np.random.default_rng(1806)
for n in [4,6,8]:
 pair,tri,inc,tp=geometry(n)
 masks=np.array([sum(1<<int(j) for j in np.flatnonzero(x)) for x in inc],dtype=np.uint64)
 for repeat in range(12):
  a=np.sort(rng.choice(len(pair),min(len(pair),3+repeat),replace=False))
  b=np.sort(rng.choice(len(tri),min(len(tri),1+repeat),replace=False))
  ids=np.r_[a,len(pair)+b]; metrics,W,k,counts=build_metrics(inc[ids])
  fast=evaluate(masks[ids],len(a))[0]
  np.testing.assert_allclose(fast,[metrics[z] for z in ['S','D_G','D_H','D_PR','f_simple','f_recip','f_dense']],atol=2e-12,equal_nan=True)
  # The bitset census excludes disconnected triples; compare all 13 connected types.
  np.testing.assert_array_equal(bit_census(W,k,.85,LUT)[:,1:],counts[:,1:])
  assert metrics['ratio_max_error']<1e-12
  assert metrics['S_reconstruction_error']<1e-12
  assert metrics['q_reconstruction_error']<1e-12
  assert np.all(np.diag(W)==0)
  if np.isfinite(metrics['D_G']):np.testing.assert_allclose(metrics['D_G_finite_size'],metrics['D_G']*len(k)/(len(k)-1))
  # Independent member relabeling preserves weights and metrics.
  perm=rng.permutation(n); mm,WW,_,_=build_metrics(inc[ids][:,perm])
  np.testing.assert_allclose(WW,W)

pair,tri,inc,tp=geometry(6)
t=np.arange(1024)/1024
x=np.tile(np.cos(2*np.pi*10*t),(6,1))
plv,coc=phase_scores(x,pair,tri)
np.testing.assert_allclose(plv,1,atol=1e-12);np.testing.assert_allclose(coc,1,atol=1e-10)
th,curve,flags=threshold_curve(np.array([[.1,.2,.7,.8],[.15,.25,.75,.85]]))
assert len(curve)==401 and curve.selected.sum()==1 and 0<=th<=1
assert np.all(np.diff(curve.density)<=0)
empty,W,_,_=build_metrics(inc[:0]);assert np.isnan(empty['S']) and np.isnan(empty['D_G'])

masks=np.array([sum(1<<int(j) for j in np.flatnonzero(x)) for x in inc],dtype=np.uint64)
c2,c3=masks[:len(pair)].astype(np.int64),masks[len(pair):].astype(np.int64)
edges=sample_count(c2,c3,8,9);assert len(np.unique(edges))==17
a=edges[:8].copy();seen=np.zeros(1<<6,np.int8);seen[a]=1;before=degree(a,6)
swaps(a,2,6,seen,1000)
np.testing.assert_array_equal(degree(a,6),before);assert len(np.unique(a))==len(a)
for model in [0,1]:
 out,_,_,_,fail=simulate(edges,8,6,c2,c3,model,8,17,2,2)
 assert fail==0 and out.shape==(8,23)
v=rng.normal(size=(14,5,6));z=standardize_feature_array(v)
np.testing.assert_allclose(z.std(axis=(0,1),ddof=0),1)
d=np.ones((4,2));p,total,extreme,mode=exact_signflip_pvalue(d,np.sqrt(2))
assert p==2/16 and total==16 and mode=='exact_enumeration'
np.testing.assert_allclose(bh_adjust(np.array([.01,.03,.9])),[.03,.045,.9])
noise=rng.normal(size=(120,3));a=np.ones((3,3))-np.eye(3)
y=kuramoto(np.array([8.,8.1,7.9]),np.zeros(3),a,noise,100,np.array([0.,2.]),burn=0,segment=1)
assert y.shape==(3,120) and np.isfinite(y).all() and np.abs(y).max()<=1
print('PASS: 36 synthetic hypergraphs; exact census parity, DMWA identities, signal scores, null preservation, inference and simulation smoke test')
