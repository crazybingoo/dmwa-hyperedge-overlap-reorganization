import itertools
import numpy as np
import pandas as pd
from scipy.signal import hilbert
from scipy.special import xlogy

def geometry(n):
    pair=np.array(list(itertools.combinations(range(n),2)),dtype=int)
    tri=np.array(list(itertools.combinations(range(n),3)),dtype=int)
    all_edges=np.vstack([np.pad(pair,((0,0),(0,1)),constant_values=-1),tri])
    inc=np.zeros((len(all_edges),n),np.int8)
    for a,e in enumerate(all_edges): inc[a,e[e>=0]]=1
    pair_lookup={tuple(v):i for i,v in enumerate(pair)}
    tri_pairs=np.array([[pair_lookup[tuple(t[[0,1]])],pair_lookup[tuple(t[[0,2]])],pair_lookup[tuple(t[[1,2]])]] for t in tri])
    return pair,tri,inc,tri_pairs


def phase_scores(x,pair,tri):
    z=hilbert(x,axis=1); ph=np.angle(z); z=np.exp(1j*ph)
    plv=np.abs(z@z.conj().T)/x.shape[1]
    mu=np.angle(z.mean(1)); s=np.sin(ph-mu[:,None]); norms=np.linalg.norm(s,axis=1)
    if (norms<1e-12).any(): raise ValueError('Undefined circular correlation: zero norm')
    sn=s/norms[:,None]; circ=np.clip(sn@sn.T,-1,1)
    mat=np.ones((len(tri),3,3))
    for a,b in ((0,1),(0,2),(1,2)): mat[:,a,b]=mat[:,b,a]=circ[tri[:,a],tri[:,b]]
    eig=np.linalg.eigvalsh(mat)
    if eig.min()< -1e-9: raise ValueError('Non-PSD circular correlation')
    eig=np.maximum(eig,0); prob=eig/eig.sum(1)[:,None]
    coc=np.clip(1+xlogy(prob,prob).sum(1)/np.log(3),0,1)
    return plv[pair[:,0],pair[:,1]],coc


def threshold_curve(values):
    if not np.isfinite(values).all():raise ValueError('Nonfinite scores')
    grid=np.arange(401,dtype=float)*.0025
    vals=np.sort(values.ravel());density=(len(vals)-np.searchsorted(vals,grid,side='left'))/len(vals)
    eps=np.finfo(float).eps
    x=(grid-grid[0])/(grid[-1]-grid[0]+eps)
    y=(density-density[-1])/(density[0]-density[-1]+eps)
    distance=np.abs((y[-1]-y[0])*x-(x[-1]-x[0])*y+x[-1]*y[0]-y[-1]*x[0])/(np.hypot(y[-1]-y[0],x[-1]-x[0])+eps)
    ix=int(np.argmax(distance))
    flags=[]
    if density[0]-density[-1]<=eps:flags.append('degenerate_density_range')
    if ix in (0,400):flags.append('endpoint_selected')
    if distance.max()<1e-12:flags.append('flat_distance')
    return float(grid[ix]),pd.DataFrame({'threshold':grid,'density':density,'normalized_density':y,'distance':distance,'selected':np.arange(401)==ix}),flags
