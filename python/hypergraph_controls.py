"""Count and order-specific participation preserving hypergraph null kernels."""
import numpy as np
from numba import njit
from kernels import evaluate

@njit(cache=True)
def degree(edges,n):
    g=np.zeros(n,np.int64)
    for e in edges:
        for j in range(n):
            if e&(1<<j):g[j]+=1
    return g


@njit(cache=True)
def swaps(e,k,n,seen,proposals):
    m=len(e);accepted=0
    if m<2:return 0
    for step in range(proposals):
        i=np.random.randint(m);j=np.random.randint(m)
        if i==j:continue
        x=e[i];y=e[j];px=np.random.randint(k);py=np.random.randint(k);vx=0;vy=0
        a=x
        for z in range(px):a=a&(a-1)
        vx=a&(-a);a=y
        for z in range(py):a=a&(a-1)
        vy=a&(-a)
        if (vx&y) or (vy&x):continue
        xx=x^vx^vy;yy=y^vy^vx
        if seen[xx] or seen[yy]:continue
        seen[x]=0;seen[y]=0;seen[xx]=1;seen[yy]=1;e[i]=xx;e[j]=yy;accepted+=1
    return accepted


@njit(cache=True)
def sample_count(c2,c3,m2,m3):
    # Partial Fisher-Yates; draws simple subsets uniformly.
    a=c2.copy();b=c3.copy()
    for i in range(m2):
        j=np.random.randint(i,len(a));a[i],a[j]=a[j],a[i]
    for i in range(m3):
        j=np.random.randint(i,len(b));b[i],b[j]=b[j],b[i]
    return np.concatenate((a[:m2],b[:m3]))


@njit(cache=True)
def simulate(edges,m2,n,c2,c3,model,B,seed,burn,thin):
    np.random.seed(seed);m=len(edges);m3=m-m2
    out=np.empty((B,23));accept=np.zeros((4,2),np.int64);proposals=np.zeros((4,2),np.int64)
    traces=np.zeros((B,2));orig=np.zeros(1<<n,np.int8)
    for e in edges:orig[e]=1
    want2=degree(edges[:m2],n);want3=degree(edges[m2:],n)
    fail=0
    for chain in range(4):
        a=edges[:m2].copy();b=edges[m2:].copy();seen2=np.zeros(1<<n,np.int8);seen3=np.zeros(1<<n,np.int8)
        for e in a:seen2[e]=1
        for e in b:seen3[e]=1
        comp2=m2>len(c2)//2;comp3=m3>len(c3)//2
        if model==1 and comp2:
            a=c2[seen2[c2]==0].copy();seen2[:]=0
            for e in a:seen2[e]=1
        if model==1 and comp3:
            b=c3[seen3[c3]==0].copy();seen3[:]=0
            for e in b:seen3[e]=1
        if model==1:
            accept[chain,0]+=swaps(a,2,n,seen2,burn*len(a));accept[chain,1]+=swaps(b,3,n,seen3,burn*len(b))
            proposals[chain,0]+=burn*len(a);proposals[chain,1]+=burn*len(b)
        for rep in range(B//4):
            ix=chain*(B//4)+rep
            if model==0:draw=sample_count(c2,c3,m2,m3)
            else:
                accept[chain,0]+=swaps(a,2,n,seen2,thin*len(a));accept[chain,1]+=swaps(b,3,n,seen3,thin*len(b))
                proposals[chain,0]+=thin*len(a);proposals[chain,1]+=thin*len(b)
                aa=c2[seen2[c2]==0] if comp2 else a;bb=c3[seen3[c3]==0] if comp3 else b
                draw=np.concatenate((aa,bb))
                if np.any(degree(aa,n)!=want2) or np.any(degree(bb,n)!=want3):fail+=1
            met,extra,_,_,_,_=evaluate(draw,m2);out[ix,:7]=met;out[ix,7:]=extra
            for e in draw[:m2]:traces[ix,0]+=orig[e]
            for e in draw[m2:]:traces[ix,1]+=orig[e]
            if m2:traces[ix,0]/=m2
            if m3:traces[ix,1]/=m3
    return out,accept,proposals,traces,fail
