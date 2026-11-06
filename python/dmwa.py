"""DMWA weights, role gaps, decomposition and exact 13-triad census."""
import itertools
import numpy as np
from scipy.special import xlogy
from numba import njit

def lut13():
    templates=[[[0,1,1],[0,0,0],[0,0,0]],[[0,0,0],[1,0,0],[1,0,0]],
    [[0,1,0],[0,0,1],[0,0,0]],[[0,1,0],[1,0,0],[1,0,0]],
    [[0,1,1],[1,0,0],[0,0,0]],[[0,1,1],[0,0,1],[0,0,0]],
    [[0,1,0],[0,0,1],[1,0,0]],[[0,1,1],[1,0,0],[1,0,0]],
    [[0,1,1],[1,0,1],[0,0,0]],[[0,1,0],[1,0,0],[1,1,0]],
    [[0,1,1],[1,0,0],[0,1,0]],[[0,1,1],[1,0,1],[1,0,0]],
    [[0,1,1],[1,0,1],[1,1,0]]]
    out=np.zeros(64,np.int64)
    for code,t in enumerate(templates,1):
        t=np.array(t)
        for perm in itertools.permutations(range(3)):
            b=t[np.ix_(perm,perm)]
            c=b[0,1]+2*b[0,2]+4*b[1,0]+8*b[1,2]+16*b[2,0]+32*b[2,1]
            assert out[c] in (0,code)
            out[c]=code
    return out


@njit(cache=True)
def motif_census(W, orders, tau, lut):
    m=len(orders)
    adj=np.zeros((m,m),np.int8)
    for i in range(m):
        for j in range(i+1,m):
            hi=max(W[i,j],W[j,i]); lo=min(W[i,j],W[j,i])
            if hi<=0: continue
            if lo/hi+1e-12>=tau:
                adj[i,j]=1; adj[j,i]=1
            elif W[i,j]>W[j,i]: adj[i,j]=1
            else: adj[j,i]=1
    # first column disconnected, then manuscript triads 1..13
    counts=np.zeros((4,14),np.int64)
    for i in range(m-2):
        for j in range(i+1,m-1):
            ij=adj[i,j]+4*adj[j,i]
            for k in range(j+1,m):
                c=ij+2*adj[i,k]+8*adj[j,k]+16*adj[k,i]+32*adj[k,j]
                typ=orders[i]+orders[j]+orders[k]-6
                counts[typ,lut[c]]+=1
    return counts


def load_gaps(q):
    m=len(q)
    def parts(x):
        total=x.sum()
        if m<2 or total<=0: return (np.nan,np.nan,np.nan)
        p=x/total; g=2*np.dot(np.arange(1,m+1),np.sort(x))/(m*total)-(m+1)/m
        return g, -xlogy(p,p).sum()/np.log(m), 1/(np.dot(p,p)*m)
    gs,hs,ps=parts(np.maximum(q,0)); gt,ht,pt=parts(np.maximum(-q,0))
    return gs-gt,ht-hs,pt-ps


def build_metrics(incidence, include_motifs=True):
    k=incidence.sum(1).astype(np.int64); m=len(k)
    O=incidence.astype(float)@incidence.T
    np.fill_diagonal(O,0)
    d=(O>0).sum(1)
    same=k[:,None]==k[None,:]
    den_d=d[:,None]+d[None,:]
    bias=np.divide(d[None,:],den_d,out=np.zeros((m,m)),where=den_d>0)
    W=(O/k[:,None])*np.where(same,bias,1) if m else O
    i,j=np.triu_indices(m,1); a=np.maximum(W[i,j],W[j,i]); b=np.minimum(W[i,j],W[j,i])
    A=a.sum(); S=b.sum()/A if A>0 else np.nan
    result={'m2':int((k==2).sum()),'m3':int((k==3).sum()),'m':m,
            'S':S,'overlap_pairs':int((a>0).sum()),'isolated_count':int((d==0).sum())}
    check=[]
    for lab,k1,k2 in [('22',2,2),('23',2,3),('33',3,3)]:
        use=((k[i]==k1)&(k[j]==k2))|((k[i]==k2)&(k[j]==k1))
        ag=a[use].sum(); bg=b[use].sum()
        result.update({f'A{lab}':ag,f'B{lab}':bg,f'P{lab}':ag/A if A else np.nan,
            f'h{lab}':bg/ag if ag else np.nan,f'T{lab}':(bg-2/3*ag)/A if A else np.nan,
            f'pairs{lab}':int((a[use]>0).sum())})
        act=use&(a>0)
        expected=np.full(act.sum(),2/3) if k1!=k2 else np.minimum(d[i[act]],d[j[act]])/np.maximum(d[i[act]],d[j[act]])
        check.append(float(np.max(np.abs(b[act]/a[act]-expected))) if act.any() else 0.)
    result['S_reconstruction_error']=abs(S-(2/3+result['T22']+result['T33'])) if A else np.nan
    result['ratio_max_error']=max(check)
    so=W.sum(1); si=W.sum(0); denom=so+si+np.finfo(float).eps
    q=(so-si)/denom
    ws=W*same; wx=W*(~same)
    qs=(ws.sum(1)-ws.sum(0))/denom; qx=(wx.sum(1)-wx.sum(0))/denom
    result['q_reconstruction_error']=float(np.max(abs(q-qs-qx))) if m else 0.
    result['D_G'],result['D_H'],result['D_PR']=load_gaps(q)
    result['D_G_finite_size']=result['D_G']*m/(m-1) if m>1 else np.nan
    for order in (2,3):
        sel=k==order
        for name,v in [('q',q),('q_same',qs),('q_cross',qx),('overlap_degree',d)]:
            result[f'{name}_{order}']=float(v[sel].mean()) if sel.any() else np.nan
        for name,v in [('positive',q>1e-12),('negative',q< -1e-12),('balanced_active',(abs(q)<=1e-12)&(d>0)),('isolated',d==0)]:
            result[f'{name}_{order}']=float(v[sel].mean()) if sel.any() else np.nan
    result['nested_pair_count']=int(np.count_nonzero((O==2)&(k[:,None]==2)&(k[None,:]==3)))
    if include_motifs:
        counts=motif_census(W,k,.85,LUT)
        result.update(motif_summary(counts))
    else: counts=None
    return result,W,k,counts


def motif_summary(counts):
    total=counts[:,1:].sum()
    groups={'simple':[1,2,3],'ffl':[6],'cycle':[7],'recip':[4,5,8,9,10,11,12,13], 'dense':[9,10,11,12,13],'complex':list(range(4,14))}
    r={'motif_connected':int(total),'motif_candidates':int(counts.sum())}
    for g,codes in groups.items(): r['f_'+g]=float(counts[:,codes].sum()/total) if total else np.nan
    for h,lab in enumerate(['222','223','233','333']):
        den=counts[h,1:].sum(); r['connected_'+lab]=int(den)
        r['Q_'+lab]=float(den/total) if total else np.nan
        for g in ('simple','recip','dense'):
            num=counts[h,groups[g]].sum()
            r[f'f_{g}_{lab}']=float(num/den) if den else np.nan
            r[f'contribution_{g}_{lab}']=float(num/total) if total else np.nan
    return r

LUT=lut13()
