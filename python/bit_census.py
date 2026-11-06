"""Exact bitset implementation of the induced three-node census."""
import numpy as np
from numba import njit, types
from numba.extending import intrinsic
from llvmlite import ir

@intrinsic
def popcount(typingctx,x):
    sig=types.uint64(types.uint64)
    def codegen(context,builder,signature,args):
        fn=builder.module.declare_intrinsic('llvm.ctpop',[ir.IntType(64)])
        return builder.call(fn,args)
    return sig,codegen


@njit(cache=True)
def bit_census(W,k,tau,lut):
    """Exact induced 3-node census by pairwise bitset intersections; no sampling."""
    m=len(k); nb=(m+63)//64
    masks=np.zeros((m,4,nb),np.uint64)
    state=np.zeros((m,m),np.int8)
    for i in range(m):
        for j in range(i+1,m):
            hi=max(W[i,j],W[j,i]);lo=min(W[i,j],W[j,i])
            if hi<=0:continue
            if lo/hi+1e-12>=tau:state[i,j]=3;state[j,i]=3
            elif W[i,j]>W[j,i]:state[i,j]=1;state[j,i]=2
            else:state[i,j]=2;state[j,i]=1
    order2=np.zeros(nb,np.uint64)
    for j in range(m):
        bit=np.uint64(1)<<np.uint64(j%64)
        if k[j]==2:order2[j//64]|=bit
        for i in range(m):masks[i,state[i,j],j//64]|=bit
    counts=np.zeros((4,14),np.int64)
    lookup=np.zeros((4,4,4),np.int64)
    for s in range(4):
        for a in range(4):
            for b in range(4):
                code=(s&1)+4*(s>>1)+2*(a&1)+16*(a>>1)+8*(b&1)+32*(b>>1)
                lookup[s,a,b]=lut[code]
    for i in range(m-2):
        for j in range(i+1,m-1):
            s=state[i,j];g=k[i]+k[j]-4;start=(j+1)//64
            firstmask=np.uint64(0xffffffffffffffff)<<np.uint64((j+1)%64)
            for a in range(4):
                for b in range(4):
                    c=lookup[s,a,b]
                    if c==0:continue
                    for block in range(start,nb):
                        bits=masks[i,a,block]&masks[j,b,block]
                        if block==start:bits&=firstmask
                        total=np.int64(popcount(bits))
                        if total:
                            n2=np.int64(popcount(bits&order2[block]))
                            counts[g,c]+=n2;counts[g+1,c]+=total-n2
    n2=np.sum(k==2);n3=m-n2
    candidates=np.array([n2*(n2-1)*(n2-2)//6,n2*(n2-1)//2*n3,n2*n3*(n3-1)//2,n3*(n3-1)*(n3-2)//6])
    for g in range(4):counts[g,0]=candidates[g]-counts[g,1:].sum()
    return counts
