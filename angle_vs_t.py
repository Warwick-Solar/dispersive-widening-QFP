import numpy as np
from solver2d import setup, run
from scipy.signal import argrelextrema
def crest_analysis(v,x,z,z0,thr=0.10,x0s=0.5):
    s=z-z0; ix=np.argmin(abs(x-x0s)); sl=v[ix,:]; mz=(s>3); vmax=np.abs(sl[mz]).max()
    cand=np.r_[argrelextrema(sl,np.greater)[0],argrelextrema(sl,np.less)[0]]
    cand=np.sort([i for i in cand if s[i]>3 and abs(sl[i])>=thr*vmax])
    def ext(i,side):
        xs = x[x>0] if side>0 else -x[x<0][::-1]
        pr = np.abs(v[x>0,i]) if side>0 else np.abs(v[x<0,i][::-1])
        j=np.argmax(pr[xs<3]); m=pr[j]; k=j+np.argmax(pr[j:]<m/np.e); return xs[k]
    zj=s[cand]; ep=np.array([ext(i,1) for i in cand]); em=np.array([ext(i,-1) for i in cand])
    if len(zj)<2: return zj,ep,em,np.nan,None,None
    cP=np.polyfit(zj,ep,1); cM=np.polyfit(zj,em,1)
    th=np.degrees(np.arctan(abs(cP[0]))+np.arctan(abs(cM[0])))
    return zj,ep,em,th,cP,cM
if __name__=='__main__':
    ts=list(range(10,46,2))
    for Lx,Nx in [(50,500),(100,1000)]:
        S=setup(R=10.0,p=10,Lx=Lx,Nx=Nx); out,_=run(S,ts)
        print('Lx=',Lx)
        for t in ts:
            zj,ep,em,th,_,_=crest_analysis(out[t],S['x'],S['z'],S['z0'])
            print(f'  t={t:3d}  ncrest={len(zj)}  z-range={zj.min():5.1f}-{zj.max():5.1f}  ext-range={ep.min():4.1f}-{ep.max():4.1f}  theta={th:5.1f}')
