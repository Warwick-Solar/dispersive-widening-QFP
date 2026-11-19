"""Data for the animation of the Fig. 2 run (rho_i/rho_e = 10, p = 10, Lx = 150): the same solver and grid as
solver2d.run, but the field is stored every NSAVE time steps (0.2 w/C_Ai) only in the displayed window
(|x| <= XW, SMIN <= z - z0 <= SMAX), and the crest analysis of the figures (angle_vs_t.crest_analysis, on the full
field) is done on the fly for t >= T_FIT.  Output: anim_data.npz (window frames, float32, ~210 MB) + anim_crests.json.
Run:  python3 anim_data.py   (~5 min; needs solver2d.py and angle_vs_t.py in the same folder), then anim_render.py.
The run is continued to t = 93; the animation stops at t = 90 (see anim_render.py)."""
import numpy as np, json, time
from solver2d import setup
from angle_vs_t import crest_analysis
TEND, NSAVE, T_FIT = 93.0, 10, 20.0
XW, SMIN, SMAX = 10.0, -5.0, 200.0
S = setup(R=10.0, p=10, Lx=150.0, Nx=1500); x, z, z0 = S['x'], S['z'], S['z0']; s = z - z0
dx, dz, dt = S['dx'], S['dz'], 0.02
c2 = S['CA2'][:, None]*dt**2
def lap(v):
    L = np.zeros_like(v)
    L[1:-1,1:-1] = ((v[2:,1:-1]-2*v[1:-1,1:-1]+v[:-2,1:-1])/dx**2 + (v[1:-1,2:]-2*v[1:-1,1:-1]+v[1:-1,:-2])/dz**2)
    return L
xm = np.abs(x) <= XW + 1e-9; zm = (s >= SMIN) & (s <= SMAX)
frames, times, crests = [], [], []
def save(v, t):
    frames.append(v[np.ix_(xm, zm)].astype(np.float32)); times.append(t)
    if t >= T_FIT - 1e-9:
        zj, ep, em, th, cP, cM = crest_analysis(v, x, z, z0, thr=0.10)
        crests.append(dict(t=t, zj=zj.tolist(), ep=ep.tolist(), em=em.tolist(), th=float(th),
                           cP=None if cP is None else cP.tolist(), cM=None if cM is None else cM.tolist()))
v_prev = S['v0'].copy(); v_prev[0,:]=v_prev[-1,:]=0; v_prev[:,0]=v_prev[:,-1]=0
save(v_prev, 0.0)
v = v_prev + 0.5*c2*lap(v_prev); n = 1; t = dt          # identical to solver2d.run
t0 = time.time(); nend = int(round(TEND/dt))
while n <= nend:
    if n % NSAVE == 0: save(v, round(n*dt, 6))
    if n % 1000 == 0: print(f'step {n}/{nend}  t={n*dt:.1f}  {time.time()-t0:.0f}s', flush=True)
    v_new = 2*v - v_prev + c2*lap(v)
    v_new[0,:]=v_new[-1,:]=0; v_new[:,0]=v_new[:,-1]=0
    v_prev, v = v, v_new; n += 1
np.savez('anim_data.npz', x=x[xm], s=s[zm], t=np.array(times), V=np.array(frames))
json.dump(crests, open('anim_crests.json', 'w'))
print('done', len(times), 'frames', time.time()-t0, 's')
