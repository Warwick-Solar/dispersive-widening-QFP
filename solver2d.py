"""2D wave equation  d2v/dt2 = C_A^2(x) (d2v/dx2 + d2v/dz2), zero-beta slab, sausage pulse.
Units: w=1, C_Ai=1. Dirichlet v=0 on all boundaries (reflective), dv/dt(0)=0.
Second-order centred differences in space, leapfrog in time."""
import numpy as np

def setup(R=10.0, p=10, Lx=50.0, Lz=400.0, Nx=500, Nz=2200, z0=100.0, wx=1.0, wz=1.0):
    x = np.linspace(-Lx, Lx, Nx+1)
    z = np.linspace(0.0, Lz, Nz+1)
    dx, dz = x[1]-x[0], z[1]-z[0]
    with np.errstate(over='ignore'):
        rho = (R-1.0)/np.cosh(np.abs(x)**p)**2 + 1.0          # rho_e=1, rho_i=R
    CA2 = R/rho                                               # C_A^2 / C_Ai^2
    X, Z = np.meshgrid(x, z, indexing='ij')
    v0 = (X/wx)*np.exp(-X**2/wx**2 - (Z-z0)**2/wz**2)
    return dict(x=x, z=z, dx=dx, dz=dz, CA2=CA2, v0=v0, R=R, p=p, z0=z0)

def run(S, t_out, dt=0.02):
    x, z, dx, dz, CA2 = S['x'], S['z'], S['dx'], S['dz'], S['CA2']
    c2 = (CA2[:, None]*dt**2)
    cfl = np.sqrt(CA2.max())*dt*np.sqrt(1/dx**2+1/dz**2)
    assert cfl < 1, f'CFL {cfl}'
    def lap(v):
        L = np.zeros_like(v)
        L[1:-1,1:-1] = ((v[2:,1:-1]-2*v[1:-1,1:-1]+v[:-2,1:-1])/dx**2 +
                        (v[1:-1,2:]-2*v[1:-1,1:-1]+v[1:-1,:-2])/dz**2)
        return L
    v_prev = S['v0'].copy(); v_prev[0,:]=v_prev[-1,:]=0; v_prev[:,0]=v_prev[:,-1]=0
    v = v_prev + 0.5*c2*lap(v_prev)                   # dv/dt(0)=0
    t = dt; out = {}
    t_out = sorted(t_out); k = 0
    while k < len(t_out):
        if abs(t - t_out[k]) < dt/2:
            out[t_out[k]] = v.copy(); k += 1
            if k == len(t_out): break
        v_new = 2*v - v_prev + c2*lap(v)
        v_new[0,:]=v_new[-1,:]=0; v_new[:,0]=v_new[:,-1]=0
        v_prev, v = v, v_new; t += dt
    return out, cfl
