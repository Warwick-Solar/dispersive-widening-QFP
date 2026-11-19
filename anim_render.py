"""Animation of the Fig. 2 wave train (rho_i/rho_e = 10, p = 10): v_x/max|v_x| in the window, colour saturated at
+-0.3, dashed slab boundaries; from t = T_FIT the crest analysis of Fig. 2 (circles = 1/e extents of the crests above
10%, lines = linear fits, full opening angle).  Same scale (inches per w) and fonts as Fig. 2, so the frame at t = 40
looks like Fig. 2, but the z-axis covers the path of the train ((z - z0)/w = -3..110).  Reads anim_data.npz and
anim_crests.json (made by anim_data.py).
Run:  python3 anim_render.py            -> wavetrain_animation.mp4 (t = 0..90, 451 frames, 25 fps, 18 s)
                                           + every frame as frames/frame_t<t>.png
      FRAME=40 python3 anim_render.py   -> single frame anim_frame_t40.png
Options (environment variables): ZLIM (default -3,110), TEND (90), EVERY (1), FPS (25), DPI (105), OUT, FRAMES_DIR.
t_end = 90: reflections of the leaky waves from x = +-150 enter the displayed band |x| <= 10 at t ~ 91."""
import os, json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
plt.rcParams.update({'font.size':12,'axes.labelsize':13,'xtick.labelsize':12,'ytick.labelsize':12})
SMIN, SMAX = [float(a) for a in os.environ.get('ZLIM', '-3,110').split(',')]
TEND = float(os.environ.get('TEND', '90')); EVERY = int(os.environ.get('EVERY', '1')); FPS = int(os.environ.get('FPS', '25'))
XLIM, CLIP, DPI = 10.0, 0.3, int(os.environ.get('DPI', '105'))
OUT = os.environ.get('OUT', 'wavetrain_animation.mp4'); FRAMES_DIR = os.environ.get('FRAMES_DIR', 'frames')
IN_PER_W = 4.1168/38.0                      # Fig. 2: axes 4.1168 in wide for (z - z0)/w = 22..60
AXH, LEFT, BOT, TOP, CBW, PAD = 2.4640, 0.80, 0.62, 0.12, 0.1232, 0.10
d = np.load('anim_data.npz'); x, s, T, V = d['x'], d['s'], d['t'], d['V']
cr = {round(c['t'], 6): c for c in json.load(open('anim_crests.json'))}
m = (s >= SMIN) & (s <= SMAX); s = s[m]; V = V[:, :, m]
idx = [i for i in range(len(T)) if T[i] <= TEND + 1e-9][::EVERY]
axw = IN_PER_W*(SMAX - SMIN); W = LEFT + axw + PAD + CBW + 0.75; H = BOT + AXH + TOP
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([LEFT/W, BOT/H, axw/W, AXH/H]); cax = fig.add_axes([(LEFT+axw+PAD)/W, (BOT+0.112)/H, CBW/W, 2.24/H])
pc = ax.pcolormesh(s, x, V[0]/np.abs(V[0]).max(), cmap='RdBu_r', vmin=-CLIP, vmax=CLIP, shading='auto', rasterized=True)
ax.axhline(1, ls='--', c='k', lw=0.6); ax.axhline(-1, ls='--', c='k', lw=0.6)
cP, = ax.plot([], [], 'o', mfc='w', mec='k', ms=5); cM, = ax.plot([], [], 'o', mfc='w', mec='k', ms=5)
fP, = ax.plot([], [], 'k', lw=1); fM, = ax.plot([], [], 'k', lw=1)
tl = ax.text(0.006, 0.95, '', transform=ax.transAxes, ha='left', va='top', fontsize=12)
al = ax.text(0.006, 0.05, '', transform=ax.transAxes, ha='left', va='bottom', fontsize=12)
ax.set_xlim(SMIN, SMAX); ax.set_ylim(-XLIM, XLIM); ax.set_xlabel('$(z-z_0)/w$'); ax.set_ylabel('$x/w$')
cb = fig.colorbar(pc, cax=cax, extend='both'); cb.set_label('$v_x/\\max|v_x|$', fontsize=11); cb.ax.tick_params(labelsize=10)
def draw(i):
    v = V[i]; pc.set_array((v/np.abs(v).max()).ravel()); t = T[i]
    tl.set_text(f'$t={t:.1f}\\,w/C_\\mathrm{{Ai}}$')
    c = cr.get(round(float(t), 6))
    if c is not None and c['cP'] is not None:
        zj = np.array(c['zj']); zz = np.linspace(zj.min(), zj.max(), 10)
        cP.set_data(zj, c['ep']); cM.set_data(zj, -np.array(c['em']))
        fP.set_data(zz, np.polyval(c['cP'], zz)); fM.set_data(zz, -np.polyval(c['cM'], zz))
        al.set_text(f'$\\theta \\approx {c["th"]:.0f}$°')
    else:
        for h in (cP, cM, fP, fM): h.set_data([], [])
        al.set_text('')
if os.environ.get('FRAME'):                 # single frame as PNG, e.g. FRAME=40
    i = int(np.argmin(abs(T - float(os.environ['FRAME'])))); draw(i)
    fig.savefig(f'anim_frame_t{T[i]:.0f}.png', dpi=DPI); print('frame', T[i])
else:
    wr = FFMpegWriter(fps=FPS, codec='libx264', bitrate=-1, extra_args=['-pix_fmt', 'yuv420p', '-crf', '18', '-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2'])
    if FRAMES_DIR: os.makedirs(FRAMES_DIR, exist_ok=True)
    with wr.saving(fig, OUT, dpi=DPI):
        for i in idx:
            draw(i); wr.grab_frame()
            if FRAMES_DIR: fig.savefig(os.path.join(FRAMES_DIR, f'frame_t{T[i]:05.1f}.png'), dpi=DPI)
    print('wrote', OUT, len(idx), 'frames,', f'{len(idx)/FPS:.1f} s', 'size', fig.get_size_inches()*DPI)
