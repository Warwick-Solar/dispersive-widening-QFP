# Dispersive widening of a fast magnetoacoustic wave train: animation

Animation accompanying the paper *Dispersive widening of fast magnetoacoustic wave trains in the solar corona*
by D. Y. Kolotkov, E. Luxton and H. Shepherd (submitted to MNRAS). It shows the evolution of the wave train of Fig. 2
of the paper.

## What is shown

`wavetrain_animation.mp4` (18 s, 25 fps, 1472 × 336 px): transverse velocity $v_x$ of an impulsively generated fast
sausage wave train in a zero-β plasma slab of constant width, with density contrast $\rho_\mathrm{i}/\rho_\mathrm{e} = 10$
and profile steepness $p = 10$, from $t = 0$ to $90\,w/C_\mathrm{Ai}$ (one frame every $0.2\,w/C_\mathrm{Ai}$).

- Colour: $v_x$ normalised to its maximum in each frame, saturated at ±0.3; dashed lines mark the slab boundaries $x = \pm w$.
- From $t = 20\,w/C_\mathrm{Ai}$, once the wave train has formed: circles mark the $1/e$ transverse extents of the crests
  whose amplitude at $x = 0.5w$ exceeds 10% of the maximum, solid lines are linear fits to them, and the full opening
  angle θ is given in the lower-left corner (method of the paper, Sections 3.1–3.2).
- The frame at $t = 40\,w/C_\mathrm{Ai}$ corresponds to Fig. 2 of the paper. The opening angle decreases from about 44° at
  $t = 20$ to 24° at $t = 40$, 18° at $t = 60$ and 12° at $t = 90\,w/C_\mathrm{Ai}$.
- The animation stops at $t = 90\,w/C_\mathrm{Ai}$, before leaky waves reflected from the side boundaries of the
  computational domain enter the displayed region.

## Files

| file | content |
|---|---|
| `wavetrain_animation.mp4` | the animation |
| `solver2d.py` | 2D wave-equation solver (second-order finite differences, leapfrog time stepping, NumPy) |
| `angle_vs_t.py` | crest identification, $1/e$ transverse extents, linear fits and opening angle |
| `anim_data.py` | runs the simulation (domain $[-150, 150]\,w \times [0, 400]\,w$, ~5 min) and stores the displayed window and the fit results |
| `anim_render.py` | draws the frames as PNG files in `frames/` |

## Reproducing the animation

Requires Python 3 with NumPy, SciPy and Matplotlib, and ffmpeg.

```bash
python3 anim_data.py      # simulation -> anim_data.npz (~210 MB) + anim_crests.json
python3 anim_render.py    # frames/frame_t000.0.png ... frame_t090.0.png
cd frames && ffmpeg -y -framerate 25 -pattern_type glob -i 'frame_t*.png' \
  -c:v libx264 -crf 18 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -movflags +faststart \
  ../wavetrain_animation.mp4
```
