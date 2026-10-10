"""Animate NEMO vs Oceananigans on CANAL with the same numerics: animate_twin.py <nemo_run_dir> <oceananigans dump prefix> <out.mp4> [still.png]
   Both vorticities are computed here from the saved u, v with the same NEMO f-point stencil."""
import glob, sys
import netCDF4, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import LinearSegmentedColormap

run, jld, out = sys.argv[1:4]
still = sys.argv[4] if len(sys.argv) > 4 else None
f0, dx = 2 * 7.292115e-5 * np.sin(np.deg2rad(38.5)), 1e4
ink, ink2, gridc, surface = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
c_nemo, c_oc = "#2a78d6", "#eb6834"
plt.rcParams.update({"font.size": 9.5, "axes.edgecolor": ink2, "axes.labelcolor": ink2, "xtick.color": ink2,
                     "ytick.color": ink2, "axes.titlecolor": ink, "figure.facecolor": surface, "axes.facecolor": surface})
dye_cmap = LinearSegmentedColormap.from_list("dye", ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
dye_cmap.set_over("#eb6834"); dye_cmap.set_under("#e87ba4")
vor_cmap = LinearSegmentedColormap.from_list("vor", ["#1c5cab", "#86b6ef", "#e4e3df", "#f2a07e", "#b8431a"])

def zeta(U, V):
    """U (t, 99, 201) on wet U rows, V (t, 100, 201) on V rows incl. both walls -> ζ/f on the 98 interior f rows"""
    return ((np.roll(V[:, 1:-1], -1, axis=2) - V[:, 1:-1]) - (U[:, 1:] - U[:, :-1])) / dx / f0

# ---- NEMO: rows 0 and 100 are land; U at the east face, V at the north face of each T cell
nc = lambda g: netCDF4.Dataset(glob.glob(f"{run}/CANAL_*_grid_{g}_0000.nc")[0])
T = nc("T")
t_nemo = T["time_counter"][:] / 86400
dye_nemo = (T["vosaline"][:, 0] / T["vovvle3t"][:, 0])[:, 1:100].astype(float)
U_nemo = nc("U")["vozocrtx"][:, 0, 1:100].astype(float)
V_nemo = nc("V")["vomecrty"][:, 0, 0:100].astype(float)

# ---- Oceananigans (JLD2 dumped by jld2npy.jl, Float32, halo 5, Fortran order (x, y, t)): u at the west face
#      (-> roll by one to NEMO's east face), v on y faces 1..100
t_oc = np.fromfile(jld + "_t.bin") / 86400
raw = lambda n, ny: np.fromfile(jld + f"_{n}.bin", dtype=np.float32).astype(float).reshape(len(t_oc), -1, 211)[:, 5:5 + ny, 5:206]
U_oc, V_oc, dye_oc = np.roll(raw("u", 99), -1, axis=2), raw("v", 100), raw("c", 99)

# match the NEMO output times (6-hourly snapshots)
idx = [int(np.argmin(abs(t_oc - t))) for t in t_nemo]
assert max(abs(t_oc[idx] - t_nemo)) < 1e-6, "output times differ"
U_oc, V_oc, dye_oc = U_oc[idx], V_oc[idx], dye_oc[idx]

z_nemo, z_oc = zeta(U_nemo, V_nemo), zeta(U_oc, V_oc)
def ke(U, V): return 0.5 * ((U**2).sum(axis=(1, 2)) + (V**2).sum(axis=(1, 2)))
def ens(z):   return 0.5 * (z**2).sum(axis=(1, 2))
ke_n, ke_o = ke(U_nemo, V_nemo), ke(U_oc, V_oc)
en_n, en_o = ens(z_nemo), ens(z_oc)
ke_n, ke_o, en_n, en_o = ke_n / ke_n[0], ke_o / ke_n[0], en_n / en_n[0], en_o / en_n[0]
days = t_nemo
for d in (10, 20, 30, 45, 60):
    n = int(np.argmin(abs(days - d)))
    print(f"day {days[n]:4.1f}  KE NEMO {ke_n[n]:.3f}  Oceananigans {ke_o[n]:.3f}   enstrophy {en_n[n]:.3f} {en_o[n]:.3f}   "
          f"dye range NEMO {np.nanmin(dye_nemo[n]):.2f}..{np.nanmax(dye_nemo[n]):.2f}  Oc {dye_oc[n].min():.2f}..{dye_oc[n].max():.2f}")

ext_t = [-1005, 1005, -495, 495]                    # T cells (km)
ext_f = [-1000, 1010, -490, 490]                    # interior f points (km), cell-centred display
fig = plt.figure(figsize=(13, 9.6))
gs = fig.add_gridspec(3, 5, height_ratios=[1, 1, 0.78], width_ratios=[1, 0.022, 0.13, 1, 0.022],
                      hspace=0.34, wspace=0.03, left=0.055, right=0.93, top=0.86, bottom=0.065)
axs = [[fig.add_subplot(gs[r, c]) for c in (0, 3)] for r in range(2)]
ims = []
for r, (lab, z, dy) in enumerate([("NEMO", z_nemo, dye_nemo), ("Oceananigans", z_oc, dye_oc)]):
    a, b = axs[r]
    ims.append((a.imshow(z[0], origin="lower", extent=ext_f, cmap=vor_cmap, vmin=-0.8, vmax=0.8, interpolation="nearest"),
                b.imshow(dy[0], origin="lower", extent=ext_t, cmap=dye_cmap, vmin=29.99, vmax=35.01, interpolation="nearest")))
    a.set_title(f"{lab}   relative vorticity ζ/f", loc="left", fontsize=10, fontweight="bold")
    b.set_title(f"{lab}   dye", loc="left", fontsize=10, fontweight="bold")
    a.set_ylabel("y (km)"); b.set_yticklabels([])
    for ax in (a, b):
        ax.set_aspect("auto"); ax.tick_params(length=2); ax.set_xticks([-1000, -500, 0, 500, 1000])
        if r == 1: ax.set_xlabel("x (km)")
        else: ax.set_xticklabels([])
cb = fig.colorbar(ims[0][0], cax=fig.add_subplot(gs[0:2, 1]), extend="both"); cb.set_label("ζ / f")
cb = fig.colorbar(ims[0][1], cax=fig.add_subplot(gs[0:2, 4]), extend="both")
cb.set_label("dye (psu)   orange > 35.01, pink < 29.99")

tk, te = fig.add_subplot(gs[2, 0:2]), fig.add_subplot(gs[2, 3:5])
for ax, (yn, yo), ttl in [(tk, (ke_n, ke_o), "Kinetic energy (relative to day 0.25) — inviscid, should not grow"),
                          (te, (en_n, en_o), "Enstrophy (relative to day 0.25)")]:
    ax.plot(days, yn, color=c_nemo, lw=1.0, alpha=0.3); ax.plot(days, yo, color=c_oc, lw=1.0, alpha=0.3, ls="--")
    ax.set_title(ttl, loc="left", fontsize=10)
    ax.set_xlim(0, 60); ax.set_xlabel("time (days)")
    ax.grid(color=gridc, lw=0.6); ax.spines[["top", "right"]].set_visible(False); ax.tick_params(length=2)
te.set_yscale("log")
tk.set_ylim(min(ke_n.min(), ke_o.min()) - 0.02, max(ke_n.max(), ke_o.max()) + 0.02)
ln_n, = tk.plot([], [], color=c_nemo, lw=2, label="NEMO"); ln_o, = tk.plot([], [], color=c_oc, lw=2, ls="--", label="Oceananigans")
le_n, = te.plot([], [], color=c_nemo, lw=2); le_o, = te.plot([], [], color=c_oc, lw=2, ls="--")
dots = [ax.plot([], [], "o", ms=5, color=c, mec=surface, mew=1.5)[0] for ax in (tk, te) for c in (c_nemo, c_oc)]
tk.legend(frameon=False, fontsize=9, loc="lower left")
fig.text(0.055, 0.965, "CANAL barotropic jet — NEMO vs Oceananigans with the same numerics", fontsize=13, color=ink, fontweight="bold")
fig.text(0.055, 0.93, "Both: WENO5 vector-invariant momentum (vorticity, KE gradient, vertical advection) · WENO5 tracer \n"
         "RK3 with split-explicit free surface, Δt = 1350 s · 10 km grid (201 periodic × 99 wet cells), f-plane, inviscid · identical initial state, including the random noise", fontsize=9.5, color=ink2, va="top")
corr = [np.corrcoef(a.ravel(), b.ravel())[0, 1] for a, b in zip(z_nemo, z_oc)]
clock = fig.text(0.93, 0.965, "", fontsize=13, color=ink, ha="right", fontweight="bold")

cortxt = fig.text(0.93, 0.938, "", fontsize=9.5, color=ink2, ha="right", va="top")

def draw(n):
    for (ia, ib), z, dy in zip(ims, (z_nemo, z_oc), (dye_nemo, dye_oc)):
        ia.set_data(z[n]); ib.set_data(dy[n])
    s = slice(0, n + 1)
    ln_n.set_data(days[s], ke_n[s]); ln_o.set_data(days[s], ke_o[s]); le_n.set_data(days[s], en_n[s]); le_o.set_data(days[s], en_o[s])
    for d, y in zip(dots, (ke_n, ke_o, en_n, en_o)): d.set_data([days[n]], [y[n]])
    clock.set_text(f"day {days[n]:5.2f}"); cortxt.set_text(f"pointwise correlation of ζ: {corr[n]:.3f}")

writer = FFMpegWriter(fps=12, bitrate=8000, extra_args=["-pix_fmt", "yuv420p"])
with writer.saving(fig, out, dpi=120):
    for n in range(len(days)):
        draw(n); writer.grab_frame()
if still:
    draw(int(np.argmin(abs(days - 6)))); fig.savefig(still, dpi=110)
print("frames:", len(days), "->", out)
