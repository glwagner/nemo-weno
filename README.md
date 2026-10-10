# nemo-weno: WENO advection for NEMO

This fork of [NEMO](https://forge.nemo-ocean.eu/nemo/nemo) adds WENO advection, ported from [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl):

- WENO3 and WENO5 tracer advection.
- A WENO vector-invariant momentum scheme, matching Oceananigans' `WENOVectorInvariant`.
- An RK3 time-stepping fix, needed to keep the momentum scheme stable.

The upstream NEMO README is [README.rst](README.rst).

The work is a stack of five pull requests:

| PR | branch → base | content |
|---|---|---|
| [#1](https://github.com/glwagner/nemo-weno/pull/1) | `nn-hls-3` → `main` | 3-point halo (`nn_hls = 3`) |
| [#2](https://github.com/glwagner/nemo-weno/pull/2) | `weno3-tracer-advection` → `main` | WENO3 tracer advection |
| [#3](https://github.com/glwagner/nemo-weno/pull/3) | `rk3-bt-forcing-ab3` → `main` | Adams–Bashforth barotropic forcing for RK3 (`nn_bt_frc`) |
| [#4](https://github.com/glwagner/nemo-weno/pull/4) | `weno5-tracer-advection` → `weno3-tracer-advection` | WENO5 tracer advection (merges #1) |
| [#5](https://github.com/glwagner/nemo-weno/pull/5) | `weno-vector-invariant` → `weno5-tracer-advection` | WENO vector-invariant momentum (draft) |

Figures and animations are on the [`pr-assets`](https://github.com/glwagner/nemo-weno/tree/pr-assets) branch.

![CANAL barotropic jet, NEMO vs Oceananigans with the same numerics](https://github.com/glwagner/nemo-weno/blob/pr-assets/canal-vi-twin/canal_nemo_vs_oceananigans.gif?raw=true)

## 1. WENO3 tracer advection (#2)

- **New module `src/OCE/TRA/traadv_weno.F90`:** a flux-form WENO-Z scheme ported from Oceananigans `WENO(order=3)`, with smoothness indicators β and ε = 1e-8. Within two cells of a wall it falls back to first-order upwind. A centred fallback undershot, so upwind replaced it.
- **Wired into `traadv.F90`**, with new entries in `&namtra_adv`:
  - `ln_traadv_weno` turns the scheme on.
  - `nn_weno_h` and `nn_weno_v` set the horizontal and vertical order (3 or 5).
- **`tests/SOLITON/EXPREF/namelist_ref`:** SOLITON keeps its own copy of `namelist_ref`, so it needed the same entries.
- **`tests/CANAL/MY_SRC/usrdef_istate.F90`:** `RANDOM_SEED` is now fixed at 1234567. The initial noise used to be unseeded, so runs weren't reproducible.
- **Tested on:** the CANAL barotropic jet (`nn_initcase = 1`, `ln_sshnoise`, passive dye 35 psu inside the jet / 30 outside, `rn_a0 = 0`). It was compared with CEN, FCT, MUS, UBS and QCK.

## 2. 3-point halo (#1)

- WENO5 needs a 3-cell stencil, so `nn_hls = 3` had to work.
- **`lib_mpp.F90`:** `n_hlsmax` raised from 2 to 3.
- **`mppini.F90`:** stops if `nn_hls > n_hlsmax`.
- **Arrays hard-coded to a 2-point halo, now sized with `nn_hls`:**
  - `zdfdrg.F90`: `A2D(2)` → `A2D(nn_hls)` dummies.
  - `zdfmxl.F90`: `nmln`, `hmlp`.
  - `icbutl.F90`: `zssh_lead_m`.
  - `traadv_fct.F90`: tridiagonal solver, `T2D(nn_hls)`.
- `sette/sette.sh` updated to match.
- **Verified with:** `run.stat` / `tracer.stat` bitwise comparisons across decompositions, using a reproducible build arch with `-ffp-contract=off`. Debug builds used `-fcheck=all` and a signalling-NaN malloc interposer.

## 3. WENO5 tracer advection (#4)

- **New shared include `src/OCE/weno_kernels.h90`:** `weno5`, `weno3`, the β and weight functions, and the reconstruction kernels.
- **Order 5** uses three Jiang–Shu smoothness indicators, as in Oceananigans. Near walls the order drops 5 → 3 → 1, controlled by masks `zm3` and `zm5`.
- **Requires** `nn_weno_h = 5` together with `nn_hls = 3`.

## 4. WENO vector-invariant momentum (#5, draft)

This ports Oceananigans `WENOVectorInvariant` with `OnlySelfUpwinding`. It needs `ln_dynadv_vec = .true.`.

- **Relative vorticity (`dynvor.F90`, new `vor_weno`):**
  - Namelist: `ln_dynvor_weno` in `&namdyn_vor`, with `nn_dynvor_weno` = 3 or 5.
  - Planetary vorticity stays with the selected EEN/ENE scheme (`ivor = ncor`). The WENO reconstruction handles only the relative vorticity, using Oceananigans' `VelocityStencil` smoothness.
- **Kinetic-energy gradient and vertical advection:**
  - Namelist: `ln_dynadv_weno` and `nn_dynadv_weno` in `&namdyn_adv`.
  - New routines `dyn_keg_weno` (`dynkeg.F90`) and `dyn_zad_weno` (`dynzad.F90`). The dispatch is in `dynadv.F90`.
  - `stp2d.F90` and `stprk3_stg.F90` exchange the vertical velocity `ww` over 2 halo points.
- **Verified with:** a kernel harness that matches Oceananigans to round-off (~1e-15) for every term. Results are bit-identical across MPI decompositions, and existing options are unchanged.

## 5. The instability, and the time-stepping fix (#3)

- **Symptom:** the vector-invariant runs blew up at both 1350 s and 675 s. Plain EEN also blew up (day 6 at 1350 s, day 15 at 675 s), not only the WENO variants. Oceananigans with the same schemes stayed stable.
- **Root cause:** four parallel investigations converged on the same answer. NEMO's RK3 (`stprk3.F90`) advances the depth-mean momentum once per step, using slow forcing frozen at time n: a forward-Euler step. Weakly dissipative vector-invariant advection injects energy through that step.
- **Fix:** an Adams–Bashforth extrapolation of the slow forcing, `nn_bt_frc` = 1, 2 or 3 in `&namdyn_spg` (`dynspg.F90`, `dom_oce.F90`).
  - `nn_bt_frc = 1` (the default) is bit-identical to `main`.
  - **First version, extrapolating all the forcing:** ORCA2 blew up at step 9.
  - **Extrapolating everything except the dissipative terms:** ORCA2 blew up at step 122. The planetary-vorticity residual couples columns inertially at depth jumps (fΔt ≈ 0.8 in ORCA2), first in a 17 m column of the Nile delta.
  - **Final version:** extrapolate only the depth mean of the momentum advection. When `nn_bt_frc > 1`, `stp2d.F90` calls `dyn_vor` twice, Coriolis (`ncor`) and then relative vorticity (`nrvm`), to isolate that part. In `dynspg_ts.F90`, `ts_frc_ab` combines it with the stored previous forcings `ufrc_b1`/`ufrc_b2`. Those are written to restarts along with `nfrc_ab`. `dynvor.F90` gates the Coriolis calls with `ll_cor`.
- **Validated on:**
  - ORCA2 with sea ice and PISCES runs 3 months stable and stays close to the default (max SSH 1.88093 vs 1.88100 m).
  - CANAL with EEN conserves kinetic energy (1.006 of its initial value at day 60), and UBS gives 0.967.
  - Restarts are exact, and results are bit-identical across decompositions for CANAL, BENCH and ORCA2.
  - The VORTEX and AGRIF tests pass.
- **Integration fix in #5:** with `nn_bt_frc > 1`, `dyn_vor_RK3` would have used EEN's relative vorticity instead of WENO's.

## 6. Final comparison against Oceananigans

- **Configuration:**
  - `ln_dynvor_weno = nn_dynvor_weno = 5`
  - `ln_dynadv_weno = nn_dynadv_weno = 5`
  - `ln_traadv_weno` with `nn_weno_h = nn_weno_v = 5`
  - `nn_hls = 3`, `nn_bt_frc = 3`, `rn_Dt = 1350`
- **Results:** stable for 60 days, and the same at 675 s. With `nn_bt_frc = 1` it blows up at step 241.
- **Matching twin ([`canal_twin5.jl`](https://github.com/glwagner/nemo-weno/blob/pr-assets/canal-vi-twin/canal_twin5.jl)):** NEMO's exact geometry (201 periodic × 99 wet cells, walls at ±495 km) and the gfortran noise field reproduced by [`noise.f90`](https://github.com/glwagner/nemo-weno/blob/pr-assets/canal-vi-twin/noise.f90).
- **Agreement:**
  - Vorticity correlation is 0.997 at 6 h, 0.95 at day 6 and 0.48 at day 30, after which the chaotic eddy field decorrelates.
  - Kinetic energy at day 60 is 0.953 (NEMO) vs 0.963 (Oceananigans).

**Still open on #5:** a multi-level test of the vertical term (BENCH or ORCA2), a cost comparison with EEN, tiling and vortex-force support, and a rebase once #1–#4 merge.
