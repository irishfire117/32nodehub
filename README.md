# 1D Steady Conduction with Spatially Varying Generation

Finite-difference solver for a one-dimensional rod with internal heat generation,
solved as a tridiagonal matrix system with Gauss–Seidel iteration.

ME452 Heat Transfer — 70-pt problem. C++17, no dependencies beyond the standard
library.

## The problem

A rod of unit length (0 <= x <= 1) in steady state:

- Left end (x = 0): held at **T = 100 K** (fixed temperature, Dirichlet)
- Right end (x = 1): exposed to ambient fluid at **T∞ = 25 K** with heat
  transfer coefficient **h = 100 W/(m²·K)** (convection, Robin)
- Electrical heater inside the rod generates **q̇(x) = 10x + 5 W/m³**
  (generation *varies with position* — this is what kills a hand solution
  worth checking against, and why we go numerical)
- Thermal conductivity **k = 1 W/(m·K)** (unit conductivity)

Governing equation (steady 1D conduction with generation):

    k · d²T/dx² + q̇(x) = 0,     0 < x < 1

with the two boundary conditions:

    T(0) = 100
    -k · dT/dx (at x=1) = h · (T_surface - T∞)

## Method: cell-centered finite differences with ghost cells

The rod is cut into **N blocks** of width Δx = L/N. The temperature unknown is
at the **center of each block**, not on the boundary:

    x_i = (i - 0.5)·Δx,     i = 1 .. N

so x_1 = Δx/2 (not 0) and x_N = L - Δx/2 (not L). Each physical face (x = 0 and
x = L) sits **halfway between a real block and a fictitious "ghost" cell** just
outside the domain.

Discretizing the heat equation at block i (central difference for the second
derivative) gives the interior equation:

    T[i-1] - 2·T[i] + T[i+1] = -q̇_i · Δx²/k

For i = 1 and i = N the neighbor is a ghost cell, which doesn't exist — so the
boundary condition is used to **express the ghost in terms of real temperatures**
and substitute it in. Every boundary condition reduces to the same form:

    T_ghost = a · T_edge + b

**Left face (Dirichlet).** The face temperature is *given* (100 K) and the face
is halfway between the ghost and block 1, so linear interpolation:
100 = (T_ghost + T_1)/2  →  **T_ghost = 2·100 − T_1**. Substituting into
block 1's equation merges the ghost's −T_1 into the diagonal:

    3·T_1 − T_2 = q̇_1·Δx²/k + 2·100

The diagonal going from 2 to **3** is the Dirichlet fingerprint.

**Right face (convection).** The face temperature is *unknown*. The face has no
volume, so whatever conducts in must convect out:

    -k · dT/dx|_L = h · (T_face − T∞)

Two first-order approximations at the face (slope from ghost to block N, face
temperature as their average) turn this into a solvable relation for the ghost:

    T_ghost = b·T_N + C

    b = (k/Δx − h/2) / (k/Δx + h/2)
    C = h·T∞ / (k/Δx + h/2)

Substituting into block N's equation gives:

    −T[N−1] + (2 − b)·T_N = q̇_N·Δx²/k + C

Sanity limits: as h → ∞, b → −1 and the row becomes the Dirichlet form pinned
at T∞ (the "3" fingerprint); as h → 0, b → 1 and the row becomes the insulated
wall (diagonal 1).

## The matrix system

All N equations assembled: a **tridiagonal** matrix — each row touches only its
own diagonal and the two off-diagonals. Everything is built by a loop plus two
corner patches:

    interior rows:    [ -1   2  -1 ]      rhs = q̇_i·Δx²/k
    row 1 (Dirichlet):[  3  -1   0 ]      rhs = q̇_1·Δx²/k + 2·Tw
    row N (convective):[ 0  -1  2-b ]     rhs = q̇_N·Δx²/k + C

For N = 10 (Δx = 0.1, so b = −2/3, 2−b = 8/3 ≈ 2.667, C = 41.667):

    row 1:   3   -1    0   0   0   0   0   0   0   0   | 200.055
    row 2:  -1    2   -1    0   0   0   0   0   0   0   |   0.065
    ...  (rows 3..9: the -1 / 2 / -1 band)  ...
    row 10:  0   0    0   0   0   0   0   0  -1  2.667  |  41.812

The matrix is almost symmetric; the only asymmetric cells are the two corners,
and they are asymmetric because the two walls have different physics.

## Solver: Gauss–Seidel iteration

Each row, solved for its own unknown, is a recipe in terms of its neighbors.
The iteration: guess all temperatures, run every recipe, repeat until nothing
moves.

    T_new[i] = (T[i-1] + T[i+1] + q̇_i·Δx²/k) / 2        (interior)
    T_new[0] = (2·Tw + T[1] + q̇_0·Δx²/k) / 3            (Dirichlet end)
    T_new[N-1] = (T[N-2] + C + q̇_N·Δx²/k) / (2 - b)     (convective end)

Gauss–Seidel (vs. Jacobi): the sweep overwrites each T[i] **in place**, so row
i's recipe already sees row i−1's *fresh* value from the same pass. Information
from the hot wall propagates rightward within one sweep instead of waiting a
full sweep per cell. Convergence is guaranteed here because the matrix is
diagonally dominant (every diagonal entry outweighs the sum of its row's
off-diagonals).

Stopping rule: track the maximum absolute update in each sweep; stop when it
drops below 1e-10 K (with a 200,000-sweep safety cap in case a broken system
never converges). When the iteration stops, every T simultaneously satisfies
its own row — which is the definition of the solution. An independent audit
re-evaluates all rows: `max |A·T − rhs|`.

Why iterative at all? For a tridiagonal system a direct O(N) Thomas algorithm
exists (and was used as the ground-truth cross-check); for large or
non-tridiagonal systems (2D/3D grids, millions of nodes) only the iterative
family is affordable.

## Results

**N = 32** (verification run, cross-checked against an independent
Thomas-algorithm solve of the identical system — both agree to solver
tolerance):

    Gauss–Seidel sweeps to 1e-10:  2346
    T_1  (x = 0.0156):  98.906 K
    T_32 (x = 0.9844):  27.051 K
    surface T(L):       25.800 K
    max |A·T − rhs|:    ~1e-10

**N = 10** (the 10-node answer):

    node   x      T (K)
    1     0.05    96.499
    2     0.15    89.441
    3     0.25    82.319
    4     0.35    75.121
    5     0.45    67.839
    6     0.55    60.461
    7     0.65    52.979
    8     0.75    45.381
    9     0.85    37.659
    10    0.95    29.802

    Gauss–Seidel sweeps:  251
    max |A·T − rhs|:      9.3e-11
    surface T(L):         25.800 K
    energy balance:  convective out = 80.025 W/m²
                     wall in + generation = 70.025 + 10.000 = 80.025 W/m²  ✓

**Reading the numbers.** T_10 = 29.80 K here vs 27.05 K at N=32 — not a
discrepancy but grid resolution: at N=10 the last block spans x = 0.9–1.0
(center 0.95), at N=32 it spans 0.97–1.0 (center 0.984). Coarser grid, center
further from the cold face, reads hotter. The *physical* quantity — the surface
temperature — agrees to four decimals (25.800) at both resolutions, so the
surface heat flux q″ = h(T_s − T∞) ≈ **80 W/m²** is converged to 4
significant figures.

The energy balance is the physics-level check: heat leaving by convection at
the right face must equal heat conducted in at the left wall plus total
generation (∫₀¹ (10x+5) dx = 10 W/m², analytic). It closes to the printed
digit.

## Build and run

    g++ -g -Wall 32nodeblock.cpp -o 32nodeblock
    ./32nodeblock

Memory-safety check (recommended while editing):

    g++ -g -fsanitize=address 32nodeblock.cpp -o 32nodeblock_asan
    ./32nodeblock_asan        # any out-of-bounds access aborts loudly

Change `N` (one line) to run other resolutions; nothing else depends on it.
VS Code: `Ctrl+Shift+B` builds (g++, -g -Wall), `Ctrl+Alt+R` builds and runs,
`F5` debugs with breakpoints.

## Files

    32nodeblock.cpp    the solver (this repo's source)
    32nodeblock        compiled binary (build artifact — rebuild, don't edit)
    .vscode/           build/run/debug tasks for VS Code
    Screenshot_*.png   the original problem statement
