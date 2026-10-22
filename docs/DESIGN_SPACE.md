# Design space (version 1)

Rule: version 1 stays inside the ranges validated against Abaqus. Every extension (new geometry
ranges, materials) needs blind Abaqus points at random positions inside the new range before it is used
in the paper.

## Omega-stiffened panel
| parameter | range / values | validated by |
|---|---|---|
| stringer pitch b_s | 150 – 285.5 mm | 90 Abaqus cases (60 symmetry, 30 periodic) |
| aspect ratio β = L / b_s | 0.75 – 3.0 | same |
| stringer: b3, h, web angle | 27.75 mm, 30 mm, 25° (fixed) | same |
| material | CFRP: E1 150, E2 9, G12 = G13 = G23 5.3 GPa, ν12 0.33 | same |
| skin / stringer | 8 / 7 plies × 0.125 mm | same |
| ply angles | {0, ±30, ±45, ±60, 90}, 50 % symmetric lay-ups | contains all validated lay-ups |
| end conditions | SS, CC | SS and CC (symmetry); SS (periodic) |
| side edges | symmetry, periodic | periodic CC not yet compared with Abaqus |

## Members (C, Z, L, T, I, RHS)
| parameter | range / values | validated by |
|---|---|---|
| web H, flange B | 50 mm, 20 mm (fixed) | 144 Abaqus cases |
| length L | 200 – 1000 mm | same |
| material | CFRP: E1 171.4, E2 9.08, G12 = G13 5.29, G23 2.645 GPa, ν12 0.32 | same |
| lay-up | 4 plies × 0.25 mm, angles from the same set | Sym1, Sym2, Unsym validated |
| end conditions | SS, CC | same |

Hat and double-hat are excluded from version 1 (GBT not compared with Abaqus, 2–5 GB memory).

## Network inputs
Scaled continuous parameters (to [-1, 1]), categorical flags, and for every laminate its 12 lamination
parameters ξ_A, ξ_B, ξ_D (Tsai–Pagano) and its thickness. The input size is therefore fixed whatever the
stacking sequence, and any lay-up — including ones never sampled — maps into the same bounded space.

## Blind Abaqus points needed for the paper
About 20–30 per family at random (Sobol) positions inside the space above, with random lay-ups.
The existing 234 cases lie on the grid of the parametric studies and stay as a second test set.
