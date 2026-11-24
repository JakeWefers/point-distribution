# Substitution data provenance

`substitution.toml` contains geometric placement data for the five marked
Schlottmann square–triangle prototiles, with inflation `2 + sqrt(3)`.

Sources:

1. D. Frettlöh, *A fractal fundamental domain with 12-fold symmetry*, Figure 1:
   https://www.math.uni-bielefeld.de/~frettloe/papers/12fdwin-arx.pdf
2. Tilings Encyclopedia, square–triangle rule diagram:
   https://tilings.math.uni-bielefeld.de/img/substitution/square-triangle/rule.png
   https://tilings.math.uni-bielefeld.de/substitution/square-triangle/

## Representation

`base` gives the centered unit-edge polygon for each of five marked types.
Each `rules` child record supplies its type, center relative to the inflated
parent's center, and a 2×2 row-major orthogonal orientation matrix. Coordinates
use screen convention (positive y downward); reflections are intentional.

For parent center `c`, parent orientation `O`, child center `q`, and child
orientation `A`, the new tile has:

- center `(2 + sqrt(3))*c + O*q`
- orientation `O*A`

Rule lengths are 13, 16, 14, 29, and 29 children. Replacement patches protrude
past the original inflated polygon, so coincident tiles must be deduplicated.
Simply pasting those patches without deduplication would overcount vertices.

## Recovery and cross-checks

The geometric coordinates were recovered from vector polygons in Figure 1.
Rounded source edges were snapped to the twelve unit directions at multiples
of 30 degrees. Integer coefficients of 1 and sqrt(3) were propagated through
each edge graph, with all cycles checked for agreement. Local parent centers
were fixed using the central child tiles in the published rules.

The orientation markings were cross-checked against the encyclopedia diagram.
The prototype-5 marker template in the PDF was reflected across its diagonal
to match the encyclopedia convention; the central white square in the
square-type-4 replacement uses a 180-degree orientation. The unmarked corner
fixes the distinguished direction of triangle type 3. Equivalent symmetric
orientations are accepted when they induce the same next replacement geometry.
These conventions are baked into the data; runtime image analysis is unnecessary.

The data are stored as Float64 values, not exact algebraic numbers. Shared
geometry is keyed at 1e-7 edge units and vertices normalized at 1e-9 edge units.
The included four-generation regression tests check the resulting coordinates
against the approved image and independently verify spacing and coverage.
This is numerical verification of an implementation of the published construction,
not a new proof of aperiodicity.

From seed type 4, generation tile counts are:

| Generation | Unique tiles | Merged duplicate child copies |
|---|---:|---:|
| 1 | 29 | 0 |
| 2 | 428 | 110 |
| 3 | 6,027 | 1,830 |
| 4 | 84,198 | 26,614 |

Generation 4 has 55,321 unique vertices. The default 10,000-point crop is
strictly inside the generated patch, has full square coverage, and has uniform
100-pixel noiseless nearest-neighbor distances within floating-point tolerance.
