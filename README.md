# Point Distribution

This project aims to produce a non-repeating pattern of homogeounsly space points
on the surface of a cylinder of given size.
P denoting the set of points whose positions we determine, we define a constellation A to be more homogonous than a constellation B, if

max_{a_A, b_A \in P u BoundarySet} dist(a_A, b_A) < max_{a_B, b_B \in P u BoundarySet} dist(a_B, b_B)

Where dist() is the euclidean distance and a_A denotes the position of point a in constellation A. The Boundary Set denotes the edges of the cylinder we are placing the dots on.

This problem consists of finding a mix between non-repeatable spacing of points so as to "not see any repetitions of constellations" and the equidistant positiong.

This attempt does it in the following way.

1. A  semi a-periodic tiling is generated to fit on the surface of the cylinder so that  the
number of nodes of the tiles estimates closely the number of points we want to place.
It is only semi a-peridic as we of course want to wrap the tiling around the cylinders circumferance. It is therefore periodic in this dimension.

2. We then apply noise the each of these nodes.



This project began as a Python experiment and was later implemented in Julia. This repository contains the Julia algorithms and command-line examples. The square–triangle constructions are based on published tiling methods; the references are credited below.

This porject began by considering the same problem but on bounded rectangluar areas on which actual aperiodic tilings can be applied.

Later on, to tile the cylinder, square trianlge tilings were used as all nodes are equidistant in these.

A spin-off of a tiing by Martin Schlottmann was used.
See example below.

![Cylinder point pattern with bounded displacement](example.png)

## Quick start

Install Julia (tested with Julia 1.10.3), then run this from the repository folder:

```sh
julia example.jl
```

The example generates the included `example.png`. Rerunning it replaces that image using the same seeds. It does not create additional coordinate or statistics files. PNG generation uses libraries bundled with Julia; no package installation or project setup is needed.

Edit the settings at the top of `example.jl` to change the cylinder dimensions, point density, noise, dot radius or image resolution.

## Files

| File | Purpose |
| --- | --- |
| `cylinder/1_tiling_coordinates.jl` | Generate square–triangle tiling vertices with a matching cylinder circumference. |
| `cylinder/2_add_noise.jl` | Apply bounded random displacement and wrap points around the circumference. |
| `cylinder/3_output_dots.jl` | Render the dots to PNG or TIFF, including dots crossing the cylinder join. |
| `example.jl` | Run the three steps with a small reproducible example. |
| `example.png` | The generated result shown above. |
| `README.md` | Explanation and usage. |

## Example settings and units

The example uses a circumference of 0.20 m, a height of 0.12 m, a requested density of 192.72 points per square inch, a maximum displacement of 7% of the tile-edge spacing, and a resolution of 300 DPI. The dot radius is 15% of the tile-edge spacing.

Cylinder dimensions are supplied in meters. Returned coordinates, spacing, dot radius and optional cutout dimensions are in millimeters. The first return value from `tiling_coordinates` is an N × 2 matrix of point coordinates.

## Noise and the cylinder join

For spacing `d` and noise percentage `p`, the displacement cap is `R = p*d/100`. Each point receives an independent displacement uniformly sampled over a disk:

```text
r = R * sqrt(u)
theta = 2*pi*v
delta = (r*cos(theta), r*sin(theta))
```

Here `u` and `v` are uniform random draws in [0, 1). Circumferential coordinates wrap across the left and right edges. Points may move beyond the top and bottom of the rectangle; the renderer clips them to the image boundary.

`add_noise` returns a conservative relative distance-change bound of `2R/d`. For 7% displacement, this bound is 14% of the original spacing. It is a geometric bound, not the measured error of a particular generated pattern.

## Image output

PNG output works with Julia alone. TIFF output additionally requires the `vips` command-line tool from libvips. To produce TIFF, install libvips and change the example output filename to end in `.tif`. Large TIFFs need a BigTIFF-capable viewer.

The output logic in `3_output_dots.jl` was written by Codex and checked by running examples and reviewing the calculations on paper.

## Background

Background on the original square–triangle construction: [Tilings Encyclopedia: square–triangle tiling](https://tilings.math.uni-bielefeld.de/substitution/square-triangle/).
