# Honeycomb V2

## Bugs and status

1. At the moment I added Evolution Operator for fix nf. Next step is to have a collection of these that evolves from initial scale to first heavy quark, from first heavy quark to the second etc until final scale. Next to next step is store these for each RungeKutta step and build interpolation table, so that I can decide on the grid once, precompute all evolution operators on the Q grid and then use them.
2. The gluon plus has some weird ridges, present in both current code and old code equally. They are present along each of the three $x_i=0$ lines, but evident only if zoomed in, if zoomed out they are very difficult to see. To be investigated

## Unsupported feature

1. It is not yet possible to obtain from the Solution struct the discretized values as vectors on the grid for $T, \Delta T, T_F^\pm$. It is only possible to rotate between evolution basis (gluon, singlet, NSi) to physical basis in terms of $\mathfrak{S}^\pm_f, \mathfrak{F}^\pm$ (the latter two for the gluon are the same between evolution and physical) 