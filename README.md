# Honeycomb V2

## Bugs and status

1. Fixed bug in PushFlavor, now there is still something wrong in the singlet on the line $x_1=0$

## Unsupported feature

1. It is not yet possible to obtain from the Solution struct the discretized values as vectors on the grid for $T, \Delta T, T_F^\pm$. It is only possible to rotate between evolution basis (gluon, singlet, NSi) to physical basis in terms of $\mathfrak{S}^\pm_f, \mathfrak{F}^\pm$ (the latter two for the gluon are the same between evolution and physical) 