# Honeycomb V2

## Bugs and status

1. Fixed bug in PushFlavor, now there is still something wrong in the singlet on the line $x_1=0$
   1. This is checked in the `convolution_check`, where already a single convolution with very trivial test functions show the problem. Fix, there was a factor 2 too much for $x_1=0 / x_3=0$ cases of Vplus13, derived from taking the sum of the two integrals from the thetas, without considering that I should use consistent $\theta$ definition (e.g. $\theta(0)=1/2$)

## Unsupported feature

1. It is not yet possible to obtain from the Solution struct the discretized values as vectors on the grid for $T, \Delta T, T_F^\pm$. It is only possible to rotate between evolution basis (gluon, singlet, NSi) to physical basis in terms of $\mathfrak{S}^\pm_f, \mathfrak{F}^\pm$ (the latter two for the gluon are the same between evolution and physical) 