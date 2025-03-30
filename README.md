# Honeycomb V2

### To Do:
- [ ] 
- [ ] At the moment Discretization has mixed convention for radius-angle ordering (specifically the flatten vector has the slower loop over the angle.) Probably must be uniformed. Everything else works with radius, angle as order for arguments/returned types.

### Notes:

1. **30.03.25** At the moment, the `Discretization` struct is just a wrapper for the various interpolation strategies and a container to store a reference to the `Grid2D` object. The discretized values of the various functions are obtainable by calling the `operator()` of the `Discretization` instance. I.e., one should have one global `Grid2D` object, one global `Discretization` object and a number of `Eigen::VectorXd` variables. The latter will be eventually wrapped in the `Solution` struct, which will take care of a number of tasks. Primarily, it will store the solution of the evolution problem and provide all the necessary operations to be able to pass it to Runge-Kutta along with the kernels. The `Solution` I think will store a pointer to the `Discretization` object, so that we may initialize the models and, depending on what I decide to do, handle directly the interpolation of the solution.