# example/umbrella_pull - OrganL 0.9b pre-release

This is an illustation of how to perform a global variable scan. 


## HOW TO RUN:

1. Copy the appropriate **organl** binary for your platform from `/binaries/` into this folder.
2. Run on commandline: `[organl] oblate.obj`


## Brief Description:

- This is an explicit example of optimizing values of various parameters.
- Uses the `DCH` (default) energy functional.


## Important Parameters:
> Since there is a high degree of non-linearity in the stimulation, a precise behavioural evaluation of parametric values is not possible. The following may be viewed as a rough baseline to work on. Read the Manual for a detailed take on running stimulations and optimising parametric values.

- `S updateNghbrsEvery XX` - defines the number of runs before the Neighbour list update. Also check out `S nghbRadius XX` which is the neighbourhood distance. Too low `updateNghbrsEvery` or too high `nghbRadius` value will result in significant computational time. Neighbourhood evaluation is a known throttlepoint for the current implementation.
- `S updateDefEvery XX` - sets the number of runs before the target `E DEF X XX` values are updated. Works together with `E INCR X YY` which allows for a gradual change in `E DEF X` value at each DEF update.
- `S writeEvery XX` - sets the number of steps after which the output is recorded as a `.vtu` file for visualisation and `.obj` files for continuing stimulation from that point.
- `S RFile [name.log]` - sets the name of the Report file. Very useful for troubleshooting, refer the manual for details on interpreting the contained data. To be used in combination with `S RFreq XX` which sets the steps after which the report line is appended to `[name.log]`.

## Other relevant Parameters:
- `E DEF 0 XX` - sets the value of $A_0$.
- `E DEF 1 XX` - sets the value of $V_0$.
- `E LAMBDA 0 XX` used to denote the constraint on area (Surface Tension).
- `E LAMBDA 1 XX` used to denote the constraint on volume (pressure).