# example/budding_couple/bud_enhance_prot_form - OrganL 0.9b pre-release

This example will form a bud by "blowing up" starting from a spherical cap.


## HOW TO RUN:

1. Copy the appropriate **organl** binary for your platform from `/binaries/` into this folder.
2. Run on commandline: `[organl] init.obj EAP LIP`


## Brief Description:

- This run uses Elastic Area Pressure - `EAP`, requiring a total of 1 energy functional parameters (set in *control.txt* under lines headed with 'E') together with active lipid mixing - `LIP`.
- The current *control.txt* blows gently and allows for a stable final structure. Use *control_fast.txt* as *control.txt* if this is taking too long.
- The starting point of the stimulation is a spherical cap (*init.obj*), with inhomogenous distribution of lipids - PC, OH and ELP (*props.txt*) on the faces. The properties of the lipids are specified in *lipids.txt*. Notice that ELP is a protein mimick as APL is -1. The boundary of the structure is ignored / frozen throughout the stimulation using the *blocks.txt*.
- The result is a bud formation. The variations of this is result is the starting point (*init.obj*) of examples `bud_enhance_prot` and `bud_enhance`.

## Important Parameters:
> Since there is a high degree of non-linearity in the stimulation, a precise behavioural evaluation of parametric values is not possible. The following may be viewed as a rough baseline to work on. Read the Manual for a detailed take on running stimulations and optimizing for good parametric values.

- `E LAMBDA 0 XX` parameter is used to denote the pressure. The choice is made in combination with the `Ka` value from *props.txt*. In general, higher the value the faster the blowup is.  