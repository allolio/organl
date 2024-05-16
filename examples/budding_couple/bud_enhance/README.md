# example/budding_couple/bud_enhance - OrganL 0.9b pre-release


These examples show the ability to manipulate, form and constrict membrane buds using mixing terms. The results are comparable to Fig. 7 of the manuscript.


## HOW TO RUN:

1. Copy the appropriate **organl** binary for your platform from `/binaries/` into this folder.
2. Run on commandline: `[organl] init.obj CHF LIP`


## Brief Explaination:

- Run uses Constrained Helfrich Functional - `CHF`, requiring a total of 4 energy functional parameters (set in *control.txt* under lines headed with 'E') together with active lipid mixing - `LIP`.
- The starting point of the stimulation is an inflated spherical cap (*init.obj*), with inhomogenous distribution of lipids - PC and OH (*props.txt*) on the faces. The properties of the lipids are specified in *lipids.txt*. The boundary of the structure is ignored / frozen throughout the stimulation using the *blocks.txt*.
- A typical result is a strong constraint formation as shown in Fig. 7 (c)-left.

## Important Parameters:
> Since there is a high degree of non-linearity in the stimulation, a precise behavioural evaluation of parametric values is not possible. The following may be viewed as a rough baseline to work on. Read the Manual for a detailed take on running stimulations and optimizing for good parametric values.

- `E DEF 0 XX` - sets the value of $A_0$. Setting it too low results in heavy distortions. Compare DCH and CHF energies to check if $A_0$ values are reasonable. Lower values is seen to produce wider neck constraint and a symmetric blob formation.
- `E DEF 1 XX` - sets the value of $V_0$. Higher values of $V_0$ is seen to result in more consistent blob heads with wider neck.
- `E LAMBDA 0 XX` used to denote the constraint on area (Surface Tension). Lower $\lambda_0$ value is seen to result in elongation in height of the blob with a tendency to bud.
- `E LAMBDA 1 XX` used to denote the constraint on volume (pressure). Low $\lambda_1$ is seen to form heavy and up-right neck constraining.