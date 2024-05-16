# example/budding_couple/bud_enhance_prot - OrganL 0.9b pre-release


These examples show the ability to manipulate, form and constrict membrane buds using mixing terms. This contains the effect of protein mimicks and the results are comparable to Fig 7 in the manuscript.


## HOW TO RUN:

1. Copy the appropriate **organl** binary for your platform from `/binaries/` into this folder.
2. Run on commandline: `[organl] init.obj CHF LIP`


## Brief Description:

- This run uses Constrained Helfrich Functional - `CHF`, requiring a total of 4 energy functional parameters (set in *control.txt* under lines headed with 'E') together with active lipid mixing - `LIP`.
- The starting point of the stimulation is an inflated spherical cap (*init.obj*), with inhomogenous distribution of lipids - PC, OH and ELP (_props.txt_) on the faces. The properties of the lipids are specified in *lipids.txt*. Notice that ELP is a protein mimick as APL is -1. The boundary of the structure is ignored / frozen throughout the stimulation using the *blocks.txt*.
- A typical result is a strong neck formation as shown in Fig. 7 (c)-right.


## Important Parameters:
> Since there is a high degree of non-linearity in the stimulation, a precise behavioural evaluation of parametric values is not possible. The following may be viewed as a rough baseline to work on. Read the Manual for a detailed take on running stimulations and optimizing for good parametric values.

- `E DEF 0 XX` - sets the value of $A_0$. Lower values suppresses the blob formation considerably.  
- `E DEF 1 XX` - sets the value of $V_0$. Lower $V_0$ values is seen to bend the blob above neck.
- `E LAMBDA 0 XX` used to denote the constraint on area (Surface Tension). Lower $\lambda_0$ value is seen to result in a heavier blob head with a short neck.
- `E LAMBDA 1 XX` used to denote the constraint on volume (pressure). Low values of $\lambda_1$ is typically seen to result in elongated blob. 