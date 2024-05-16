# example/ade - OrganL 0.9b pre-release

This example on Area Difference Elasticity (ADE) and features the use of contraints on area, volume, and the integrated mean curvature for Monte Carlo Stimulation of the Membrane. Result is comparable to Fig 7. (b) of the manuscript.

Refer to Supplementary Material to the manuscript for details on implementation of ADE.


## HOW TO RUN:

1. Copy the appropriate **organl** binary for your platform from `/binaries/` into this folder.
2. Run on commandline: `[organl] init.obj ADE`


## Brief Description:

- This run features the Area Difference Elasticity - `ADE` where the area difference between leaflets and the bilayer thickness is accounted for in addition to DCH functional. This requires a total of 6 energy functional parameters as set in _control.txt_ under lines headed with `E`. 
- The starting point of the stimulation is a almost prolate mesh (_init.obj_), with uniform kappa (_props.txt_) on all faces. 
- A typical result is a heavy neck formation in the middle as shown in Fig. 7 (b)-right.

## Important Parameters:
> Since there is a high degree of non-linearity in the stimulation, a precise behavioural evaluation of parametric values is not possible. The following may be viewed as a rough baseline to work on. Read the Manual for a detailed take on running stimulations and finding good values for the parameters.

- `E DEF 0 XX` - sets the value of $A_0$. Compare DCH and CHF energies to check if $A_0$ values are reasonable. Considerable smaller values results in heavy distortions.
- `E DEF 1 XX` - sets the value of $V_0$. Considerable smaller values is seen to result in abnormal spike formations and deformations.
- `E DEF 2 XX` - sets the value of $M_0$. High values result in a smoother initial structure.
- `E LAMBDA 0 XX` used to denote the constraint on area (Surface Tension). Lower $\lambda_0$ value is seen to commonly result in more branched budding.
- `E LAMBDA 1 XX` used to denote the constraint on volume (pressure). Lower $\lambda_1$ seems to result in heavy elongation along length of the prolate structure.
- `E LAMBDA 2 XX` used to denote the constraint on integrated mean curvature. Low values of $\lambda_2$ is seen to give only slight changes from the initial structure. Could possibly be used to smoothen out the mesh.