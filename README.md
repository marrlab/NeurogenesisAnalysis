# NeurogenesisAnalysis
contains code and data accompanying
### Subtle Changes in Clonal Dynamics Underlie the Age-related Decline in Neurogenesis
by Lisa Bast<sup>1,2,\*</sup>, Filippo Calzolari<sup>3,4,5,\*,†</sup>, Michael Strasser<sup>1,6</sup>, Jan Hasenauer<sup>1,2</sup>, Fabian Theis<sup>1,2</sup>, Jovica Ninkovic<sup>3,4,†</sup> and Carsten Marr<sup>1,†</sup>

<sup>\*</sup> Equal contribution
<sup>†</sup> Correspondence (shared senior authorship)
 
- download the required data from https://hmgubox.helmholtz-muenchen.de/d/dedd78d4e8c6463abfd0/
- required software: MATLAB (R2017a)
- usage of Toolboxes:
  - CERENA (https://github.com/CERENADevelopers/CERENA/)
  - PESTO (https://github.com/ICB-DCM/PESTO/)
  - AMICI (https://github.com/ICB-DCM/AMICI)
  which are already included in folder 'tools'

###  Pre-analysis: 
RUN_N_activationSystem.m provides an ODE stem cell compartment fit using toolboxes AMICI and PESTO. 


###  Neurogenesis model 
#### 1. Modelselection
RUN_N_modelselection.m defines and fits all 64 models to both data sets (young and old mice) using toolboxes CERENA 	and PESTO, and calculates BIC scores. RUN_N_getModelselectionResults.m provides and visualises model selection results and applies model averaging
#### 2. Evaluation
RUN_N_evaluation.m builds age-independent and age-dependent models (from Hill-function fits) and compares it to data of 	Shook et al. and Daynac et al.
#### 3. Simulation
RUN_N_simulation.m provides model prediction from SSA simulations (using CERENA toolbox) and tree simulations.



