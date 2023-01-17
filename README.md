# PROVER-M
PROVER-M is a near-field model for projecting disposals of fine sediments in coastal and estuarine environments.


### Introduction

PROVER-M ("Prozessverständnis von Verklappungsvorgängen von Feinsediment in tidebeeinflussten Gewässern") is a near-field model for projecting disposals of fine sediments in coastal and estuarine environments.


### Installation

Two options are provided for running PROVER-M:

**Option A**: This option is designed for users that want to apply PROVER-M without any adaptions to the source code. A stand-alone _PROVER.exe_ file for Windows is provided that installs PROVER-M, including all necessary files for running PROVER-M, is provided. By double-clicking on _PROVER.exe_ and following the installation instructions, the model is installed. After installation, the user can call PROVER-M by executing the installed program from your system. Now the GUI opens and the user can insert his simulation input for simulating disposals.

**Option B**: Is designed for users that want to use the PROVER-M source code and possibly adapt or expand parts of it. Here, no further installation is needed, but having MATLAB installed is a prerequisite. In the subdirectory _scr_ all scripts necessary to run PROVER-M are provided. If input should be given through the GUI, the user needs to execute _main_appDesigner.mlapp_ and the GUI should open within MATLAB. All other scripts can simply be opened and edited within the MATLAB editor. A list and short description of the scripts can be found below.


### Use of program

**How to use**
* To begin a simulation, first, the GUI App needs to be started 
	* either by executing the installed PROVER-M application.
	* or by opening the scr-directory and double klick "PROVER_M.mlapp" (requires MATLAB).
	* or opening and playing the "PROVER_M.mlapp" within MATLAB (requires MATLAB). 
* The model input can be set by changing individual parameters or by selecting an existing input text file via the "Input case" dropdown menu or the "Load.."-Button.
* After setting the desired input parameters, the configuration can be saved as a text file by clicking the "Save as.."-Button. 
* An exemplary simulation case based on the study of Delo & Burt (1987) is configured through "Default". Please check the publication for a more detailed case description.
* By clicking the "START"-Button, the simulation is initiated.
* In the right half of the GUI a live feed of the cloud propagation through the water column and main parameters is presented.


### Input

*(including order of magnitude)*

---

**Ambient parameters**
- Water depth   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*(10<sup>1</sup> to 10<sup>2</sup>)*
- Ambient velocities   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*(-10<sup>0</sup> to 10<sup>0</sup>)*
- Ambient densities   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*(10<sup>3</sup>)*
---	
**Hopper settings**
- Disposal volume   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*(10<sup>2</sup> to 10<sup>4</sup>)*
- Hopper draft   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*(10<sup>0</sup> to 10<sup>1</sup>)*
- Dumping instances &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; *(divides a disposal into n-equal intervals that are disposed consecutively)*
---	
**Settling**			
* -1 to 1 ; pure reduction factor for positive values. For negative values the reduction factor is based on the critical shear stress.
---
**Coefficients** *(All coefficients should be changed carefully and within the order of magnitude)*
- Entrainment Phase 1 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  *(10<sup>-1</sup>)*
- Entrainment Phase 2 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  *(10<sup>-1</sup>)*
- Mass  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  *(10<sup>0</sup>)*
- Drag Phase 1  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  *(10<sup>0</sup>)*
- Drag Phase 2 	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  *(10<sup>-2</sup> to 10<sup>-1</sup>)*
- Friction  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  *(10<sup>-2</sup>)*
- Stripping  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  *(10<sup>-3</sup>)*
---
**Sediment characteristics**
- Type
- Density
- Volumetric fraction
- Fall velocity
- Void ratio
- Critical shear stress
- Cohesiveness
- Allow the material to be stripped (true for fine sediments, including fine sand)

--- 


### Output
---
**GUI live feed**
- Time
- Cloud radius
- Stripping volume
- Cloud width
- Cloud height
- Settling volume
---
**Output file**
- Cloud parameters as textfile
- Stripped and settled sediments as textfile
- Cloud and sediment variables as .mat files
---


### Functionality (Flow chart)
- An overview of the functionality of the program code in the form of a flow chart can be found \[here](link to paper).


### Necessary packages
 For executing the PROVER.exe on a Windows system no necessary packages or programs are needed.
 
 Additionally, the scr-directory includes the following files, if you want to adapt or change the source code:
- *PROVER_M.mlapp*
	- A MATLAB App that starts the graphical user interface (GUI), for input parameters, simulation start and live feed.
- *prover_m_main.m*
	- Main program code, where the bookkeeping of the cloud and ambient parameters occurs.
- *prover_m_rk4.m*
	- Parameter gradients are approximated utilizing the Runge-Kutta 4th order method in this numerical solver function.
- *prover_m_phase1.m*
	- The addressed variables in the phase of convective descent are calculated using the conservation equations.
- *prover_m_phase2.m*
	- The addressed variables in the phase of dynamic collapse are calculated using the energy concept equations. 

For running the source code files, an installed and licened version of MATLAB by Mathworks is needed.


### License

* GNU GPL License

PROVER-M has been developed in MATLAB and is provided via the Application Compiler provided by MATLAB®. © 1984 - 2022a The MathWorks, Inc.


### More information

A full account of software functionalities and implemented methods can be found \[here](link to paper).


### References

* Delo, E. and Burt, T.N. (1987) Dispersal of dredged material - Tees field study September 1986. Technical Report. Hydraulics Research Wallingford. 
* Koh, R. C., & Chang, Y. C. (1973). Mathematical model for barged ocean disposal of wastes (Vol. 1). Office of Research and Development, US Environmental Protection Agency.
* Brandsma, M. G., & Divoky, D. J. (1976). Development of models for prediction of short-term fate of dredged material discharged in the estuarine environment. TETRA TECH INC PASADENA CALIF.
* Johnson, B. H.; Fong, M. T. (1995): Development and verification of numerical models for predicting the initial fate of dredged material disposed in open water. Report 2 - Theoretical developments and verification results (Technical Report DRP-93-1).
* Escudier, M. P., & Maxworthy, T. (1973). On the motion of turbulent thermals. Journal of Fluid Mechanics, 61(3), 541-552.
* Turner, J. S. (1960). A comparison between buoyant vortex rings and vortex pairs. Journal of Fluid Mechanics, 7(3), 419-432.
* Partheniades, E. (1965). Erosion and deposition of cohesive soils. Journal of the Hydraulics Division, 91(1), 105-139.

  			
