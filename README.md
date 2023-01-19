# PROVER-M
PROVER-M is a near-field model for projecting disposals of fine sediments in coastal and estuarine environments.


### Introduction

PROVER-M ("Prozessverständnis von Verklappungsvorgängen von Feinsediment in tidebeeinflussten Gewässern") is a near-field model for projecting disposals of fine sediments in coastal and estuarine environments.


### Installation

Two options are provided for running PROVER-M:

**Option A**: This option is designed for users that want to apply PROVER-M without any adaptions to the source code. A stand-alone *PROVER.exe* file for Windows is provided that installs PROVER-M, including all necessary files for running the model. By double-clicking on *PROVER.exe* and following the installation instructions, the model is installed. After installation, users can call PROVER-M by executing the installed program from their system. Input for simulating disposals may be inserted in the graphical user interface (GUI).

**Option B**: This option is designed for users that want to use the PROVER-M source code and possibly adapt or expand parts thereof. Here, no further installation is needed, but having MATLAB installed is a prerequisite. All scripts necessary to run PROVER-M are located in the *scr* directory. If the GUI shall be used to provide input, the user needs to execute the *PROVER_M.mlapp* file to open the GUI in MATLAB. All other scripts can simply be opened and edited in the MATLAB editor. A list and short description of the scripts can be found below.


### Use of program

**How to use**
* To begin a simulation, the GUI App needs to be started either by: 
	1. Executing the installed PROVER-M application.
	2. Opening the scr-directory and double-clicking on the file *PROVER_M.mlapp* (requires MATLAB).
	3. Opening the *PROVER_M.mlapp* in MATLAB (requires MATLAB). 
* The model input can be set by changing individual parameters or by selecting an existing input text file via the "Input case" dropdown menu or the "Load.."-button.
* After setting the desired input parameters, the configuration can be saved as a text file by clicking the "Save as.."-button. 
* The input can be reset to the last saved version of the selected input case in two ways. While "Reset sediment data" resets the sediment characteristics, "Reset input" resets the entire input configuration.
* By clicking the "START"-button, the simulation is initiated.
* In the right half of the GUI, a live feed of the cloud propagation through the water column and main parameters is presented.


### Input

*(including suggested order of magnitude)*

---

**Ambient parameters** *(incl. suggested order of magnitude)*
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
* -1 to 1 (For positive values [0-1], the reduction factor determines the fraction of the sediment to be settled. For a reduction factor of -1, the settling is only based on the critical shear stress.)
---
**Coefficients** *(All coefficients should be selected carefully and within the limits of the given order of magnitude)*
- Entrainment phase 1 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  *(10<sup>-1</sup>)*
- Entrainment phase 2 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  *(10<sup>-1</sup>)*
- Mass  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  *(10<sup>0</sup>)*
- Drag phase 1  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  *(10<sup>0</sup>)*
- Drag phase 2 	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  *(10<sup>-2</sup> to 10<sup>-1</sup>)*
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
- Allowing the material to be stripped (applies to fine sediments, including fine sand)

--- 


### Output
---
**Output settings**
- Output time step interval (sets the interval for the output variables for stripped and settled sediments)
- Time step for plotting 
	* 'Yes': Uses the output time step interval for updating the live feed of the cloud propagation through the water column within the GUI.
	* 'No': Updates the live feed of the cloud propagation through the water column within the GUI using a preset interval of 25 time steps.
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
- Cloud and sediment variables as MAT-files
---
*Output files are stored in a directory "output". If the PROVER-M stand alone .exe has been used, the data will be stored under C:\Users\username\AppData\Local\Temp\username\mcrCache9.12\PROVER0\PROVER_M\output. The accessability will be changed with the next update.*


### Functionality (Flow chart)
- An overview of the functionality of the program code in the form of a flow chart can be found \[here](link to paper).


### Prerequisites
For executing the *PROVER.exe* on a Windows system no necessary packages or programs are needed.
 
 
### Source Code 
If a user wants to adapt or change the source code, the scr-directory includes the following files:
- *PROVER_M.mlapp*
	- A MATLAB App that starts the GUI, which may be used for selecting input parameters, starting simulations and accessing the live feed.
- *prover_m_main.m*
	- Main program code, where the bookkeeping of the cloud and ambient parameters occurs.
- *prover_m_rk4.m*
	- This numerical solver function utilizes the Runge-Kutta 4th order method to approximate parameter gradients
- *prover_m_phase1.m*
	- In this function, the addressed variables in the phase of convective descent are calculated using the conservation equations.
- *prover_m_phase2.m*
	- In this function, the addressed variables in the phase of dynamic collapse are calculated using the energy concept equations. 

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

  			
