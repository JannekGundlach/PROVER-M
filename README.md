# PROVER-M
PROVER-M is a near-field model for projecting disposals of fine sediments in coastal and estuarine environments.

### Introduction

PROVER-M (Prozessverständnis von Verklappungsvorgängen von Feinsediment in tidebeeinflussten Gewässern) is a near-field model for projecting disposals of fine sediments in coastal and estuarine environments.

### Installation

Two options are provided for running PROVER-M:

**Option A**: Is designed for users that want to apply PROVER-M without any adaptions to the source code. A stand-alone _PROVER.exe_ file, that installs PROVER-M on your system (Windows), including all necessary files for running PROVER-M, is provided. By double-clicking on _PROVER.exe_ and following the installation instructions, the model is installed. After installation, the user can call PROVER-M by executing the installed program from your system. Now the GUI opens and the user can insert his simulation input for simulating disposals.

**Option B**: Is designed for users that want to use the PROVER-M source code and possibly adapt or expand parts of it. Here, no further installation is needed, but having MATLAB installed is a prerequisite. In the subdirectory _scr_ all scripts necessary to run PROVER-M are provided. If input should be given through the GUI, the user needs to execute _main_appDesigner.mlapp_ and the GUI should open within MATLAB. All other scripts can simply be opened and edited within the MATLAB editor. A list and short description of the scripts can be found below.

### Use of program

**How to use**
* To begin a simulation, first, the GUI App needs to be started 
	* either by executing the installed PROVER-M application.
	* or by opening the scr-directory and double klick "PROVER_M.mlapp".
	* or opening and playing the "PROVER_M.mlapp" within MATLAB. 
* The model input can be set by changing individual parameters or by selecting an existing input text file via the "Input case" dropdown menu or the "Load.."-Button.
* After setting the desired input parameters, the configuration can be saved as a text file by clicking the "Save as.."-Button. 
* An exemplary simulation case based on the study of Delo & Burt (1987) is configured through "Default". Please check the publication for a more detailed case description.
* By clicking the "START"-Button, the simulation is initiated.
* In the right half of the GUI a live feed of the cloud propagation through the water column and main parameters is presented.

### Input

*(including order of magnitude)*

---

**Ambient parameters**
	- Water depth 			*(10<sup>1</sup> to 10<sup>2</sup>)*
	- Ambient velocities		*(-2.5 to 2.5)*
	- Ambient densities		*(around 10<sup>3</sup>)*
---	
**Hopper settings**
	- Disposal volume		*(10<sup>2</sup> to 10<sup>4</sup>)*
	- Hopper draft			*(10<sup>0</sup> to 10<sup>1</sup>)*
	- Dumping instances *(divides a disposal into n-equal intervals that are disposed consecutively)*
---	
**Settling**			*(-1 to 1 ; pure reduction factor for positive values. For negative values the reduction factor is based on the critical shear stress.)*
---
**Coefficients** *(All coefficients should be changed carefully and within an order of magnitude)*
	- Entrainment 	(Phase 1)	*(10<sup>-1</sup>)*
	- Entrainment 	(Phase 2)	*(10<sup>-1</sup>)*
	- Mass				*(10<sup>0</sup>)*
	- Drag 		(Phase 1)	*(10<sup>0</sup>)*
	- Drag 		(Phase 2)	*(10<sup>-2</sup> to 10<sup>-1</sup>)*
	- Friction			*(10<sup>-2</sup>)*
	- Stripping			*(10<sup>-3</sup>)*
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

###Output
**GUI live feed**
	- Time
	- Cloud radius
	- Stripping volume
	- Cloud width
	- Cloud height
	- Settling volume
**Output file**
	- Cloud parameters as textfile
	- Stripped and settled sediments as textfile
	- Cloud and sediment variables as .mat files



###Functionality (Flow chart)
	\- An overview of the functionality of the program code in the form of a flow chart can be found \[here](link to paper).

### Necessary packages

Included in the download are the following files:
	\- main_appDesigner.mlapp
		\- A MATLAB App that starts the graphical user interface (GUI), for input parameters, simulation start and live feed.
	\- prover_m_main.m
		\- Main program code, where the bookkeeping of the cloud and ambient parameters occurs.
	\- prover_m_rk4.m
		\- Parameter gradients are approximated utilizing the Runge-Kutta 4th order method in this numerical solver function.
	\- prover_m_phase1.m
		\- The addressed variables in the phase of convective descent are calculated using the conservation equations.
	\- prover_m_phase2.m
		\- The addressed variables in the phase of dynamic collapse are calculated using the energy concept equations. 

### License

* GNU GPL License

PROVER-M has been developed in MATLAB and is provided via the Application Compiler provided by MATLABÂ®. Â© 1984 - 2022a The MathWorks, Inc.

### More information

A full account of software functionalities and implemented methods can be found \[here](link to paper).

??? Please file any comments or concerns through the github issue tracker \[here](link to github).
	\- wollen wir das ???

### References

* Bowers, G. W., & Goldenblatt, M. K. (1978). Calibration of a predctive model for instantanneously discharged dedged material. US Environmental Protecion Agency, Corvallis, OR. EPA-699/3-78-089.
* Brandsma, M. G., & Divoky, D. J. (1976). Development of models for prediction of short-term fate of dredged material discharged in the estuarine environment. TETRA TECH INC PASADENA CALIF.
* ??? Cheng 1999 --> siehe "func_phase2_stfate.m", Zeile 122
  		\- Cheng, N. -., & Chiew, Y. -. (1999). Analysis of initiation of sediment suspension from bed load. Journal of Hydraulic Engineering, 125(8), 855-861. doi:10.1061/(ASCE)0733-9429(1999)125:8(855)			- ist das das richtige???
  Escudier, M. P., & Maxworthy, T. (1973). On the motion of turbulent thermals. Journal of Fluid Mechanics, 61(3), 541-552.
* Johnson, B. H.; Fong, M. T. (1995): Development and verification of numerical models for predicting the initial fate of dredged material disposed in open water. Report 2 - Theoretical developments and verification results (Technical Report DRP-93-1).
* ??? Turner, J. S. (1960). A comparison between buoyant vortex rings and vortex pairs. Journal of Fluid Mechanics, 7(3), 419-432.
  		\- weg? siehe "func_phase1_stfate.m" Zeile 121
  			
