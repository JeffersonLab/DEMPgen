# DEMPgen

Event generator for Deep Exclusive Meson Production.

## Building

- To build the event generator, the compiler needs acces to both a compiled installation of CERN ROOT and its source code. ROOT version 6.08.06 or later is supported, and must be installed with the MathMore package enabled. Be sure not to change the location of either the ROOT source or compiled files after installation, as this will interfere with ROOT's built in CMake configurators.

- CMake is also required. CMake 2.8 is the minimum supported version, and CMake 3 has been tested as well.

- Before attempting to build, run the setup script -
  - source setup.sh
  - OR
  - source setup.csh

Depending upon shell.

- Next, create a build directory and cd to it. Take note of the location of the source directory (where CMakeLists.txt should be stored) and run the commands:

  - mdkir build
  - cd build
  - cmake ..
  - make -j8

- As a one liner - 
  - mkdir build && cd build && cmake ../ && make -j8

- The event generator can now be run using the following command from the main (root) directory.

  - ./build/DEMPgen Config.json
  - Data will be saved under data/OutputFiles/
  - As long as the path to the generator and config file are correct, DEMPgen should be executable from anywhere
  
### Building on the JLab iFarm

- Building on the JLab iFarm requires you to set up some software versions beforehand, to build successfully, I did the following - 
  - Comment out any CUE or other initialisation in your .login/.cshrc scripts
  - Login to an ifarm node
  - module use /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/modulefiles
  - module load root/6.30.06-gcc11.4.0
  - module load cmake/3.19.4
 
- Following this, build DEMGen using the one liner above without any issues, you will need to load these modules when running the generator subsequently (this is done by default in the farm job scripts)

## File/Directory Structure

The files within the top level directory and the directory structure of the generator are outlined below. Note that files within subdirectories are not outlined in full. Please refer to comments within each source/header files and README files within subdirectories (where applicable) for details on individual files within subdirectories.

### Top Level Directory Files

- [Batch_Submission_EIC.sh](Batch_Submission_EIC.sh) - Shell script to generate and submit EIC event generation jobs to a local batch queueing system. See file comments and description below for details on usage.
- [JLab_Batch_Submission.sh](JLab_Batch_Submission.sh) - Shell script to generate and submit EIC event generation jobs to the JLab farm batch queueing system. See file comments and description below for details on usage.
- [Process_EIC.csh](Process_EIC.csh) - .csh script that creates a `.json` config file based upon inputs and runs DEMPgen with produced configuration file. See file comments and description below for details on usage.
- [Process_EIC_iFarm.csh](Process_EIC_iFarm.csh) - .csh script that creates a `.json` config file based upon inputs and runs DEMPgen with produced configuration file. This version is configured for use on the JLab iFarm system. See file comments and description below for details on usage.
- [Config_EIC.json](Config_EIC.json) - An example `.json` config file for EIC event generation. Used by shell scripts below. See comments within file for details on flags.
- [Config_SoLID.json](Config_SoLID.json) - An example `.json` config file for SoLID event generation. See comments within file for details on flags.
- [CMakeLists.txt](CMakeLists.txt)- CMakeLists.txt file.
- [LICENSE](LICENSE) - Copyright and license information.
- [EvGenFlowChart.xml](EvGenFlowChart.xml) - Event generation flow chart.
- [Test_EIC.json](Test_EIC.json) - `.json` file with all configuration options to gereate a test EIC event sample. See `reference_README.md` within the `reference_output` directory for details.
- [README.md](README.md) - The README file.

### Directory Structure

- [data/](data/) - Stores output files from DEMPgen.
  - [data/input/](data/input) - Contains input cross-section files.
- [debug/](debug) - Contains debug files for the SoLID module.
- [include/](include) - Contains header files for the SoLID module and the following subdirectory:
  - [json/](include/json/) - Contains `.json` files associated with the header files for the SoLID module.
- [src/](src/) - Contains source files for DEMPgen and the following subdirectories:
  - [eic_evgen/](src/eic_evgen/) - Contains source files, header files, and pion cross-section parameterization files for the EIC module. Also contains the following subdirectories:
    - [CrossSection_Params/](src/eic_evgen/CrossSection_Params/) - Contains kaon cross-section parameterization files for the EIC module.
    - [process_routine/](src/eic_evgen/process_routine/) - Contains the main DEMPgen processing routine.
- [Dec2023_ePIC_Simulation_Campaign/](Dec2023_ePIC_Simulation_Campaign/) - Contain information about the files submitted to the ePIC simulation campaign in December 2023 respectively, see documentation within directory for more details.
- [Jul2024_ePIC_Simulation_Campaign/](Jul2024_ePIC_Simulation_Campaign/) - Contain information about the files submitted to the ePIC simulation campaign in July 2024, see documentation within directory for more details.
- [reference_output/](reference_output/) - Contains test sample files and the instructions on how to reproduce them, as detailed in the `reference_README.md` file within this directory.

### Instructions for Comprehensive Test Run

Instructions for a comprehensive test run, and a reference file to compare to (using for example, diff) are included in -

- [reference_output/reference_README.md](reference_output/reference_README.md)

## Configuration

The file Config_EIC.json and Config_SoLID.json contain all the configuration options for EIC and SoLID event generation respectively. Use these as a template for other configuration files, which may be given as an argument to the event generator.
- Note that Config_EIC.json is used in shell scripts (see below) that run large numbers of jobs on batch queueing systems. *Do not* modify this file unless you know what you're doing with it or you really need to edit something.
  - Copying the template to a new file and editing that is strongly recommended in all cases.


### Ejectile calculation methods

- The EIC module of DEMPgen has two different calculation methods that may be used to calculate the ejectile properties.
- These are referred to as the "Analytical" and "Solve" methods. Either can be used in a simulation run, just change the config file to switch between them.
- In your config .json file, select between the two versions using the "calc_method" argument, this can be set to Analytical or Solve, e.g.
  - "calc_method": "Analytical"
  - OR
  - "calc_method": "Solve"
  - If the method is not entered exactly as above, it will not be recognised and the generator will default to the analytical method.

- If you are using the shell scripts provided to run a lot of simulations, modify Config_EIC.json BEFORE running the shell scripts.

## Output

The event data is output to the configured file location relative to the build directory. The TTree in this file contains all kinematic data for all particles in the laboratory rest frame. Variables with the prefix "Vert" represent values read at the interaction vertex. Values with the prefix "Lab" represent values read after all correcting effects (multiple scattering, ionization, etc.) have been applied. 

## Sources and Classes

### Particle Class

The Particle class inherits from the root TLorentzVector class so that its functions are immediately accessible through the Particle objects. It's member variables include its mass, charge, and its GEMC compatible id.

### DEMPEvent Class

The Event class contains all particles of a single event, in one frame of reference. It includes methods enabling transformation to other frames of reference and coordinate systems.

### Asymmetry Class

The Asymmetry class calculates asymmetry amplitudes based on Monte-Carlo data by Goloskokov and Kroll. The asymmetry objects are persistent between events.

### Customrand Class

A class to hold various distribution functions for randomly generated variables, such as scattered electron energy. These functions should be persistent between events, and will provide fast random numbers.

### ScatteredParticleGen

Stores the kinematic ranges for the scattered electron and generates them with random energy and direction within the range, using sphere point picking.

### TargetGen

Generates the target neutron with Fermi momentum (when enabled) and proton for FSI (when enabled, may also have Fermi momentum).

### ProductGen

This class reads in the kinematic variables for the incident, target, and scattered particles and uses conservation laws to solve for the remaining two particles. The pion is first given a random direction. A fast root finding algorithm then calculates the pion's momentum magnitude, which is then used to find the proton's momentum. Pion direction may also be passed as an argument for debugging purposes.

### SigmaCalc

This class returns the cross sections and weights for the event. It acts as an interface between the current version of the event generator and header files from the old event generator (seen under branch "original"), which contain the parameterization of the cross sections.

### TreeBuilder

Manages the ROOT TTree object to be stored in the output .root file. Contains methods to easily add all kinematics data stored in a particle or DEMPEvent object. 

### Matter Effects

Transforms particle kinematics based on three effects caused by transition through matter: Ionization, Bremsstrahlung, and Multiple Scattering. 

### FSI

Computes the effect of final state interaction between the produced pion and one of the two protons of the target nucleus. The momentum of the outgoing pion and recoiled proton, as well as  the cross section of the interaction is calculated based on elastic pion-nucleon scattering.

### JsonCpp

This project uses [JsonCpp](https://github.com/open-source-parsers/jsoncpp "JsonCpp Github") to read in configuration options. The amalgamated sources for JsonCpp are redistributed with this project in compliance with the MIT license.

## Processing Scripts and json examples

- There are several shell scripts in the main directory that can be used to run the generator and produce some output, information on these scripts is provided below.
- To remove clutter from the main directory, example .json scripts were moved to a new folder -
  - json_examples
- Remember to set the ejectile calculation method in the Config_EIC.json file *before* running these scripts if you do not want to use the default analytical method.

### Process_EIC.csh

!!! NOTICE !!!  
This script copies Config_EIC.json and formats a new file based upon this, DO NOT MODIFY Config_EIC.json (other than the calculation method) if you want to use this script!  
!!! NOTICE !!!  

- To facilitate the submission of batch jobs, I created a csh script to automatically construct .json config files and run them. This script can also be utilised to run the generator manually, without the need to go and edit a json file. 
- The script requires 8 arguments (which is a lot, I know), but in the K+ case, it expects 9. They are as follows -  

  - Arg 1 - FileNum -> For batch running, we typically run X files of Y events, this argument is just X, if you're running manually as a test, just input 1 or whatever you fancy.
  - Arg 2 - NumEvents -> The number of events thrown for this file, set this to whatever you want to run. For reference, with the Pi+/K+ generator, 1B files takes ~1 hour.
  - Arg 3 - EBeamE -> The electron beam energy, set this to whatever you want, typically, we use 5, 10 or 18 (the nominal max for the EIC). 
  - Arg 4 - HBeamE -> The hadron beam energy, again, set this to whatevr you want. Typically we use 41, 100 or 275 (41 and 275 being the nominal min/max).
  - Arg 5 - RandomSeed -> The random seed, self explanatory. Set this however you like, the batch submission job randomly generates a random seed to feed in here.
  - Arg 6 - OutputType -> The format of the output file, select from LUND, Pythia6 (for ECCE/Fun4All) or HEPMC3 (for ePIC), the default is HEPMC3 if your choice is invalid.
  - Arg 7 - InteractionPoint -> The interaction point, choose from ip6 or ip8. The default is ip6 if your choice is invalid.
  - Arg 8 - Ejectile -> The produced ejectile (meson) in the reaction, choose from omega, pi+, pi0 or K+.
  - Arg 9 - RecoilHadron -> OPTIONAL - This only matters if you select K+ as the ejectile, in this case, choose from Lambda or Sigma0 here. If your choice is invalid (or you don't specify arg9), the default is Lambda.
  
- So as an example if you executed the following -  
  
  - ./Process_EIC.csh 1 100000 18 275 24432 HEPMC3 ip6 K+ Lambda
  
- You would run the generator for 18 GeV e- on 275 protons for ip6, throwing 100000 events with the K+/Lambda generator.
  
### Batch_Submission_EIC.sh

- This script creates and submits batch jobs. It is designed for use with the torque queueing system on Lark at the University of Regina.
The jobs the script creates and submits all execute the Process_EIC.csh script described above. This script requries a very similar set of arguments -  

  - Arg 1 - NumFiles -> The batch script is designed to run X jobs of Y events, this number is just X, the number of files you want to run.
  - Arg 2 - NumEvents -> The number of events thrown for this file, set this to whatever you want to run. For reference, with the Pi+/K+ generator, 1B files takes ~1 hour.
  - Arg 3 - EBeamE -> The electron beam energy, set this to whatever you want, typically, we use 5, 10 or 18 (the nominal max for the EIC).
  - Arg 4 - HBeamE -> The hadron beam energy, again, set this to whatevr you want. Typically we use 41, 100 or 275 (41 and 275 being the nominal min/max).
  - Arg 5 - OutputType -> The format of the output file, select from LUND, Pythia6 (for ECCE/Fun4All) or HEPMC3 (for ePIC), the default is HEPMC3 if your choice is invalid. 
  - Arg 6 - InteractionPoint -> The interaction point, choose from ip6 or ip8. The default is ip6 if your choice is invalif.
  - Arg 7 - Ejectile -> The produced ejectile (meson) in the reaction, choose from omega, pi+, pi0 or K+.
  - Arg 8 - RecoilHadron -> OPTIONAL - This only matters if you select K+ as the ejectile, in this case, choose from Lambda or Sigma0 here. If your choice is invalid (or you don't specify arg9), the default is Lambda.

- The script automatically generates a random seed itself using the /dev/urandom function.

### Process_EIC_iFarm.csh

- This version is for use on the JLab iFarm/Farm.
- This script uses the same arguments as Process_EIC.csh.
- If you are processing interactively (i.e. on the iFarm), you will need to make sure that you have executed -
  - module use /group/halla/modulefiles
  - module load root
- You should also check the path set at the top looks OK.
    - By default, /eic/users/${USER}/DEMPgen is assumed.

### JLab_Batch_Submission.sh

- This version should be used to submit jobs to the Farm batch queueing system (swif2).
- It uses the same arguments as Batch_Submission_EIC.sh.
- You should also check the paths set throughout the script look ok (for example, in the COMMAND: ... line).
    - By default, /eic/users/${USER}/DEMPgen is assumed.

## Acknowledgments

## License
DEMPgen is licensed under the GNU General Public License v3.0.