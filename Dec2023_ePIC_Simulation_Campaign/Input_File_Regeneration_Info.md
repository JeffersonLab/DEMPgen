# Love Preet - 14/12/2023 - University of Regina

To reproduce the files used for the December simulation campaign from DEMPgen, follow the provided steps:

- Download the generator from https://github.com/JeffersonLab/DEMPgen/tree/develop.

- After successfully compliling the genertor by following the instructions in the README.md file, the files can be generted using a shell script, namely Process_EIC.sh
  - To produce large numbers of files, it is easier to use a batch computing system. A script to create and generate files is provided - Batch_Submission_EIC.sh
  - This script was used to produce the files for the December 2023 simulation campaign. An equivalent script to produce jobs on the JLab farm - JLab_Batch_Submission.sh can also be used.

- The Batch_Submission_EIC.sh script was run with the following arguments -
  - ./Batch_Submission_EIC.sh 60 100000000 5 41 HEPMC3 ip6 pi+
  - ./Batch_Submission_EIC.sh 60 1000000000 5 100 HEPMC3 ip6 pi+
  - ./Batch_Submission_EIC.sh 150 1000000000 10 100 HEPMC3 ip6 pi+ 
  - ./Batch_Submission_EIC.sh 60 100000000 5 41 HEPMC3 ip6 K+ Lambda
  - ./Batch_Submission_EIC.sh 60 1000000000 5 100 HEPMC3 ip6 K+ Lambda
  - ./Batch_Submission_EIC.sh 150 1000000000 10 100 HEPMC3 ip6 K+ Lambda
  - ./Batch_Submission_EIC.sh 60 100000000 5 41 HEPMC3 ip6 K+ Sigma0
  - ./Batch_Submission_EIC.sh 60 1000000000 5 100 HEPMC3 ip6 K+ Sigma0
  - ./Batch_Submission_EIC.sh 150 1000000000 10 100 HEPMC3 ip6 K+ Sigma0 

- The batch submission scripts submits these arguments to Process_EIC.sh (or Process_EIC_iFarm.csh in the JLab farm), see these scripts for further details on the arguments provided.

- When the batch job is created, a random seed is generated using the - od -An -N3 -i /dev/urandom - command, such that each job has a "random" random seed.

- The random seed for each job is recorded in a .txt file, along with information on the event generation such as the number of events cut at various stages. This is produced by the generator along with a corresponding .dat file which contains the generated events.
  - The .txt files generated from the execution of the command above are recorded in three sub directories within this folder -
    - pion - pi+ commands, subdivided into each beam energy combo
    - K+Lambda - K+ Lambda commands, subdivided into each beam energy combo
    - K+Sigma0 - K+ Sigma0 commands, subdivided into each beam energy combo 

-  A specific set of events could be recreated by running Process_EIC.sh (or Process_EIC_iFarm.csh) with the corresponding random seed. For example, to recreate file 41 of the 5 on 100 pi+ events, we could run 
  - ./Process_EIC.csh 41 1000000000 5 100 12034245 HEPMC3 ip6 pi+
    - The first argument is the file number and the 5th is the random seed
  - This should reproduce the input .dat file exactly

- Note that for the 5on41 energy combination for both pions and kaons (lambda and sigma), to achieve full distribution coverage, the ejectile angle was changed to 100. 
  - This is denoted as "EjectileX_Theta_High": 100 in the Config_EIC.json file.
  - Edit Config_EIC.json with this change BEFORE running the 5on41 events 
- For the other energy combinations(or 5on100 and 10on100), this parameter was left as the default - "EjectileX_Theta_High": 50.
  - Set back to this value if running after a 5on41 set of files

- Once the files are generted with DEMPgen, they are processed through the Monte Carlo Afterburner (https://github.com/eic/afterburner) to incorporate crossing angle, beam  effects, and veterx spread for the EIC.

- This is done using the - Process_Files_Afterburner.sh - script included in this directory.
** This script is run within the ePIC container/eic-shell. See the script for further details.
** The script takes file lists as arguments, these are included in this directory for reference/convenience.
