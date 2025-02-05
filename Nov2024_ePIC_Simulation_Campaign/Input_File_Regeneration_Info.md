# Love Preet (navisaharan3@gmail.com)- 19/11/2024 - University of Regina

To reproduce the files used for the November simulation campaign from DEMPgen, follow the steps below:

- Download the version 1.2.0 release of the generator from https://github.com/JeffersonLab/DEMPgen.

- After successfully compliling the genertor by following the instructions in the README.md file, the files can be generted using a shell script, namely Process_EIC.sh

- The Process_EIC.sh script was run with the following arguments -
  - ./Process_EIC.csh 1 14500000 5 41 1975090 HEPMC3 ip6 pi+
  - ./Process_EIC.csh 1 23000000 5 41 1975090 HEPMC3 ip6 pi+
  - ./Process_EIC.csh 1 95000000 5 41 1975090 HEPMC3 ip6 pi+
  - ./Process_EIC.csh 1 33000000 10 100 7187711 HEPMC3 ip6 pi+
  - ./Process_EIC.csh 1 30000000 10 100 7187711 HEPMC3 ip6 pi+
  - ./Process_EIC.csh 1 128000000 10 100 7187711 HEPMC3 ip6 pi+
  - ./Process_EIC.csh 1 53000000 18 275 6978096 HEPMC3 ip6 pi+
  - ./Process_EIC.csh 1 56000000 18 275 6978096 HEPMC3 ip6 pi+
  - ./Process_EIC.csh 1 200000000 18 275 6978096 HEPMC3 ip6 pi+
  - ./Process_EIC.csh 1 10300000 5 41 8681159 HEPMC3 ip6 K+ Lambda
  - ./Process_EIC.csh 1 18000000 5 41 8681159 HEPMC3 ip6 K+ Lambda
  - ./Process_EIC.csh 1 72800000 5 41 8681159 HEPMC3 ip6 K+ Lambda

- Note that DEMPgen was run three times for each energy combination by varying the scattered electron angles, defined as "e_Theta_Low" and "e_Theta_High" in the Config_EIC.json file, and with a different number of events, specified as a second argument, to ensure sufficient statistics in each Q2 bin. The required statistics are calculated based on the reconstructed events and efficiency values obtained during the physics analysis. Once the files are geneated, their naming scheme is modified based on the Q2 bin. The values of "e_Theta_Low" and "e_Theta_High" are recorded in a .txt file, along with information on the event generation such as the number of events cut at various stages. This is produced by the generator along with a corresponding .dat file which contains the generated events.
  - The .txt files generated from the execution of the command above are recorded in two sub directories within this folder -
    - pion - pi+ commands, subdivided into each beam energy combo
    - K+Lambda - K+ Lambda commands, subdivided into only one beam energy combo 
   
-  A specific set of events could be recreated by running Process_EIC.sh (or Process_EIC_iFarm.csh) with the corresponding random seed and the scattered electron ranges provided in the .txt file. For example, to recreate first file of the 5 on 41 pi+ events, we could run 
  - ./Process_EIC.csh 1 14500000 5 41 1975090 HEPMC3 ip6 pi+
    - The first argument is the file number and the 5th is the random seed
  - This should reproduce the input .dat file exactly

- The commands above also produce the output in an optional root output format. These files are produced and backed up on multiple systems, but are not included on GitHub.
  - Contact Stephen Kay (stephen.kay@york.ac.uk) or Love Preet (navisaharan3@gmail.com) for access to these root files

- Once the files are generted with DEMPgen, they are processed through the Monte Carlo Afterburner (https://github.com/eic/afterburner) to incorporate crossing angle, beam  effects, and veterx spread for the EIC.

- This is done by running the following command within the ePIC container/eic-shell. 
  - abconv /work/eic/users/preet/HEPMC3_Files/Nov2024_Files/pion/10on100/eic_DEMPgen_10on100_ip6_pi+_q2_3_10.dat -o eic_DEMPgen_10on100_ip6_pi+_q2_3_10
    - The first argument is the input file path. The corresponding files are provided in the Afterburner_Filelists directory for each beam energy combination within this folder.
    - The last argument is the output file name, based on the given input file.
