# Stephen JD Kay (stephen.kay@york.ac.uk)- 18/07/2025 - University of York

To reproduce the files used for the July 2025 simulation campaign from DEMPgen, follow the steps below:

- Download the version 1.2.4 release of the generator from https://github.com/JeffersonLab/DEMPgen.

- After successfully compliling the genertor by following the instructions in the README.md file, remember to source setup.csh before running DEMPgen. Then, execute DEMPgen on the .json config files included in this directory. These were executed via -
  - ./build/DEMPgen Config_EIC_10on130_ip6_Pi+_q2_3_10_11000000.json
  - ./build/DEMPgen Config_EIC_10on130_ip6_Pi+_q2_10_20_170600000.json
  - ./build/DEMPgen Config_EIC_10on130_ip6_Pi+_q2_20_35_350000000.json
  - ./build/DEMPgen Config_EIC_10on250_ip6_pi+_q2_3_10_100600000.json
  - ./build/DEMPgen Config_EIC_10on250_ip6_pi+_q2_10_20_250000000.json
  - ./build/DEMPgen Config_EIC_10on250_ip6_pi+_q2_20_35_586000000.json

 - The .txt files generated from the execution of the commands above are recorded in a sub directory within this folder -
    - pion - pi+ commands, subdivided into each beam energy combo

- The commands above also produce the output in an optional root output format. These files are produced and backed up on multiple systems, but are not included on GitHub.
  - Contact Stephen Kay (stephen.kay@york.ac.uk) or Love Preet (navisaharan3@gmail.com) for access to these root files

- Once the files are generted with DEMPgen, they are processed through the Monte Carlo Afterburner (https://github.com/eic/afterburner) to incorporate crossing angle, beam  effects, and veterx spread for the EIC.

- This is done by running the following command within the ePIC container/eic-shell. 
  - 'abconv' runs the afterburner
    - The first argument is the input file path.
    - The last argument is the output file name, based on the given input file.
    - In this case, a specific beam energy combination configuration was enabled using the -p flag
    - The pion files were processed from the /work/eic/users/sjdkay/Jul2025_Campaign_Input directory (within eic-shell) -
      - abconv -p ip6_ep_130x10 Generator_Output/pion/10on130/eic_DEMPgen_10on130_ip6_Pi+_73000000_q2_3_10.hepmc3 -o Afterburner_Output/pion/10on130/DEMPgen_10on130_ip6_Pi+_q2_3_10_10x130
      - abconv -p ip6_ep_130x10 Generator_Output/pion/10on130/eic_DEMPgen_10on130_ip6_Pi+_170600000_q2_10_20.hepmc3 -o Afterburner_Output/pion/10on130/DEMPgen_10on130_ip6_Pi+_q2_10_20_10x130
      - abconv -p ip6_ep_130x10 Generator_Output/pion/10on130/eic_DEMPgen_10on130_ip6_Pi+_350000000_q2_20_35.hepmc3 -o Afterburner_Output/pion/10on130/DEMPgen_10on130_ip6_Pi+_q2_20_35_10x130
      - abconv -p ip6_ep_250x10 Generator_Output/pion/10on250/eic_DEMPgen_10on250_ip6_pi+_100600000_q2_3_10.hepmc3 -o Afterburner_Output/pion/10on250/DEMPgen_10on250_ip6_Pi+_q2_3_10_10x250
      - abconv -p ip6_ep_250x10 Generator_Output/pion/10on250/eic_DEMPgen_10on250_ip6_pi+_250000000_q2_10_20.hepmc3 -o Afterburner_Output/pion/10on250/DEMPgen_10on250_ip6_Pi+_q2_10_20_10x250
      - abconv -p ip6_ep_250x10 Generator_Output/pion/10on250/eic_DEMPgen_10on250_ip6_pi+_586000000_q2_20_35.hepmc3 -o Afterburner_Output/pion/10on250/DEMPgen_10on250_ip6_Pi+_q2_20_35_10x250
