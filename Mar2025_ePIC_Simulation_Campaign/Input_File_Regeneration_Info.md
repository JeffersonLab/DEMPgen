# Stephen JD Kay (stephen.kay@york.ac.uk)- 04/03/2025 - University of York

To reproduce the files used for the March 2025 simulation campaign from DEMPgen, follow the steps below:

- Download the version 1.2.3 release of the generator from https://github.com/JeffersonLab/DEMPgen.

- After successfully compliling the genertor by following the instructions in the README.md file, DEMPgen was executed on the .json config files included in this directory. These were executed via -
  - ./build/DEMPgen Config_EIC_10on130_ip6_Pi+_q2_3_10_11000000.json
  - ./build/DEMPgen Config_EIC_10on130_ip6_Pi+_q2_10_20_18000000.json
  - ./build/DEMPgen Config_EIC_10on130_ip6_Pi+_q2_20_35_34000000.json

 - The .txt files generated from the execution of the commands above are recorded in a sub directory within this folder -
    - pion - pi+ commands, subdivided into each beam energy combo (only 10x130 was run this time)

- The commands above also produce the output in an optional root output format. These files are produced and backed up on multiple systems, but are not included on GitHub.
  - Contact Stephen Kay (stephen.kay@york.ac.uk) or Love Preet (navisaharan3@gmail.com) for access to these root files

- The output files are initially in .dat format, this was renamed to .hepmc3 for input into afterburner.

- Once the files are generted with DEMPgen, they are processed through the Monte Carlo Afterburner (https://github.com/eic/afterburner) to incorporate crossing angle, beam  effects, and veterx spread for the EIC.

- This is done by running the following command within the ePIC container/eic-shell. 
  - abconv 
    - The first argument is the input file path. The corresponding files are provided in the Afterburner_Filelists directory for each beam energy combination within this folder.
    - The last argument is the output file name, based on the given input file.
    - In this case, a specific configuration was enabled using the -p flag
    - The files were processed with the 10x100 hidiv setting, from /work/eic/users/sjdkay/Mar2025_Campaign_Input (within eic-shell) -
      - abconv -p ip6_hidiv_100x10 Generator_Output/pion/10on130/eic_input_DEMPgen_10on130_ip6_Pi+_11000000_q2_3_10.hepmc3 -o Afterburner_Output/pion/10on130/DEMPgen_10on130_ip6_Pi+_q2_3_10_10x100_hidiv.root
      - abconv -p ip6_hidiv_100x10 Generator_Output/pion/10on130/eic_input_DEMPgen_10on130_ip6_Pi+_18000000_q2_10_20.hepmc3 -o Afterburner_Output/pion/10on130/DEMPgen_10on130_ip6_Pi+_q2_10_20_10x100_hidiv.root
      - abconv -p ip6_hidiv_100x10 Generator_Output/pion/10on130/eic_input_DEMPgen_10on130_ip6_Pi+_34000000_q2_20_35.hepmc3 -o Afterburner_Output/pion/10on130/DEMPgen_10on130_ip6_Pi+_q2_20_35_10x100_hidiv.root
    - The files were also processed with the 10x100 hiacc setting, from /work/eic/users/sjdkay/Mar2025_Campaign_Input (within eic-shell) -
      - abconv -p ip6_hiacc_100x10 Generator_Output/pion/10on130/eic_input_DEMPgen_10on130_ip6_Pi+_11000000_q2_3_10.hepmc3 -o Afterburner_Output/pion/10on130/DEMPgen_10on130_ip6_Pi+_q2_3_10_10x100_hiacc.root
      - abconv -p ip6_hiacc_100x10 Generator_Output/pion/10on130/eic_input_DEMPgen_10on130_ip6_Pi+_18000000_q2_10_20.hepmc3 -o Afterburner_Output/pion/10on130/DEMPgen_10on130_ip6_Pi+_q2_10_20_10x100_hiacc.root
      - abconv -p ip6_hiacc_100x10 Generator_Output/pion/10on130/eic_input_DEMPgen_10on130_ip6_Pi+_34000000_q2_20_35.hepmc3 -o Afterburner_Output/pion/10on130/DEMPgen_10on130_ip6_Pi+_q2_20_35_10x100_hiacc.root

    
