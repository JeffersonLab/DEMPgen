# Stephen JD Kay (stephen.kay@york.ac.uk) - 06/05/2025 - University of York

To reproduce the files used for the May 2026 simulation campaign from DEMPgen, follow the steps below:

- Download the version 1.2.5 release of the generator from https://github.com/JeffersonLab/DEMPgen.

- After successfully compliling the genertor by following the instructions in the README.md file, remember to source setup.csh before running DEMPgen. Then, execute DEMPgen on the .json config files included in this directory. These were executed via -
  - ./build/DEMPgen Config_EIC_9on130_ip6_pi+_q2_3_10_76600000.json
  - ./build/DEMPgen Config_EIC_9on130_ip6_pi+_q2_10_20_170000000.json
  - ./build/DEMPgen Config_EIC_9on130_ip6_pi+_q2_20_35_365000000.json
  - ./build/DEMPgen Config_EIC_9on275_ip6_pi+_q2_3_10_120000000.json
  - ./build/DEMPgen Config_EIC_9on275_ip6_pi+_q2_10_20_321000000.json
  - ./build/DEMPgen Config_EIC_9on275_ip6_pi+_q2_20_35_726000000.json

 - The .txt files generated from the execution of the commands above are recorded in a sub directory within this folder -
    - pion - pi+ commands, subdivided into each beam energy combo

- The commands above also produce the output in an optional root output format. These files are produced and backed up on multiple systems, but are not included on GitHub.
  - Contact Stephen Kay (stephen.kay@york.ac.uk) for access to these root files

- Once the files are generted with DEMPgen, they are processed through the Monte Carlo Afterburner (https://github.com/eic/afterburner) to incorporate crossing angle, beam  effects, and veterx spread for the EIC.

- This is done by running the following command within the ePIC container/eic-shell. 
  - 'abconv' runs the afterburner
    - The first argument is the input file path.
    - The last argument is the output file name, based on the given input file.
    - In this case, a specific beam energy combination configuration was enabled using the -p flag. Note that the configurations below were specified as no dedicated flag for the new beam energy combinations in use here existed
    - The pion files were processed from the /work/eic/users/sjdkay/May2026_Campaign_Input directory (within eic-shell) -
      - abconv -p ip6_ep_130x9 Generator_Output/pion/9on130/eic_DEMPgen_9on130_ip6_pi+_q2_3_10_76600000.hepmc3 -o Afterburner_Output/pion/9on130/DEMPgen_v1.2.5_DEMP_Pi+_9x130_q2_3to10
      - abconv -p ip6_ep_130x9 Generator_Output/pion/9on130/eic_DEMPgen_9on130_ip6_pi+_q2_10_20_170000000.hepmc3 -o Afterburner_Output/pion/9on130/DEMPgen_v1.2.5_DEMP_Pi+_9x130_q2_10to20
      - abconv -p ip6_ep_130x9 Generator_Output/pion/9on130/eic_DEMPgen_9on130_ip6_pi+_q2_20_35_365000000.hepmc3 -o Afterburner_Output/pion/9on130/DEMPgen_v1.2.5_DEMP_Pi+_9x130_q2_20to35
      - abconv -p ip6_ep_275x9 Generator_Output/pion/9on275/eic_DEMPgen_9on275_ip6_pi+_q2_3_10_120000000.hepmc3 -o Afterburner_Output/pion/9on275/DEMPgen_v1.2.5_DEMP_Pi+_9x275_q2_3to10
      - abconv -p ip6_ep_275x9 Generator_Output/pion/9on275/eic_DEMPgen_9on275_ip6_pi+_q2_10_20_321000000.hepmc3 -o Afterburner_Output/pion/9on275/DEMPgen_v1.2.5_DEMP_Pi+_9x275_q2_10to20
      - abconv -p ip6_ep_275x9 Generator_Output/pion/9on275/eic_DEMPgen_9on275_ip6_pi+_q2_20_35_726000000.hepmc3 -o Afterburner_Output/pion/9on275/DEMPgen_v1.2.5_DEMP_Pi+_9x275_q2_20to35
- The latest version of afterburner at time of use (08/05/26) already produces output in the correct hepmc3.tree.root file format
