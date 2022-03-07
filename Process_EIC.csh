#!/bin/csh

# SJDK - 09/02/22 - New script which takes in a whole bunch of inputs to create jobs/config files. Note that I'm not 100% happy with the pathing in this file so it should be tweaked and optimised at some point

echo""
echo "This file is intended to be run as part of a batch job submission, however, you can also run it on its own."
echo "Expected input is - FileNumber NumberOfEvents ElectronBeamEnergy HadronBeamEnergy RandomSeed OutputType IP Particle Hadron(Optional, for K+)"
echo "Please see the README for more info."
echo ""

if ($#argv != 8 && $#argv != 9) then
    echo "! ERROR !"
    echo "! ERROR ! - Expected 8 or 9 arguments, please see the opening information -! ERROR !"
    echo "! ERROR !"
    exit 0
endif

set FileNum=$1
set NumEvents=$2
set EBeamE=$3
set HBeamE=$4
set RandomSeed=$5
set OutputType=$6
set InteractionPoint=$7
set Particle=$8

if ($Particle == "K+" && $#argv == 8 ) then
    echo "! WARNING !"
    echo "! WARNING ! - For K+ production expect a hadron specified, defaulting to Lambda - ! WARNING !"
    echo "! WARNING !"
    set Hadron="Lambda"
else if ($Particle == "K+" && $#argv == 9 ) then
    set Hadron=$9
else 
    set Hadron=""
endif

echo "Running target polarisation up, FF setting for file $FileNum with $NumEvents events per file for $EBeamE GeV e- on $HBeamE GeV p using random seed $RandomSeed, using $OutputType format output for $Particle $Hadron events."
    
# Set the config file name based upon inputs
set ConfigFilename = 'Config_EIC_'$EBeamE'on'$HBeamE'_'$InteractionPoint'_'$Particle$Hadron'_'$NumEvents'_'$FileNum'.json'

# Copy the default config file to our constructed filename
cp Config_EIC.json $ConfigFilename

# Use sed commands to change our config file based upon inputs
sed -i 's/"file_name" \:.*/"file_name" \: "DEMPGen_'$EBeamE'on'$HBeamE'_'$InteractionPoint'_'$Particle$Hadron'_'$NumEvents'_'$FileNum'",/' $ConfigFilename
sed -i 's/"n_events" \:.*/"n_events" \: '$NumEvents',/' $ConfigFilename
sed -i 's/"generator_seed"\:.*/"generator_seed" \: '$RandomSeed',/' $ConfigFilename
sed -i 's/"ebeam"\:.*/"ebeam" \: '$EBeamE',/' $ConfigFilename
sed -i 's/"hbeam"\:.*/"hbeam" \: '$HBeamE',/' $ConfigFilename
sed -i 's/"particle"\:.*/"particle" \: "'$Particle'",/' $ConfigFilename
sed -i 's/"hadron"\:.*/"hadron" \: "'$Hadron'",/' $ConfigFilename
sed -i 's/"det_location"\:.*/"det_location" \: "'$InteractionPoint'",/' $ConfigFilename
sed -i 's/"OutputType"\:.*/"OutputType"\: "'$OutputType'",/' $ConfigFilename

# Run our new config file
cd data/
./../build/DEMPgen ../$ConfigFilename
sleep 5

# Filename as it's created is a bit odd, so rename it
set OriginalOutput = 'eic_input_DEMPGen_'$EBeamE'on'$HBeamE'_'$InteractionPoint'_'$Particle$Hadron'_'$NumEvents'_'$FileNum'.dat'
set RenamedOutput = 'eic_DEMPGen_'$EBeamE'on'$HBeamE'_'$InteractionPoint'_'$Particle$Hadron'_'$NumEvents'_'$FileNum'.dat'
mv "LundFiles/"$OriginalOutput "LundFiles/"$RenamedOutput

rm -rf ../$ConfigFilename
