#!/bin/csh

# SJDK - 09/02/22 - New script which takes in a whole bunch of inputs to create jobs/config files. Note that I'm not 100% happy with the pathing in this file so it should be tweaked and optimised at some point

source setup.csh

echo""
echo "This file is intended to be run as part of a batch job submission, however, you can also run it on its own."
echo "Expected input is - FileNumber NumberOfEvents ElectronBeamEnergy HadronBeamEnergy RandomSeed OutputType IP Ejectile RecoilHadron(Optional, for K+)"
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
set Ejectile=$8

if ($Ejectile == "K+" && $#argv == 8 ) then
    echo "! WARNING !"
    echo "! WARNING ! - For K+ production expect a hadron specified, defaulting to Lambda - ! WARNING !"
    echo "! WARNING !"
    set RecoilHadron="Lambda"
else if ($Ejectile == "K+" && $#argv == 9 ) then
    set RecoilHadron=$9
else 
    set RecoilHadron=""
endif

echo "Running file $FileNum with $NumEvents events per file for $EBeamE GeV e- on $HBeamE GeV p using random seed $RandomSeed, using $OutputType format output for $Ejectile $RecoilHadron events."
    
# Set the config file name based upon inputs
set ConfigFilename = 'Config_EIC_'$EBeamE'on'$HBeamE'_'$InteractionPoint'_'$Ejectile$RecoilHadron'_'$NumEvents'_'$FileNum'.json'

# Copy the default config file to our constructed filename
cp Config_EIC.json $ConfigFilename

# Use sed commands to change our config file based upon inputs
sed -i 's/"file_name" \:.*/"file_name" \: "DEMPgen_'$EBeamE'on'$HBeamE'_'$InteractionPoint'_'$Ejectile$RecoilHadron'_'$NumEvents'_'$FileNum'",/' $ConfigFilename
sed -i 's/"n_events" \:.*/"n_events" \: '$NumEvents',/' $ConfigFilename
sed -i 's/"generator_seed"\:.*/"generator_seed" \: '$RandomSeed',/' $ConfigFilename
sed -i 's/"ebeam"\:.*/"ebeam" \: '$EBeamE',/' $ConfigFilename
sed -i 's/"hbeam"\:.*/"hbeam" \: '$HBeamE',/' $ConfigFilename
sed -i 's/"ejectile"\:.*/"ejectile" \: "'$Ejectile'",/' $ConfigFilename
sed -i 's/"recoil_hadron"\:.*/"recoil_hadron" \: "'$RecoilHadron'",/' $ConfigFilename
sed -i 's/"det_location"\:.*/"det_location" \: "'$InteractionPoint'",/' $ConfigFilename
sed -i 's/"OutputType"\:.*/"OutputType"\: "'$OutputType'",/' $ConfigFilename

# Run our new config file
./build/DEMPgen $ConfigFilename
sleep 5

# Filename as it's created is a bit odd, so rename it
set OriginalOutput = 'eic_input_DEMPgen_'$EBeamE'on'$HBeamE'_'$InteractionPoint'_'$Ejectile$RecoilHadron'_'$NumEvents'_'$FileNum'.dat'
set RenamedOutput = 'eic_DEMPgen_'$EBeamE'on'$HBeamE'_'$InteractionPoint'_'$Ejectile$RecoilHadron'_'$NumEvents'_'$FileNum'.dat'
mv "data/output/"$OriginalOutput "data/output/"$RenamedOutput

rm -rf ../$ConfigFilename
