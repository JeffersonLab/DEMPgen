#!/bin/csh

set FileNum=$1
set NumEvents=$2
set RandomSeed=$3

if ( ! $?4 ) then
    echo "Output type not specified, defaulting to Pythia6"
    set OutputType="Pythia6"
else
    set OutputType=$4
endif

echo "Running target polarisation up, FF setting for file $FileNum with $NumEvents events per file using random seed $RandomSeed using $OutputType output format."
cp Config_EIC.json Config_EIC_10on100_$FileNum.json
sed -i 's/"file_name" \:.*/"file_name" \: "DEMPGen_10on100_'$NumEvents'_'$FileNum'",/' Config_EIC_10on100_$FileNum.json
sed -i 's/"n_events" \:.*/"n_events" \: '$NumEvents',/' Config_EIC_10on100_$FileNum.json
sed -i 's/"generator_seed"\:.*/"generator_seed" \: '$RandomSeed',/' Config_EIC_10on100_$FileNum.json
sed -i 's/"ebeam"\:.*/"ebeam" \: '10.0',/' Config_EIC_10on100_$FileNum.json
sed -i 's/"hbeam"\:.*/"hbeam" \: '100.0',/' Config_EIC_10on100_$FileNum.json
sed -i 's/"OutputType"\:.*/"OutputType"\: "'$OutputType'",/' Config_EIC_10on100_$FileNum.json
cd data/
./../build/DEMPgen ../Config_EIC_10on100_$FileNum.json
sleep 5
mv "LundFiles/eic_input_DEMPGen_10on100_"$NumEvents"_"$FileNum".dat" "LundFiles/eic_DEMPGen_10on100_"$NumEvents"_"$FileNum".dat"
rm -rf ../Config_EIC_10on100_$FileNum.json
