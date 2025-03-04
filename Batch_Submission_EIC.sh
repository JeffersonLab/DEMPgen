#! /bin/bash                            
                                                                            
# SJDK - 09/02/22 - An updated version of the batch submission script. This script is designed for the torque queueing system on Lark at the University of Regina. However, it could quickly be adapted for use on the JLab iFarm for example.
# This script creates batch job files and submits them, these jobs run the Process_EIC.csh script

echo "Running as ${USER}" # Checks who you're running this as
# If output directory for .out/.err files doesn't exist, make it
if [ ! -d  "/home/${USER}/trq_output" ]; then
    echo "/home/${USER}/trq_output directory doesn't exist, making this directory for you..."
    mkdir "/home/${USER}/trq_output"
    echo "Directory created, check there for output and error logs from your job."
fi
# If 7 or 8 arguments not given, complain
if [[ "$#" -ne 7 && "$#" -ne 8 ]]; then
    echo ""
    echo "!!! ERROR !!! - Expected 7 or 8 arguments - !!! ERROR !!!"
    echo "Expect - NumFiles NumEvents EBeamE HBeamE OutputType InteractionPoint Ejectile RecoilHadron(optional)"
    echo "See the Config_EIC.json file or the README for options and try again, exiting"
    echo "!!! ERROR !!! - Expected 7 or 8 arguments - !!! ERROR !!!"
    echo ""
    exit 0
fi

# Set variables equal to arguments provided
NumFiles=$1
NumEvents=$2
EBeamE=$3
HBeamE=$4
OutputType=$5
InteractionPoint=$6
Ejectile=$7

# If K+ specified, check the 8th argument, expect this to exist for K+, if it does NOT (case 1), set a default
if [[ $Ejectile == "K+" && -z "$8" ]]; then
    echo "!!! WARNING !!! - For K+ production expect a hadron specified, defaulting to Lambda - !!! WARNING !!!"
    RecoilHadron="Lambda"
elif [[ $Ejectile == "K+" && ! -z "$8" ]]; then # If 8th argument is not a blank string (i.e. it exists), set the RecoilHadron to this
    RecoilHadron=$8
else # Any other case (non K+), set RecoilHadron to be a blank string. We don't actually care for Pi+, Pi0 production etc.
    RecoilHadron=""
fi

i=1
while [[ $i -le $NumFiles ]]; do
    # This is the name of the job submission script the shell script creates
    batch="${USER}_EICDempGen_${EBeamE}on${HBeamE}_${Ejectile}${RecoilHadron}_${InteractionPoint}_${NumEvents}_${i}_Job.txt" # The name of the job submission script it'll create each time
    echo "Running ${batch} for file ${i}"
    cp /dev/null ${batch}
    RandomSeed=$(od -An -N3 -i /dev/urandom)
    echo "#!/bin/csh" >> ${batch} # Tells your job which shell to run in
    echo "#PBS -N DEMPgen_${EBeamE}on${HBeamE}_${Ejectile}${RecoilHadron}_${InteractionPoint}_${NumEvents}_${i}" >> ${batch} # Name your job                     
    echo "#PBS -m abe" >> ${batch} # Email you on job start, end or error
    #echo "#PBS -M ${USER}@jlab.org" >>${batch} # Your email address, change it to be what you like
    echo "#PBS -r n" >> ${batch} # Don't re-run if it crashes
    echo "#PBS -o  /home/${USER}/trq_output/${EBeamE}on${HBeamE}_${Ejectile}${RecoilHadron}_${InteractionPoint}_${NumEvents}_${i}.out" >> ${batch} # Output directory and file name, set to what you like
    echo "#PBS -e  /home/${USER}/trq_output/${EBeamE}on${HBeamE}_${Ejectile}${RecoilHadron}_${InteractionPoint}_${NumEvents}_${i}.err" >> ${batch} # Error output directory and file name
    echo "date" >> ${batch} 
    echo "cd /home/apps/DEMPgen/" >> ${batch} # Tell your job to go to the directory with the script you want to run
    echo "./Process_EIC.csh ${i} ${NumEvents} ${EBeamE} ${HBeamE} ${RandomSeed} ${OutputType} ${InteractionPoint} ${Ejectile} ${RecoilHadron}" >> ${batch} # Run your script, change this to what you like
    echo "date">>${batch}
    echo "exit">>${batch} # End of your job script
    echo "Submitting batch"
    eval "qsub ${batch} 2>/dev/null" # Use qsub to actually submit your job
    echo " "
    i=$(( $i + 1 ))
    sleep 2
    rm ${batch}
done

