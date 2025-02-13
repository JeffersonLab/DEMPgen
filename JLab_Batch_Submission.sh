#! /bin/bash                            
                                                                            
# SJDK - 07/03/22 - A version of the batch submission script that is intended to run on the JLab iFarm
# This script creates batch job files and submits them, these jobs run the Process_EIC.csh script

# This has been tested and works successfully

echo "Running as ${USER}" # Checks who you're running this as

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

Timestamp=$(date +'%d_%m_%Y')
Workflow="EIC_DEMPgen_${USER}_${Timestamp}" # Change this as desired

while true; do
    read -p "Do you wish to begin a new batch submission? (Please answer yes or no) " yn
    case $yn in
        [Yy]* )
            i=1
            (
		while [[ $i -le $NumFiles ]]; do
		    # This is the name of the job submission script the shell script creates
		    batch="${USER}_EICDempGen_${EBeamE}on${HBeamE}_${Ejectile}${RecoilHadron}_${InteractionPoint}_${NumEvents}_${i}_Job.txt" # The name of the job submission script it'll create each time
		    echo "Running ${batch} for file ${i}"
		    cp /dev/null ${batch}
		    RandomSeed=$(od -An -N3 -i /dev/urandom)
		    echo "PROJECT: c-kaonlt"  >> ${batch} # Is eic a valid project?
		    echo "TRACK: analysis" >> ${batch}
		    echo "JOBNAME: DEMPGen_${EBeamE}on${HBeamE}_${Ejectile}${RecoilHadron}_${InteractionPoint}_${NumEvents}_${i}" >> ${batch}
                    echo "MEMORY: 2000 MB" >> ${batch} # Request 2GB RAM - probably too much
		    echo "CPU: 1" >> ${batch} # Request 1 CPU core per job
		    echo "COMMAND:/group/eic/users/${USER}/DEMPGen/Process_EIC_iFarm.csh ${i} ${NumEvents} ${EBeamE} ${HBeamE} ${RandomSeed} ${OutputType} ${InteractionPoint} ${Ejectile} ${RecoilHadron}" >> ${batch}
                    echo "MAIL: ${USER}@jlab.org" >> ${batch}
		    echo "Submitting batch"
		    eval "swif2 add-jsub ${Workflow} -script ${batch} 2>/dev/null" # Swif2 job submission, uses old jsub scripts
		    echo " "
		    i=$(( $i + 1 ))
		    sleep 2
		    rm ${batch}
		    if [ $i -gt $NumFiles ]; then
			echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
			echo " "
			echo "###############################################################################################################"
			echo "############################################ END OF JOB SUBMISSIONS ###########################################"
			echo "###############################################################################################################"
			echo " "
		    fi
		done
	    )
	    eval 'swif2 run ${Workflow}'
	    break;;
        [Nn]* ) 
	    exit;;
        * ) echo "Please answer yes or no.";;
    esac
done

