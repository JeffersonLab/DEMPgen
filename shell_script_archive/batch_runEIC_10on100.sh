#! /bin/bash                                                                                                        
# Batch submission job script, executes run_EIC_batch_10on100.csh script

echo "Running as ${USER}" # Checks who you're running this as
NumFiles=$1
NumEvents=$2
if [[ -z "$3" ]]; then
    echo "No output format specified, script will default to Pythia6 (Fun4All)"
    OutputType="Pythia6"
else OutputType=$3
fi
##Output history file##                                                                               
historyfile=hist.$( date "+%Y-%m-%d_%H-%M-%S" ).log # Creates a log file

i=1
while [[ $i -le $NumFiles ]]; do
    batch="${USER}_EICDempGen_10on100_${i}_Job.txt" # The name of the job submission script it'll create each time
    echo "Running ${batch} for file ${i}"
    cp /dev/null ${batch}
    RandomSeed=$(od -An -N3 -i /dev/random)
    echo "#!/bin/csh" >> ${batch} # Tells your job which shell to run in
    echo "#PBS -N DEMPGen_10on100_${NumEvents}_${i}" >> ${batch} # Name your job                     
    echo "#PBS -m abe" >> ${batch} # Email you on job start, end or error
    #echo "#PBS -M ${USER}@jlab.org" >>${batch} # Your email address, change it to be what you like
    echo "#PBS -r n" >> ${batch} # Don't re-run if it crashes
    echo "#PBS -o  /home/${USER}/trq_output/DEMPGen_10on100_${NumEvents}_${i}.out" >> ${batch} # Output directory and file name, set to what you like
    echo "#PBS -e  /home/${USER}/trq_output/DEMPGen_10on100_${NumEvents}_${i}.err" >> ${batch} # Error output directory and file name
    echo "date" >> ${batch} 
    echo "cd /home/apps/DEMPgen/" >> ${batch} # Tell your job to go to the directory with the script you want to run
    echo "./run_EIC_batch_10on100.csh ${i} ${NumEvents} ${RandomSeed} ${OutputType}" >> ${batch} # Run your script, change this to what you like
    echo "date">>${batch}
    echo "exit">>${batch} # End of your job script
    echo "Submitting batch"
    eval "qsub ${batch} 2>/dev/null" # Use qsub to actually submit your job
    echo " "
    i=$(( $i + 1 ))
    sleep 2
    rm ${batch}
done
