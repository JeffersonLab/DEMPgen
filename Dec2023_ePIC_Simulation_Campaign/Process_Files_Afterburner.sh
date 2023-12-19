#! /bin/bash

### Stephen Kay, University of Regina
### 07/03/22
### stephen.kay@uregina.ca
### A script to process a bunch of HEPMC3 DEMPGEN files through afterburner
### Required input is - $FileList $OutputDir $OutputNaming $IP
### Must be run in eic-shell/container for ATHENA/ePIC (or otherwise have afterburner set up)

# 18/10/22 - Updates, added some lines to check start/end of output naming to reduce weird file naming (__ in the output etc)
# Added logic statement to prevent processing the file if it already exists

echo "Running as ${USER}"
FileList=$1 # First arg is a file list, this should be a plain text file where each line is a path to a file
if [[ -z "$1" ]]; then
    echo "I need a file list to process!"
    echo "Please provide a file list as the first argument"
    exit 1
fi
OutputDir=$2 # Second argument is an output directory, this should be a path to where you want to output your files
if [[ -z "$2" ]]; then
    echo "I need an output directory!"
    echo "Please provide an output directory as the second argument"
    exit 2
fi

OutputNaming=$3 # Third argument is the output naming convention, this should be a string that the script can iteratively write files as
if [[ -z "$3" ]]; then
    echo "I need an output naming convention!"
    echo "Please provide an output naming convention as the third argument!"
    exit 3
fi
# Check if an argument was provided for the IP, if not, assume IP6
if [[ -z "$4" ]]; then
    IP="IP6"
    echo "IP argument not provided, assuming IP6"
else
    IP=$4
fi

# Check the file list and output directory exist, warn and exit if they do not
if [ ! -f $FileList ]; then
    echo "!!! WARNING !!!"
    echo "!!! $FileList - Does not exist - Double check pathing and try again !!!"
    echo "!!! WARNNING !!!"
    exit 4
fi
if [ ! -d $OutputDir ]; then
    echo "!!! WARNING !!!"
    echo "!!! $OutputDir - Does not exist - Double check pathing and try again !!!"
    echo "!!! WARNNING !!!"
    exit 4
fi

if [[ $IP != "IP6" && $IP != "IP8" ]]; then
    echo "!!! NOTICE !!!"
    echo "Expected IP6 or IP8 to be specified as IP value, defaulting to IP6"
    echo "!!! NOTICE !!!"
    IP="IP6"
fi

# Do some quick editing of our output directories and output naming schemes
if [ "${OutputDir: -1}" == "/" ]; then # If output directory ends in a /, snip it off
    OutputDir=${OutputDir::len-1}
fi
if [ "${OutputNaming: -1}" == "_" ] || [ "${OutputNaming: -1}" == "/" ]; then # If the output naming convention provided ends in a _ or /, snip it off
    OutputNaming=${OutputNaming::len-1} 
fi
if [[ $OutputNaming = /* ]]; then # If output naming convention begins with a /, snip it off
    OutputNaming=${OutputNaming:1} 
fi

i=1 # An iterator we'll use to name files and count our loop in a moment
while IFS='' read -r line || [[ -n "$line" ]]; do # Read the file list, line by line
    File=$line
    if [ ! -f $File ]; then # Check file exists, if it doesn't, skip to the next line in the file
	echo "Entry on line $i - $File - Not found, skipping"
	continue
    fi
    # Set up the output file name
    OutputFilename="${OutputDir}/${OutputNaming}_${i}"
    if [ ! -f "${OutputFilename}.hist.root" ]; then # Only process if the file doesn't already exist
	# Process our files through afterburner
	if [ $IP == "IP6" ]; then # If IP6, process with IP6 (default) preset
	    eval abconv $File -o $OutputFilename
	else # If IP8 specified, process with IP8 preset
	    eval abconv $File -o $OutputFilename -p 3
	fi
	mv "${OutputFilename}.hist.root" "${OutputFilename}.hepmc.tree.root"
    else 
	echo "${OutputFilename}.hepmc and .root already exist, processing skipped"
    fi
    i=$(( $i + 1)) # Add to our iterator
done < "$FileList"

exit 0
