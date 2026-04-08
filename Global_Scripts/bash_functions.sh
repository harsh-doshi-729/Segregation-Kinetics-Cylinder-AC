#!/bin/bash
# This shell file contains functions that may be used in other bash scripts

# Reading confinement diameter functions from database:

# A function to read the values of the confinement cylinder diameters for various architectures
# and set them to a global dictionary named diameters
read_diameters(){
    filePath=$1
    while IFS="," read -r column1 column2 column3 # IFS is the Input Field Separator for our input; read command reads from input line by line
    # the -r option prevents any backslashes in the input to escape characters while reading
    do
        diameters["$column1"]="$column3" # skipping column 2 since it is the R_g
    done < <(tail -n +2 $filePath) # The input is taken from the process tail which returns the file from the second line; header skipped
    # using "process substitution to feed the modified file into the while loop and read command"
}

# Prints the diameters stored in the diameters dictionary
print_diameters(){
    echo "Confinement diameters for various architectures:"
    for key in ${!diameters[@]} # '!' allows iteration over the keys
    do
        echo "$key: ${diameters[$key]}"
    done
}

# A function to read the seeds for a particular architecture 
# Sets the read seeds in an array called "seeds"
read_seeds(){
    fileName=$1
    architecture=$2
    # Checking for arguments:
    if [ -z ${fileName} ] || [ -z ${architecture} ]; then
        echo "Not enough arguments passed to read_seeds(). Please pass the filePath and architecure name."
        return 1
    fi

    # Checking if file exists:
    if ! [ -e ${fileName} ]; then
        echo "The file ${fileName} could not be found! Terminating"
        return 1
    fi

    architectureFound=0 #  A flag to indicate whether the architecture has been found or not
    while read -r line # IFS is the Input Field Separator for our input; read command reads from input line by line
    # the -r option prevents any backslashes in the input to escape characters while reading
    do
        # Parsing the file line by line:
        # Ignoring comments:
        first_char=${line:0:1} # extracting first character
        # echo $first_char
        if [ "${first_char}" = "#" ] || [ "${first_char}" = "\n" ]; then # A comment string
            continue
        fi
        # else:
        # Checking if architecture was found in the previous line:
        if [ "${architectureFound}" == 1 ]; then
            # echo "Seeds: ${line}"
            # assign seeds:
            # Set IFS?
            seeds=( ${line} ) # tokenizes the string with the delimiter in IFS
            break
        fi
            
        if [ "${line}" = "${architecture}" ]; then
            # Architecture found!
            echo "Architecture ${architecture} found!"
            architectureFound=1
        fi

    done < "$fileName"

    if [ "${architectureFound}" == 0 ]; then
        echo "Architecture ${architecture} was not found in ${fileName}!"
        return 1
    fi

    return 0
}

# A function to calculate the ceil of a float
ceil() {                                              
    echo "define ceil (x) {if (x/1 == x) {return x/1} \
    else { if (x<0) {return x/1 -1} \
    else {return x/1 + 1 }}} ; ceil($1)" | bc
}