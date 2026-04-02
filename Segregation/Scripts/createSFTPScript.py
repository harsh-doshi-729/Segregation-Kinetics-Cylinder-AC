import numpy as np
import sys
import io
# importing the scripts that contains all the local and cluster paths
# The path to the directory containing this module is set as the environment variable PYTHONPATH
import system_file_paths as sysPaths
import cluster_file_paths as clusPaths
import Segregation_Parameters as segParam

filePath = f"SFTP_script.txt" 
architecture = "Arc0"
numberOfRuns = 50
numberOfMonomers = 200

modes = ["initial", "mixed", "distributions", "post", "concatenation", "segregation", "put-positions", "backup", "custom"]
chosenMode = "segregation"
special_simulation = ""

def SetConstants() -> None:
    """Reads the argument(s) passed while invoking the script from the command line and sets it to the architecture and number of Monomers"""
    global architecture
    global numberOfMonomers
    global chosenMode
    global special_simulation

    numberOfMandatoryArguments = 3
    if len(sys.argv) <= numberOfMandatoryArguments:
        print("Not enough arguments passed while invoking the script! Please pass the number of monomers, the name of the architecture, and the write mode in the following format:")
        print("python /<path>/plotCoMDistribution.py <numberOfMonomers> <architecture> <write-mode>")
        print(f"Accepted modes: {modes}")
        # print(f"An optional argument for the special simulation name can be passed to work in the corresponding named directory.")
        quit() # terminating the script
    
    architecture = sys.argv[2]
    numberOfMonomers = sys.argv[1]
    if numberOfMonomers.isdigit():
        numberOfMonomers = int(numberOfMonomers)
    else:
        print("The first argument passed was not an integer! Please pass either 200 or 500 for the number of monomers.")

    chosenMode = sys.argv[3]
    if not chosenMode in modes:
        print(f"The passed mode is not recognized. Please pass one of the following: {modes}. Terminating")
        quit()

    if len(sys.argv) > numberOfMandatoryArguments + 1: # Optional argument passed
        special_simulation = sys.argv[numberOfMandatoryArguments+1]

def GetFolder(baseFolder: str, numberOfMonomers: int, architecture: str, special_simulation: str|None = None) -> str:
    """Returns the local or cluster folder by incorporating the special_simulation, number of monomers, and architecture"""
    if special_simulation is None:
        special_simulation = sysPaths.SPECIAL_SIMULATION
    return f"{sysPaths.GetFolder(baseFolder, numberOfMonomers, special_simulation)}{architecture}/"

def WriteSingleRunCommands(file: io.TextIOWrapper, runIndex: int, commands: str) -> None:
    """Writes the common commands to enter and exit a run folder incorporated with the specialised commands passed"""
    file.write(f"cd run{runIndex}\n")
    file.write(f"lmkdir run{runIndex}\n") # if the folder doesn't already exist
    file.write(f"lcd run{runIndex}\n")
    file.write(f"{commands}")
    file.write("cd ..\n")
    file.write("lcd ..\n")

def WritePutInitialConfigurations() -> None:
    """Writes the commands to put the random initial mixed states for different runs to the corresponding run folder on the cluster"""
    spath = GetFolder(clusPaths.CREATE_INITIAL_STATES, numberOfMonomers, architecture) # remote path
    localPath = GetFolder(sysPaths.CREATE_INITIAL_STATES, numberOfMonomers, architecture)
    with open(filePath, "w") as file:
        file.write(f"cd {spath}\n")
        file.write(f"lcd {localPath}\n")
        for i in range(1, numberOfRuns+1):
            file.write(f"mkdir run{i}\n")
            file.write(f"cd run{i}/\n")
            file.write(f"put initial_configuration{i}.txt initial_configuration.txt\n")
            file.write("cd ..\n")
        file.write("exit\n")

def WriteGetMixedStates() -> None:
    """Writes commands to get the initial mixed states and the CoM files from each run folder"""
    spath = GetFolder(clusPaths.CREATE_INITIAL_STATES, numberOfMonomers, architecture) # remote path
    with open(filePath, "w") as file:
        file.write(f"cd {spath}\n")
        # file.write(f"lcd {architecture}/ ") # going to local directory
        file.write("get mixed_state_*\n")
        for i in range(1, numberOfRuns+1):
            WriteSingleRunCommands(file, i, "get com*\n")
        file.write("exit\n")

def WriteGetShrinkRelaxData(singleRun: bool = True) -> None:
    """Writes commands to get all the shrink-relax data"""
    global numberOfRuns
    spath = GetFolder(clusPaths.CREATE_INITIAL_STATES, numberOfMonomers, architecture) # remote path
    localPath = GetFolder(sysPaths.CREATE_INITIAL_STATES, numberOfMonomers, architecture) # local path
    with open(filePath, "w") as file:
        file.write(f"cd {spath}\n")
        file.write(f"lcd {localPath}\n") # going to local directory
        # if singleRun:
        #     numberOfRuns = 1
        # commands = "get -r *\n"
        # for i in range(1, numberOfRuns + 1):
        #     WriteSingleRunCommands(file, i, commands)
        file.write("get -r *\n")
        file.write("exit\n")


def WriteGetPostCoM(singleRun: bool = True) -> None:
    """Writes commands to get the post mixing CoM files"""
    spath = GetFolder(clusPaths.CREATE_INITIAL_STATES, numberOfMonomers, architecture) # remote path
    localPath = GetFolder(sysPaths.CREATE_INITIAL_STATES, numberOfMonomers, architecture) # local path
    with open(filePath, "w") as file:
        file.write(f"cd {spath}\n")
        file.write(f"lcd {localPath}\n") # going to local directory
        if singleRun:
            numberOfRuns = 1
        commands = """get post_com*
        get com_reg*\n"""
        for i in range(1, numberOfRuns + 1):
            WriteSingleRunCommands(file, i, commands)
        file.write("exit\n")

def WriteGetConcatenatedStates() -> None:
    """Writes commands to get the concatenation verifying simulation dump"""
    spath = GetFolder(clusPaths.CHECK_CONCATENATION, numberOfMonomers, architecture) # remote path
    localPath = GetFolder(sysPaths.CHECK_CONCATENATION, numberOfMonomers, architecture) # local path
    with open(filePath, "w") as file:
        file.write(f"cd {spath}\n")
        file.write(f"lcd {localPath}\n") # going to local directory
        commands = """get both_log.lammps
        get com*
        get separate_log.lammps\n"""
        for i in range(1, numberOfRuns+1): # getting the dumps for all states
            WriteSingleRunCommands(file, i, commands)
        file.write("exit\n")

def WriteGetSegregationData() -> None:
    """Writes commands to get all the data generated during the segregation runs"""
    spath = GetFolder(clusPaths.NEW_SEGREGATION, numberOfMonomers, architecture) # remote path
    localPath = GetFolder(sysPaths.NEW_SEGREGATION, numberOfMonomers, architecture) # local path
    with open(filePath, "w") as file:
        file.write(f"cd {spath}\n")
        file.write(f"lcd {localPath}\n") # going to local directory
        commands = """get *\n"""
        for i in range(1, numberOfRuns + 1):
            WriteSingleRunCommands(file, i, commands)
        file.write("exit\n")

def WritePutPositionData() -> None:
    """Writes commands to put the Create_Initial_States distribution_positions in the cluster. This is useful for performing Analysis on the cluster."""
    spath = GetFolder(clusPaths.CREATE_INITIAL_STATES, numberOfMonomers, architecture) # remote path
    localPath = GetFolder(sysPaths.CREATE_INITIAL_STATES, numberOfMonomers, architecture)
    with open(filePath, "w") as file:
        file.write(f"cd {spath}\n")
        file.write(f"lcd {localPath}\n")
        for i in range(1, numberOfRuns+1):
            file.write(f"mkdir run{i}\n")
            file.write(f"cd run{i}/\n")
            file.write(f"lcd run{i}/\n")
            file.write(f"put distribution_positions.dump\n")
            file.write("cd ..\n")
            file.write("lcd ..\n")
        file.write("exit\n")

def WriteBackupScript() -> None:
    """Writes commands to get the .zip files from a set of special simulations and architectures defined in Segregation_Parameters.py"""
    # Printing the special simulations and architectures that will be used:
    print(f"Getting files from the following special simulations: {segParam.BACKUP_SPECIAL_SIMULATIONS}")
    print(f"Getting files from the following architectures: {segParam.BACKUP_ARCHITECTURES}")
    with open(filePath, "w") as file:
        for special_simulation in segParam.BACKUP_SPECIAL_SIMULATIONS:
            file.write(f"lcd {segParam.SFTP_LOCAL_BASE_FOLDER}b{numberOfMonomers}/Previous_Attempts/\n")
            file.write(f"lmkdir {special_simulation}\n") # Making folder if it doesn't already exist
            file.write(f"lcd {special_simulation}\n")
            for arc in segParam.BACKUP_ARCHITECTURES:
                file.write(f"lmkdir {arc}\n") # Making folder if it doesn't already exist
                file.write(f"lcd {arc}\n")
                spath = GetFolder(segParam.SFTP_S_BASE_FOLDER, numberOfMonomers, arc, special_simulation)
                file.write(f"cd {spath}\n")
                for i in range(1, numberOfRuns + 1):
                    commands = "get *.zip\n"
                    WriteSingleRunCommands(file, i, commands)
                file.write("lcd ..\n")
                file.write("cd ..\n")
        file.write("exit\n")

            

def WriteCustomScript() -> None:
    """Writes a custom set of commands as specified in the Segregation_Parameters.py"""
    spath = GetFolder(segParam.SFTP_S_BASE_FOLDER, numberOfMonomers, architecture, special_simulation) # remote path
    localPath = GetFolder(segParam.SFTP_LOCAL_BASE_FOLDER, numberOfMonomers, architecture, special_simulation)
    with open(filePath, "w") as file:
        file.write(f"cd {spath}\n")
        file.write(f"lcd {localPath}\n")
        for i in range(1, numberOfRuns + 1):
            file.write(segParam.Get_SFTP_Commands(i))
        file.write("exit\n")

# Scripting:
if __name__ == "__main__":
    SetConstants() # setting architecture and number of monomers
    if chosenMode == "initial":
        WritePutInitialConfigurations()
    elif chosenMode == "mixed":
        WriteGetMixedStates()
    elif chosenMode == "distributions":
        WriteGetShrinkRelaxData(False)
    elif chosenMode == "post":
        WriteGetPostCoM()
    elif chosenMode == "concatenation":
        WriteGetConcatenatedStates()
    elif chosenMode == "segregation":
        WriteGetSegregationData()
    elif chosenMode == "put-positions":
        WritePutPositionData()
    elif chosenMode == "backup":
        WriteBackupScript()
    elif chosenMode == "custom":
        WriteCustomScript()
    else:
        print("The corresponding function for the passed node could not be found!")