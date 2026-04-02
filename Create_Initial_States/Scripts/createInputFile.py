# creating the LAMMPS data file for Arc0 or Arc_Lin architecture:

from matplotlib.offsetbox import PaddedBox
import numpy as np
import math
import sys
import io
from typing import Tuple
from typing import List

# indexOfPolymer = 2 # the index of which Arc0 polymer is being created
numberOfPolymers = 2
numberofAtomTypes = 2
numberOfMonomers = 200 # In one polymer
architecture = "Arc2"
numberOfBonds = numberOfMonomers
offsetAngles = []
numberOfRuns = 1 # The number of initial files to be written; these file will only be different if atoms are placed randomly

# box information: cylinder
WALL = 0 # width of the box wall; in LJ units; half of sigma
box_length = 25 
diameter = 6.24 # in units of sigma; This value alone determines the dimensions of the box; set it to the desired value
side_length = 2 * (2.5) # in units of sigma; the side length of the cuboid box housing the cyclinder (includes the wall width)
box_dimensions = []
isBoxCentredAtOrigin = True # A flag to indicate whether the box is centred at the origin or has a corner at the origin; this will determine the range of coordinates for the atoms

# cross link:
insertTerCrossLink = False
crossLinks = [] # TODO: implement a function to read cross links from the database

def SetCylinderBoxDimensions() -> None:
    """Sets the box_length and side_length variables for a cylindrical box based on the number of monomers"""
    global box_length
    global side_length
    global box_dimensions
    #cuboid box: square prism that houses the cylinder
    radius = diameter / 2
    box_length = radius * 10 + 2 * WALL # in LJ units; in terms of eq. bond length
    side_length = 2 * radius + 2 * WALL # diameter of the cylinder 

def SetCubicalBoxDimensions() -> None:
    """Sets the box_length and side_length variables for a big cubical box to house relaxed rings"""
    global box_length
    global side_length
    
    PADDING = 5
    radius = numberOfMonomers / 2 / math.pi + PADDING
    box_length = round(radius * 2) 
    side_length = box_length

def SetConstants() -> None:
    """Reads the command line arguments and sets the global variables"""
    global numberOfMonomers
    global architecture
    global numberOfBonds
    global offsetAngles
    global insertTerCrossLink
    global box_dimensions
    global numberOfRuns

    inputs = sys.argv
    numberOfMandatoryArguments = 2
    if len(inputs) < numberOfMandatoryArguments + 1:
        print("Not enough arguments passed from the command line! Please pass the number of monomers (integer) and the architecture (string).")
        print("Optionally, additional arguments for offset angles for each polymer can be passed.")
        exit(1)

    if inputs[1].isdigit():
        numberOfMonomers = int(inputs[1])
    else:
        print(f"The passed number of monomers ({inputs[1]}) could not be converted to an integer! Terminating.")
        exit(1)

    architecture = inputs[2]

    # reading optional arguments:
    if len(inputs) > numberOfMandatoryArguments + 1:
        for i in range(numberOfMandatoryArguments + 1, len(inputs)):
            try:
                offsetAngles.append(float(inputs[i]))
            except ValueError:
                print(f"The optionl argument {inputs[i]} could not be converted to a float! Terminating.")
                exit(1)
    else:
        offsetAngles = [0, 0] # both polymer's angles set to 0

    if numberOfPolymers == 1:
        insertTerCrossLink = False

    numberOfBonds = numberOfPolymers * (numberOfMonomers + len(crossLinks)) # Total in all polymers
    if insertTerCrossLink:
        numberOfBonds += 1

    SetCubicalBoxDimensions()
    if isBoxCentredAtOrigin:
        box_dimensions = [[-side_length/2, side_length/2], [-side_length/2, side_length/2], [-box_length/2, box_length/2]]
    else:
        box_dimensions = [[0,side_length], [0, side_length], [0, box_length]]

def GetRingCoordinates(radius: float, index: int, indexOfPolymer: int, offsetAngle: float = 0) -> Tuple[float, float, float]:
    """Returns the 3D coordinates of a monomer with the passed index place along a ring of the passed radius of a given polymer.
    An offset angle (in degrees) can be passed to rotate the ring about the y axis"""
    angle = index/numberOfMonomers * 2 * math.pi + math.radians(offsetAngle) # in radian
    # calculating the coordinates for the indexed monomer place on a ring of the given radius in the x-z plane
    y = indexOfPolymer * side_length / (numberOfPolymers + 1) - side_length/2
    x = radius * math.sin(angle)
    z = radius * math.cos(angle)
    if not isBoxCentredAtOrigin: # The box has a corner at the origin such that all coordinates are positive
        y = indexOfPolymer * side_length / (numberOfPolymers + 1)
        x += side_length / 2
        z += box_length / 2
    return (x, y, z)

def GetRelaxedRingCoordinates(index: int, indexOfPolymer: int, offsetAngle: float = 0) -> Tuple[float, float, float]:
    """Returns the 3D coordinates for one monomer with the passed index on a ring such that it is unit distance away from its neighbours
    An offset angle can be passed to rotate the ring of monomers about the y axis"""
    radius = numberOfMonomers / 2 / math.pi
    # verifying that the ring will fit in the box:
    if (2*radius > box_length) or (2*radius > side_length):
        print(f"A ring of radius {radius} will not fit in the specified box! Terminating\n")
        exit(1)
    
    return GetRingCoordinates(radius, index, indexOfPolymer, offsetAngle = offsetAngle)

def GetRelaxedSemiRingCoordinates(index: int, indexOfPolymer: int, offsetAngle: float = 0) -> Tuple[float, float, float]:
    """
    Returns the 3D coordinates for one monomer placed on a semicircle, meant for linear polymers. 
    The radius of this semicircle is the same as that of the relaxed ring used for ring polymers.
    Args:
        index: The index of the monomer (starting from zero)
        indexOfPolymer: The index of the polymer (starting from 1)
        offsetAngle: An angle in radians to rotate the whole semi circle with respect to the box. Default value = 0
    """
    radius = numberOfMonomers / 2 / math.pi
    # verifying that the ring will fit in the box:
    if (2*radius > box_length) or (2*radius > side_length):
        print(f"A semi circle of radius {radius} will not fit in the specified box! Terminating\n")
        exit(1)

    angle = index/numberOfMonomers * math.pi + math.radians(offsetAngle) # in radian
    # calculating the coordinates for the indexed monomer place on a ring of the given radius in the x-z plane
    y = indexOfPolymer * side_length / (numberOfPolymers + 1) - side_length/2
    x = radius * math.sin(angle)
    z = radius * math.cos(angle)
    if not isBoxCentredAtOrigin: # The box has a corner at the origin such that all coordinates are positive
        y = indexOfPolymer * side_length / (numberOfPolymers + 1)
        x += side_length / 2
        z += box_length / 2
    return (x, y, z)

def GetRandomCoordinates(index: int) -> Tuple[float, float, float]:
    """Returns random 3D coordinates in the cuboidal box set by box_dimensions"""
    coords = [0, 0, 0]
    for index in range(3):
        coords[index] = np.random.uniform(box_dimensions[index][0], box_dimensions[index][1])
    return tuple(coords)

# to get random coordinates within the cylinder
def GetRandomCoordinatesInCylinder() -> Tuple[float, float, float]:
    """Returns random 3D coordinates in a cylindrical box with dimensions given in box_dimensions"""
    coords = [0, 0, 0] # cylindrical coordinate system: r, phi, z
    coords[0] = np.random.uniform(0, side_length / 2) # the radial coordinate from 0 to max radius
    coords[1] = np.random.uniform(0, 2 * math.pi) # the polar angle
    coords[2] = np.random.uniform(box_dimensions[2][0], box_dimensions[2][1]) # the z coordinate

    # converting to cartesian coordinates:
    z = coords[2]
    x = coords[0] * math.cos(coords[1])
    y = coords[0] * math.sin(coords[1])
    if not isBoxCentredAtOrigin:
        x += side_length / 2 # offset the coordinates to lie within the bounds of the box
        y += side_length / 2

    return (x, y, z)

def GetCoordinatesAlongALine(index: int, indexOfPolymer: int) -> Tuple[float, float, float]:
    """Returns 3D coordinates for the passed polymer along a line parallel to the z axis"""
    # The y coord depends on index of polymer
    y = indexOfPolymer * side_length / (numberOfPolymers + 1) - side_length/2
    x = 0 # The line is in the x = 0 plane
    PADDING = 0.5
    z = (box_length - PADDING)/numberOfMonomers * index - box_length / 2
    if not isBoxCentredAtOrigin:
        x += side_length/2
        y += side_length/2
        z += box_length/2
    
    return (x, y, z)

    
def PrintHeader(file: io.TextIOWrapper) -> None:
    """Writes the header section for a LAMMPS input file to the file passed"""
    # initial comment:
    file.write(f"# input data file for {numberOfPolymers} {architecture} polymer(s) with {numberOfMonomers} monomers\n# Offset angles:")
    for i in range(numberOfPolymers): # writing the offset angles used for the polymers
        file.write(f" {offsetAngles[i]:.2f}")
    file.write("\n\n")
    # number of atoms:
    file.write(f'{numberOfPolymers * numberOfMonomers} atoms\n')
    file.write(f'{numberofAtomTypes} atom types\n')
    # number of bonds:
    file.write(f'{numberOfBonds} bonds\n')
    file.write('1 bond types\n\n')
    # dimensions of box:
    directions = ['x', 'y', 'z']
    for index in range(3):
        file.write(f'{box_dimensions[index][0]} {box_dimensions[index][1]} {directions[index]}lo {directions[index]}hi \n')

def PrintAtomsSection(file: io.TextIOWrapper, offsetAngles: Tuple[float] = (0, 0)) -> None:
    """Writes the atoms section for a LAMMPS input file to the file passed
    The atoms are place along a big ring with an offset angle from the z axis.
    There is an offset angle for each polymer; the angles must be passed as a tuple"""

    if len(offsetAngles) < numberOfPolymers:
        print("Not enough offset angles passed! Please pass one for each polymer. Terminating.")
        exit(1)
    # header:
    file.write("\nAtoms # bond\n\n")
    # content: atom-id, molecule-id, atom type, x, y, z
    for indexOfPolymer in range(1, numberOfPolymers+1):
        for index in range(numberOfMonomers):
            # computing positions:
            x, y, z = GetRelaxedRingCoordinates(index, indexOfPolymer, offsetAngles[indexOfPolymer - 1])
            atomID = (indexOfPolymer - 1) * numberOfMonomers + (index + 1)
            file.write("%i %i %i %.6lf %.6lf %.6lf \n" % (atomID, 0, indexOfPolymer, x, y, z)) # atom ID, moleculed ID set to 0, atom type, coords


def PrintCrossLinks(file: io.TextIOWrapper) -> None:
    """Writes the cross link (bonds) to the file as part of the bonds section of a LAMMPS input file"""
    counter = 1
    for indexOfPolymer in range(numberOfPolymers):
        for links in crossLinks:
            atom1_ID = indexOfPolymer * numberOfMonomers + links[0]
            atom2_ID = indexOfPolymer * numberOfMonomers + links[1]
            file.write(f"{numberOfPolymers * numberOfMonomers + counter} 1 {atom1_ID} {atom2_ID}\n") # bond id, bond type, atom1ID, atom2ID
            counter += 1

def InsertTerCrossLink(file: io.TextIOWrapper) -> None:
    """Writes the ter cross link bond to the passed file as part of the bonds section of a LAMMPS input file"""
    atom1_ID = numberOfMonomers
    atom2_ID = numberOfPolymers * numberOfMonomers
    file.write(f"{numberOfBonds} 1 {atom1_ID} {atom2_ID}") # bond id, bond type, atom1ID, atom2ID

def PrintBondsSection(file: io.TextIOWrapper) -> None:
    """Writes the bond section for a LAMMPS input file to the file passed"""
    #header:
    file.write("\nBonds\n\n")
    for indexOfPolymer in range(numberOfPolymers):
        for index  in range(indexOfPolymer * numberOfMonomers + 1, (indexOfPolymer + 1) * numberOfMonomers + 1):
            nextIndex = index + 1
            if index == (indexOfPolymer + 1) * numberOfMonomers: # last atom
                nextIndex = indexOfPolymer * numberOfMonomers + 1
            file.write(f"{index} 1 {index} {nextIndex}\n") # bond id, bond type, atom1ID, atom2ID
    PrintCrossLinks(file)
    if insertTerCrossLink:
        InsertTerCrossLink(file)

if __name__ == "__main__":
    SetConstants()
    filePath = f"initial_configuration.txt" # relative path
    filePrefix = filePath.split('.')[0]
    # SetCubicalBoxDimensions()
    # offsetAngles = (0, 90) # in degrees
    for i in range(1, numberOfRuns + 1):
        # Setting random seed
        np.random.seed(i)
        if i != 1:
            filePath = f"{filePrefix}{i}.txt"
        with open(filePath, "w") as file:
            PrintHeader(file)
            PrintAtomsSection(file, offsetAngles)
            PrintBondsSection(file)
    
