# Only required to get the example working without
# making this package pip installable
import sys
sys.path.append("..")

from DesiSequencer import readRawFile
from DesiSequencer.sequence import sequenceFromDesiFile


if __name__ == "__main__":
    # Before beginning, we must first check the utils.py script and change the monomers
    # to whatever we used. Once that is done, we can begin.

    # Set some information for sequencing
    endcap_mass=262.084
    endcap_name='Tyr(OMe)'

    # Convert rawfile 
    converted_file = readRawFile('./examples/test.RAW')

    # Sequence the converted file
    sequenceFromDesiFile(converted_file, endcap_mass = endcap_mass, endcap_name = endcap_name)

    