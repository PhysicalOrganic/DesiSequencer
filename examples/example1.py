# Only required to get the example working without
# making this package pip installable
#TODO Make package pip installable
import sys
sys.path.append("..")

from pathlib import Path

from DesiSequencer import readRawFile
from DesiSequencer.sequence import sequenceFromDesiFile


if __name__ == "__main__":
    # Before beginning, we must first check the utils.py script and change the monomers
    # to whatever we used. Once that is done, we can begin.

    # Set endcap information and debug settings
    endcap_mass = 262.084
    endcap_name = 'Tyr(OMe)'
    debug = True

    # Convert rawfile, this is done using the python version, but the resulting
    # file can be made using the C# code as well (stated in README)
    #converted_file = readRawFile('./examples/test.RAW')
    converted_file = Path('/Users/admin-js99925/Documents/ANSLYN GROUP/DESI april 2024/04172024_oligomers_CVprofiling_A-F_1sttry.txt')

    # Sequence the converted file
    sequenceFromDesiFile(converted_file,
                         endcap_mass = endcap_mass,
                         endcap_name = endcap_name,
                         debug=debug)

