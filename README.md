# PeakPicker

# Workflow

    1. Use Visual Studio and Program.cs in the RawFileReaderExample directory to decode raw files into textfiles.

    2. Change the directory information around line 149 to the appropriate location of the raw files.

    3. Change the run and debug preferences to input the CLI arguments.

    4. Run the program. The output will be a text file with the scan data which can be processed, plottd, and used to determine sequence information.
    
    5. Use the remaining Python scripts of `DesiSequencer` to parse the converted data files.

    6. In `sequence.py`, scroll down the the "if __name__ == "__main__": portion of the script. You can enter the path to the converted datafile here.

    7. Run the sequence.py script.





