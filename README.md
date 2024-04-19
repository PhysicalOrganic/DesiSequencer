# PeakPicker

Converts raw DESI mass spectrometer data (.RAW) into plain text files which can be analyzed
by the Python code to extract sequence information.
## Workflow

    1. Use Visual Studio and Program.cs in the RawFileReaderExample directory to decode raw files
       into textfiles. Change the directory information around line 149 to the appropriate location of the raw files.

       The readRawFile.py script can also be used to extract the data from the .RAW
       file types. However, there are several instances in which the *ion count* or *intensity* of a
       signal was incorrectly read as 1 unit off (i.e., signal that was 187,004 became 187,003). The
       source of this difference is unknown.

    2. Change the run and debug preferences to input any CLI arguments. In the project properties, change working dir /Users/jameshoward/Documents/Programming/PeakPicker/DesiSequencer/data

    3. Run the program. The output will be a text file with the scan data which can be processed,
       plotted, and used to determine sequence information.

    4. Use the remaining Python scripts of `DesiSequencer` to parse the converted data files.

    5. In `sequence.py`, scroll down the the "if __name__ == "__main__": portion of the script.
       You can enter the path to the converted datafile here.

    6. Run the sequence.py script.

    Recommended: look at the example1.py file in DesiSequenceer/examples/


# Install
install this https://www.mono-project.com/download/stable/#download-mac
install this https://dotnet.microsoft.com/en-us/download/dotnet/8.0 (the runtime installer and the sdk)
https://dotnet.microsoft.com/en-us/download/dotnet/thank-you/runtime-8.0.4-macos-arm64-installer
https://visualstudio.microsoft.com/thank-you-downloading-visual-studio-mac/?sku=communitymac&rel=17

pip install pythonnet
pip install numpy
conda install cmake

git clone https://github.com/dotnet/runtime.git
cd


$ brew install pkg-config
$ brew cask install mono-mdk

export DOTNET_ROOT=/usr/local/share/dotnet/

pip install pandas
pip install matplotlib
pip install scikit-learn

install C# in vscode

pip install xlsxwriter