
import numpy as np

from pathlib import Path

import clr

# Add references to the DLLs
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parents[0]))
clr.AddReference('ThermoFisher.CommonCore.Data')
clr.AddReference('ThermoFisher.CommonCore.RawFileReader')
clr.AddReference('ThermoFisher.CommonCore.BackgroundSubtraction')
clr.AddReference('ThermoFisher.CommonCore.MassPrecisionEstimator')
clr.AddReference('System.Collections')

from ThermoFisher.CommonCore.RawFileReader import RawFileReaderAdapter
from ThermoFisher.CommonCore.Data.Business import Device

# imports I added
from ThermoFisher.CommonCore.Data.Interfaces import IRawDataExtended

def getSpectrum(
        rawFile: IRawDataExtended, 
        scanNumber: int, 
        scanFilter = None, 
        outputData: bool = False, 
        verbose: bool = False):
    '''Gets the spectrum from the RAW file. Modified by James.

    Args:
        rawFile (IRawDataPlus): the RAW file being read.
        scanNumber (int): the scan number being read.
        scanFilter (str): the scan filter for that scan. Get with
            rawFile.GetFilterForScanNumber(scanNumber)
        outputData (bool): the output data flag.
    '''

    # Check for a valid scan filter
    if scanFilter is None:
        scanFilter = rawFile.GetFilterForScanNumber(scanNumber)

    # Get the scan statistics from the RAW file for this scan number
    scanStatistics = rawFile.GetScanStatsForScanNumber(scanNumber)

    # Check to see if the scan has centroid data or profile data.  Depending upon the
    # type of data, different methods will be used to read the data.  While the ReadAllSpectra
    # method demonstrates reading the data using the Scan.FromFile method, generating the
    # Scan object takes more time and memory to do, so that method isn't optimum.
    
    if scanStatistics.IsCentroidScan:
        # Get the centroid (label) data from the RAW file for this
        # scan
        centroidStream = rawFile.GetCentroidStream(scanNumber, False)

        if verbose:
            print('Spectrum (centroid/label) {} - {} points'.format(scanNumber, centroidStream.Length))

        # Print the spectral data (mass, intensity, charge values).
        # Not all of the information in the high resolution centroid
        # (label data) object is reported in this example.  Please
        # check the documentation for more information about what is
        # available in high resolution centroid (label) data.
        if outputData:
            for i in range(centroidStream.Length):
                print('  {} - {:.4f}, {:.0f}, {:.0f}'.format(
                    i, centroidStream.Masses[i], centroidStream.Intensities[i], centroidStream.Charges[i]))
            print()

    else:
        # Get the segmented (low res and profile) scan data
        segmentedScan = rawFile.GetSegmentedScanFromScanNumber(scanNumber, scanStatistics)

        if verbose:
            print('Spectrum (normal data) {} - {} points'.format(scanNumber, segmentedScan.Positions.Length))

        # Print the spectral data (mass, intensity values)
        if outputData:
            for i in range(segmentedScan.Positions.Length):
                print('  {} - {:.4f}, {:.0f}'.format(
                    i, segmentedScan.Positions[i], segmentedScan.Intensities[i]))
            print()

def getSpectrum(
        rawFile: IRawDataExtended, 
        scanNumber: int, 
        scanFilter = None, 
        verbose: bool = False):
    '''Gets the spectrum from the RAW file. Modified by James.

    Args:
        rawFile (IRawDataPlus): the RAW file being read.
        scanNumber (int): the scan number being read.
        scanFilter (str): the scan filter for that scan. Get with
            rawFile.GetFilterForScanNumber(scanNumber)
        outputData (bool): the output data flag.
    '''

    # Check for a valid scan filter
    if scanFilter is None:
        scanFilter = rawFile.GetFilterForScanNumber(scanNumber)

    # Get the scan statistics from the RAW file for this scan number
    scanStatistics = rawFile.GetScanStatsForScanNumber(scanNumber)

    # Check to see if the scan has centroid data or profile data.  Depending upon the
    # type of data, different methods will be used to read the data.  While the ReadAllSpectra
    # method demonstrates reading the data using the Scan.FromFile method, generating the
    # Scan object takes more time and memory to do, so that method isn't optimum.
    
    if scanStatistics.IsCentroidScan:
        # Get the centroid (label) data from the RAW file for this
        # scan
        centroidStream = rawFile.GetCentroidStream(scanNumber, False)

        if verbose:
            print('Spectrum (centroid/label) {} - {} points'.format(scanNumber, centroidStream.Length))

        # Print the spectral data (mass, intensity, charge values).
        # Not all of the information in the high resolution centroid
        # (label data) object is reported in this example.  Please
        # check the documentation for more information about what is
        # available in high resolution centroid (label) data.
        # list of position, intensity pairs
        pairs = []
        for i in range(centroidStream.Length):
            pairs.append([centroidStream.Masses[i], centroidStream.Intensities[i]])
            #print('  {} - {:.4f}, {:.0f}, {:.0f}'.format(
            #    i, centroidStream.Masses[i], centroidStream.Intensities[i], centroidStream.Charges[i]))
        return np.array(pairs)

    else:
        # Get the segmented (low res and profile) scan data
        segmentedScan = rawFile.GetSegmentedScanFromScanNumber(scanNumber, scanStatistics)

        if verbose:
            print('Spectrum (normal data) {} - {} points'.format(scanNumber, segmentedScan.Positions.Length))

        # Print the spectral data (mass, intensity values)
        for i in range(segmentedScan.Positions.Length):
            print('  {} - {:.4f}, {:.0f}'.format(
                i, segmentedScan.Positions[i], segmentedScan.Intensities[i]))
            print()

def readRawFile(filename: str, 
                outPath: Path = None) -> Path:
    '''
    Python wrapper for writing spectral data from the DESI
    to a file path.

    Parameters
    ----------
    filename: str
        Relative directory path to the raw file.

    outPath: Path
        Path which the data will be written to. Defaults to
        the same directory and name with the extension '.out'

    Returns
    ----------
    outPath: Path
        This is just returned as a convenience in your script. 

    '''

    # Create the IRawDataPlus object for accessing the RAW file
    rawFile = RawFileReaderAdapter.FileFactory(filename)

    # Get the number of instruments (controllers) present in the RAW file
    # and set the selected instrument to the MS instrument, first instance
    # of it. We must select the instrument before anything else.
    rawFile.SelectInstrument(Device.MS, 1)

    if outPath is None:
        filename = Path(filename)
        name = filename.stem
        outPath = Path(filename.parent / str(name + '.out'))

    n_scans = rawFile.RunHeaderEx.SpectraCount

    with open(outPath, 'w') as o:
        for i in range(n_scans):
            s = getSpectrum(rawFile, i + 1, rawFile.GetFilterForScanNumber(i + 1), False)
            for pair in s:
                mass = round(pair[0], 4)
                if len(str(mass).split('.')[1]) == 4:
                    pass
                elif len(str(mass).split('.')[1]) == 3:
                    mass = str(mass) + '0'
                elif len(str(mass).split('.')[1]) == 2:
                    mass = str(mass) + '00'
                elif len(str(mass).split('.')[1]) == 1:
                    mass = str(mass) + '000'
                else:
                    raise Exception(mass)
                o.write(f'{mass}, {int(round(pair[1], 0))};')
            o.write('|')

    return outPath