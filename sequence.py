import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

import numpy as np

from sklearn.neighbors import KernelDensity
from scipy.signal import argrelextrema

from utils import MonomerMasses
from plotting import plotTwoMassSpectra, plotMassSpectrum, plot_spots
from massCluster import findLargestMassPeak


def getArgs():
    pass

def get_spots(intensities: list[float],
              threshold: float,
              spot_width: int) -> list[int]:
    # identify the spots based on intensity
    IsolateSpots = []
    TempSpots = []
    #for i, intensity in enumerate(intensities): # Go through all avg intensities
    #    if intensity > threshold and len(TempSpots) < spot_width: # If it is above the threshold and #TODO less than the 9th scan intensity?
    #        TempSpots.append(i) # Add it to the temporary spots
    #    elif len(TempSpots) >= spot_width: # If the count is above or equal to the spot width (6 currently)
    #        IsolateSpots.append(TempSpots) # Add this to the isolate spots
    #        TempSpots = [] # Reset the temp spots?
    count = 0
    for i, intensity in enumerate(intensities): # Go through all avg intensities
        if intensity > threshold and count < 9: # If it is above the threshold and #TODO less than the 9th scan intensity?
            count += 1
            TempSpots.append(i) # Add it to the temporary spots
        elif count >= spot_width: # If the count is above or equal to the spot width (6 currently)
            IsolateSpots.append(TempSpots) # Add this to the isolate spots
            TempSpots = [] # Reset the temp spots?
            count = 0


    #for i, intensity in enumerate(AverageScanIntensities): # Go through all avg intensities
    #    if intensity > IntensityThreshold and count < 9: # If it is above the threshold and #TODO less than the 9th scan intensity?
    #        count += 1
    #        TempSpots.append(i) # Add it to the temporary spots
    #    elif count >= SpotWidth: # If the count is above or equal to the spot width (6 currently)
    #        IsolateSpots.append(TempSpots) # Add this to the isolate spots
    #        TempSpots = [] # Reset the temp spots?
    #        count = 0

    return IsolateSpots

def CleanUpScan(
        scan: str,
        mz_round: int = 2) -> dict:
    '''
    Takes the raw scan data formatted as a string and splits it
    into a dictionary of m/z:intensity pairs.

    Parameters
    ----------
    scan: str
        Original scan data formatted like:
        ...m/z, intensity;m/z, intensity;m/z, intensity;...

    mz_round: int
        Number of decimal points to round each m/z

    Returns
    ----------
        dict
    '''

    NewDict = {}

    scan = scan.split(';')
    for pair in scan:
        if len(pair) > 0:
            pair = pair.strip(' ')
            pair = pair.split(',')
            NewDict[round(float(pair[0]),2)] = round(float(pair[1]),1)

    return NewDict

def getAverageIntensity(
        scan: dict,
        mass_thresh: float = 500) -> float:
    '''
    Gets the average intensity of all the masses of a scan
    which is formatted as a dictionary of mass:intensity key:value
    pairs. Returns the average intensity.

    scan: dict
        Dictionary of mass:intensity pairs which represents
        a mass spectrum scan

    mass_thresh: float
        Threshold overwhich masses are considered
    '''
    return np.average([scan[x] for x in scan.keys() if x >= mass_thresh])

def findPeaks(
    scanData: pd.DataFrame,
    parent_mass: float,
    intensity_threshold: float = 0.10,
    minimum_mass: int = 250,
    debug = True): #110000
    '''
    Identifies the signals associated with the loss of a monomer
    from the parent oligomer. The parent oligomer is identified
    by the mass

    Parameters
    ----------
    ScanData: pd.DataFrame
        Dataframe containing two columns: Mass, Intensity. Typically the averaged
        mass spectrum of a particular sample (i.e., the average of multiple scans
        which are associated with a spot on a DESI plate)

    parent_mass: float
        Mass of the parent oligomer of length n which is used to find the
        n - x oligomers resulting from sequencing. x represents the sequence
        iteration.

    intensity_threshold: float (default = 0.05)
        The percentage of intensity of the largest intensity signal
        within the scanData["Intensity"] column above which signals
        are considered for selection.

        Default is 5% of maximum intensity signal

    minimum_mass: int (default = )
        The mass below which signals are considered insignificant.

    Returns
    ----------
    monomers: dict
        Mass:intensity pairs of the selected "peaks" which will be used
        to subtract and thus get the mass loss. This mass loss will eventually
        be matched to a monomer identity.
    '''

    # Convert intensity threshold to count
    intensity_percent = intensity_threshold
    intensity_threshold = intensity_threshold * scanData['Intensity'].max()

    # Get the largest mass in the spectrum
    # This is the parent peak of the entire spectrum
    LargestMassPeakIntensity = float(scanData[scanData['Mass'] == parent_mass]['Intensity'].iloc[0])

    if debug:
        print(f'[DEBUG] In findPeaks: intensity_percent = {intensity_percent}')
        print(f'[DEBUG] In findPeaks: intensity_threshold = {intensity_threshold}')
        print(f'[DEBUG] In findPeaks: LargestMassPeakIntensity = {LargestMassPeakIntensity} counts.')

    # dictionary containing identity of monomers and the mass at which they are found
    monomers = {}

    # First "monomer" is not actually a monomer, it is the mass:intensity pair
    # of the oligo which will lose a monomer once sequenced.
    # i.e., this is the parent peak
    monomers[parent_mass] = LargestMassPeakIntensity

    # Get the lighest and heaviest monomers
    heaviest_monomer = max(MonomerMasses.keys())
    lightest_monomer = min(MonomerMasses.keys())

    current_peak = parent_mass
    while True:

        # Get the window of signals that we are looking at in tmp_df
        lower_bound = current_peak - (heaviest_monomer * 1.05)
        upper_bound = current_peak - (lightest_monomer * 0.95)

        # Primary exit condition
        if upper_bound <= minimum_mass:
            if debug:
                print(f'[DEBUG] Minimum mass ({minimum_mass}) was surpassed with lower_bound of {lower_bound}')
            break

        tmp_df = scanData[(scanData['Mass'] >= lower_bound) & (scanData['Mass'] <= upper_bound)]

        # Filter out low intensity signals in the window
        #print(f'WARNING: for debugging: only considering masses within 50% of max signal in window')
        tmp_df = tmp_df[tmp_df['Intensity'] >= (tmp_df['Intensity'].max() * intensity_percent)] # Change 0.5 here to intensity_percent to return to normal
        #print(tmp_df.shape, ScanData.shape, 'lower bound',current_peak - heaviest_monomer, 'upper bound', current_peak - lightest_monomer)

        # Add the range to the figure
        if len(monomers) % 2 == 0:
            color = 'purple'
        else:
            color = 'green'

        h = scanData[(scanData['Mass'] >= lower_bound) & (scanData['Mass'] <= upper_bound)].Intensity.max() * 1.05
        rect = plt.Rectangle((lower_bound, 0), width = upper_bound - lower_bound, height = h, color=color, alpha = 0.3)
        plt.gca().add_patch(rect)

        # If the df is empty, stop sequencing
        if tmp_df.empty:
            print(f'tmp_df is empty')
            break

        # Go through each of the masses and monomers
        for mass, monomer_identity in MonomerMasses.items():
            # The scores are the difference between
            tmp_df[f'delta_{monomer_identity}'] = abs((current_peak  - tmp_df['Mass']) - mass)

        #print(f'Minimum values along axes 0: \n{tmp_df.min(axis=0)}')
        #print(f'Minimum values along axes 1: \n{tmp_df.min(axis=1)}')
        # Get the minimum of each row
        tmp_df['min'] = tmp_df.min(axis=1)
        tmp_df['score'] = tmp_df['min'] * (1 / (tmp_df['Intensity'] * 2)) #
        tmp_df['score'] = tmp_df['score'] * (1 / tmp_df['score'].max()) # Normalize the score to 1

        if debug:
            print(f'[DEBUG] Scoring DataFrame for monomer {len(monomers)}')
            print(f'[DEBUG] Selecting loss of monomer from previous mass (i.e., current_peak) {current_peak}')

        # Select the column that has th lowest delta value # OLD, didn't work because there were some that were close
        #TODO Warn users if there are similar or other deltas close to the absolute minimum
        #global_min_column = tmp_df[[x for x in tmp_df.columns if 'delta' in x]].min().idxmin()
        #current_peak = float(tmp_df.loc[tmp_df[global_min_column] == tmp_df[global_min_column].min(), 'Mass'].iloc[0])

        # Get the monomer identity? Might not be necessary
        global_min_column = tmp_df[[x for x in tmp_df.columns if 'delta' in x]].min().idxmin()
        # Reset current peak
        current_peak = float(tmp_df.loc[tmp_df['score'] == tmp_df['score'].min(), 'Mass'].iloc[0])
        monomers[current_peak] = float(tmp_df.loc[tmp_df['Mass'] == current_peak, 'Intensity'].iloc[0])

        if debug:
            pass
            #print(tmp_df[[x for x in tmp_df.columns if 'Ser' in x or 'Phe' in x or 'Leu' in x or 'Cha' in x or 'Tyr' in x or x in ['Mass', 'Intensity', 'min','score']]])
            #print(f'Global min column: {MatchMasstoMonomer()}')
            #TODO Match a single mass to a monomer

    if debug:
        print(f'[DEBUG] In findPeaks: Selected masses: {list(monomers.keys())}')

    return monomers

def MatchMasstoMonomer(
        MonomerMasses: dict,
        SpectraMasses,
        monomer_tolerance: float = 1,
        endcap_mass: float = 262.084,
        endcap_name: str = 'Tyr(OMe)',
        debug: bool = False):
    '''
    Function which assigns a monomer identity to a given mass

    Monomer masses was a dictionary of mass:intensity pairs
    '''
    # MonomerMasses has the keys as the mass

    # Is there a reason Mass:intensity pairs are shipped together into this function?
    masses_identified_as_oligos = list(SpectraMasses.keys())

    # List of monomer names that were identified
    monomer_names = []

    for x in range(len(masses_identified_as_oligos)-1):

        # Test for endcap + last monomer
        # Difference in mass between the current oligomer
        # and the endcap (which would be the last monomer)
        endcap_difference = masses_identified_as_oligos[x] - endcap_mass
        # Differences between the endcap difference and each monomer.
        difs = {monomer_symbol:k - endcap_difference for k, monomer_symbol in MonomerMasses.items()}
        # Test if any of these are less than the tolerance
        for monomer_symbol, m in difs.items():
            if abs(m) <= monomer_tolerance:
                monomer_names.append(monomer_symbol)
                monomer_names.append(endcap_name)
                return monomer_names

        mass_loss = masses_identified_as_oligos[x] - masses_identified_as_oligos[x+1]

        # Iterate over all the possible monomers and their masses
        # Dictionary for storing the difference between the mass loss and a potential monomer.
        # The lowest of this dictionary's values is the most likely monomer
        differences = {}
        for mass, monomer_symbol in MonomerMasses.items():
            differences[monomer_symbol] = abs(round(mass - mass_loss, 2))

        if debug:
            from pprint import pprint
            print('Differences between mass loss and possible monomers')
            pprint(differences)
            print()

        lowest_match_difference = min(differences.values())
        match = list(differences.keys())[list(differences.values()).index(min(differences.values()))]

        if lowest_match_difference >= 1:
            print(f'WARNING: Mass difference for {masses_identified_as_oligos[x]} - {masses_identified_as_oligos[x+1]} = {round(mass_loss, 2)}')
            print(f'was larger than the acceptable tolerance for monomer identification\nTolerance: {monomer_tolerance}\tDiff: {lowest_match_difference}\n')

        monomer_names.append(match)

    return monomer_names

def sequenceFromDesiFile(
        file: Path,
        endcap_mass=262.084,
        endcap_name = None,
        debug: bool = False
        ):
    '''
    Takes a converted file (.RAW to .txt) and
    sequences it.
    '''

    if endcap_name is None:
        endcap_name = 'END CAP (UNSPECIFIED)'

    name = str(file.stem)

    if debug:
        print(f"Processing File: {name}")

    with open(file, 'r', encoding='utf-8') as f:
        raw_data = f.read()

    # split combined scans into individual scan events
    raw_data = raw_data.split("|")

    #### CLEANING UP RAW FILE DATA ####
    # cleaning up scans to remove brackets and commas and whatnot and format them into a dictionary
    CombinedScans = []
    for Scan in raw_data:
        NewScan = CleanUpScan(Scan)
        if (len(NewScan) != 0):
            CombinedScans.append(NewScan)

    if debug:
        print("Number of Scans in the entire data file:", len(CombinedScans))

    # Get the average scan intensities of each scan to identify the spots
    AverageScanIntensities = []

    # plot the average intensity of each scan to visualize spots
    for i in CombinedScans:
        AverageScanIntensities.append(getAverageIntensity(i, mass_thresh=500))

    # this is the scan intensity that you check to see what is noise and what is oligos
    IntensityThreshold = 10000

    # this is how many scans should make up one spot to get rid of the spikes there sometimes are
    SpotWidth = 6
    IsolateSpots = get_spots(intensities=AverageScanIntensities,
                             threshold=IntensityThreshold,
                             spot_width=SpotWidth)

    #print(f'WARNING: Debugging the scan with manual spot picking.')
    #IsolateSpots = [[0, 1, 2, 3, 4, 5], [11, 12, 13, 14, 15, 16, 17], [23, 24, 25, 26, 27, 28, 29, 30], [37, 38, 39, 40, 41, 42, 43, 44, 45], [49, 50, 51, 52, 53, 54, 55, 56, 57], [62, 63, 64, 65, 66, 67, 68]]

    if len(IsolateSpots) != 1:
        print(f'Found {len(IsolateSpots)} different spots in the data file.\n')

    if debug:
        print(f"Spots found at scan Numbers: {IsolateSpots}\tNumber of scans: {len(CombinedScans)}\n")

    # List of dictionaries which have mass:intensity pairs
    # Is the length of the number of spots
    CombinedSpotScans = []

    # Iterate through the groups of
    # For every list of scans (dictionaries of mass:intensity pairs), these are combined because they represent a single sample or spot
    for i, scan_group in enumerate(IsolateSpots):

        # This gets the mass scan of the first peak (i[0] is the
        # first index of the peak)
        BasePeak = CombinedScans[scan_group[0]]

       # Each peak has multiple indices over which it exists. Iterate over them
       # For every scan that corresponds to a peak in the total ion count spectrum
        for j in scan_group:

            # For every mass: intensity pair in that scan
            for k in CombinedScans[j].keys():

                # If the mass is already in the "BasePeak", increase it's intensity
                if k in BasePeak.keys():
                    BasePeak[k] += CombinedScans[j][k] # I think we are adding the intensity of those masses from each scan to a single
                else:
                    BasePeak[k] = CombinedScans[j][k]
        CombinedSpotScans.append(BasePeak)

    if debug:
        plot_spots(average_scan_intensities=AverageScanIntensities,
                   isolated_spots=IsolateSpots,
                   save=Path(file.parent / f'{file.stem}_identified_spots.png'))

    # This contains all of the dataframes for each spot
    # which we will use to output the results later
    #raise NotImplementedError('This code above should be implemented.')
    dataframes_for_spots = []
    ms_spectra_for_spots = []

    # Iterate over every spot (which is the combined intensities over all scans which were done for that spot)
    for i, x in enumerate(CombinedSpotScans):

        # Initiate plotting here so we can plot within other functions
        fig, ax = plt.subplots(1,1, figsize = (12,6))

        ScansTemp = pd.DataFrame(list(x.items()), columns=['Mass', 'Intensity'])

        # Find the parent oligomer mass
        parent_mass = findLargestMassPeak(scanData=ScansTemp, threshold = 0.07, debug=debug)
        parent_mass_intensity = float(ScansTemp[ScansTemp["Mass"] == parent_mass]['Intensity'].iloc[0])

        #if i == 1:
        #    exit('DEBUG EXITING')
        print(f'SPOT {i}: Parent mass in spectrum: {parent_mass}\tIntensity: {parent_mass_intensity}\n')

        # Find the peaks
        FoundSpectraMassesIntensity = findPeaks(scanData=ScansTemp, parent_mass=parent_mass, debug=debug)

        # Convert the mass:intensity pairs into list of masses
        FoundSpectraMasses = list(FoundSpectraMassesIntensity.keys())

        # Match the list of masses to a particular monomer loss
        MassMatches = MatchMasstoMonomer(MonomerMasses, FoundSpectraMassesIntensity, endcap_mass=endcap_mass, endcap_name=endcap_name)

        print(f'[DEBUG] MassMatches=')
        print(MassMatches)

        # This is just stored for printing
        results = []
        # Iterate over the masses
        for j, mass_found_in_spectrum in enumerate(FoundSpectraMasses):
            print(f'[DEBUG] j={j} mass_found_in_spectrum={mass_found_in_spectrum}')
            if j == len(MassMatches):
                print(f'[WARNING] There were {j} or more masses identified in spectrum while we have {len(MassMatches)} mass matches')
                break
            if MassMatches[j] == endcap_name:
                results.append([j, endcap_mass, endcap_name + '(endcap)'])
            else:
                results.append([j, mass_found_in_spectrum, MassMatches[j]])

        result = pd.DataFrame(results, columns=['Signal No.', 'M/Z', 'Monomer lost']).set_index('Signal No.')

        # Append the result dataframe
        dataframes_for_spots.append(result)

        # Plot the full spectrum
        x, y = zip(*sorted(x.items()))
        ax.stem(x, y, linefmt='black', markerfmt='')
        ax.set_title(f"Spot {i + 1} of {len(IsolateSpots)} ({file.stem})")
        ax.set_ylabel('Counts')
        ax.set_xlabel(r'$m/z$')
        ax.stem(FoundSpectraMassesIntensity.keys(), FoundSpectraMassesIntensity.values(), linefmt='red', markerfmt='') #, width=1.5

        ax.text(parent_mass*1.005, parent_mass_intensity * 1.005, f'Parent peak {parent_mass}', color='red')

        tmp_image_path = Path(file.parent / str(name + f'_spot_{i + 1}' + ".png"))
        plt.savefig(tmp_image_path, format='png', dpi = 600)
        ms_spectra_for_spots.append(tmp_image_path)

    # Write results to the excel file
    savePath = Path(file.parent / str(name + ".xlsx"))
    writer = pd.ExcelWriter(savePath)
    wb = writer.book
    for i, d in enumerate(dataframes_for_spots):
        sheet_name = f"Spot {i}"
        d.to_excel(writer, sheet_name=sheet_name)

        # get the worksheet object
        ws = writer.sheets[sheet_name]
        ws.insert_image('E1', ms_spectra_for_spots[i])

    writer.close()

    for p in ms_spectra_for_spots:
        p.unlink()

    print(f'\nFile saved at: {savePath.absolute()}')

if __name__ == "__main__":
    # Define a file to analyze
    f = Path('./DesiSequencer/data/Oligomers March 2023/manualProfiling_mayaAngelou_spot1.txt')

    # Set some parameters
    endcap_mass=262.084
    endcap_name='Tyr(OMe)'
    debug = True

    sequenceFromDesiFile(f, debug=debug, endcap_mass=262.084, endcap_name='Tyr(OMe)')



