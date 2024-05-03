

from pathlib import Path

from scipy.signal import argrelextrema, savgol_filter, find_peaks

import matplotlib.pyplot as plt

from sklearn.neighbors import KernelDensity

import numpy as np
import pandas as pd

def read_desi_data_file(file: Path,
                        mz_round: int = 2) -> dict:
    '''
    Reads a converted (.RAW to .txt) desi file
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

    with open(file, 'r', encoding='utf-8') as f:
        raw_data = f.read()

    # Split combined scans into individual scan events
    # This is a list of strings
    raw_data = raw_data.split("|")

    cleaned_scans = []

    # Go through each scan
    for scan in raw_data:

        mass_intensity_pairs = {}

        # Get the mass:intensity pairs
        scan = scan.split(';')

        if scan == ['']:
            continue

        for pair in scan:
            if len(pair) > 0:
                pair = pair.strip(' ')
                pair = pair.split(',')
                mass_intensity_pairs[round(float(pair[0]), mz_round)] = round(float(pair[1]),1)

        cleaned_scans.append(mass_intensity_pairs)

    return cleaned_scans


def get_average_intensity_of_scan(scan: dict,
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

def get_spots(intensities: list[float],
              threshold: float,
              spot_width: int,
              look_ahead_bias: float = 0.8,
              polyorder: int = 2,
              plot_smoothed: bool = False,
              debug: bool = False) -> list[int]:

    if look_ahead_bias <= 0 or look_ahead_bias > 1:
        raise ValueError(f'Looking ahead must be greater than 0 and less than or equal to 1, not {look_ahead_bias}.')

    def movingaverage(interval, window_size):
        window = np.ones(int(window_size))/float(window_size)
        return np.convolve(interval, window, 'same')

    # Smooth the spectrum scan
    smoothed = savgol_filter(intensities, window_length=spot_width, polyorder=polyorder)
    smoothed = movingaverage(smoothed, window_size=spot_width)

    if threshold == 'auto' and not isinstance(threshold, float):
        threshold = np.average(smoothed) * 0.8
        print(f'Setting threshold to {threshold} counts')

    # Find the extreme
    extrema = list(find_peaks(smoothed)[0])

    print(f'[INFO] Found {len(extrema)} extrema in the spot scan. This means there are six peaks/samples.')

    # Make a list for holding the spots which are the extrema plus surrounding scans if they exceed threshold
    spots = []

    for extreme_point in extrema:
        scans_associated_with_spot = []
        scans_associated_with_spot.append(extreme_point)

        # Look ahead
        look_ahead_list = list([idx for idx, _ in enumerate(smoothed) if idx > extreme_point])
        #print(f'extrema: {extreme_point} looking head at {look_ahead_list}')
        for idx in look_ahead_list:
            if idx <= extreme_point:
                continue
            if smoothed[idx] < threshold * look_ahead_bias: # Look ahead since the trailing tails might be better?
                break

            #print(f'Scan {idx} is greater than {extreme_point} and has intensity greater than threshold')
            scans_associated_with_spot.append(idx)

        # Look behind
        look_behind_list = list([idx for idx, _ in enumerate(smoothed) if idx < extreme_point])
        look_behind_list.reverse()
        #print(f'extrema: {extreme_point} looking behind at {look_behind_list}')
        for idx in look_behind_list:
            if idx >= extreme_point or idx in scans_associated_with_spot:
                continue
            if smoothed[idx] < threshold:
                break
            #print(f'Scan {idx} is less than {extreme_point} and has intensity greater than threshold')
            scans_associated_with_spot.append(idx)

        spots.append(sorted(list(set(scans_associated_with_spot))))

    fig, ax = plt.subplots(1,1, figsize = (8,6))
    ax.plot(np.arange(0, len(intensities)), intensities, color='black')
    if plot_smoothed:
        ax.plot(np.arange(0, len(smoothed)), smoothed, color='purple', label='smoothed')
        ax.scatter(extrema, [smoothed[z] for z in extrema], marker='x', color='red')
    else:
        ax.scatter(extrema, [intensities[z] for z in extrema], marker='x', color='red')
    ax.set_xlabel('scan number')
    ax.set_ylabel('average scan intensity')
    ax.hlines(threshold, xmin=0, xmax=len(intensities), color='red', label='intensity threshold')
    ax.legend()
    for extreme_point, spot_group in zip(extrema, spots):
        if debug:
            print(f'Extreme point {extreme_point} has scans {spot_group}')
        ax.plot(spot_group, [intensities[z] for z in spot_group])
    plt.show()

    return spots

def combine_scans(scans: list[dict], debug: bool = True) -> dict:
    '''
    Combines scans together
    '''
    # This gets the mass scan of the first peak (i[0] is the first index of the peak)
    base_peak = scans[0]

    # Each peak has multiple indices over which it exists. Iterate over them
    # For every scan that corresponds to a peak in the total ion count spectrum
    for scan in scans:

        # For every mass: intensity pair in that scan
        for mass, intensity in scan.items():

            # If the mass is already in the "BasePeak", increase it's intensity
            if mass in base_peak.keys():
                base_peak[mass] += intensity
            else:
                base_peak[mass] = intensity

    return base_peak


def mass_to_monomer(mass: int,
                    mass_monomer_definitions: dict,
                    monomer_tolerance: float = 1,
                    endcap_mass: float = 262.084,
                    endcap_name: str = 'Tyr(OMe)',
                    debug: bool = False):
    '''
    Function which assigns a monomer identity to a given mass

    Monomer masses was a dictionary of mass:intensity pairs
    '''
    # MonomerMasses has the keys as the mass

    # Get the lighest and heaviest monomers
    heaviest_monomer = max(mass_monomer_definitions.keys())
    lightest_monomer = min(mass_monomer_definitions.keys())

    endcap_difference = mass - endcap_mass

    # Get the window of signals that we are looking at in tmp_df
    lower_bound = mass - (heaviest_monomer * 1.05)
    upper_bound = mass - (lightest_monomer * 0.95)

    #print(f'Bounds for mass: {mass}: {lower_bound} - {upper_bound}: endcap dif: {endcap_difference}')

    # Return the endcap if found
    if abs(endcap_difference) <= 0.1 * endcap_mass:
        if endcap_difference >= monomer_tolerance:
            print(f'[WARNING] In range of endcap but found {mass} instead of {endcap_mass}')
        return endcap_name

    difs = {monomer_symbol:abs(monomer_mass_defn - mass) for monomer_mass_defn, monomer_symbol in mass_monomer_definitions.items()}
    df = pd.DataFrame([difs])
    #display(df)
    return df.min().idxmin()


def identify_masses_from_sequencing(scan_data: pd.DataFrame,
                                    mass_monomer_definitions: dict,
                                    parent_mass: float,
                                    intensity_threshold: float = 0.10,
                                    endcap_mass: float = 262.084,
                                    debug = True) -> tuple[dict, list[tuple]]: #110000
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

    bounds: list[tuple]
        List of tuples of the lower and upper bounds for the windows of M/Z
        values within which the peaks were identified (useful for plotting)
    '''

    # Convert intensity threshold to count
    intensity_percent = intensity_threshold
    intensity_threshold = intensity_threshold * scan_data['Intensity'].max()

    # Get the largest mass in the spectrum
    # This is the parent peak of the entire spectrum
    LargestMassPeakIntensity = float(scan_data[scan_data['Mass'] == parent_mass]['Intensity'].iloc[0])

    if debug:
        print(f'[DEBUG] intensity_percent: {intensity_percent}\tfirst intensity_threshold: {intensity_threshold}\tLargestMassPeakIntensity: {LargestMassPeakIntensity} counts.')

    # dictionary containing identity of monomers and the mass at which they are found
    monomers = {}

    # First "monomer" is not actually a monomer, it is the mass:intensity pair
    # of the oligo which will lose a monomer once sequenced (i.e., this is the parent peak).
    monomers[parent_mass] = LargestMassPeakIntensity
    bounds = []

    # Get the lighest and heaviest monomers
    heaviest_monomer = max(mass_monomer_definitions.keys())
    lightest_monomer = min(mass_monomer_definitions.keys())

    current_peak = parent_mass
    while True:

        # Get the window of signals that we are looking at in tmp_df
        lower_bound = current_peak - (heaviest_monomer * 1.02)
        upper_bound = current_peak - (lightest_monomer * 0.98)
        bounds.append((lower_bound, upper_bound))

        # Primary exit condition
        if upper_bound <= endcap_mass:
            if debug:
                print(f'[DEBUG] Upperbound mass {upper_bound} was lower than the minimum mass {endcap_mass}')
            break

        tmp_df = scan_data[(scan_data['Mass'] >= lower_bound) & (scan_data['Mass'] <= upper_bound)].copy(deep=True)

        # Filter out low intensity signals in the window
        #tmp_df = tmp_df[tmp_df['Intensity'] >= (tmp_df['Intensity'].max() * intensity_percent)] # Change 0.5 here to intensity_percent to return to normal

        # If the df is empty, stop sequencing
        if tmp_df.empty:
            print(f'[WARNING] Sequencing terminated because range of signals was empty!')
            break

        # Go through each of the masses and monomers
        for mass, monomer_identity in mass_monomer_definitions.items():
            # The scores are the difference between
            tmp_df[f'delta_{monomer_identity}'] = abs((current_peak  - tmp_df['Mass']) - mass)

        #print(f'Minimum values along axes 0: \n{tmp_df.min(axis=0)}')
        #print(f'Minimum values along axes 1: \n{tmp_df.min(axis=1)}')
        # Get the minimum of each row
        tmp_df['min'] = tmp_df.min(axis=1)
        tmp_df['score'] = tmp_df['min']
        tmp_df['intensity_score'] = tmp_df['Intensity'] / np.sum(tmp_df['Intensity'].to_list())
        tmp_df['intensity_score'] = 1 / tmp_df['intensity_score']
        tmp_df['intensity_score'] = tmp_df['intensity_score'] * (1 / tmp_df['intensity_score'].max())
        tmp_df['score'] = tmp_df['score'] * (1 / tmp_df['score'].max()) # Normalize the score to 1
        tmp_df['score'] = tmp_df['score'] + tmp_df['intensity_score']

        # Get the monomer identity? Might not be necessary
        global_min_column = tmp_df[[x for x in tmp_df.columns if 'delta' in x]].min().idxmin()

        # If we're looking at the first monomer loss, it has to be Phe
        if len(monomers) == 1:
            current_peak = tmp_df.loc[tmp_df['delta_Phe'] == tmp_df[['delta_Phe']].min().iloc[0], 'Mass'].iloc[0]
            monomers[current_peak] = float(tmp_df.loc[tmp_df['Mass'] == current_peak, 'Intensity'].iloc[0])
            continue

        # Check if the endcap is in this range
        elif lower_bound < endcap_mass and upper_bound > endcap_mass:
            if debug:
                print(f'[DEBUG] In endcap range {lower_bound} - {upper_bound}')
            tmp_df['score'] = abs(tmp_df['Mass'] - endcap_mass)
            #display(tmp_df)


        # Reset current peak
        new_peak_to_use = float(tmp_df.loc[tmp_df['score'] == tmp_df['score'].min(), 'Mass'].iloc[0])
        #if debug:
        #    print(f'[DEBUG] Selecting loss of monomer from previous mass (i.e., current_peak) {current_peak}. Global min: {global_min_column}. Min score: {tmp_df["score"].min()}. Selected {new_peak_to_use} ')
        current_peak = new_peak_to_use
        monomers[current_peak] = float(tmp_df.loc[tmp_df['Mass'] == current_peak, 'Intensity'].iloc[0])

    if debug:
        print(f'[DEBUG] Selected masses: {list(monomers.keys())}')

    return monomers, bounds


def chunks(lst, n):
    '''Yield successive n-sized chunks from lst.'''
    for i in range(0, len(lst), n):
        yield lst[i:i + n]



def find_parent_mass(scanData: pd.DataFrame,
                     threshold: float = 0.01,
                     bandwidth: float = 3.5,
                     debug = False) -> float:
    '''
    Finds the largest intensity signal in the group of mass spectra
    signals which are associated with the largest masses measured.

    Parameters
    ----------
    ScanData: pd.DataFrame
        Dataframe containing two columns: Mass, Intensity. Typically the averaged
        mass spectrum of a particular sample (i.e., the average of multiple scans
        which are associated with a spot on a DESI plate)

    intensity_threshold: float (default = 0.02)
        The percentage of intensity of the largest intensity signal
        within the scanData["Intensity"] column above which signals
        are considered for selection.

        Default is 2% of maximum intensity signal

    bandwidth: float
        Width of the kernel used for estimating peak density.

    debug: bool (default = False)
        Turns on various debug printing and plotting

    Returns
    ----------
    parent_peak: float
        The m/z value of the parent peak selected by this function
    '''
    y = scanData['Intensity'].to_numpy()

    # Get the threshold value
    #threshold = threshold * max(y)

    if threshold >= 0.01:
        print(f'[WARNING] Small parent peaks can cause unusual results with thresholds larger than 0.01.')

    threshold = threshold * max(y)

    # Get a copy of the df above the intensity threshold
    tmp_df = scanData[scanData['Intensity'] > threshold].copy(deep=True)

    # 1D array of masses
    a = tmp_df['Mass'].to_numpy()

    # Calculate the density function
    kde = KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(a.reshape(-1,1))

    # Generate linepoints for evaulating with density function
    # we probably wont look at anything larger than 2000 m/z
    s = np.linspace(0, 2000, num=1000)

    # Density function? (less negative values (greater values) indicate greater density of points)
    e = kde.score_samples(s.reshape(-1,1))

    # Get the min and maxima indices of the evaluated line points
    # These are not the indices of the original data
    mi, ma = argrelextrema(e, np.less)[0], argrelextrema(e, np.greater)[0]

    # Determine in the first point in both ma and mi
    # is either a minimum or a maximum. I don't think this is necessary
    if s[ma[0]] < s[mi[0]]:
        first_point_identity = 'max'
    elif s[ma[0]] > s[mi[0]]:
        print('Min point comes first')
        first_point_identity = 'min'
    else:
        raise ValueError('Could not determine if maximum or minimum value came first')

    # The intervals are going to be pairs of minima
    if first_point_identity == 'max':

        dfs_to_plot = []

        splits = chunks(s[mi], 2)
        for i in range(len(s[mi]) - 1):

            # Get the split
            split = s[mi][i: i + 2]

            # If it is the last split
            # this will get the group of signals with the highest mass
            if i == len(s[mi]) - 2:
                temp = scanData[(scanData['Mass'] >= split[1]) & (scanData['Intensity'] >= threshold)]
                parent_peak = float(temp[temp['Intensity'].astype(float) == temp['Intensity'].astype(float).max()]['Mass'].iloc[0])
                parent_peak_intensity = float(temp[temp['Intensity'] == temp['Intensity'].max()]['Intensity'].iloc[0])

    else:
        raise Exception('Code was not made to properly handle min first.')

    #if debug:
    #    fig, axes = plt.subplots(1,2)
    #
    #    axes[0].plot(scanData['Mass'], y, lw=0.5)
    #    axes[0].scatter([parent_peak], [parent_peak_intensity], color='red')
    #    axes[1].plot(s, e, color = 'blue')
    #    plt.savefig('./debug-parent-peak_ident.png', dpi=600)
    #    plt.clf()

    return parent_peak