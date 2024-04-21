import pandas as pd
import numpy as np
import sklearn
from sklearn.neighbors import KernelDensity
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt


def getGaussian(xdata: np.ndarray, ydata: np.ndarray, sigma: float, line_x: np.ndarray) -> np.ndarray:
    '''
    Generates a set of points for plotting a Gaussian
    broadened signal of a mass spectrum.

    Parameters
    ----------
    xdata: np.ndarray

    ydata: np.ndarray

    sigma: float
        Value which controls the broadness of the
        gaussian curve

    line_x: np.ndarray
        list of integers which represent the span over which
        the guassian is calculated.

    Returns
    ----------
    np.ndarray
        y-values for the calculated function
    '''
    lp=[]
    for line_point in line_x:
        tot=0
        for x,y in zip(xdata,ydata):
            tot+=y*np.exp(-((((x-line_point)/sigma)**2)))
        lp.append(tot)
    return np.array(lp)

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def playGround(
        scanData: pd.DataFrame,
        sigma = 0.8,
        threshold: float = 0.02,
        debug=True) -> pd.DataFrame:
    '''
    Add
    '''

    x, y = scanData['Mass'].to_numpy(), scanData['Intensity'].to_numpy()
    threshold = threshold * max(y)
    print(threshold)
    #ax.stem(x, y, linefmt='blue', markerfmt='', basefmt=" ")

    # Plot Gaussian curves
    #line_x=np.linspace(min(ax.get_xlim()) * 1.1, max(ax.get_xlim()) * 0.9, num=1000, endpoint=True)

    #line_x=np.linspace(min(x), max(x), num=400, endpoint=True)
    #broaden = getGaussian(x, y, sigma=sigma,line_x=line_x)
    #ax.plot(line_x,broaden,"--k", linewidth=0.2)

    # Horizontal line showing the average intensity of the Gaussian function
    #avg_intensity = np.average(y)
    #ax.hlines(y=0.05, xmin=min(x), xmax=max(x), colors='red', label='Avg intensity (gaussian)')

    # Get location of the clusters of peaks on x axis (these are indices)
    #locs = np.where(y>avg_intensity)

    # Get a df for determining density based clustering
    # This does not reset the index
    ######################### HARD CODED LIMIT 0.05 #####################
    tmp_df = scanData[scanData['Intensity'] >= threshold].copy(deep=True)

    # 1D array of masses
    a = tmp_df['Mass'].to_numpy()

    # Calculate the density function
    kde = KernelDensity(kernel='gaussian', bandwidth=3.5).fit(a.reshape(-1,1))

    # Generate linepoints for evaulating with density function
    # This parameter is hard coded because we probably wont look
    # at anything larger than 2000 m/z
    s = np.linspace(0, 2000, num=1000)

    # Density function? (less negative values (greater values) indicate greater density of points)
    e = kde.score_samples(s.reshape(-1,1))
    #e = (((e - min(e)) / (np.average(e)))) * -1
    #e = e - 14


    # Get the min and maxima indices of the evaluated line points
    # These are not the indices of the original data
    mi, ma = argrelextrema(e, np.less)[0], argrelextrema(e, np.greater)[0]

    # Determine in the first point in both ma and mi
    # is either a minimum or a maximum. I don't think
    # this is necessary
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
            #print(i, len(s[mi]) - 1)
            split = s[mi][i: i + 2]
            #print(split)
            if i == 0: # this will color the group of signals with the lowest mass (doesn't seem to work correctly)
                temp = scanData[(scanData['Mass'] >= 0) & (scanData['Mass'] < split[0]) & (scanData['Intensity'] >= threshold)]
            elif i == len(s[mi]) - 2: # this will color the group of signals with the highest mass
                temp = scanData[(scanData['Mass'] >= split[1]) & (scanData['Intensity'] >= threshold)]
            else: # this will get everything else
                temp = scanData[(scanData['Mass'] >= split[0]) & (scanData['Mass'] < split[1]) & (scanData['Intensity'] >= threshold)]

            # This is for coloring EVERYTHING because no intensity threshold
            #dfs_to_plot.append(scanData[(scanData['Mass'] >= split[0]) & (scanData['Mass'] < split[1])]) # For everything

            if not temp.empty:
                dfs_to_plot.append(temp)
            else:
                print(temp.empty, 'a')
    else:
        raise Exception('Code was not made to properly handle min first.')

    fig, ax = plt.subplots(1,1, figsize=(12,8))
    cmap = plt.cm.get_cmap('viridis', len(dfs_to_plot))
    cmap = [z for z in cmap.colors]

    ax.stem(scanData['Mass'], scanData['Intensity'], linefmt='black', markerfmt='', basefmt=" ")
    for i, df in enumerate(dfs_to_plot):
        # Convert to 255
        color = cmap[i] * 255
        color = [int(round(f, 0)) for f in color]
        color = '#{:02x}{:02x}{:02x}'.format(*color)
        if i % 2 == 0:
            color = 'r'
        else:
            color = 'blue'
        ax.stem(df['Mass'], df['Intensity'], linefmt=color, markerfmt='', basefmt=" ")


    plt.show()



    plt.plot(s[:mi[0]+1], e[:mi[0]+1], 'black', # Black line to the first minimum
     s[mi[0]:mi[1]+1], e[mi[0]:mi[1]+1], 'purple', # Purple line extending from first min to next
     s[mi[1]:], e[mi[1]:], 'b', # The remainder is colored blue
     s[ma], e[ma], 'go',
     s[mi], e[mi], 'ro')
    plt.show()
    exit()
    return scanData

def findLargestMassPeak(
        scanData: pd.DataFrame,
        threshold: float = 0.03,
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

def findClusters(
        scanData: pd.DataFrame,
        sigma = 0.8,
        threshold: float = 0.02) -> pd.DataFrame:
    '''
    Splits the scan data into smaller dataframes of the
    clusters of peaks
    '''

    x, y = scanData['Mass'].to_numpy(), scanData['Intensity'].to_numpy()
    threshold = threshold * max(y)

    tmp_df = scanData[scanData['Intensity'] >= threshold].copy(deep=True)

    # 1D array of masses
    a = tmp_df['Mass'].to_numpy()

    # Calculate the density function
    kde = KernelDensity(kernel='gaussian', bandwidth=3.5).fit(a.reshape(-1,1))

    # Generate linepoints for evaulating with density function
    # This parameter is hard coded because we probably wont look
    # at anything larger than 2000 m/z
    s = np.linspace(0, 2000, num=1000)

    # Density function? (less negative values (greater values) indicate greater density of points)
    e = kde.score_samples(s.reshape(-1,1))

    # Get the min and maxima indices of the evaluated line points
    # These are not the indices of the original data
    mi, ma = argrelextrema(e, np.less)[0], argrelextrema(e, np.greater)[0]

    # Determine in the first point in both ma and mi
    # is either a minimum or a maximum. I don't think
    # this is necessary
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
            #print(i, len(s[mi]) - 1)
            split = s[mi][i: i + 2]

            dfs_to_plot.append(scanData[(scanData['Mass'] >= split[0]) & (scanData['Mass'] < split[1])].copy(deep=True)) # For everything
    else:
        raise Exception('Code was not made to properly handle min first.')

    return dfs_to_plot