from pathlib import Path

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

def plotTwoMassSpectra(
        ms1: pd.DataFrame,
        ms2: pd.DataFrame) -> None:
    '''
    Plots two mass spectra which are represented by two
    pd.DataFrames

    Parameters
    ----------
    ms1: pd.DataFrame
        Dataframe with a "Mass" and "Intensity" columns

    ms2: pd.DataFrame
        Dataframe with a "Mass" and "Intensity" columns

    Returns
    ----------
        None
    '''

    fig, ax = plt.subplots(1,1)
    #ax.set_title(f"Overlay")

    x, y = ms1['Mass'], ms1['Intensity']
    ax.stem(x, y, linefmt='blue', markerfmt='')

    x, y = ms2['Mass'], ms2['Intensity']
    ax.stem(x, y, linefmt='red', markerfmt='')

    plt.show()

def plotMassSpectrum(
        ms1: pd.DataFrame) -> None:
    '''
    Plots  mass spectra which is represented as a
    pd.DataFrames

    Parameters
    ----------
    ms1: pd.DataFrame
        Dataframe with a "Mass" and "Intensity" columns

    Returns
    ----------
        None
    '''

    fig, ax = plt.subplots(1,1, figsize=(12,8))
    #ax.set_title(f"Overlay")

    x, y = ms1['Mass'], ms1['Intensity']
    ax.stem(x, y, linefmt='blue', markerfmt='', basefmt=" ")

    plt.show()


def plot_spots(average_scan_intensities: list[float],
               isolated_spots: list[list[int]],
               save: Path,
               ):

    '''
    Parameters
    ----------
    average_scan_intensities: list[float]
        List of the ion intensity (counts) of a
        series of scans.

    isolated_spots: list[list[int]]
        List of list of ints where the inner ints
        are the indices of the scan intensities that
        are associated with a spot. The inner lists
        are the number of spots.

    save: Path
        Path to the saved image

    Returns
    ----------
        None
    '''
    fig, ax = plt.subplots(1,1)
    ax.plot(np.arange(0,len(average_scan_intensities)), average_scan_intensities, color='lightgray')
    for _i, scan_group in enumerate(isolated_spots):
        ax.plot(scan_group, [average_scan_intensities[z] for z in scan_group])
    ax.set_title('DESI scan of total intensity')
    ax.set_xlabel('Scan number (roughly equal to time)')
    ax.set_ylabel('Total scan intensity (ion count??)')
    plt.savefig(save, dpi=600)