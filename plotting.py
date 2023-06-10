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

def plotSpots(n_scans: int, AvgScanIntensity: list[float]) -> None:
    plt.plot(np.arange(n_scans), AvgScanIntensity, color='black')
    plt.title('Spot Diagram')
    plt.xlabel('Scan Number')
    plt.ylabel('Average Intensity of All Massess')
    plt.show()