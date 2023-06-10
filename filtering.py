import pandas as pd

def filterLowMasses(ScanData: pd.DataFrame) -> pd.DataFrame:
    '''
    Filters a ScanData pd.DataFrame for masses less than 
    1,254.6 m/z units for masses that are even possible with
    a particular set of monomers

    Parameters
    ----------
    ScanData: pd.DataFrame
        hass a Mass column and Intensity column

    Returns
    ----------
        pd.DataFrame
    '''
    with open('./possible_masses.txt', 'r') as infile:
        possibleMasses = infile.readlines()
        possibleMasses = [float(x.strip('\n')) for x in possibleMasses]

    ScanData['Possible'] = ScanData['Mass'].isin(possibleMasses)
    return ScanData[(ScanData['Possible'] == True) | (ScanData['Mass'] >= 1254)]