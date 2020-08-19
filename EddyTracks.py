"""
EddyTracks.py

Author: Mark Yamane

Functions created for eddy data processing and statistics in
Yamane MT, Claret M, Lelong P
"Near-inertial wave interactions with cyclonic and anticyclonic eddies in the northwestern Mediterranean"
"""

import h5py
import numpy as np

def EddyTracks(fpath='/'):
    """
    Used to process MATLAB structure data from from the AMEDA eddy tracking algorithm into a Python list using dictionary
    comprehensions and NumPy.
    
    Input:
        fpath     file path for eddy_tracks.mat file, default is current directory
    Output:
        tracks    list holding variables for each eddy track (e.g., tracks[0] is a dictionary of variables for the first eddy)
    """
    f = h5py.File(fpath + 'eddy_tracks.mat')
    data = f['tracks']            # relevant variables
    varnames = list(f['tracks'])  # variable names
    nTracks = len(data['x1'])     # number of tracks (could work with any variable)

    # add dictionaries for each variable to the corresponding eddy index (0 to N-1)
    tracks = []
    for itrack in range(nTracks):
        variables = {var:[] for var in varnames}  # initialize a dictionary with each variable name
        for var in varnames:
            # go through each variable
            if var[0:6] == 'shapes':
                # translate shapes data to list of numpy arrays (one array for each time step)
                coordinates = []                             # initialize list of coordinates
                numSteps = len(f[data[var][itrack,0]][(0)])  # number of timesteps
                # each timestep has its own hdf5 object reference, so append each step to a list
                for step in range(numSteps):
                    coords = np.array(f[f[data[var][itrack,0]][(0, step)]])
                    coordinates.append(coords)
                variables[var] = coordinates
            else:
                # translate data from hdf5 object reference to numpy array
                variables[var] = np.array(f[data[var][itrack,0]][(0)])
        tracks.append(variables)
    return tracks

def DifferentiateEddies(tracks, lifetime=0):
    """
    Differentiate eddies by rotation. 'lifetime' can be used to filter for eddies with a minimum lifetime (default is 0).
    
    Input:
        tracks    data for each track (output from EddyTracks)
    Output:
        cTracks   list of cyclones from tracks
        aTracks   list of anticyclones from tracks
    """
    cTracks = []  # initialize list for cyclones
    aTracks = []  # initialize list for anticyclones
    
    for i, track in enumerate(tracks):
        # iterate through all tracks
        timesteps = track['step'][-1] - track['step'][0]+1
        if timesteps >= lifetime and track['type'][0] == 1.:
            # add to list of cyclonic eddies
            cTracks.append(track)
        if timesteps >= lifetime and track['type'][0] == -1.:
            # add to list of anticyclonic eddies
            aTracks.append(track)
    return cTracks, aTracks

def fillCoords(trackSteps, trackLons, trackLats, fill='linear'):
    """
    Find gaps in timesteps, longitudes, and latitudes and fills coordinates using the specified interpolation
    (default fills points linearly).
    
    Input:
        trackSteps    steps from tracks
        trackLons     longitudes corresponding to trackSteps
        trackLats     latitudes corresponding to trackSteps
        fill          type of interpolation for lons and lats (midpoint, linear, begin, end)
    Output:
        steps         filled steps
        lons          filled longitudes
        lats          filled latitudes
    """
    prev = trackSteps[0]
    steps = np.array([])
    lons = []
    lats = []
    for i, step in enumerate(trackSteps):
        lon = trackLons[i]
        lat = trackLats[i]
        if step - prev > 1:
            # there is a gap
            stepFill = np.arange(prev+1, step, 1.)
            numFill = len(stepFill)
            if fill == 'midpoint':
                # fill using mid-point of gap
                lonFill = np.ones(len(stepFill))*((lon + trackLons[i-1])/2)
                latFill = np.ones(len(stepFill))*((lat + trackLats[i-1])/2)
            elif fill == 'linear':
                # fill with linearly spaced positions between gap
                lonFill = np.linspace(trackLons[i-1], lon, num=numFill)
                latFill = np.linspace(trackLats[i-1], lat, num=numFill)
            elif fill == 'begin':
                # fill with beginning position of gap
                lonFill = np.ones(len(stepFill))*(trackLons[i-1])
                latFill = np.ones(len(stepFill))*(trackLats[i-1])
            elif fill == 'end':
                # fill with end position of gap
                lonFill = np.ones(len(stepFill))*lon
                latFill = np.ones(len(stepFill))*lat
            else:
                raise ValueError('Invalid fill type.')
            steps = np.append(steps, stepFill)
            lons = np.append(lons, lonFill)
            lats = np.append(lats, latFill)
        steps = np.append(steps, step)
        lons = np.append(lons, lon)
        lats = np.append(lats, lat)
        prev = step
    steps.flatten()
    lons.flatten()
    lats.flatten()
    return steps, lons, lats


############################################################
"""Functions for eddy statistics"""

def haversine(lon1, lat1, lon2, lat2):
    """Calculates distance between two geographic coordinates (pos1, pos2) using the haversine formula"""
    R = 6373.                                      # radius of Earth (in km)
    coords = np.float64([lon1, lat1, lon2, lat2])  # ensure coords are not integers
    lon1, lat1, lon2, lat2 = np.radians(coords)    # convert to radians
    
    dlon = lon2 - lon1  # distance btwn longitudes
    dlat = lat2 - lat1  # distance btwn latitudes
    
    # haversine formula
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    return R*c

def geodist(lons, lats):
    """Calculates total distance travelled from lists of longitudes and latitudes"""
    if len(lons) != len(lats):
        raise ValueError('Coordinate lists must have the same length.')
        
    dist = 0.
    for i in range(len(lons)):
        # sum distances
        currLon = lons[i]
        currLat = lats[i]
        if i > 0:
            dist += haversine(prevLon, prevLat, currLon, currLat)
        prevLon = currLon
        prevLat = currLat
    return dist