#Module containing functions used in the finding of vortex center and calculation of tangential winds

#Imports 
from math import radians, cos, sin, asin, sqrt
import numpy as np

# Function to find index of element closest to specified value:
def find_nearest(array, value):
    ''' Written by Michael Fischer, HRD'''
    array = np.asarray(array)
    X = np.abs(array - value)
    idx = np.where( X == np.nanmin(X) )
    return array[idx]

def find_nearest_ind(array, value):
    ''' Written by Michael Fischer, HRD'''
    array = np.asarray(array)
    X = np.abs(array - value)
    idx = np.where( X == np.nanmin(X) )
    return idx

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)

    Output is  in km

    Taken from: https://stackoverflow.com/questions/4913349/haversine-formula-in-python-bearing-and-distance-between-two-gps-points
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles
    return c * r

def calc_distance_from_point(lons,lats,lon0,lat0):

   ''' Calculate great circle distance of all lat/lon pairs from a reference point 
 
       lons and lats need to be 2d for this to work
   '''

   [I,J] = np.shape(lons)
   dist = np.zeros_like(lons)
   for i in range(I):
      for j in range(J):
         dist[i,j] = haversine(lons[i,j],lats[i,j],lon0,lat0)

   return dist

def cut_latlons(lons,lats,latlon_dim,lon0,lat0,npoints):

   '''
    This function cuts a grid of lats and lons to be smaller for use when calculating ideal_angle variable when doing the vortex center finding. The reason this needs to be done is that if the grid is too large, python will run out of memory when creating ideal_angle since it is 4d.

    Parameters
    ----------
    lons : numpy.ndarray
       1 or 2d array of lons
    lats : numpy.ndarray
       1 or 2d array lf lats
    latlon_dim : int
       dimensions of lats/lons
    lon0 : float (can be None)
       first guess longitude of tc center 
    lat0 : float (can be None)
       first guess latitude of tc center
    npoints : int
       number of gridpoints wide/2 in each direction to cut the grid

    Returns
    -------
    lons_bool : numpy.ndarray of bools
       array of bools size of longitude direction representing whether gridpoint will be kept or not
    lats_bool : numpy.ndarray of bools
       array of bools size of latitude direction representing whether gridpoints will be kept or not

   '''

   #If lat/lon dimension is 1, we need to make a meshgrid
   if latlon_dim == 1:
      lats1d = lats
      lons1d = lons
      lons2d,lats2d = np.meshgrid(lons,lats)
   #lats and lons are already in a mesh
   # this is probably a lambert conformal grid, so were going to make a simplification and create a new lat/lon mesh using the the first column of lats and the first row of lons. "Good enough for government work - Howie Bluestein"
   else:
      lats1d = lats[:,0]
      lons1d = lons[0,:]
      lons2d,lats2d = np.meshgrid(lons1d,lats1d)

   #Set these flags, have to be turned on for code to return properly
   lats_done=False
   lons_done=False

   #First see if the grid dimensions are smaller or the same as npoints*2
   [nlat,nlon] = np.shape(lons2d)
   if nlat < npoints*2:
      lats_bool = np.ones_like(lats1d,dtype=bool)
      lats_done=True
   if nlon < npoints*2:
      lons_bool = np.ones_like(lons1d,dtype=bool)
      lons_done=True 
   
   #lon0 and lat0 could be None, if so the grid will just be cut around the center
   if lat0 is None:
      if not lats_done:
         lats_bool = np.zeros_like(lats1d,dtype=bool)
         lat_center_indx = int(len(lats1d)/2)
         lats_bool[lat_center_indx-npoints:lat_center_indx+npoints] = 1
         lats_done=True

   if lon0 is None:
      if not lons_done:
         lons_bool = np.zeros_like(lons1d,dtype=bool)
         lon_center_indx = int(len(lons1d)/2)
         lons_bool[lon_center_indx-npoints:lon_center_indx+npoints] = 1
         lons_done=True
       
   #If cutting is done, we can return now 
   if lats_done and lons_done:
      return lons_bool,lats_bool
 
   #Calculate distance of all lat/lon pairs from given center
   r = calc_distance_from_point(lons2d,lats2d,lon0,lat0)
   [I,J] = np.where(r == np.min(r))
   I=I[0]
   J=J[0] 

   #If lats or lons aren't done, cut around the girdpoint closest to the first guess lat/lon point
   if not lats_done:
      lats_bool = np.zeros_like(lats1d,dtype=bool)
      lats_bool[I-npoints:I+npoints]=1
      lats_done = True
   
   if not lons_done:
      lons_bool = np.zeros_like(lons1d,dtype=bool)
      lons_bool[J-npoints:J+npoints] = 1
      lons_done = True

   #Sanity check that both lats and lons were cut
   if not lats_done and not lons_done:
      print("lats or lons were not cut right, exitting...")
      sys.exit()

   return lons_bool,lats_bool


def cut_uv_grids(lons,lats,latlon_dim,u3d,v3d,lon_guess,lat_guess,grid_cut):

   '''
    This function cuts the u and v grids to a smaller size so that they can be processed by the vortex recentering routine.

    Parameters
    ----------
    lons : numpy.ndarray
       1 or 2d array of lons
    lats : numpy.ndarray
       1 or 2d array lf lats
    latlon_dim : int
       dimensions of lats/lons
    u3d : numpy.ndarray
       3d gird of u winds (vertical,lat,lon)
    v3d : numpy.ndarray
       3d grid of u winds (vertical,lat,lon)
    lon_guess : float
       first guess of TC lon center (can be None)
    lat_guess : float
       first guess of TC lat center (can be None)
    grid_cut : int
       width of grid in both dimensions/2 to cut 

    Returns
    -------
    lons2d : numpy.ndarray
       2d mesh of longitudes (original grid) 
    lats2d : numpy.ndarray
       2d mesh of latitudes (original grid)
    lonscut : numpy.ndarray
       2d mesh of cut longitudes
    latscut : numpy.ndarray
       2d mesh of cut latitudes
    ucut : numpy.ndarray
       3d grid of cut u winds (vert,latscut,lonscut)
    vcut : numpy.ndarray
       3d grid of cut v winds (vert,latscut,lonscut)

   '''
   #Get bools for lats and lons for cutting grid
   lons_bool,lats_bool = cut_latlons(lons,lats,latlon_dim,lon_guess,lat_guess,grid_cut)

   #Lats and lons need to be in a mesh, if not already
   if latlon_dim == 1:
      lons2d,lats2d = np.meshgrid(lons,lats)
   else:
      lons2d = lons
      lats2d = lats

   #Cut lats and lons
   lonscut = lons2d[np.ix_(lats_bool,lons_bool)]
   latscut = lats2d[np.ix_(lats_bool,lons_bool)]
 
   #Cut u and v grids
   numvert = np.shape(u3d)[0]
   ucut = u3d[np.ix_(np.ones(numvert,dtype=bool),lats_bool,lons_bool)]
   vcut = v3d[np.ix_(np.ones(numvert,dtype=bool),lats_bool,lons_bool)]

   return lons2d,lats2d,lonscut,latscut,ucut,vcut
