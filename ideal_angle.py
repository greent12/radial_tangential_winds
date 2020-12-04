import numpy as np

def calc_ideal_angle(lons,lats):

   """
    This function calculates an array of angles for every gridpoint in a 2d grid that represent the flow of an idealized, axisymetric vortex. The output of this function is used in the routine "recenter_tc" and should be called before that routine is called.

    Inputs variables are the following:
    
    1) lons is a two-dimensional array [lat,lon] of the longitudinal coordinates
    2) lats is a two-dimensional array [lat,lon] of the latitudinal coordinates

    Output variables are the following:
   
    1) ideal_angle is a 4d array of holding the angles representing the winds of an idealized vortex. For each grid point in the 2d lat/lon grid, another 2d lat/lon grid of the angles are stored, making it 4d.

    Written by Michael Fischer, HRD
    Modified by Tyler Green, 11/19/20
   """

   ideal_angle = np.full((lons.shape[0],lons.shape[1],lons.shape[0],lons.shape[1]),np.nan,dtype='f4')

   # Loop over meridional grid points:
   for ybi in range(ideal_angle.shape[0]):

      
      # Loop over zonal grid points:
      for xbi in range(ideal_angle.shape[1]):

         #Tyler calculate XD and YD meshgrid using great circle distance from current lat and lon point
         XD = lons - lons[ybi,xbi]
         YD = lats - lats[ybi,xbi]

         # Compute the angle of wind in the idealized vortex:
         ideal_angle[ybi,xbi,:,:] = np.arctan2(YD,XD) + np.pi/2.

         # Correct for instances where the angle exceeds pi:
         highy = np.where(ideal_angle[ybi,xbi,:,:] > np.pi)[0]
         highx = np.where(ideal_angle[ybi,xbi,:,:] > np.pi)[1]

         ideal_angle[ybi,xbi,highy,highx] = ideal_angle[ybi,xbi,highy,highx] - 2.*np.pi
         #Tyler add: at center point, angle should be a singularity, make angle nan
         ideal_angle[ybi,xbi,ybi,xbi] = np.nan

   return ideal_angle


