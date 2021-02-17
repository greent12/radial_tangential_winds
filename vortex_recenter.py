import numpy as np
import sys
import warnings
warnings.filterwarnings("ignore")

from uv_to_rt_winds import uv_to_rt
from recenter_utils import find_nearest,find_nearest_ind
import recenter_utils
from rxpad_check import check_rxpad

#fortran import, could be used but still needs work
#from ideal_angle import create_ideal_angle

def recenter_tc(uwind,vwind,ideal_angle,lons,lats,num_sectors,dist_coeff,wind_coeff,rxpad,spad,num_iterations,olon=None,olat=None):
    
   """
    This function is optimized to re-center TDR merged analyses and/or (hopefully) model analyses
    
    ***
    NOTE: For this function to work, the "ideal_angle" variable should have already been created.
    Here, we assume ideal_angle is a four-dimensional array, shaped proportionally to uwind/vwind.
    
    Furthermore, this function assumes data is on an evenly-spaced grid of value delta_grid.
    ***
    
    Inputs variables are the following:
    
    1) uwind is the zonal component of the TC flow. Should be two-dimensional array ([lat,lon]).
    2) vwind is the meridional component of the TC flow. Should be two-dimensional array ([lat,lon]).
    3) ideal_angle is the angle of the flow associated with a symmetric vortex. 4-D array ([lat,lon,lat,lon]).
    4) lons is a two-dimensional array of the longitudinal coordinates
    5) lats is a two-dimensional array of the latitudinal coordinates
    6) num_sectors is the number (integer) of azimuthal sectors used to compute the mean error (e.g., quadrants = 4).
    7) dist_coeff is the coefficient to weight the distance errors (recommend this to be > than wind_coeff).
    8) wind_coeff is the coefficient to weight the wind speed errors (recommend this to be < than dist_coeff).
    9) rxpad is the radial distance in km past the RMW to weight winds for error calculations 
    10) spad is the number of gridpoints in each direction from the center/2 to search for maximizing weighted 
        tangential wind difference
    11) num_iterations is the amount of times to loop to try to find a center that converges
    12) olon is the first-guess center longitude (float). Default is None, and midpoint of grid will be used.
    13) olat is the first-guess center latitude (float). Default is None, and midpoint of grid will be used.

   Written by Michael Fischer, HRD
   Moddified by Tyler Green:
   First: 11/18/20    
   Last : 12/21/20
   """
    
   #Sectors to process 
   angle_thresh = np.linspace(-1.*np.pi,np.pi,num_sectors+1) # Bounds of azimuthal sectors to average angle errors

   #Tyler changed rad_gaussian creation, look for RMW every 1km
   gr0=1.0
   grf=175
   curr_delta=1.0 # Annulus width used to search for RMW
   rad_gaussian = np.arange(gr0,grf+curr_delta,curr_delta) # Array of radii to use for radius of maximum wind computations

   #Tyler add converged flag for determining if lat/lon needs to be output to None
   converged=False

   ####################################
   #Start recentering
   ####################################

   # Prepare to iterate to find array indices of TC center; initialize as NaNs:
   yloc = np.nan # Meridional index of where TC center was found
   xloc = np.nan # Zonal index of where TC center was found

   # Determine grid indices closest to first-guess TC center, if provided:
   if olon != None:
      try:
         # Find closest lat/lon index to real-time center (olon_np/olat_np):
         lat_close = find_nearest(lats,olat)[0]
         pnyi = np.where(lats == lat_close)[0][0] # meridional index of prior guess center
         lon_close = find_nearest(lons[pnyi,:],olon)[0]
         pnxi = np.where(lons[pnyi,:] == lon_close)[0][0] # zonal index of prior guess center
      except IndexError:
         # If I cannot identify the nearest grid point for whatever reason, just use midpoint of grid:
         pnyi = int(round(lons.shape[0]/2.,0))
         pnxi = int(round(lons.shape[1]/2.,0))

   # Use midpoint of grid if no first-guess center is provided:
   if olon == None:
      pnyi = int(round(lons.shape[0]/2.,0))
      pnxi = int(round(lons.shape[1]/2.,0))

   # Compute wind speed:
   ws = np.sqrt(uwind**2 + vwind**2)

   #Tyler 12/21/20
   #Removed the following variables, as they were unused: angle_dif,weighted_dif, obs_angle

   # Create an array to be filled with mean errors:
   sector_mean_error = np.full((num_sectors,ideal_angle.shape[0],ideal_angle.shape[1]),np.nan,dtype='f4')

   # Set assumed previous error difference to infinity; will be replaced in loop:
   prev_mean_dif =  np.inf

   #Tyler 02/01/21
   # make sure that wind speed within searching area is not all nan's. If so, have to bump up rxpad
   rxpad=check_rxpad(lons,lats,uwind,vwind,pnyi,pnxi,rad_gaussian,curr_delta,rxpad)

   ### Begin TC center search ###

   for n in range(num_iterations):

      print('  recenter_tc: current iteration is:',n)

      # If first iteration has been completed, copy center estimate from previous iteration:
      if n >= 1:
         if prev_mean_dif < np.inf:
            pnyi = np.copy(yloc)
            pnxi = np.copy(xloc)
         else:
            tc_center_lon = np.nan
            tc_center_lat = np.nan
            break
      
      #Tyler 12/21/20 removed the conditional 'if n>=0...' as this will always be true 
      # Loop over meridional grid points:
      for ybi in np.arange(pnyi - spad,pnyi + spad+1,1):
            
         # Ensure the algorithm isn't searching outside the bounds of the domain:
         if ybi >= int(ideal_angle.shape[-2]):
             continue
         if ybi < 0:
             continue

         # Loop over zonal grid points:
         for xbi in np.arange(pnxi - spad,pnxi + spad+1,1):
            # Ensure the algorithm isn't searching outside the bounds of the domain:
            if xbi >= int(ideal_angle.shape[-1]):
               continue
            if xbi < 0:
               continue

            # Ensure the mean error at this grid point hasn't been computed already:
            if n > 1:
               if np.isfinite(np.nanmean(sector_mean_error[:,ybi,xbi])):
                  continue

            # Compute zonal displacement:
            XD = lons - lons[ybi,xbi]
            YD = lats - lats[ybi,xbi]
            curr_dist = recenter_utils.haversine(lons,lats,lons[ybi,xbi],lats[ybi,xbi],deg=True) #Tyler switched the distance calculating function to just haversine, should be faster

            # Compute relative angle:
            angle = np.arctan2(YD,XD) #angle (radians) of each grid point from TC center

            # Calculate tangential wind:
            curr_rw,curr_vt = uv_to_rt(lons,lats,uwind,vwind,lons[ybi,xbi],lats[ybi,xbi])

            ### Determine RMW ###

            # Here we will determine the gaussian distance to be equal to the approximate RMW #
            # Identify the annulus where the tangential wind is maximized:
            vt_ann_azi = np.full((np.size(rad_gaussian)),np.nan,dtype='f4')

            for ri in range(np.size(rad_gaussian)):
               #Tyler 12/21/20, took out the two np.where(np.logical_and statements
               #Only need one logical_and statment to slice up the curr_vt array
               rmw_ann_xy = np.logical_and(curr_dist >= rad_gaussian[ri] - 0.5*curr_delta,\
                                                   curr_dist < rad_gaussian[ri] + 0.5*curr_delta)

               # Compute the maximum value from all annuli:
               vt_ann_azi[ri] = np.nanmean(curr_vt[rmw_ann_xy])

            vt_azi_max = np.nanmax(vt_ann_azi)
              
            try:
               curr_rmw = rad_gaussian[np.where(vt_ann_azi == vt_azi_max)[0][0]]
            except IndexError:
               #print('Could not identify RMW... Looking for vt_azi_max of:',vt_azi_max)
               curr_rmw = 250. # Use a default, broad RMW that should encompass whole domain
                     #continue

            # Use a Gaussian weighted distance:
            # https://en.wikipedia.org/wiki/Gaussian_function

            # Establish the radial distance weighting for errors:
            dist_weight = np.exp(-1.*(((curr_dist - 0.)**2)/(2.*(curr_rmw**2))))

            # Compute sum of weights (use default coefficients here):
            curr_weight = dist_coeff*dist_weight + wind_coeff*(ws/np.nanmax(ws))

            # Compute the angle difference between idealized vortex and observed flow:
            curr_angle_dif = np.arctan2(vwind,uwind) - ideal_angle[ybi,xbi,:,:]

            # Correct for angle differences outside of specified range:
            low_angle_y = np.where(curr_angle_dif < -1.*np.pi)[0]
            low_angle_x = np.where(curr_angle_dif < -1.*np.pi)[1]

            curr_angle_dif[low_angle_y,low_angle_x] = 2.*np.pi + curr_angle_dif[low_angle_y,low_angle_x]

            high_angle_y = np.where(curr_angle_dif > np.pi)[0]
            high_angle_x = np.where(curr_angle_dif > np.pi)[1]

            curr_angle_dif[high_angle_y,high_angle_x] = curr_angle_dif[high_angle_y,high_angle_x] - 2.*np.pi

            # Compute the finalized errors at each grid point:
            curr_weighted_dif = curr_weight*curr_angle_dif

            # Focus on inner core:
            outer_y = np.where(curr_dist > curr_rmw + rxpad)[0]
            outer_x = np.where(curr_dist > curr_rmw + rxpad)[1]

            curr_weighted_dif[outer_y,outer_x] = np.nan
               
            # Separate grid into azimuthal sectors:  
            for thi in range(np.size(angle_thresh) - 1):
               angle_upper = angle_thresh[thi+1]
               angle_lower = angle_thresh[thi]
               xythi = np.logical_and(angle >= angle_lower, angle <= angle_upper)

               sector_mean_error[thi,ybi,xbi] = np.nanmean(np.abs(curr_weighted_dif[xythi]))

            # Compute mean error of each quadrant:
            curr_mean_dif = np.nanmean(sector_mean_error[:,ybi,xbi])

            # Check to see if the current location yields a better center estimate than the previous iteration:
            if curr_mean_dif < prev_mean_dif:

               # Set previous error equal to current error for next iteration:
               prev_mean_dif = np.copy(curr_mean_dif)

               #Tyler 12/21/20 removed the following variables here:
               # angle_dif,weighted_dif,obs_angle, as they were unused to begin with
               # Store TC location estimate:
               tc_center_lon = lons[ybi,xbi]
               tc_center_lat = lats[ybi,xbi]

               yloc = ybi
               xloc = xbi

            # Delete old variables to conserve memory:
            del YD
            del XD
            del curr_dist
            del angle

      # If TC center estimate converges, break from loop and keep estimate:
      if n >= 1:
         if np.logical_and(yloc == pnyi, xloc == pnxi):
            converged=True
            break

   # Delete old variables to conserve memory:
   del sector_mean_error
   del yloc
   del xloc
   del pnyi
   del pnxi

   #If the estimate for tc center did not converge, return None for both
   # This will only happen if this routine is called with num_iterations parameter
   # greater than 1
   if num_iterations >1 and not converged:
      return None,None
 
   return tc_center_lon,tc_center_lat
