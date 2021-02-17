import numpy as np
import recenter_utils
from uv_to_rt_winds import uv_to_rt

def check_rxpad(lons,lats,uwind,vwind,ybi,xbi,rad_gaussian,curr_delta,rxpad_org):

   thresh=0.75

   #Get current distance from center
   curr_dist = recenter_utils.haversine(lons,lats,lons[ybi,xbi],lats[ybi,xbi],deg=True)

   #Get tangential wind
   curr_rw,curr_vt = uv_to_rt(lons,lats,uwind,vwind,lons[ybi,xbi],lats[ybi,xbi])

   # Identify the annulus where the tangential wind is maximized
   vt_ann_azi = np.full((np.size(rad_gaussian)),np.nan,dtype='f4')

   for ri in range(np.size(rad_gaussian)):
      rmw_ann_xy = np.logical_and(curr_dist >= rad_gaussian[ri] - 0.5*curr_delta,\
                                  curr_dist < rad_gaussian[ri] + 0.5*curr_delta)

      # Compute the maximum value from all annuli:
      vt_ann_azi[ri] = np.nanmean(curr_vt[rmw_ann_xy])

   vt_azi_max = np.nanmax(vt_ann_azi)
   curr_rmw = rad_gaussian[np.where(vt_ann_azi == vt_azi_max)[0][0]]

   #Find rxpad that gives at least 25% of valid grid point that are not nan's
   rxpad=rxpad_org
   valid_vt = curr_vt[curr_dist<=curr_rmw+rxpad]
   npoints_nan = np.sum(np.isnan(valid_vt))
   percent_points_nan = npoints_nan / np.size(valid_vt)

   while percent_points_nan>=thresh:
      if rxpad+curr_rmw > np.nanmax(curr_dist):
         break

      rxpad+=5.0
      valid_vt = curr_vt[curr_dist<=curr_rmw+rxpad]
      npoints_nan = np.sum(np.isnan(valid_vt))
      percent_points_nan = npoints_nan / np.size(valid_vt)

   #Done
   if rxpad != rxpad_org:
      print("rxpad_check: new rxpad: {} does not match original: {}".format(rxpad,rxpad_org))

   return rxpad

