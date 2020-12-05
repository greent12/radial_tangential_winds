import numpy as np
import sys
import warnings
warnings.filterwarnings("ignore")
import xarray as xr
import os
import matplotlib.pyplot as plt 

#My imports
from uv_to_rt_winds import uv_to_rt
from ideal_angle import calc_ideal_angle
from vortex_recenter import recenter_tc
from recenter_utils import cut_latlons,cut_uv_grids

from main_tools \
   import read_namelist,initialize_read_data, extract_vertical_coords, \
          get_vert_indicies,create_output_directory,extract_uv_winds, \
          plot_winds_pressure_heights,save_data

###############################################################################
# START CODE EXECUTION
###############################################################################

#Read in namelist
inputs,vortex_parms = read_namelist("namelist.input")

#Initialize
data,vertical_name = initialize_read_data(inputs['input_file'],\
                                          inputs['vertical_option'])
#Extract coordinate values from data
verts,lats,lons,latlon_dim = extract_vertical_coords(data,vertical_name)

#Find indicies of vertical levels wanted
vert_indicies = get_vert_indicies(verts,inputs['all_levels'],inputs['levels'],\
                                  vertical_name)

#Create output directory if it does not exist
create_output_directory(inputs['output_dir'])

#Extract u and v winds (3d) from data
u3d,v3d = extract_uv_winds(data)

#Cut u and v grid so that they can be processed in the vortex finding techinque
#extract lats,lons,u, and v and cut grid
lons2d,lats2d,lonscut,latscut,ucut,vcut=cut_uv_grids(lons,lats,latlon_dim,u3d,v3d,inputs['lon_guess'],inputs['lat_guess'],vortex_parms['grid_cut'])

#Create array of ideal vortex wind angles
print("Creating ideal_angle...")
ideal_angle = calc_ideal_angle(lonscut,latscut)
print("Done creating ideal angle")

#Initialize arrays to hold vertical levels processed, found tc lat and lons, and twind/rwind
numvert = len(vert_indicies)
[K,J,I] = np.shape(u3d)
verts_processed=np.full(numvert,np.nan,dtype=float)
tc_lats = np.full(numvert,np.nan,dtype=float)
tc_lons = np.full(numvert,np.nan,dtype=float)
twind = np.full((numvert,J,I),np.nan,dtype=float)
rwind = np.full((numvert,J,I),np.nan,dtype=float)

#Loop through vertical levels
print("Start processing of input file")
for zcounter,zi in enumerate(vert_indicies):
   
   print("Working on level: {} {}".format(verts[zi],vertical_name))

   #Slice winds for this level
   uwind = ucut[zi,:,:]
   vwind = vcut[zi,:,:]

   #Find tc center
   tc_lon,tc_lat= recenter_tc(uwind,vwind,ideal_angle,lonscut,latscut,\
                              vortex_parms['num_sectors'],\
                              vortex_parms['dist_coeff'],\
                              vortex_parms['wind_coeff'],\
                              vortex_parms['rxpad'],\
                              inputs['lon_guess'],\
                              inputs['lat_guess'])

   #Append lat, lon, and vertical level
   tc_lats[zcounter] = tc_lat
   tc_lons[zcounter] = tc_lon
   verts_processed[zcounter] = verts[zi]

   #Calculate tangential wind of orginal horizontal field if tc_lat/lon was found
   if tc_lon is not None and tc_lat is not None:
      rwind_zi,twind_zi =  uv_to_rt(lons2d,lats2d,u3d[zi,:,:],v3d[zi,:,:],tc_lon,tc_lat,deg=True)
      rwind[zcounter,:,:] = rwind_zi
      twind[zcounter,:,:] = twind_zi

   #Make plot
   if inputs['plot']:
      plot_winds_pressure_heights(lons2d,lats2d,verts,zi,u3d,v3d,data,vertical_name,tc_lon,tc_lat,inputs['output_dir'])
 
#Save data to nc file
save_data(verts_processed,vertical_name,lats,lons,rwind,twind,tc_lons,tc_lats,inputs['output_dir'],inputs['output_file'])

print("Done")
