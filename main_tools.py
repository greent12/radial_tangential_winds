import sys
import os
from vertical_options import vertical_dict
import xarray as xr
import numpy as np
import f90nml
from netCDF4 import Dataset

def read_namelist(namelist):
   '''
    This function reads a fortran namelist for parameters used in the main program

    Parameters
    ----------
    namelist : str
      fortran namelist file

    Returns
    -------
    main_inputs : dict
      dictionary of namelist inputs used for main code
    vortex_inputs : dict
      dictionary of namelist inputs uses specifically for the vortex recentering routine

    Notes
    -----
    There are 3 mandatory namelist parameters for the  "&inputs" section, which are 'input_file', 'all_levels', and 'vertical_option'
    There are 5 mandatory namelist parameters for the  "&vortex_finder_input" section which are 'grid_cut', 'num_sectors','dist_coeff','wind_coeff', and 'rxpad'

    List of possible namelist parameters:
    -------------------------------------

    &inputs
     1. input_file (string) : full path to input grib file to read 
     2. all_levels (bool)   : whether to process all vertical levels in input_file or not, if this is ".false.", then the namelist parmeter: 'levels' needs to be set
     3. vertical_option (string) : the name of the vertical coordinate, the options for this are located in vertical_options.py. One of the dictonary keys should be chosen, which then maps to what that vertical coorindate is called in grib format. If needed, options can be added to this dictionary for vertical coordinates that are not currently present
     4. levels (float or list of floats) : comma seperated list of vertical levels to process. These should correspond to whatever vertical coordinate chosen. These will only be used if 'all_levels' namelist parameter is ".false."
     5. output_dir (string) : path to directory for outputing plots and nc file. If it is not provided, it will take the directory that 'input_file' is in
     6. lat_guess (float) : first guess for TC center latitude, if you don't want to make a guess, just leave it out of the namelist and it will be set to 'None'
     7. lon_guess (float) : same as above but for longitude
     8. plot (bool) : whether to make plots or not for tangential wind at each vertical level, if not in namelist, it is set to 'False'
     9. output_file (string) : filename that will be used to write tangential/radial & lat/lon centers to. It nees to have a .nc extension. If a name is not provided, the input_file will be appended with a "_rtwind.nc" extension. This file is written to the 'output_dir'

   &vortex_finder_input
    1. grid_cut (integer) : width of the grid to cut/2 (in both dimensions lat/lon) this needs to be relatively small to avoid running out of RAM, 75 seems to work well
    2. num_sectors (integer) : see 'vortex_recenter.py' docstring
    3. dist_coeff  (float) : ' '
    4. wind_coeff (float) : ' '
    5. rxpad (float) : ' '

    A full list of namelist option are available in "sample_namelist"

   '''

   print("Reading namelist...")

   #Establish list of mandatory namelist parameters, if one of these is missing code will exit
   mandatory_main = ["input_file","all_levels","vertical_option"]
   mandatory_vortex = [ "grid_cut","num_sectors","dist_coeff","wind_coeff","rxpad" ]

   #Try reading the namelist, if there's an exception, code will exit
   try:
      nml = f90nml.read(namelist)
   except: 
      print("  Cannot read namelist file, make sure it exists...")
      sys.exit()
   
   #Two namelist sections, one for main inputs and one for vortex finding parameters
   main_inputs = dict(nml['inputs'])
   vortex_inputs = dict(nml['vortex_finder_input'])

   #Scan for mandatory arguments
   for arg in mandatory_main:
      if not arg in main_inputs.keys():
         print("  Mandatory argument: {} missing in {}, exiting...".format(arg,namelist))
         sys.exit()

   for arg in mandatory_vortex:
      if not arg in vortex_inputs.keys():
         print("  Vortex argument: {} missing in {}, exiting...".format(arg,namelist))
         sys.exit()

   #If "all_levels" is False and no "levels" parameter is given, this is a problem
   if not main_inputs['all_levels'] and "levels" not in main_inputs:
      print("  'all_levels' parameter is False but no levels specified in {} ... exiting".format(namelist))
      sys.exit()

   #Add "levels" to inputs dictionary if not there, also typecast this to a list if there is only one level given
   if "levels" not in main_inputs:
      main_inputs.update({"levels":[]})
   else:
      if type(main_inputs['levels']) is not list:
         main_inputs['levels'] = [ main_inputs['levels'] ]
      else:
         main_inputs['levels'] = sorted(main_inputs['levels'])

   #If output directory is not given, it will use directory that the input file is in
   if "output_dir" not in main_inputs:
      mysplit=main_inputs['input_file'].split("/")
      mysplit.pop(-1)
      output_dir = ""
      for dir in mysplit:
         output_dir+=dir+"/"
      main_inputs.update({"output_dir" : output_dir}) 
      print("  No 'output_dir' specified, using 'input_file' directory")

   #If no initial lat guess is specified, set to None
   if "lat_guess" not in main_inputs:
      main_inputs.update({"lat_guess" : None})
      print("  No 'lat_guess' provided, setting to None")

   #If no initial lon guess is specified, set to None
   if "lon_guess" not in main_inputs:
      main_inputs.update({"lon_guess" : None})
      print("  No 'lon_guess' provided, setting to None")

   #If no "plot" parameter is specified, set to False
   if "plot" not in main_inputs:
      main_inputs.update({"plot" : False})
      print("  No 'plot' provided, setting to False")

   #If no name for output nc file is given, just use input file name with "_trwind.nc" appended to the end 
   if "output_file" not in main_inputs:
      output_file=main_inputs['input_file'].split("/")[-1].split(".grib")[0]+"_trwind.nc"      
      print("  No 'output_file' provided, setting to: {}".format(output_file))
      main_inputs.update({"output_file" : output_file})
   else:
      #Check if there is a ".nc" extension
      output_file = main_inputs['output_file']
      if not output_file.endswith(".nc"):
         print("  Make sure 'output_file' ends with '.nc' ... exiting...")
         sys.exit()

   return main_inputs,vortex_inputs


def initialize_read_data(file,vertical_opt):

   '''
    This function does some initialization before calculations are done in finding t and r winds. Specifically, it makes sure that the input file exists, the chosen vertical coordinate is available (see vertical_options.py), and then reads the input file with xarray to find variables with the chosen vertical coordinate.

    Parameters
    ----------
    file : str
           path to input grib file to read from
    vertical_opt : str
           choice for vertical coordinate for variables, see dict keys in vertical_options.py

    Returns
    -------
    data_array : xarray.core.dataset.Dataset
           xarray dataset containing grib information from input file
    vertical_grib_name : str
           vertical coordinate name as it would be stored in the grib file

   '''

   #Check to make sure file exists
   if not os.path.exists(file):
      raise FileNotFoundError("{} could not be found".format(file))
   else:
      print("File found: {}, continuing...".format(file))

   #Check if vertical level option is available
   if vertical_opt in vertical_dict:
      vertical_grib_name=vertical_dict[vertical_opt]
   else:
      print("Selected vertical level: {} is not available, please check ""vertical_options.py"" and retry/modify this code accordinly".format(vertical_opt))
      sys.exit()

   #Open file with xarray
   print("Opening file for reading...")
   data_array=xr.open_dataset(file,engine='cfgrib',
                        backend_kwargs={'filter_by_keys':
                                       {'typeOfLevel':vertical_grib_name}})
   print("Done reading file")

   return data_array,vertical_grib_name


def extract_vertical_coords(data_array,vertical_grib_name):
     
   '''
    This function extracts the the coordinate values for vertical coordinate, lats and lons from an input data_array.

    Parameters
    ----------
    data_array : xarray.core.dataset.Dataset
           data array read from grib file 
    vertical_grib_name : str
           vertical coordinate choice in grib file

    Returns
    -------
    vert : numpy.ndarray 
           1d array of vertical levels
    lats : numpy.ndarray
           1/2 d array of latitudes
    lons : numpy.ndarray
           1/2 d array of longitudes
    lats.ndim : int
           dimension of lats/lons

   '''

   #Exctract vertical levels,lats, and lons 
   vert = data_array.coords[vertical_grib_name].values
   vert = vert.astype(np.float32)
   lats = data_array.latitude.values
   lons = data_array.longitude.values

   #Check to make sure dimensions of lats and lons exist
   if lats.ndim != lons.ndim:
      print("dimension of lats does not match dimension of lons, there seems to be a problem, exitting...")
      sys.exit()

   return vert,lats,lons,lats.ndim


def get_vert_indicies(verts,all_levels,levels,vertical_grib_name):

   '''
    This function gets the indicies of wanted vertical levels to process.

    Parameters
    ----------
    verts : numpy.ndarray
            array of vertical levels to look in 
    all_levels : bool
            whether to process all vertical levels or not
    levels : list/numpy 1d array
            vertical levels to try to find in "verts"
    vertical_grib_name : str
            name of vertical coordinate in grib file

    Returns
    -------
    vert_indxs : list
            indicies of vertical levels to process

    Notes
    -----
    if no vertical levels wanted could be found, this script will exit
   '''

   if all_levels:
      print("Processing all vertical levels in file")
      vert_indxs = list(range(len(verts)))
      return vert_indxs
   else:
      vert_indxs = []
      for lev in levels:
         try: 
            vert_indxs.append(np.where(verts == lev)[0][0])
            print("Found level: {} {} in input file".format(lev,vertical_grib_name))
         except:
            print("Could not find level: {} {} in input file".format(lev,vertical_grib_name))
   if not vert_indxs:
      print("Could not find any levels to process, exitting...")
      sys.exit()

   return vert_indxs
      
def create_output_directory(output_dir):

   ''' Creates directory if it does not already exist. Argument is a str filepath '''
   if not os.path.exists(output_dir):
      try:
         os.mkdir(output_dir)
         print("Created output directory: {}".format(output_dir))
      except:
         print("Output directory cannot be made, likely that upward directory does not exist...exitting")
         sys.exit() 
   else:
      print("Output directory already exists")

def extract_uv_winds(data_array):
   uwind_3d = data_array.u.values
   vwind_3d = data_array.v.values

   return uwind_3d,vwind_3d

def save_data(verts,vertname,lats,lons,rwind,twind,tc_lons,tc_lats,output_dir,filename):
   
   ''' 
    This function saves t/r wind and found tc lats and lons into a netcdf file
   
    Parameters
    ----------
    verts : list
       list of vertical levels that were processed (will be used as coordinate in nc file)
    vertname : str
       name of the vertical coordiante used
    lats : numpy.ndarray
       1 or 2 D array of latitudes (will be used as coordinate in nc file)
    lons : numpy.ndarray
       1 or 2 D array of longitudes (will be used as coordinate in nc file)
    rwind : numpy.ndarray
       3D array of radial wind that was processed, vertical coordinate should be the same dimension of 'verts'
    twind : numpy.ndarray
       Same as above but for tangential wind
    tc_lons : numpy.ndarray
       1D array of found TC longitude centers, if one wasn't found for a specific level it will be 'nan'
    tc_lats : numpy.ndarray
       Same as above but for lats
    output_dir : str
       Path to directory to write nc file to
    filenmae : str
       Name of nc file to write to
    
   Returns
   -------
   nothing

   Notes
   -----
   all missing values for each variable stored in the nc file will be np.nan

  '''

   #Put filepath together    
   filepath = output_dir+"/"+filename
   ncfile = Dataset(filepath,mode="w",format='NETCDF4')

   #nc atributes
   ncfile.title='Tangential and Radial Wind Components'

   #Create dimensions and corresponding variables
   vert_dim = ncfile.createDimension(vertname,len(verts))
   vert_var = ncfile.createVariable(vertname,np.float32,(vertname))
   vert_var[:] = verts

   if lats.ndim==1 and lons.ndim==1:
      lat_dim = ncfile.createDimension('lat',len(lats))
      lat_var = ncfile.createVariable('lats',np.float32,('lat'))
      lat_var[:] = lats

      lon_dim = ncfile.createDimension('lon',len(lons))
      lon_var = ncfile.createVariable('lons',np.float32,('lon'))
      lon_var[:] = lons
   else:
      [nlat,nlon] = np.shape(lats)

      lat_dim = ncfile.createDimension('lat',nlat)
      lon_dim = ncfile.createDimension('lon',nlon)

      lat_var = ncfile.createVariable('lats',np.float32,('lat','lon'))
      lat_var[:,:] = lats

      lon_var = ncfile.createVariable('lons',np.float32,('lat','lon'))
      lon_var[:,:] = lons

   lat_var.units = "degrees_north"
   lat_var.long_name = 'latitude'
   lon_var.units = "degrees_east"
   lon_var.long_name = "longitude"
   
   #Write variables
   twind_var = ncfile.createVariable('twind',np.float32,(vertname,'lat','lon'))
   twind_var[:,:,:] = twind
   twind_var.units = "meters_per_second"
   twind_var.long_name = "tangential wind component"
  
   rwind_var = ncfile.createVariable('rwind',np.float32,(vertname,'lat','lon'))
   rwind_var[:,:,:] = rwind
   rwind_var.units = "meters_per_second"
   rwind_var.long_name = "radial wind component"

   tc_lons_var = ncfile.createVariable('tc_center_lons',np.float32,(vertname))
   tc_lons_var[:] = tc_lons
   tc_lons_var.units = "degrees_west"
   tc_lons_var.long_name = "TC center longitude"
      
   tc_lats_var = ncfile.createVariable('tc_center_lats',np.float32,(vertname))
   tc_lats_var[:] = tc_lats
   tc_lats_var.units = "degrees_north"
   tc_lats_var.long_name = "TC center latitude"

   #Close nc file
   ncfile.close()

