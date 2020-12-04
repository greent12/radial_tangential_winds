#This module holds a dictionary of vertical level type pairs
# The keys are the vertical level types that the user would pick in namelist.input and the values are what they correspond to in grib files
#Add to this dictionary if needed

vertical_dict={"msl"                 : "meanSea",
               "pressure"            : "isobaricInhPa",
               "height_above_ground" : "heightAboveGround",
               "height_above_sea"    : "heightAboveSea"
               }
