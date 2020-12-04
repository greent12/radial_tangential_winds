import numpy as np

def uv_to_rt(lons,lats,uwind,vwind,lon0,lat0,deg=True):
   """
   Convert u and v winds to radial (r) and tangential (t) winds.
   """
   #Convert all things that may be in degrees to radians for numpy
   if deg:
      lon0=np.radians(lon0)
      lat0=np.radians(lat0)
      lons=np.radians(lons)
      lats=np.radians(lats)

   X=np.cos(lats)*np.sin(lons-lon0)
   Y=np.cos(lat0)*np.sin(lats) - np.sin(lat0)*np.cos(lats)*np.cos(lons-lon0)
   theta=np.arctan2(Y,X)
   rwind=np.add(uwind*np.cos(theta),vwind*np.sin(theta))
   twind=np.subtract(vwind*np.cos(theta),uwind*np.sin(theta))

   return rwind,twind
