MODULE calc_distance_from_point

USE kinds
USE constants

IMPLICIT NONE

CONTAINS

  SUBROUTINE calc_distance_g2p(lonpoint,latpoint,nlon,nlat,lons,lats,dist)
  REAL(KIND=single),INTENT(IN) :: lonpoint,latpoint
  REAL(KIND=single),DIMENSION(nlat,nlon),INTENT(IN) :: lons,lats
  INTEGER,INTENT(IN) :: nlon,nlat
  REAL(KIND=single),DIMENSION(nlat,nlon),INTENT(OUT) ::dist
  REAL(KIND=single),DIMENSION(nlat,nlon) ::dlon,dlat,a,lon2_rad,lat2_rad
  REAL(KIND=single) :: lon1_rad,lat1_rad

  !Initialize distance array
  dist=0.0
 
  !Convert to radians
  lon1_rad = lonpoint*PI/180.
  lat1_rad = latpoint*PI/180.
  lon2_rad = lons*PI/180.
  lat2_rad = lats*PI/180.

  !Calculate distance using Haversine formula
  dlon = lon2_rad - lon1_rad
  dlat = lat2_rad - lat1_rad
  a = SIN(dlat/2.)**2 + COS(lat1_rad)*COS(lat2_rad)*SIN(dlon/2.)**2
  dist = 2.*ASIN(SQRT(a))*RE
 
  END SUBROUTINE

END MODULE
