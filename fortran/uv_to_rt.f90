MODULE uv_to_rt

!Module holding routine to convert uv winds to tangential/radial winds
USE constants
USE kinds

IMPLICIT NONE

CONTAINS

SUBROUTINE calc_uv_to_rt(nlon,nlat,lons,lats,lon0,lat0,uwind,vwind,rwind,twind)

!INPUTS
! nlon : number of longitudes in grid
! nlat : number of latitudes in grid
! lons : 2d grid (nlat,nlon) of longitudes on the grid
! lats : 2d grid (nlat,nlon) of latitudes on the gird
! lon0 : longitude of center 
! lat0 : latitude of center
! uwind : 2d grid (nlat,nlon) of zonal wind
! vwind : 2d grid (nlat,nlon) of meridional wind
!
!OUTPUTS
! rwind : 2d grid (nlat,nlon) of radial winds
! twind : 2d grid (nlat,nlon) of tangential winds

INTEGER,INTENT(IN) :: nlon,nlat
REAL(KIND=single),DIMENSION(nlat,nlon),INTENT(IN) :: lons,lats,uwind,vwind
REAL(KIND=single),INTENT(IN) :: lon0,lat0
REAL(KIND=single),DIMENSION(nlat,nlon),INTENT(OUT) :: rwind,twind

REAL(KIND=single) :: lon0_rad,lat0_rad
REAL(KIND=single),DIMENSION(nlat,nlon) :: lons_rad,lats_rad,X,Y,theta

lon0_rad = lon0*PI/180.
lat0_rad = lat0*PI/180.
lons_rad = lons*PI/180.
lats_rad = lats*PI/180.

X=COS(lats_rad)*SIN(lons_rad-lon0_rad)
Y=COS(lat0_rad)*SIN(lats_rad) - SIN(lat0_rad)*COS(lats_rad)*COS(lons_rad-lon0_rad)
theta = ATAN2(Y,X)

rwind = uwind*COS(theta) + vwind*SIN(theta)
twind = vwind*cos(theta) - uwind*sin(theta)

END SUBROUTINE

END MODULE
