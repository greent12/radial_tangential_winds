MODULE ideal_angle

!This module contains the subroutine to calculate the 4D array holding ideal angles
USE read_test_data
USE kinds
USE constants
USE omp_lib

IMPLICIT NONE

CONTAINS

SUBROUTINE calc_ideal_angle(nlat,nlon,lats,lons,ideal_angle)
! Calculate ideal angle 4D array
! INPUTS
!  nlat : number of latitudes in grid, used to specify shape of lats/lons
!  nlon : number of longitudes in grid, used to specify shape of lats/lons
!  lats : 2D array (nlat,nlon) holding latitudes
!  lons : 2D array (nlat,nlon) holding longitudes
!
! OUTPUTS
!  ideal_angle : 4D array (nlat,nlon,nlat,nlon) of ideal angles for idealized vortex, values will be between -pi and pi

  INTEGER,INTENT(IN) :: nlat,nlon
  INTEGER :: ybi,xbi
  REAL(KIND=single),DIMENSION(nlat,nlon),INTENT(IN) :: lats,lons
  REAL(KIND=single),DIMENSION(nlat,nlon) :: ideal_angle_i
  REAL(KIND=single),DIMENSION(nlat,nlon,nlat,nlon),INTENT(OUT) :: ideal_angle
  double precision :: timei,timef

  WRITE(*,*) " ","IN CALC_IDEAL_ANGLE ROUTINE"
  !$omp parallel do shared(ideal_angle) private(ideal_angle_i)
  DO xbi = 1,nlon
    DO ybi=1,nlat
      ideal_angle_i=atan2((lats-lats(ybi,xbi)),(lons-lons(ybi,xbi))) + pi/2.
      WHERE(ideal_angle_i>pi)
        ideal_angle_i = ideal_angle_i - 2.*pi
      END WHERE
      ideal_angle(ybi,xbi,:,:) = ideal_angle_i
      ideal_angle(ybi,xbi,ybi,xbi) = fill_val
    ENDDO
  ENDDO
  !$omp end parallel do
  WRITE(*,*) " ","LEAVING CALC_IDEAL_ANGLE ROUTINE"

END SUBROUTINE

END MODULE
