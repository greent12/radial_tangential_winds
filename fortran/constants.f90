MODULE constants

!Put all constant parameters in this module
!The fill value parameter needs to be positive for things to work correctly in the recenter.f90 script
USE kinds

IMPLICIT NONE
SAVE
REAL(KIND=double),PARAMETER :: pi=4.0_double*DATAN(1.0_double)
REAL(KIND=single),PARAMETER :: fill_val=99999.
INTEGER,PARAMETER :: fill_val_int=99999
REAL(KIND=single),PARAMETER :: RE=6371.
INTEGER,PARAMETER :: MAXVERTLEVS=100

END MODULE
