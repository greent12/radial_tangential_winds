program recenter

use read_test_data  !Has subroutine that reads in test data

implicit none
integer, parameter :: GRIDDIM=150
real, dimension(GRIDDIM,GRIDDIM) :: lons,lats,uwind,vwind

!Read in a simple case 150 x 150 grids from Hurricane Matthew 19z analysis HWRF run
call load_data(lats,lons,uwind,vwind)



end program
