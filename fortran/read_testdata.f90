module read_test_data

implicit none

contains

subroutine load_data(lats,lons,u,v)
integer :: i,j
real, dimension(150,150),intent(out) :: lats,lons,u,v


open(14,file="lats.txt",action="read")
open(15,file="lons.txt",action="read")
open(16,file="uwind.txt",action="read")
open(17,file="vwind.txt",action="read")

do i =1,150
   read(14,*) (lats(i,j), j=1,150)
   read(15,*) (lons(i,j), j=1,150)
   read(16,*) (u(i,j), j=1,150)
   read(17,*) (v(i,j), j=1,150)
enddo

close(14)
close(15)
close(16)
close(17)

end subroutine 

end module 
