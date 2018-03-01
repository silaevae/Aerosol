subroutine turbulence(Nscreen,time_num)

use a_parameters
use a_variables


    implicit none

    ! Variables
    
    real(4):: B(1:5), C(1:8)
    integer(4) :: len, m1, m2, nwork, isp, iex, m, k, t, time_num, Nscreen
    real(4):: amin, amax, r0, disp, am2, am
    real(4), allocatable:: phi(:,:)
    character(len=255):: str
    
    allocate (phi(0:grid.num_x-1, 0:grid.num_x-1))
    
    if (Nscreen<10)then
      write(str,'(a19,i1,a3)'),'TurbScreens\1a16000',Nscreen,'.fi' !Преобразуем число в строку str
    else
      write(str,'(a18,i2,a3)'),'TurbScreens\1a1600',Nscreen,'.fi'
    end if
    
    open (NTurb,file=str, status='OLD', form='BINARY')
    
    read(NTurb) len
    read(NTurb) m1
    read(NTurb) m2
    read(NTurb) B(:)
    read(NTurb) isp
    read(NTurb) nwork
    read(NTurb) r0
    read(NTurb) iex
    read(NTurb) C(:)
    read(NTurb) phi(:,:)
             
 do t=0,time_num-1	
  do m=0,grid.num_x-1
   do k=0,grid.num_x-1
   
     beam.e(k,m,t)=beam.e(k,m,t)*cdexp(cmplx(0.0d0,phi(k, m),complex_kind))
     
   end do
  end do
 end do
    

   close(NTurb)
    

    
   deallocate (phi)

end subroutine turbulence
