module a_functions

use a_types
use a_variables
use a_parameters
use dflib

implicit none 

contains
 
!subroutine calculateFluence()
 !  use a_parameters
!   use a_variables
!   
!   integer(4) :: n,kkk,k,nnn,t_index
!
!   if(.not. associated(beam.fluence)) then
!      allocate(beam.fluence(0:grid.num_x-1,0:grid.num_x-1))
!   end if
!
!   nnn=0
! !do n=grid.num_x/2-num_2,grid.num_x/2+num_2,sdvig
!   do n=0,grid.num_x-1
!   kkk=0
!   !do k=grid.num_x/2-num_2,grid.num_x/2+num_2,sdvig
!   do k=0,grid.num_x-1
!     beam.fluence(kkk,nnn)=0.0
!     do t_index=0,grid.time_number-1	
!       beam.fluence(kkk,nnn)=beam.fluence(kkk,nnn)+grid.delta_time(t_index)*(cabs(beam.e(k,n,t_index))**2)
!     end do
!     kkk=kkk+1
!  end do
! nnn=nnn+1
! end do
 
!end subroutine calculateFluence
 
!**************************************************************************
!Вычисляет случайную величину по распределению пуассона с параметром = a_param
!**************************************************************************
 function puasson_random(a_param)
  real(real_kind) a_param
  integer puasson_random
  real(8) calc_value, rand_value
  puasson_random=0
  if (a_param>100) then
	 puasson_random=nint(random_normal(a_param,sqrt(a_param)))
  else  
   calc_value=exp(-a_param)
   call random_number(rand_value)
   do while (rand_value>calc_value)
	rand_value=rand_value-calc_value   
	puasson_random=puasson_random+1
    calc_value=calc_value*a_param/puasson_random
   end do
  end if
   
 end function puasson_random

!**************************************************************************
!
!**************************************************************************
 function calc_energy(field,num,num_t,rad,max_value,n_point,k_point,tt_point)
  complex(complex_kind),pointer :: field(:,:,:) 
  integer(integer_kind) num, num_t, n_point, k_point, tt_point
  real(real_kind) calc_energy, max_value, rad
  real(real_kind) n_r, num_2
  real(real_kind) tmp_energy, tmp_int
  integer(integer_kind) n, k, nt, k2 
   
  n_r=(rad*num)**2
  
  num_2=num/2-0.5
  calc_energy=0
  max_value=0
  tt_point=0
  n_point=0
  k_point=0
  do nt=0,num_t-1
   do k=0,num-1
	!k2=(k-num_2)**2
	do n=0,num-1
	 if (((k-num_2)**2+(n-num_2)**2)<=n_r) then
	 tmp_int=(cdabs(field(n,k,nt)))**2 
	 tmp_energy=tmp_int*grid.delta_time(nt)
	 calc_energy=calc_energy+tmp_energy
	 if (tmp_int>max_value) then 
	   max_value=tmp_int
	   n_point=n
	   k_point=k
	   tt_point=nt
	 end if
	 end if  
	end do
   end do
  end do

 end function calc_energy
 
 
 
 !***********************************************************************
 !                     Integral Electron Density
!***********************************************************************/
function calc_plasma(density,num,rad,max_value)
  real(real_kind),pointer :: density(:,:) 
  integer(integer_kind) num
  real(real_kind) calc_plasma, rad, max_value
  real(real_kind) tmp_plasma
  integer(integer_kind) n, k, nt, n_r, k2, num_2 
   
  n_r=nint((rad*num)**2)
  
  num_2=num/2
  calc_plasma=0
  max_value=0
  
   do k=0,num-1
	!k2=(k-num_2)**2
	do n=0,num-1
	 !if ((k2+(n-num_2)**2)<=n_r) then 
	 tmp_plasma=density(n,k)
	 calc_plasma=calc_plasma+tmp_plasma
	end do
   end do
  
  max_value=maxval(density)  !Максимальное значение концентрации электронов
 end function calc_plasma
!***********************************************************************
 !                     Normal distribution
!***********************************************************************/
  
 function random_normal(m,s)
  ! normal distribution with mean m and standard deviation s
    real(real_kind) :: random_normal
	real(real_kind) :: m,s
	real(real_kind) :: x1, x2, w
	    
  do 
    call random_number(x1)
	x1 = 2. * x1 - 1.
    call random_number(x2)
	x2 = 2. * x2 - 1.
    w = x1**2 + x2**2
    if (w > 0.0  .and. w < 1.0) exit
  end do
  w = sqrt(-2.*log(w)/w)
  x1 =x1* w
  ! x2 *= w;  // a second normally distributed result not used
  !!Now make the Box-Muller transformation to
  !!get two normal deviates. Return one and
  !!save the other for next time
  random_normal= x1 * s + m
 end function random_normal

 
!**************************************************************************
!                        ReturnDate(date)
!
!
!**************************************************************************
subroutine ReturnDate(date)

  character(len=*) :: date
  character(len=2) :: day,mon
  character(len=4) :: yr
  integer(2) :: iyr,imon,iday

  call getdat(iyr,imon,iday)

  write(yr,'(I4)') iyr

  write(mon,'(I2)') imon
  if(imon<10) mon(1:1)='0'
  
  write(day,'(I2)') iday
  if(iday<10) day(1:1)='0'

  date=day//'.'//mon//'.'//yr

end subroutine ReturnDate

!**************************************************************************
!                         ReturnTime(time)
!
!
!**************************************************************************
subroutine ReturnTime(time)

 character(len=*) :: time
 integer(2) :: ihr,imin,isec,i100th
 character(len=2) :: chr,cmin,csec,c100th

 call gettim(ihr,imin,isec,i100th)

 write(chr,'(I2)') ihr
 if(ihr<10) chr(1:1)='0'

 write(cmin,'(I2)') imin
 if(imin<10) cmin(1:1)='0'

 write(csec,'(I2)') isec
 if(isec<10) csec(1:1)='0'

 write(c100th,'(I2)') i100th
 if(i100th<10) c100th(1:1)='0'

 time=chr//':'//cmin//':'//csec//':'//c100th

end subroutine ReturnTime





end module a_functions


