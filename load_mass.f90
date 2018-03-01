! загрузка массивов функций рассеяния для капель размером 2..15 мкм
subroutine load_scattering_mass

use a_parameters
use a_types
use a_variables

implicit none

!**************** локальные переменные ****************************
 type tscattering_mass
    integer :: maxScPoint ! всего точек в файле maxScPoint+1
    real :: maxScR, &     ! максимальное значение R (значение в последней точке) [м]
	        areaR         ! область в которой рассматриваем интерференцию, если ==0, 
			              ! то по первому минимуму
	real, pointer :: scat(:,:)
 end type tscattering_mass

 type(tscattering_mass), pointer :: scattering_mass (:)

 integer,&
  parameter :: dropR_min   = 2,  & ! 2 Область размера капель
               dropR_max   = 15, & ! 15
			   TEST_load   = 0     !если >0, рисуем поверхности интерференции в \ScatFunc

 integer(integer_kind):: j_loc,i_loc,m,k,d
 character(3) :: str_loc
 character(2) :: str_loc2


!**************** загружаем файлы рассеяния ****************************
   write(*,'(a)'), ' --- LOAD SCATTERING FILES -----------------------------'
   write(*,*)
   allocate(scattering_mass(dropR_min:dropR_max)) 
   do j_loc=dropR_min,dropR_max
  	    if (j_loc>100) then 
		  write (str_loc,'(i3)') j_loc 
	      write(*,'(i3,a)'), j_loc, '  load scattering_fortran'//str_loc//'.dat'
          open (1,file='ScatFunc\scattering_fortran'//str_loc//'.dat',STATUS='OLD')
        else
		  if (j_loc>=10) then
		   write (str_loc2,'(i2)') j_loc
	      else
	       write (str_loc2,'(a1,i1)') '0',j_loc
	      end if 
	      write(*,'(i3,a)'), j_loc, '  load scattering_fortran'//str_loc2//'.dat'
          open (1,file='ScatFunc\scattering_fortran'//str_loc2//'.dat',STATUS='OLD')
        end if 

 
	  read(1,*), scattering_mass(j_loc).maxScPoint
	  read(1,*), scattering_mass(j_loc).maxScR
	  read(1,*), scattering_mass(j_loc).areaR
	  
      allocate(scattering_mass(j_loc).scat(0:1,0:scattering_mass(j_loc).maxScPoint)) 
                
	  do m=0, scattering_mass(j_loc).maxScPoint
		!Координата по r
!	    scattering_mass(j_loc).scat(index_loc,0)=
!       scattering_mass(j_loc).maxScR*index_loc/scattering_mass(j_loc).maxScPoint  
		!амплитуда и фаза рассеяния
		read(1,*), scattering_mass(j_loc).scat(0,m),scattering_mass(j_loc).scat(1,m) 
      end do
      close(1)
!----------------- TEST -----------------
!	    write(*,'(a,i4)') 'Number of point = ' , scattering_mass(j_loc).maxScPoint
!	    write(*,*) 'R_max =' , scattering_mass(j_loc).maxScR
!	    write(*,*) 'areaR =' , scattering_mass(j_loc).areaR
!		m=m-1
!		write (*,*) 'm = ', m 
!		!write (*,*) 'r = ', scattering_mass(j_loc).scat(index_loc,0) 
!		write (*,*) '1 ', scattering_mass(j_loc).scat(0,m) 
!		write (*,*) '2 ', scattering_mass(j_loc).scat(1,m) 

   end do

!   write (*,*) 'h =', grid.h*beam.abs_a0, ' beam.abs_a0', beam.abs_a0

!**************** инициализация массивов рассеяния ****************************
!   write(*,'(a)'), ' -------------------------------------------------------'
   allocate(scatter(dropR_min:dropR_max)) 
   do j_loc=dropR_min,dropR_max
	 scatter(j_loc).r = 1.e-6*j_loc ! размер капли в метрах
	 
	 if (scattering_mass(j_loc).areaR<=0) then
	   i_loc=0
	   do while (scattering_mass(j_loc).scat(0,i_loc)>scattering_mass(j_loc).scat(0,i_loc+1))
	     i_loc=i_loc+1
	     if (i_loc==scattering_mass(j_loc).maxScPoint) exit
	   end do
	   ! область интерференции в шагах сетки
	   scatter(j_loc).area=nint(scattering_mass(j_loc).maxScR*i_loc/(scattering_mass(j_loc).maxScPoint*grid.h*beam.abs_a0))       
     else
       ! область интерференции в шагах сетки
	   scatter(j_loc).area=nint(scattering_mass(j_loc).areaR/(grid.h*beam.abs_a0))
	 end if
	 if (scatter(j_loc).area>grid.num_x/10) scatter(j_loc).area=grid.num_x/10 
	 
	 d=scatter(j_loc).area
	 allocate(scatter(j_loc).e(-d:d,-d:d)) 
     do m=-d,d
       do k=-d,d
         if ( (k**2+m**2)<=d**2 ) then
		   i_loc=nint( grid.h*beam.abs_a0/scattering_mass(j_loc).maxScR*sqrt(1.0*(k**2+m**2))*scattering_mass(j_loc).maxScPoint)
		   scatter(j_loc).e(k,m)=scattering_mass(j_loc).scat(0,i_loc)*cdexp(i*scattering_mass(j_loc).scat(1,i_loc))
		 else
		   scatter(j_loc).e(k,m)=0
		 end if
	   end do
	 end do
   
!     write(*,'(a,i3,a)'), ' droplet ', j_loc,' uploaded' 
!----------------- TEST -----------------
	  write (*,*) 'droplet',j_loc, 'area',scatter(j_loc).area
!	  write (*,*) 'r0', scatter(j_loc).area*grid.h*beam.abs_a0
  if (TEST_load>0) then
	 if (j_loc<10) then 
	   write (str_loc,'(a1,i1)') '0',j_loc
	 else
  	    if (j_loc>100) then 
		  write (str_loc,'(i3)') j_loc 
        else
		  write (str_loc,'(i2)') j_loc
        end if 
	 end if
     open (1,file='ScatFunc\scattering_fortran'//str_loc//'.grd')
     write (1,110)
     write (1,111) 2*d+1, 2*d+1
     write (1,112) -d*grid.h, d*grid.h !от чего до чего по оси Х
     write (1,112) -d*grid.h, d*grid.h !от чего до чего по оси Y
     write (1,112) 0.9, 1.1 !от чего до чего по оси Z
 
     do m=-d,d
	     write(1,113) (cdabs(scatter(j_loc).e(k,m) +1.0*cdexp(i*0))**2,k=-d,d)
	 end do
     close(1)
  end if
!-----------------------------------------
   end do

   do j_loc=dropR_min,dropR_max
     deallocate(scattering_mass(j_loc).scat) 
   end do
   deallocate(scattering_mass) 

   write(*,'(a)'), ' -------------------------------------------------------'
   write(*,*) 
 
if (TEST_load>0) call exit

110    FORMAT('DSAA')
111    FORMAT(2I8)
112    FORMAT(2F15.7)
113    FORMAT(35F15.7)

end subroutine load_scattering_mass