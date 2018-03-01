!subroutine test_self_focusing(dksi,Max_Intensity)
!use a_variables

!real :: dksi, Max_Intensity

!********** проверка шага по z на соответствие с условием для самофокусировки
!********** проверка фазового набега за дифракционный шаг
!	do while (0.5*3.77*beam.r*Max_Intensity*dksi>pi/5)
!      dksi=dksi/2
!      write(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
!      write(*,'(a,e9.3)') 'Decrease step z, now skin = ', grid.dz
!      write(*,'(a,f7.2)') 'I/Icr=', (beam.r*Max_Intensity)
!      write(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
!    end do!
!
!end subroutine test_self_focusing


subroutine self_focusing(dksi,Max_Intensity,time_num)
use a_parameters
use a_variables
implicit none
integer(integer_kind) :: k,m,t,time_num
real(real_kind) :: dksi, Max_Intensity, tmp_v

tmp_v=0.5*beam.RcrAB*beam.r*dksi
 do t=0,time_num-1	
  do m=0,grid.num_x-1
   do k=0,grid.num_x-1
     beam.e(k,m,t)=beam.e(k,m,t)*cdexp(cmplx(0.0d0,tmp_v*(cdabs(beam.e(k,m,t))**2),complex_kind ))
   end do
  end do
 end do
  

 
!********** проверка шага по z на соответствие с условием для самофокусировки
!********** проверка фазового набега за дифракционный шаг

!!if (0.5*beam.RcrAB*beam.r*Max_Intensity*dksi>AdaptiveCriterium) DecreseStepSelfFocus=1
 
phaseKerr=0.5*beam.RcrAB*beam.r*Max_Intensity !*dz - набег фазы при самофокусировке
  
!   dksi=dksi/2
!   write(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
!   write(*,'(a,e9.3)') 'Decrease step z, now skin = ', grid.dz
!   write(*,'(a,f7.2)') 'I/Icr=', (beam.r*Max_Intensity)
!   write(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
! end do

!*********************************  
  	 !ЛИНЕЙНОЕ ПОГЛОЩЕНИЕ
if (beam.linear_abs==1) then 	 
 if (z_old>zk_aerosol ) then 
  !curr_ae_screen=curr_ae_screen+1
  
  do t=0,time_num-1	
   do m=0,grid.num_x-1
    do k=0,grid.num_x-1
      beam.e(k,m,t)=beam.e(k,m,t)*cdexp(cmplx(-2.0d0*pi*(media.abs_drop_size)**2*media.abs_drop_concentration*grid.dz*beam.abs_ld/2,0.0d0,complex_kind))
    end do
   end do
  end do
  
  !zk_aerosol=media.start_z_aerosol+curr_ae_screen*media.delta_z_aerosol
 end if
end if 
!*********************************

end subroutine self_focusing
