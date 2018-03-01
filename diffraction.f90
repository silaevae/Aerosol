module a_difraction

!use a_parameters
 use a_variables
 use MKL_DFTI

 implicit none

!************************************************************
 contains

!************************************************************
 subroutine do_difraction(delta_z,num_x, num_t0, num_t1)
!! дифрагирует поле field (указатель на двумерный массив num_x:num_x) на расстояние delta_z

   integer(4) n,k,num_2,num_x, tn, num_t0, num_t1
   real(real_kind) delta_z, k2
   real(real_kind) time_2,time_1
   complex(complex_kind), allocatable :: field(:)
   integer(4):: L(1:2)
   type(DFTI_DESCRIPTOR),pointer :: Dfti_Handle 
   integer(4):: Dfti_Status
   
   
   L(1)=num_x
   L(2)=num_x

!*** c промежуточным массивом f
!  complex, allocatable :: f(:,:)
!   integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM



!time_2=time_1
call cpu_time(time_1)
!write(*,'(a,F15.10,a)'), ' Time of init_difraction is ',time_1 - time_2,' seconds'

!$OMP PARALLEL PRIVATE (time_1, time_2, Dfti_Handle, Dfti_Status, n,k,num_2,k2, field) default(shared)
allocate (field(0:num_x*num_x-1))
!allocate(f(0:num_x-1,0:num_x-1))

 
     Dfti_Status = DftiCreateDescriptor(Dfti_Handle, DFTI_DOUBLE, DFTI_COMPLEX, 2, L)
     Dfti_Status = DftiSetValue(Dfti_Handle, DFTI_BACKWARD_SCALE, 1./(L(1)*L(2)))
     Dfti_Status = DftiCommitDescriptor(Dfti_Handle)

   
!$OMP DO SCHEDULE (DYNAMIC, 4)

   do tn=num_t0,num_t1-1
     !print*, 't_index=',tn,', proc', OMP_GET_THREAD_NUM(),'of', OMP_GET_NUM_THREADS()
     
    !do k=0,num_x-1
    !  do n=0,num_x-1
    !    f(n,k)=beam.e(n,k,tn)
    !  end do
    ! end do

!time_2=time_1
!call cpu_time(time_1)
!write(*,'(a,F15.10,a)'), ' Time of f(n,k)=field(n,k,tn) is ',time_1 - time_2,' seconds'


!*********  делаем бпф

     field = pack (beam.e(:,:,tn), .true.) !преобразует поле е в одномерный массив
     
     Dfti_Status = DftiComputeForward(Dfti_Handle, field)!делаем прямое преобразование Фурье для массива field
     
     beam.e(:,:,tn) = reshape (field, L) !преобразует field обратно в двумерный массив
     
     !call cfft2d (beam.e(0,0,tn), num_x, num_x, -1 )
     !call cfft2d (f(0,0), num_x, num_x, -1 )
       !  получаем спектр, у которого абсолютный шаг по частоте 2*pi/grid.abs_size_x
       !  tmp_val=-i*delta_z*beam.abs_ld*4*(pi**2)/(2*beam.abs_k0*(grid.abs_size_x**2))
 !time_2=time_1
!call cpu_time(time_1)
!write(*,'(a,F15.10,a)'), ' Time of cfft2d is ',time_1 - time_2,' seconds'

     num_2=num_x/2

   
!********  добавляем фазовые набеги в приближении kx,ky<<k0


     do k=0,num_2
      k2=k**2
	  do n=0,num_2
      beam.e(n,k,tn)=beam.e(n,k,tn)*cdexp(cmplx(0.,-delta_z*(n**2+k2),complex_kind))
      end do
      do n=num_2+1,num_x-1
       beam.e(n,k,tn)=beam.e(n,k,tn)*cdexp(cmplx(0.,-delta_z*((n-num_x)**2+k2),complex_kind))
      end do
     end do
   
     do k=num_2+1,num_x-1
      k2=(k-num_x)**2
	  do n=0,num_2
       beam.e(n,k,tn)=beam.e(n,k,tn)*cdexp(cmplx(0.,-delta_z*(n**2+k2), complex_kind))
      end do
      do n=num_2+1,num_x-1
       beam.e(n,k,tn)=beam.e(n,k,tn)*cdexp(cmplx(0.,-delta_z*((n-num_x)**2+k2),complex_kind))
      end do
     end do
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !do k=0,num_2
     ! k2=k**2
	  !do n=0,num_2
      ! f(n,k)=f(n,k)*cexp(cmplx(0.,-delta_z*(n**2+k2)))
      !end do
      !do n=num_2+1,num_x-1
      ! f(n,k)=f(n,k)*cexp(cmplx(0.,-delta_z*((n-num_x)**2+k2)))
      !end do
     !end do
   
     !do k=num_2+1,num_x-1
      !k2=(k-num_x)**2
	  !do n=0,num_2
      ! f(n,k)=f(n,k)*cexp(cmplx(0.,-delta_z*(n**2+k2)))
      !end do
      !do n=num_2+1,num_x-1
       !f(n,k)=f(n,k)*cexp(cmplx(0.,-delta_z*((n-num_x)**2+k2)))
      !end do
     !end do
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!time_2=time_1
!call cpu_time(time_1)
!write(*,'(a,F15.10,a)'), ' Time of phase increase is ',time_1 - time_2,' seconds'

!*********  обратное бпф
     field = pack (beam.e(:,:,tn), .true.) !преобразует поле е в одномерный массив
     
     Dfti_Status = DftiComputeBackward(Dfti_Handle, field)!делаем обратное преобразование Фурье для массива field
     
     beam.e(:,:,tn) = reshape (field, L) !преобразует field обратно в двумерный массив

  !call cfft2d (beam.e(0,0,tn), num_x, num_x, 1)
   
  !call cfft2d (f(0,0), num_x, num_x, 1)

!time_2=time_1
!call cpu_time(time_1)
!write(*,'(a,F15.10,a)'), ' Time of cfft2d -1 is ',time_1 - time_2,' seconds'


 !    do k=0,num_x-1
 !     do n=0,num_x-1
 !      beam.e(n,k,tn)=f(n,k)
 !     end do
 !    end do

   end do

!$OMP END DO  

!deallocate(f)
deallocate (field)
Dfti_Status = DftiFreeDescriptor(Dfti_Handle)
  
!$OMP END PARALLEL

time_2=time_1
call cpu_time(time_1)
write(*,'(a,F15.10,a)'), ' Time of difraction is ',time_1 - time_2,' seconds'
  
 end subroutine do_difraction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! subroutine do_difraction_thread(arg)
!   type(tArgument) :: arg
!	 
!   call do_difraction(arg.delta_z_tmp,arg.grid.num_x,arg.t0,arg.t1)
!   
! end subroutine do_difraction_thread(arg)


end module a_difraction


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! без промежуточного массива не работает, если beam.e указатель 
!!!! и размер выделенной памяти больше ~256x256
!! subroutine do_difraction(delta_z)
!!   integer(integer_kind) n,k,num_2
!!   real(real_kind) delta_z
!!   complex(complex_kind) tmp_val
!!
!!   call cfft2d ( beam.e, grid.num_x, grid.num_x, -1 )
!!!  делаем бпф
!!!  получаем спектр, у которого абсолютный шаг по частоте 2*pi/grid.abs_size_x
!!   tmp_val=-i*delta_z*beam.abs_ld*4*(pi**2)/(2*beam.abs_k0*(grid.abs_size_x**2))
!!
!!   num_2=grid.num_x/2
!!   
!!   do n=0,num_2
!!    do k=0,num_2
!!     beam.e(n,k)=beam.e(n,k)*exp(tmp_val*(n**2+k**2))
!!    end do
!!    do k=num_2+1,grid.num_x-1
!!     beam.e(n,k)=beam.e(n,k)*exp(tmp_val*(n**2+(k-grid.num_x)**2))
!!    end do
!!   end do
!! 
!!   do n=num_2+1,grid.num_x-1
!!    do k=0,num_2
!!     beam.e(n,k)=beam.e(n,k)*exp(tmp_val*((n-grid.num_x)**2+k**2))
!!    end do
!!    do k=num_2+1,grid.num_x-1
!!     beam.e(n,k)=beam.e(n,k)*exp(tmp_val*((n-grid.num_x)**2+(k-grid.num_x)**2))
!!    end do
!!   end do
!!
!!   call cfft2d ( beam.e, grid.num_x, grid.num_x, 1 )
!!
!! end subroutine do_difraction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! end module a_difraction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!  для тестирования   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   complex, allocatable :: f1d(:)
!   complex, allocatable :: f1d_temp(:)
!   allocate(f1d(0:grid.num_x-1))
!   allocate(f1d_temp(0:2*(grid.num_x-1)))
!   do k=0,grid.num_x-1
!     f1d(k)=beam.e(512,k+1)
!   end do
! 
!  open(NIntFNum,FILE=IntGraphName//'signal_test_'//IntFName,access='sequential',status='unknown')  
!  do k=0,grid.num_x-1
!   write(NIntFNum,*),(k-grid.num_x/2+1)*grid.h,' ',cabs(f1d(k))
!  end do
!  close(NIntFNum)
!
!  call cfft1d ( f1d, grid.num_x, 0, f1d_temp)  
!  call cfft1d ( f1d, grid.num_x, 1, f1d_temp)
!    
!  open(NIntFNum,FILE=IntGraphName//'spectr_test_'//IntFName,access='sequential',status='unknown')  
!  do k=0,grid.num_x-1
!    write(NIntFNum,*),(k-grid.num_x/2+1)/grid.size_x,' ',cabs(f1d(k))
!  end do
!  close(NIntFNum)  
!  deallocate(f1d)  
!  deallocate(f1d_temp)  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
