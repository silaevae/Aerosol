module a_dispersion

!use a_parameters
 use a_variables
 use MKL_DFTI

 implicit none

!************************************************************
 contains

!**************************************************************************
!                             FFTDispersion
!
!           Calculates 2nd and 3d order dispersion of the pulse 
!                     using fast fourier transform
!**************************************************************************
subroutine do_dispersion

   integer(4),automatic            :: k,n,j
   real(real_kind),automatic       :: w,dw
   complex(complex_kind),automatic :: factor2, factor3

   type(DFTI_DESCRIPTOR), pointer  :: Dfti_Handel
   integer(4) :: Dfti_Status

     dw      = 2.*pi/grid.time_number/grid.dt
     factor2 = 1/2.*i*beam.abs_ld/beam.abs_lds*grid.dz
     !factor3 = -1/3.*i/Pulse.Ldisp3*grid.dz*0.
          
   
    !$OMP PARALLEL PRIVATE(Dfti_Handel, Dfti_Status, w)   
     Dfti_Status = DftiCreateDescriptor(Dfti_Handel, DFTI_DOUBLE, DFTI_COMPLEX, 1, grid.time_number)
     Dfti_Status = DftiSetValue(Dfti_Handel, DFTI_BACKWARD_SCALE, 1./(grid.time_number))
     Dfti_Status = DftiCommitDescriptor(Dfti_Handel)
     
     !$OMP DO
     do n=0,grid.num_x-1
      do k=0,grid.num_x-1
         
        Dfti_Status = DftiComputeForward(Dfti_Handel, beam.e(k,n,:))     

!       Multiplying by phase factor
        do j=0,grid.time_number/2-1
          w=DFLOTJ(j)*dw
          beam.e(k,n,j)=beam.e(k,n,j)*zexp(factor2*w**2) !+factor3*w**3)
          
        end do
        do j=grid.time_number/2,grid.time_number-1
          w=DFLOTJ(j-(grid.time_number))*dw
          beam.e(k,n,j)=beam.e(k,n,j)*zexp(factor2*w**2) !+factor3*w**3)
          
        end do
 
       
        Dfti_Status = DftiComputeBackward(Dfti_Handel, beam.e(k,n,:))
         
      end do
     end do
    !$OMP END DO
     
     Dfti_Status = DftiFreeDescriptor(Dfti_Handel)
    
    !$OMP END PARALLEL   
   
end subroutine do_dispersion

end module a_dispersion