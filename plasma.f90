module a_plasma
!************************************************************
! ���� �� ������ M_Physics �������� �.
!
! Types:
!
! Contains:
!
!   subroutine TabFuncCreation
!--   subroutine DifractionX(j)
!--   subroutine DifractionY(j)
!--   subroutine ElectronKerr(j)
!--   subroutine VKRKerr(j)
!   subroutine Plasma(j)
!   subroutine EnergyAbsorption(j)
!************************************************************

 use a_medium




 implicit none

 !���������� ���������� ������� ��������� ���������
 integer(4),parameter :: maxRair     =9999 ,& !���������� ��������� � ���������� ������� ��� ��������� ��������� � �������
                         maxRmethanol=10000   !���������� ��������� � ���������� ������� ��� ��������� ��������� � ��������

 real(real_kind) :: IImin ,& !����������� ������������� � ���������� �������
                    IImax    !������������ ������������� � ���������� �������

 real(real_kind) :: Rmethanol(1:maxRmethanol) ,& !��������� ��������� ��������
                    RN(1:maxRair)             ,& !��������� ��������� �����
                    RO(1:maxRair)                !��������� ��������� ���������

 real(real_kind),pointer           :: integM(:,:) !��������������� ������ ��� ���������� ���
 complex(complex_kind),pointer :: VKRfunc(:)  !������� ������� ���
 complex(complex_kind),pointer :: FFTcoef(:)  !������������ ���
 
 real(real_kind) :: TimeCPU_Difraction   ,&  !������������ ����� ������������� �� ������ ��������� ���������
                    TimeCPU_ElectronKerr ,&  !������������ ����� ������������� �� ������ ��������� ������������ ������� �����
                    TimeCPU_VKRKerr      ,&  !������������ ����� ������������� �� ������ ��������� ������������ ��������������� ���������
                    TimeCPU_Plasma       ,&  !������������ ����� ������������� �� ������ ��������� ���������� ������������
                    TimeCPU_EnergyAbsorption !������������ ����� ������������� �� ������ ��������� ���������� ������� ��� ���������

!**************************************************************************
                              CONTAINS
!**************************************************************************
!**************************************************************************
!                           TabFuncCreation
!
!  ��������� �� ����� � ������ ���������� ������� ��� ��������� ��������� �
!  ������� ������� ��� ���
!**************************************************************************
 subroutine TabFuncCreation

  integer(4)      :: unit,k,j
  real(real_kind) :: Om,Gam,Lam,t

  if(TimeFlag /= 0 .and. PlasmaFlag/=0) then

     unit=NIntFNumPlasma

     SelectCase(MediumFlag)
     Case(1) !===   ������   ================================================

        open (unit,file='TabFunc\N2.dat') !,recordtype='stream_cr')
        read (unit,*) RN(:)
        close(unit)
 
        open (unit,file='TabFunc\O2.dat') !,recordtype='stream_cr')
        read (unit,*) RO(:)
        close(unit)

        !��������� t0/21.e-15 �������� �������������� �������� ���������,
        !����������� � ����� ��� ������ t0=21 ��, �� ������ ������������
        !������������ ��������, ���������� ���������� t0
        RN=RN*beam.pulse_duration/21.d-15
        RO=RO*beam.pulse_duration/21.d-15 

        IImin=1.d15/beam.abs_i0 !��������� ��������� �������� ��� ��������������, ������� � I=1.e15 ��/�**2
        
          
        IImax=maxRair*IImin

!write (*,*) 'IImin', IImin
!write (*,*) 'IImax', IImax

     Case(2) !===   �������   ===============================================

        open (unit,file='TabFunc\Methanol.dat') !,recordtype='stream_cr')
        read (unit,*) Rmethanol(:)
        close(unit)

        Rmethanol=Rmethanol*beam.pulse_duration !����������

        IImin=1.d14/beam.abs_i0 !��������� ��������� �������� ��� ��������������, ������� � I=1.e14 ��/�**2
        IImax=maxRmethanol*IImin

     End Select !============================================================

  end if

!  !===   ������� ������� ��� ���   =========================================
!  if(TimeFlag/=0 .and. VKRKerrFlag/=0) then
!
!     Om=20.6e12  !����������� ������� [��]
!     Gam=26.0e12 !����������� ������� [��]
!     Lam=sqrt(Om*Om-Gam*Gam/4.)
!
!     !������ ������� ������� ���
!     t=0
!     do k=0,mfur
!        t=k*grid.delta_time
!        VKRfunc(k)=cmplx( Om*Om*exp(-0.5*Gam*beam.pulse_duration*t)*sin(Lam*beam.pulse_duration*t)/Lam*beam.pulse_duration, 0. )
!     end do
!
!     !��������� ������ ������� ������� ���
!     call cfft1d(VKRfunc(:),mfur+1,0,FFTcoef(:))
!     call cfft1d(VKRfunc(:),mfur+1,1,FFTcoef(:))
!
!  end if

 end subroutine TabFuncCreation

!**************************************************************************
!!                         VKRKerr(istart,iend)
!!
!!                      ��� �� ������������ ���������
!!
!! istart - ����� ������� ����� �������
!! iend - ������ ������� ����� �������
!!**************************************************************************
! subroutine VKRKerr(istart,iend)
!
!  integer(4) :: istart,iend
!  integer(4),automatic :: k,n,j
!  real(real_kind),automatic :: tux
!  complex(complex_kind),automatic :: IntDum(0:mfur),dn(0:Grid.Nt),dnDum(0:mfur)
!
!  do n=istart,iend
!  do k=0,Grid.Nx
!
!     !�������� ������������� �������� � ������ ������� �����
!     IntDum(:)=0.
!     do j=0,Grid.Nt
!        IntDum((mfur+1)/2-(Grid.Nt+1)/2+j)=cmplx(cabs(beam.e(k,n,j))**2,0.)
!     end do
!
!     !��������� ������ �������������
!     call cfft1d(IntDum(:),mfur+1,0,FFTcoef(:))
!     call cfft1d(IntDum(:),mfur+1,1,FFTcoef(:))
!      
!     !����������� ������ ������������� �� �������� ������� ������� ��� ��� ���������� ������� ������ ����
!     do j=0,mfur
!        dnDum(j)=IntDum(j)*VKRfunc(j)
!     end do
!
!     !�� ������� ��������������� ����� ����
!     call cfft1d(dnDum(:),mfur+1,0,FFTcoef(:))
!     call cfft1d(dnDum(:),mfur+1,-1,FFTcoef(:))
!
!     !�� ������� ������� ����� ����� ������ ��� �����
!     do j=0,Grid.Nt
!        dn(j)=dnDum((mfur+1)/2-(Grid.Nt+1)/2+j)
!        dn(j)=dn(j)*mfur*grid.delta_time !�� ������� ������, �� ��� ���������� ����������� ���� �������
!     end do
!
!     do j=0,Grid.Nt
!
!        tux=-0.5*Medium.g*Medium.R*real(dn(j))
!        beam.e(k,n,j)=beam.e(k,n,j)*cexp(cmplx(0.,tux*Grid.dz))
!
!        !���� ����� ���������� ���� �� ������� ��� �� z ��������
!        !AdaptiveCriterium, �� �� ��������� ���� dz ���������� � 2 ����
!        if(abs(tux*Grid.dz) >= AdaptiveCriterium) adapt2=1
!
!     end do
!
!  end do
!  end do
!
! end subroutine VKRKerr



!!**************************************************************************
!!                       ElectronKerr(istart,iend,j)
!!
!!                       ����������� ������ �����
!!
!! istart - ����� ������� ����� �������
!! iend - ������ ������� ����� �������
!! j - ����� ���������� ���� ��������
!!**************************************************************************
! subroutine ElectronKerr(istart,iend,j)
!
!  integer(4) :: j,istart,iend
!  integer(4),automatic :: k,n
!  real(real_kind),automatic :: tux
!
!  IF(ElectronKerrFlag /=0 .and. Medium.g /=1.) THEN
!
!     do n=istart,iend
!        do k=0,Grid.Nx
!
!           tux=-0.5*(1.-Medium.g)*Medium.R*cabs(beam.e(k,n,j))**2
!           beam.e(k,n,j)=beam.e(k,n,j)*cexp(cmplx(0.,tux*Grid.dz))
!
!           !���� ����� ���������� ���� �� ������� ��� �� z ��������
!           !AdaptiveCriterium, �� �� ��������� ���� dz ���������� � 2 ����
!           if(abs(tux*Grid.dz) >= AdaptiveCriterium) adapt1=1
!
!        end do
!     end do
!
!  END IF
!
! end subroutine ElectronKerr

!**************************************************************************
!                         Plasma(istart,iend,j)
!
!                       ���������� ������������
!
! istart - ����� ������� ����� �������
! iend - ������ ������� ����� �������
! tn - ����� ���������� ���� ��������
!**************************************************************************
 subroutine Plasma(istart,iend,tn)

  integer(integer_kind) :: tn,istart,iend
  integer(integer_kind),automatic :: k,n,m
  real(real_kind),automatic :: Int,Rm,RmN,RmO,tux,tux2 !,absE
  real(real_kind),automatic :: tmp_rm_r
  integer :: tmp_rm_i
  
  IF(PlasmaFlag /= 0) THEN

     !��������� ������ �������� ��� ������� �����������

!$OMP PARALLEL PRIVATE (Int,Rm,RmN,RmO,tux,tux2,k,n,tmp_rm_i,tmp_rm_r)

!$OMP DO SCHEDULE (DYNAMIC, 4)
     do n=istart,iend
       do k=0,grid.num_x-1
         Medium.Ne_old(k,n)=Medium.Ne(k,n)
       end do
     end do
!$OMP END DO  

!     if(MediumFlag == 1) then
!$OMP DO SCHEDULE (DYNAMIC, 4)
       do n=istart,iend
         do k=0,Grid.num_x-1
           Medium.NeN_old(k,n)=Medium.NeN(k,n)
           Medium.NeO_old(k,n)=Medium.NeO(k,n)
         end do
       end do
!$OMP END DO  

!     end if
	 
!$OMP DO SCHEDULE (DYNAMIC, 4)
	 do n=istart,iend
        do k=0,grid.num_x-1
           !if (tn==0) then
             !Int=cdabs(beam.e(k,n,tn))**2
           !else
           Int=cdabs(beam.e(k,n,tn))**2
             !Int=(cdabs(beam.e(k,n,tn-1))**2+cdabs(beam.e(k,n,tn))**2)/2
           !end if

!           if(MediumFlag==2) then
!              absE=cabs(beam.e(k,n,tn))
!           end if

!           Select Case(MediumFlag)
!           Case(1) !===   ������   ================================================

              !����� �������� ��������� �� ���������� �������
              if(IImin>Int) then
                 Medium.NeN(k,n)=Medium.NeN_old(k,n)
                 Medium.NeO(k,n)=Medium.NeO_old(k,n)
				 
                 Medium.Ne(k,n) =Medium.NeO(k,n)+Medium.NeN(k,n)
              else 
!write (*,*) 'Int/IImin', Int/IImin
               if(IImax<Int) then
                 RmN=RN(maxRair)
                 RmO=RO(maxRair)
              else
                 !��������� ������ II �������� ������� � ���������� ����� IImin, ��
                 !������ � ��� ����� ��������� �������. ���� Int-IImin - ���������� ��
                 !������ �������, �� (Int-IImin)/IImin - ���������� ����� �� ����
                 !����������. � ����� ������ ��� ����� �������� � �������.

!                 RmN=RN(nint(Int/IImin-1.)+1) 
!                 RmO=RO(nint(Int/IImin-1.)+1)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 !��� ������� �����������
                 tmp_rm_i=idint(Int/IImin-1.d0)+1
                 !tmp_rm_i=nint(Int/IImin-1.)+1
                 tmp_rm_r=Int/IImin-1.d0*tmp_rm_i
                 RmN=RN(tmp_rm_i)*(1.d0-tmp_rm_r)+RN(tmp_rm_i+1)*tmp_rm_r
                 RmO=RO(tmp_rm_i)*(1.d0-tmp_rm_r)+RO(tmp_rm_i+1)*tmp_rm_r
! ����� ����� �������� �� �������������, �.� ���������������� ������ �� 
! ������� grid.delta_time ����� ������, ���� ������������� �������� ������������� 
! � ������������� RmN,
! �������� �� �����-�����
! � ���-�� ������� y=ln(0.78-N), ������� �������������� N=0.78-exp(y)
! k1 ��� j-1 , �.� I~beam.e(k,n,j-1)
! k2 I ������� j-1 � j, �� y + k1/2
! k3 I ������� j-1 � j, �� y + k2/2
! k4 I ��� j � y + k3

  
 ! write(TMP_File,*) 'nint(Int/IImin-1.)+1', nint(Int/IImin-1.)+1
  
              end if
  !write(*,*) 'RmN*grid.delta_time', RmN, grid.delta_time, RmN*grid.delta_time
  !write(*,*) 'RmO*grid.delta_time', RmO, grid.delta_time, RmO*grid.delta_time

              Medium.NeN(k,n)=0.78-(0.78-Medium.NeN_old(k,n))*exp(-1.*RmN*grid.delta_time(tn)) !0.78 �.�. 78% ���� ������� ������� ��� ����
              Medium.NeO(k,n)=0.21-(0.21-Medium.NeO_old(k,n))*exp(-1.*RmO*grid.delta_time(tn)) !0.21 �.�. 21% ���� ������� ������� ��� ��������
              
              Medium.Ne(k,n) =Medium.NeO(k,n)+Medium.NeN(k,n)!+Medium.NeN_old(k,n)+Medium.NeO_old(k,n))/2
              end if  


              tux=-0.5*Medium.Rpl*Medium.Ne(k,n)

              beam.e(k,n,tn)=beam.e(k,n,tn)*cdexp(cmplx(0.0d0,tux*grid.dz,complex_kind))

!           Case(2) !===   �������   ===============================================
!
!              !1. ������������� ��������� -------------------------------------------
!
!              !����� �������� ��������� �� ���������� �������
!              !(*IonSpeed - ���������� ������ ���������)
!              if(IImin>Int*IonSpeed) then
!                 Rm=0.
!              else if (IImax<Int*IonSpeed) then
!                 Rm=Rmethanol(maxRmethanol)
!              else
!                 !��������� ������ II �������� ������� � ���������� ����� IImin, ��
!                 !������ � ��� ����� ��������� �������. ���� Int-IImin - ���������� ��
!                 !������ �������, �� (Int-IImin)/IImin - ���������� ����� �� ����
!                 !����������. � ����� ������ ��� ����� �������� � �������.
!                 Rm=Rmethanol(nint(Int*IonSpeed/IImin-1.)+1) 
!              end if
!
!              tux=exp(-1.*Rm*grid.delta_time) !��������� ���������� ��������, �.�. ����� ���������� Ne
!              Medium.Ne(k,n)=1.-(1.-Medium.Ne_old(k,n))*tux
!
!              !2. ������� ��������� -------------------------------------------------
!
!              Medium.Ne(k,n)=Medium.Ne(k,n)*exp( Medium.Kbl*absE**3/(1.+Medium.Kpl**2*absE**2)*grid.delta_time )
!
!              !----------------------------------------------------------------------
!
!              tux2=-0.5*Medium.Rpl*Medium.Ne(k,n)*Medium.Kpl*absE/(1.+Medium.Kpl**2*absE**2)
!
!              tux=0.5*Medium.Rpl*Medium.Ne(k,n)/(1.+Medium.Kpl**2*absE**2)
!    
!              beam.e(k,n,j)=beam.e(k,n,j)*exp(tux2*Grid.dz)*cexp(i*tux*Grid.dz)
! 
!           End Select !============================================================

           if (Medium.plasma/=1) then
            if(Medium.Ne(k,n)>=NePlasma) then
             Medium.plasma=1
             NeZXlen_start=length-1.e-6
             NeZXlen=NeZXlen_start
                         
             ! ���������� ����������� ������ ������� �������� ����� ���������� ������������ NePlasma
            end if  
           end if  

           !���� ����� ���������� ���� �� ������� ��� �� z ��������
           !AdaptiveCriterium, �� �� ��������� ���� dz ���������� � 2 ����
           !if(abs(tux*Grid.dz)>=AdaptiveCriterium) DecreseStepPlasma=1 !dzPlasma=abs(tux*grid.dz)

        end do
     end do

!$OMP END DO  
!$OMP END PARALLEL
     

  END IF

 end subroutine Plasma

!**************************************************************************
!                      EnergyAbsorption(istart,iend,j)
!
!                  ���������� ������� ��� ���������
!
! istart - ����� ������� ����� �������
! iend - ������ ������� ����� �������
! j - ����� ���������� ���� ��������
!**************************************************************************
 subroutine EnergyAbsorption(istart,iend,j)

  integer(integer_kind) :: j,istart,iend
  integer(integer_kind),automatic :: k,n
  real(real_kind),automatic :: dum,tux

  IF(EnergyAbsorptionFlag /= 0) THEN

!$OMP PARALLEL PRIVATE (n,k,dum,tux) FIRSTPRIVATE (Isat, j) 
!$OMP DO SCHEDULE (DYNAMIC, 4)
     do n=istart,iend
        do k=0,grid.num_x-1

           Select Case(MediumFlag)
           Case(1) !===   ������   =================================================

              dum=cdabs(beam.e(k,n,j))**2

              !����� �������� ���, ��� ������������� ������ �� ������ ��������� �
              !����. ����� ����� ���� tuxM ��������� � �������������. ����� �����
              !��������, ������ ������� ���������� ������� � ������������� Isat
              if(dum >= Isat) then

                 !���������� ������� ��� ��������� �����
                 tux=-0.5*Medium.RdN2/( dum )*(Medium.NeN(k,n)-Medium.NeN_old(k,n))/grid.delta_time(j)*Grid.dz
                 beam.e(k,n,j)=beam.e(k,n,j)*cdexp(cmplx(tux,0.0d0,complex_kind))

                 !���������� ������� ��� ��������� ���������
                 tux=-0.5*Medium.RdO2/( dum )*(Medium.NeO(k,n)-Medium.NeO_old(k,n))/grid.delta_time(j)*Grid.dz
                 beam.e(k,n,j)=beam.e(k,n,j)*cdexp(cmplx(tux,0.0d0, complex_kind))

              end if

           Case(2) !===   �������   ================================================
   
              dum=cdabs(beam.e(k,n,j))**2

              !����� ��������� ���, ��� ������������� ������ �� ������ ��������� �
              !����. ����� ����� ���� tuxM ��������� � �������������. ����� �����
              !�������� ������ ������� ���������� ������� � ������������� Isat
              if(dum >= Isat) then
                 tux=-0.5*Medium.Rd/( dum )*(Medium.Ne(k,n)-Medium.Ne_old(k,n))/grid.delta_time(j)*Grid.dz
                 beam.e(k,n,j)=beam.e(k,n,j)*cdexp(cmplx(tux,0.0d0, complex_kind))
              end if

           End Select !=============================================================

        end do
     end do
!$OMP END DO  
!$OMP END PARALLEL

  END IF

 end subroutine EnergyAbsorption

end module a_plasma
