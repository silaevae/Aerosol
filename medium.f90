module a_medium
!************************************************************
! ���� �� ������ M_Medium �������� �.
!
! Types:
!
!   tMedium
!--   tPlasmaChannel
!
! Contains:
!
!  subroutine MediumAlloc
!  subroutine MediumDealloc
!  subroutine MediumCalculation
!--  subroutine PlasmaAnalysis
!--  subroutine FindPlasmaChannels
!--  function   CheckNe(k,n,direction)
!--  subroutine MoveToDirection(k,n,direction)
!--  subroutine FindPCBorder(nowk,nown,PlasmaChannel)
!--  subroutine FindPCpoints(PlasmaChannel)
!--  subroutine PlasmaChannelsCalculation
!************************************************************

! use M_Pulse
! use M_Plot
! use msflib
 use a_parameters
 use a_variables
 
 implicit none

 integer(4) :: MediumFlag=1 !���� ������ �����(=1 - ������, =2 - �������)

 !=== ��������� ������� ===
 real(real_kind),&
 parameter :: Air_n   = 1.000273      ,& !���������� ����������� �������
              !Air_n2  = 5.6e-23       ,& !����������� ���������� ������������ ��� ������� �� ����� ����� 800�� [�**2/��]
              Air_n2  = 1.7d-23       ,& !����������� ���������� ������������ ��� ������� �� ����� ����� 800�� [�**2/��]
              Air_k2  = 16.           ,& !�����������, ����������� ��������� ��������� ��������, ��� ������� �� ����� ����� 800�� [��**2/�]
              Air_N0  = 2.68d25       ,& !����� ����������� ������� � ������� [1/�**3]
              Air_WO2 = 19.386499d-19 ,& !��������� ��������� ��������� [��]
              Air_WN2 = 24.994164d-19 ,& !��������� ��������� ����� [��]
              Air_g   = 0.5           ,& !����. ����������� ����������� �������������� � ����������� ���������� ������������ (=0 - ������ ��������������, =1- ������ �����������)
              Air_Pcr = 4.d9 ,& !6.e9          ,& !���������� �������� � ������� [��]
              Air_Ne0 = 0.               !���������� ��������� ���������� � ���������� ��������� ��������� [1/�**3]

 !=== ��������� �������� ===
 real(real_kind),&
 parameter :: Methanol_n   = 1.3286       ,& !���������� ����������� ��������
              Methanol_N0  = 1.5e28       ,& !����� ����������� ������� � �������� [1/�**3]
              Methanol_Wg  = 9.933578e-19 ,& !������ ����������� ���� � �������� [��]
              Methanol_g   = 0.           ,& !����. ����������� ����������� �������������� � ����������� ���������� ������������ (=0 - ������ ��������������, =1- ������ �����������)
              Methanol_Pcr = 8.e6         ,& !���������� �������� � �������� [��]
              Methanol_Sc  = 1.e-19       ,& !������� ������� ������������ [�**2]
              Methanol_Ne0 = 1.e16           !���������� ���������� ���� ������������ � ���������� ��������� ��������� [1/�**3]

!**************************************************************************
!                               MEDIUM
!**************************************************************************
 type tMedium

  real(real_kind) :: n   ,& !���������� �����������
                     n2  ,& !����������� ���������� ������������
                     N0  ,& !����� ����������� ������ (�������) �����
                     Wp0 ,& !������������� ���������� �������
                     R   ,& !�������� ���������� ������������
                     Rpl ,& !�������� ���������� ������������
                     g   ,& !����. ����������� ����������� �������������� � ����������� ���������� ������������ (=0 - ������ ��������������, =1- ������ �����������)
                     Pcr ,&
                     Ne0    !���������� ��������� ���������� ��� ���������� ���� ������������ � ���������� ��������� ���������

  real(real_kind) :: Wc0,& !������������� ������� ������������
                     Kpl,& !Medium.Wc0/beam.abs_w0
                     Kbl   !����. ������� ���������

  real(real_kind),pointer :: Ne(:,:)     ,& !������������ ��������� ���������� ��� ���������� ���� ������������
                             Ne_old(:,:)    !������������ ��������� ���������� ��� ���������� ���� ������������ � ����������� ���� �� �������

  !--- ���������� ��� ������� ---------------------------------------------
  real(real_kind) :: lO2  ,& !������� �������������� �������� ��� ��������� ���������
                     lN2  ,& !������� �������������� �������� ��� ��������� �����
                     RdO2 ,& !�������� ������������ ���������� ���������
                     RdN2    !�������� ������������ ���������� �����

  real(real_kind),pointer :: NeN(:,:)    ,& !������������ ��������� ���������� �����
                             NeN_old(:,:),& !������������ ��������� ���������� ����� � ����������� ���� �� �������
                             NeO(:,:)    ,& !������������ ��������� ���������� ���������
                             NeO_old(:,:)   !������������ ��������� ���������� ��������� � ����������� ���� �� �������

  !--- ���������� ��� �������� --------------------------------------------
  real(real_kind) :: Sc ,& !������� ������� ������������
                     Wg ,& !������ ����������� ����
                     l  ,& !������� �������������� �������� ��� ���������
                     Rd    !�������� ������������ ����������
  !------------------------------------------------------------------------

  integer(4) :: plasma=0 !���� ����������� ������ (=1 - ������ ������������)

 end type tMedium

 type(tMedium) :: Medium

!**************************************************************************
!                              PlasmaChannel
!**************************************************************************
! type tPlasmaChannel
!
!  real(real_kind) :: length,& !����� ����������� ������
!                     RhoNe    !��������� ������ � ���������� ������
!
!  real(real_kind) :: x,y !���������� ������ ������� ����������� ������
!
!  integer(4) :: NumBor !����������� Border
!  type(tGridPoint),pointer :: Border(:) !������� ����������� ������
!
!  integer(4) :: NumPoints !����������� PCpoints
!  type(tGridPoint),pointer :: PCpoints(:) !����� ����������� ������
!
! end type tPlasmaChannel
!
! integer(4) :: NumPlasmaChannel=0 !���������� ���������� �������
!
! type(tPlasmaChannel),pointer :: PlasmaChannel(:)
!
! integer(4),pointer :: PCMap(:,:) !����� ���������� �������
!
!**************************************************************************
                                 CONTAINS
!**************************************************************************
!**************************************************************************
!                               MediumAlloc
!
!                    �������� ������ ��� ������� Medium
!**************************************************************************
 subroutine MediumAlloc
  integer(4),automatic :: k,n
    
  if(TimeFlag/=0 .and. PlasmaFlag/=0) then

    allocate(Medium.Ne(0:grid.num_x,0:grid.num_x))
    allocate(Medium.Ne_old(0:grid.num_x,0:grid.num_x))

    if(MediumFlag==1) then !������
      allocate(Medium.NeN(0:grid.num_x,0:grid.num_x))
      allocate(Medium.NeO(0:grid.num_x,0:grid.num_x))
      allocate(Medium.NeN_old(0:grid.num_x,0:grid.num_x))
      allocate(Medium.NeO_old(0:grid.num_x,0:grid.num_x))
    end if

    do n=0,grid.num_x-1
      do k=0,grid.num_x-1
         Medium.Ne(k,n)=Medium.Ne0
      end do
    end do
    if(MediumFlag == 1) then
      do n=0,grid.num_x-1
        do k=0,grid.num_x-1
          Medium.NeN(k,n)=Medium.Ne0
          Medium.NeO(k,n)=Medium.Ne0
        end do
      end do
    end if

!���� ������� ������
!    if(PAflag /= 0) then
!      allocate(PCMap(0:grid.num_x,0:grid.num_x))
!    end if

  end if

 end subroutine MediumAlloc

!**************************************************************************
!                              MediumDealloc
!
!                ����������� ������ ������� ��������� Medium
!**************************************************************************
 subroutine MediumDealloc

  if(TimeFlag/=0 .and. PlasmaFlag/=0) then

    deallocate(Medium.Ne)
    deallocate(Medium.Ne_old)

    if(MediumFlag==1) then !������
      deallocate(Medium.NeN,Medium.NeO)
      deallocate(Medium.NeN_old,Medium.NeO_old)
    end if


!���� ������� ������
!    if(PAflag /= 0) then
!      deallocate(PCMap)
!    end if

  end if

 end subroutine MediumDealloc

!**************************************************************************
!                            MediumCalculation
!
!                       ��������� ���������� Medium
!**************************************************************************
 subroutine MediumCalculation

  !������� ���������� �����
  
  select case(MediumFlag)
  case(1) !=== ������ =====================================================

    Medium.Pcr = Air_Pcr
    Medium.n   = Air_n
    Medium.n2  = beam.RcrAB*beam.abs_lambda0**2/(8.*pi*Medium.n*Medium.Pcr) !Air_n2
    Medium.N0  = Air_N0
    Medium.g   = Air_g
    Medium.Ne0 = Air_Ne0
    Medium.lO2 = Air_WO2 / beam.Wph !������� �������������� �������� ��� ��������� ���������
    Medium.lN2 = Air_WN2 / beam.Wph !������� �������������� �������� ��� ��������� �����
    Medium.RdO2= 1./beam.abs_i0*beam.abs_ld*Medium.N0/beam.pulse_duration*Medium.lO2*beam.Wph !�������� ������������ ���������� ���������
    Medium.RdN2= 1./beam.abs_i0*beam.abs_ld*Medium.N0/beam.pulse_duration*Medium.lN2*beam.Wph !�������� ������������ ���������� �����

    Medium.Wp0 = sqrt(el**2*Medium.N0/eps0/me) !������������� ���������� �������
    Medium.R   = 2.*beam.abs_k0**2/Medium.n*Medium.n2*beam.abs_a0*beam.abs_a0*beam.abs_i0 !�������� ���������� ������������
    Medium.Rpl = beam.abs_k0*beam.abs_ld/Medium.n**2*Medium.Wp0**2/beam.abs_w0**2 !�������� ���������� ������������
    

  case(2) !=== ������� ====================================================

    Medium.Pcr= Methanol_Pcr
    Medium.n  = Methanol_n
    Medium.n2 = beam.RcrAB*beam.abs_lambda0**2/(8.*pi*Medium.n*Medium.Pcr)


!   marb=0.367/sqrt((sqrt(Medium.R/beam.RcrAB)-0.852)**2-0.0219)
!write(unit,'(A,f18.7)') ' ���������� ����������� ��������� � �������� P/Pcr=',Medium.R/beam.RcrAB
   
    Medium.N0 = Methanol_N0
    Medium.g  = Methanol_g
    Medium.Ne0= Methanol_Ne0/Medium.N0
    Medium.Sc = Methanol_Sc
    Medium.Wg = Methanol_Wg
    Medium.l  = Medium.Wg / beam.Wph !������� �������������� �������� ��� ���������
    Medium.Rd = 1./beam.abs_i0*beam.abs_ld*Medium.N0/beam.pulse_duration*Medium.l*beam.Wph !�������� ������������ ����������

    Medium.Wp0 = sqrt(el**2*Medium.N0/eps0/me) !������������� ���������� �������
    Medium.R   = 2.*beam.abs_k0**2/Medium.n*Medium.n2*beam.abs_a0**2*beam.abs_i0 !�������� ���������� ������������
    Medium.Rpl = beam.abs_k0**2*beam.abs_a0**2/Medium.n**2*Medium.Wp0**2/beam.abs_w0**2 !�������� ���������� ������������

  end select
  !========================================================================

 end subroutine MediumCalculation

!!**************************************************************************
!!                           PlasmaAnalysis
!!
!!                         ����������� ������
!!**************************************************************************
! subroutine PlasmaAnalysis
!
!  NumPlasmaChannel=0
!
!  call FindPlasmaChannels !���� ���������� ������
!
!  call PlasmaChannelsCalculation !��������� ��������� ���������� �������
!  
! end subroutine PlasmaAnalysis
!
!
!!**************************************************************************
!!                            FindPlasmaChannels
!!
!!            ���� ���� ����� ������������� ���������� �������
!!**************************************************************************
! subroutine FindPlasmaChannels
!
!  integer(4) :: i,k,n,finish,posn,posl
!
!  real(real_kind) :: maximum
!
!  type(tPlasmaChannel),allocatable :: PCold(:)
!
!  type(tGridPoint) :: MaxPoint
!
!  !������������� ����� ���������� �������
!  do n=0,grid.num_x
!  do k=0,grid.num_x
!    if(Medium.Ne(k,n) >= NePlasma) then
!      PCMap(k,n)=-1 !�������� ����� ��� ������ ����
!    else
!      PCMap(k,n)=0  !�������� ����� ��� ��� ������
!    end  if
!  end do
!  end do
!
!  finish=0
!
!  DO WHILE(finish==0)
!
!    !���� �������� �� ������������� ������ �������
!    maximum=0.
!    do n=0,grid.num_x
!    do k=0,grid.num_x
!      if(PCMap(k,n)==-1 .and. Medium.Ne(k,n)>maximum) then
!        maximum=Medium.Ne(k,n)
!        MaxPoint.k=k
!        Maxpoint.n=n
!      end if
!    end do
!    end do
!
!    if(maximum /= 0.) then !������ ����� ���������� �����
!
!       if(NumPlasmaChannel/=0) then !��������� ������ ���������� ������
!          allocate(PCold(NumPlasmaChannel))
!          PCold=PlasmaChannel
!          deallocate(PlasmaChannel)
!       end if
!
!       NumPlasmaChannel=NumPlasmaChannel+1 !����������� ���������� ���������� ������� �� �������
!       allocate(PlasmaChannel(NumPlasmaChannel))
!
!       if(allocated(PCold)) then
!          PlasmaChannel(1:NumPlasmaChannel-1)=PCold
!          deallocate(PCold)
!       end if
!
!       !���� ������� ������ ����������� ������
!       call FindPCBorder(MaxPoint.k,MaxPoint.n,PlasmaChannel(NumPlasmaChannel))
!
!       !���� ���������� ����� ������ ����������� ������
!       call FindPCpoints(PlasmaChannel(NumPlasmaChannel))
!
!       !�������� ����� ������������� ������ ����������� ������
!       do i=1,PlasmaChannel(NumPlasmaChannel).NumPoints
!          PCMap(PlasmaChannel(NumPlasmaChannel).PCpoints(i).k,PlasmaChannel(NumPlasmaChannel).PCpoints(i).n)=1
!       end do
! 
!    else !������ ��� ���������� �������
!
!       finish=1
!
!    end if
!
!  END DO
!
!  !��������� ������������ ���������� ���������� �������
!  i=0
!  do n=0,grid.num_x
!  do k=0,grid.num_x
!     if(PCMap(k,n)==1) i=i+1
!  end do
!  end do
!
!  if(i/=sum(PlasmaChannel.NumPoints)) then !������� ��������� � �������������� ����
!     call WriteHistory('������ � ��������� FindPlasmaChannel',&
!                       '���������� ����� ����� ��� ������������'// &
!                       ' ������ �� ��������� � ����������� ����� ������������� ���������� �������')
!  end if
!
! end subroutine FindPlasmaChannels
!
!!**************************************************************************
!!                           CheckNe(k,n,direction)
!!
!! ��������� ������� ������ � ���� ����� �������� � (k,n) � �����������
!! ���������� ���������� direction. CheckNe=0 ���� � ��� ��� ������ � =1
!! ���� ������ ��� ����.
!!
!! k,n - ���������� ���� ����� ����� � ������� ����������� ������� ������
!! direction - ����������� ������ (1-����, 2-�����, 3-���, 4-����)
!!**************************************************************************
! function CheckNe(k,n,direction)
!
!  integer(4) :: k,n,direction,CheckNe
!
!  CheckNe=0
!
!  select case(direction)
!  case(1) !�����
!    if(Medium.Ne(k,n+1)>=NePlasma) CheckNe=1
!  case(2) !������
!    if(Medium.Ne(k+1,n)>=NePlasma) CheckNe=1
!  case(3) !����
!    if(Medium.Ne(k,n-1)>=NePlasma) CheckNe=1
!  case(4) !�����
!    if(Medium.Ne(k-1,n)>=NePlasma) CheckNe=1 
!  end select
!
! end function CheckNe
!
!!**************************************************************************
!!                  MoveToDirection(k,n,direction)
!!
!!       ������ ������� ���� ����� � ����������� direction
!!
!! k,n - ���������� ���� ����� �� �������� ������������
!! direction - ����������� �������� (1-����, 2-�����, 3-���, 4-����)
!!**************************************************************************
! subroutine MoveToDirection(k,n,direction)
!
!  integer(4) :: k,n,direction
!
!  select case(direction) !��������� � �������� �����������
!  case(1)
!    n=n+1
!  case(2)
!    k=k+1
!  case(3)
!    n=n-1
!  case(4)
!    k=k-1
!  end select
!
! end subroutine MoveToDirection
!
!!**************************************************************************
!!               FindPCBorder(nowk,nown,PlasmaChannel)
!!
!!                ������� ������� ����������� ������
!!
!! nowk,nown - ���������� ���� ����� �� �������� ���������� �����
!! PlasmaChannel - ���������� ����� ��� ������� ������
!!**************************************************************************
! subroutine FindPCBorder(nowk,nown,PlasmaChannel)
!
!  integer(4) :: nowk,nown,bor,direction,nstop,checkresult,&
!                           left
!  type(tGridPoint) :: Border((grid.num_x+1)*(grid.num_x+1)),StartPoint
!  type(tPlasmaChannel) :: PlasmaChannel
!
!  Border.k=-1
!  Border.n=-1
!  bor=0
!
!  StartPoint.k=-1
!  StartPoint.n=-1
!
!  !=== ������� �� ������� =================================================
!  nstop=0
!  direction=1 !�������� � �������� ���������� ����������� �������� �������� �����
!  do while(nstop==0)
!
!    checkresult=CheckNe(nowk,nown,direction) !��������� �������� (nowk,nown) ����� � ����������� direction
!
!    if(checkresult==1) then
!       call MoveToDirection(nowk,nown,direction) !��������� � ��������� �����������
!    else
!       StartPoint.k=nowk !���������� ����� � ������� ������ ����� �������
!       StartPoint.n=nown
!       bor=bor+1
!       Border(bor)=StartPoint
!       direction=2 !��� ���������� ������� ������������ �������, ����� ������� ������ ���������� �����
!       nstop=1
!    end if
!
!  end do
!
!  !=== ������� ������� ====================================================
!  nstop=0
!  DO WHILE(nstop==0)
!
!    checkresult=CheckNe(nowk,nown,direction) !��������� �������� (nowk,nown) ����� � ����������� direction
!
!    if(checkresult==1) then !��������� ����� ����������� ������
!
!      call MoveToDirection(nowk,nown,direction) !��������� � ��������� �����������
!
!      bor=bor+1
!      Border(bor).k=nowk
!      Border(bor).n=nown
!
!      !������� ��� �����
!      left=direction-1
!      if(left==0) left=4
!
!      checkresult=CheckNe(nowk,nown,left)
!
!      if(checkresult==1) then
!        call MoveToDirection(nowk,nown,left)
!        bor=bor+1
!        Border(bor).k=nowk
!        Border(bor).n=nown
!      end if
!
!    else !�������� ����������� �������� �� ������� �������
!
!      direction=direction+1 !�������� ����������� ��������
!      if(direction==5) direction=1
!
!      checkresult=CheckNe(nowk,nown,direction) !��������� �������� ����� � ����� �����������
!    
!      if(checkresult==1) then !�������� ����� ����������� ����������� ������
!
!         call MoveToDirection(nowk,nown,direction) !��������� � ��������� �����������
!
!         bor=bor+1
!         Border(bor).k=nowk
!         Border(bor).n=nown
!
!         direction=direction-1 !������������ � ��������� �����������
!         if(direction==0) direction=4
!
!      end if
!
!    end if
!
!    !���� ������ � ����� � ������� �������� ����� �������, �� ������� �������
!    if(nowk==StartPoint.k .and. nown==StartPoint.n) nstop=1
!
!  END DO
!
!  PlasmaChannel.NumBor=FindInArray(-1,Border(:).k,grid.num_x*grid.num_x)-1
!  allocate(PlasmaChannel.Border(PlasmaChannel.NumBor))
!  PlasmaChannel.Border=Border(1:PlasmaChannel.NumBor)
!
! end subroutine FindPCBorder
!
!!**************************************************************************
!!                    FindPCpoints(PlasmaChannel)
!!
!!    ������� ���� ����� ������������� ������� ����������� ������
!!
!! PlasmaChannel - ���������� ����� ��� ���� ������
!!**************************************************************************
! subroutine FindPCpoints(PlasmaChannel)
!
!  type(tPlasmaChannel) :: PlasmaChannel
!  integer(4) :: k,n,posy,pos,i,CrossingPoints(PlasmaChannel.NumBor),&
!                           numcross,pc,count,pcX,pcY
!  type(tGridPoint) :: PCpoints(grid.num_x*grid.num_x),PCpointsY(grid.num_x*grid.num_x),&
!                      PCpointsX(grid.num_x*grid.num_x)
!
!  pc=0
!  PCpoints.k=-1
!  PCpoints.n=-1
!
!  !��� ������� � ���� ����� ������� ����� ��������� ������
!  do n=0,grid.num_x
!
!    posy=FindInArray(n,PlasmaChannel.Border(:).n,PlasmaChannel.NumBor)
!
!    if(posy/=-1) then !���������� ����� � ����� n ������������� �������
!
!       !������� ����� ����������� � �������� ��� ������� n
!       CrossingPoints=grid.num_x*grid.num_x
!       numcross=0
!       do i=1,PlasmaChannel.NumBor
!          if(PlasmaChannel.Border(i).n==n) then
!             numcross=numcross+1
!             CrossingPoints(numcross)=PlasmaChannel.Border(i).k
!          end if
!       end do
!
!       !��������� ������ ����� �����������
!       count=numcross
!       call sortqq(loc(CrossingPoints(1:numcross)),count,SRT$INTEGER4)
!       if(count /= numcross) stop 'Sorting Error in procedure FindPCpoints'
!
!       do k=0,grid.num_x
!       
!          if(PCMap(k,n)/=-1) cycle !�� ������������� ����� ������������� ������ ������� � ����� ��� ��� ������
!
!          do i=1,numcross
!
!             if(k==CrossingPoints(i)) then !����� ����������� �������
!                pc=pc+1
!                PCpoints(pc).k=k
!                PCpoints(pc).n=n
!                exit
!             end if
!
!             if(i==numcross) exit !������ ��� ����������� ����� ������� �����������
!
!             if(k>CrossingPoints(i) .and. k<CrossingPoints(i+1)) then !����� ����� ����� ������� ����������� � ��������
!                pc=pc+1
!                PCpoints(pc).k=k
!                PCpoints(pc).n=n
!             end if
!
!          end do
!
!       end do !k=0,grid.num_x
!
!    end if
!
!  end do !n=0,grid.num_x
!
!  PlasmaChannel.NumPoints=FindInArray(-1,PCpoints.k,grid.num_x*grid.num_x)-1
!  allocate(PlasmaChannel.PCpoints(PlasmaChannel.NumPoints))
!  PlasmaChannel.PCpoints=PCpoints(1:PlasmaChannel.NumPoints)
!
! end subroutine FindPCpoints
!
!!**************************************************************************
!!                     PlasmaChannelsCalculation
!!
!! ������������ �������������� ���������� �������
!!**************************************************************************
! subroutine PlasmaChannelsCalculation
!
!  integer(4) :: i,j
!  integer(4),pointer :: k(:),n(:)
!  real(real_kind) :: sum,sumx,sumy
!
!  DO i=1,NumPlasmaChannel
!
!    k=>PlasmaChannel(i).PCpoints.k
!    n=>PlasmaChannel(i).PCpoints.n
!   
!    !���� ���������� ������ ������� ����������� ������
!    sum=0.
!    do j=1,PlasmaChannel(i).NumPoints
!      sum=sum+Medium.Ne(k(j),n(j))
!    end do
!
!    sumx=0.
!    sumy=0.
!    do j=1,PlasmaChannel(i).NumPoints
!      sumx=sumx+Grid.x(k(j))*Medium.Ne(k(j),n(j))
!      sumy=sumy+Grid.y(n(j))*Medium.Ne(k(j),n(j))
!    end do
!
!    PlasmaChannel(i).x=sumx/sum
!    PlasmaChannel(i).y=sumy/sum
!
!    !���� ��������� ������ � ���������� ������
!    sum=0.
!    do j=1,PlasmaChannel(i).NumPoints
!      sum=sum + (Grid.x(k(j)+1)-Grid.x(k(j)-1))/2.* &
!                (Grid.y(n(j)+1)-Grid.y(n(j)-1))/2.* Medium.Ne(k(j),n(j))
!    end do
!
!    PlasmaChannel(i).RhoNe=sum*beam.abs_a0*beam.abs_a0*Medium.N0 !��������� � [1/�]
!
!  END DO
!
! end subroutine PlasmaChannelsCalculation
!
end module a_medium
