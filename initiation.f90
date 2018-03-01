!////////////////////////////////////////////////////////////////////////// 
!				 ������� ��������� ����������		 
!/////////////////////////////////////////////////////////////////////////
subroutine AerosolInitDialog(res)
  integer::res
  character :: str
  
    write(*,'(a)') 'Project AEROSOL is running'
    write(*,'(a)') '<><><><><><><><><><><><><>'
    write(*,*)
    write(*,'(a)') 'Select one of next options:'
    write(*,'(a)') '(N): Start calculation'
    write(*,'(a)') '(S): Start calculation with intensity Spectrum output'
    write(*,*)
    read(*,*) str
    write(*,*)

    select case(str)
    case('n')    
	   res=0
	case('s') 
       res=1
    case default
       res=-1
	   write(*,'(a)') 'Program aborted'	   
	   call exit
	end select

end subroutine AerosolInitDialog


subroutine initiation

 use a_parameters
 use a_types
 use a_variables
 use a_medium
 
 implicit none
 real(real_kind) :: time_step_param, time_step_param_1, t_s_1, t_s_2, t_s_3, t_s_4, dt_s
 integer(integer_kind) :: t_index, tmp_t 
 

 beam.pulse_duration = 140.d-15 !140.e-15 !�������� ������������ �������� e-1 �� ������������� [s]

 
 grid.time_number = 500 !����� ��c������������ ��������� �����
 grid.time_start = -2.5 !-2.0 !������ ���������������� ��������� ������� [������������ ��������] 
 grid.time_end = 2.0 !����� ���������������� ��������� ������� [������������ ��������]
 
 allocate(grid.current_time(0:grid.time_number)) !��� ������������ ����� �� �������
 allocate(grid.delta_time(0:grid.time_number)) 
 if (grid.time_number>1) then 
   TimeFlag=1 !=0 - ������������ ������; =1 - �������.
   !grid.delta_time = (grid.time_end-grid.time_start)/(grid.time_number-1) !�������� �������
   grid.dt = (grid.time_end-grid.time_start)/(grid.time_number-1) !�������� ������� ��� ����������� ����
   time_step_param=4. !��������� ���������� � ������������ ����
   time_step_param_1=sqrt(time_step_param) !��������� �������� �����
   t_s_1=-1.
   t_s_2=-0.7
   t_s_3=0.45
   t_s_4=0.75
   
   dt_s=(grid.time_end-grid.time_start+(time_step_param_1-1)*(t_s_4-t_s_1)+(time_step_param-time_step_param_1)*(t_s_3-t_s_2))/((grid.time_number-1)*time_step_param)
 
   grid.time_0=int((t_s_1-grid.time_start)/dt_s/time_step_param+(t_s_2-t_s_1)/dt_s/time_step_param_1+(0-t_s_2)/dt_s)
   
   grid.current_time(grid.time_0)=0
   
   t_index=grid.time_0
   tmp_t=t_index
   do while ((grid.current_time(t_index)>t_s_2).and.(t_index>0))
     t_index=t_index-1
     grid.current_time(t_index)=grid.current_time(tmp_t)+dt_s*(t_index-tmp_t)
   end do
   tmp_t=t_index
   do while ((grid.current_time(t_index)>t_s_1).and.(t_index>0))
     t_index=t_index-1
     grid.current_time(t_index)=grid.current_time(tmp_t)+dt_s*time_step_param_1*(t_index-tmp_t)
   end do
   tmp_t=t_index
   do while (t_index>0)
     t_index=t_index-1
     grid.current_time(t_index)=grid.current_time(tmp_t)+dt_s*time_step_param*(t_index-tmp_t)
   end do
   
   t_index=grid.time_0
   tmp_t=t_index
   do while ((grid.current_time(t_index)<t_s_3).and.(t_index<grid.time_number))
     t_index=t_index+1
     grid.current_time(t_index)=grid.current_time(tmp_t)+dt_s*(t_index-tmp_t)
   end do
   tmp_t=t_index
   do while ((grid.current_time(t_index)<t_s_4).and.(t_index<grid.time_number))
     t_index=t_index+1
     grid.current_time(t_index)=grid.current_time(tmp_t)+dt_s*time_step_param_1*(t_index-tmp_t)
   end do
   tmp_t=t_index
   do while (t_index<grid.time_number)
     t_index=t_index+1
     grid.current_time(t_index)=grid.current_time(tmp_t)+dt_s*time_step_param*(t_index-tmp_t)
   end do

   do t_index=1,grid.time_number-1	
     grid.delta_time(t_index)=grid.current_time(t_index)-grid.current_time(t_index-1)
   end do
   grid.delta_time(0)=grid.delta_time(1)


   open(NIntFNum1,FILE=ResultDirName//'time_scale.dat',access='sequential',status='unknown') 
   write(NIntFNum1,*) 'index    time    time_step'
   do t_index=0,grid.time_number-1	
     write(NIntFNum1,*) t_index, grid.current_time(t_index), grid.delta_time(t_index) 
   end do
   close(NIntFNum1)

   
 else 
   TimeFlag=0 !=0 - ������������ ������; =1 - �������.
   grid.delta_time(0)=0
   grid.dt=0 
 end if 
 
 beam.abs_a0 = 2.5d-3 !1.0e-2 !��������� ������ ����� �� ��� X [�]
 beam.abs_r0 = 0 !-0.48     !�������� ��������� ������ [�]
 beam.abs_lambda0 = 0.8d-6 !����� ����� ���������   [�]
 beam.abs_e0 = 1.d0     !��������� ���������     [�/�]
 beam.abs_i0 = 0.5*eps0*sqrt(n0)*(beam.abs_e0**2)  !��������� ������������� [��/�**2]
 
 beam.r = 50 !1.0e-1 !50 !30 ! ��������� �������� ����� � ����������� (�������� ���������������)

 beam.abs_i0 = beam.r*Air_Pcr/(pi*(beam.abs_a0**2))    !��������� ������������� [��/�**2]
 !write(*,*) 'beam.abs_i0',beam.abs_i0

 beam.TheorEnergy = pi*sqrt(pi)*(beam.abs_a0**2)*beam.pulse_duration*beam.abs_i0
! Pulse.AbsEnergy   = 0.
 beam.TheorPower  = pi*(beam.abs_a0**2)*beam.abs_i0

! beam.I0 = 
 beam.RcrAB =3.77 !��� ��������� �������� ���������������� �����
  
 AdaptiveCriterium= pi/15.d0 !���� ���������� ����� ���� �� ��� �������� ���� ��������, �� dz ����������
 Isat=1.d14/beam.abs_i0  !����� �������������, ������� � �������� ��������� ���������� ��� ��������� [��/�**2]
 Max_Int_Const=1000 !����� ���������� �������������, ����� ��������������� ����
 !*********** ��������
 !TimeFlag=0 


 PlasmaFlag=1           !���� ���������� ������������
 EnergyAbsorptionFlag=1  !���� ���������� ��� ���������
 NePlasma=1.d-7          !������������ ����������, ��� ���������� ������� ���������, ��� ������������ ������ [Medium.N0]

 beam.test=1		!=0, ���� ������ ���� ������� ����� � 1-� ������ � ������
                    !=-1 ������ ������ ���
 beam.linear_abs=0  !=1, ����� ��������� �������� ����������                     

 media.droptest=-1  !=0 ����� ��������� ��������
                    !=-1 ����� ��������� ��������, ������ ������������ � ���� (����� ������, ������, x, y)
                    !=1 ��������� � ������ ������ �������� �� �����
                    
 media.turbulence=0 !=1 �������������� ����, �������� ������� ������
                    !=0 �������������� ���  
 SpFilterFlag=0    !���� ��������� ������������� ����
 FilterWidth=0.5   !������ ����� ����, �� ������� ������������ ������
 
 media.droptest_filename='ScatFunc\droplet_position.txt'   !'_Graph\droplet_position.txt' !��� ����� ���� �������� ������ � �������� ������                    
 
 media.abs_drop_concentration = 10.d6 !������������ ������ [�-3] 
 media.abs_drop_size = 15.d-6 ! ������ ������ [�]	� ������ ��������������� ��������
 media.area_x = 3  !������ ������� ������������ ������ [������ �����]
 counter_loc_drop = 250 !����� ����, � ������ �������� ��������� �����, ���� beam.test=0
 
!****************** �������������� �������� **************************************
 media.poly=0			!1-������������ �������������� ��������, 0-��������������
 media.droplet_type=nint(media.abs_drop_size*1000000) !�� ��������� ��� ��������������� ��������
!���������� ������, ������������ ���������, ��������� - 9 ����� ������ (2mkm - 10mkm)
!*********************************************************************************

!*************** ���������********************************************************
 media.dispersion=0	!1-��������� ����, 0-���
!*********************************************************************************

!*************** ��������������� *************************************************
 media.self_focusing=1	!1-��������������� ����, 0-���
!*********************************************************************************

 grid.num_x = 1024 !����� ����� �� Ox (= Oy)
 grid.size_x = 8 !������ ���������� ������ [������ �����]
 grid.abs_size_z = 7 !������ ���������� ������ [�]
 grid.abs_dz_diffraction_step= 0.015 !���������� ��� ������� ��������� [�] 
 grid.num_z = nint(grid.abs_size_z/grid.abs_dz_diffraction_step) !����� ����� �� Oz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !grid.abs_dz_step = 5.e-2 !2.e-3 !���������� ��� ������� ��������� / ������ ������������ ������[�] 

 beam.abs_k0 = 2*pi/beam.abs_lambda0               !�������� �����          [�-1]
 beam.abs_w0 = c0*beam.abs_k0                      !������� ���������       [��]
 beam.Wph = hpl*beam.abs_w0               !������� ������ ������ ��������� [��]
  
 beam.abs_ld = beam.abs_k0*(beam.abs_a0**2)    !������������� ����� ����� [�]
 beam.abs_lds = (beam.pulse_duration)**2*(1.e15)**2/Air_k2 !������������� ����� [�] 
 
 !marb=0.367/sqrt((sqrt(beam.r)-0.852)**2-0.0219)*beam.abs_ld !���������� ��������������� [�]
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!********************** ������� ����� ������� ������� ***************************
! �� ����� ����� ��� ����  
! grid.abs_limit_z=49.06 !49.06			    !������� ������� ������� ���� � ����� [�]
! grid.limit_z=grid.abs_limit_z/beam.abs_ld	!������� ������� ������� ���� � ����� [������������� ����� �����]
! grid.limit_z=grid.abs_size_z/beam.abs_ld	!������� ������� ������� ���� � ����� [������������� ����� �����]
!********************************************************************************

 grid.size_z = grid.abs_size_z/beam.abs_ld !������ ���������� ������ [������������� ����� �����]

 grid.h = grid.size_x/(grid.num_x-1)     !���������� ��� ����� [������ �����]
 grid.dz = grid.abs_dz_diffraction_step/beam.abs_ld    !���������� ��� ����� (�������������) [������������� ����� �����]
 grid.abs_size_x = grid.size_x*beam.abs_a0 !������ ���������� ������ [�]
! grid.abs_size_z = grid.size_z*beam.abs_ld !������ ���������� ������ [�]

!������� ����� ������ � ��������� ����
! media.avr_skin_drop_number = media.abs_drop_concentration*4*((media.area_x*beam.abs_a0)**2)*grid.abs_size_z/grid.num_z 
! ������� ����� ������ �� ������������� �����
! media.avr_skin_drop_number = media.abs_drop_concentration*4*((media.area_x*beam.abs_a0)**2)*grid.abs_size_z 

! ����������, ����� ������� ���������� ��������� ������������ ������ delta_z_aerosol
 media.delta_z_aerosol= grid.dz !0.009d0/beam.abs_ld ! ����� ���������� dz
 media.start_z_aerosol= -1.e-5 !0.819d0/beam.abs_ld  !��������� ������� ������������ ������
 media.length_aerosol = 0.05d0/beam.abs_ld !������ ������������ ���� [������������� �����]

 media.end_z_aerosol= 200./beam.abs_ld !media.start_z_aerosol+media.length_aerosol !����� ������������ ���� 
 !������� ����� ������ � ��������� ����
 media.avr_skin_drop_number = media.abs_drop_concentration*4*((media.area_x*beam.abs_a0)**2)*beam.abs_ld*media.delta_z_aerosol

!��������� ������� ������� ��������������
media.delta_z_turb= 0.5d0/beam.abs_ld !���������� ����� ��������
media.start_z_turb= -1.e-5  !��������� ������� ������������� ������

! �������������� �����������
 media.z_focus_lens= 200/beam.abs_ld          ! ������������ ������������ �����  
 beam.abs_r1 = 0 !-0.6     !�������������� ����������� [�]

!**************** ����� � ���� ***************************************************
 write_spectrum_flag=0 ! ���� ������ ������� �������������
 
 media.build_field=1 !=1, ���� ����� �� ������ i-� ���� ������� ���� ��� Surfer'�; =0, ���� ���� ������� �� ����
 media.build_field_dz=grid.dz !5*grid.dz ! ����������� ����� ����� i 
 media.build_field_start_z= grid.dz-1.e-5 ! ���������� ���������� ������� �������


 media.build_field_dz=(1.-1.e-5)*media.build_field_dz !����� �������������� ������ ����������

 
 NeZXlen_start= 100 !2.67               !���������� �� ������� �������� ������������ ���������� � ���������� �������[Ld]
 !FlZXYlen_start= 50*grid.dz            !����������, �� ������� �������� ������������ ���������� � ������� [Ld] 
 PdzAfterPlasma= 0.015d0/beam.abs_ld !media.delta_z_aerosol ! 0.00025 !0.0025   !��� ������ �������� ����� ����������� ������ [Ld]

!***************************************************

 call init_gauss(grid.num_x, grid.time_number)

! lenslet �� �����  
! call phase_modulation(grid.size_x,grid.num_x-1,24,0.3,beam.e,2.1e-2,-1)

 
! subroutine phase_modulation(a_0,md1,N_d,k,W,R_f,sign)
!real a_0 - ������ ������� [������ �����] = grid.size_x
!integer md1 - ����� ���������� ����� - 1; grid.num_x-1 
!integer N_d - ���������� (��������) ����� ���� (������) = 32
!real k - ������������ ����������� ����������� (h=k*d) ��������� ������� ����������� � ������� ����� 
!complex W(0:md1,0:md1) - ����
!real R_f - �������� ��������� ���� 2.1e-2
!integer sign - ���� � ���������� : -1 ��� +1
 
 contains

!�������� ������ ��� beam.e(0:num_,0:num_) � �������� � ��� �����
 subroutine init_gauss(num_,tnum_)
   integer(integer_kind) num_, tnum_, n,k,tn
   real(real_kind) :: num_2, tnum_2, tix, koef,tkoef, koef_r, time_k
   allocate(beam.e(0:num_-1,0:num_-1,0:tnum_-1))

!   allocate(beam.etmp(0:num_-1,0:num_-1,0:129))
!
!  do tn=0,129
!   do n=0,num_-1
!   do k=0,num_-1
!	   beam.etmp(k,n,tn)=k+n+tn
! end do    
! end do   
! end do   


! ��������� �����
!   if (media.self_focusing==1) then
!      allocate(beam.e2(0:num_-1,0:num_-1))
!   end if

! num_2=num_/2		!�������� ����c� ���������� �� ���� num_2
  num_2=num_/2-0.5 !�������� ����c� ���������� ���������� ����� ������ num_2 � num_2-1
  tnum_2=tnum_/2-0.5
   !koef=(grid.h**2)/2. 
   tkoef=(grid.dt**2)/2.
  if (beam.abs_r0/=0) then 
    koef_r=beam.abs_k0*(beam.abs_a0**2)/beam.abs_r0
  else
    koef_r=0
  end if  
 
  do tn=0,tnum_-1
   time_k = exp(-((grid.current_time(tn))**2)/2)
    !time_k = exp(-((tn-tnum_2)**2)*tkoef) !��� ����������� ����� �� �������
   do n=0,num_-1
   do k=0,num_-1
	 !tix=((n-num_2)**2+(k-num_2)**2)*koef 
	 tix=((n-num_2)**2+(k-num_2)**2)*(grid.size_x)**2/(2.*(num_-1)**2)
     beam.e(k,n,tn)=cmplx(time_k*exp(-tix),0.0d0, complex_kind)
!	 if (beam.abs_r0/=0) then 
	 beam.e(k,n,tn)=beam.e(k,n,tn)*exp(i*tix*koef_r) 
!	 end if
     !beam.e(k,n,tn)=0
	 !if (2*tix<1) then 
	 ! beam.e(k,n,tn)=1
	 !end if
   end do
   end do
  end do
 end subroutine init_gauss

end subroutine initiation

subroutine add_focus(num_,tnum_)
 use a_parameters
! use a_types
 use a_variables
! use a_medium

   integer(integer_kind) num_, tnum_,num_2,n,k,tn
   real(real_kind) :: tix, koef, koef_r, time_k

   num_2=num_/2		!�������� ����c� ���������� �� ���� num_2
!  num_2=num_/2-0.5 !�������� ����c� ���������� ���������� ����� ������ num_2 � num_2-1
   koef=(grid.h**2)/2 
  
  if (beam.abs_r1/=0) then 
    koef_r=beam.abs_k0*(beam.abs_a0**2)/beam.abs_r1
  else
    koef_r=0
  end if  
  
  do tn=0,tnum_-1
   do n=0,num_-1
   do k=0,num_-1
	 tix=((n-num_2)**2+(k-num_2)**2)*koef 
	   beam.e(k,n,tn)=beam.e(k,n,tn)*exp(i*tix*koef_r) 
   end do
   end do
  end do


end subroutine add_focus
