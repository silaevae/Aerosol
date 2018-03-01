!subroutine scattering(droplet_radius,r,sc_ampl,sc_ph)
!! r-���������� � �������� �����
!! droplet_radius - ������ ����� � ���
!use a_parameters
!use a_types
!use a_variables
!use a_functions
!
!implicit none
!
!real(real_kind) :: r,sc_ampl,sc_ph
!integer(integer_kind) :: droplet_radius
!
!real(real_kind) :: teta_min,teta_find,teta_current !,phi0
!integer(integer_kind) :: index_loc,index_find, delta_i
!
!   r=r*beam.abs_a0		!�������� r � ������� ������� ��
!                        !beam.abs_a0 - ��������� ������ ����� �� ��� X [�]
!   teta_find=atan2(r,grid.abs_dz_step) !�������� �������� �������� �������
!                                       !grid.abs_dz_step - ���������� ��� ������� ��������� [�]
!   teta_min=1                  !������������ ������ � ����������� ����������
!   
!   !���������������� ����� (������ :))
!   index_loc=1000 !maxScPoint ! ������ �������
!   index_find=(index_loc+1)/2
!   delta_i=(index_find+1)/2
!
!   if (teta_find>scattering_mass(droplet_radius).scat(index_loc-1,0)) then
!     index_find=index_loc-1
!     delta_i=0
!   end if 
!   do while (delta_i>1)
!     if (teta_find>scattering_mass(droplet_radius).scat(index_find,0)) then 
!		index_find=index_find+delta_i
!	 else
!		index_find=index_find-delta_i
!	 end if
!	 if (index_find>index_loc-1) then 
!	   index_find=index_loc-1
!	 end if
!	 if (index_find<0) then 
!	   index_find=0
!	 end if
!	 delta_i=(delta_i+1)/2
!   end do
!
!   if (delta_i==0) then
!     sc_ampl=0
!     sc_ph=0
!   else
!	 if (index_find>1) then 
!       if (index_find<index_loc-2) then 
!   	     delta_i=index_find
!	   else
!	     delta_i=index_loc-3
!       end if
!	 else 
!	   delta_i=2
!	 end if
!
!	 do index_loc=delta_i-2,delta_i+2
!      teta_current=scattering_mass(droplet_radius).scat(index_loc,0); !���������� �� r
!      if (abs(teta_current-teta_find)<teta_min) then
!         teta_min=abs(teta_current-teta_find)
!         index_find=index_loc
!      end if
!     end do
!    
!	 sc_ampl=media.scatter*scattering_mass(droplet_radius).scat(index_find,1)
!     sc_ph=scattering_mass(droplet_radius).scat(index_find,2)
!   end if
!
!!   write(*,*) 'index_find+2', scattering_mass(droplet_radius,index_find+2,0)-teta_find
!!   write(*,*) 'index_find+1', scattering_mass(droplet_radius,index_find+1,0)-teta_find
!!   write(*,*) '0', scattering_mass(droplet_radius,index_find,0)-teta_find
!!   write(*,*) '-1', scattering_mass(droplet_radius,index_find-1,0)-teta_find
!!   write(*,*) '-2', scattering_mass(droplet_radius,index_find-2,0)-teta_find
!!   do index_loc=0,999             !����� ����������� �������� ���� � �����
!!      teta_current=scattering_mass(droplet_radius,index_loc,0); !���������� �� r
!!      if (abs(teta_current-teta_find)<teta_min) then
!!         teta_min=abs(teta_current-teta_find)
!!         index_find=index_loc
!!      end if
!!   end do
!!   write(*,*) 'delta_i', delta_i , 'index_find', index_find
!!   write(*,*) 'grid.abs_dz_step', grid.abs_dz_step , 'grid.h', grid.h
!!   write(*,*) 'grid.h*beam.abs_a0*media.scat_num(media.droplet_type)', beam.abs_a0*grid.h*media.scat_num(media.droplet_type)
!!   write(*,*) 'r', r , 'droplet_radius', droplet_radius
!!   write(*,*) 'teta_find', teta_find , 'index_find', index_find
!!   write(*,*) 'sc_ampl', sc_ampl , 'sc_ph', sc_ph
!!	read(*), index_loc
!
!end subroutine scattering

subroutine define_droplet_type(kind)

use a_parameters
use a_types
use a_variables
use a_functions

integer(integer_kind) :: kind		!������ ����� � ���
real(real_kind) :: random_loc
   call random_number(random_loc)	!��������� �������� ������������ ���������� �� ������� [0,,1]
   
   if (random_loc>=0.00 .and. random_loc<=0.087) then !���� ������� ��� ������ �����-������������� ����� � ��������� ��������, ��
   kind=2	!������ ����� = 2 ���
   end if

   if (random_loc>0.087 .and. random_loc<=0.279) then !���� ������� ��� ������ �����-������������� ����� � ��������� ��������, ��
   kind=3	!������ ����� = 3 ���
   end if

   if (random_loc>0.279 .and. random_loc<=0.514) then !���� ������� ��� ������ �����-������������� ����� � ��������� ��������, ��
   kind=4	!������ ����� = 4 ���
   end if

   if (random_loc>0.514 .and. random_loc<=0.718) then !���� ������� ��� ������ �����-������������� ����� � ��������� ��������, ��
   kind=5	!������ ����� = 5 ���
   end if

   if (random_loc>0.718 .and. random_loc<=0.855) then !���� ������� ��� ������ �����-������������� ����� � ��������� ��������, ��
   kind=6	!������ ����� = 6 ���
   end if

   if (random_loc>0.855 .and. random_loc<=0.932) then !���� ������� ��� ������ �����-������������� ����� � ��������� ��������, ��
   kind=7	!������ ����� = 7 ���
   end if
   
   if (random_loc>0.932 .and. random_loc<=0.970) then !���� ������� ��� ������ �����-������������� ����� � ��������� ��������, ��
   kind=8	!������ ����� = 8 ���
   end if

   if (random_loc>0.970 .and. random_loc<=0.988) then !���� ������� ��� ������ �����-������������� ����� � ��������� ��������, ��
   kind=9	!������ ����� = 9 ���
   end if

   if (random_loc>0.988 .and. random_loc<=1.000) then !���� ������� ��� ������ �����-������������� ����� � ��������� ��������, ��
   kind=10	!������ ����� = 10 ���
   end if


end subroutine define_droplet_type
