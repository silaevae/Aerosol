subroutine a_field_oper(dr_type_in,x_coord,y_coord,t_coord)

use a_parameters
use a_types
use a_variables
use a_functions


implicit none

!����������� ����������� ���� "beam.e" ��� �������� �������, ��������� ������
!�������� ������� �������� Droplet_x,Droplet_y - ���������� �� ���� Ox � Oy ��������������
	!integer :: dr_type
	integer(integer_kind) :: dr_type_in, x_coord,y_coord,t_coord
	integer(integer_kind) dr_type
    character(len=255) :: str1
	real(real_kind) :: sum1,sum2,alpha0,beta !,alphah,phi_loc
	complex (complex_kind) :: Amp0
!	real(real_kind) :: var_h
!,var_h2,c_re, c_im
!	real(real_kind) :: tangens_k_re,tangens_k_im,tangens_m_re,tangens_m_im
!	real(real_kind) :: alpha_h,alpha_0,amplitude_loc,energy_tough,beta

	integer(integer_kind) :: k,m,p, d_area, d_area_2 !, smoth
!	real(real_kind),pointer :: phase(:,:), dphase(:,:), e0(:,:), a_scat(:,:)
!	complex(complex_kind),pointer :: e_field(:,:)		!"��������� ����", �.� ���� � ������� �������������
    dr_type=dr_type_in
    if (dr_type==0) then
     if (media.poly==0) then 
	  !�������������� ��������
	  !write(*,*) 'media.droplet_type', media.droplet_type
	  !read(*,*)str1
	  
	  !dr_type = 0
	  dr_type = media.droplet_type 
	   
     else
	  !�������������� ��������
	  call define_droplet_type(dr_type)	!��������� ������ ������� � ������������ � �����-��������������
	 end if
	end if 

	d_area=scatter(dr_type).area

!    d_area=media.scat_num(media.droplet_type)

!			write (*,*) 'x_coord', x_coord, 'y_coord', y_coord
!			write (*,*) 'd_area', d_area
			!write (*,*) 'grid.num_x', grid.num_x

	
    if ((x_coord<=d_area).or.(y_coord<=d_area).or.(grid.num_x<=d_area+x_coord+1).or.(grid.num_x<=d_area+y_coord+1)) return
	
	media.different_types(dr_type)=media.different_types(dr_type)+1   
	
!	allocate(e_field(0:2*d_area,0:2*d_area))
!	allocate(phase(0:2*d_area,0:2*d_area))
!	allocate(dphase(0:2*d_area,0:2*d_area))
!	allocate(e0(0:2*d_area,0:2*d_area))
!    allocate(a_scat(0:2*d_area,0:2*d_area))
          
!    energy_tough=0            
    
!	do m=0,2*d_area			!�������������� ������� ��������/"������"
!     do k=0,2*d_area
!       e_field(k,m)=beam.e(x_coord-d_area+k,y_coord-d_area+m,t_coord)		!�������� ���� � "��������� ����", � ����� ������� ����� ��������� ������� �������������
!       
!	   phase(k,m)=atan2(imag(e_field(k,m)),real(e_field(k,m)))
!!	 
!!	   if (real(e_field(k,m))<0) then										!������ �������������� ���� ����
!!	      phase(k,m)=pi+atan(imag(e_field(k,m))/real(e_field(k,m)))		!� ������� �������������
!!       else
!!          if (real(e_field(k,m))/=0) then
!!		     phase(k,m)=atan( imag(e_field(k,m))/real(e_field(k,m)) )	
!!		  else
!!		     
!!			 if (imag(e_field(k,m))>0) then
!!			  phase(k,m)=pi/2
!!			 else
!!			  phase(k,m)=-pi/2
!!			 end if
!!		  end if
!!       end if															!.......................
!	   e0(k,m)=cabs(e_field(k,m))										!������ �������������� ��������� ����
! 
!!
!!���� ��������, ��� ��� ������ ����� �������� ����������, ������� �������������� ������
!!	   call scattering( media.droplet_type, sqrt( ((k-d_area)*grid.h)**2+((m-d_area)*grid.h)**2 ),a_scat(k,m),dphase(k,m) )	
!
!	   !�������� ��������� � ���� ������� ���������
!
!!	   var_h2=(k-d_area)**2+(m-d_area)**2
!!	   var_h=sqrt( var_h2 )
!!      if ( var_h<=d_area) then                         
!!	      energy_tough=energy_tough+a_scat(k,m)**2
!!      end if
!       
!	   end do
!    end do
!    
!!    beta=sqrt( 1.73E-10/(energy_tough*((grid.h*beam.abs_a0)**2)) )
!!   do k=0, 2*d_area
!!	   do m=0, 2*d_area
!!      var_h2=(k-d_area)**2+(m-d_area)**2
!!	   var_h=sqrt( var_h2 )
!!       if ( var_h<=d_area) then                         
!!	      a_scat(k,m)=a_scat(k,m)*beta   
!!       end if	   	   
!!	   end do
!!	end do

    select case (dr_type)
! beta ������ ���������� ������ ���������� ��������� ��� ������ ��������, ��� ��������������� �������
! ���� ������� �������, �� ������ ������� ��������� 2*pi*a^2
! ���� ���������, �� ���� ���������� ������� �������������
! �� ������ ������� ���������� �������� ��������� ��� ��������������� ������ ������ ���������,
! ��� �����-�� � ��-�� \Kirkhgof\scatter_kirkhgof.exe � ���� "�������� ���������" (���������� ���� ����� �������� ��� �������������)

!���� ����� ��� ���� �������� ���������������� ������������� � ��������������� �������, 
!����� ���������� ���� ����� ������ �������������� �-�� + ������� �� ���� ��������� ������
    case(2)
       !ro=15,7  sigma_poln~2.8 �������� ��� dz=2 � angle~0.3 =1,5007 (������ ���.)  ��� angle=1.4 =2,264
      beta=1.3
    case(4)
       !ro=31.4  sigma_poln~2.0 �������� ��� dz=2 � angle~0,2275=0,8519 (������ ���.) ��� angle=1.4 =1,7766
      beta=1.15 
    case(5)
       !ro=39.2 sigma_poln~2.1 �������� ��� dz=2 � angle~0,197 =0,895 (������ ���.) ��� angle=1.4 =1,826
      beta=1.205 
    case(6)
       !ro=47  sigma_poln~2.1 �������� ��� dz=2 � angle~0,16=0,985 (������ ���.) ��� angle=1.4 =1,907 angle~0.08 =0,8389 (������ ���.)
      beta=1.015 !- ��� ���� ���  !beta=1.16 !- ��� ������� ���
    case(8)
       !ro=62.8  sigma_poln~2.1 �������� ��� dz=2 � angle~0.104 =0,9905 (������ ���.) ��� angle=1.4 =2,015
      beta=1.11 
    case(10)
       !ro=78.5  sigma_poln~2.-2.1 �������� ��� dz=2 � angle~0.09 =0,8767 (������ ���.) ��� angle=1.4 =1,843
      beta=1.123 
    case(12)
       !ro=94.2  sigma_poln~2.-2.1 �������� ��� dz=2 � angle~0.075 =0,92305 (������ ���.) ��� angle=1.4 =__
      beta=1.077 
    case(14)
       !ro=110  sigma_poln~2.-2.1 �������� ��� dz=2 � angle~0,0615 =0,951573 (������ ���.) ��� angle=1.4 =__
      beta=1.048 
    case(15)
       !ro=118  sigma_poln~2.-2.1 �������� ��� dz=2 � angle~0.06 =0,9197 (������ ���.) ��� angle=1.4 =1,922
      beta=1.081 !��� dz=2mm 2-0,919
      !beta=1.083 !0.083 !��� dz=3mm 
      !���� ��������, ��� ���� ��������� ������ ����������, �� ������? 
      !���� ���������������� ������ ���������� ����, �� ���� ����� 0,97614; ���������� � ���� ���� 0,4445
      !����� beta=1,063 ��� ���� ��� �������
    case default
      beta=1.1
    end select 
!      beta=1
!   Sum1=Sum1-beta*2*pi*(cabs(Amp0)*scatter(dr_type).r/(grid.h*beam.abs_a0))**2 !� ����������
!0.16 ��� 6 ���
!0.08 ��� 15 ���   
!write(*,*) beta

    Sum1=0				!������ ������������ ������������ � ������������������ ����
    Sum2=0
	Amp0=beam.e(x_coord,y_coord,t_coord)
	
	d_area_2=d_area**2
	
	do m=-d_area, d_area
  	   p=m**2
       do k=-d_area, d_area
!       dphase(k,m)=dphase(k,m) !+phase(d_area,d_area)
!	   var_h=sqrt( 1.0*(k**2+m**2) )	   
       if ( (p+k**2)<=d_area_2) then                         
!		 write (*,*) k,m
         !Sum2=Sum2+(e0(k,m))**2+(a_scat(k,m)*e0(d_area,d_area))**2+2*e0(k,m)*a_scat(k,m)*e0(d_area,d_area)*cos(phase(k,m)-dphase(k,m))
!         Sum1=Sum1+(e0(k,m))**2
         Sum1=Sum1+(cdabs( beam.e(x_coord+k,y_coord+m,t_coord) ))**2
	
!if ((k.eq.0).and.(m.eq.0)) then 
!scatter(dr_type).e(k,m)=0         
!else
!	scatter(dr_type).e(k,m)=0
!end if 

		 Sum2=Sum2+(cdabs( beam.e(x_coord+k,y_coord+m,t_coord) + Amp0*scatter(dr_type).e(k,m) ))**2    
		 
		 !!!!! ����� � ����
		 !!!!! ��� ����� ������� ���������� ������ ������� ��������� Amp0 �� �������������� �������
		 !!!!! ����������� ��������� � ���������� ������ ������� � ��������� ��������, ��� 
		 !!!!! ������ ������������ �����
		 !!!!! ������, ���������� �� ���� ����� ���������� ������� � ���� ��� ��� ������� ���� ������
		     
	   end if
	   end do
    end do

! open(NIntFNum,FILE=IntGraphName//'tmp_a.grd',access='sequential',status='unknown')  
! WRITE (NIntFNum,110)
! 110    FORMAT('DSAA')
! WRITE (NIntFNum,'(2I8)') 2*d_area+1, 2*d_area+1   !������� ����� �� ������� �����
! WRITE (NIntFNum,'(2F15.7)') -d_area*grid.h,d_area*grid.h !�� ���� �� ���� �� ��� �
! WRITE (NIntFNum,'(2F15.7)') -d_area*grid.h,d_area*grid.h !�� ���� �� ���� �� ��� Y
! WRITE (NIntFNum,'(2F15.7)') 0, cabs(scatter(dr_type).e(0,0)) !�� ���� �� ���� �� ��� Z
! do m=-d_area, d_area
!    write(NIntFNum,'(35F15.7)') (cabs(scatter(dr_type).e(k,m)),k=-d_area, d_area)
!  ����  
!    write(NIntFNum,'(35F15.7)') (atan2(AIMAG(scatter(dr_type).e(k,m)),REAL(scatter(dr_type).e(k,m))),k=-d_area, d_area)
!  end do
! close(NIntFNum)
    

!********** TEST ******************
! beta=2

   Sum1=Sum1-beta*pi*(cdabs(Amp0)*scatter(dr_type).r/(grid.h*beam.abs_a0))**2 !� ����������

!!	write(*,*), Sum1,Sum2
!!    if (media.abs_drop_size == 6e-6) then
!   alpha0=sum1
!	   Sum1=Sum1-0.1*pi*(Amp0*media.abs_drop_size/(grid.h*beam.abs_a0))**2
!       Sum1=Sum1-0,16*2*pi*(Amp0*scatter(dr_type).r/(grid.h*beam.abs_a0))**2  !0.16 ��� 6 ���
!	   Sum1=Sum1-0.084*2*pi*(Amp0*scatter(dr_type).r/(grid.h*beam.abs_a0))**2 !0.08 ��� 15 ���   

!!! 0.1 - ��������� � ��������, ���� ������� ��������� ������ �� ���� ������� ���������
!!    end if
!	!� ��������� ��� ������� ������ ����������� �����������, ������������ ����������� ������� � ����������� � ������� ��������� �����. ���������
!!	if (media.abs_drop_size == 10e-6) then
!!	   Sum1=Sum1-0.208*pi*(e0(d_area,d_area)*media.abs_drop_size/(grid.h*beam.abs_a0))**2
!!    end if
!write(*,*) 'Sum-sigma*A2/Sum', sum1/alpha0

	alpha0=sqrt(Sum1/Sum2)
!	write(*,*),'alpha',alpha0
!	write(*,*),'Amp0',Amp0
!   write(*,*),'dr_type',dr_type

!!    energy_old=0


!����������� �� 'smoth' ������� ������
!   ����������� ��������, �.�. �������� ���������� �������������� ������, �������� ���
!   ������� ������� � ������� ���������, 
!   ������������ �������� ���������� ���� �� ������� ~3% ��� 15 ��� � <1% ��� 6
!   �������� alpha0 ������ ���������� �� 1 ������ ��� 1%
!    smoth=1
!    write(*,*) 'Sum1/Sum2', Sum1/Sum2

!    Sum1=0
!    Sum2=0
    do m=-d_area, d_area
  	   p=m**2
       do k=-d_area, d_area
!	   var_h2=d_area
!	   var_h=sqrt( (k-var_h2)**2+(m-var_h2)**2 )
!	   if ( (k**2+m**2)<=(d_area-smoth)**2 ) then                         
       if ( (p+k**2)<=d_area_2 ) then                         
!         Sum1=Sum1+(cabs( beam.e(x_coord+k,y_coord+m,t_coord) ))**2
		   beam.e(x_coord+k,y_coord+m,t_coord)=alpha0*( beam.e(x_coord+k,y_coord+m,t_coord) + Amp0*scatter(dr_type).e(k,m) )	      
!         Sum2=Sum2+(cabs( beam.e(x_coord+k,y_coord+m,t_coord) ))**2
!	   else
!  	     if ( (k**2+m**2)<=d_area**2 ) then                         
!		   var_h=(1.0*d_area-sqrt(1.0*(k**2+m**2)))/smoth
!		   beam.e(x_coord+k,y_coord+m,t_coord)=(1+var_h*(alpha0-1))*beam.e(x_coord+k,y_coord+m,t_coord) + Amp0*alpha0*var_h*scatter(dr_type).e(k,m)	      
!            write(*,*) 'field', cabs(Amp0*scatter(dr_type).e(k,m))/cabs(beam.e(x_coord+k,y_coord+m,t_coord))
!		 end if
	   end if

       end do
    end do		!������ ���� beam.e ��������

!    write(*,*) 'Sum2/Sum1 after', Sum2/Sum1
!    write(*,*) 'scatter(dr_type).e(k,m)', cabs(scatter(dr_type).e(0,0))



!!   write (*,*) 'Energy Old / Energy New',energy_old,energy_new 
!    write (*,*), alpha0  
   
!   deallocate(e_field)
!   deallocate(phase)
!   deallocate(dphase)
!   deallocate(e0)
!   deallocate(a_scat)          
   

end subroutine a_field_oper




subroutine kill_limit_points(points_number, t_coord)

use a_parameters
use a_types
use a_variables
use a_functions


implicit none

integer(integer_kind) :: k,m,points_number, t_coord, num_
real(real_kind) :: lam_,phi_loc,amplitude_loc
complex(complex_kind) :: tmp_cmpl

										!��������� �����, �� ������� ���������� ���������

  num_=grid.num_x-1
  lam_=pi/points_number

  do m=0,points_number
   tmp_cmpl=cmplx(0.5*(1.- cos(lam_*m) ),0.0d0, real_kind)
   do k=0,num_
	beam.e(k,m,t_coord)=beam.e(k,m,t_coord)*tmp_cmpl         
   end do
  end do
 
  do m=num_-points_number, num_
   tmp_cmpl=cmplx(0.5d0*(1.0d0- cos(lam_*(num_-m)) ),0.0d0, real_kind)
   do k=0,num_
	beam.e(k,m,t_coord)=beam.e(k,m,t_coord)*tmp_cmpl
   end do
  end do

  do m=0, num_
   do k=0,points_number
   	beam.e(k,m,t_coord)=beam.e(k,m,t_coord)*cmplx(0.5d0*(1.- cos(lam_*k)) ,0.)
    end do
   do k=num_-points_number, num_
   	beam.e(k,m,t_coord)=beam.e(k,m,t_coord)*cmplx(0.5d0*(1.- cos(lam_*(num_-k)) ),0.)
   end do
  end do

end subroutine kill_limit_points


!��������� ����������� ������ ������� ��������
subroutine kill_limit_points_sin(points_number, t_coord)
use a_variables
implicit none

integer(integer_kind) :: k,m,points_number, t_coord, num_
real(real_kind) :: lam_,phi_loc,amplitude_loc
complex(complex_kind) :: tmp_cmpl, pi2

										!��������� �����, �� ������� ���������� ���������
!��������� ����������� ������ ������� ��������
 num_=grid.num_x-1
 lam_=pi/points_number/2
 pi2=pi/2
 
  do m=0,points_number
   tmp_cmpl=cmplx(sin(lam_*m),0.)
   do k=0,num_
	beam.e(k,m,t_coord)=beam.e(k,m,t_coord)*tmp_cmpl         
   end do
  end do
 
  do m=num_-points_number, num_
   tmp_cmpl=cmplx(sin(lam_*(num_-m)),0.)
   do k=0,num_
	beam.e(k,m,t_coord)=beam.e(k,m,t_coord)*tmp_cmpl
   end do
  end do

  do m=0, num_
   do k=0,points_number
   	beam.e(k,m,t_coord)=beam.e(k,m,t_coord)*cmplx(sin(lam_*k),0.0d0,complex_kind)
    end do
   do k=num_-points_number, num_
   	beam.e(k,m,t_coord)=beam.e(k,m,t_coord)*cmplx(sin(lam_*(num_-k)),0.0d0,complex_kind)
   end do
  end do

end subroutine kill_limit_points_sin
