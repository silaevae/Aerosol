module a_parameters
 
 implicit none

 integer, parameter :: integer_kind=4,real_kind=8,complex_kind=8
                       !real(8) is specified as REAL*8.

 integer(integer_kind) :: TimeFlag=0             ,& !=0 - ������������ ������; =1 - �������.
						  PlasmaFlag=0           ,& !���� ���������� ������������					   
                          EnergyAbsorptionFlag=0 ,& !���� ���������� ��� ���������
                          SpFilterFlag=0            !���� ��������� ������������� ���� 
 
 integer(integer_kind) :: DecreseStepPlasma ,& !���� ����������� ���� ������������ ���������� �������������
						  DecreseStepSelfFocus !���� ����������� ���� ������������ ����������� �������� �����
						  !adapt1=0 ,& !���� ����������� ���� ������������ ����������� �������� �����
						  !adapt2=0 ,& !���� ����������� ���� ������������ ���
						  !adapt3=0    !���� ����������� ���� ������������ ���������� �������������

 integer(integer_kind) :: Flag_Write_NeZXFNum  ,& !���� ������ ���������� ���� ������
                          Flag_Write_FlZXYFNum ,& !���� ������ ���������� �������
                          counter_loc_drop     ,& !����� ����, � ������� ��������� ���� �����
		                  curr_ae_screen       ,& !����� ������������ ����
		                  init_turb_screen     ,& !��������� ����� ����� � ������������ �������
		                  curr_turb_screen     ,& !����� ������������� ������
		                  Counter_NeZXFNum     ,& !������� ���������� ����� ������
		                  Counter_FlZXFNum     ,&
		                  Counter_FlZXYFNum       

 real(real_kind) ::		AdaptiveCriterium ,& !���� ���������� ����� ���� �� ��� �������� ���� ��������, �� dz ���������� � ��� ����
                        Max_Int_Const     ,& !����� ���������� �������������, ����� ��������������� ����
                        Isat              ,& !����� �������������, ������� � �������� ��������� ���������� ��� ��������� [��/�**2]
						NePlasma          ,& !������������ ����������, ��� ���������� ������� ���������, ��� ������������ ������ [Medium.N0]
                        NeZXlen           ,& !���������� �� ������� ������������ ���������� � ���������� �������
			            NeZXlen_start     ,& !���������� �� ������� �������� ������������ ���������� � ���������� �������
			            FlZXYlen          ,& !����������, �� ������� ������������ ���������� � ���������� ������� �������
			            FlZXYlen_start    ,& !����������, �� ������� �������� ������������ ���������� � ���������� ������� �������
						zk_aerosol        ,& !����������, �� ������� ��������� ����������� �����
						zk_turb           ,& !����������, �� ������� ��������� ������������ �����
						z                 ,& !����������
						z_old             ,& !���������� �� ���������� ����
						PdzAfterPlasma    ,& !��� ������ �������� ����� ����������� ������ [Ld]
						FilterWidth       ,& !������ ����������� ������������� ����
						IonSpeed=1           !������������ ��� ���������� ������ ��������� (=1 - ����� �� ��������)

!���������
! integer(integer_kind),parameter :: maxScPoint     =1000 !���������� ��������� � ���������� ������� ���������


 real(real_kind),&
 parameter :: pi  = 3.141592653589793238462643d0 ,& !����� ��
			  el  = 1.60219d-19   ,& !����� ��������� [��]
			  me  = 9.1095d-31    ,& !����� ��������� [��]
			  n0  = 1.000273d0      ,& !���������� ����������� �������
              n_drop = 1.33d0       ,& !����������� ����������� ������
			  eps0= 8.854d-12     ,& !������������� ���������� [�/�]
			  hpl = 1.055d-34     ,& !���������� ������ [��*�]
			  c0  = 2.9979d8      !,& !�������� ����� � ������� [�/�]
!			  Nl  = 2.68e25       ,& !����� �������� [1/�**3]
!			  hpl = 1.055e-34     ,& !���������� ������ [��*�]
!			  WO2 = 19.386499e-19 ,& !��������� ��������� ��������� [��]
!			  WN2 = 24.994164e-19 ,& !��������� ��������� ����� [��]
!              Om  = 20.6e12       ,& !����������� �������� [��]
!              Gam = 26.0e12       ,& !����������� �������� [��]
!              nl2 = 3.e-23        ,& !����������� ���������� ������������ ��� ������� [�**2/��]
!              g   = 0.5              !����. ����������� ����������� �������������� � ����������� ���������� ������������


 complex(complex_kind),parameter :: i=cmplx(0.0d0,1.0d0,complex_kind) !������ �������

 integer(integer_kind) :: write_spectrum_flag=0 ! ���� ������� ������� �������������

character(len=*),&
  parameter :: &
			   !*****************
			   !����� ����������:
			   !*****************
			   ResultDirName  = '_Graph\'                    ,& !��� �����������
			   IntGraphName   = '_Graph\intencity\'          ,& !��� �������� �������������
			   PlasmGraphName   = '_Graph\plasma\'          ,& !��� �������� ������
			   IntGraphName_time = '_Graph\intencity_time\' ,&   ! ��� ������������� ��������� �����
			   FluenceGraphName = '_Graph\fluence\' ,&
			   CrossGraphName   = '_Graph\crossection\',& !��� �������� ���������� �������
        !��� ���������� ������ ��������� ������ "\" �� �����, ����� ����� ��������� ������ ��� �������� ����������
		!���������� ���������� ������ �������������� ������� �� �
		!������ DirNames (��. ����)
			   !*****************
			   !����� ������:
			   !*****************
			   HistoryFName   = 'history.txt'        ,& !��� ��������������� �����
			   IntFName		  = 'int.dat'               !��� �������� ����������� ������������� I �� ��������� x � y
 
 integer(integer_kind),&
  parameter :: &
			   !*********************
			   !������ ��������� �/�:
			   !*********************
			   HistoryFNum   = 1  ,&  !��� ��������������� �����
			   NIntFNum		 = 2  ,&  !��� ����� ����� ������ �������� �������������
               NIntFNum1	 = 3  ,&
               NIntFNum2	 = 4  ,&
               NIntFNum3	 = 5  ,&
               FNDropTest	 = 6  ,&
               NIntFNumPlasma= 8  ,&  !��� ���������� ������� ��� ��������� ��������� � ..
	           NeZXFNum		 = 9  ,&    !��� ������ ������
  	           NeTimeFNum    = 10 ,&   !��� ������ ������� ������ � ����������� �� �������
  	           FluenceFNum   = 11 ,&   !��� ������ �������
  	           TMP_File      = 15 ,&
  	           NeXFNum       = 16 ,&  !��� ������ ������ � ���������� �������
  	           FluenceXFNum  = 17 ,&    !��� ������ ������� � ���������� ������� �� ����� ����������
  	           NImFNum	     = 18 ,&
  	           NIntFNum4	 = 19 ,&
  	           NIntFNum5	 = 20 ,&
  	           NFlFNum1      = 21 ,&
  	           NPlFNum1      = 22 ,&
  	           NeZXFNum1	 = 23 ,&
  	           NIntTFNum	 = 24 ,&
  	           NImTFNum	     = 25 ,&
  	           NReTFNum   	 = 26 ,&
  	           NReFNum   	 = 27 ,&
  	           NTurb         = 28    !��� ������ ������� ������� ��������������

end module a_parameters
