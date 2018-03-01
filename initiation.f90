!////////////////////////////////////////////////////////////////////////// 
!				 Задание начальных параметров		 
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
 

 beam.pulse_duration = 140.d-15 !140.e-15 !половина длительности импульса e-1 по интенсивности [s]

 
 grid.time_number = 500 !число раcсматриваемых временных слоев
 grid.time_start = -2.5 !-2.0 !начало рассматриваемого интервала времени [длительность импульса] 
 grid.time_end = 2.0 !конец рассматриваемого интервала времени [длительность импульса]
 
 allocate(grid.current_time(0:grid.time_number)) !для неоднородной сетки по времени
 allocate(grid.delta_time(0:grid.time_number)) 
 if (grid.time_number>1) then 
   TimeFlag=1 !=0 - стационарная задача; =1 - импульс.
   !grid.delta_time = (grid.time_end-grid.time_start)/(grid.time_number-1) !интервал времени
   grid.dt = (grid.time_end-grid.time_start)/(grid.time_number-1) !интервал времени при равномерном шаге
   time_step_param=4. !отношение начального и минимального шага
   time_step_param_1=sqrt(time_step_param) !отношение соседних шагов
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
   TimeFlag=0 !=0 - стационарная задача; =1 - импульс.
   grid.delta_time(0)=0
   grid.dt=0 
 end if 
 
 beam.abs_a0 = 2.5d-3 !1.0e-2 !Начальный радиус пучка по оси X [м]
 beam.abs_r0 = 0 !-0.48     !Кривизна волнового фронта [м]
 beam.abs_lambda0 = 0.8d-6 !Длина волны излучения   [м]
 beam.abs_e0 = 1.d0     !Начальная амплитуда     [В/м]
 beam.abs_i0 = 0.5*eps0*sqrt(n0)*(beam.abs_e0**2)  !Начальная интенсивность [Вт/м**2]
 
 beam.r = 50 !1.0e-1 !50 !30 ! Отношение мощности пучка к критической (параметр самофокусировки)

 beam.abs_i0 = beam.r*Air_Pcr/(pi*(beam.abs_a0**2))    !Начальная интенсивность [Вт/м**2]
 !write(*,*) 'beam.abs_i0',beam.abs_i0

 beam.TheorEnergy = pi*sqrt(pi)*(beam.abs_a0**2)*beam.pulse_duration*beam.abs_i0
! Pulse.AbsEnergy   = 0.
 beam.TheorPower  = pi*(beam.abs_a0**2)*beam.abs_i0

! beam.I0 = 
 beam.RcrAB =3.77 !для кругового гауссова коллимированного пучка
  
 AdaptiveCriterium= pi/15.d0 !Если нелинейный набег фазы за шаг превысит этот параметр, то dz уменьшится
 Isat=1.d14/beam.abs_i0  !Порог интенсивности, начиная с которого считается поглощение при ионизации [Вт/м**2]
 Max_Int_Const=1000 !порог превышения интенсивности, когда останавливается счет
 !*********** временно
 !TimeFlag=0 


 PlasmaFlag=1           !Флаг плазменной нелинейности
 EnergyAbsorptionFlag=1  !Флаг поглощения при ионизации
 NePlasma=1.d-7          !Концентрация электронов, при достижении которой считается, что образовалась плазма [Medium.N0]

 beam.test=1		!=0, если только один фазовый экран с 1-й каплей в центре
                    !=-1 капель вообще нет
 beam.linear_abs=0  !=1, тогда считается линейное ослабление                     

 media.droptest=-1  !=0 капли генерятся случайно
                    !=-1 капли генерятся случайно, запись расположения в файл (номер экрана, размер, x, y)
                    !=1 положение и размер капель задается из файла
                    
 media.turbulence=0 !=1 турбулентность есть, строятся фазовые экраны
                    !=0 турбулентности нет  
 SpFilterFlag=0    !Флаг наложения спектрального окна
 FilterWidth=0.5   !Ширина краев окна, по которым сглаживается спектр
 
 media.droptest_filename='ScatFunc\droplet_position.txt'   !'_Graph\droplet_position.txt' !имя файла куда кидаются данные о разбросе капель                    
 
 media.abs_drop_concentration = 10.d6 !концентрация частиц [м-3] 
 media.abs_drop_size = 15.d-6 ! размер частиц [м]	в случае монодисперсного аэрозоля
 media.area_x = 3  !радиус области раскидывания капель [радиус пучка]
 counter_loc_drop = 250 !номер слоя, в центре которого находится капля, если beam.test=0
 
!****************** Полидисперсный аэрозоль **************************************
 media.poly=0			!1-моделируется полидисперсный аэрозоль, 0-монодисперсный
 media.droplet_type=nint(media.abs_drop_size*1000000) !По умолчанию для монодисперсного аэрозоля
!Количество частиц, отличающихся размерами, неизменно - 9 типов частиц (2mkm - 10mkm)
!*********************************************************************************

!*************** Дисперсия********************************************************
 media.dispersion=0	!1-дисперсия есть, 0-нет
!*********************************************************************************

!*************** Самофокусировка *************************************************
 media.self_focusing=1	!1-самофокусировка есть, 0-нет
!*********************************************************************************

 grid.num_x = 1024 !число точек по Ox (= Oy)
 grid.size_x = 8 !полный поперечный размер [радиус пучка]
 grid.abs_size_z = 7 !полный продольный размер [м]
 grid.abs_dz_diffraction_step= 0.015 !продольный шаг расчета дифракции [м] 
 grid.num_z = nint(grid.abs_size_z/grid.abs_dz_diffraction_step) !число точек по Oz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !grid.abs_dz_step = 5.e-2 !2.e-3 !продольный шаг расчета рассеяния / ширина аэрозольного экрана[м] 

 beam.abs_k0 = 2*pi/beam.abs_lambda0               !Волновое число          [м-1]
 beam.abs_w0 = c0*beam.abs_k0                      !Частота излучения       [Гц]
 beam.Wph = hpl*beam.abs_w0               !Энергия одного фотона излучения [Дж]
  
 beam.abs_ld = beam.abs_k0*(beam.abs_a0**2)    !Дифракционная длина пучка [м]
 beam.abs_lds = (beam.pulse_duration)**2*(1.e15)**2/Air_k2 !дисперсионная длина [м] 
 
 !marb=0.367/sqrt((sqrt(beam.r)-0.852)**2-0.0219)*beam.abs_ld !расстояние самофокусировки [м]
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!********************** Задание длины области расчета ***************************
! не понял зачем это надо  
! grid.abs_limit_z=49.06 !49.06			    !верхняя граница расчета поля в длину [м]
! grid.limit_z=grid.abs_limit_z/beam.abs_ld	!верхняя граница расчета поля в длину [дифракционная длина пучка]
! grid.limit_z=grid.abs_size_z/beam.abs_ld	!верхняя граница расчета поля в длину [дифракционная длина пучка]
!********************************************************************************

 grid.size_z = grid.abs_size_z/beam.abs_ld !полный продольный размер [дифракционная длина пучка]

 grid.h = grid.size_x/(grid.num_x-1)     !поперечный шаг сетки [радиус пучка]
 grid.dz = grid.abs_dz_diffraction_step/beam.abs_ld    !продольный шаг сетки (дифракционный) [дифракционная длина пучка]
 grid.abs_size_x = grid.size_x*beam.abs_a0 !полный поперечный размер [м]
! grid.abs_size_z = grid.size_z*beam.abs_ld !полный продольный размер [м]

!Среднее число частиц в расчетном слое
! media.avr_skin_drop_number = media.abs_drop_concentration*4*((media.area_x*beam.abs_a0)**2)*grid.abs_size_z/grid.num_z 
! Среднее число частиц на дифракционной длине
! media.avr_skin_drop_number = media.abs_drop_concentration*4*((media.area_x*beam.abs_a0)**2)*grid.abs_size_z 

! расстояние, через которое происходит генерация аэрозольного экрана delta_z_aerosol
 media.delta_z_aerosol= grid.dz !0.009d0/beam.abs_ld ! равно начальному dz
 media.start_z_aerosol= -1.e-5 !0.819d0/beam.abs_ld  !положение первого аэрозольного экрана
 media.length_aerosol = 0.05d0/beam.abs_ld !ширина аэрозольного слоя [дифракционные длины]

 media.end_z_aerosol= 200./beam.abs_ld !media.start_z_aerosol+media.length_aerosol !конец аэрозольного слоя 
 !Среднее число частиц в расчетном слое
 media.avr_skin_drop_number = media.abs_drop_concentration*4*((media.area_x*beam.abs_a0)**2)*beam.abs_ld*media.delta_z_aerosol

!генерация фазовых экранов турбулентности
media.delta_z_turb= 0.5d0/beam.abs_ld !расстояние между экранами
media.start_z_turb= -1.e-5  !положение первого турбулентного экрана

! дополнительная фокусировка
 media.z_focus_lens= 200/beam.abs_ld          ! расположение фокусирующей линзы  
 beam.abs_r1 = 0 !-0.6     !дополнительная фокусировка [м]

!**************** Вывод в файл ***************************************************
 write_spectrum_flag=0 ! флаг записи спектра интенсивности
 
 media.build_field=1 !=1, если нужно на каждом i-м слое строить поле для Surfer'а; =0, если поля строить не надо
 media.build_field_dz=grid.dz !5*grid.dz ! указывается какое имеем i 
 media.build_field_start_z= grid.dz-1.e-5 ! расстояние построения первого графика


 media.build_field_dz=(1.-1.e-5)*media.build_field_dz !чтобы нейтрализовать ошибку округления

 
 NeZXlen_start= 100 !2.67               !Расстояние на котором начинает записываться информация о плазменных каналах[Ld]
 !FlZXYlen_start= 50*grid.dz            !Расстояние, на котором начинает записываться информация о флюенсе [Ld] 
 PdzAfterPlasma= 0.015d0/beam.abs_ld !media.delta_z_aerosol ! 0.00025 !0.0025   !Шаг записи графиков после образования плазмы [Ld]

!***************************************************

 call init_gauss(grid.num_x, grid.time_number)

! lenslet на входе  
! call phase_modulation(grid.size_x,grid.num_x-1,24,0.3,beam.e,2.1e-2,-1)

 
! subroutine phase_modulation(a_0,md1,N_d,k,W,R_f,sign)
!real a_0 - размер области [радиус пучка] = grid.size_x
!integer md1 - число поперечных шагов - 1; grid.num_x-1 
!integer N_d - поперечное (линейное) число линз (четное) = 32
!real k - безразмерный коэффициент сглаживания (h=k*d) отношение области сглаживания к размеру линзы 
!complex W(0:md1,0:md1) - поле
!real R_f - фокусное рассояние линз 2.1e-2
!integer sign - знак в экспоненте : -1 или +1
 
 contains

!выделяет память под beam.e(0:num_,0:num_) и загоняет в нее гаусс
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


! непонятно зачем
!   if (media.self_focusing==1) then
!      allocate(beam.e2(0:num_-1,0:num_-1))
!   end if

! num_2=num_/2		!максимум гаусcа приходится на узел num_2
  num_2=num_/2-0.5 !максимум гаусcа приходится посередине между узлами num_2 и num_2-1
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
    !time_k = exp(-((tn-tnum_2)**2)*tkoef) !для равномерной сетки по времени
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

   num_2=num_/2		!максимум гаусcа приходится на узел num_2
!  num_2=num_/2-0.5 !максимум гаусcа приходится посередине между узлами num_2 и num_2-1
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
