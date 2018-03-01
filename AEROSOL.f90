!  AEROSOLE.f90 
!
!  FUNCTIONS:
!	AEROSOL      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: AEROSOL
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

program aerosole

  use a_variables
  use a_parameters
  use a_functions
  !use a_field_oper
  use a_difraction
  use a_dispersion
  use a_plasma
  use sp_calc_and_write
  use dflib
  

  !use dfmt
    
  implicit none 

  real(real_kind) :: time_begin,time_end,new,kn
  real(real_kind) :: time_1,time_2
  real(real_kind) :: delta_z_tmp,rand_value, x, y, zk,tix, dz0
  integer(integer_kind) :: local1, local2, break_count, index_z_point,k,num_,num_2,n,m
  integer(integer_kind) :: x_point,y_point,t_point
 !real(real_kind) :: zk_aerosol
  integer(integer_kind) :: droplet_number,m_local,Droplet_x,Droplet_y,k_loc,m_loc,counter_loc,k11,m11
  integer(integer_kind) :: t_index,kt,kt1
  real(real_kind) :: tk,tk1
  character(len=255) :: string_name
  character(len=255) :: str1,str2,str3,str4,str5,str6,str7,str8,str9,str10,str11, strT, strNT
  integer(integer_kind) :: exitresult
  real(real_kind) :: tmp_real,tmp_realR,tmp_realF,tmp_realF2,tmp_realF10, tmp_after, Energy_Init, plasma_real, max0
!  complex :: drop(1:5)
!  real :: drop_z(1:5)
  !integer(integer_kind) :: curr_ae_screen
  integer(integer_kind) :: drop_type, lastopenskin, time_layer1, time_layer2
  !CONJG(X) - X комплексное сопряжение с Х;
  !CMPLX(X1,X2) - преобразование к комплексному числу; АВS(Х) - IXI модуль Х;
  !CABS(X) - модуль комплексного числа;
  !AIMAG(X) - мнимая часть комплексного числа Х;
  !REAL(X) - вещественная часть комплексного Х;
  
  
  
!Задали переменные
!******************************************************


!OpenMP Test

  integer :: myid, nthreads
  integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
  integer :: MY_NUM_OF_THREADS
  
  MY_NUM_OF_THREADS=4
  
  call OMP_SET_NUM_THREADS(MY_NUM_OF_THREADS)
  !$OMP PARALLEL default(none) private(myid) &
  !$OMP shared(nthreads)
  ! Determine the number of threads and their id

  myid = OMP_GET_THREAD_NUM()
  nthreads = OMP_GET_NUM_THREADS();
  !$OMP BARRIER
  if (myid==0) then 
   write(*,'(a)') 'OpenMP TEST'
   write(*,'(a,i1)') 'Number of Threads = ', nthreads
   write(*,'(a)') '*******************************'
  end if 
  !$OMP END PARALLEL

  exitresult=0
  !call AerosolInitDialog(exitresult)
  select case(exitresult)
  case(0)    
	  write(*,'(a)') 'Calculation initiated'
	  write(*,'(a)') '*******************************'
	  write(*,*)
  case(1) 
	  write(*,'(a)') 'Calculation initiated with Spectrum Output'
	  write(*,'(a)') '*******************************'
	  write(*,*)
      write_spectrum_flag=1
!      call exit
  end select

  call CreateDir
  
  call cpu_time(time_begin)

  call initiation      !Задание начальных условий
  
  call load_scattering_mass

  call MediumCalculation !Вычисляет переменные Medium
  call MediumAlloc      !Выделение памяти под массивы среды
  call TabFuncCreation
  
  !do n=1,5	
  !   drop_z(n)=drop_z(n)*media.delta_z_aerosol-1.e-7
  !end do 
 
  call random_seed()

  write(*,'(a40, e8.2)'), ' -- Average droplet number in skin   -- ', media.avr_skin_drop_number
  write(*,*)

  z=0
  counter_loc=1
  break_count=0

  zk_aerosol=media.start_z_aerosol  !положение текущего аэрозольного экрана (распределены по трассе равномерно, через расстояние media.delta_z_aerosol)
  zk_turb=media.start_z_turb        !положение текущего турбулентного экрана
  curr_ae_screen=0    !счетчик аэрозольных экранов
  init_turb_screen=8  !начальный номер файла с турбулентным экраном-1
  curr_turb_screen=0  !счетчик турбулентных экранов
  drop_type=0
  lastopenskin=0
  
  zk=media.build_field_start_z

  DecreseStepPlasma=0
  DecreseStepSelfFocus=0
  NeZXlen=NeZXlen_start
  FlZXYlen=FlZXYlen_start

  Flag_Write_NeZXFNum=0
  Counter_NeZXFNum=0 !число записанных слоев плазмы
  Counter_FlZXFNum=0 !число записанных профилей флюенса F(x)   
  
  local1=grid.num_x/(2*grid.size_x)
  !local1=grid.num_x

  !open(TMP_File,FILE=ResultDirName//'tmp.dat' ,access='sequential',status='unknown')  						!Начало расчета

  open(NIntFNum1,FILE=ResultDirName//'I(z)_int.dat' ,access='sequential',status='unknown')
  open(NFlFNum1,FILE=ResultDirName//'Fl(z).dat' ,access='sequential',status='unknown')
  open(NPlFNum1,FILE=ResultDirName//'Pl(z).dat' ,access='sequential',status='unknown') 						!Начало расчета
  open(NIntFNum2,FILE=ResultDirName//'Length.txt',access='sequential',status='unknown')
  if(PlasmaFlag/=0) then
     open(NeZXFNum,FILE=ResultDirName//'Ne_zx.grd',access='sequential',status='unknown')
     
     WRITE(NeZXFNum,110)
     WRITE(NeZXFNum,111) grid.num_x/2+1, 100   !сколько точек на сколько точек
     WRITE(NeZXFNum,112) -grid.h*grid.num_x/4,grid.h*grid.num_x/4 !от чего до чего по оси Х
     WRITE(NeZXFNum,112) NeZXlen_start,NeZXlen_start+0.1 !от чего до чего по оси Z
     WRITE(NeZXFNum,112) 0, 0.01 !от чего до чего по концентрации
     
     !Для флюенса
     !open(FluenceXFNum,FILE=ResultDirName//'Fluence_zx.grd',access='sequential',status='unknown')
     !WRITE(FluenceXFNum,110)
     !WRITE(FluenceXFNum,111) grid.num_x/2+1, 100   !сколько точек на сколько точек
     !WRITE(FluenceXFNum,112) -grid.h*grid.num_x/4,grid.h*grid.num_x/4 !от чего до чего по оси Х
     !WRITE(FluenceXFNum,112) NeZXlen_start,NeZXlen_start+0.1 !+0.83 !от чего до чего по оси Z
     !WRITE(FluenceXFNum,112) 0, 25.0 !от чего до чего по флюенсу
     
  end if
  
  if (media.droptest/=1) then
    if (media.droptest==-1) then
      open(FNDropTest,FILE=media.droptest_filename,access='sequential',status='unknown')
    end if
  else
    write(*,*), ' -- Type filename with drolet distribution (default - d:', media.droptest_filename,')' 
    write(*,*)
    read(*,*) str1
    if (str1=='d') then
      str1=media.droptest_filename
    end if
    open(FNDropTest,FILE=str1,access='sequential',status='unknown')
    if (eof(FNDropTest)) then
     write(*,'(a40, e8.2)'), 'File is empty' 
     call exit
    else !'(i,a1,i,a1,i,a1,i)'
     read (FNDropTest,*) lastopenskin,drop_type,Droplet_x,Droplet_y
     !write(*,*), 'read',lastopenskin,drop_type,Droplet_x,Droplet_y
     !read(*,*)
    end if 
  end if 
 
 
 
 beam.Energy=calc_energy(beam.e,grid.num_x,grid.time_number,1.d0,Max_Int,x_point,y_point,t_point)
 Energy_Init=beam.Energy
 call WriteParameters
 
 !Флюенс в самом начале до дифракции, аэрозоля и т.д.
  write(str7,'(a9,i,a4)'),'Fluence_x',0,'.dat'           
  call write_Fluence_x(1,FluenceGraphName//str7,Max_Fl)
  Fluence_Init=Max_Fl
  call write_Fluence_x(1,FluenceGraphName//str7,Max_Fl)
  
 if (write_spectrum_flag/=0) then
   open(NIntFNum3,FILE=IntGraphName//'init_sp.dat',access='sequential',status='unknown') 
   call sp_int_calc(beam.e,grid.num_x,t_point)
   close(NIntFNum3)
 end if 

     						!Начало расчета

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!********* Основной цикл по z **************************

  write(*,'(a)'), '---- START MAIN CYCLE ---------------------------------'
  write(*,*)

!		call write_xy_s(2,IntGraphName//'init.grd', t_point, grid.num_x/4)
!call exit

  do while (z<grid.size_z)

!call cpu_time(time_1)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!********** поглощение на краю расчетной сетки ***************
!local1 ширина полосы поглощения с края сетки

		!kn=2*pi*grid.num_x/(beam.abs_a0*grid.size_x)
		!local2=nint(    ( ( grid.dz*beam.abs_ld*kn )/sqrt( (2*pi/beam.abs_lambda0)**2-(kn)**2 ) ) / (grid.h*beam.abs_a0)   )
		!if (local2<local1) then
          !local1=local2
		!end if
		!if (z==0) then 
		   !write(*,'(a40, i4)'), ' -- Theoretical absorbtion zone width -- ',local1
		!end if

  		!if (local1<grid.num_x/64) then 
     	  !local1=grid.num_x/64
		  !write(*,'(a,i4)'), ' -- Low Limit of abs. zone width achieves, new value -- ',local1
		!end if 

  		!if (local1>grid.num_x/16) then 
     	  !local1=grid.num_x/16
		  !write(*,'(a,i4)'), ' -- High Limit of abs. zone width achieves, new value -- ',local1
		!end if 
		write(*,*),  '------------------------------------------------------------'

! вызов функции поглощения на краю расчетной сетки

!$OMP PARALLEL FIRSTPRIVATE (local1)
!$OMP DO SCHEDULE (DYNAMIC, 4)
	  do t_index=0,grid.time_number-1	
!!		call kill_limit_points(local1,t_index)
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!поступило предложение резать четвертью периода
        call kill_limit_points_sin(local1,t_index)
	  end do
!$OMP END DO  
!$OMP END PARALLEL

	
		!Проверка kill_limit_points
		         !call calc_energy(beam.e,grid.num_x,grid.time_number,Max_Int,x_point,y_point,t_point)
         	     !write(*,'(a35,i4,a1,i4,a1,i2,a8,F7.3)') ' *** Maximum intensity ***  Point=(',x_point,',',y_point,',', t_point,') Value=',Max_Int
				 !call write_xy(1,IntGraphName//'Field_TMP.grd', 0)
		         !exit

!time_2=time_1
!call cpu_time(time_1)
!write(*,'(a,F15.10,a)'), ' Time of "LIMIT POINT" operation is ',time_1 - time_2,' seconds'
!время прохода 0.05 ue



!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      if (z>media.end_z_aerosol) then
        zk_aerosol=grid.size_z+media.delta_z_aerosol
        media.end_z_aerosol=zk_aerosol
      end if 

      if (z>media.z_focus_lens) then
        media.z_focus_lens=2*grid.size_z
        call add_focus(grid.num_x, grid.time_number)
      end if 


! определение случайного числа раскидываемых капель, при равномерном 
! прохождении в аэрозоле расстояния delta_z_aerosol
	  if (beam.test<0) then
	    droplet_number=0
	  else if (beam.test==0) then
	    droplet_number=1	! тестовый случай
	  else
       if (z>zk_aerosol ) then
         curr_ae_screen=curr_ae_screen+1
         droplet_number=puasson_random(media.avr_skin_drop_number)
       !  if (droptest>0) then
       !      write(*,*), 'z>drop_z(droptest)',z,drop_z(droptest)

        !   if (z>drop_z(droptest)) then 
	     !    Droplet_x=int(grid.num_x/2+real(drop(droptest)))	!Сгенерим координаты по Оx, Oy   
		 !    Droplet_y=int(grid.num_x/2+aimag(drop(droptest)))	!Сгенерим координаты по Оx, Oy   
          !   write(*,*), 'Test droplet in the screen (x,y)=',Droplet_x,Droplet_y
	       !  droptest=droptest-1
          !   if ( (Droplet_x-grid.num_x/2)**2+(Droplet_y-grid.num_x/2)**2 <(media.area_x/grid.h)**2     ) then
		   !     do t_index=0,grid.time_number-1	
			!        call a_field_oper(Droplet_x,Droplet_y,t_index)	!Опишем изменения светового поля, вносимые каплей воды
		     !   end do
    		! !end if
		  ! end if     
         !end if
         zk_aerosol=media.start_z_aerosol+curr_ae_screen*media.delta_z_aerosol
         !zk_aerosol=zk_aerosol+media.delta_z_aerosol !сделано чтобы не накапливалась ошибка определения расстояния z
	   else
	      droplet_number=0
	   end if
	  end if
!    read(*), counter_loc


!********** Обработка поля в присутствии капель (интерференция рассеянной и прошедшей компонент поля) ***********
!time_2=time_1
!call cpu_time(time_1)

     if (media.droptest==1) then 
       if (lastopenskin==curr_ae_screen) then
        droplet_number=1
       else
        droplet_number=0
       end if 
     end if
     
     m_local=droplet_number
    
!((((((((((((((((((((((((((((((((absorb test 			  
!Tmp_drop_absorb=0
     
     do while (m_local>0)
       m_local=m_local-1
        
       if (media.droptest/=1) then
        
        if (media.poly==0) then 
	      !Монодисперсный аэрозоль
	      drop_type = media.droplet_type 
	    else
	      !Полидисперсный аэрозоль
	      call define_droplet_type(drop_type)	!Определим радиус частицы в соответствии с гамма-распределением
	    end if
  
        call random_number(rand_value)
		! возможно это и более правильно, но есть ощущение, что возиться с media.scat_num(2) не надо
		!		   Droplet_x=int(grid.num_x/2-media.area_x/grid.h)+int(rand_value*( 2*media.area_x/grid.size_x)*(grid.num_x-2*media.scat_num(2)) )	!Сгенерим координаты по Оx, Oy   
		!		   Droplet_x=Droplet_x+media.scat_num(2)
	    Droplet_x=int(grid.num_x/2+(2*rand_value-1)*media.area_x/grid.h)	!Сгенерим координаты по Оx, Oy   
        
		call random_number(rand_value)

		Droplet_y=int(grid.num_x/2+(2*rand_value-1)*media.area_x/grid.h)	!Сгенерим координаты по Оx, Oy   

        if(media.droptest==-1) then 
           write(FNDropTest,'(i3,a1,i3,a1,i4,a1,i4)') curr_ae_screen,' ',drop_type,' ',Droplet_x,' ',Droplet_y
        end if   
       end if
   if (beam.test>0) then 
        if ( (Droplet_x-grid.num_x/2)**2+(Droplet_y-grid.num_x/2)**2 <(media.area_x/grid.h)**2     ) then
		!		write(*,*), Droplet_x,Droplet_y
		!	    write(*,*), '*** Energy before droplet ***   ', calc_energy(beam.e,grid.num_x)

!$OMP PARALLEL FIRSTPRIVATE (drop_type,Droplet_x,Droplet_y)
!$OMP DO SCHEDULE (DYNAMIC, 4)
		    do t_index=0,grid.time_number-1	
!((((((((((((((((((((((((((((((((absorb test 			  
			  !Tmp_drop_absorb=Tmp_drop_absorb+(cabs(beam.e(Droplet_x,Droplet_y,t_index)))**2 !в поглощение
			  call a_field_oper(drop_type,Droplet_x,Droplet_y,t_index)	!Опишем изменения светового поля, вносимые каплей воды
 	    
		    end do
		     
!$OMP END DO  
!$OMP END PARALLEL
         !     write(*,*), '*** Energy after droplet ***   ', calc_energy(beam.e,grid.num_x)
		end if
		


        if(media.droptest==1) then 
          if (eof(FNDropTest)) then
            lastopenskin=-1
            write (*,*) 'Droplet distribution file ended'
          else
            read (FNDropTest,*) lastopenskin,drop_type,Droplet_x,Droplet_y
            ! write (*,*) drop_type,Droplet_x,Droplet_y 
          end if 
        end if
  
        if (lastopenskin==curr_ae_screen) then 
          droplet_number=droplet_number+1
          m_local=1
        end if
   end if    
	 end do
	 
	 		 !Запись флюенса сразу за слоем с каплями
             if (counter_loc==1)then
               !write(str7,'(a27)'),'Fluence_x_after_droplet.dat' !Преобразуем число в строку str7
               !call write_Fluence_x(1,CrossGraphName//str7)!Одномерный
               write(str8,'(a28)'),'Fluence_xy_after_droplet.grd'
               call write_Fluence_xy(1,FluenceGraphName//str8, int(1.d0*grid.num_x/grid.size_x))
              !call write_Fluence_xy(1,FluenceGraphName//str8,FluenceGraphName//str7 int(1.*grid.num_x/grid.size_x))!двумерный
             end if
             
             !Запись начального флюенса
              !if (counter_loc==0)then
                !write(str7,'(a27)'),'Fluence_x_init.dat' !Преобразуем число в строку str7
                !call write_Fluence_x(1,CrossGraphName//str7)!Одномерный
                !write(str8,'(a28)'),'Fluence_xy_init.grd'
                !call write_Fluence_xy(1,FluenceGraphName//str8,FluenceGraphName//str7, int(1.*grid.num_x/grid.size_x))!двумерный
              !end if

!((((((((((((((((((((((((((((((((absorb test 			  
!      beam.Energy=calc_energy(beam.e,grid.num_x,grid.time_number,1.*media.area_x/grid.size_x,Max_Int,x_point,y_point,t_point)	   
!      Tmp_drop_absorb=1-2*pi*(Tmp_drop_absorb/beam.Energy)*((scatter(drop_type).r/(grid.h*beam.abs_a0))**2)
!	  Tmp_drop_absorb=sqrt(Tmp_drop_absorb)
!     
!       do t_index=0,grid.time_number-1	
!		do k=0,grid.num_x-1
!         do n=0,grid.num_x-1
!          beam.e(n,k,t_index)=Tmp_drop_absorb*beam.e(n,k,t_index)
!         end do
!        end do
!       end do
	  
write(*,'(a, i5)'), ' Droplet quantity on the current skin = ', droplet_number
!********** Обработка поля в присутствии капель закончена***********************


!time_2=time_1
!call cpu_time(time_1)
!write(*,'(a,F15.10,a)'), ' Time of 2 operation is ',time_1 - time_2,' seconds'

!время прохода 
!beam.test<0 -- 0.00 ue
!beam.test=0 -- 0.00 ue
!beam.test>0 -- 300. ue при 33'000 каплях
!				310. ue при 33'600 каплях
!				316. ue при 34'000 каплях
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!тестовые выводы в файл; одна капля находится в центре слоя counter_loc
	!if (z==0)then
	if (counter_loc==counter_loc_drop)then
     if (beam.test==0) then
		write(*,*), 'Test droplet on the current skin'
		    do t_index=0,grid.time_number-1	
			  call a_field_oper(0,grid.num_x/2,grid.num_x/2,t_index)!Опишем изменения светового поля, вносимые каплей воды
		    end do
		  		  
       if (media.build_field==1) then
    	write(*,*), '*** Draw Field after droplet ***   '
	     beam.Energy=calc_energy(beam.e,grid.num_x,grid.time_number,1.d0,Max_Int,x_point,y_point,t_point)
		 call write_xy_s(1,IntGraphName//'Field_after_droplet.grd',0,grid.num_x/16)
       end if  
	 end if
    end if     

!time_2=time_1
!call cpu_time(time_1)
!write(*,'(a,F15.10,a)'), ' Time of 3 operation is ',time_1 - time_2,' seconds'
!время прохода 
!beam.test<0 -- 5.00 ue

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!ТУРБУЛЕНТНОСТЬ

if (media.turbulence==1) then
   strT=' '
   if (z>zk_turb) then
         strT='turb'
         curr_turb_screen=curr_turb_screen+1
         call turbulence(curr_turb_screen+init_turb_screen,grid.time_number)
   		 zk_turb=media.start_z_turb+curr_turb_screen*media.delta_z_turb
   end if
end if
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!ДИФРАКЦИЯ
     delta_z_tmp=grid.dz*beam.abs_ld*4*(pi**2)/(2*beam.abs_k0*(grid.abs_size_x**2))		! deta_z - набег фазы в параксиальном приближении для первой пространственной гармоники kx
		!	 write(*,*), '*** Energy before difraction ***   ', calc_energy(beam.e,grid.num_x)


!     difraction_thread(1).Arg
     
!     Handle1=CreateThread(0,0,do_difraction_thread,loc(difraction_thread(1).Arg),0,0)
!     Handle2=CreateThread(0,0,do_difraction_thread,loc(difraction_thread(2).Arg),0,0)
 	 
    call do_difraction(delta_z_tmp,grid.num_x,0,grid.time_number)!/2-1)
	 
!	 call do_difraction(delta_z_tmp,beam.e,grid.num_x,grid.time_0, grid.time_number)
		!	 write(*,*), '*** Energy after difraction ***   ', calc_energy(beam.e,grid.num_x)

	 !z=z+grid.dz
	 z_old=z !предыдущее расстояние
	 !z=nint(z/grid.dz+1)*grid.dz  !сделано чтобы не накапливалась ошибка определения расстояния z
	 z=z+grid.dz
	 length=z !Пройденная длина [Ld]

!		 call write_xy_s(1,IntGraphName//'Field_after_droplet.grd',0)
!time_2=time_1
!call cpu_time(time_1)
!write(*,'(a,F15.10,a)'), ' Time of DIFRACTION is ',time_1 - time_2,' seconds'
!время прохода 5.2 ue при 2048
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    beam.Energy=calc_energy(beam.e,grid.num_x,grid.time_number,1.d0,Max_Int,x_point,y_point,t_point)
	write(*,'(a,e9.2)') ' *** Energy *** ', beam.Energy !beam.abs_i0
		!	write(*,*) x_point,y_point

	write(*,'(a35,i4,a1,i4,a1,i4,a8,F7.3)') ' *** Maximum intensity ***  Point=(',x_point,',',y_point,',', t_point,') Value=',Max_Int

!time_2=time_1
!call cpu_time(time_1)
!write(*,'(a,F15.10,a)'), ' Time of "calc_energy" operation is ',time_1 - time_2,' seconds'
!время прохода 0.2 ue при 2048
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!ДИСПЕРСИЯ
      if (media.dispersion==1) then
         call do_dispersion
      end if

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxtime_2=time_1

!САМОФОКУСИРОВКА

     if (media.self_focusing==1) then 
   		 call self_focusing(grid.dz,Max_Int,grid.time_number)
	 end if
!time_2=time_1
!call cpu_time(time_1)
!write(*,'(a,F15.10,a)'), ' Time of SELF-FOCUSING operation is ',time_1 - time_2,' seconds'
!read(*,*) str1
!время прохода 0.6 ue при 2048
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!ПЛАЗМА
    
    kt1=1
    tk1=-0.75
    kt=1
    tk=-0.75
   
  if(TimeFlag /= 0 .and. PlasmaFlag/=0) then
    do m=0,grid.num_x-1
       do n=0,grid.num_x-1
         Medium.Ne(n,m)=0.
         Medium.NeN(n,m)=0.
         Medium.NeO(n,m)=0.
       end do
    end do  
    
    do t_index=0,grid.time_number-1	
       
       call plasma(0,grid.num_x-1,t_index)
       call EnergyAbsorption(0,grid.num_x-1,t_index)  
   
  if (length>=NeZXlen) then 
  
      if (grid.current_time(t_index)>=tk .and. kt<3) then
       kt=kt+1
       
      select case(kt)
       case(1)    
	    tk=0
	   case(2)    
        tk=grid.current_time(grid.time_number-1)-1.e-5
       case default
        tk=100
	  end select
       !tk1=(0.25*kt1-1.)
       
       write(str2,'(a5,i,i,a4)'),'Ne_xy',t_index,Counter_NeZXFNum+1,'.grd' !Преобразуем число в строку str2
       
       !write(str3,'(a6,i,i,a4)'),'Int_x',t_index,Counter_NeZXFNum+1,'.grd'!Преобразуем число в строку str3
       !call write_x_s(1,IntGraphName_time//str3, t_index, int(1.*grid.num_x/grid.size_x))
       
       !Опишем поле и плазму в поперечном сечении пучка на разных слоях импульса
         
        
       open(NeTimeFNum,FILE=PlasmGraphName//str2,access='sequential',status='unknown')  

     
       WRITE(NeTimeFNum,110)
   
!       WRITE(NeTimeFNum,111) grid.num_x, grid.num_x   !сколько точек на сколько точек
!       WRITE(NeTimeFNum,112) -grid.num_x/2*grid.h,(grid.num_x/2-1)*grid.h !от чего до чего по оси Х
!       WRITE(NeTimeFNum,112) -grid.num_x/2*grid.h,(grid.num_x/2-1)*grid.h !от чего до чего по оси y
!       WRITE(NeTimeFNum,112) 0, 1. !от чего до чего по концентрации
!         do n=0,grid.num_x-1
!              write(NeTimeFNum,113) (Medium.Ne(m,n), m=0, grid.num_x-1)
!              !(m-grid.num_x/2)*grid.h,' ',(n-grid.num_x/2)*grid.h,' ',Medium.Ne(n,m)
!         end do

       num_2=int(1.*grid.num_x/grid.size_x)       !area_ !media.scat_num(media.droplet_type)
       WRITE (NeTimeFNum,111) 2*num_2+1, 2*num_2+1   !сколько точек на сколько точек
       WRITE (NeTimeFNum,112) -num_2*grid.h,num_2*grid.h !от чего до чего по оси Х
       WRITE (NeTimeFNum,112) -num_2*grid.h,num_2*grid.h !от чего до чего по оси Y
       WRITE (NeTimeFNum,112) 0, 0.003 !от чего до чего по оси Z
 
       do n=grid.num_x/2-num_2,grid.num_x/2+num_2
          write(NeTimeFNum,113) (Medium.Ne(k,n),k=grid.num_x/2-num_2,grid.num_x/2+num_2)
       end do

       close(NeTimeFNum)
      end if
      
       time_layer1=grid.time_number+10
       time_layer2=grid.time_number+10
       if (grid.current_time(t_index)> 0.5 .and. grid.current_time(t_index)< 0.51)then       
           time_layer1=t_index
       end if
       if (grid.current_time(t_index)> 0.25 .and. grid.current_time(t_index)< 0.26)then       
           time_layer2=t_index
       end if
       
       if (t_index == t_point .or. t_index == time_layer1 .or. t_index == time_layer2 .or. t_index == grid.time_0) then
           write(str11,'(a6,i,i,a4)'),'Int_x',t_index,counter_loc,'.dat'!Преобразуем число в строку str3
           !write(str4,'(a6,i,i,a4)'),'Re_x',t_index,counter_loc,'.dat'!Преобразуем число в строку str3
           !write(str5,'(a6,i,i,a4)'),'Im_x',t_index,counter_loc,'.dat'!Преобразуем число в строку str3
           call write_x_s(1,IntGraphName_time//str11, t_index, grid.num_x/2) !int(1.*grid.num_x/grid.size_x))  !IntGraphName_time//str4,IntGraphName_time//str5 для мнимой и действительной частей
       end if
  
     end if  

 
      
    !Запись двумерных графиков на разных слоях импульса -0.75,-0.5,-0.25,0,0.25,0.5,0.75 длительности
    
 
   !if (grid.current_time(t_index)>=tk .and. kt<8) then
       !kt=kt+1
       !tk=(0.25*kt-1.)
    
     !if (counter_loc>=counter_loc_drop)then
          !write(str4,'(a4,i,i,a4)'),'Ne_x',t_index,counter_loc,'.dat' !Преобразуем число в строку str4
          !write(str5,'(a5,i,i,a4)'),'Int_x',t_index,counter_loc,'.dat' !Преобразуем число в строку str5
         
          !call write_x(1,CrossGraphName//str5,t_index)!Запись двухмерного графика интенсивност I(x)
         
     !Запись двухмерного графика плазмы Ne(x)     
        !open(NeXFNum,FILE=CrossGraphName//str4,access='sequential',status='unknown')
          !n=grid.num_x/2
          !num_2=n
          !do k=0,grid.num_x-1
             !write(NeXFNum,'(f,a,e)'),(k-num_2)*grid.h,' ',Medium.Ne(k,n)
          !end do
        !close(NeXFNum)
     !end if
     
   !end if
        
 end do
    

  !Запись графика флюенса после слоя с каплями
    !if (counter_loc>=counter_loc_drop)then
      !write(str2,'(a11,i4,a4)'),'Fluence_xy',counter_loc,'.grd'
      !call write_Fluence_xy(1,FluenceGraphName//str2,int(1.*grid.num_x/grid.size_x))!двумерный fl(x,y)
    !end if
    
     
    plasma_real=calc_plasma(Medium.Ne,grid.num_x,1.0/grid.size_x,Max_Pl) 
    phasePlasma =0.5*Medium.Rpl*Max_Pl !*dz - нелинейный набег фазы при плазменной нелинейности
                
    if(length>=NeZXlen) then
   
!*******************************************************************************
!              Наложение спектрального фильтра на верхние частоты
!*******************************************************************************	            
 	   if (SpFilterFlag == 1 .and. FilterWidth > 0)then
 	   
          beam.energyF1=calc_energy(beam.e,grid.num_x,grid.time_number,1.d0,Max_Int,x_point,y_point,t_point)
          !write(*,*)'energy before filter', beam.energyF1/Energy_Init	   
 	      open(NIntFNum4,FILE=IntGraphName//'sp_time.dat',access='sequential',status='unknown')
          open(NIntFNum5,FILE=IntGraphName//'sp_okno_time.dat',access='sequential',status='unknown')
 	      call sp_filter(grid.time_number)
 	      close(NIntFNum4)
 	      close(NIntFNum5)
 	      
 	      beam.energyF2=calc_energy(beam.e,grid.num_x,grid.time_number,1.d0,Max_Int,x_point,y_point,t_point)
 	      !write(*,*)'energy after filter', beam.energyF2/Energy_Init
 	       	      
 	      call energy_correction(grid.num_x, grid.time_number, beam.EnergyF1, beam.EnergyF2)
 	       
 	      !tmp_after=calc_energy(beam.e,grid.num_x,grid.time_number,1.,Max_Int,x_point,y_point,t_point)
 	      !write(*,*)'energy after correction', tmp_after/Energy_Init
 	      
 	   end if   
 	   
  !*******************************************************************************
              
    !write(str10,'(a5,i4,a4)'),'Spect',Counter_NeZXFNum+1,'.dat' !Преобразуем число в строку str10
      
     n=grid.num_x/2
     write(*,*), '*** Draw PLASMA Field ***   '
     !do m=grid.num_x/4, 3*grid.num_x/4-1 !0,grid.num_x-1
     !     write(NeZXFNum1,'(f,a,f,a,e)') length*100,' ',(m-grid.num_x/2)*grid.h,' ',Medium.Ne(m,n)
     !end do
       write(NeZXFNum,113) (Medium.Ne(m,n), m=grid.num_x/4, 3*grid.num_x/4)
       Flag_Write_NeZXFNum=1
       Counter_NeZXFNum=Counter_NeZXFNum+1 !число записанных слоев плазмы
         NeZXlen=NeZXlen_start+Counter_NeZXFNum*PdzAfterPlasma !NeZXlen+PdzAfterPlasma
     
     
   end if  
   
  !write(str7,'(a9,i,a4)'),'Fluence_x',Counter_FlZXFNum+1,'.dat'           
  !call write_Fluence_x(1,FluenceGraphName//str7,Max_Fl) !для записи F(x,z)
  !Counter_FlZXFNum=Counter_FlZXFNum+1
  
  write(*,*), '*** Write FLUENCE Field ***'
 if(length>=FlZXYlen) then
    
    write(str9,'(a11,i,a4)'),'Fluence_xy',Counter_FlZXYFNum+1,'.grd' !Преобразуем число в строку str2
    call write_Fluence_xy(1,FluenceGraphName//str9, int(1.d0*grid.num_x/grid.size_x))
    Flag_Write_FlZXYFNum=1
    Counter_FlZXYFNum=Counter_FlZXYFNum+1
    FlZXYlen=FlZXYlen_start+Counter_FlZXYFNum*PdzAfterPlasma
 end if        
end if


!time_2=time_1
!call cpu_time(time_1)
!write(*,'(a,F15.10,a)'), ' Time of PLASMA operation is ',time_1 - time_2,' seconds'
!read(*,*) str1
!время прохода 1.7 ue при 2048
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   

!time_2=time_1
!call cpu_time(time_1)
!write(*,'(a,F15.10,a)'), ' Time of EnergyAbsorption operation is ',time_1 - time_2,' seconds'
!read(*,*) str1
!время прохода 0.3 ue при 2048
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


 !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 !                        Адаптивное изменение шага при увеличении
 !                              нелинейного набега фазы
 !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx   
	!if (DecreseStepPlasma==1.or.DecreseStepSelfFocus==1) then
    !    DecreseStepPlasma=0
	!    DecreseStepSelfFocus=0
    ! if (grid.dz>(0.001/beam.abs_ld)) then
	!  grid.dz=grid.dz/2
	
	  dzDiffr=grid.abs_dz_diffraction_step/beam.abs_ld !начальный дифракционный шаг
	  
	  if(media.self_focusing==1) then
	  	  dzKerr=AdaptiveCriterium/phaseKerr !шаг в зависимости от фазы при самофокусировке
	  	  write(*,*) phaseKerr, dzKerr
	  else
	      dzKerr=1.e5	  
	  end if
	  
	  if(PlasmaFlag==1) then
	      dzPlasma=AdaptiveCriterium/phasePlasma !шаг в зависимости от плазменной нелинейности
	      write(*,*) phasePlasma
	  else
	      dzPlasma=1.e5
	  end if
	  
	  dz0=grid.dz	  
	  grid.dz=dmin1(dzDiffr, dzKerr, dzPlasma)
	  
	  if (grid.dz<dz0) then
          write(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
          write(*,'(a,e9.3)') 'Decrease step z, now skin = ', grid.dz
          write(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      end if
	! else
	!  grid.dz=0.0005/beam.abs_ld    
    !  write(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    !  write(*,'(a,e9.3)') 'Now skin = ', grid.dz
    !  write(*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	! end if
    !end if
		

	str1='-'

	if (media.build_field==1) then
		!  	write(*,*), '*** z>zk ***   ', z, zk, media.build_field_dz 
 
	  if (z>=zk) then
    	write(*,*), '*** Draw Field ***   '
    	
        if (counter_loc>=1000) then
		  !write(str1,'(a5,i4,a4)'),'Field',counter_loc,'.grd' !Преобразуем число в строку str1
		  write(str10,'(a5,i4,a4)'),'Spect',counter_loc,'.dat' !Преобразуем число в строку str10
		  !write(str3,'(a5,i4,a4)'),'Int_x',counter_loc,'.dat'!Преобразуем число в строку str3
          write(str3,'(a4,i4,a4)'),'Re_x',counter_loc,'.dat'
          write(str8,'(a4,i4,a4)'),'Im_x',counter_loc,'.dat'
          write(str6,'(a11,i4,a4)'),'Field_time',counter_loc,'.dat'
          
  		  !write(str2,'(a11,i4,a4)'),'Fluence_xy',counter_loc,'.grd' !Преобразуем число в строку str2
		else
	      if (counter_loc>=100) then 
		    !write(str1,'(a5,i3,a4)'),'Field',counter_loc,'.grd' !Преобразуем число в строку str1
		    write(str10,'(a5,i3,a4)'),'Spect',counter_loc,'.dat' !Преобразуем число в строку str10
		    !write(str3,'(a5,i3,a4)'),'Int_x',counter_loc,'.dat'!Преобразуем число в строку str3
   		    write(str3,'(a4,i3,a4)'),'Im_t',counter_loc,'.dat'
            write(str8,'(a4,i3,a4)'),'Re_t',counter_loc,'.dat'
            write(str6,'(a11,i3,a4)'),'Field_time',counter_loc,'.dat'
          !write(str2,'(a11,i3,a4)'),'Fluence_xy',counter_loc,'.grd' !Преобразуем число в строку str2
		  else
            if (counter_loc>=10) then 
		      !write(str1,'(a6,i2,a4)'),'Field0',counter_loc,'.grd' !Преобразуем число в строку str1
		      write(str10,'(a6,i2,a4)'),'Spect0',counter_loc,'.dat' !Преобразуем число в строку str10
		      !write(str3,'(a5,i2,a4)'),'Int_x',counter_loc,'.dat'!Преобразуем число в строку str3
   		      write(str3,'(a4,i2,a4)'),'Im_t',counter_loc,'.dat'
              write(str8,'(a4,i2,a4)'),'Re_t',counter_loc,'.dat'
              write(str6,'(a11,i2,a4)'),'Field_time',counter_loc,'.dat'
          !write(str2,'(a12,i2,a4)'),'Fluence_xy0',counter_loc,'.grd' !Преобразуем число в строку str2
		    else
		  	  !write(str1,'(a7,i1,a4)'),'Field00',counter_loc,'.grd' !Преобразуем число в строку str1
		  	  write(str10,'(a7,i1,a4)'),'Spect00',counter_loc,'.dat' !Преобразуем число в строку str10
		  	  !write(str3,'(a5,i1,a4)'),'Int_x',counter_loc,'.dat'!Преобразуем число в строку str3
  		  	  write(str3,'(a4,i1,a4)'),'Im_t',counter_loc,'.dat'
              write(str8,'(a4,i1,a4)'),'Re_t',counter_loc,'.dat'
              write(str6,'(a11,i1,a4)'),'Field_time',counter_loc,'.dat'
  		  !write(str2,'(a13,i1,a4)'),'Fluence_xy00',counter_loc,'.grd' !Преобразуем число в строку str2
		    end if
		  end if
		end if  	
		
		!call write_xy_s(2,IntGraphName//str1, t_point, grid.num_x/2-1)
		call write_t(1,IntGraphName//str6,IntGraphName//str3,IntGraphName//str8, grid.num_x/2-1, grid.num_x/2-1)!запись поперечного профиля интенсивности от времени для x,y = 0
		
		!call write_x_s(1, IntGraphName//str6,IntGraphName//str3,IntGraphName//str8, grid.time_0, grid.num_x/2-1) !CrossGraphName//str4, CrossGraphName//str5, grid.time_0, grid.num_x/2) !int(1.*grid.num_x/grid.size_x))
        !call write_Fluence_xy(1,FluenceGraphName//str2, int(1.*grid.num_x/grid.size_x))

        !if (beam.test<=0) then
		  !write(*,*), '*** Write Ox profile ***'
          !call write_x_s(1,IntGraphName//'Ox_'//str2,t_point, nint(1.*grid.num_x/grid.size_x))  ! запись поля на первом слое сразу за каплями перед дифракцией
          !call write_x_s(1,IntGraphName//'Ox_'//str2,t_point, grid.num_x/2)  ! запись поля на первом слое сразу за каплями перед дифракцией	    
	    !end if




		if (write_spectrum_flag/=0) then
                  open(NIntFNum3,FILE=IntGraphName//str10,access='sequential',status='unknown') 
                  call sp_int_calc(beam.e,grid.num_x,t_point)
                  close(NIntFNum3)
 	    end if 
 	    
 	          
 	    zk=z+media.build_field_dz
	  end if
    end if
    tmp_real=calc_energy(beam.e,grid.num_x,grid.time_number,(grid.size_x/2-1)/grid.size_x,max0,x_point,y_point,t_point) !энергия в области на один радиус пучка меньшей
    tmp_realR=calc_energy(beam.e,grid.num_x,grid.time_number,1.d0/grid.size_x,max0,x_point,y_point,t_point) !энергия в радиусе пучка
	tmp_realF=calc_energy(beam.e,grid.num_x,grid.time_number,1.d0/(30.d0*grid.size_x),max0,x_point,y_point,t_point) !энергия в размере филамента
	tmp_realF2=calc_energy(beam.e,grid.num_x,grid.time_number,1.d0/(15.d0*grid.size_x),max0,x_point,y_point,t_point)!энергия в 2х размере филамента
	tmp_realF10=calc_energy(beam.e,grid.num_x,grid.time_number,1.d0/(20.d0*grid.size_x),max0,x_point,y_point,t_point)!энергия в 2х размере филамента
		!	 write(*,'(i3,F15.7,F15.7,a)'),counter_loc,z,z*beam.abs_ld,'	'//str1

    str2=''
    if (Flag_Write_NeZXFNum>0) then
      Flag_Write_NeZXFNum=0
      write(str2,'(a9,i4)') '  plasma_', Counter_NeZXFNum
    end if
    
    str5=''
    if (Flag_Write_FlZXYFNum>0) then
      Flag_Write_FLZXYFNum=0
      write(str5,'(a10,i4)') '  fluence_', Counter_FlZXYFNum
    end if 
    
    write(NIntFNum2,'(i4,F15.7,F15.7,a15,i,F15.7,F15.7,F15.7,F15.7,F15.7,F15.7,F15.7,a15,a15, a15)') counter_loc,z,z*beam.abs_ld, '	'//str1, droplet_number, beam.Energy/Energy_Init, tmp_real/Energy_Init,tmp_realR/Energy_Init, tmp_realF2/Energy_Init,tmp_realF/Energy_Init,tmp_realF10/Energy_Init, plasma_real*10000, str5, str2, strT
      
	 
		 ! write(str11,'(a4,i3,a4)'),'Skin',counter_loc,'.dat' !Преобразуем число в строку str1
		 !call write_x(1,IntGraphName//str11 )!IntFName)  ! запись файла после дифракции    

	 write(NIntFNum1,'(F15.7,F15.7,F15.7,F15.7)'),z, grid.current_time(t_point),cdabs(beam.e(grid.num_x/2,grid.num_x/2,grid.time_0))**2,Max_Int
	 write(NFlFNum1,'(F15.7,F15.7)'),z, Max_Fl/Fluence_Init
	 write(NPlFNum1,'(F15.7,F15.7)'),z, Max_Pl
	 
!time_2=time_1
!call cpu_time(time_1)
!write(*,'(a,F15.10,a)'), ' Time of 6 operation is ',time_1 - time_2,' seconds'
!время прохода 1.1 ue при write_xy(8,... N=2048
!			   6.2 ue при write_xy(1,... N=2048
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

     time_1=time_end
	 call cpu_time(time_end)
     write(*,'(a,F5.2,a2,F7.2,a)'), ' Time of operation is ',time_end-time_1,' (',time_end - time_begin,') seconds'
	 write(*,'(a12,e9.3,a)'), ' Z_LENGTH = ', z,' diffraction length'

if (Max_Int>Max_Int_Const) then

    	write(*,*), '****** Draw  Final Field ********   '
       do t_index=0,grid.time_number-1	
         write(str1,'(a11,i2,a4)'),'Field_Last_',t_index,'.grd'
 		 call write_xy(4,IntGraphName//str1, t_index)
	   end do		

  	   call write_t(1,IntGraphName//'Field_Last_Time.dat', x_point, y_point);
	    
		z=grid.size_z
end if

	counter_loc=counter_loc+1  

end do
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 
  if(PlasmaFlag/=0) then
     close(NeZXFNum)
    !close(NeZXFNum1)
     !     close(FluenceXFNum)
  end if
    close(NIntFNum1)
    close(NFlFNum1)
    close(NPlFNum1)
  close(NFlFNum1)
  close(NIntFNum2)
  if (media.droptest/=0) close(FNDropTest)
  call MediumDealloc
!  if (media.poly==1) then
	 do m_local=2,10
	    write(*,*), 'Droplet number with radius', m_local, ' = ', media.different_types(m_local)
     end do   
!  end if
  call save_statistics  
    
  call beepqq(500,200)
  deallocate(beam.e)
!  deallocate(beam.fluence)

  write(*,*), 'Program ended with z=',length   
  stop
  
110    FORMAT('DSAA')
111    FORMAT(2I8)
112    FORMAT(2F15.7)
113    FORMAT(35F15.7)  

end program aerosole