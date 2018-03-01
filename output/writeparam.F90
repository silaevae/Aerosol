!**************************************************************************
!                              WriteParameters
!
!                    Записывает начальные параметры в файл
!**************************************************************************
 subroutine WriteParameters

 use a_parameters
 use a_variables
 use a_medium
 use a_functions 
 
  implicit none

  integer(4) :: unit
  real(real_kind) :: size, realEnergy
  character(len=12) NowDate,NowTime

  unit=HistoryFNum
  
  open(unit,FILE=ResultDirName//'Parameters.txt',access='sequential',status='unknown')

  call ReturnDate(NowDate)
  call ReturnTime(NowTime)
       
     write(unit,'(A)') '                                     ПАРАМЕТРЫ ЗАДАЧИ'
         write(unit,*)
     write(unit,'(A)') ' Дата записи  '//NowDate
     write(unit,'(A)') ' Время записи '//NowTime
         write(unit,*)
         write(unit,*)
!**************************************************************************
!                                GRID
!**************************************************************************
    write(unit,'(A)') '===   ПАРАМЕТРЫ РАСЧЕТНОЙ СЕТКИ   ==='
    write(unit,'(A,I)') ' Количество равномерных шагов по X: nx=',Grid.num_x
!    write(unit,'(A,I)') ' Количество равномерных шагов по Y: ny=',Grid.NyConst+1
    write(unit,'(A,I)') ' Количество шагов по T: nt=',Grid.time_number
!    write(unit,'(A,I)') ' Число узлов на сетке с учетом нелинейности шага (по оси X): ',Grid.Nx+1
!    write(unit,'(A,I)') ' Число узлов на сетке с учетом нелинейности шага (по оси Y): ',Grid.Ny+1
!    write(unit,'(A,I)') ' Число сэкономленных узлов (по оси X): ', Grid.NxConst-Grid.Nx
!    write(unit,'(A,I)') ' Число сэкономленных узлов (по оси Y): ', Grid.NyConst-Grid.Ny
!    write(unit,'(A,f18.7)') ' Величина последнего шага по Ox [radx]: ',(Grid.x(Grid.Nx)-Grid.x(Grid.Nx-1))/Pulse.exc
!    write(unit,'(A,f18.7)') ' Величина последнего шага по Oy [radx]: ',(Grid.y(Grid.Ny)-Grid.y(Grid.Ny-1))/Pulse.exc
    write(unit,'(A,f18.7,A)') ' Пространственное разрешение сетки: Grid.h*beam.abs_a0=',Grid.h*beam.abs_a0*1.e6,' [мкм]'
    !write(unit,'(A,f18.7,A)') ' Временное разрешение сетки: grid.dt=',grid.dt,' [в длительностях импульса]'
          write(unit,*)
    write(unit,'(A,f18.7,A)') ' Начальный шаг по оси Z: grid.dz=', grid.dz,' [Ld]'
    write(unit,'(A,f18.7,A)') ' Начальный шаг по оси Z: grid.dz*beam.abs_ld=', grid.dz*beam.abs_ld,' [м]'
          write(unit,*)
          write(unit,*)
!**************************************************************************
!                                PULSE
!**************************************************************************
     write(unit,'(A)') '===   ПАРАМЕТРЫ ИМПУЛЬСА   ==='
 write(unit,'(A,f18.7,A)') ' beam.pulse_duration=',beam.pulse_duration*1.d15, '!Длительность импульса [фс]'    
 write(unit,'(A,f18.7,A)') ' beam.abs_a0=',beam.abs_a0*1.e2,' !Начальный радиус пучка [см]'
 write(unit,'(A,f18.7,A)') ' beam.abs_r0=',beam.abs_r0,' !Начальная фокусировка [м]'
 write(unit,'(A,f18.7,A)') ' Продольный масштаб: Ld=', beam.abs_ld,' [м]'
 write(unit,'(A,f18.7,A)') ' Дисперсионная длина: Lds=', beam.abs_lds,' [м]'
 write(unit,'(A,e18.7,A)') ' Частота излучения: w0=',beam.abs_w0,' [1/c]'
 write(unit,'(A,e18.7,A)') ' Волновое число: k0=',beam.abs_k0,' [1/м]'
 write(unit,'(A,e18.7,A)') ' Энергия одного фотона излучения: Wph=',beam.Wph,' [Дж]'

 select case(TimeFlag)
 case(0) !Стационарная задача
    write(unit,'(A,e18.7,A)') ' beam.TheorEnergy=', beam.TheorEnergy,' !Энергия импульса (рассчитывается теоретически для импульса не пропущенного через сетку) [Дж]'																	  
    write(unit,'(A,e18.7,A)') ' beam.TheorPower=',beam.TheorPower,' !Мощность пучка (рассчитывается теоретически) [Вт]'

!write(unit,'(A,e18.7,A)') ' beam.Energy=',beam.Energy*,' Энергия импульса [Дж]'
 case(1) !Импульс
 
    write(unit,'(A,e18.7,A)') ' beam.TheorEnergy=',beam.TheorEnergy,&
                              ' !Энергия импульса (рассчитывается теоретически для импульса не пропущенного через сетку) [Дж]'
    realEnergy = beam.Energy*(grid.h)**2*beam.pulse_duration*(beam.abs_a0**2)*beam.abs_i0
    write(unit,'(A,e18.7,A)') ' beam.Energy=',realEnergy,' !Энергия импульса [Дж]' !beam.Energy
    write(unit,'(A,e18.7,A)') ' beam.abs_i0=',beam.abs_i0/10000,' !Начальная пиковая интенсивность импульса (рассчитывается теоретически) [Вт/см2]'

 end select

         write(unit,*)
   if (media.self_focusing==1) then 
   write(unit,'(A,f18.7)') ' Количество критических мощностей в импульсе P/Pcr=',beam.r   !Medium.R/Pulse.RcrAB
   else
   write(unit,'(A)') ' Самофокусировки нет '
   end if 
         write(unit,*)
         write(unit,*)
!**************************************************************************
!                                MEDIUM
!**************************************************************************
     write(unit,'(A)') '===   ПАРАМЕТРЫ СРЕДЫ   ==='
   write(unit,'(A,f18.7)') ' Коэффициент керровской нелинейности: n2=',Medium.n2  
   write(unit,'(A,f18.7)') ' Параметр керровской нелинейности: R=',Medium.R
   write(unit,'(A,e18.7)') ' Параметр плазменной нелинейности: Rpl=',Medium.Rpl
Select Case(MediumFlag)
Case(1) !Воздух
   write(unit,'(A,e18.7)') ' Параметр нелинейности поглощения кислорода: RdO2=',Medium.RdO2
   write(unit,'(A,e18.7)') ' Параметр нелинейности поглощения азота: RdN2=',Medium.RdN2
   write(unit,'(A,f18.7)') ' Порядок многофотонного процесса при ионизации кислорода: lO2=',Medium.lO2
   write(unit,'(A,f18.7)') ' Порядок многофотонного процесса при ионизации азота: lN2=',Medium.lN2
Case(2) !Метанол
   write(unit,'(A,e18.7)') ' Параметр нелинейности поглощения: Rd=',Medium.Rd
   write(unit,'(A,f18.7)') ' Порядок многофотонного процесса при ионизации: l=',Medium.l
End Select
         write(unit,*)
         write(unit,*)
!**************************************************************************
!                                AEROSOL
!**************************************************************************
     write(unit,'(A)') '===   ПАРАМЕТРЫ АЭРОЗОЛЯ   ==='
     
   
   if (beam.linear_abs==1) then
    write(unit,'(A)') '   Линейное ослабление  '
   else if (beam.test==0)then
    write(unit,'(A,i)') ' Одна капля в слое', counter_loc_drop
   else if (beam.test==-1) then
    write(unit,'(A)') '   Капель нет  '
   end if
          
   write(unit,'(A,f18.7)') ' Концентрация аэрозоля [m-3]: n=',media.abs_drop_concentration
   write(unit,'(A,f18.7)') ' Радиус области раскидывания капель [радиус пучка]: area=',media.area_x
   write(unit,'(A,f18.7)') ' Расcтояние между аэрозольными экранами [m]: delta_z_aerosol=', media.delta_z_aerosol*beam.abs_ld
      
   if (media.poly==0) then	
     write(unit,'(A)') ' Тип аэрозоля: монодисперсный'
   else
     write(unit,'(A)') ' Тип аэрозоля: полидисперсный'
   end if
   write(unit,'(A,f18.7)') ' Среднее число частиц в расчетном слое: avr_skin_drop_number=', media.avr_skin_drop_number
   write(unit,'(A,f18.7)') ' Размер частиц [м] в случае монодисперсного аэрозоля: R=', media.abs_drop_size
         write(unit,*)
         write(unit,*)
!**************************************************************************
!                                NET
!**************************************************************************
!if(NetFlag /= 0) then
!      write(unit,'(A)') '===   ПАРАМЕТРЫ СЕТКИ   ==='
!  write(unit,'(A,f18.7,A)') ' Net.Zt=',Net.Zt,' !Расстояние Тальбо [м]'
!          write(unit,*)
!          write(unit,*)
!end if
!**************************************************************************
!                             RandomField
!**************************************************************************
!if(RandomFlag /= 0) then
!     write(unit,'(A)') '===   ПАРАМЕТРЫ СЛУЧАЙНОГО ПОЛЯ   ==='
! write(unit,'(A,f18.7,A)') ' RandomField.RInt=',RandomField.RInt,' !Дисперсия шума вычисленная численно'
! write(unit,'(A,f18.7,A)') ' RandomField.aver=',RandomField.aver,' !Среднее значение случайного поля'
!         write(unit,*)
!         write(unit,*)
!end if
!**************************************************************************
!                                PLOT
!**************************************************************************
!      write(unit,'(A)') '===   ПАРАМЕТРЫ ГРАФИКОВ   ==='
!    write(unit,'(A,e18.7)') ' Единица измерения флюенса: Pulse.FluenceOne=',Pulse.FluenceOne
!          write(unit,*)
!          write(unit,*)
!**************************************************************************
!                                SPECIAL
!**************************************************************************
 if(TimeFlag==0) then
   size=(Grid.num_x)*(Grid.num_x)*2.*complex_kind
 else
   size=(Grid.num_x)*(Grid.num_x)*(Grid.time_number)*2.*complex_kind
 end if
write(unit,'(A,e18.7,A)') ' Величина массива beam.E: ',size,'  [b]'
write(unit,'(A,e18.7,A)') ' Величина массива beam.E: ',size/1024,'  [Kb]'
write(unit,'(A,f18.7,A)') ' Величина массива beam.E: ',size/1024/1024,'  [Mb]'
write(unit,'(A,f18.7,A)') ' Величина массива beam.E: ',size/1024/1024/1024,'  [Gb]'
            write(unit,*)
            write(unit,*)

  close(unit)

 end subroutine WriteParameters
