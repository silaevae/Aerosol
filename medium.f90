module a_medium
!************************************************************
! взял из модуля M_Medium Федорова В.
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

 integer(4) :: MediumFlag=1 !Флаг выбора среды(=1 - воздух, =2 - метанол)

 !=== Параметры воздуха ===
 real(real_kind),&
 parameter :: Air_n   = 1.000273      ,& !Показатель преломления воздуха
              !Air_n2  = 5.6e-23       ,& !Коэффициент керровской нелинейности для воздуха на длине волны 800нм [м**2/Вт]
              Air_n2  = 1.7d-23       ,& !Коэффициент керровской нелинейности для воздуха на длине волны 800нм [м**2/Вт]
              Air_k2  = 16.           ,& !Коэффициент, описывающий дисперсию групповой скорости, для воздуха на длине волны 800нм [фс**2/м]
              Air_N0  = 2.68d25       ,& !Число нейтральных молекул в воздухе [1/м**3]
              Air_WO2 = 19.386499d-19 ,& !Потенциал ионизации кислорода [Дж]
              Air_WN2 = 24.994164d-19 ,& !Потенциал ионизации азота [Дж]
              Air_g   = 0.5           ,& !Коэф. описывающий соотношение безынерционной и инерционной керровской нелинейности (=0 - только безынерционная, =1- только инерционная)
              Air_Pcr = 4.d9 ,& !6.e9          ,& !Критичская мощность в воздухе [Вт]
              Air_Ne0 = 0.               !Количество свободных электронов в отсутствие лазерного излучения [1/м**3]

 !=== Параметры метанола ===
 real(real_kind),&
 parameter :: Methanol_n   = 1.3286       ,& !Показатель преломления метанола
              Methanol_N0  = 1.5e28       ,& !Число нейтральных молекул в метаноле [1/м**3]
              Methanol_Wg  = 9.933578e-19 ,& !Ширина запрещенной зоны в метаноле [Дж]
              Methanol_g   = 0.           ,& !Коэф. описывающий соотношение безынерционной и инерционной керровской нелинейности (=0 - только безынерционная, =1- только инерционная)
              Methanol_Pcr = 8.e6         ,& !Критичская мощность в метаноле [Вт]
              Methanol_Sc  = 1.e-19       ,& !Сечение упругих столкновений [м**2]
              Methanol_Ne0 = 1.e16           !Количество электронов зоны проводимости в отсутствие лазерного излучения [1/м**3]

!**************************************************************************
!                               MEDIUM
!**************************************************************************
 type tMedium

  real(real_kind) :: n   ,& !Показатель преломления
                     n2  ,& !Коэффициент керровской нелинейности
                     N0  ,& !Число нейтральных атомов (молекул) среды
                     Wp0 ,& !Нормированная плазменная частота
                     R   ,& !Параметр керровской нелинейности
                     Rpl ,& !Параметр плазменной нелинейности
                     g   ,& !Коэф. описывающий соотношение безынерционной и инерционной керровской нелинейности (=0 - только безынерционная, =1- только инерционная)
                     Pcr ,&
                     Ne0    !Количество свободных электронов или электронов зоны проводимости в отсутствие лазерного излучения

  real(real_kind) :: Wc0,& !Нормированная частота столкновений
                     Kpl,& !Medium.Wc0/beam.abs_w0
                     Kbl   !Коэф. ударной ионизации

  real(real_kind),pointer :: Ne(:,:)     ,& !Концентрация свободных электронов или электронов зоны проводимости
                             Ne_old(:,:)    !Концентрация свободных электронов или электронов зоны проводимости с предыдущего шага по времени

  !--- Переменные для воздуха ---------------------------------------------
  real(real_kind) :: lO2  ,& !Порядок многофотонного процесса при ионизации кислорода
                     lN2  ,& !Порядок многофотонного процесса при ионизации азота
                     RdO2 ,& !Параметр нелинейности поглощения кислорода
                     RdN2    !Параметр нелинейности поглощения азота

  real(real_kind),pointer :: NeN(:,:)    ,& !Концентрация свободных электронов азота
                             NeN_old(:,:),& !Концентрация свободных электронов азота с предыдущего шага по времени
                             NeO(:,:)    ,& !Концентрация свободных электронов кислорода
                             NeO_old(:,:)   !Концентрация свободных электронов кислорода с предыдущего шага по времени

  !--- Переменные для метанола --------------------------------------------
  real(real_kind) :: Sc ,& !Сечение упругих столкновений
                     Wg ,& !Ширина запрещенной зоны
                     l  ,& !Порядок многофотонного процесса при ионизации
                     Rd    !Параметр нелинейности поглощения
  !------------------------------------------------------------------------

  integer(4) :: plasma=0 !Флаг образования плазмы (=1 - плазма образовалась)

 end type tMedium

 type(tMedium) :: Medium

!**************************************************************************
!                              PlasmaChannel
!**************************************************************************
! type tPlasmaChannel
!
!  real(real_kind) :: length,& !Длина плазменного канала
!                     RhoNe    !Плотность плазмы в плазменном канале
!
!  real(real_kind) :: x,y !Координаты центра тяжести плазменного канала
!
!  integer(4) :: NumBor !Размерность Border
!  type(tGridPoint),pointer :: Border(:) !Граница плазменного канала
!
!  integer(4) :: NumPoints !Размерность PCpoints
!  type(tGridPoint),pointer :: PCpoints(:) !Точки плазменного канала
!
! end type tPlasmaChannel
!
! integer(4) :: NumPlasmaChannel=0 !Количество плазменных каналов
!
! type(tPlasmaChannel),pointer :: PlasmaChannel(:)
!
! integer(4),pointer :: PCMap(:,:) !Карта плазменных каналов
!
!**************************************************************************
                                 CONTAINS
!**************************************************************************
!**************************************************************************
!                               MediumAlloc
!
!                    Выделяет память под массивы Medium
!**************************************************************************
 subroutine MediumAlloc
  integer(4),automatic :: k,n
    
  if(TimeFlag/=0 .and. PlasmaFlag/=0) then

    allocate(Medium.Ne(0:grid.num_x,0:grid.num_x))
    allocate(Medium.Ne_old(0:grid.num_x,0:grid.num_x))

    if(MediumFlag==1) then !Воздух
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

!Флаг анализа плазмы
!    if(PAflag /= 0) then
!      allocate(PCMap(0:grid.num_x,0:grid.num_x))
!    end if

  end if

 end subroutine MediumAlloc

!**************************************************************************
!                              MediumDealloc
!
!                Освобождает память занятую массивами Medium
!**************************************************************************
 subroutine MediumDealloc

  if(TimeFlag/=0 .and. PlasmaFlag/=0) then

    deallocate(Medium.Ne)
    deallocate(Medium.Ne_old)

    if(MediumFlag==1) then !Воздух
      deallocate(Medium.NeN,Medium.NeO)
      deallocate(Medium.NeN_old,Medium.NeO_old)
    end if


!Флаг анализа плазмы
!    if(PAflag /= 0) then
!      deallocate(PCMap)
!    end if

  end if

 end subroutine MediumDealloc

!**************************************************************************
!                            MediumCalculation
!
!                       Вычисляет переменные Medium
!**************************************************************************
 subroutine MediumCalculation

  !Задание параметров среды
  
  select case(MediumFlag)
  case(1) !=== Воздух =====================================================

    Medium.Pcr = Air_Pcr
    Medium.n   = Air_n
    Medium.n2  = beam.RcrAB*beam.abs_lambda0**2/(8.*pi*Medium.n*Medium.Pcr) !Air_n2
    Medium.N0  = Air_N0
    Medium.g   = Air_g
    Medium.Ne0 = Air_Ne0
    Medium.lO2 = Air_WO2 / beam.Wph !Порядок многофотонного процесса при ионизации кислорода
    Medium.lN2 = Air_WN2 / beam.Wph !Порядок многофотонного процесса при ионизации азота
    Medium.RdO2= 1./beam.abs_i0*beam.abs_ld*Medium.N0/beam.pulse_duration*Medium.lO2*beam.Wph !Параметр нелинейности поглощения кислорода
    Medium.RdN2= 1./beam.abs_i0*beam.abs_ld*Medium.N0/beam.pulse_duration*Medium.lN2*beam.Wph !Параметр нелинейности поглощения азота

    Medium.Wp0 = sqrt(el**2*Medium.N0/eps0/me) !Нормированная плазменная частота
    Medium.R   = 2.*beam.abs_k0**2/Medium.n*Medium.n2*beam.abs_a0*beam.abs_a0*beam.abs_i0 !Параметр керровской нелинейности
    Medium.Rpl = beam.abs_k0*beam.abs_ld/Medium.n**2*Medium.Wp0**2/beam.abs_w0**2 !Параметр плазменной нелинейности
    

  case(2) !=== Метанол ====================================================

    Medium.Pcr= Methanol_Pcr
    Medium.n  = Methanol_n
    Medium.n2 = beam.RcrAB*beam.abs_lambda0**2/(8.*pi*Medium.n*Medium.Pcr)


!   marb=0.367/sqrt((sqrt(Medium.R/beam.RcrAB)-0.852)**2-0.0219)
!write(unit,'(A,f18.7)') ' Количество критических мощностей в импульсе P/Pcr=',Medium.R/beam.RcrAB
   
    Medium.N0 = Methanol_N0
    Medium.g  = Methanol_g
    Medium.Ne0= Methanol_Ne0/Medium.N0
    Medium.Sc = Methanol_Sc
    Medium.Wg = Methanol_Wg
    Medium.l  = Medium.Wg / beam.Wph !Порядок многофотонного процесса при ионизации
    Medium.Rd = 1./beam.abs_i0*beam.abs_ld*Medium.N0/beam.pulse_duration*Medium.l*beam.Wph !Параметр нелинейности поглощения

    Medium.Wp0 = sqrt(el**2*Medium.N0/eps0/me) !Нормированная плазменная частота
    Medium.R   = 2.*beam.abs_k0**2/Medium.n*Medium.n2*beam.abs_a0**2*beam.abs_i0 !Параметр керровской нелинейности
    Medium.Rpl = beam.abs_k0**2*beam.abs_a0**2/Medium.n**2*Medium.Wp0**2/beam.abs_w0**2 !Параметр плазменной нелинейности

  end select
  !========================================================================

 end subroutine MediumCalculation

!!**************************************************************************
!!                           PlasmaAnalysis
!!
!!                         Анализирует плазму
!!**************************************************************************
! subroutine PlasmaAnalysis
!
!  NumPlasmaChannel=0
!
!  call FindPlasmaChannels !Ищем плазменные каналы
!
!  call PlasmaChannelsCalculation !Вычисляем параметры плазменных каналов
!  
! end subroutine PlasmaAnalysis
!
!
!!**************************************************************************
!!                            FindPlasmaChannels
!!
!!            Ищет узлы сетки принадлежащие плазменным каналам
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
!  !Инициализация карты плазменных каналов
!  do n=0,grid.num_x
!  do k=0,grid.num_x
!    if(Medium.Ne(k,n) >= NePlasma) then
!      PCMap(k,n)=-1 !Отмечаем точки где плазма есть
!    else
!      PCMap(k,n)=0  !Отмечаем точки где нет плазмы
!    end  if
!  end do
!  end do
!
!  finish=0
!
!  DO WHILE(finish==0)
!
!    !Ищем максимум не принадлежаций другим каналам
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
!    if(maximum /= 0.) then !Найден новый плазменный канал
!
!       if(NumPlasmaChannel/=0) then !Сохраняем старые плазменные каналы
!          allocate(PCold(NumPlasmaChannel))
!          PCold=PlasmaChannel
!          deallocate(PlasmaChannel)
!       end if
!
!       NumPlasmaChannel=NumPlasmaChannel+1 !Увеличиваем количество плазменных каналов на единицу
!       allocate(PlasmaChannel(NumPlasmaChannel))
!
!       if(allocated(PCold)) then
!          PlasmaChannel(1:NumPlasmaChannel-1)=PCold
!          deallocate(PCold)
!       end if
!
!       !Ищем границу нового плазменного канала
!       call FindPCBorder(MaxPoint.k,MaxPoint.n,PlasmaChannel(NumPlasmaChannel))
!
!       !Ищем координаты точек нового плазменного канала
!       call FindPCpoints(PlasmaChannel(NumPlasmaChannel))
!
!       !Отмечаем точки принадлежащие новому плазменному каналу
!       do i=1,PlasmaChannel(NumPlasmaChannel).NumPoints
!          PCMap(PlasmaChannel(NumPlasmaChannel).PCpoints(i).k,PlasmaChannel(NumPlasmaChannel).PCpoints(i).n)=1
!       end do
! 
!    else !Больше нет плазменных каналов
!
!       finish=1
!
!    end if
!
!  END DO
!
!  !Проверяем правильность нахождения плазменных каналов
!  i=0
!  do n=0,grid.num_x
!  do k=0,grid.num_x
!     if(PCMap(k,n)==1) i=i+1
!  end do
!  end do
!
!  if(i/=sum(PlasmaChannel.NumPoints)) then !Выводим сообщение в информационный файл
!     call WriteHistory('Ошибка в процедуре FindPlasmaChannel',&
!                       'Количество узлов сетки где образовалась'// &
!                       ' плазма не совпадает с количеством узлов принадлежащих плазменным каналам')
!  end if
!
! end subroutine FindPlasmaChannels
!
!!**************************************************************************
!!                           CheckNe(k,n,direction)
!!
!! Проверяет наличие плазмы в узле сетки соседним с (k,n) в направлении
!! задаваемом переменной direction. CheckNe=0 если в нем нет плазмы и =1
!! если плазма там есть.
!!
!! k,n - координаты узла сетки рядом с которым проверяется наличие плазмы
!! direction - направление поиска (1-верх, 2-право, 3-низ, 4-лево)
!!**************************************************************************
! function CheckNe(k,n,direction)
!
!  integer(4) :: k,n,direction,CheckNe
!
!  CheckNe=0
!
!  select case(direction)
!  case(1) !Вверх
!    if(Medium.Ne(k,n+1)>=NePlasma) CheckNe=1
!  case(2) !Вправо
!    if(Medium.Ne(k+1,n)>=NePlasma) CheckNe=1
!  case(3) !Вниз
!    if(Medium.Ne(k,n-1)>=NePlasma) CheckNe=1
!  case(4) !Влево
!    if(Medium.Ne(k-1,n)>=NePlasma) CheckNe=1 
!  end select
!
! end function CheckNe
!
!!**************************************************************************
!!                  MoveToDirection(k,n,direction)
!!
!!       Делает текущим узел сетки в направлении direction
!!
!! k,n - координаты узла сетки из которого перемещаемся
!! direction - направление движения (1-верх, 2-право, 3-низ, 4-лево)
!!**************************************************************************
! subroutine MoveToDirection(k,n,direction)
!
!  integer(4) :: k,n,direction
!
!  select case(direction) !Двигаемся в заданном направлении
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
!!                Находит границу плазменного канала
!!
!! nowk,nown - координаты узла сетки из которого начинается поиск
!! PlasmaChannel - плазменный канал чьи границы ищутся
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
!  !=== Доходим до границы =================================================
!  nstop=0
!  direction=1 !Выбираем в качестве начального направление движения движение вверх
!  do while(nstop==0)
!
!    checkresult=CheckNe(nowk,nown,direction) !Проверяем соседнюю (nowk,nown) точку в направлении direction
!
!    if(checkresult==1) then
!       call MoveToDirection(nowk,nown,direction) !Двигаемся в выбранном направлении
!    else
!       StartPoint.k=nowk !Запоминаем точку с которой начали обход границы
!       StartPoint.n=nown
!       bor=bor+1
!       Border(bor)=StartPoint
!       direction=2 !При достижении границы поворачиваем направо, чтобы граница всегда оставалась слева
!       nstop=1
!    end if
!
!  end do
!
!  !=== Обходим границу ====================================================
!  nstop=0
!  DO WHILE(nstop==0)
!
!    checkresult=CheckNe(nowk,nown,direction) !Проверяем соседнюю (nowk,nown) точку в направлении direction
!
!    if(checkresult==1) then !Следующая точка принадлежит каналу
!
!      call MoveToDirection(nowk,nown,direction) !Двигаемся в выбранном направлении
!
!      bor=bor+1
!      Border(bor).k=nowk
!      Border(bor).n=nown
!
!      !Смотрим что слева
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
!    else !Изменяем направление движения по часовой стрелке
!
!      direction=direction+1 !Изменили направления движения
!      if(direction==5) direction=1
!
!      checkresult=CheckNe(nowk,nown,direction) !Проверили соседнюю точку в новом направлении
!    
!      if(checkresult==1) then !Соседняя точка принадлежит плазменному каналу
!
!         call MoveToDirection(nowk,nown,direction) !Двигаемся в выбранном направлении
!
!         bor=bor+1
!         Border(bor).k=nowk
!         Border(bor).n=nown
!
!         direction=direction-1 !Возвращаемся к исходному направлению
!         if(direction==0) direction=4
!
!      end if
!
!    end if
!
!    !Если пришли в точку с которой начинали обход границы, то граница найдена
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
!!    Находит узлы сетки принадлежащие данному плазменному каналу
!!
!! PlasmaChannel - плазменный канал чьи узлы ищутся
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
!  !Для каждого у ищем точки лежащие между границами канала
!  do n=0,grid.num_x
!
!    posy=FindInArray(n,PlasmaChannel.Border(:).n,PlasmaChannel.NumBor)
!
!    if(posy/=-1) then !Существует точка с таким n принадлежащая границе
!
!       !Находим точки пересечения с границей для данного n
!       CrossingPoints=grid.num_x*grid.num_x
!       numcross=0
!       do i=1,PlasmaChannel.NumBor
!          if(PlasmaChannel.Border(i).n==n) then
!             numcross=numcross+1
!             CrossingPoints(numcross)=PlasmaChannel.Border(i).k
!          end if
!       end do
!
!       !Сортируем массив точек пересечения
!       count=numcross
!       call sortqq(loc(CrossingPoints(1:numcross)),count,SRT$INTEGER4)
!       if(count /= numcross) stop 'Sorting Error in procedure FindPCpoints'
!
!       do k=0,grid.num_x
!       
!          if(PCMap(k,n)/=-1) cycle !Не рассматриваем точки принадлежащие другим каналам и точки где нет плазмы
!
!          do i=1,numcross
!
!             if(k==CrossingPoints(i)) then !Точка принадлежит границе
!                pc=pc+1
!                PCpoints(pc).k=k
!                PCpoints(pc).n=n
!                exit
!             end if
!
!             if(i==numcross) exit !Дальше нет промежутков между точками пересечения
!
!             if(k>CrossingPoints(i) .and. k<CrossingPoints(i+1)) then !Точка лежит между точками пересечения с границей
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
!! Рассчитывает характеристики плазменных каналов
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
!    !Ищем координату центра тяжести плазменного канала
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
!    !Ищем плотность плазмы в плазменном канале
!    sum=0.
!    do j=1,PlasmaChannel(i).NumPoints
!      sum=sum + (Grid.x(k(j)+1)-Grid.x(k(j)-1))/2.* &
!                (Grid.y(n(j)+1)-Grid.y(n(j)-1))/2.* Medium.Ne(k(j),n(j))
!    end do
!
!    PlasmaChannel(i).RhoNe=sum*beam.abs_a0*beam.abs_a0*Medium.N0 !Переводим в [1/м]
!
!  END DO
!
! end subroutine PlasmaChannelsCalculation
!
end module a_medium
