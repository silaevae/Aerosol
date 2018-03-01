module a_variables

 use a_parameters
 use a_types

 implicit none

 real(real_kind) a, Max_Int, Max_Fl, Max_Pl, Fluence_Init
 real(real_kind) dzDiffr, dzKerr, dzPlasma, phaseKerr, phasePlasma
 real(real_kind) :: length=0. !Пройденная длина [Ld]
 
!**************** Задание массива рассеяния ****************************
 type(tscatter), pointer :: scatter(:)
!***********************************************************************

 type(tbeam) :: beam
 type(tmedia) :: media
 type(tgrid) :: grid
 
 complex, pointer :: sp(:,:)
 real, pointer :: sp_r(:),sp_i(:) 

end module a_variables