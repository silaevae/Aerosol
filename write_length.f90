!**************************************************************************
!
!	Записывает соответствующую номеру слоя координату по Oz
!						 
!**************************************************************************
subroutine write_length(z_n,z_loc,string)
 use a_parameters
 use a_variables

  
  integer(integer_kind) :: z_n 
  real(real_kind) :: z_loc
 character(len=*) :: string
 
 open(NIntFNum,FILE=string,access='sequential',status='unknown')  
! endfile(NIntFNum)
  
  write(NIntFNum,*),z_n ,' ', z_loc

 close(NIntFNum)

end subroutine write_length
