subroutine save_statistics

use a_parameters
use a_types
use a_variables
use a_functions

implicit none

integer(integer_kind)::j_loc
   open (1,file=ResultDirName//'droplet_distribution.txt',access='sequential',status='unknown')
   write(1,*), 'Радиус частицы        ', 'Количество частиц'
   do j_loc=2,10
      write(1,*), j_loc, '        ', media.different_types(j_loc)      
   end do
   close(1)
end subroutine save_statistics
