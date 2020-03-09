!> Global variables
module module_globals
  implicit none

  ! parameters
  complex(8),parameter :: zero=cmplx(0.d0,0.d0,kind(1.d0)) !< 0
  complex(8),parameter :: ione=cmplx(0.d0,1.d0,kind(1.d0)) !< I
  complex(8),parameter :: one=cmplx(1.d0,0.d0,kind(1.d0)) !< 1
  real(8),parameter :: tiny=epsilon(1.d0)
  real(8),parameter :: pi=acos(-1.d0) !< Pi

  ! AutoDiff
  integer,parameter :: ntaylor=12 !< 最大微分回数
  
end module module_globals
