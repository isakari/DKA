!> main
program main
  use module_globals
  use module_autodiff
  use module_dka
  implicit none

  type(iCj) :: choose(0:ntaylor)
  type(DurandKernerAberth) :: dka
  
  call init_autodiff(choose)
  call init_dka(dka)
  call find_roots_dka
  call uninit_dka(dka)
  call uninit_autodiff(choose)
  
end program main
