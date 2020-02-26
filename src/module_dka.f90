module module_dka
  use module_globals
  use module_autodiff
  implicit none

  private
  public DurandKernerAberth, init_dka, uninit_dka, find_roots_dka

  type DurandKernerAberth
     integer :: m !< m次方程式を解く
     real(8) :: x0 !< sum_[i=1->M] q(i)(x-x0)^i = 0を解く
     real(8),allocatable :: q(:) !< sum_[i=1->M] q(i)(x-x0)^i = 0を解く
     complex(8) :: beta !< 解はbetaを中心とする
     real(8) :: zeta !< 半径zetaの円の内部に存在する
     real(8),allocatable :: c(:) !< sum_[i=1->M] q(i)(x-x0)^i = sum_[i=1->M] c(i)(x-beta)^i
     complex(8),allocatable :: alpha(:) !< 答え
  end type DurandKernerAberth

  integer,pointer :: m !< m次方程式を解く
  real(8),pointer :: x0 !< sum_[i=1->M] q(i)(x-x0)^i = 0を解く
  real(8),pointer :: q(:) !< sum_[i=1->M] q(i)(x-x0)^i = 0を解く
  complex(8),pointer :: beta !< 解はbetaを中心とする
  real(8),pointer :: zeta !< 半径zetaの円の内部に存在する
  real(8),pointer :: c(:)!< sum_[i=1->M] q(i)(x-x0)^i = sum_[i=1->M] c(i)(x-beta)^i
  complex(8),pointer :: alpha(:) !< 答え
  
contains

  subroutine init_dka(dka_)
    type(DurandKernerAberth),target,intent(inout) :: dka_

    ! dka_%m=3; m=>dka_%m
    ! dka_%x0=2.d0; x0=>dka_%x0
    ! allocate(dka_%q(0:m)); dka_%q=(/-5.d0,4.d0,3.d0,2.d0/); q=>dka_%q

    ! dka_%m=7; m=>dka_%m
    ! dka_%x0=0.d0; x0=>dka_%x0
    ! allocate(dka_%q(0:m)); dka_%q=(/1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0/); q=>dka_%q
    
    ! dka_%m=7; m=>dka_%m
    ! dka_%x0=0.d0; x0=>dka_%x0
    ! allocate(dka_%q(0:m)); dka_%q=(/-1.d0,-3.d0,-3.d0,-1.d0,1.d0,3.d0,3.d0,1.d0/); q=>dka_%q

    dka_%m=7; m=>dka_%m
    dka_%x0=0.d0; x0=>dka_%x0
    allocate(dka_%q(0:m)); dka_%q=(/-5040.d0,13068.d0,-13132.d0,6769.d0,-1960.d0,322.d0,-28.d0,1.d0/); q=>dka_%q

    dka_%beta=0.d0 ;beta=>dka_%beta
    dka_%zeta=0.d0 ;zeta=>dka_%zeta
    allocate(dka_%c(0:m)); c=>dka_%c
    allocate(dka_%alpha(m)); alpha=>dka_%alpha

  end subroutine init_dka

  subroutine uninit_dka(dka_)
    type(DurandKernerAberth),target,intent(inout) :: dka_

    deallocate(dka_%alpha)
    deallocate(dka_%c)
    deallocate(dka_%q)
    nullify(m,x0,q,beta,zeta,alpha)

  end subroutine uninit_dka

  subroutine find_roots_dka
    type(Dcomplex8) :: Dbeta, Dg, Dc, Dalphaj
    integer :: i, j, k
    real(8) :: tmp
    complex(8) :: ztmp

    beta=-1.d0/dble(m)*(q(m-1)-m*x0*q(m))/q(m)
    
    Dbeta=ad(beta)
    Dg=q(m)
    do i=m-1,0,-1
       Dg=Dg*(Dbeta-x0)+q(i)
    end do
    Dc=d2t(Dg)
    
    zeta=-1.d0
    do i=1,m
       tmp=((m-1)*abs(Dc%f(m-i)/Dc%f(m)))**(1.d0/dble(i))
       if(tmp.gt.zeta) zeta=tmp
    end do

    open(1,file="range.res")
    do i=1,101
       tmp=2.d0*pi*i/100.d0
       ztmp=beta+zeta*exp(ione*tmp)
       write(1,*) real(ztmp), aimag(ztmp)
    end do
    close(1)


    do i=1,m
       alpha(i)=beta+zeta*exp((2.d0*dble(i-1)*pi+1.5d0)*ione/dble(m))
       write(10+i,*) real(alpha(i)), aimag(alpha(i))
    end do
    
    do k=1,1000
       do j=1,m
          Dalphaj=ad(alpha(j))
          Dg=q(m)
          do i=m-1,0,-1
             Dg=Dg*(Dalphaj-x0)+q(i)
          end do
          alpha(j)=alpha(j)-Dg%f(0)/Dg%f(1)*0.2d0
          write(10+j,*) real(alpha(j)),aimag(alpha(j))
       end do
    end do
    do i=1,m
       write(*,*) alpha(i)
    end do
    
  end subroutine find_roots_dka
  
end module module_dka
