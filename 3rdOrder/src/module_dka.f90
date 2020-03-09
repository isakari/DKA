module module_dka
  use module_globals
  use module_autodiff
  implicit none

  private
  public DurandKernerAberth, init_dka, uninit_dka, find_roots_dka

  type DurandKernerAberth
     integer :: m !< m次方程式を解く
     real(8) :: x0 !< sum_[i=1->M] q(i)(x-x0)^i = 0を解く
     complex(8),allocatable :: q(:) !< sum_[i=1->M] q(i)(x-x0)^i = 0を解く
     complex(8) :: beta !< 解はbetaを中心とする
     real(8) :: zeta !< 半径zetaの円の内部に存在する
     complex(8),allocatable :: c(:) !< sum_[i=1->M] q(i)(x-x0)^i = sum_[i=1->M] c(i)(x-beta)^i
     complex(8),allocatable :: alpha(:) !< 答え
  end type DurandKernerAberth

  integer,pointer :: m !< m次方程式を解く
  real(8),pointer :: x0 !< sum_[i=1->M] q(i)(x-x0)^i = 0を解く
  complex(8),pointer :: q(:) !< sum_[i=1->M] q(i)(x-x0)^i = 0を解く
  complex(8),pointer :: beta !< 解はbetaを中心とする
  real(8),pointer :: zeta !< 半径zetaの円の内部に存在する
  complex(8),pointer :: c(:)!< sum_[i=1->M] q(i)(x-x0)^i = sum_[i=1->M] c(i)(x-beta)^i
  complex(8),pointer :: alpha(:) !< 答え
  
contains

  subroutine init_dka(dka_)
    type(DurandKernerAberth),target,intent(inout) :: dka_

    ! dka_%m=3; m=>dka_%m
    ! dka_%x0=2.d0; x0=>dka_%x0
    ! allocate(dka_%q(0:m)); dka_%q=(/-5.d0,4.d0,3.d0,2.d0/); q=>dka_%q

!!$    dka_%m=7; m=>dka_%m
!!$    dka_%x0=0.d0; x0=>dka_%x0
!!$    allocate(dka_%q(0:m)); dka_%q=(/1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0/); q=>dka_%q
    
    ! dka_%m=7; m=>dka_%m
    ! dka_%x0=0.d0; x0=>dka_%x0
    ! allocate(dka_%q(0:m)); dka_%q=(/-1.d0,-3.d0,-3.d0,-1.d0,1.d0,3.d0,3.d0,1.d0/); q=>dka_%q

    ! dka_%m=7; m=>dka_%m
    ! dka_%x0=0.d0; x0=>dka_%x0
    ! allocate(dka_%q(0:m)); dka_%q=(/-5040.d0,13068.d0,-13132.d0,6769.d0,-1960.d0,322.d0,-28.d0,1.d0/); q=>dka_%q

    dka_%m=10; m=>dka_%m !10次方程式
    dka_%x0=0.d0; x0=>dka_%x0
    allocate(dka_%q(0:m)); q=>dka_%q
    q( 0)=-ione
    q( 1)=1.d0+3.d0*ione
    q( 2)=-3.d0-4.d0*ione
    q( 3)=4.d0+4.d0*ione
    q( 4)=-4.d0-2.d0*ione
    q( 5)=2.d0-2.d0*ione
    q( 6)=2.d0+4.d0*ione
    q( 7)=-4.d0-4*ione
    q( 8)=4.d0+3.d0*ione
    q( 9)=-3.d0-ione
    q(10)=one
    
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
    type(Dcomplex8) :: Dbeta, Dg, Dc, Dalphaj, Dzeta
    integer :: i, j, k, iflag(m)
    real(8) :: tmp, tmp2, brake
    complex(8) :: ztmp, zaccel(m)

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
    do k=1,10 !解の存在範囲を絞るためのNewton法
       Dzeta=ad(zeta)
       Dg=abs(Dc%f(m))*Dzeta**m
       do i=0,m-1
          Dg=Dg-abs(Dc%f(i))*Dzeta**i
       end do
       zeta=zeta-real(Dg%f(0))/real(Dg%f(1))
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
    
    tmp=tiny
    tmp2=tiny
    ! Newton法で sum_[0->npad] q_i(alp(j)-x0)^i=0を解く
    k=0
    iflag(:)=0
    brake=1.d0
    do
       if(k.gt.1000000) exit
       tmp2=tmp
       tmp=tiny
       zaccel(:)=zero
       do j=1,m
          do i=1,m
             if(i.ne.j) zaccel(j)=zaccel(j)+one/(alpha(j)-alpha(i))
          end do
       end do
       do j=1,m
          if(iflag(j).eq.1) cycle
          Dalphaj=ad(alpha(j))
          Dg=q(m)
          do i=m-1,0,-1
             Dg=Dg*(Dalphaj-x0)+q(i)
          end do
          alpha(j)=alpha(j)-Dg%f(0)/Dg%f(1)/(one-Dg%f(0)/Dg%f(1)*zaccel(j))*brake
          write(10+j,*) real(alpha(j)),aimag(alpha(j))
          ztmp=q(m)
          do i=m-1,0,-1
             ztmp=ztmp*(alpha(j)-x0)+q(i)
          end do
          if(abs(ztmp).gt.tmp) then
             tmp=abs(ztmp)
          end if
          if(abs(ztmp).lt.tiny) iflag(j)=1
       end do
       if(sum(iflag(:)).eq.m) exit
       k=k+1
    end do
    write(*,*) "# DKA is converged after", k, "iterations."

    do i=1,m
       write(*,*) alpha(i)
       write(9,*) real(alpha(i)), aimag(alpha(i))
    end do
  end subroutine find_roots_dka
  
end module module_dka
