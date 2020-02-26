!> 自動微分
module module_autodiff
  use module_globals
  implicit none

  integer,parameter :: n=ntaylor

  private
  public iCj, init_autodiff, uninit_autodiff, ad, d2t

  !-------------------------------------------------------------------------------------------------
  type iCj
     real(8),allocatable :: a(:)
  end type iCj
  type(iCj),pointer :: choose(:) !< choose(i)%a(j)=iCjを(integer overflowを防ぐため)real(8)で保持

  real(8) :: fact(0:n) !< fact(n)=n!

  integer :: maxfdb
  !> Faa Di Brunoの公式
  type FaaDiBruno
     ! \sum i x Mi =nを満たすリストと(i,Mi)に依って決まるFaa di Brunoの公式に現れる係数を保存
     integer :: n !< n階微分に寄与する
     integer :: s !< 非零のMiの数
     integer,allocatable :: part(:) !< 整数nの分割
     integer,allocatable :: m(:) !< 非零のMiのリスト
     integer :: sum_m !< sum(m(:))
     integer,allocatable :: i(:) !< m(:)に対応するiのリスト
     real(8),allocatable :: c !< [n!/(M1!(1!)^M1,...,Mn!(n!)^Mn]を(integer overflowを防ぐため)real(8)で保持
  end type FaaDiBruno
  type(FaaDiBruno),allocatable :: fdb(:)

  ! workspaces for integer partitions
  integer :: idx, knd
  integer,allocatable :: parts(:)
  !-------------------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------------------
  !> Dcomplex8 type
  public Dcomplex8
  type Dcomplex8
     complex(8) :: f(0:n) !f(i)=(d^i)f/d(x^i) (for i=0...n)
  end type Dcomplex8
  !-------------------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------------------  
  !> Equal assignment
  public assignment (=)
  interface assignment (=)
     module procedure EqualDR  !Dcomplex=real
     module procedure EqualDC  !Dcomplex=complex
     module procedure EqualDD !Dcomplex=Dcomplex
  end interface assignment (=)

  !> unary operator +
  public operator (+)
  interface operator (+)
     module procedure PlusD
  end interface operator (+)

  !> Addition operator
  interface operator (+)
     module procedure AddRD
     module procedure AddCD
     module procedure AddDR
     module procedure AddDC
     module procedure AddDD
  end interface operator (+)

  !> Unary operator -
  public operator (-)  
  interface operator (-)
     module procedure MinusD
  end interface operator (-)

  !> Subtraction operator
  interface operator (-)
     module procedure SubtractRD
     module procedure SubtractCD
     module procedure SubtractDR
     module procedure SubtractDC
     module procedure SubtractDD
  end interface operator (-)

  !> Multiplication operator
  public operator (*)
  interface operator (*)
     module procedure MultiplyRD
     module procedure MultiplyCD
     module procedure MultiplyDR
     module procedure MultiplyDC
     module procedure MultiplyDD
  end interface operator (*)

  !> Division operator
  public operator (/)
  interface operator (/)
     module procedure DivideDD
     module procedure DivideDR
     module procedure DivideRD
     module procedure DivideDC
     module procedure DivideCD
  end interface operator (/)

  !> Power operator
  public operator (**)
  interface operator (**)
     module procedure PowerDI
     module procedure PowerDR
     module procedure PowerDC
     module procedure PowerRD
     module procedure PowerCD     
  end interface operator (**)

  ! !-------------------
  ! ! Logical operators.
  ! !-------------------
  ! 使う気がしないので定義しない
  ! ! Equal operator.
  ! interface operator (.eq.)  ! or (==)
  !    module procedure eq_dd
  !    module procedure eq_dr
  !    module procedure eq_rd
  !    module procedure eq_di
  !    module procedure eq_id
  ! end interface operator (.eq.)
  ! ! Not equal operator.
  ! interface operator (.ne.)  ! or (/=)
  !    module procedure ne_dd
  !    module procedure ne_dr
  !    module procedure ne_rd
  !    module procedure ne_di
  !    module procedure ne_id
  ! end interface operator (.ne.)
  ! ! Less than operator.
  ! interface operator (.lt.)  ! or (<)
  !    module procedure lt_dd
  !    module procedure lt_dr
  !    module procedure lt_rd
  !    module procedure lt_di
  !    module procedure lt_id
  ! end interface operator (.lt.)
  ! ! Less than or equal operator.
  ! interface operator (.le.)  ! or (<=)
  !    module procedure le_dd
  !    module procedure le_dr
  !    module procedure le_rd
  !    module procedure le_di
  !    module procedure le_id
  ! end interface operator (.le.)
  ! ! Greater than operator.
  ! interface operator (.gt.)  ! or (>)
  !    module procedure gt_dd
  !    module procedure gt_dr
  !    module procedure gt_rd
  !    module procedure gt_di
  !    module procedure gt_id
  ! end interface operator (.gt.)
  ! ! Greater than or equal operator.
  ! interface operator (.ge.)  ! or (>=)
  !    module procedure ge_dd
  !    module procedure ge_dr
  !    module procedure ge_rd
  !    module procedure ge_di
  !    module procedure ge_id
  ! end interface operator (.ge.)

  !----------------
  ! Math functions.
  !----------------

  ! ! Absolute value function
  ! |z|は正則関数でないので、Dcomplex8に対するabsは定義しない  
  ! public abs
  ! interface abs
  !    module procedure absD
  ! end interface abs

  ! ! Integer function
  ! 複素関数に対しては定義しない
  ! interface int
  !    module procedure intDual
  ! end interface int

  ! ! Nearest integer function
  ! 複素関数に対しては定義しない  
  ! interface nint
  !    module procedure nintDual
  ! end interface nint

  !> Real function
  public real
  interface real
     module procedure realD
  end interface real

  !> Imag function
  public aimag
  interface aimag
     module procedure aimagD
  end interface aimag
 
  !> Congjugate function
  public conjg
  interface conjg
     module procedure conjgD
  end interface conjg

  ! ! Sign function
  ! 使う気がしないので定義しない
  ! interface sign
  !    module procedure sign_dd
  !    module procedure sign_dr
  !    module procedure sign_rd
  ! end interface sign

  !> Sine function
  public sin
  interface sin
     module procedure sinD
  end interface sin

  !> Cosine function
  public cos
  interface cos
     module procedure cosD
  end interface cos

  !> Tangent function
  public tan
  interface tan
     module procedure tanD
  end interface tan

  !> Sqrt function
  public sqrt
  interface sqrt
     module procedure sqrtD
  end interface sqrt

  !> Log function
  public log
  interface log
     module procedure logD
  end interface log

  !> Log10 function
  public log10
  interface log10
     module procedure log10C
     module procedure log10D
  end interface log10

  !> Exp function
  public exp
  interface exp
     module procedure expD
  end interface exp

  !> Sinh function
  public sinh
  interface sinh
     module procedure sinhD
  end interface sinh

  !> Cosh function
  public cosh
  interface cosh
     module procedure coshD
  end interface cosh

  !> Tanh function
  public tanh
  interface tanh
     module procedure tanhD
  end interface tanh

  !> Acos function
  public acos  
  interface acos
     module procedure acosD
  end interface acos

  !> Asin function
  public asin
  interface asin
     module procedure asinD
  end interface asin

  !> Atan function
  public atan
  interface atan
     module procedure atanD
  end interface atan

  ! ! Atan2 function
  ! atan2はそもそも複素数に対して定義されていないのでまぁいいか...。
  ! public atan2
  ! interface atan2
  !    module procedure atan2DD
  !    module procedure atan2DR
  !    module procedure atan2RD
  !    module procedure atan2DC
  !    module procedure atan2CD     
  ! end interface atan2

  !> Bessel (Jn) function
  public bessel_jn
  interface bessel_jn
     module procedure bessel_jnD
  end interface bessel_jn

  !> Bessel (J0) function  
  public bessel_j0
  interface bessel_j0
     module procedure bessel_j0D
  end interface bessel_j0

  !> Bessel (J1) function
  public bessel_j1
  interface bessel_j1
     module procedure bessel_j1D
  end interface bessel_j1

  !> Bessel (Yn) function
  public bessel_yn
  interface bessel_yn
     module procedure bessel_ynD
  end interface bessel_yn

  !> Bessel (Y0) function
  public bessel_y0
  interface bessel_y0
     module procedure bessel_y0D
  end interface bessel_y0

  !> Bessel (Y1) function
  public bessel_y1
  interface bessel_y1
     module procedure bessel_y1D
  end interface bessel_y1

  !> Hankel (H^(m)_n) function
  public hankel
  interface hankel
     module procedure hankelD
  end interface hankel
  
  ! ! Max function (limited to combinations below, but that
  ! ! can be extended)
  ! 使わない気がする
  ! interface max
  !    module procedure max_dd
  !    module procedure max_ddd
  !    module procedure max_dr
  !    module procedure max_rd
  ! end interface max

  ! ! Min function (limited for now to 2 arguments, but that
  ! ! can be extended)
  ! 使わない気がする  
  ! interface min
  !    module procedure min_dd
  !    module procedure min_dr
  !    module procedure min_rd
  ! end interface min

  !=================================================================

contains

  !> convrt f^(n) to f^(n)/n!
  type(Dcomplex8) function d2t(z) result(out)
    type(Dcomplex8),intent(in) :: z
    integer :: i
    out%f(0)=z%f(0)
    out%f(1)=z%f(1)
    do i=2,ntaylor
       out%f(i)=z%f(i)/fact(i)
    end do
  end function d2t
  
  !> assign complex(8) independent variable to type(Dcomplex8)
  type(Dcomplex8) function ad(z) result(out)
    complex(8),intent(in) :: z

    out%f(0)=z
    out%f(1)=one
    out%f(2:ntaylor)=zero

  end function ad
  
  !=================================================================================================
  !> compute binomial coefficient
  subroutine gen_choose(choose_)
    type(iCj),intent(inout),target :: choose_(0:n)
    integer :: i, j
    real(8) :: tmp(-1:n+1)
    tmp(:)=0.d0
    tmp(0)=1.d0
    tmp(1)=1.d0
    allocate(choose_(0)%a(0:0))
    choose_(0)%a(0)=1.d0
    do i=1,n
       allocate(choose_(i)%a(0:i))
       do j=0,i
          choose_(i)%a(j)=tmp(j)
       end do
       do j=n,0,-1
          tmp(j)=tmp(j)+tmp(j-1)
       end do
    end do
    ! do i=0,n
    !    write(*,*) i,int(choose(i)%a(:))
    ! end do
    choose=>choose_
  end subroutine gen_choose

  !-------------------------------------------------------------------------------------------------

  !> compute and store factorial
  subroutine gen_factorial
    integer :: i

    fact(0)=1.d0
    fact(1)=1.d0
    do i=2,n
       fact(i)=fact(i-1)*dble(i)
    end do

  end subroutine gen_factorial

  !-------------------------------------------------------------------------------------------------

  !> the number of partition of an integer
  recursive subroutine count_partitions(p, px)
    implicit none
    integer, intent(in) :: p, px
    integer :: k, ist

    if(p.eq.0) then
       knd=knd+1;
    else
       do k=px,1,-1
          ist=idx
          idx=idx+1
          call count_partitions(p-k,min(p-k,k))
          idx=ist
       end do
    end if

  end subroutine count_partitions

  !-------------------------------------------------------------------------------------------------

  !> partition of an integer  
  recursive subroutine gen_partition(p, px)
    implicit none
    integer, intent(in) :: p, px
    integer :: k, ist

    if(p==0) then
       knd=knd+1;
       allocate(fdb(knd)%part(idx))
       fdb(knd)%part(:)=parts(1:idx)
    else
       do k=px,1,-1
          ist=idx
          idx=idx+1
          parts(idx)=k
          call gen_partition(p-k,min(p-k,k))
          idx=ist
       end do
    end if

  end subroutine gen_partition

  !-------------------------------------------------------------------------------------------------

  !> compute coefficients for Faa di Bruno formula
  subroutine gen_fdb
    integer :: i, j, cnt

    knd=0
    do j=1,n
       idx=0
       call count_partitions(j,j)
    end do
    maxfdb=knd
    ! write(*,*) "ntaylor,maxfdb",ntaylor, maxfdb
    allocate(fdb(maxfdb))
    knd=0
    do j=1,n
       allocate(parts(j))
       idx=0
       parts=0
       call gen_partition(j,j)
       deallocate(parts)
    end do

    do j=1,maxfdb
       fdb(j)%n=sum(fdb(n)%part(:))
       cnt=1
       do i=2,size(fdb(j)%part)
          if(.not.(fdb(j)%part(i)==fdb(j)%part(i-1))) then
             cnt=cnt+1
          end if
       end do
       fdb(j)%s=cnt
       allocate(fdb(j)%i(cnt))
       allocate(fdb(j)%m(cnt))
    end do

    do j=1,maxfdb
       fdb(j)%n=sum(fdb(j)%part(:))
       fdb(j)%m(:)=0
       cnt=1
       fdb(j)%i(cnt)=fdb(j)%part(cnt)
       fdb(j)%m(cnt)=1
       do i=2,size(fdb(j)%part)
          if(.not.(fdb(j)%part(i)==fdb(j)%part(i-1))) then
             cnt=cnt+1
             fdb(j)%i(cnt)=fdb(j)%part(i)
          else
          end if
          fdb(j)%m(cnt)=fdb(j)%m(cnt)+1
       end do
       fdb(j)%sum_m=sum(fdb(j)%m)
       fdb(j)%c=fact(fdb(j)%n)
       !fdb(j)%c=1.d0
       do i=1,fdb(j)%s
          fdb(j)%c=fdb(j)%c/(fact(fdb(j)%m(i))*fact(fdb(j)%i(i))**fdb(j)%m(i))
       end do
    end do

    ! do j=1,maxfdb
    !    write(*,*) j,fdb(j)%n,fdb(j)%sum_m,fdb(j)%s,fdb(j)%m,fdb(j)%i,fdb(j)%c
    ! end do

  end subroutine gen_fdb

  !-------------------------------------------------------------------------------------------------  
  !> initialise
  subroutine init_autodiff(choose_)
    type(iCj),intent(inout) :: choose_(0:n)
    call gen_choose(choose_)
    call gen_factorial
    call gen_fdb
  end subroutine init_autodiff

  !-------------------------------------------------------------------------------------------------
  !> finalise
  subroutine uninit_autodiff(choose_)
    type(iCj),intent(inout) :: choose_(0:n)
    integer :: i
    do i=0,n
       deallocate(choose_(i)%a)
    end do
    deallocate(fdb)
  end subroutine uninit_autodiff

  !=================================================================================================  

  !> Functions for the equal assignment.
  elemental subroutine EqualDR(res,inp)
    type(Dcomplex8),intent(out) :: res
    real(8),intent(in) :: inp
    res%f(0)=cmplx(inp,0.d0,kind(1.d0))
    res%f(1:n)=0.d0
  end subroutine EqualDR

  !> Functions for the equal assignment.  
  elemental subroutine EqualDC(res,inp)
    type(Dcomplex8),intent(out) :: res
    complex(8),intent(in) :: inp
    res%f(0)=inp
    res%f(1:n)=0.d0
  end subroutine EqualDC

  !> Functions for the equal assignment.  
  elemental subroutine EqualDD(res,inp)
    type(Dcomplex8),intent(out) :: res
    type(Dcomplex8),intent(in) :: inp
    res%f(:)=inp%f(:)
  end subroutine EqualDD

  !> Function for the unary operator +.
  elemental function PlusD(v1) result(v2)
    type(Dcomplex8),intent(in) :: v1
    type(Dcomplex8):: v2
    v2%f(:)=v1%f(:)
  end function PlusD

  !> Functions for the addition operator.
  elemental function AddRD(v1,v2) result(v3)
    real(8),intent(in) :: v1
    type(Dcomplex8),intent(in) :: v2
    type(Dcomplex8) :: v3
    v3%f(0)=cmplx(v1+real(v2%f(0)),aimag(v2%f(0)),kind(1.d0))
    v3%f(1:n)=v2%f(1:n)
  end function AddRD

  !> Functions for the addition operator.  
  elemental function AddDR(v1,v2) result(v3)
    type(Dcomplex8),intent(in) :: v1
    real(8),intent(in) :: v2
    type(Dcomplex8) :: v3
    v3%f(0)=cmplx(real(v1%f(0))+v2,aimag(v1%f(0)),kind(1.d0))
    v3%f(1:n)=v1%f(1:n)
  end function AddDR

  !> Functions for the addition operator.  
  elemental function AddCD(v1,v2) result(v3)
    complex(8),intent(in) :: v1
    type(Dcomplex8),intent(in) :: v2
    type(Dcomplex8) :: v3
    v3%f(0)=v1+v2%f(0)
    v3%f(1:n)=v2%f(1:n)
  end function AddCD

  !> Functions for the addition operator.
  elemental function AddDC(v1,v2) result(v3)
    type(Dcomplex8),intent(in) :: v1
    complex(8),intent(in) :: v2
    type(Dcomplex8) :: v3
    v3%f(0)=v1%f(0)+v2
    v3%f(1:n)=v1%f(1:n)
  end function AddDC
  
  !> Functions for the addition operator.
  elemental function AddDD(v1,v2) result(v3)
    type(Dcomplex8),intent(in) :: v1, v2
    type(Dcomplex8) :: v3
    v3%f(:)=v1%f(:)+v2%f(:)
  end function AddDD

  !> Function for the unary operator -.
  elemental function MinusD(v1) result(v2)
    type(Dcomplex8),intent(in) :: v1
    type(Dcomplex8):: v2
    v2%f(:)=-v1%f(:)
  end function MinusD

  !> Functions for the subtraction operator.
  elemental function SubtractRD(v1,v2) result(v3)
    real(8),intent(in) :: v1
    type(Dcomplex8),intent(in) :: v2
    type(Dcomplex8) :: v3
    v3%f(0)=cmplx(v1-real(v2%f(0)),-aimag(v2%f(0)),kind(1.d0))
    v3%f(1:n)=-v2%f(1:n)
  end function SubtractRD

  !> Functions for the subtraction operator.  
  elemental function SubtractDR(v1,v2) result(v3)
    type(Dcomplex8),intent(in) :: v1
    real(8),intent(in) :: v2
    type(Dcomplex8) :: v3
    v3%f(0)=cmplx(real(v1%f(0))-v2,aimag(v1%f(0)),kind(1.d0))
    v3%f(1:n)=v1%f(1:n)
  end function SubtractDR

  !> Functions for the subtraction operator.  
  elemental function SubtractCD(v1,v2) result(v3)
    complex(8),intent(in) :: v1
    type(Dcomplex8),intent(in) :: v2
    type(Dcomplex8) :: v3
    v3%f(0)=v1-v2%f(0)
    v3%f(1:n)=-v2%f(1:n)
  end function SubtractCD

  !> Functions for the subtraction operator.  
  elemental function SubtractDC(v1,v2) result(v3)
    type(Dcomplex8),intent(in) :: v1
    complex(8),intent(in) :: v2
    type(Dcomplex8) :: v3
    v3%f(0)=v1%f(0)-v2
    v3%f(1:n)=v1%f(1:n)
  end function SubtractDC

  !> Functions for the subtraction operator.  
  elemental function SubtractDD(v1,v2) result(v3)
    type(Dcomplex8),intent(in) :: v1, v2
    type(Dcomplex8) :: v3
    v3%f(:)=v1%f(:)-v2%f(:)
  end function SubtractDD
  
  !> Functions for the multiplication operator.
  elemental function MultiplyRD(v1,v2) result(v3)
    real(8),intent(in) :: v1
    type(Dcomplex8),intent(in) :: v2
    type(Dcomplex8) :: v3
    v3%f(:)=v1*v2%f(:)
  end function MultiplyRD

  !> Functions for the multiplication operator.  
  elemental function MultiplyDR(v1,v2) result(v3)
    type(Dcomplex8),intent(in) :: v1
    real(8),intent(in) :: v2
    type(Dcomplex8) :: v3
    v3%f(:)=v1%f(:)*v2
  end function MultiplyDR

  !> Functions for the multiplication operator.  
  elemental function MultiplyCD(v1,v2) result(v3)
    complex(8),intent(in) :: v1
    type(Dcomplex8),intent(in) :: v2
    type(Dcomplex8) :: v3
    v3%f(:)=v1*v2%f(:)
  end function MultiplyCD

  !> Functions for the multiplication operator.
  elemental function MultiplyDC(v1,v2) result(v3)
    type(Dcomplex8),intent(in) :: v1
    complex(8),intent(in) :: v2
    type(Dcomplex8) :: v3
    v3%f(:)=v1%f(:)*v2
  end function MultiplyDC

  !> Functions for the multiplication operator.  
  elemental function MultiplyDD(v1,v2) result(v3)
    type(Dcomplex8),intent(in) :: v1
    type(Dcomplex8),intent(in) :: v2
    type(Dcomplex8) :: v3
    integer :: i, j
    v3%f(0)=v1%f(0)*v2%f(0)
    do i=1,n
       v3%f(i)=zero
       do j=0,i
          v3%f(i)=v3%f(i)+choose(i)%a(j)*v1%f(j)*v2%f(i-j)
       end do
    end do
  end function MultiplyDD

  !> Function for the inversion operator.
  elemental function invertD(inp) result(res)
    type(Dcomplex8),intent(in) :: inp
    type(Dcomplex8) :: res
    integer :: i
    real(8) :: pm
    complex(8) :: ztmp, fext(0:n)
    pm=1.d0
    ztmp=one/inp%f(0)
    do i=0,n
       fext(i)=pm*fact(i)*ztmp
       ztmp=ztmp/inp%f(0)
       pm=-pm
    end do
    res=polyfdb(fext,inp)
  end function invertD

  !> Functions for the division operator.
  elemental function divideDD(v1,v2) result(v3)
    type(Dcomplex8), intent(in) :: v1, v2
    type(Dcomplex8) :: v3
    v3=v1*invertD(v2)
  end function divideDD

  !> Functions for the division operator.  
  elemental function divideDR(v1,v2) result(v3)
    type(Dcomplex8),intent(in) :: v1
    real(8),intent(in) :: v2
    type(Dcomplex8) :: v3
    v3%f(:)=v1%f(:)/v2
  end function divideDR

  !> Functions for the division operator.  
  elemental function divideRD(v1,v2) result(v3)
    real(8),intent(in) :: v1
    type(Dcomplex8),intent(in) :: v2
    type(Dcomplex8) :: v3
    v3=v1*invertD(v2)
  end function divideRD

  !> Functions for the division operator.  
  elemental function divideDC(v1,v2) result(v3)
    type(Dcomplex8),intent(in) :: v1
    complex(8),intent(in) :: v2
    type(Dcomplex8) :: v3
    v3%f(:)=v1%f(:)/v2
  end function divideDC

  !> Functions for the division operator.  
  elemental function divideCD(v1,v2) result(v3)
    complex(8),intent(in) :: v1
    type(Dcomplex8),intent(in) :: v2
    type(Dcomplex8) :: v3
    v3=v1*invertD(v2)
  end function divideCD
  
  !> Functions for the power operator.
  elemental function PowerDI(v1,v2) result(v3)
    type(Dcomplex8),intent(in) :: v1
    integer,intent(in) :: v2
    type(Dcomplex8) :: v3
    integer :: i
    real(8) :: tmp
    complex(8) :: ztmp
    complex(8) :: fext(0:n)
    tmp=1.d0
    ztmp=v1%f(0)**v2
    do i=0,n
       fext(i)=tmp*ztmp
       tmp=tmp*dble(v2-i)
       ztmp=ztmp/v1%f(0)
    end do 
    v3=polyfdb(fext,v1)
  end function PowerDI

  !> Functions for the power operator.  
  elemental function PowerDR(v1,v2) result(v3)
    type(Dcomplex8),intent(in) :: v1
    real,intent(in) :: v2
    type(Dcomplex8) :: v3
    integer :: i
    real(8) :: tmp
    complex(8) :: ztmp
    complex(8) :: fext(0:n)
    tmp=1.d0
    ztmp=v1%f(0)**v2
    do i=0,n
       fext(i)=tmp*ztmp
       tmp=tmp*(v2-dble(i))
       ztmp=ztmp/v1%f(0)
    end do 
    v3=polyfdb(fext,v1)
  end function PowerDR

  !> Functions for the power operator.  
  elemental function PowerDC(v1,v2) result(v3)
    type(Dcomplex8),intent(in) :: v1
    complex(8),intent(in) :: v2
    type(Dcomplex8) :: v3
    integer :: i
    complex(8) :: tmp, ztmp
    complex(8) :: fext(0:n)
    tmp=one
    ztmp=v1%f(0)**v2
    do i=0,n
       fext(i)=tmp*ztmp
       tmp=tmp*(v2-cmplx(dble(i),0.d0,kind(1.d0)))
       ztmp=ztmp/v1%f(0)
    end do 
    v3=polyfdb(fext,v1)
  end function PowerDC

  !> Functions for the power operator.  
  elemental function PowerRD(v1,v2) result(v3)
    real(8),intent(in) :: v1
    type(Dcomplex8),intent(in) :: v2
    type(Dcomplex8) :: v3
    integer :: i
    complex(8) :: fext(0:n)
    fext(0)=v1**v2%f(0)
    do i=1,n
       fext(i)=fext(i-1)*log(v1)
    end do
    v3=polyfdb(fext,v2)
  end function PowerRD

  !> Functions for the power operator.  
  elemental function PowerCD(v1,v2) result(v3)
    complex(8),intent(in) :: v1
    type(Dcomplex8),intent(in) :: v2
    type(Dcomplex8) :: v3
    integer :: i
    complex(8) :: fext(0:n)
    fext(0)=v1**v2%f(0)
    do i=1,n
       fext(i)=fext(i-1)*log(v1)
    end do
    v3=polyfdb(fext,v2)
  end function PowerCD
  
  ! !----------------------------------
  ! ! Functions for the equal operator.
  ! !----------------------------------

  ! logical function eq_dd(lhs, rhs)
  !   type (DualNumber), intent(in) :: lhs, rhs
  !   eq_dd = lhs%f0 == rhs%f0
  ! end function eq_dd

  ! logical function eq_dr(lhs, rhs)
  !   type (DualNumber), intent(in) :: lhs
  !   real, intent(in) :: rhs
  !   eq_dr = lhs%f0 == rhs
  ! end function eq_dr

  ! logical function eq_rd(lhs, rhs)
  !   real, intent(in) :: lhs
  !   type (DualNumber), intent(in) :: rhs
  !   eq_rd = lhs == rhs%f0
  ! end function eq_rd

  ! logical function eq_di(lhs, rhs)
  !   type (DualNumber), intent(in) :: lhs
  !   integer, intent(in) :: rhs
  !   eq_di = lhs%f0 == rhs
  ! end function eq_di

  ! logical function eq_id(lhs, rhs)
  !   integer, intent(in) :: lhs
  !   type (DualNumber), intent(in) :: rhs
  !   eq_id = lhs == rhs%f0
  ! end function eq_id

  ! !--------------------------------------
  ! ! Functions for the not equal operator.
  ! !--------------------------------------

  ! logical function ne_dd(lhs, rhs)
  !   type (DualNumber), intent(in) :: lhs, rhs
  !   ne_dd = lhs%f0 /= rhs%f0
  ! end function ne_dd

  ! logical function ne_dr(lhs, rhs)
  !   type (DualNumber), intent(in) :: lhs
  !   real, intent(in) :: rhs
  !   ne_dr = lhs%f0 /= rhs
  ! end function ne_dr

  ! logical function ne_rd(lhs, rhs)
  !   real, intent(in) :: lhs
  !   type (DualNumber), intent(in) :: rhs
  !   ne_rd = lhs /= rhs%f0
  ! end function ne_rd

  ! logical function ne_di(lhs, rhs)
  !   type (DualNumber), intent(in) :: lhs
  !   integer, intent(in) :: rhs
  !   ne_di = lhs%f0 /= rhs
  ! end function ne_di

  ! logical function ne_id(lhs, rhs)
  !   integer, intent(in) :: lhs
  !   type (DualNumber), intent(in) :: rhs
  !   ne_id = lhs /= rhs%f0
  ! end function ne_id

  ! !--------------------------------------
  ! ! Functions for the less than operator.
  ! !--------------------------------------

  ! logical function lt_dd(lhs, rhs)
  !   type (DualNumber), intent(in) :: lhs, rhs
  !   lt_dd = lhs%f0 < rhs%f0
  ! end function lt_dd

  ! logical function lt_dr(lhs, rhs)
  !   type (DualNumber), intent(in) :: lhs
  !   real, intent(in) :: rhs
  !   lt_dr = lhs%f0 < rhs
  ! end function lt_dr

  ! logical function lt_rd(lhs, rhs)
  !   real, intent(in) :: lhs
  !   type (DualNumber), intent(in) :: rhs
  !   lt_rd = lhs < rhs%f0
  ! end function lt_rd

  ! logical function lt_di(lhs, rhs)
  !   type (DualNumber), intent(in) :: lhs
  !   integer, intent(in) :: rhs
  !   lt_di = lhs%f0 < rhs
  ! end function lt_di

  ! logical function lt_id(lhs, rhs)
  !   integer, intent(in) :: lhs
  !   type (DualNumber), intent(in) :: rhs
  !   lt_id = lhs < rhs%f0
  ! end function lt_id

  ! !-----------------------------------------------
  ! ! Functions for the less than or equal operator.
  ! !-----------------------------------------------

  ! logical function le_dd(lhs, rhs)
  !   type (DualNumber), intent(in) :: lhs, rhs
  !   le_dd = lhs%f0 <= rhs%f0
  ! end function le_dd

  ! logical function le_dr(lhs, rhs)
  !   type (DualNumber), intent(in) :: lhs
  !   real, intent(in) :: rhs
  !   le_dr = lhs%f0 <= rhs
  ! end function le_dr

  ! logical function le_rd(lhs, rhs)
  !   real, intent(in) :: lhs
  !   type (DualNumber), intent(in) :: rhs
  !   le_rd = lhs <= rhs%f0
  ! end function le_rd

  ! logical function le_di(lhs, rhs)
  !   type (DualNumber), intent(in) :: lhs
  !   integer, intent(in) :: rhs
  !   le_di = lhs%f0 <= rhs
  ! end function le_di

  ! logical function le_id(lhs, rhs)
  !   integer, intent(in) :: lhs
  !   type (DualNumber), intent(in) :: rhs
  !   le_id = lhs <= rhs%f0
  ! end function le_id

  ! !-----------------------------------------
  ! ! Functions for the greater than operator.
  ! !-----------------------------------------

  ! logical function gt_dd(lhs, rhs)
  !   type (DualNumber), intent(in) :: lhs, rhs
  !   gt_dd = lhs%f0 > rhs%f0
  ! end function gt_dd

  ! logical function gt_dr(lhs, rhs)
  !   type (DualNumber), intent(in) :: lhs
  !   real, intent(in) :: rhs
  !   gt_dr = lhs%f0 > rhs
  ! end function gt_dr

  ! logical function gt_rd(lhs, rhs)
  !   real, intent(in) :: lhs
  !   type (DualNumber), intent(in) :: rhs
  !   gt_rd = lhs > rhs%f0
  ! end function gt_rd

  ! logical function gt_di(lhs, rhs)
  !   type (DualNumber), intent(in) :: lhs
  !   integer, intent(in) :: rhs
  !   gt_di = lhs%f0 > rhs
  ! end function gt_di

  ! logical function gt_id(lhs, rhs)
  !   integer, intent(in) :: lhs
  !   type (DualNumber), intent(in) :: rhs
  !   gt_id = lhs > rhs%f0
  ! end function gt_id

  ! !--------------------------------------------------
  ! ! Functions for the greater than or equal operator.
  ! !--------------------------------------------------

  ! logical function ge_dd(lhs, rhs)
  !   type (DualNumber), intent(in) :: lhs, rhs
  !   ge_dd = lhs%f0 >= rhs%f0
  ! end function ge_dd

  ! logical function ge_dr(lhs, rhs)
  !   type (DualNumber), intent(in) :: lhs
  !   real, intent(in) :: rhs
  !   ge_dr = lhs%f0 >= rhs
  ! end function ge_dr

  ! logical function ge_rd(lhs, rhs)
  !   real, intent(in) :: lhs
  !   type (DualNumber), intent(in) :: rhs
  !   ge_rd = lhs >= rhs%f0
  ! end function ge_rd

  ! logical function ge_di(lhs, rhs)
  !   type (DualNumber), intent(in) :: lhs
  !   integer, intent(in) :: rhs
  !   ge_di = lhs%f0 >= rhs
  ! end function ge_di

  ! logical function ge_id(lhs, rhs)
  !   integer, intent(in) :: lhs
  !   type (DualNumber), intent(in) :: rhs
  !   ge_id = lhs >= rhs%f0
  ! end function ge_id

  !----------------
  ! Math functions.
  !----------------
  !> Faa di Bruno
  pure function polyfdb(f,g) result(v)
    complex(8),intent(in) :: f(0:n)
    type(Dcomplex8), intent(in) :: g
    type(Dcomplex8) :: v
    integer :: j
    v%f(:)=zero
    v%f(0)=f(0)
    do j=1,maxfdb
       v%f(fdb(j)%n)=v%f(fdb(j)%n)+f(fdb(j)%sum_m)*fdb(j)%c*product(g%f(fdb(j)%i(:))**fdb(j)%m)
    end do
  end function polyfdb

  ! ! Absolute value function.
  ! |z|は正則関数でないので、Dcomplex8に対するabsは定義しない.
  ! elemental function absD(v1) result(v2)
  !   type(Dcomplex8),intent(in) :: v1
  !   type(Dcomplex8) :: v2
  ! end function absD

  ! ! Integer function.
  ! 複素関数に対しては定義しない
  ! function intDual(v1) result (v2)
  !   type (DualNumber), intent(in) :: v1
  !   integer :: v2
  !   v2 = int(v1%f0)
  ! end function intDual

  ! ! Nearest integer function.
  ! 複素関数に対しては定義しない  
  ! function nintDual(v1) result (v2)
  !   type (DualNumber), intent(in) :: v1
  !   integer :: v2
  !   v2 = nint(v1%f0)
  ! end function nintDual

  !> Real function.
  elemental function realD(v1) result (v2)
    type(Dcomplex8), intent(in) :: v1
    type(Dcomplex8) :: v2
    v2%f(:)=cmplx(real(v1%f(:)),0.d0,kind(1.d0))
  end function realD

  !> Imag function.
  elemental function aimagD(v1) result (v2)
    type(Dcomplex8), intent(in) :: v1
    type(Dcomplex8) :: v2
    v2%f(:)=cmplx(aimag(v1%f(:)),0.d0,kind(1.d0))
  end function aimagD
  
  !> Conjugate function.
  elemental function conjgD(v1) result(v2)
    type(Dcomplex8),intent(in) :: v1
    type(Dcomplex8) :: v2
    v2%f(:)=conjg(v1%f(:))
  end function conjgD

  ! ! Functions for the sign function.
  ! 使わない気がする
  ! function sign_dd(v1,v2) result (v3)
  !   type (DualNumber), intent(in) :: v1, v2
  !   type (DualNumber) :: v3
  !   real :: ssign
  !   if(v2%f0 < 0.0) then
  !      ssign = -1.0
  !   else
  !      ssign =  1.0
  !   endif
  !   v3 = ssign*v1
  ! end function sign_dd
  ! function sign_dr(v1,v2) result (v3)
  !   type (DualNumber), intent(in) :: v1
  !   real,              intent(in) :: v2
  !   type (DualNumber) :: v3
  !   real :: ssign
  !   if(v2 < 0.0) then
  !      ssign = -1.0
  !   else
  !      ssign =  1.0
  !   endif
  !   v3 = ssign*v1
  ! end function sign_dr
  ! function sign_rd(v1,v2) result (v3)
  !   real,              intent(in) :: v1
  !   type (DualNumber), intent(in) :: v2
  !   type (DualNumber) :: v3
  !   real :: ssign
  !   if(v2%f0 < 0.0) then
  !      ssign = -1.0
  !   else
  !      ssign =  1.0
  !   endif
  !   v3 = ssign*v1
  ! end function sign_rd

  !> Sine function.
  elemental function sinD(v1) result(v2)
    TYPE(Dcomplex8),intent(in) :: v1
    TYPE(Dcomplex8) :: v2
    integer :: i
    complex(8) :: fext(0:n)
    complex(8) :: tmp(0:3)
    tmp(0)=sin(v1%f(0))
    tmp(1)=cos(v1%f(0))
    tmp(2)=-tmp(0)
    tmp(3)=-tmp(1)
    do i=0,n
       fext(i)=tmp(mod(i,4))
    end do
    v2=polyfdb(fext,v1)
  end function sinD

  !> Cosine function.
  elemental function cosD(v1) result(v2)
    TYPE(Dcomplex8),intent(in) :: v1
    TYPE(Dcomplex8) :: v2
    integer :: i
    complex(8) :: fext(0:n)
    complex(8) :: tmp(0:3)
    tmp(0)=cos(v1%f(0))
    tmp(1)=-sin(v1%f(0))
    tmp(2)=-tmp(0)
    tmp(3)=-tmp(1)
    do i=0,n
       fext(i)=tmp(mod(i,4))
    end do
    v2=polyfdb(fext,v1)
  end function cosD
    
  !> Tangent function.
  elemental function tanD(v1) result (v2)
    type(Dcomplex8), intent(in) :: v1
    type(Dcomplex8) :: v2
    v2=sin(v1)/cos(v1)
  end function tanD

  !> Sqrt function
  elemental function sqrtD(v1) result(v2)
    type(Dcomplex8),intent(in) :: v1
    type(Dcomplex8) :: v2
    integer :: i
    complex(8) :: fext(0:n)
    fext(0)=sqrt(v1%f(0))
    fext(1)=0.5d0/sqrt(v1%f(0))
    do i=1,n-1
       fext(i+1)=-fext(i)*(dble(i)-0.5d0)/v1%f(0)
    end do
    v2=polyfdb(fext,v1)
  end function sqrtD

  !> Log function
  elemental function logD(v1) result(v2)
    type(Dcomplex8),intent(in) :: v1
    type(Dcomplex8) :: v2
    integer :: i
    real(8) :: pm
    complex(8) :: ztmp, fext(0:n)
    fext(0)=log(v1%f(0)) !returns the principal value
    pm=1.d0
    ztmp=one/v1%f(0)
    do i=1,n
       fext(i)=pm*fact(i-1)*ztmp
       ztmp=ztmp/v1%f(0)
       pm=-pm
    end do
    v2=polyfdb(fext,v1)
  end function logD

  !> Log10 for complex variable
  elemental function log10C(v1) result(v2)
    complex(8),intent(in) :: v1
    complex(8) :: v2
    v2=log(v1)/log(10.d0) !returns the principal value
  end function log10C
  
  !> Log function
  elemental function log10D(v1) result(v2)
    type(Dcomplex8),intent(in) :: v1
    type(Dcomplex8) :: v2
    integer :: i
    real(8) :: pm, inv_log10
    complex(8) :: ztmp, fext(0:n)
    inv_log10=1.d0/log(10.d0)
    fext(0)=log(v1%f(0))*inv_log10 !returns the principal value
    pm=1.d0
    ztmp=one/v1%f(0)
    do i=1,n
       fext(i)=pm*fact(i-1)*ztmp*inv_log10
       ztmp=ztmp/v1%f(0)
       pm=-pm
    end do
    v2=polyfdb(fext,v1)
  end function log10D

  !> Exp function
  elemental function expD(v1) result(v2)
    type(Dcomplex8),intent(in) :: v1
    type(Dcomplex8) :: v2
    complex(8) :: fext(0:n)
    fext(:)=exp(v1%f(0))
    v2=polyfdb(fext,v1)
  end function expD

  !> Sinh function
  elemental function sinhD(v1) result(v2)
    type(Dcomplex8),intent(in) :: v1
    type(Dcomplex8) :: v2
    integer :: i
    real(8) :: pm
    complex(8) :: fext(0:n)
    pm=1.d0
    do i=0,n
       fext(i)=0.5d0*(sinh(v1%f(0))*(1.d0+pm)+cosh(v1%f(0))*(1.d0-pm))
       pm=-pm
    end do
    v2=polyfdb(fext,v1)
  end function sinhD

  !> Cosh function
  elemental function coshD(v1) result(v2)
    type(Dcomplex8),intent(in) :: v1
    type(Dcomplex8) :: v2
    integer :: i
    real(8) :: pm
    complex(8) :: fext(0:n)
    pm=-1.d0
    do i=0,n
       fext(i)=0.5d0*(sinh(v1%f(0))*(1.d0+pm)+cosh(v1%f(0))*(1.d0-pm))
       pm=-pm
    end do
    v2=polyfdb(fext,v1)
  end function coshD

  !> Tanh function
  elemental function tanhD(v1) result(v2)
    type(Dcomplex8),intent(in) :: v1
    type(Dcomplex8) :: v2
    v2=sinh(v1)/cosh(v1)
  end function tanhD

  !> Acos function
  !> branch cutは実軸から[-1,1]を除いた部分  
  elemental function acosD(v1) result(v2)
    type(Dcomplex8),intent(in) :: v1
    type(Dcomplex8) :: v2
    type(Dcomplex8) :: tmp
    integer :: i
    complex(8) :: fext(0:n)
    tmp%f(:)=zero
    tmp%f(0)=v1%f(0)
    tmp%f(1)=one
    tmp=-1.d0/sqrt(1.d0-tmp*tmp)
    fext(0)=acos(v1%f(0))
    do i=1,n
       fext(i)=tmp%f(i-1)
    end do
    v2=polyfdb(fext,v1)
  end function acosD

  !> Asin function
  !> branch cutは実軸から[-1,1]を除いた部分
  elemental function asinD(v1) result(v2)
    type(Dcomplex8),intent(in) :: v1
    type(Dcomplex8) :: v2
    type(Dcomplex8) :: tmp
    integer :: i
    complex(8) :: fext(0:n)
    tmp%f(:)=zero
    tmp%f(0)=v1%f(0)
    tmp%f(1)=one
    tmp=1.d0/sqrt(1.d0-tmp*tmp)
    fext(0)=asin(v1%f(0))
    do i=1,n
       fext(i)=tmp%f(i-1)
    end do
    v2=polyfdb(fext,v1)
  end function asinD

  !> Atan function
  !> branch cutは虚軸から[-i,i]を除いた部分  
  elemental function atanD(v1) result(v2)
    type(Dcomplex8),intent(in) :: v1
    type(Dcomplex8) :: v2
    type(Dcomplex8) :: tmp    
    integer :: i
    complex(8) :: fext(0:n)
    tmp%f(:)=zero
    tmp%f(0)=v1%f(0)
    tmp%f(1)=one
    tmp=1.d0/(1.d0+tmp*tmp)
    fext(0)=atan(v1%f(0))
    do i=1,n
       fext(i)=tmp%f(i-1)
    end do
    v2=polyfdb(fext,v1)
  end function atanD

  ! ! Atan2 function
  ! 使わない気がする
  ! elemental function atan2DD(v1,v2) result(v3)
  !   type(Dcomplex8),intent(in) :: v1, v2
  !   type(Dcomplex8) :: v3

  !   a = v1%f0; b = v1%f1
  !   c = v2%f0; d = v2%f1

  !   v3%f0 = atan2(a,c)
  !   v3%f1 = (c*b - a*d)/(a*a + c*c)
  ! end function atan2DD

  !> Bessel functions
  !> They are no longer elemental because zbesj is not pure.
  function bessel_jnD(v1,v2) result(v3)
    integer,intent(in) :: v1
    type(Dcomplex8),intent(in) :: v2
    type(Dcomplex8) :: v3
    integer :: i, k
    real(8) :: pm, coef
    complex(8) :: fext(0:n)
    ! for slatec
    integer :: nz, ierr, kode, nn
    real(8) :: zi, zr, fnu, cyr(max(v1-n,0):max(v1+n,n-v1)), cyi(max(v1-n,0):max(v1+n,n-v1))
    complex(8) :: besj(v1-n:v1+n)
    nn=max(v1+n,n-v1)-max(v1-n,0)+1
    kode=1
    zr=real(v2%f(0))
    zi=aimag(v2%f(0))
    fnu=max(v1-n,0)
    call zbesj(zr,zi,fnu,kode,nn,cyr(max(v1-n,0)),cyi(max(v1-n,0)),nz,ierr)
    do i=v1-n,min(v1+n,0)
       besj(i)=sign(1.d0,dble(i))**(abs(i))*cmplx(cyr(abs(i)),cyi(abs(i)),kind(1.d0))
    end do
    do i=max(0,v1-n),v1+n
       besj(i)=cmplx(cyr(i),cyi(i),kind(1.d0))
    end do
    fext(:)=zero
    fext(0)=besj(v1)
    coef=0.5d0
    do i=1,n
       pm=1.d0
       do k=0,i
          fext(i)=fext(i)+pm*choose(i)%a(k)*besj(v1-i+2*k)
          pm=-pm
       end do
       fext(i)=fext(i)*coef
       coef=coef*0.5d0
    end do
    v3=polyfdb(fext,v2)
  end function bessel_jnD

  function bessel_j0D(v1) result(v2)
    type(Dcomplex8),intent(in) :: v1
    type(Dcomplex8) :: v2
    v2=bessel_jn(0,v1)
  end function bessel_j0D

  function bessel_j1D(v1) result(v2)
    type(Dcomplex8),intent(in) :: v1
    type(Dcomplex8) :: v2
    v2=bessel_jn(1,v1)
  end function bessel_j1D

  !> Bessel functions
  !> They are no longer elemental because zbesj is not pure.
  function bessel_ynD(v1,v2) result(v3)
    integer,intent(in) :: v1
    type(Dcomplex8),intent(in) :: v2
    type(Dcomplex8) :: v3
    integer :: i, k
    real(8) :: pm, coef
    complex(8) :: fext(0:n)
    ! for slatec
    integer :: nz, ierr, kode, nn
    real(8) :: zi, zr, fnu, cyr(max(v1-n,0):max(v1+n,n-v1)), cyi(max(v1-n,0):max(v1+n,n-v1))
    real(8) :: wkr(max(v1-n,0):max(v1+n,n-v1)), wki(max(v1-n,0):max(v1+n,n-v1))
    complex(8) :: besy(v1-n:v1+n)
    nn=max(v1+n,n-v1)-max(v1-n,0)+1
    kode=1
    zr=real(v2%f(0))
    zi=aimag(v2%f(0))
    fnu=max(v1-n,0)
    call zbesy(zr,zi,fnu,kode,nn,cyr(max(v1-n,0)),cyi(max(v1-n,0)),nz,wkr(max(v1-n,0)),wki(max(v1-n,0)),ierr)
    do i=v1-n,min(v1+n,0)
       besy(i)=sign(1.d0,dble(i))**(abs(i))*cmplx(cyr(abs(i)),cyi(abs(i)),kind(1.d0))
    end do
    do i=max(v1-n,0),v1+n
       besy(i)=cmplx(cyr(i),cyi(i),kind(1.d0))
    end do
    fext(:)=zero
    fext(0)=besy(v1)
    coef=0.5d0
    do i=1,n
       pm=1.d0
       do k=0,i
          fext(i)=fext(i)+pm*choose(i)%a(k)*besy(v1-i+2*k)
          pm=-pm
       end do
       fext(i)=fext(i)*coef
       coef=coef*0.5d0
    end do
    v3=polyfdb(fext,v2)
  end function bessel_ynD

  !> Bessel function
  function bessel_y0D(v1) result(v2)
    type(Dcomplex8),intent(in) :: v1
    type(Dcomplex8) :: v2
    v2=bessel_yn(0,v1)
  end function bessel_y0D

  !> Bessel function  
  function bessel_y1D(v1) result(v2)
    type(Dcomplex8),intent(in) :: v1
    type(Dcomplex8) :: v2
    v2=bessel_yn(1,v1)
  end function bessel_y1D

  !> Hankel function  
  function hankelD(m,v1,v2) result(v3)
    integer,intent(in) :: m, v1
    type(Dcomplex8),intent(in) :: v2
    type(Dcomplex8) :: v3
    integer :: i, k
    real(8) :: pm, coef
    complex(8) :: fext(0:n)
    ! for slatec
    integer :: nz, ierr, kode, nn
    real(8) :: zi, zr, fnu
    real(8) :: cyr(max(v1-n,0):max(v1+n,n-v1)), cyi(max(v1-n,0):max(v1+n,n-v1))
    real(8) :: cyr2(max(v1-n,0):max(v1+n,n-v1)), cyi2(max(v1-n,0):max(v1+n,n-v1))
    real(8) :: wkr(max(v1-n,0):max(v1+n,n-v1)), wki(max(v1-n,0):max(v1+n,n-v1))
    complex(8) :: besh(v1-n:v1+n)
    real(8) :: pm2
    pm2=sign(1.d0,1.5d0-m)
    
    nn=max(v1+n,n-v1)-max(v1-n,0)+1
    kode=1
    zr=real(v2%f(0))
    zi=aimag(v2%f(0))
    fnu=max(v1-n,0)
    call zbesj(zr,zi,fnu,kode,nn,cyr(max(v1-n,0)),cyi(max(v1-n,0)),nz,ierr)
    call zbesy(zr,zi,fnu,kode,nn,cyr2(max(v1-n,0)),cyi2(max(v1-n,0)),nz,wkr(max(v1-n,0)),wki(max(v1-n,0)),ierr)
    do i=v1-n,min(v1+n,0)
       besh(i)=sign(1.d0,dble(i))**(abs(i))&
            *cmplx(cyr(abs(i))-pm2*cyi2(abs(i)),cyi(abs(i))+pm2*cyr2(abs(i)),kind(1.d0))
    end do
    do i=max(v1-n,0),v1+n
       besh(i)=cmplx(cyr(i)-pm2*cyi2(i),cyi(i)+pm2*cyr2(i),kind(1.d0))
    end do
    fext(:)=zero
    fext(0)=besh(v1)
    coef=0.5d0
    do i=1,n
       pm=1.d0
       do k=0,i
          fext(i)=fext(i)+pm*choose(i)%a(k)*besh(v1-i+2*k)
          pm=-pm
       end do
       fext(i)=fext(i)*coef
       coef=coef*0.5d0
    end do
    v3=polyfdb(fext,v2)

    ! v3=bessel_jn(v1,v2)+pm*ione*bessel_yn(v1,v2)

  end function hankelD
  
  ! 以下、使わない気がする
  ! ! Max functions
  ! function max_dd(v1, v2) result (v3)
  !   type (DualNumber), intent(in) :: v1, v2
  !   type (DualNumber) :: v3

  !   if(v1%f0 > v2%f0) then
  !      v3 = v1
  !   else
  !      v3 = v2
  !   endif
  ! end function max_dd

  ! function max_ddd(v1, v2, v3) result (v4)
  !   type (DualNumber), intent(in) :: v1, v2, v3
  !   type (DualNumber) :: v4

  !   if(v1%f0 > v2%f0) then
  !      v4 = v1
  !   else
  !      v4 = v2
  !   endif

  !   if(v3%f0 > v4%f0) v4 = v3
  ! end function max_ddd

  ! function max_dr(v1, v2) result (v3)
  !   type (DualNumber), intent(in) :: v1
  !   real,              intent(in) :: v2
  !   type (DualNumber) :: v3

  !   if(v1%f0 > v2) then
  !      v3 = v1
  !   else
  !      v3 = v2
  !   endif
  ! end function max_dr

  ! function max_rd(v1, v2) result (v3)
  !   real,              intent(in) :: v1
  !   type (DualNumber), intent(in) :: v2
  !   type (DualNumber) :: v3

  !   if(v1 > v2%f0) then
  !      v3 = v1
  !   else
  !      v3 = v2
  !   endif
  ! end function max_rd

  ! ! Min functions
  ! function min_dd(v1, v2) result (v3)
  !   type (DualNumber), intent(in) :: v1, v2
  !   type (DualNumber) :: v3

  !   if(v1%f0 < v2%f0) then
  !      v3 = v1
  !   else
  !      v3 = v2
  !   endif
  ! end function min_dd

  ! function min_dr(v1, v2) result (v3)
  !   type (DualNumber), intent(in) :: v1
  !   real,              intent(in) :: v2
  !   type (DualNumber) :: v3

  !   if(v1%f0 < v2) then
  !      v3 = v1
  !   else
  !      v3 = v2
  !   endif
  ! end function min_dr

  ! function min_rd(v1, v2) result (v3)
  !   real,              intent(in) :: v1
  !   type (DualNumber), intent(in) :: v2
  !   type (DualNumber) :: v3

  !   if(v1 < v2%f0) then
  !      v3 = v1
  !   else
  !      v3 = v2
  !   endif
  ! end function min_rd

end module Module_Autodiff
