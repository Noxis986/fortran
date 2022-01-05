module decomposition
!--------------------------------
!License:GNU Fortran (Ubuntu 9.3.0-17ubuntu1~20.04) 9.3.0
!Authors: Jorge Alvarez Espiño
!         Juan Manuel Fornos Andrade
!Last updated:29/12/21
!-------------------------------
!Module with LDLt decomposition and solving procedure for diagonal, upper and 
!lower matrices
!-------------------------------
!To Do:
!-------------------------------
use iso_fortran_env, only: real64
implicit none
contains
!-------------------------------
subroutine ldlt_decomp(a, l, d, lt)
!LDLt decomposition
real(real64), intent(in) :: a(:,:)
real(real64), intent(inout) :: l(:,:), d(:,:), lt(:,:)
real(real64) :: error(size(a,2),size(a,2)), error2(size(a,2),size(a,2)) 
real(real64) :: s, m
integer i, j, k
!Calculate L and D
do i=1, size(a,2)
    l(i,i)=1
    s=0._real64
    do k=1,i-1
        s=s+(l(i,k)**2)*d(k,k)
    end do
    d(i,i)=a(i,i) - s
    do j=i+1, size(a,2)
        m=0._real64
        do k=1,i-1
            m=m+l(j,k)*l(i,k)*d(k,k)
        end do    
        l(j,i)=(a(j,i) - m)/d(i,i)
    end do
end do

!Transpose L
do i=1,size(a,2)
    do j=1,size(a,2)
        lt(i,j)=l(j,i)
    end do
end do

!print L,D and Lt matrices
print*, 'Matrix L: '
print*,' '
do i=1,size(a,2) 
    print*, l(i,:)
end do
print*, 'Matrix D: '
print*,' '
do i=1,size(a,2) 
    print*, d(i,:)
end do
print*, 'Matrix Lt: '
print*,' '
do i=1,size(a,2) 
    print*, lt(i,:)
end do

!Simple error check 
error=matmul(matmul(l,d),lt)
error2= a-error
do i=1, size(a,2)
    do j=1, size(a,2)
        if (a(i,j) /= error(i,j)) then
            print*, 'ERROR: a - L*D*Lt = '
            do k= 1,size(a,2)
                print *, error2(i,:)
            end do
            if (abs(error2(i,j))>1E10*epsilon(error2(i,j))) then
                stop 'Error too big to ignore'
            else
                continue
            end if 
        end if
    end do
end do
print*,' ' 
print*, 'LDLt decomposition done successfully'
end subroutine
!-----------------------------        
subroutine diagonal_resolution(a,b,x)
!Solution for a diagonal matrix
real(real64), intent(in) :: a(:,:)
real(real64), intent(in) :: b(:)
real(real64), allocatable :: x(:)
integer i

do i=1, size(a,2)
    x(i)=b(i)/a(i,i)
end do
end subroutine
!----------------------------
subroutine backward_substitution(a,b,x)
!Backward substitution procedure for solving a linear system where 'a' is an
!upper triangular matrix
real(real64), intent(in) :: a(:,:)
real(real64), intent(in) :: b(:)
real(real64), allocatable :: x(:)
integer :: i, n
n=size(a,2)

x(n)=b(n)/a(n,n)
do i= n-1, 1, -1
 x(i)= (b(i)-dot_product(a(i,i+1:n),x(i+1:n)))/a(i,i)
end do
end subroutine
!-----------------------------
subroutine forward_substitution(a,b,x)
!Forward substitution for a lower matrix
real(real64), intent(in) :: a(:,:)
real(real64), intent(in) :: b(:)
real(real64), allocatable :: x(:)
integer :: i, n
n=size(a,2)

x(1)=b(1)/a(1,1)
do i= 2, n
    x(i) = (b(i)-dot_product(a(i,1:i-1),x(1:i-1)))/a(i,i)
end do
end subroutine
!----------------------------
end module
