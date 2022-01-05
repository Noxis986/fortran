module checks
!--------------------------------
!License:GNU Fortran (Ubuntu 9.3.0-17ubuntu1~20.04) 9.3.0
!Authors: Jorge Alvarez Espiño
!         Juan Manuel Fornos Andrade
!Last updated:27/12/21
!-------------------------------
!Check if 'a' matrix is a symmetric and positive definite matrix 
!-------------------------------
!To Do
!-------------------------------
use iso_fortran_env, only: real64
implicit none
contains
!-------------------------------
subroutine symmetric(a)
! Checks if 'a' is a symmetric matrix
real(real64), intent(in):: a(:,:)
integer :: i,j
logical :: sym

do i=1,size(a,2)
    do j=1,size(a,2)
        sym=a(i,j)==a(j,i)
        if (sym .eqv. .false.) then
            print*, 'ERROR: The a matrix is not symmetrical'
            stop
        end if
    end do
end do
print*, 'Matrix is symmetric'

end subroutine
!--------------------------------
subroutine lu(a, n, l, u)
!LU method (done as an example in a lecture) for evaluating determinants ahead
! as det(a)=det(l)*det(u).
real(real64), intent(in) :: a(:,:)
real(real64), intent(inout) :: l(:,:), u(:,:)
integer :: i,j,n

do j=1,n
    do i= 1, n
        u(i,j)=a(i,j) - dot_product(l(i,1:i-1),u(1:i-1,j))
    end do
    l(j,j)=1._real64
    do i=j+1, n
        l(i,j)=(a(i,j)-dot_product(l(i,1:j-1),u(1:j-1,j)))/u(j,j)
    end do
end do

end subroutine
!--------------------------------
subroutine positive(a)
!Checks if 'a' is a positive definite matrix
real(real64), intent(in):: a(:,:)
real(real64), allocatable :: l(:,:), u(:,:), b(:,:)
integer :: i,j,n
real(real64) :: k

if (a(1,1)<0) then
    print *, ' ERROR: The a matrix is not a positive definite matrix, a(1,1)<0'
    stop
end if
!Determinants of the upper left submatrices, evaluation via LU decomposition
!det(a)=det(l)*det(u)
!The determinant of a triangular matrix ('l' and 'u') is equal to the product
!of its diagonal elements.
n=size(a,2)
do i=2,n
    allocate(b(i,i),l(i,i),u(i,i))
    b=a(1:i,1:i)
    call lu(b,i,l,u)
    k=1._real64
    do j=1,i
        k=k*u(j,j)*l(j,j)
        if (k<0) then
            print*, 'ERROR: a is not positive definite, determinant of upper &
                left submatrix n=',i,' <0'
            stop
        end if    
    end do
    deallocate(b,l,u)
end do
print*, 'Matrix is positive definite'

end subroutine
!-------------------------------
end module
