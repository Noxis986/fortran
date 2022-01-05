program LDLt
!--------------------------------                                               
!License:GNU Fortran (Ubuntu 9.3.0-17ubuntu1~20.04) 9.3.0                       
!Authors: Jorge Alvarez Espiño                                                  
!         Juan Manuel Fornos Andrade                                            
!Last updated:30/12/21                                                          
!-------------------------------                                                
!Resolution of a linear system via LDLt method                
!-------------------------------                                                
!To Do:                                                                          
!-------------------------------   
use iso_fortran_env, only: real64
use checks
use decomposition
implicit none
real(real64), allocatable :: a(:,:), l(:,:), d(:,:), lt(:,:), x(:), b(:)
integer :: i, n
character(1) :: iter
!-------------------------
!Structure of data.dat:
!n (dimension)
!matrix a (coefficient matrix, row for line)
!column vector b (ax=b)
!------------------------
!Read dimension and initialize variables
print*,'-------------------------'
print*, 'The matrix dimension n is:'
read*, n
print*, n
allocate(a(n,n), l(n,n), d(n,n), lt(n,n), x(n), b(n))
!Read a
do i= 1,n
    read*, a(i,:)
end do
!Read vector b
do i= 1,n
    read*, b(i)
end do

!Print matrix a
print*, 'The coefficient matrix a is:'
do i=1,n
    print*,a(i,:)
end do
!Print vector b
print*, 'Vector b:'
do i=1,n
    print*,b(i)
end do
print*,'-------------------------'

!Requirements for doing LDLt method, if they are not met,
!program LDLt will stop.
call symmetric(a)
call positive(a)
print*,'-------------------------'

!LDLt decomposition, if there is an error 1E6 times epsilon, 
!program LDLt will stop.
call ldlt_decomp(a,l,d,lt)

!L*x=b, forward substitution(module decomposition).
call forward_substitution(l,b,x)
!D*b=x, diagonal resolution (module decomposition).
call diagonal_resolution(d,x,b)
!Lt*x=b, back substitution (module decomposition). x is the solution vector
call backward_substitution(lt,b,x)

print*,'-------------------------'
print*, 'Solution vector x:'
do i=1,n
    write(iter, '(I1)') i
    print *, 'x'//iter//' =', x(i)
end do 
print*, '-------------------------'
end program
