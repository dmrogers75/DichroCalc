program coupling
!* The diabatisation scheme works well only in the case where the ref states are very
!  similar to the diabatic state of interest. 
!  Also, the choice of the chemical properties is also crucial.
!*
  ! Units have to be consistent. Stay with Atomic Units? Seems to work fine with
  ! cm^-1 and debye.
  ! Scale the length of iso to match adia? Seems to work fine.
!*
  implicit none
  integer :: ndim   ! number of vector elements (number of rows)
  integer :: ntrans ! number of columns (one trans per column)
  double precision, allocatable :: mu_iso(:,:), mu_a(:,:), mu_d(:,:), mu_iso_save(:)
  !---transition dipole moments for isolated, adiabatic and diabatic states
  double precision, allocatable :: c(:,:), ct(:,:), h_a(:,:), h_d(:,:), evec(:,:), eval(:)
  !=====================================================================
  !---variables for Singular value decomposition
  integer, parameter :: lwmax=1000
  ! ndim is actually the number of rows. For TDM, this is normally 3 (x, y , z components)
  integer :: m, n, lda, ldu, ldvt, info, lwork
  double precision :: work(lwmax)
  double precision, allocatable :: a(:,:), u(:,:), vt(:,:), s(:)
  !=====================================================================
  double precision :: scal_fact, theta1, theta2, r1i, r1j, r2i, r2j, norm_iso, norm_a, error 
  double precision :: dr, sumdiff, sumdiff_save, acc_rate
  integer :: i, j, n_arg, clen, step, nstep, nacc, ntrial
  character (len=32) :: arg
  character(len=:), allocatable :: arg1
 !
  !---Read all command line arguments
  n_arg = 0
  do
    call get_command_argument(n_arg, arg)
    if (len_trim(arg) == 0) exit
    if (n_arg == 1) then
      clen=len_trim(arg)
      allocate(character(len=clen) :: arg1)
      arg1=arg(1:clen)
      !write(6,*) 'length of argument ', clen, arg1
    end if
    !write(6,*) n_arg, trim(arg)
    n_arg = n_arg + 1
  end do
 !
  !---check if an input is given (also check if more than 1 input is given)
  if (n_arg <= 1) then
    write(6,*) 'Need an input?'
    write(6,*) './name.exe input'
    stop
  else
    if (n_arg > 2) then
      write(6,*) 'Warning: ', n_arg-1, ' inputs are given. Only need one!'
      WRITE(6,*) 'Warning: First argument will be used!'
    end if
  end if
 !
  !---read input (mu_iso, mu_a and h_a)
  open(1,file=arg1)
  read(1,*) ! title (first line)
  read(1,*) ndim !(second line)
  read(1,*) ! title (third line)
  read(1,*) ntrans !(fourth line)
  !write(6,*) ntrans
  !---Once ntrans is read we have to allocate memory for arrays before carrying
  !on reading input
  allocate(mu_iso(ndim,ntrans), mu_a(ndim,ntrans), mu_d(ndim,ntrans))
  allocate(c(ntrans,ntrans), ct(ntrans,ntrans))
  allocate(h_a(ntrans,ntrans), h_d(ntrans,ntrans))
  allocate(evec(ntrans,ntrans), eval(ntrans))
  !---now we can carry on reading things in
  read(1,*) ! title !
  do i=1,ndim ! x y z
    read(1,*) mu_iso(i,:) ! (lines 4-6)
    !write(6,*) i, mu_iso(i,:)
  end do
  read(1,*) ! title (line 7)
  do i=1,ndim ! x y z
    read(1,*) mu_a(i,:) ! (line 8-10)
    !write(6,*) i, mu_a(i,:)
  end do
  read(1,*) ! title ! (line 11)
  do i=1,ntrans
    read(1,*) h_a(i,:) ! (line 12-13)
    !write(6,*) i, h_a(i,:)
  end do
  !---convert from nanometer to wave number
  do i=1,ntrans
    !h_a(i,i)=10d6/h_a(i,i) ! cm^-1
    !h_a(i,i)=h_a(i,i)/27.211 ! a.u
  end do
  !write(6,*)
  !do i=1,ntrans
  !  write(6,*) h_a(i,:)
  !end do
 !
  !write(6,*) 'mu_iso '
  !do i=1,ndim
  !  write(6,'(4(F10.4))') mu_iso(i,:)
  !end do
  !write(6,*) 'norm iso ', get_norm(mu_iso)
 !
  !write(6,*) 'mu_a '
  !do i=1,ndim
  !  write(6,'(4(F10.4))') mu_a(i,:)
  !end do
  !write(6,*) 'norm adi ', get_norm(mu_a)
  !write(6,*)
 !
  !---scale mu_a to match the length of m_iso
  norm_iso=get_norm(mu_iso)
  norm_a=get_norm(mu_a)
  scal_fact=sqrt(norm_iso)/sqrt(norm_a)
  mu_a=mu_a*scal_fact
  !write(6,*) 'norm adi new', get_norm(mu_a)  
 !
  !---allocate memory for arrays
  m=ntrans; n=ntrans
  lda=m; ldu=m; ldvt=n
  allocate(a(lda,n), u(ldu,m), vt(ldvt,n), s(min(m,n)) )
  !---make matrix a
  a=matmul(transpose(mu_a),mu_iso)
  call svd()
  ct=matmul(u,vt)
  c=transpose(ct)
  mu_d=matmul(mu_a,ct)
 !
  !---write out mu_iso and mu_d for comparison. In the ideal case they should be
  !very simular.
  !write(6,*) '===============write out mu_iso and mu_d for comparison=================='
  !write(6,*) 'mu_iso '
  !do i=1,ndim
  !  write(6,'(4(F10.4))') mu_iso(i,:)
  !end do
  !write(6,*) 'mu_d '
  !do i=1,ndim
  !  write(6,'(4(F10.4))') mu_d(i,:)
  !end do
  !write(6,*) '===============write out mu_iso and mu_d for comparison==================='
  !write(6,*)
 !
  !---unitary matrix
  !write(6,*) 'unitary matrix c'
  !do i=1,ntrans
  !  write(6,'(6(f18.4))') c(i,:)
  !end do
  !write(6,*)
  !---calculate the angles between final results and those from iso. In the
  !ideal case these angles should be zero.
  do i=1,ntrans
    theta1=get_costheta(mu_iso(:,i),mu_d(:,i))
    !write(6,*) 'angle ', i, theta1, acos(theta1)/acos(-1d0)*180d0
  end do
  !write(6,*)
  !
  !---unitary transformation the adiabatic hamiltonion to the diabatic one
  h_d=matmul(c,matmul(h_a,ct)) ! au
   write(6,*) 'diabatic hamiltonian (with couplings)'
   do i=1,4
     write(6,'(6(f18.4))') h_d(i,:)
   end do
  !write(6,*)
  !write(6,'(A10, F10.4)') "J_1xy", h_d(1,2)*1000d0
  !write(6,'(A10, F10.4)') "J_xx", h_d(1,3)*1000d0
  !write(6,'(A10, F10.4)') "J_xy", h_d(1,4)*1000d0
  !write(6,'(A10, F10.4)') "J_yx", h_d(2,3)*1000d0
  !write(6,'(A10, F10.4)') "J_yy", h_d(2,4)*1000d0
  !write(6,'(A10, F10.4)') "J_2xy", h_d(3,4)*1000d0
  write(6,'(F10.4)') h_d(1,2)
 !
  !---diagonalize the diabatic to get back to the adiabatic state (mainly for
  !housekeeping)
  evec=h_d
  call diamat(evec,eval,ntrans)
  !---eigenvector (or the unitary matrix)
  !write(6,*) 'eigenvectors '
  !do i=1,ntrans
  !  write(6,'(6(f18.4))') evec(i,:)
  !end do
  !write(6,*)
 !
  !write(6,*) 'adiabatic hamiltonian after diagonalize the diabatic one'
  h_a=0d0
  do i=1,ntrans
    h_a(i,i)=eval(i)
  end do
  !---convert h_d from au to cm^-1
  !do i=1,m
  !  write(6,'(6(f18.4))') h_a(i,:)
  !end do
!*
contains
!*
  subroutine svd()
    implicit none
    !---query the optimal workspace.
    lwork = -1
    call dgesvd( 'all', 'all', m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
    lwork = min( lwmax, int( work( 1 ) ) )
    !---compute svd.
    call dgesvd( 'all', 'all', m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info )
    !---check for convergence.
    if( info.gt.0 ) then
      write(*,*)'the algorithm computing svd failed to converge.'
      stop
    end if
  end subroutine svd
 !
!*
  subroutine diamat(a,eig,n)                                                                                                                
    implicit none                                                                                                                           
    integer :: n,l,inf                                                                                                                      
    double precision ::  a(n,n),eig(n),work(n*(3+n/2))                                                                                      
    l=n*(3+n/2)                                                                                                                             
    !---this routine is from lapack,                                                                                                        
    !so remember to compile with -llapack                                                                                                   
    call dsyev('V','U',n,a,n,eig,work,l,inf)                                                                                                
  end subroutine diamat                                                                                                                     
!*
  double precision function get_norm(mat)
    double precision, intent(in) :: mat(:,:)
    get_norm=0d0
    do i=1,ndim
      do j=1,ntrans
        get_norm=get_norm+mat(i,j)**2
      end do
    end do
  end function get_norm
!*
  double precision function get_length(v)
    double precision, intent(in) :: v(ndim)
    get_length=sqrt( dot_product(v,v) )
  end function get_length
!*
  double precision function get_costheta(v1,v2)
    double precision, intent(in) :: v1(ndim), v2(ndim)
    get_costheta=dot_product(v1,v2)
    get_costheta=get_costheta/get_length(v1)
    get_costheta=get_costheta/get_length(v2)
  end  function get_costheta
!*
end program coupling
