program Hubbard

!  use mpi
#include <slepc/finclude/slepceps.h>
  use slepceps
  use, intrinsic :: iso_fortran_env, only : iostat_end

  implicit none

  integer :: num_sites
  integer :: num_upspin
  integer :: num_dnspin
  double precision :: t1
  double precision :: t2
  double precision :: u1
  double precision :: v1

  double precision :: input_val
  character(len=30):: str1, str2
  integer :: stat

  integer , parameter :: mem_factor = 60 
  double precision , parameter :: small_val = 1E-8
  integer*8 :: num_states_up, num_states_dn, num_states_tmp
!  integer :: total_states
  integer*8 :: total_states  ! for 18 sites

  integer*8 :: m_max
  integer, allocatable :: comb(:), states_tmp(:), states_up(:), states_dn(:)
  double precision, allocatable :: T_matrix(:,:)
  double precision, allocatable :: V_matrix(:,:)
  double precision :: num_sum, sum_tmp, error
  integer*8 :: ii, jj, i, j, k, sign_tmp, tmp_j, tmp_i
  integer*8 :: counter       ! for 18 sites
  integer*8 :: tmp
  integer*8 :: key, pos
  integer*8 :: bool_up, bool_dn, bool_sum_i
!  integer :: up_index, dn_index
  integer*8 :: up_index, dn_index  ! for 18 sites

  ! for mpi, non-zero elements
  integer :: num_procs
!  integer :: l_a, l_start, l_end, l_rows
  integer*8 :: l_a, l_start, l_end, l_rows ! for 18 sites
  integer*8, allocatable :: accumu_num_nz(:)

  ! for PETSc/SLEPSc
  Vec            xr            ! for eigenvector
!  PetscScalar    kr, kr_tmp    ! for eigenvalue
  PetscComplex   kr, kr_tmp    ! for eigenvalue
  Mat            P_H_matrix
  EPSType        tname
  EPS            eps
  PetscInt       P_total_states, P_a, P_start, P_end, num_values
  PetscInt       P_tmp, loc_a
  PetscMPIInt    rank
  PetscErrorCode ierr
  PetscBool      terse
  PetscInt       nev,ncv, mpd, nconv, maxits
  PetscReal      tol, P_norm
  PetscInt, allocatable :: d_nnz(:), o_nnz(:)
  PetscInt, allocatable :: col(:), vec_indices(:)
!  PetscScalar, allocatable :: non_zeros(:), vec_val(:)
  PetscComplex, allocatable :: non_zeros(:), vec_val(:)

  ! for timer
  real :: start, finish, time1, time2

  ! for file
  character :: file1*100

  open(unit=66, file='inputs', status = 'old')
  do
     read( 66, *, iostat = stat ) str1, str2, input_val

     select case( stat )
     case( 0 )
       if(str1 == 'num_sites') num_sites = int(input_val)
       if(str1 == 'num_upspin') num_upspin = int(input_val)
       if(str1 == 'num_dnspin') num_dnspin = int(input_val)
       if(str1 == 't1') t1 = dble(input_val)
       if(str1 == 't2') t2 = dble(input_val)
       if(str1 == 'u1') u1 = dble(input_val)
       if(str1 == 'v1') v1 = dble(input_val)
     case( iostat_end )
       exit
     case default
       write( *, * ) 'Error in reading file'
       stop
     end select
  end do

  m_max = num_upspin
  if (num_dnspin > num_upspin) then
    m_max = num_dnspin
  end if

  num_states_up = Combination(num_sites,num_upspin)
  num_states_dn = Combination(num_sites,num_dnspin)

  allocate(comb(0:m_max-1)); comb = 0
  allocate(states_up(0:num_states_up-1)); states_up = 0
  allocate(states_dn(0:num_states_dn-1)); states_dn = 0

  ! Gen is used to generate basis states recorded by integer ( | 0101 > = 5 )
  num_states_tmp = num_states_up
  allocate(states_tmp(0:num_states_up-1)); states_tmp = 0
  ii = 0
  call Gen (0, num_upspin-1)
  do ii = 0 , num_states_up-1
    states_up(ii)=states_tmp(ii)
  end do
  deallocate(states_tmp)

  num_states_tmp = num_states_dn
  allocate(states_tmp(0:num_states_dn-1)); states_tmp = 0
  ii = 0
  call Gen (0, num_dnspin-1)
  do ii = 0 , num_states_dn-1
    states_dn(ii)=states_tmp(ii)
  end do
  deallocate(states_tmp)
  deallocate(comb)

  total_states = num_states_up*num_states_dn

! Generating T_matrix and V_matrix
  allocate(T_matrix(0:num_sites-1, 0:num_sites-1)); T_matrix = 0.0
  allocate(V_matrix(0:num_sites-1, 0:num_sites-1)); V_matrix = 0.0


  call MPI_Init ( ierr )
  call MPI_Comm_rank ( MPI_COMM_WORLD, rank, ierr ) ! Determine this process's rank.
  call MPI_Comm_size ( MPI_COMM_WORLD, num_procs, ierr ) ! Find out the number of processes available.

  include 'num_sites_16.f90' ! H_matrix's interaction terms has been written

!!!  call MPI_Finalize ( ierr )

  if( rank .eq. 0) then
    write(*,*)
    write(*,"(A,I5)") 'num_sites  = ', num_sites
    write(*,"(A,I5)") 'num_upspin = ', num_upspin
    write(*,"(A,I5)") 'num_dnspin = ', num_dnspin
    write(*,"(A,F10.2)") 't1 = ', t1
    write(*,"(A,F10.2)") 't2 = ', t2
    write(*,"(A,F10.2)") 'u1 = ', u1
    write(*,"(A,F10.2)") 'v1 = ', v1
    write(*,*)
  end if
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  if (rank < mod(total_states,num_procs)) then
    l_start = rank*(total_states/num_procs+1)
    l_end = l_start + total_states/num_procs+1
  else
    l_start = mod(total_states,num_procs)*(total_states/num_procs+1)+(rank-mod(total_states,num_procs))*(total_states/num_procs)
    l_end = l_start + total_states/num_procs
  end if

  tmp = total_states/num_procs
  tmp = tmp*mem_factor
  allocate(col(0:tmp-1)); col = 0
  allocate(non_zeros(0:tmp-1)); non_zeros = (0.0,0.0)

  l_rows = l_end-l_start
  allocate(d_nnz(0:l_rows-1)); d_nnz = 0
  allocate(o_nnz(0:l_rows-1)); o_nnz = 0
  allocate(vec_indices(0:l_rows-1)); vec_indices = 0
  allocate(vec_val(0:l_rows-1)); vec_val = 0
  allocate(accumu_num_nz(0:l_rows)); accumu_num_nz = 0


  call cpu_time(start)


  counter = 0
  do l_a=l_start, l_end-1 ! start from 0
    up_index = int(l_a/num_states_dn)
    dn_index = mod(l_a,num_states_dn)
    num_sum = 0.0
    ! H_matrix : for interaction terms
    do i=0, num_sites-1

      bool_up = 0
      if ( btest(states_up(up_index),i) .eqv. .true. ) bool_up = 1
      bool_dn = 0
      if ( btest(states_dn(dn_index),i) .eqv. .true. ) bool_dn = 1

      num_sum = num_sum + u1*bool_up*bool_dn  ! for hubbard u
      bool_sum_i =  bool_up + bool_dn
      sum_tmp = 0.0
      do j=i+1, num_sites-1
        if ( btest(states_up(up_index),j) .eqv. .true. ) then
          bool_up = 1
        else
          bool_up = 0
        end if
        if ( btest(states_dn(dn_index),j) .eqv. .true. ) then
          bool_dn = 1
        else
          bool_dn = 0
        end if
        sum_tmp = sum_tmp + V_matrix(i,j) * (bool_up + bool_dn)
      end do
      num_sum = num_sum + bool_sum_i*sum_tmp
    end do
    if (abs(num_sum) .gt. small_val) then
     col(counter) = l_a
     non_zeros(counter) = num_sum
     d_nnz(l_a-l_start) = d_nnz(l_a-l_start)+1
     counter = counter + 1
    end if
    ! H_matrix : hopping term
    do i=0, num_sites-1
      do j=0, num_sites-1
        ! spin_dn hopping
        if ((abs(T_matrix(i,j)) > small_val) .and. (btest(states_dn(dn_index),j) .eqv. .true.) .and. (btest(states_dn(dn_index),i) .eqv. .false.)) then
          ! determine sign (+ or -)
          sign_tmp = 1
          tmp_j = j
          tmp_i = i
          if ( j > i) then
            tmp_j = i
            tmp_i = j
          end if
          do k = tmp_j+1, tmp_i-1
            if (btest(states_dn(dn_index),k) .eqv. .true.) then
              sign_tmp = -sign_tmp
            end if
          end do
          ! generate the new state by operation and then search the position of new state
          key = states_dn(dn_index)
          key = ibclr(key, j)
          key = ibset(key, i)
          pos = Binary_search(states_dn, num_states_dn, key)
          col(counter) = up_index*num_states_dn+pos    ! start from 0
          non_zeros(counter) = sign_tmp*T_matrix(i,j)  ! Again, one should use complex conjugate for complex H.
          if( ( l_start <= col(counter) ).and. ( col(counter) <= l_end-1) )then
            d_nnz(l_a-l_start) = d_nnz(l_a-l_start)+1
          else
            o_nnz(l_a-l_start) = o_nnz(l_a-l_start)+1
          end if
          counter = counter + 1
        end if
        ! spin_up hopping
        if ((abs(T_matrix(i,j)) > small_val) .and. (btest(states_up(up_index),j) .eqv. .true.) .and. (btest(states_up(up_index),i) .eqv. .false.)) then
          ! determine sign (+ or -)
          sign_tmp = 1
          tmp_j = j
          tmp_i = i
          if ( j > i) then
            tmp_j = i
            tmp_i = j
          end if
          do k = tmp_j+1, tmp_i-1
            if (btest(states_up(up_index),k) .eqv. .true.) then
              sign_tmp = -sign_tmp
            end if
          end do
          ! generate the new state by operations and then search the position of new state
          key = states_up(up_index)
          key = ibclr(key, j)
          key = ibset(key, i)
          pos = Binary_search(states_up, num_states_up, key)
          col(counter) = pos*num_states_dn+dn_index    ! start from 0
          non_zeros(counter) = sign_tmp*T_matrix(i,j)  ! Again, one should use complex conjugate for complex H.
          if( ( l_start <= col(counter) ).and. ( col(counter) <= l_end-1) )then
            d_nnz(l_a-l_start) = d_nnz(l_a-l_start)+1
          else
            o_nnz(l_a-l_start) = o_nnz(l_a-l_start)+1
          end if
          counter = counter + 1
        end if
      end do
    end do
    accumu_num_nz(l_a-l_start+1) = counter
  end do


  deallocate(states_up)
  deallocate(states_dn)
  deallocate(T_matrix)
  deallocate(V_matrix)

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if(rank .eq. 0) write(*,*) 'Matrix construction finished!'

  call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
!  call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
!  call MPI_Comm_size ( MPI_COMM_WORLD, num_procs, ierr ) ! Find out the number of processes available.


  call MatCreate(PETSC_COMM_WORLD,P_H_matrix,ierr);CHKERRA(ierr)
  P_total_states = total_states
  call MatSetSizes(P_H_matrix,PETSC_DECIDE,PETSC_DECIDE,P_total_states,P_total_states,ierr);CHKERRA(ierr)
  call MatSetType(P_H_matrix, MATMPIAIJ, ierr);CHKERRA(ierr)

!  call MatSetUp(P_H_matrix,ierr);CHKERRA(ierr)
  P_tmp = 0
  call MatMPIAIJSetPreallocation(P_H_matrix, P_tmp, d_nnz, P_tmp, o_nnz, ierr);CHKERRA(ierr)
  call MatGetOwnershipRange(P_H_matrix,P_start,P_end,ierr);CHKERRA(ierr)

  deallocate(d_nnz)
  deallocate(o_nnz)

  if ((P_start .ne. l_start) .or. (P_end .ne. l_end)) write(*,*) 'Inconsistency!'

  P_tmp = 1
  do P_a = P_start, P_end-1
    loc_a = P_a-P_start
    num_values = accumu_num_nz(loc_a+1)-accumu_num_nz(loc_a)
    call MatSetValues(P_H_matrix, P_tmp, P_a, num_values, col(accumu_num_nz(loc_a):accumu_num_nz(loc_a)-1), &
                    & non_zeros(accumu_num_nz(loc_a):accumu_num_nz(loc_a)-1), INSERT_VALUES, ierr);CHKERRA(ierr)
  end do

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if(rank .eq. 0) write(*,*) 'Petsc: Matrix elements population finished!'

  deallocate(col)
  deallocate(non_zeros)
  deallocate(accumu_num_nz)

  call MatAssemblyBegin(P_H_matrix, MAT_FINAL_ASSEMBLY, ierr);CHKERRA(ierr)
  call MatAssemblyEND(P_H_matrix, MAT_FINAL_ASSEMBLY, ierr);CHKERRA(ierr)

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if(rank .eq. 0) write(*,*) 'Petsc: Matrix assembly finished!'

  call cpu_time(finish)
  time1 = finish-start

! diagonalization

  call cpu_time(start)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create the eigensolver and display info
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** Create eigensolver context
      call EPSCreate(PETSC_COMM_WORLD,eps,ierr);CHKERRA(ierr)

!     ** Set operators. In this case, it is a standard eigenvalue problem
      call EPSSetOperators(eps,P_H_matrix,PETSC_NULL_MAT,ierr);CHKERRA(ierr)
      call EPSSetProblemType(eps,EPS_HEP,ierr);CHKERRA(ierr)
      call EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL, ierr);CHKERRA(ierr)

      call EPSSetType(eps,EPSKRYLOVSCHUR, ierr);CHKERRA(ierr)
!      call EPSSetType(eps,EPSGD, ierr);CHKERRA(ierr)
!      call EPSSetType(eps,EPSARNOLDI, ierr);CHKERRA(ierr)
!      call EPSSetType(eps,EPSLANCZOS, ierr);CHKERRA(ierr)
!      call EPSSetType(eps,EPSPOWER , ierr);CHKERRA(ierr)
!      call EPSSetType(eps,EPSSUBSPACE , ierr);CHKERRA(ierr)
!      call EPSSetType(eps,EPSJD , ierr);CHKERRA(ierr)
!      call EPSSetType(eps,EPSRQCG , ierr);CHKERRA(ierr)
!      call EPSSetType(eps,EPSLOBPCG , ierr);CHKERRA(ierr)
!      call EPSSetType(eps,EPSCISS , ierr);CHKERRA(ierr)
!The parameters ncv and mpd are intimately related, so that the user is advised
!to set one of them at most. Normal usage is that (a) in cases where nev is
!small, the user sets ncv (a reasonable default is 2*nev); and (b) in cases
!where nev is large, the user sets mpd.

!The value of ncv should always be between nev and (nev+mpd), typically
!ncv=nev+mpd. If nev is not too large, mpd=nev is a reasonable choice, otherwise
!a smaller value should be used.

      mpd=30 
      nev=20
      ncv=50
      call EPSSetDimensions(eps, nev, ncv, mpd, ierr)

      tol=1.0d-8
      maxits=500 
      call EPSSetTolerances(eps, tol, maxits, ierr)

!     ** Set solver parameters at runtime
      call EPSSetFromOptions(eps,ierr);CHKERRA(ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Solve the eigensystem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call EPSSolve(eps,ierr);CHKERRA(ierr)

!     ** Optional: Get some information from the solver and display it
      call EPSGetType(eps,tname,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,120) tname
      endif
 120  format (' Solution method: ',A)
      call EPSGetDimensions(eps,nev,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,130) nev
      endif
 130  format (' Number of requested eigenvalues:',I4)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** show detailed info unless -terse option is given by user
      call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-terse',terse,ierr);CHKERRA(ierr)
  if (terse) then
    call EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_NULL_VIEWER,ierr);CHKERRA(ierr)  ! show eigenvalues and their errors
  else
    call PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL,ierr);CHKERRA(ierr)
    call EPSReasonView(eps,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRA(ierr)
    call EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRA(ierr)
    call PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRA(ierr)

!==================  For printing eigenvector  ==============================================================================!

    P_tmp = 0  ! 0 means the ground state
    counter = 0 ! for skipping results with large error
    
    call MatCreateVecs( P_H_matrix, xr, PETSC_NULL_VEC, ierr )    !  Preparing a vec for an eigenvector
    call EPSGetEigenpair(eps, counter, kr, PETSC_NULL_SCALAR, xr, PETSC_NULL_VEC, ierr)
    call EPSComputeError(eps, counter, EPS_ERROR_RELATIVE, error, ierr)
      if( rank .eq. 0) write(*,*) 'eigenval: ', kr, 'error: ', error
    do while( error .gt. small_val)
      counter = counter + 1
      call EPSComputeError(eps, counter, EPS_ERROR_RELATIVE, error, ierr)
      if( rank .eq. 0) write(*,*) 'error: ', error
    end do

    kr_tmp = kr

    do while( abs(kr_tmp-kr) .lt. small_val )
!        call VecView(xr, PETSC_VIEWER_STDOUT_WORLD, ierr)                                 ! show eigenvector
      !  Preparing a vec for an eigenvector

!      write(*,*) P_tmp, kr

      do P_a = P_start, P_end-1
        vec_indices(P_a-P_start) = P_a
      end do

      call EPSGetEigenpair(eps, counter, kr, PETSC_NULL_SCALAR, xr, PETSC_NULL_VEC, ierr)
      call VecGetValues(xr, P_end-P_start, vec_indices, vec_val, ierr)
      call VecNorm(xr, NORM_2, P_norm, ierr)
      if (rank .eq. 0)    write(*,*)  kr, 'norm2 =', P_norm
      write(unit=file1, fmt="('vec_',i2.2,'_',i5.5,'.txt')") P_tmp, rank
      open(unit=111+rank,file=file1)
      do jj=0, P_end-P_start-1
        write(111+rank, *) vec_val(jj)
      enddo
      close(111+rank)

      kr_tmp = kr
      P_tmp = P_tmp + 1
      counter = counter + 1
      call EPSGetEigenpair(eps, counter, kr, PETSC_NULL_SCALAR, xr, PETSC_NULL_VEC, ierr)

    enddo

  endif


  call EPSDestroy(eps,ierr);CHKERRA(ierr)
  call MatDestroy(P_H_matrix,ierr);CHKERRA(ierr)

  call SlepcFinalize(ierr)

  call cpu_time(finish)
  time2 = finish-start


  if (rank .eq. 0) then

    write(*,9) "Time for H_matrix        ", int(time1/60/60), " hr ", mod(int(time1/60.0),60), " min ", mod(int(time1),60), " sec "
    write(*,9) "Time for diagonalization ", int(time2/60/60), " hr ", mod(int(time2/60.0),60), " min ", mod(int(time2),60), " sec "
9 FORMAT(A, I5, A, I2, A, I2, A)  

  endif



  call MPI_Finalize ( ierr )

  include 'contains.f90'

end program Hubbard

