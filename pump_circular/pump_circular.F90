program pump_circular

!  use mpi
#include <slepc/finclude/slepceps.h>
#include <slepc/finclude/slepcmfn.h>
  use slepceps
  use slepcmfn
  use, intrinsic :: iso_fortran_env, only : iostat_end

  implicit none

  integer :: num_sites
  integer :: num_upspin
  integer :: num_dnspin
  double precision :: t1
  double precision :: t2
  double precision :: u1
  double precision :: v1
  integer :: num_files_per_vec
  integer :: degeneracy

  double precision :: input_val
  character(len=30):: str1, str2
  integer :: stat

  integer , parameter :: mem_factor = 60 
  double precision , parameter :: small_val = 1D-8
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
  integer   :: f_start, f_end
  integer*8, allocatable :: accumu_num_nz(:)

  integer, allocatable :: rx(:), ry(:), pair_sites(:,:)
  double precision :: pi, qx, qy, sum_square, Cdw, Szz
  double precision :: tmp_norm, del_s, del_s_star, del_d_square, del_px, del_py, del_s_s_star
  double precision :: del_px_py, del_px_m_py, del_px_ipy, del_px_m_ipy

  ! for reading files
  character :: file1*100, file2*100
  integer, allocatable :: num_lines(:), num_files(:)
  complex*8 :: eye, sum_i

  ! for PETSc/SLEPSc
  Vec            xr_0, xr_1     ! for eigenvector
!  PetscScalar    kr, kr_tmp    ! for eigenvalue
  PetscComplex   kr, kr_tmp     ! for eigenvalue
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
  PetscComplex, allocatable :: non_zeros(:), vec_val(:)
  PetscComplex, allocatable :: coeff(:)

  ! for timer
  real :: start, finish, time1, time2

  ! for time evolution
  MFN              mfn            ! MFN solver context 
  FN               f              ! the function, exp() in this example 
  double precision dt, t, t_max, t_min   ! alpha = -i*dt
  PetscScalar      alpha          ! to compute beta*exp(alpha*H) 
  PetscScalar      beta
  PetscScalar      P_de

  ! For new states
  integer :: num_upspin_1, num_dnspin_1
  integer*8 :: total_states_1
  integer*8 :: num_states_up_1, num_states_dn_1
  integer, allocatable :: states_up_1(:), states_dn_1(:)
  double complex, allocatable :: coeff_1(:), del_tmp(:), del(:,:)
!  double precision, allocatable :: coeff_1(:), del_tmp(:), del(:,:)
  integer*8 :: key_up, key_dn, pos_up, pos_dn, new_index
  integer :: sign1, sign_count

  ! For pump
  double precision :: A0, t0, sigma, omega
  double precision :: Ax, Ay, phi
  integer :: x_ij, y_ij
  integer :: deg_tmp

  t_min = -1D1
  t_max = 5D1
  dt = 5D-2
  alpha = (0.0, -5D-2)
  beta = (1.0, 0.0)

  open(unit=66, file='inputs_circular', status = 'old')
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
       if(str1 == 'num_files_per_vec') num_files_per_vec = int(input_val)
       if(str1 == 'degeneracy') degeneracy = int(input_val)
       if(str1 == 'A0') A0 = dble(input_val)
       if(str1 == 't0') t0 = dble(input_val)
       if(str1 == 'sigma') sigma = dble(input_val)
       if(str1 == 'omega') omega = dble(input_val)
     case( iostat_end )
       exit
     case default
       write( *, * ) 'Error in reading file'
       stop
     end select
  end do

  pi = 3.141592653589793
  eye = ( 0.0, 1D0)

  allocate(rx(0:num_sites-1)); rx = 0.0
  allocate(ry(0:num_sites-1)); ry = 0.0
  include 'lattice_16.f90'

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
    write(*,"(A,F10.2)") 'A0 = ', A0
    write(*,"(A,F10.2)") 't0 = ', t0
    write(*,"(A,F10.2)") 'sigma = ', sigma
    write(*,"(A,F10.2)") 'omega = ', omega
    write(*,*)
  end if

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


!  deallocate(states_up)
!  deallocate(states_dn)
!  deallocate(T_matrix)
!  deallocate(V_matrix)

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

!  deallocate(d_nnz)
!  deallocate(o_nnz)

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

!  deallocate(col)
!  deallocate(non_zeros)
!  deallocate(accumu_num_nz)

  call MatAssemblyBegin(P_H_matrix, MAT_FINAL_ASSEMBLY, ierr);CHKERRA(ierr)
  call MatAssemblyEND(P_H_matrix, MAT_FINAL_ASSEMBLY, ierr);CHKERRA(ierr)

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if(rank .eq. 0) write(*,*) 'Petsc: Matrix assembly finished!'

  call cpu_time(finish)
  time1 = finish-start

  ! For pairings
  num_upspin_1 = num_upspin - 1 
  num_dnspin_1 = num_dnspin - 1

  m_max = num_upspin_1
  if (num_dnspin_1 > num_upspin_1) then
    m_max = num_dnspin_1
  end if

  num_states_up_1 = Combination(num_sites,num_upspin_1)
  num_states_dn_1 = Combination(num_sites,num_dnspin_1)
  total_states_1 = num_states_up_1 * num_states_dn_1
!  allocate(coeff_1(0:total_states_1-1)); coeff_1 = 0.0
  allocate(coeff_1(0:total_states_1-1)); coeff_1 = (0.0, 0.0)
!  if ( rank .eq. 0) allocate(del(0:4, 0:total_states_1-1)); del = 0.0
  if ( rank .eq. 0) allocate(del(0:4, 0:total_states_1-1)); del = (0.0, 0.0)

  allocate(comb(0:m_max-1)); comb = 0
  allocate(states_up_1(0:num_states_up_1-1)); states_up_1 = 0
  allocate(states_dn_1(0:num_states_dn_1-1)); states_dn_1 = 0

  ! Gen is used to generate basis states recorded by integer ( | 0101 > = 5 )
  num_states_tmp = num_states_up_1
  allocate(states_tmp(0:num_states_tmp-1)); states_tmp = 0
  ii = 0
  call Gen (0, num_upspin_1-1)
  do ii = 0 , num_states_up_1-1
    states_up_1(ii)=states_tmp(ii)
  end do
  deallocate(states_tmp)

  num_states_tmp = num_states_dn_1
  allocate(states_tmp(0:num_states_tmp-1)); states_tmp = 0
  ii = 0
  call Gen (0, num_dnspin_1-1)
  do ii = 0 , num_states_dn_1-1
    states_dn_1(ii)=states_tmp(ii)
  end do
  deallocate(states_tmp)
  deallocate(comb)


  allocate(pair_sites(0:num_sites-1, 0:4)); pair_sites = 0
  !include 'pairing_sites_16.f90'

  call MPI_Barrier(MPI_COMM_WORLD, ierr)


! read initial vector

  if (rank < mod(num_files_per_vec,num_procs)) then
    f_start = rank*(num_files_per_vec/num_procs+1)
    f_end = f_start + num_files_per_vec/num_procs+1
  else
    f_start =mod(num_files_per_vec,num_procs)*(num_files_per_vec/num_procs+1)+(rank-mod(num_files_per_vec,num_procs))*(num_files_per_vec/num_procs)
    f_end = f_start + num_files_per_vec/num_procs
  end if

  P_tmp = total_states/num_files_per_vec
  allocate(num_lines(0:num_files_per_vec-1)); num_lines = P_tmp
  P_tmp = mod(total_states,num_files_per_vec)
  do ii = 0, P_tmp-1
    num_lines(ii) = num_lines(ii) + 1
  end do

  l_start = 0
  do ii=0, f_start-1
    l_start = l_start + num_lines(ii)
  end do

  l_end = 0
  do ii=0, f_end-1
    l_end = l_end + num_lines(ii)
  end do

  allocate(coeff(0:l_end-l_start-1)); coeff = (0.0, 0.0)
  if(rank .eq. 0) allocate(del_tmp(0:total_states_1-1)); del_tmp = (0.0, 0.0)

! prepare vectors for time-evolution


  call MatCreateVecs( P_H_matrix, xr_0, PETSC_NULL_VEC, ierr )    !  Preparing a vec for an eigenvector
  call MatCreateVecs( P_H_matrix, xr_1, PETSC_NULL_VEC, ierr )  

  do l_a = l_start, l_end-1
    vec_indices(l_a-l_start) = l_a
  end do

  if( rank .eq. 0 ) write(*,*) 'alpha = ', alpha

  do deg_tmp = 0, degeneracy-1
    counter = 0
    if(rank .eq. 0) then
      write(unit=file2, fmt="('Time_evolution_vec_',i2.2'')") deg_tmp
      open(unit=1,file=file2,status="new")  
      write(1,2002,advance="no")
2002 format ('    t      Cdw         Sdw        s        s_star    d_square      px         py        px+py      px-py     px+ipy     px-ipy'/)
!    t      Cdw         Sdw        s        s_star    d_square      px         py        px+py      px-py     px+ipy     px-ipy
!-10.00   0.317379   0.148253   0.631937   6.962271   4.467057   1.834872   2.372745   2.103808   2.372745   2.372745   2.372745
    end if

    do i = f_start, f_end-1 
      write(unit=file1, fmt="('vec_',i2.2,'_',i5.5,'.txt')") deg_tmp, i
      open(unit=111+rank,file=file1,status="old")
      do ii=0, num_lines(i)-1
        read(111+rank,*) coeff(counter)
        counter = counter+1
      end do
      close(111+rank)
    end do

    call MPI_AllReduce( dot_product(coeff,coeff), tmp_norm, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
    tmp_norm = sqrt(tmp_norm)
    if( rank .eq. 0 ) write(*, fmt="('vec ', i2, ' norm2 = ', f9.6)") deg_tmp, tmp_norm
    if( abs(tmp_norm-1.0) .gt. small_val) stop

!  here, do the initial measurements

    include 'for_measurements_final.f90'


    if(rank .eq. 0) then

      write(1,2003,advance="no") t_min, Cdw, Szz, del_s, del_s_star, del_d_square, del_px, del_py, del_px_py, del_px_m_py, del_px_ipy, del_px_m_ipy
2003 format (f6.2,2x,11(f9.6,2x)/)

    end if

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    ! finish initial measurements

    ! prepare vectors for time evolution
    call VecSetValues(xr_0, l_rows, vec_indices, coeff, INSERT_VALUES, ierr)
    call VecAssemblyBegin(xr_0, ierr)
    call VecAssemblyEnd(xr_0, ierr)

    call VecNorm(xr_0, NORM_2, P_norm, ierr)
    if( rank .eq. 0 ) write(*, fmt="('vec ', i2, '  xr_0_i  norm2 = ', f9.6)") deg_tmp, P_norm


    ! use MFN
    call MFNCreate( PETSC_COMM_WORLD, mfn, ierr)
    call MFNSetOperator( mfn, P_H_matrix, ierr)
    call MFNGetFN( mfn, f, ierr)
    call FNSetType( f, FNEXP, ierr)
    call FNSetScale( f, alpha, beta, ierr)
    call MFNSetFromOptions( mfn, ierr)



    do t = t_min+dt, t_max, dt

! =========================================================== time-evolved  H_Matrix   START
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
          ! For hopping phase
          Ax = A0*exp(-(t-t0)**2/(2*sigma**2))*cos(omega*(t-t0))
          Ay = A0*exp(-(t-t0)**2/(2*sigma**2))*sin(omega*(t-t0))
          x_ij = rx(i)-rx(j)
          y_ij = ry(i)-ry(j)
          if(x_ij == -3) x_ij =  1
          if(x_ij ==  3) x_ij = -1
          if(y_ij == -3) y_ij =  1
          if(y_ij ==  3) y_ij = -1
          phi = Ax*x_ij + Ay*y_ij
          non_zeros(counter) = conjg( sign_tmp*T_matrix(i,j)*exp(eye*phi) )! Again, one should use complex conjugate for complex H. (because row-wise settings in PETSC)
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
          ! For hopping phase
          Ax = A0*exp(-(t-t0)**2/(2*sigma**2))*cos(omega*(t-t0))
          Ay = A0*exp(-(t-t0)**2/(2*sigma**2))*sin(omega*(t-t0))
          x_ij = rx(i)-rx(j)
          y_ij = ry(i)-ry(j)
          if(x_ij == -3) x_ij =  1
          if(x_ij ==  3) x_ij = -1
          if(y_ij == -3) y_ij =  1
          if(y_ij ==  3) y_ij = -1
          phi = Ax*x_ij + Ay*y_ij
          non_zeros(counter) = conjg( sign_tmp*T_matrix(i,j)*exp(eye*phi) )! Again, one should use complex conjugate for complex H. (because row-wise settings in PETSC)
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


!  deallocate(states_up)
!  deallocate(states_dn)
!  deallocate(T_matrix)
!  deallocate(V_matrix)

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
!  if(rank .eq. 0) write(*,*) 'Matrix construction finished!'

!  call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
!  call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
!  call MPI_Comm_size ( MPI_COMM_WORLD, num_procs, ierr ) ! Find out the number of processes available.


!  call MatCreate(PETSC_COMM_WORLD,P_H_matrix,ierr);CHKERRA(ierr)
!  P_total_states = total_states
!  call MatSetSizes(P_H_matrix,PETSC_DECIDE,PETSC_DECIDE,P_total_states,P_total_states,ierr);CHKERRA(ierr)
!  call MatSetType(P_H_matrix, MATMPIAIJ, ierr);CHKERRA(ierr)

!  call MatSetUp(P_H_matrix,ierr);CHKERRA(ierr)
!  P_tmp = 0
!  call MatMPIAIJSetPreallocation(P_H_matrix, P_tmp, d_nnz, P_tmp, o_nnz, ierr);CHKERRA(ierr)
!  call MatGetOwnershipRange(P_H_matrix,P_start,P_end,ierr);CHKERRA(ierr)

!  deallocate(d_nnz)
!  deallocate(o_nnz)

!  if ((P_start .ne. l_start) .or. (P_end .ne. l_end)) write(*,*) 'Inconsistency!'

  P_tmp = 1
  do P_a = P_start, P_end-1
    loc_a = P_a-P_start
    num_values = accumu_num_nz(loc_a+1)-accumu_num_nz(loc_a)
    call MatSetValues(P_H_matrix, P_tmp, P_a, num_values, col(accumu_num_nz(loc_a):accumu_num_nz(loc_a)-1), &
                    & non_zeros(accumu_num_nz(loc_a):accumu_num_nz(loc_a)-1), INSERT_VALUES, ierr);CHKERRA(ierr)
  end do

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
!  if(rank .eq. 0) write(*,*) 'Petsc: Matrix elements population finished!'

!  deallocate(col)
!  deallocate(non_zeros)
!  deallocate(accumu_num_nz)

  call MatAssemblyBegin(P_H_matrix, MAT_FINAL_ASSEMBLY, ierr);CHKERRA(ierr)
  call MatAssemblyEND(P_H_matrix, MAT_FINAL_ASSEMBLY, ierr);CHKERRA(ierr)

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
!  if(rank .eq. 0) write(*,*) 'Petsc: Matrix assembly finished!'



! =========================================================== time-evolved  H_Matrix   END










    ! vector evolves with one timestep
      call MFNSolve( mfn, xr_0, xr_1, ierr)  !  input: xr_0  ;  output: xr_1

    ! here for time evolution measurements
    ! get coeffs from xr_1
      call VecGetValues(xr_1, l_rows, vec_indices, coeff, ierr)


      include 'for_measurements_final.f90'

    if(rank .eq. 0) then

        write(1,2003,advance="no") t, Cdw, Szz, del_s, del_s_star, del_d_square, del_px, del_py, del_px_py, del_px_m_py, del_px_ipy, del_px_m_ipy
!2003 format (f6.2,2x,10(f9.6,2x)/)

      end if

      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      !  finish one timestep measuement

      call VecCopy( xr_1, xr_0, ierr)        !  xr_0 <= xr_1
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end do

    call MFNDestroy( mfn, ierr)
    close(1)


  end do


  call MPI_Finalize ( ierr )

  include 'contains.f90'

end program pump_circular

