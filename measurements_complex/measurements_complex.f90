program measurements

  use mpi
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

  double precision , parameter :: small_val = 1D-8
  integer*8 :: m_max
  integer, allocatable :: comb(:), states_tmp(:), states_up(:), states_dn(:)
  integer*8 :: num_states_up, num_states_dn, num_states_tmp
  integer*8 :: total_states  ! for 18 sites
  integer*8 :: ii, i, k, P_tmp
  integer*8 :: l_a, l_start, l_end, l_rows, counter ! for 18 sites
  integer   :: f_start, f_end
  integer*8 :: bool_up, bool_dn
  integer*8 :: up_index, dn_index  ! for 18 sites
  integer :: ierr, rank, num_procs
  integer, allocatable :: rx(:), ry(:), pair_sites(:,:)
  double precision :: pi, qx, qy, sum_square, Cdw, Szz, Cdw_tmp, Szz_tmp
  double precision :: tmp, del_s, del_s_star, del_d_square, del_px, del_py, del_s_s_star, del_px_ipy, del_d_xy
  double precision :: del_s_tmp, del_s_star_tmp, del_d_square_tmp, del_px_tmp, del_py_tmp, del_s_s_star_tmp, del_px_ipy_tmp
  double precision :: del_d_xy_tmp
  double complex, allocatable :: coeff(:)
!  double precision, allocatable :: coeff(:)
  character :: file1*100
  integer, allocatable :: num_lines(:), num_files(:)
  complex*8 :: eye, sum_i

  ! For new states
  integer :: num_upspin_1, num_dnspin_1
  integer*8 :: total_states_1
  integer*8 :: num_states_up_1, num_states_dn_1
  integer, allocatable :: states_up_1(:), states_dn_1(:)
  double complex, allocatable :: coeff_1(:), del_tmp(:), del(:,:)
!  double precision, allocatable :: coeff_1(:), del_tmp(:), del(:,:)
  integer*8 :: key_up, key_dn, pos_up, pos_dn, new_index
  integer :: sign1, sign_count

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
       if(str1 == 'num_files_per_vec') num_files_per_vec = int(input_val)
       if(str1 == 'degeneracy') degeneracy = int(input_val)
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
  total_states = num_states_up*num_states_dn

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

  call MPI_Init ( ierr )
  call MPI_Comm_rank ( MPI_COMM_WORLD, rank, ierr ) ! Determine this process's rank.
  call MPI_Comm_size ( MPI_COMM_WORLD, num_procs, ierr ) ! Find out the number of processes available.

!!!  call MPI_Finalize ( ierr )

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

!  write(*,*) rank, f_start, f_end, l_start, l_end

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

!  allocate(coeff(0:l_end-l_start-1)); coeff = 0.0
  allocate(coeff(0:l_end-l_start-1)); coeff = (0.0, 0.0)
!  if(rank .eq. 0) allocate(del_tmp(0:total_states_1-1)); del_tmp = 0.0
  if(rank .eq. 0) allocate(del_tmp(0:total_states_1-1)); del_tmp = (0.0, 0.0)

  Cdw_tmp = 0.0
  Szz_tmp = 0.0
  del_s_tmp = 0.0
  del_s_star_tmp = 0.0
  del_d_square_tmp = 0.0
  del_px_tmp = 0.0
  del_py_tmp = 0.0
  del_s_s_star_tmp = 0.0
  del_px_ipy_tmp = 0.0
  del_d_xy = 0.0

  do P_tmp = 0, degeneracy-1
    counter = 0
    do i = f_start, f_end-1 
      write(unit=file1, fmt="('vec_',i2.2,'_',i5.5,'.txt')") P_tmp, i
      open(unit=111+rank,file=file1,status="old")
      do ii=0, num_lines(i)-1
        read(111+rank,*) coeff(counter)
        counter = counter+1
      end do
      close(111+rank)
    end do

    call MPI_AllReduce( dot_product(coeff,coeff), tmp, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
!    call MPI_AllReduce( dot_product(coeff,coeff), tmp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!    call MPI_Reduce( dot_product(coeff,coeff), tmp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    tmp = sqrt(tmp)
    if( rank .eq. 0 ) write(*, fmt="('vec ', i2, ' norm2 = ', f9.6)") P_tmp, tmp
    if( abs(tmp-1.0) .gt. small_val) stop


    ! For Cdw(q) (Szz)
    qx = pi
    qy = pi
    sum_square = 0.0
    do l_a=l_start, l_end-1 ! start from 0
      up_index = int(l_a/num_states_dn)
      dn_index = mod(l_a,num_states_dn)
      sum_i = (0.0, 0.0)
      do i=0, num_sites-1
        bool_up = btest(states_up(up_index),i)
        bool_dn = btest(states_dn(dn_index),i)
        sum_i = sum_i + exp(eye*(qx*rx(i)+qy*ry(i)))*(bool_up + bool_dn)
      end do
      sum_square = sum_square + (abs(coeff(l_a-l_start))*cabs(sum_i))**2
    end do
    call MPI_Reduce(sum_square, Cdw, 1, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    sum_square = 0.0
    do l_a=l_start, l_end-1 ! start from 0
      up_index = int(l_a/num_states_dn)
      dn_index = mod(l_a,num_states_dn)
      sum_i = (0.0, 0.0)
      do i=0, num_sites-1
        bool_up = btest(states_up(up_index),i)
        bool_dn = btest(states_dn(dn_index),i)
        sum_i = sum_i + exp(eye*(qx*rx(i)+qy*ry(i)))*(bool_up - bool_dn)
      end do
      sum_square = sum_square + (abs(coeff(l_a-l_start))*cabs(sum_i))**2
    end do
    call MPI_Reduce(sum_square, Szz, 1, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    if(rank .eq. 0) then
      Cdw = Cdw/num_sites
      Cdw_tmp = Cdw_tmp + Cdw
      Szz = Szz/4/num_sites
      Szz_tmp = Szz_tmp + Szz
      write(*,"('vec ',i4)") P_tmp 
      write(*,"('Cdw( pi, pi)      = ', f9.6)") Cdw
      write(*,"('Szz( pi, pi)      = ', f9.6)") Szz
    end if

    include 'pairing_sites_16.f90'
    do ii = 0, 4
!      coeff_1 = 0.0
      coeff_1 = (0.0, 0.0)
      do l_a=l_start, l_end-1 ! start from 0
        up_index = int(l_a/num_states_dn)
        dn_index = mod(l_a,num_states_dn)
        do i=0, num_sites-1
          ! For delta
          if ( (btest(states_up(up_index),i) .eqv. .true.) .and. (btest(states_dn(dn_index),pair_sites(i,ii)) .eqv. .true.) ) then
            key_up = states_up(up_index)
            key_dn = states_dn(dn_index)
            key_dn = ibclr(key_dn, pair_sites(i,ii))
            key_up = ibclr(key_up, i)
            pos_up = Binary_search(states_up_1, num_states_up_1, key_up)
            pos_dn = Binary_search(states_dn_1, num_states_dn_1, key_dn)
            sign_count = num_upspin
            do k = num_sites-1, pair_sites(i,ii)+1, -1
              if (btest(states_dn(dn_index),k) .eqv. .true.) sign_count = sign_count+1
            end do
            do k = num_sites-1, i+1, -1
              if (btest(states_up(up_index),k) .eqv. .true.) sign_count = sign_count+1
            end do
            sign1 = 1
            if (mod(sign_count,2) .eq. 1) sign1 = -1
            new_index = pos_up*num_states_dn_1+pos_dn
            coeff_1(new_index) = coeff_1(new_index) + sign1*coeff(l_a-l_start)
          end if
        end do
      end do
!      call MPI_Reduce(coeff_1, del(ii,:), total_states_1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(coeff_1, del(ii,:), total_states_1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    end do

!    if(rank .eq. 0) del_tmp = 0.0
    if(rank .eq. 0) del_tmp = (0.0, 0.0)

    if(rank .eq. 0) then

      del_tmp = del(0,:)
      del_s =  dot_product(del_tmp, del_tmp)/num_sites
      del_s_tmp = del_s_tmp + del_s
      del_tmp = (del(1,:) + del(2,:) + del(3,:) + del(4,:))/2
      del_s_star =  dot_product(del_tmp, del_tmp)/num_sites
      del_s_star_tmp = del_s_star_tmp + del_s_star
      del_tmp = (del(1,:) - del(3,:))/sqrt(2.0)
      del_px =  dot_product(del_tmp, del_tmp)/num_sites
      del_px_tmp = del_px_tmp + del_px
      del_tmp = (del(2,:) - del(4,:))/sqrt(2.0)
      del_py =  dot_product(del_tmp, del_tmp)/num_sites
      del_py_tmp = del_py_tmp + del_py
      del_tmp = (del(1,:) + del(3,:) - del(2,:) - del(4,:))/2
      del_d_square =  dot_product(del_tmp, del_tmp)/num_sites
      del_d_square_tmp = del_d_square_tmp + del_d_square


!      del_tmp = del(0,:) + (del(1,:) + del(2,:) + del(3,:) + del(4,:))/2
!      del_s_s_star =  dot_product(del_tmp, del_tmp)/num_sites
!      del_s_s_star_tmp = del_s_s_star_tmp + del_s_s_star

      del_tmp = (del(1,:) - del(3,:) + eye*(del(2,:) - del(4,:) ))/2
      del_px_ipy =  dot_product(del_tmp, del_tmp)/num_sites
      del_px_ipy_tmp = del_px_ipy_tmp + del_px_ipy


!      write(*,"('vec ',i4)") P_tmp 
      write(*,"('Delta_s           = ', f9.6)") del_s
      write(*,"('Delta_s_star      = ', f9.6)") del_s_star
      write(*,"('Delta_px          = ', f9.6)") del_px
      write(*,"('Delta_py          = ', f9.6)") del_py
      write(*,"('Delta_d_(x^2-y^2) = ', f9.6)") del_d_square
!      write(*,"('Delta_s+s_star    = ', f9.6)") del_s_s_star
      write(*,"('Delta_px+ipy      = ', f9.6)") del_px_ipy

    end if

  call MPI_Barrier(MPI_COMM_WORLD, ierr)

    pair_sites = 0
    include 'pairing_sites_16_dxy.f90'

!    del = 0.0
    del = (0.0, 0.0)

    do ii = 1, 4
!      coeff_1 = 0.0
      coeff_1 = (0.0, 0.0)
      do l_a=l_start, l_end-1 ! start from 0
        up_index = int(l_a/num_states_dn)
        dn_index = mod(l_a,num_states_dn)
        do i=0, num_sites-1
          ! For delta
          if ( (btest(states_up(up_index),i) .eqv. .true.) .and. (btest(states_dn(dn_index),pair_sites(i,ii)) .eqv. .true.) ) then
            key_up = states_up(up_index)
            key_dn = states_dn(dn_index)
            key_dn = ibclr(key_dn, pair_sites(i,ii))
            key_up = ibclr(key_up, i)
            pos_up = Binary_search(states_up_1, num_states_up_1, key_up)
            pos_dn = Binary_search(states_dn_1, num_states_dn_1, key_dn)
            sign_count = num_upspin
            do k = num_sites-1, pair_sites(i,ii)+1, -1
              if (btest(states_dn(dn_index),k) .eqv. .true.) sign_count = sign_count+1
            end do
            do k = num_sites-1, i+1, -1
              if (btest(states_up(up_index),k) .eqv. .true.) sign_count = sign_count+1
            end do
            sign1 = 1
            if (mod(sign_count,2) .eq. 1) sign1 = -1
            new_index = pos_up*num_states_dn_1+pos_dn
            coeff_1(new_index) = coeff_1(new_index) + sign1*coeff(l_a-l_start)
          end if
        end do
      end do
!      call MPI_Reduce(coeff_1, del(ii,:), total_states_1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(coeff_1, del(ii,:), total_states_1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    end do


!    if(rank .eq. 0) del_tmp = 0.0
    if(rank .eq. 0) del_tmp = (0.0, 0.0)

    if(rank .eq. 0) then

      del_tmp = (del(1,:) + del(3,:) - del(2,:) - del(4,:))/2
      del_d_xy =  dot_product(del_tmp, del_tmp)/num_sites
      del_d_xy_tmp = del_d_xy_tmp + del_d_xy

      write(*,"('Delta_d_xy        = ', f9.6)") del_d_xy
      write(*,*)

    end if
  end do

  if(rank .eq. 0) then
    write(*,*)
    write(*,"('Average Cdw( pi, pi)      = ', f9.6)") Cdw_tmp / degeneracy
    write(*,"('Average Szz( pi, pi)      = ', f9.6)") Szz_tmp / degeneracy
    write(*,"('Average Delta_s           = ', f9.6)") del_s_tmp/degeneracy
    write(*,"('Average Delta_s_star      = ', f9.6)") del_s_star_tmp/degeneracy
    write(*,"('Average Delta_d_square    = ', f9.6)") del_d_square_tmp/degeneracy
    write(*,"('Average Delta_px          = ', f9.6)") del_px_tmp/degeneracy
    write(*,"('Average Delta_py          = ', f9.6)") del_py_tmp/degeneracy
!    write(*,"('Average Delta_s_s_star    = ', f9.6)") del_s_s_star_tmp/degeneracy
    write(*,"('Average Delta_px_ipy      = ', f9.6)") del_px_ipy_tmp/degeneracy
    write(*,"('Average Delta_d_xy        = ', f9.6)") del_d_xy_tmp/degeneracy
    write(*,*)
  end if


  call MPI_Finalize (ierr)

  include 'contains.f90'

end program measurements

