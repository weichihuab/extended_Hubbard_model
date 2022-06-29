
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
          Ax = A0*exp(-(t-t0)**2/(2*sigma**2))*cos(omega*(t-t0))*pol_x/sqrt(pol_x**2+pol_y**2) 
          Ay = A0*exp(-(t-t0)**2/(2*sigma**2))*cos(omega*(t-t0))*pol_y/sqrt(pol_x**2+pol_y**2) 
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
          Ax = A0*exp(-(t-t0)**2/(2*sigma**2))*cos(omega*(t-t0))*pol_x/sqrt(pol_x**2+pol_y**2) 
          Ay = A0*exp(-(t-t0)**2/(2*sigma**2))*cos(omega*(t-t0))*pol_y/sqrt(pol_x**2+pol_y**2) 
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
