
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
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

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
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    if(rank .eq. 0) then
      Cdw = Cdw/num_sites
      Szz = Szz/4/num_sites
!      write(*,"('vec ',i4)") P_tmp 
!      write(*,"('Cdw( pi, pi)      = ', f9.6)") Cdw
!      write(*,"('Szz( pi, pi)      = ', f9.6)") Szz
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
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end do

!    if(rank .eq. 0) del_tmp = 0.0
    if(rank .eq. 0) del_tmp = (0.0, 0.0)

    if(rank .eq. 0) then

      del_tmp = del(0,:)
      del_s =  dot_product(del_tmp, del_tmp)/num_sites
      del_tmp = (del(1,:) + del(2,:) + del(3,:) + del(4,:))/2
      del_s_star =  dot_product(del_tmp, del_tmp)/num_sites

      del_tmp = (del(1,:) + del(3,:) - del(2,:) - del(4,:))/2
      del_d_square =  dot_product(del_tmp, del_tmp)/num_sites

      del_tmp = (del(1,:) - del(3,:))/sqrt(2.0)
      del_px =  dot_product(del_tmp, del_tmp)/num_sites
      del_tmp = (del(2,:) - del(4,:))/sqrt(2.0)
      del_py =  dot_product(del_tmp, del_tmp)/num_sites


!      del_tmp = del(0,:) + (del(1,:) + del(2,:) + del(3,:) + del(4,:))/4
!      del_s_s_star =  dot_product(del_tmp, del_tmp)/num_sites

      del_tmp = (del(1,:) - del(3,:) + (del(2,:) - del(4,:) ))/2
      del_px_py =  dot_product(del_tmp, del_tmp)/num_sites

      del_tmp = (del(1,:) - del(3,:) - (del(2,:) - del(4,:) ))/2
      del_px_m_py =  dot_product(del_tmp, del_tmp)/num_sites

      del_tmp = (del(1,:) - del(3,:) + eye*(del(2,:) - del(4,:) ))/2
      del_px_ipy =  dot_product(del_tmp, del_tmp)/num_sites

      del_tmp = (del(1,:) - del(3,:) - eye*(del(2,:) - del(4,:) ))/2
      del_px_m_ipy =  dot_product(del_tmp, del_tmp)/num_sites


!      write(*,"('vec ',i4)") P_tmp 
!      write(*,"('Delta_s           = ', f9.6)") del_s
!      write(*,"('Delta_s_star      = ', f9.6)") del_s_star
!      write(*,"('Delta_px          = ', f9.6)") del_px
!      write(*,"('Delta_py          = ', f9.6)") del_py
!      write(*,"('Delta_d_(x^2-y^2) = ', f9.6)") del_d_square
!      write(*,"('Delta_s+s_star    = ', f9.6)") del_s_s_star
!      write(*,"('Delta_px+ipy      = ', f9.6)") del_px_ipy

    end if

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

