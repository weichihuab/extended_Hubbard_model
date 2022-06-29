
  if (num_sites /= 16) then
   write(*,*) 'Number of sites is not 16.'
   write(*,*) 'Stop.'
   stop
  end if

  T_matrix( 0, 1) = -t1 
  T_matrix( 0, 4) = -t1
  T_matrix( 0, 3) = -t1
  T_matrix( 0,12) = -t1

  T_matrix( 1, 2) = -t1 
  T_matrix( 1, 5) = -t1
  T_matrix( 1, 0) = -t1
  T_matrix( 1,13) = -t1

  T_matrix( 2, 3) = -t1 
  T_matrix( 2, 6) = -t1
  T_matrix( 2, 1) = -t1
  T_matrix( 2,14) = -t1

  T_matrix( 3, 0) = -t1 
  T_matrix( 3, 7) = -t1
  T_matrix( 3, 2) = -t1
  T_matrix( 3,15) = -t1

  T_matrix( 4, 5) = -t1 
  T_matrix( 4, 8) = -t1
  T_matrix( 4, 7) = -t1
  T_matrix( 4, 0) = -t1

  T_matrix( 5, 6) = -t1 
  T_matrix( 5, 9) = -t1
  T_matrix( 5, 4) = -t1
  T_matrix( 5, 1) = -t1

  T_matrix( 6, 7) = -t1 
  T_matrix( 6,10) = -t1
  T_matrix( 6, 5) = -t1
  T_matrix( 6, 2) = -t1

  T_matrix( 7, 4) = -t1 
  T_matrix( 7,11) = -t1
  T_matrix( 7, 6) = -t1
  T_matrix( 7, 3) = -t1

  T_matrix( 8, 9) = -t1 
  T_matrix( 8,12) = -t1
  T_matrix( 8,11) = -t1
  T_matrix( 8, 4) = -t1

  T_matrix( 9,10) = -t1 
  T_matrix( 9,13) = -t1
  T_matrix( 9, 8) = -t1
  T_matrix( 9, 5) = -t1

  T_matrix(10,11) = -t1 
  T_matrix(10,14) = -t1
  T_matrix(10, 9) = -t1
  T_matrix(10, 6) = -t1

  T_matrix(11, 8) = -t1 
  T_matrix(11,15) = -t1
  T_matrix(11,10) = -t1
  T_matrix(11, 7) = -t1

  T_matrix(12,13) = -t1 
  T_matrix(12, 0) = -t1
  T_matrix(12,15) = -t1
  T_matrix(12, 8) = -t1

  T_matrix(13,14) = -t1 
  T_matrix(13, 1) = -t1
  T_matrix(13,12) = -t1
  T_matrix(13, 9) = -t1

  T_matrix(14,15) = -t1 
  T_matrix(14, 2) = -t1
  T_matrix(14,13) = -t1
  T_matrix(14,10) = -t1

  T_matrix(15,12) = -t1 
  T_matrix(15, 3) = -t1
  T_matrix(15,14) = -t1
  T_matrix(15,11) = -t1

!---------------------------------!

  V_matrix( 0, 1) = v1
  V_matrix( 0, 4) = v1
  V_matrix( 0, 3) = v1
  V_matrix( 0,12) = v1

  V_matrix( 1, 2) = v1 
  V_matrix( 1, 5) = v1
  V_matrix( 1, 0) = v1
  V_matrix( 1,13) = v1

  V_matrix( 2, 3) = v1 
  V_matrix( 2, 6) = v1
  V_matrix( 2, 1) = v1
  V_matrix( 2,14) = v1

  V_matrix( 3, 0) = v1 
  V_matrix( 3, 7) = v1
  V_matrix( 3, 2) = v1
  V_matrix( 3,15) = v1

  V_matrix( 4, 5) = v1 
  V_matrix( 4, 8) = v1
  V_matrix( 4, 7) = v1
  V_matrix( 4, 0) = v1

  V_matrix( 5, 6) = v1 
  V_matrix( 5, 9) = v1
  V_matrix( 5, 4) = v1
  V_matrix( 5, 1) = v1

  V_matrix( 6, 7) = v1 
  V_matrix( 6,10) = v1
  V_matrix( 6, 5) = v1
  V_matrix( 6, 2) = v1

  V_matrix( 7, 4) = v1 
  V_matrix( 7,11) = v1
  V_matrix( 7, 6) = v1
  V_matrix( 7, 3) = v1

  V_matrix( 8, 9) = v1 
  V_matrix( 8,12) = v1
  V_matrix( 8,11) = v1
  V_matrix( 8, 4) = v1

  V_matrix( 9,10) = v1 
  V_matrix( 9,13) = v1
  V_matrix( 9, 8) = v1
  V_matrix( 9, 5) = v1

  V_matrix(10,11) = v1 
  V_matrix(10,14) = v1
  V_matrix(10, 9) = v1
  V_matrix(10, 6) = v1

  V_matrix(11, 8) = v1 
  V_matrix(11,15) = v1
  V_matrix(11,10) = v1
  V_matrix(11, 7) = v1

  V_matrix(12,13) = v1 
  V_matrix(12, 0) = v1
  V_matrix(12,15) = v1
  V_matrix(12, 8) = v1

  V_matrix(13,14) = v1 
  V_matrix(13, 1) = v1
  V_matrix(13,12) = v1
  V_matrix(13, 9) = v1

  V_matrix(14,15) = v1 
  V_matrix(14, 2) = v1
  V_matrix(14,13) = v1
  V_matrix(14,10) = v1

  V_matrix(15,12) = v1 
  V_matrix(15, 3) = v1
  V_matrix(15,14) = v1
  V_matrix(15,11) = v1

!----------------------------

  T_matrix( 0, 5) = -t2 
  T_matrix( 0, 7) = -t2
  T_matrix( 0,15) = -t2
  T_matrix( 0,13) = -t2

  T_matrix( 1, 6) = -t2 
  T_matrix( 1, 4) = -t2
  T_matrix( 1,12) = -t2
  T_matrix( 1,14) = -t2

  T_matrix( 2, 7) = -t2 
  T_matrix( 2, 5) = -t2
  T_matrix( 2,13) = -t2
  T_matrix( 2,15) = -t2

  T_matrix( 3, 4) = -t2 
  T_matrix( 3, 6) = -t2
  T_matrix( 3,14) = -t2
  T_matrix( 3,12) = -t2

  T_matrix( 4, 9) = -t2 
  T_matrix( 4,11) = -t2
  T_matrix( 4, 3) = -t2
  T_matrix( 4, 1) = -t2

  T_matrix( 5,10) = -t2 
  T_matrix( 5, 8) = -t2
  T_matrix( 5, 0) = -t2
  T_matrix( 5, 2) = -t2

  T_matrix( 6,11) = -t2 
  T_matrix( 6, 9) = -t2
  T_matrix( 6, 1) = -t2
  T_matrix( 6, 3) = -t2

  T_matrix( 7, 8) = -t2 
  T_matrix( 7,10) = -t2
  T_matrix( 7, 2) = -t2
  T_matrix( 7, 0) = -t2

  T_matrix( 8,13) = -t2 
  T_matrix( 8,15) = -t2
  T_matrix( 8, 7) = -t2
  T_matrix( 8, 5) = -t2

  T_matrix( 9,14) = -t2 
  T_matrix( 9,12) = -t2
  T_matrix( 9, 4) = -t2
  T_matrix( 9, 6) = -t2

  T_matrix(10,15) = -t2 
  T_matrix(10,13) = -t2
  T_matrix(10, 5) = -t2
  T_matrix(10, 7) = -t2

  T_matrix(11,12) = -t2 
  T_matrix(11,14) = -t2
  T_matrix(11, 6) = -t2
  T_matrix(11, 4) = -t2

  T_matrix(12, 1) = -t2 
  T_matrix(12, 3) = -t2
  T_matrix(12,11) = -t2
  T_matrix(12, 9) = -t2

  T_matrix(13, 2) = -t2 
  T_matrix(13, 0) = -t2
  T_matrix(13, 8) = -t2
  T_matrix(13,10) = -t2

  T_matrix(14, 3) = -t2 
  T_matrix(14, 1) = -t2
  T_matrix(14, 9) = -t2
  T_matrix(14,11) = -t2

  T_matrix(15, 0) = -t2 
  T_matrix(15, 2) = -t2
  T_matrix(15,10) = -t2
  T_matrix(15, 8) = -t2

  do ii=0, num_sites-1
     do jj=0, num_sites-1
        if((abs(T_matrix(ii,jj)-T_matrix(jj,ii)) .gt. 1E-8) .or. (abs(V_matrix(ii,jj)-V_matrix(jj,ii)) .gt. 1E-8)) then
          write(*,*) "T_matrix or V_matrix is not Hermitian"
          stop
        endif
     enddo
  enddo

  if (rank .eq. 0) write(*,*) "Pass Hermitian test" 

