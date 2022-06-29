
contains
  recursive function Combination(m,n) result(a)
    integer*8 :: a
    integer, intent(in) :: m,n

    if (n == 1) then
      a = m
    else
      a = nint(real(m+1-n) / n * Combination(m,n-1))
    end if

  end function Combination

  recursive subroutine Gen(m,m_max)

    implicit none
    integer, intent (in) :: m
    integer, intent (in) :: m_max
    integer :: n
    integer :: i
    integer :: num

    if (m > m_max) then
      num = 0
      do i= 0, m_max
        num = num + 2**comb(i)
      end do
      states_tmp(num_states_tmp-1-ii) = num
      ii = ii+1
    else
      do n = num_sites-1, 0, -1
        if ((m == 0) .or. (n < comb (m - 1))) then
          comb (m) = n
          call gen (m + 1,m_max)
        end if
      end do
    end if

  end subroutine Gen


  integer function Binary_search(array, n, key)

    implicit none
    integer*8, intent (in) :: n, key
    integer, intent (in) :: array(0:n-1)
    integer :: l, r, m  ! left, right, and middle

    l = 0
    r = n-1
    m = (l+r)/2
    if ( (key < array(l) ) .or. (key > array(r)) ) then
      Binary_search = -1
      return
    end if
    
    do while ( l <= r )
      if ( key > array(m) ) then
        l = m+1
        m = (l+r)/2
      else if (key < array(m)) then
        r = m-1
        m = (l+r)/2
      else if (key == array(m)) then
        Binary_search = m
        return
      end if
    end do

  Binary_search = -1
  return

  end function


