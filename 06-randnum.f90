Module ran_mod
    use mo_math
    use mo_config
  Implicit None
! ran return a uniform random number between 0-1  
! norma return a normal distribution  
contains 

    
 subroutine init_random_seed()
    use iso_fortran_env, only: int64
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid
    integer(int64) :: t
          
    call random_seed(size = n)
    allocate(seed(n))
            ! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/urandom", access="stream", &
        form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
        read(un) seed
       close(un)
    else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
        call system_clock(t)
        if (t == 0) then
            call date_and_time(values=dt)
            t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24_int64 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
        end if
        pid = getpid()
        t = ieor(t, int(pid, kind(t)))
        do i = 1, n
            seed(i) = lcg(t)
        end do
    end if
    call random_seed(put=seed)
    contains
            ! This simple PRNG might not be good enough for real work, but is
            ! sufficient for seeding a better PRNG.
        function lcg(s)
            integer :: lcg
            integer(int64) :: s
            if (s == 0) then
                s = 104729
            else
                s = mod(s, 4294967296_int64)
            end if
            s = mod(s * 279470273_int64, 4294967291_int64)
            lcg = int(mod(s, int(huge(0), int64)), kind(0))
        end function lcg
    end subroutine init_random_seed
   
    
    
    function ran()   !returns random number between 0 - 1  
    implicit none 
    integer , save :: flag = 0
    double precision :: ran 
    if(flag==0) then 
      call init_random_seed()
      flag = 1 
    endif
    call random_number(ran)     ! built in fortran 90 random number function  
  end function ran
  
  function normal(mean,sigma) 
    implicit none 
    integer :: flag 
    double precision, parameter :: pi = 3.141592653589793239  
    double precision :: u1, u2, y1, y2, normal, mean, sigma 
    save flag 
    data flag /0/ 
 
    u1 = ran(); u2 = ran() 
    do while(u1==0.d0)
        u1 = ran()
    enddo
    do while(u2==0.d0)
        u2 = ran()
    enddo
   
    if (flag.eq.0) then
      y1 = sqrt(-2.0d0*log(u1))*cos(2.0d0*pi*u2) 
      normal = mean + sigma*y1 
      flag = 1 
    else 
      y2 = sqrt(-2.0d0*log(u1))*sin(2.0d0*pi*u2) 
      normal = mean + sigma*y2 
      flag = 0 
    endif  
  end function normal 
End Module ran_mod

