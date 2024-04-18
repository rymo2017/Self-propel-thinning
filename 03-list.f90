module mo_list
    !
    !  verlet list, see more at page 147 of <computer simulation of liquids>
    !
    use mo_syst
    use mo_config
    implicit none

    ! global constants.
    !! skin
    real(8), private, parameter :: set_nlcut = 0.70d0
    !! for sake of memory saving, we consider listmax neighbor of one particle at most
    !! enlarge this if you study 3D system or high volume fraction system
    integer, private, parameter :: listmax = 64

    ! neighbor list of one particle
    type tplistone
        ! neighbor number
        integer    :: nbsum
        ! neighbor index
        integer    :: nblist(listmax)
        ! precalculated relative relation, used for distance calculation in perodical cell
        !integer(1) :: iround(free,listmax)
        !integer(1) :: cory(listmax)
        !
!        real(8)    :: con0(free)
        real(8)     :: rx0,ry0
    end type

    ! neighbor list struct of system
    type tplist
        integer                                    :: natom
        type(tplistone), allocatable, dimension(:) :: list
        ! contact number
        integer, allocatable, dimension(:)         :: nbi
        ! tag particle is rattler or not
        integer, allocatable, dimension(:)         :: rattlerflag
        real(8)    :: nlcut
    end type

    type(tplist) :: nb


contains

    subroutine init_list( tnb,tcon )
        implicit none

        ! para list
        type(tplist), intent(inout) :: tnb
        type(tpcon),  intent(in)    :: tcon

        ! local
        integer :: tnatom

        tnatom    = tcon%natom
        tnb%natom = tnatom
        tnb%nlcut = set_nlcut

        if ( allocated(tnb%list) ) then
            if ( size(tnb%list) /= tnatom ) then
                deallocate( tnb%list )
                allocate( tnb%list(tnatom) )
            end if
        else
            allocate( tnb%list(tnatom) )
        end if
    end subroutine

    subroutine make_list( tnb,tcon )
        implicit none
         
        type(tpcon), intent(in) :: tcon
        type(tplist), intent(inout) :: tnb
   
        real(8) :: lainv(free)
        real(8) :: rxij, ryij, dij, rij
        integer :: i, j ,k
        integer :: cory, itemp
   
        associate(          &
           natom => tcon%natom, &
           rx    => tcon%rx,&
           ry    => tcon%ry,&
           r     => tcon%r,&
           la    => tcon%la,&
           strain => tcon%strain,&
           list   => tnb%list, &
           nlcut => tnb%nlcut &
           )
   
           nlcut = set_nlcut
           lainv = 1.d0/la
   
           list(:)%nbsum = 0
   
           do i=1,natom
               list(i)%rx0 = rx(i)
               list(i)%ry0 = ry(i)
   
               do j=i+1,natom
                   rxij = rx(i) - rx(j)
                   ryij = ry(i) - ry(j)
                   cory = nint( ryij * lainv(2) )
                   rxij = rxij - cory * strain * la(2)
                   rxij = rxij - nint( rxij * lainv(1) ) * la(1)
                   ryij = ryij - cory * la(2)
                   rij = dsqrt( rxij**2 + ryij**2 )
                   dij = r(i) + r(j)
                   if( rij > (dij+nlcut) ) cycle
                   if( list(i)%nbsum < listmax ) then
                       list(i)%nbsum = list(i)%nbsum + 1
                       list(i)%nblist(list(i)%nbsum) = j
                   endif
               enddo
           enddo
       end associate
end subroutine




    function check_list( tnb, tcon ) result(flag)
        !
        !  determine remake list or not
        !
        implicit none

        ! para list
        type(tplist), intent(in) :: tnb
        type(tpcon),  intent(in) :: tcon

        ! result
        logical :: flag

        ! local
        real(8) :: maxdis, drax,dray, dr2
        integer :: i

        associate(               &
            natom => tcon%natom, &
            rx    => tcon%rx,    &
            ry    => tcon%ry,     &
            nlcut => tnb%nlcut   &
            )

            maxdis = 0.d0
            do i=1, tcon%natom
                drax = tcon%rx(i) - tnb%list(i)%rx0
                dray = tcon%ry(i) - tnb%list(i)%ry0
                dr2 = drax**2 + dray**2
                if ( maxdis < dr2 ) maxdis = dr2
            end do

        flag = .false.
        if ( maxdis > 0.25 * nlcut**2 ) flag = .true.

        end associate
    end function


end module



