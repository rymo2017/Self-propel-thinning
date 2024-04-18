module mo_config
    !
    !  base struct and method of configuration generating, storing...
    !
    use mo_syst
    use mo_math
    implicit none

    type tpcon

        integer :: natom
        real(8), allocatable, dimension(:) :: rx,ry,vx,vy,fx,fy
        real(8), allocatable, dimension(:) :: sfx,sfy
        real(8), allocatable, dimension(:) ::  r, theta, deta_theta
        real(8) :: la(free)
        real(8) :: phi
        real(8) :: Ea

        real(8) :: strain,stress
    
    end type

    type(tpcon) :: con, con0, contemp


contains

    subroutine init_system( tcon, tnatom, tphi )
        implicit none
            
        type(tpcon), intent(inout) :: tcon
        integer,     intent(in)    :: tnatom
        real(8), optional, intent(in) ::tphi

        tcon%natom = tnatom
        if(present(tphi)) tcon%phi = tphi

        allocate( tcon%theta(tnatom) )
        allocate( tcon%deta_theta(tnatom) )
        allocate( tcon%rx(tnatom) )
        allocate( tcon%ry(tnatom) )
        allocate( tcon%vx(tnatom) )
        allocate( tcon%vy(tnatom) )
        allocate( tcon%fx(tnatom) )
        allocate( tcon%fy(tnatom) )
        allocate( tcon%r(tnatom) )
        allocate( tcon%sfx(tnatom) )
        allocate( tcon%sfy(tnatom) )

    end subroutine


    subroutine gen_rand_config( tcon, tseed, tabsf, tphi )
        implicit none
            
        ! para list
        type(tpcon), intent(inout)        :: tcon
        integer,     intent(inout)        :: tseed
        real(8),      intent(in)          :: tabsf
        real(8),     intent(in), optional :: tphi
    

        ! local
        integer :: i, j, k
        real(8) :: sumx,sumy
        real(8), allocatable, dimension(:,:) :: locsfa

        ! initialized rand
        call init_rand(tseed)
        tseed = 0

        if ( present( tphi ) ) tcon%phi = tphi
        
        allocate( locsfa(free,sets%natom) )

        associate(                &
            natom  => tcon%natom, &
            rx     => tcon%rx,    &
            ry     => tcon%ry,    &
            fx     => tcon%fx,    &
            sfx    => tcon%sfx,   &
            fy     => tcon%fy,    &
            sfy    => tcon%sfy,   &
            vx     => tcon%vx,    &
            vy     => tcon%vy,    &
            r      => tcon%r,     &
            theta => tcon%theta,&
            deta_theta => tcon%deta_theta,&
            la     => tcon%la,    &
            strain => tcon%strain &
            )

            ! radius
            r(1:natom/2)       = 0.5d0
            r(natom/2+1:natom) = 0.5d0 * ratio

            ! box length
            la     = calc_box_length(tcon)
            strain = 0.d0

            ! config
            call random_number(rx)
            call random_number(ry)
            do i=1, natom
                rx(i) = ( rx(i) - 0.5d0 ) * la(1)
                ry(i) = ( ry(i) - 0.5d0 ) * la(2)
            end do


            call random_number( theta )
            do i=1,natom
                theta(i) = 2.d0*pi*theta(i)
                locsfa(1,i) = tabsf * dcos(theta(i))
                locsfa(2,i) = tabsf * dsin(theta(i))
            enddo
            do i=1,natom
                sfx(i) = locsfa(1,i)
                sfy(i) = locsfa(2,i)
            enddo

            
        
            ! f v
            deta_theta = 0.d0
            vx = 0.d0
            vy = 0.d0
            fx = 0.d0
            fy = 0.d0
        

        end associate
    end subroutine

    subroutine trim_config( tcon )
        implicit none
            
        type(tpcon), intent(inout) :: tcon

        integer :: i,cory

        associate(          &
            natom => tcon%natom,&
            rx => tcon%rx,&
            ry => tcon%ry,&
            la => tcon%la,&
            strain => tcon%strain&
            )
            
            do i=1,natom
                cory = nint( ry(i)/la(2) )
                rx(i) = rx(i) - strain * cory * la(2)
                rx(i) = rx(i) - nint( rx(i)/la(1) ) * la(1)
                ry(i) = ry(i) - cory * la(2)
            enddo
        end associate
    end subroutine


            



        
     pure function calc_box_length(tcon) result(l)
        implicit none

        ! para list
        type(tpcon), intent(in) :: tcon

        ! result
        real(8) :: l

        ! local
        real(8) :: phi
        real(8) :: sdisk, volume

        phi = tcon%phi

        ! V_n(r) = pi^(n/2) / Gamma( n/2 + 1 ) * r^n
        sdisk = sqrt(pi**free) / gamma(dble(free)/2.d0+1) * sum(tcon%r**free)

        ! box length
        volume = sdisk / phi
        l      = volume ** ( 1.d0 / dble(free) )
    end function

end module

