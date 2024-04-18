module mo_force
    use mo_syst
    use mo_config
    use mo_list
    use ran_mod
!    use mo_network
    implicit none

contains
    subroutine calc_deta_theta( tcon,rotate_mean,rotate_v,tabsf,dt )
        implicit none
            
        type(tpcon), intent(inout) :: tcon
        real(8), intent(in) :: rotate_mean
        real(8), intent(in) :: rotate_v
        real(8), intent(in) :: tabsf
        real(8), intent(in) :: dt
        
        real(8) :: double_rotatev
        
    

        integer :: i

        double_rotatev = dsqrt(2.d0*rotate_v/dt)

        do i=1,tcon%natom
            tcon%deta_theta(i) = normal( rotate_mean,double_rotatev )
        enddo
        do i=1,tcon%natom
            tcon%theta(i) = tcon%theta(i) + tcon%deta_theta(i)*dt
            tcon%sfx(i) = tabsf*dcos( tcon%theta(i) )
            tcon%sfy(i) = tabsf*dsin( tcon%theta(i) )
        enddo


    end subroutine



     subroutine calc_attra_force_vis( tcon,tnb,tb,shearrate )
        implicit none

        type(tpcon),  intent(inout) :: tcon
        type(tplist), intent(in) :: tnb
        real(8), intent(in) :: tb
        real(8), intent(in) :: shearrate


        real(8) :: rxij, ryij, dij, rij, fr, fxij, fyij, dijcut, pxij, pyij,pij
        integer :: i, j, k ,jj, cory

        associate(                     &
            natom    => tcon%natom,    &
            r        => tcon%r,        &
            rx       => tcon%rx,       &
            ry       => tcon%ry,       &
            fx       => tcon%fx,       &
            fy       => tcon%fy,       &
            sfx      => tcon%sfx,      &
            sfy      => tcon%sfy,      &
            vx       => tcon%vx,       &
            vy       => tcon%vy,       &
            la       => tcon%la,       &
            Ea       => tcon%Ea,       &
            strain   => tcon%strain,   &
            stress   => tcon%stress,   &
            list     => tnb%list       &
            )

        do i=1,natom
            fx(i) = 0.d0
            fy(i) = 0.d0
        enddo
       
        stress = 0.d0

    


        do i=1,natom
            do jj=1,list(i)%nbsum
                j = list(i)%nblist(jj)
                rxij = rx(i) - rx(j)
                ryij = ry(i) - ry(j)
                cory = nint( ryij/la(2) )
                rxij = rxij - cory * strain * la(2)
                rxij = rxij - nint(rxij/la(1))*la(1)
                ryij = ryij - cory * la(2)
                rij = dsqrt( rxij**2 + ryij**2 )
                dij = r(i) + r(j)
                dijcut = dij * ( 1.d0 + 2.d0*attra )
        

                if( rij > dijcut ) cycle

                if( rij .le. dij*(1.d0+attra) ) then
                    fr = 2.d0 * (1.d0 - rij/dij)**(alpha-1) / rij / dij
                else
                    fr = - 2.d0 * ( 1.d0 +2.d0*attra - rij/dij )**(alpha-1)/rij/dij
                endif

                pxij = vx(i)  - vx(j)
                pyij = vy(i)  - vy(j)
                pxij = pxij - cory*shearrate*la(2)
                

                fxij = fr * rxij - tb*pxij
                fyij = fr * ryij - tb*pyij
                

                fx(i) = fx(i) + fxij
                fx(j) = fx(j) - fxij
                fy(i) = fy(i) + fyij
                fy(j) = fy(j) - fyij

                stress = stress - ryij*fxij - rxij*fyij
                
            enddo
        enddo

        stress = stress / product(la) / free

        end associate
    end subroutine


        
        
end module            
                






    
