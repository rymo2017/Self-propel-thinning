program main
    use mo_syst
    use mo_var
    use mo_config
    use mo_list
    use mo_force

    implicit none


    real(8) :: shearrate, dt, dt2,dt22,tb,tabsf,rotate_mean,rotate_v
    real(8) :: strain0,strain1,sumstress,stress,stress0,stress1,sumstress0,sumstress1,stressxy,sumstressxy
    integer :: nframe, iprint,stepmax,seed,nfile
    character(30) ::tempp,filetemp,initfile,finalconfig,temppp,filetempp
    type(tpcon) :: confire



    !call readvar
    call testvar
    sets%natom = 2048
    sets%phi = 0.75d0
    rotate_mean = 0.d0
    rotate_v = 1d-3

    nframe = 0
    nfile = 0
    sumstress = 0.d0
    sumstress0 = 0.d0
    sumstress1 = 0.d0
    sumstressxy = 0.d0
    dt2 = 0.5d0 * dt
    dt22 = 0.5d0 * dt**2
    seed = sets%seed

    write( tempp,'(i6)' ) seed
    filetemp = 'strain-stress'//trim(adjustl(tempp))//'.txt'
    finalconfig = 'finalcon'//trim(adjustl(tempp))//'.txt'
    initfile = 'con'//trim(adjustl(tempp))//'.txt'


    open(22,file=filetemp)


    call init_system( con,sets%natom )
    call gen_rand_config( con,sets%seed,tabsf,sets%phi )


    open(1,file=initfile)
    do i=1,sets%natom
        read(1,'(3es26.16)') con%rx(i),con%ry(i),con%r(i)
    enddo
    close(1)

    call init_list(nb,con)
    call make_list(nb,con)

    call calc_attra_force_vis(con,nb,tb,shearrate)


    do while( con%strain .lt. strain1 )
        step = step + 1
        con%strain = con%strain + shearrate * dt

        

        call calc_deta_theta( con,rotate_mean,rotate_v,tabsf,dt ) 
        
        con%rx = con%rx + con%vx * dt + (con%fx + con%sfx) * dt22
        con%ry = con%ry + con%vy * dt + (con%fy + con%sfy) * dt22
        con%vx = con%vx + (con%fx + con%sfx) * dt2
        con%vy = con%vy + (con%fy + con%sfy) * dt2
      
        !call predic( con,dt )
        
        if( check_list( nb, con ) ) then
            call make_list( nb, con )
        end if
                
        call calc_attra_force_vis( con, nb,tb,shearrate )


        confire = con
        call trim_config(confire)
        stressxy = 0.d0
        do i=1,sets%natom
            stressxy = stressxy - confire%sfx(i)*confire%ry(i) - confire%sfy(i)*confire%rx(i)
        enddo

        

        stressxy = stressxy / product(con%la) / 2.d0

        con%vx = con%vx + (con%sfx + con%fx) * dt2
        con%vy = con%vy + (con%sfy + con%fy) * dt2

        !call correc( con,dt)!,shearrate )
        


        if( con%strain .gt. strain0  ) then
            nframe = nframe + 1
            sumstress1 = sumstress1 + con%stress
            sumstressxy = sumstressxy + stressxy
            


            if(mod(step,100*iprint)==0) then
                nfile=nfile+1
                write( temppp,'(i6)' ) nfile
                filetempp = 'con'//trim(adjustl(temppp))//'_'//trim(adjustl(tempp))//'.txt'
                open(2,file=filetempp)
                write(2,*) con%la, attra, shearrate, con%strain, tb
                do i=1,sets%natom
                    write(2,'(9es26.16)') con%rx(i),con%ry(i),con%r(i),con%vx(i),con%vy(i),con%sfx(i),con%sfy(i),con%fx(i),con%fy(i)
                enddo
                close(2)
            endif

        
        endif

        if( mod(step,iprint) ==0 ) then
            write(22,'(3es26.16)') con%strain, (sumstress1+sumstressxy)/nframe, (sumstress1+sumstressxy)/nframe/shearrate
        endif
    enddo


    call trim_config(con)

    open(1,file=finalconfig)
    do i=1,sets%natom
        write(1,'(7es26.16)') con%rx(i), con%ry(i),con%r(i),con%vx(i),con%vy(i),con%fx(i),con%fy(i)
    enddo

contains
    subroutine readvar
        implicit none
            
        read *, tabsf
        read *, sets%seed
        read *, shearrate
        read *, dt
        read *, strain0
        read *, strain1
        read *, iprint
        read *, tb
    end subroutine

    subroutine testvar
        implicit none
            
        tabsf = 1d-7
        sets%seed = 1
        shearrate = 2d-5
        dt = 5d-3
        strain0 = 10
        strain1 = 100
        iprint = 10000
        tb = 1d-2
    end subroutine

end program

