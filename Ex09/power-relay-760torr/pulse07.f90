!
! CO2 REACTION TEST
! ZDPLASKIN
!
! 2020/04/17
!

program test_co2
    !
    !
    use ZDPlasKin
    use option !mass
    ! use input
    implicit none
    ! double precision, parameter :: 	gas_temperature  = 400.0d0, & ! gas temperature, K
                            	    ! density_ini_CO2   = 0d0,  & ! initial Ar density, cm-3
                                 ! frequency        = 13.56d6, & ! voltage frequency, Hz
                                    ! density_ini_elec = 1.0d9      ! initial electron density, cm-3
    double precision            :: 	time  = 0.0d0, time_end = 7.0d-6, dtime = 1.0d-15, & ! times, s
								    powertime = 3.0d-8, cooltime = 1.0d-6, ptime = 0.0d0, &
                                    reduced_field = 0d0,  & ! reduced electric field, Td

                                    e_dens = 0.0d0, & ! electron density, cm-3
                                    CO2_dens = 0.0d0, & ! cm-3
                                    O2_dens = 0.0d0, & ! cm-3
                                    CO_dens = 0.0d0, & ! cm-3
                                    CO2va_dens = 0.0d0, & ! cm-3
                                    CO2vb_dens = 0.0d0, & ! cm-3
                                    CO2v1_dens = 0.0d0, & ! cm-3
                                    CO2vd_dens = 0.0d0, & ! cm-3
                                    O_dens = 0.0d0, & ! cm-3


                                    pressure = 760, & ! Torr
                                    torr2pa = 101325/760, &
                                    ! kb = 1.38d-23, &
                            	    density_ini_CO2 = 0d0,  & ! initial Ar density, cm-3
                                    density_ini_elec = 1.0d9, &

                                    Va = 200d0, &
                                    Prf = 0d0, &
                                    Ppeak = 8d7, &
                                    frequency = 13.56d6, &
                                    gap = 1d-3, & ! m 
                                    omega = 0d0, &
                                    
                                    ele_charge = 1.602176487d-19, &
                                    m_ele = 9.10938215d-31, &  
                                    epsilon_CO2 = 1.0d0, &
                                    epsilon_quartz = 4d0, &
                                    ! v_m = 3d10, &

                                    v_m = 7.64d9, & ! s^-1
                                    v_dr = 0d0, &
                                    ! n0 = 2.4144d21, &
                                    ! n0 = 3.22d17, & ! cm^-3
                                    sigma = 0d0, &
                                    fpart = 0d0, &
                                    all_neutral = 0d0, &
                                    ele_mobility_n = 0d0, &
                                    current_dens = 0d0, &
                                    ele_temperature = 0d0, & 
                                    gas_temperature = 400.0d0, &
                                    dens = 0d0, &
									density_old(species_max), varstep


	double precision, parameter :: 	density_min = 1.0d5, & ! Min density solved
									absvar = 0.1, relvar = 0.1, &    ! Timestep adjusters, knowing the variation in densities
									max_field_err = 1.0d-3,&     ! max allowed variation of the field
									dtime_min = 1.0d-15     ! min value for the timestep


    integer                     ::  i, n, ios
    character( len = 20 )       ::  name

    call ZDPlasKin_init()
    call ZDPlaskin_set_config(QTPLASKIN_SAVE=.true.)

    ! call set_default
    call reactionAnalysis
    call getind()

    ! density_ini_CO2 = pressure * torr2pa /(kb*gas_temperature) * 1e-6
    density_ini_CO2 = 1.8d19

    open( 1, file = "baton-6" )

    read( 1, '( 1a30, 1pe30.8 )' )  name, time
    read( 1, '( 1a30, 1pe30.8 )' )  name, ptime
    read( 1, '( 1a30, 1pe30.8 )' )  name, ele_temperature
    read( 1, '( 1a30, 1pe30.8 )' )  name, gas_temperature
    read( 1, '( 1a30, 1pe30.8 )' )  name, reduced_field

    do while( .true. )
        read( 1, '( 1i5, 1a30, 1pe30.8 )', iostat = ios ) i, name, dens 
            density(i) = dens
        if( ios /= 0 ) exit 
    end do 

    close(1)

    ! fpart=0d0
    ! do i=1,nCO2-2
    !     fpart=fpart+CO2deg(i)*dexp(-CO2evib(i)*ec/kb/gas_temperature)
    ! enddo
    ! do i=1,nCO2-2
    !     density(indCO2(i))=density_ini_CO2*CO2deg(i)*dexp(-CO2evib(i)*ec/kb/gas_temperature)/fpart
    ! enddo

    ! call ZDPlasKin_set_density( 'e', density_ini_elec)
    ! call ZDPlasKin_set_density('CO2^+', density_ini_elec)

    ! reduced_field = 0d0

    ! omega = 2 * 3.1415926 * frequency
    ! sigma = ele_charge**2 * density_ini_elec / ( m_ele * v_m ) * v_m**2 / ( v_m**2 + omega**2 )
    ! reduced_field = sqrt( 2*Prf / sigma ) / density_ini_CO2 * 1e17

    call ZDPlasKin_set_conditions(GAS_TEMPERATURE=gas_temperature, REDUCED_FIELD=reduced_field)
    ! call ZDPlasKin_set_conditions( ELEC_TEMPERATURE =ele_temperature )

    write(*,'(4(1pe12.4))') time, dtime, reduced_field

    open( 10, file = "data-conditions.txt", status='OLD', position='APPEND')
    open( 11, file = "data-maindens.txt", status='OLD', position='APPEND' )

    do while(time .lt. time_end)
    	! call ZDPlasKin_timestep(time,dtime)

        call ZDPlasKin_get_density("e", e_dens)
        call ZDPlasKin_get_density_total( ALL_NEUTRAL = all_neutral )
        call ZDPlasKin_get_conditions( ELEC_MOBILITY_N = ele_mobility_n )
        current_dens = ele_charge * e_dens * ele_mobility_n * reduced_field * 1d-17

        if(ptime < ( powertime / 2 ) ) then
            ! reduced_field = 200.0
	        ! call ZDPlasKin_get_density("e", e_dens)
            ! call ZDPlasKin_get_density_total(ALL_NEUTRAL = all_neutral)
            ! call ZDPlasKin_get_conditions(ELEC_DRIFT_VELOCITY = v_dr)
            ! sigma = ele_charge * v_dr * e_dens / ( reduced_field*1d-17 * ( all_neutral ))
            ! reduced_field = sqrt( 2*Prf / sigma ) / ( all_neutral ) * 1d17

            Prf = ptime * ( Ppeak / ( powertime / 2 ) )
            ! current_dens = cond(4)
            ! reduced_field = Prf / ( all_neutral * current_dens ) * 1d17 
            reduced_field = sqrt(Prf / ( all_neutral * ele_charge * e_dens * ele_mobility_n ) ) * 1d17 

            ! call ZDPlasKin_set_conditions(REDUCED_FIELD=reduced_field)

        else if( ( ptime > ( powertime / 2 ) ) .and. ( ptime < powertime ) ) then 
            ! call ZDPlasKin_get_density_total(ALL_NEUTRAL = all_neutral)
            Prf = ptime * -1 *( Ppeak / ( powertime / 2 ) ) + 2 * Ppeak 
            ! current_dens = cond(4)
            ! reduced_field = Prf / ( all_neutral * current_dens ) * 1d17 
            reduced_field = sqrt(Prf / ( all_neutral * ele_charge * e_dens * ele_mobility_n ) ) * 1d17 

        else if( ( ptime > powertime ) .and. ( ptime < cooltime ) ) then 
            reduced_field = 0d0
            ! call ZDPlasKin_get_conditions( ELEC_TEMPERATURE = ele_temperature )
            ! if (ele_temperature < gas_temperature) then
            !     ele_temperature = gas_temperature
            !     call ZDPlasKin_set_conditions( ELEC_TEMPERATURE = ele_temperature )
            

        else
            ptime = 0
	        ! call ZDPlasKin_get_density("e", e_dens)
            ! call ZDPlasKin_get_density_total(ALL_NEUTRAL = all_neutral)
            ! call ZDPlasKin_get_conditions(ELEC_DRIFT_VELOCITY = v_dr)
            ! ! sigma = ele_charge * v_dr * e_dens / ( reduced_field*1d-17 * ( all_neutral ))
            ! sigma = ele_charge**2 * e_dens / ( m_ele * v_m ) * v_m**2 / ( v_m**2 + omega**2 )
            ! reduced_field = sqrt( 2*Prf / sigma ) / ( all_neutral ) * 1d17
            Prf = ptime * ( Ppeak / ( powertime / 2 ) )
            reduced_field = sqrt(Prf / ( all_neutral * ele_charge * e_dens * ele_mobility_n ) ) * 1d17 
            dtime = 1.0d-15
        endif

        call ZDPlasKin_set_conditions(REDUCED_FIELD=reduced_field)
        ! write(*,'(4(1pe12.4))') time, dtime, reduced_field

    	call ZDPlasKin_timestep(time, dtime)


        call ZDPlasKin_get_conditions(REDUCED_FIELD=reduced_field, ELEC_TEMPERATURE = ele_temperature )
        call ZDPlasKin_get_density("E", e_dens )
        call ZDPlasKin_get_density("CO2", CO2_dens)
        call ZDPlasKin_get_density("CO", CO_dens)
        call ZDPlasKin_get_density("O2", O2_dens)
        call ZDPlasKin_get_density("CO2(Va)", CO2va_dens )
        call ZDPlasKin_get_density("CO2(Vb)", CO2vb_dens  )
        call ZDPlasKin_get_density("CO2(V1)", CO2v1_dens)
        call ZDPlasKin_get_density("CO2(Vd)", CO2vd_dens)
        call ZDPlasKin_get_density("O", O_dens)

        write(10, '(6(1pe30.8))' ) time, reduced_field, ele_temperature, ele_temperature/11605, e_dens, Prf
        write(11, '(9(1pe30.8))' ) time, CO2_dens, CO2va_dens, CO2vb_dens, CO2v1_dens, CO2vd_dens, &
                                    CO_dens, O2_dens, O_dens


        ! time = time + dtime
        ! ptime = ptime + dtime

        ! call ZDPlasKin_timestep(time,dtime)
        write(*,'(4(1pe12.4))') time, dtime, reduced_field

        if ((ptime > 5.0d-10) .and. (ptime < 1.0d-8) ) then 
            dtime = 1.0d-12
        else if ((ptime > 1.0d-8) .and. (ptime < 1.0d-6) ) then 
            dtime = 1.0d-10
        else if ((ptime > 1.0d-6) .and. (ptime < 1.0d-4) ) then 
            dtime = 1.0d-8
        else if ( (ptime > 1.0d-4).and. (ptime < 1.0d-1) )  then 
            dtime = 1.0d-6
        else if (ptime > 1.0d-1)  then 
            ptime = 0d0
            dtime = 1.0d-15
        endif 

        time = time + dtime
        ptime = ptime + dtime

		! varstep = maxval(abs(density(:)-density_old(:))/(0.5d0*(density(:)+density_old(:))+ density_min))
		! dtime = dtime/(varstep/relvar + absvar)
		! dtime = min(dtime,varstep*max_field_err)
		! dtime = max(dtime,dtime_min)		

    end do

    open( 2, file = "baton-7" )
    call ZDPlasKin_get_conditions( ELEC_TEMPERATURE = ele_temperature, &
        GAS_TEMPERATURE = gas_temperature, REDUCED_FIELD = reduced_field )
    write( 2, '( (1a30), (1pe30.8) )' ) "TIME", time
    write( 2, '( (1a30), (1pe30.8) )' ) "PTIME", ptime
    write( 2, '( (1a30), (1pe30.8) )' ) "ELEC_TEMPERATURE", ele_temperature
    ! write( 1, '(2(1pe30.8))' ) "GAS_TEMPERATURE", gas_temperature
    ! write( 1, '(2(1pe30.8))' ) "REDUCED_FIELD", reduced_field
    write( 2, '( (1a30), (1pe30.8) )' ) "GAS_TEMPERATURE", gas_temperature
    write( 2, '( (1a30), (1pe30.8) )' ) "REDUCED_FIELD", reduced_field

    do i = 1, species_max
      write( 2, '( (1i5), (1a30), (1pe30.8) )') i, species_name(i),density(i)
    enddo

    close(2)

    close(10)
    close(11)

    ! write(*,'(/,A,$)') 'PRESS ENTER TO EXIT ...'
    ! read(*,*)
    
end program test_co2