!
! CO2 REACTION TEST
! ZDPLASKIN
!
! 2015/06/16
!

program co_2reac
    !
    ! declare variables and modules
    !
    use ZDPlasKin
    use option !mass
    use input
    implicit none
    !double precision, parameter :: reduced_field = 65d0

    double precision :: ivel, area,massflow

    double precision :: mass(species_max),omega
    double precision :: mavg
    double precision :: time = 0.0d0, dtime = 1.0d-12 ! times, s
    integer :: i,icount=0,icount_save=100,iflag
    !
    double precision :: emob, Te, all_neutral, pow, eden
    double precision :: EN
    !
    double precision :: energydep, sedep, powertotal, sei, jouletotal, particle, seinew
    double precision :: den0, vel0, vel1, pos, pres1, den1, rho0, rho1
    !
    double precision :: CO2sum,COsum,O2sum,CO2avge,COavge,CO2vta,CO2vt1,COvtemp
    double precision :: conversion = 0d0

    double precision :: density_old(species_max)
    double precision :: varstep
    double precision, parameter :: density_min = 1.0d5, & ! Min density solved
                                    absvar = 0.01, relvar = 0.01, & ! Timestep adjusters, knowing the variation in densities
                                    max_field_err = 1.0d-3,& ! max allowed variation of the field
                                     dtime_min = 1.0d-15 ! min value for the timestep
    double precision :: fpart
    double precision :: J,V,Vdr,species_source(species_max),source_terms(reactions_max)
    real(kind= 8) :: area_J, conversion_CO2 =0d0, conversion_CO =0d0, conversion_O2 =0d0
    real(kind= 8) :: sccm,Th
    
    integer                      :: k,intgeid_k,intgeii_k,intgeidi_k,intgeida_k,intgnr_k
    double precision             :: intgeid=0d0,&   !integerate ground electron impact dissociation
                                    intgeii=0d0,&   !integerate ground electron impact ionization
                                    intgeidi=0d0,&  !integerate ground electron impact dissociative ionization
                                    intgeida=0d0,&  !integerate ground electron impact dissociative attachment
                                    intgnr=0d0,&      !integerate ground neutral reactions
                                    intvseid=0d0,intvseii=0d0,intvseidi=0d0,intvseida=0d0,intvsnr=0d0,&      !symmetric vibration
                                    intvaeid=0d0,intvaeii=0d0,intvaeidi=0d0,intvaeida=0d0,intvanr=0d0      !symmetric vibration

    double precision :: powtotal,powelas,powinelas
    real(kind= 8) :: epower=0,powerelas=0,powerinelas=0
    
    integer :: flag(9)=1
    real(kind= 8) :: eedf(2,128)
    
    
    double precision :: t1_CPU, t2_CPU, tc_CPU

    ! double precision :: Va = 200d0, &
    !                     Prf = 100d0, &
    !                     frequency = 13.56d6, &
    !                     gap = 1d-3, & ! m 
    !                     epsilon_CO2 = 1.0d0, &
    !                     epsilon_quartz = 4d0, &
    !                     omega = 0d0, &
    !                     ele_charge = 1.602176487d-19, &
    !                     m_ele = 9.10938215d-31, &  
    !                     ! v_m = 3d10, &
    !                     v_m = 7.64d9, & ! s^-1
    !                     v_dr = 0d0, &
    !                     ! n0 = 2.4144d21, &
    !                     n0 = 3.22d17, & ! cm^-3
    !                     sigma = 0d0, &
    !                     ! n_heavy = 0d0
    !                     all_neutral = 0d0
    ! print
    !
    write(*,'(/,A)') 'CO2 REACTION TEST'
    !
    call CPU_TIME(t1_CPU) ! to get the processing time at the end
    !
    !
    ! initialization of ZDPlasKin
    !
    call ZDPlasKin_init()
    call ZDPlaskin_set_config(QTPLASKIN_SAVE=.true.)
    call set_default
    call reactionAnalysis
    call getind()
    !

    !print*,"begin?"
    !pause
    !
    ! mass initialization
    !
    mass = 0d0
    do i=1,species_max
        call getmass(species_name(i),mass(i)) !relative atomic mass
    enddo
    ! calculate sum of densities

    call sum_densities(CO2sum,COsum,O2sum,CO2avge,COavge,CO2vta,CO2vt1,COvtemp) 
    !
    ! set initial conditions and densities
    !
    den0 = pressure/kb/igastemp*1d-6 !cm-3
    area = pi*radius** 2 !cm2
    area_J = pi*(radius/2.)**2
    ivel = gasflow/area/60. !cm/s

    omega=2.0*pi*efreq

    fpart=0d0
    do i=1,nCO2-2
        fpart=fpart+CO2deg(i)*dexp(-CO2evib(i)*ec/kb/igastemp)
    enddo
    do i=1,nCO2-2
        density(indCO2(i))=den0*CO2deg(i)*dexp(-CO2evib(i)*ec/kb/igastemp)/fpart
    enddo

    call ZDPlasKin_set_density( 'e',ieden)
    call ZDPlasKin_set_density('CO2^+',ieden)

    rho0 = sum(density*mass*mc) !kg*cm-3
    massflow = rho0*gasflow/60. !massflow is constant in simulation, kg*s-1

    EN = voltage/length/den0*1d17
    vel0 = ivel

    pos = 0d0
    energydep = 0d0
    sedep = 0d0
    powertotal = 0d0
    sei = 0d0
    jouletotal = 0d0
    particle = 0d0
    seinew = 0d0
    vel1 = 0d0
    den1 = 0d0
    pres1 = 0d0
    rho1 = 0d0
    intrates = 0d0
    
    !
    !set condition of ZDPlasKin
    !
    call ZDPlasKin_set_conditions(GAS_TEMPERATURE=igastemp,REDUCED_FIELD=EN)
    !
    ! print column headers and initial values
    !

    ! screen, initiation
    write(*,'(4(A12))') 'Time_s', ( trim(species_name(i)), i = 1, species_max )
    write(*,'(4(1pe12.4))') time, density(:)
    !
    ! time integration
    !
    do while(time .lt. timeend)

        density_old(:) = density(:)
        vel0 = vel1

        if(time < 1d-8) then
            call ZDPlasKin_set_conditions(REDUCED_FIELD=EN)
        else
            call ZDPlasKin_set_conditions(ELEC_TEMPERATURE=igastemp,GAS_TEMPERATURE=igastemp)
        endif

        call ZDPlasKin_timestep(time,dtime)
        time = time + dtime


        call ZDPlasKin_get_density_total(ALL_NEUTRAL = all_neutral)
        call ZDPlasKin_get_conditions(ELEC_MOBILITY_N = emob, ELEC_TEMPERATURE = te,&
                                       ELEC_FREQUENCY_N = efreq, ELEC_POWER_N = pow,&
                                ELEC_POWER_ELASTIC_N = powelas,ELEC_DRIFT_VELOCITY = Vdr,&
                        elec_power_inelastic_N = powinelas,ELEC_EEDF = eedf,GAS_TEMPERATURE = Th)
        
        !
        ! velocity and position
        !                
        eden=density(species_electrons)

        ! ! -------------------                
        
        den1 = sum(density)
        pres1 = den1*kb*igastemp*1.0d6
        density = density/pres1*pressure !correct the density,是否考虑影响很小（无，1.24%；有，1.25%）
        rho1 = sum(density*mass*mc) !new mass density, kg*cm-3
        vel1 = massflow/rho1/area !new velocity, cm*s-1
        ! pos = pos+0.5*(vel0+vel1)*dtime !new position, cm
        ! sccm = pres1*1e-6/60/igastemp/kb
        ! !
        ! ! SEI
        ! ! method II:
        ! powertotal = powertotal+pow*eden*all_neutral*0.5*(vel0+vel1)*area_J*dtime ! eV*s-1
        ! powerelas = powerelas+powelas*eden*all_neutral*0.5*(vel0+vel1)*area_J*dtime
        ! powerinelas = powerinelas+powinelas*eden*all_neutral*0.5*(vel0+vel1)*area_J*dtime        
        ! sei = powertotal/(gasflow*sccm) !
        ! ! method III:
        ! jouletotal = jouletotal+pow*eden*all_neutral*0.5*(vel0+vel1)*area_J*dtime*dtime
        ! particle = particle+gasflow*sccm*dtime
        ! seinew = jouletotal/particle
        
        ! ! e-power
        ! epower = epower + eden*ec*EN*1d-17*all_neutral*Vdr*dtime

        ! ! conversion
        ! call sum_densities(CO2sum,COsum,O2sum,CO2avge,COavge,CO2vta,CO2vt1,COvtemp) 
        ! conversion = (1.0-CO2sum*vel1/(den0*ivel))*100d0
        ! conversion_CO2 = (CO2sum*vel1*area_J+den0*ivel*(area-area_J))/(den1*vel1*area_J+den0*ivel*(area-area_J))*100d0
        ! conversion_CO = COsum*vel1*area_J/(den1*vel1*area_J+den0*ivel*(area-area_J))*100d0
        ! conversion_O2 = O2sum*vel1*area_J/(den1*vel1*area_J+den0*ivel*(area-area_J))*100d0

        ! !------------------------------------
        ! ! circuit part
        ! !------------------------------------
        ! ! Current calculation (considering that the electrons are the only current carriers)

        ! ! --------------

        J = 1.6d-19 * area_J * eden * Vdr
        V = voltage - resistance * J

        EN = Voltage / ( length + resistance * J / ( EN*all_neutral/1.0d17 + 1.0d-99 )) / all_neutral * 1.0d17
        EN = 0.5d0 * ( EN + abs(EN) ) ! non-negative
        ! en=min(en,500.)
        !------------------------------------
        
        intrates=intrates+rrt*area_J*0.5*(vel0+vel1)*dtime
        call update_energy
        
        ! ! --------------------

        ! !------------------------------------------------------------------------------------
        ! ! Time-Integrated Net Rates
        ! !------------------------------------------------------------------------------------
        ! do k=i,reactions_max
        !     if(reaction_sign(k)=='bolsig:CO2->CO+O') then
        !         intgeid=intgeid+rrt(k)*dtime
        !         intvseid=intvseid+sum(rrt(k+1:k+4))*dtime
        !         intvaeid=intvaeid+sum(rrt(k+5:k+25))*dtime
        !     endif
        !     if(reaction_sign(k)=='bolsig:CO2->CO2^+') then
        !         intgeii=intgeii+rrt(k)*dtime
        !         intvseii=intvseii+sum(rrt(k+1:k+4))*dtime
        !         intvaeii=intvaeii+sum(rrt(k+5:k+25))*dtime
        !     endif
        !     if(reaction_sign(k)=='bolsig:CO2->CO^+' .or. &
        !             reaction_sign(k)=='bolsig:CO2->C^+' .or. &
        !                 reaction_sign(k)=='bolsig:CO2->O^+') then
        !         intgeidi=intgeidi+rrt(k)*dtime
        !         intvseidi=intvseidi+sum(rrt(k+1:k+4))*dtime
        !         intvaeidi=intvaeidi+sum(rrt(k+5:k+25))*dtime
        !     endif
        !     if(reaction_sign(k)=='bolsig:CO2->CO+O^-') then
        !         intgeida=intgeida+rrt(k)*dtime 
        !         intvseida=intvseida+sum(rrt(k+1:k+4))*dtime 
        !         intvaeida=intvaeida+sum(rrt(k+5:k+25))*dtime 
        !     endif
        !     if(reaction_sign(k)=='CO2+ANY_NEUTRAL+5.520_EV=>CO+O+ANY_NEUTRAL' .or. &
        !             reaction_sign(k)=='CO2+O+0.350_EV=>CO+O2' .or. &
        !                 reaction_sign(k)=='CO2+C=>2CO+5.640_EV') then
        !         intgnr=intgnr+rrt(k)*dtime
        !         intvsnr=intvsnr+sum(rrt(k+1:k+4))*dtime
        !         intvanr=intvanr+sum(rrt(k+5:k+25))*dtime                
        !     endif
        ! enddo

        ! ! -----------------
            
        ! output
        ! if(mod(icount,icount_save) .eq. 0) then
            ! screen print
            write(*,'(99(1pe14.6))') time, dtime, EN, eden, Te, Th, sei, conversion_CO 
            ! write in the output file
            
            ! ! --------

            ! !time temporal
            ! write(10,'(8(1pe12.4))') time, density(species_electrons), Te, 3377/log10(density(1)/density(6)),&
            !                             sei, conversion_CO2, conversion_CO, conversion_O2
            ! write(11,'(3(1pe12.4))') time, sei, conversion,EN
            ! write(12,'(8(1pe12.4))') time, te, CO2vt1, COvtemp, 6.23/log(density(1)/density(27)), CO2vta
            ! write(13,'(99(1pe12.4))') time, CO2sum, COsum, O2sum, density(O3_ind), density(O_ind),&
            !                         density(CO2_ind),density(CO_ind),density(O2_ind),&
            !                         density(C_ind),density(C2O_ind)
            
            ! !VT,VV',VV,VVco,Dm,Do,ed
            ! write(14,'(10(1pe12.4))') time,e_evx/(sccm*gasflow), (v_vt+v_vtx)/e_evx*100d0, v_vvt/e_evx*100d0, &
            !                  v_vv/e_evx*100d0, v_vvi1/e_evx*100d0, v_vdm/e_evx*100d0, v_vdo/e_evx*100d0,&
            !                     (ve_ed+ve_ei+ve_edi+ve_eda)/e_evx*100d0
            ! write(15,'(3(1pe12.4))') time, V, J 

            ! !electron energy
            ! write(16,'(99(1pe12.4))') time, e_evs/powerinelas*100d0, e_eva/powerinelas*100d0, e_eex/powerinelas*100d0,&
            !             e_ed/powerinelas*100d0,e_ei/powerinelas*100d0,e_edi/powerinelas*100d0,e_eda/powerinelas*100d0,&
            !                         CO_evx/powerinelas*100d0,CO_ed/powerinelas*100d0,CO_ei/powerinelas*100d0!,powerinelas/powerinelas*100d0
            ! write(17,'(99(1pe12.4))') time, powertotal/(gasflow*sccm), e_evs/(gasflow*sccm),e_eva/(gasflow*sccm), &
            !                         e_eex/(gasflow*sccm),e_ed/(gasflow*sccm),e_ei/(gasflow*sccm),e_edi/(gasflow*sccm),&
            !                         e_eda/(gasflow*sccm),(powertotal-e_evs-e_eva-e_eex)/(gasflow*sccm),&
            !                         (powertotal-e_evs-e_eva-e_eex-e_ed-e_ei-e_edi-e_eda-epower)/(gasflow*sccm)
                        
            ! !last
            ! rewind(18)
            ! write(18,'(99(1pe12.4))') time, powertotal/(gasflow*sccm), e_evs/(gasflow*sccm),e_eva/(gasflow*sccm), &
            !                         e_eex/(gasflow*sccm),e_ed/(gasflow*sccm),e_ei/(gasflow*sccm),e_edi/(gasflow*sccm),&
            !                         e_eda/(gasflow*sccm),(powertotal-e_evs-e_eva-e_eex)/(gasflow*sccm),&
            !                         (powertotal-e_evs-e_eva-e_eex-e_ed-e_ei-e_edi-e_eda)/(gasflow*sccm)
            
            ! !CO2 Va-Vb V1-V21
            ! write(19,'(99(1pe12.4))') time, density(1:28)
            
            ! !
            ! rewind(20)
            ! write(20,'(A)') "Time-Integrated Net Rates"
            ! write(20,'(8(1pe12.4))') intgeid,intgeii,intgeidi,intgeida,intgnr
            ! write(20,'(8(1pe12.4))') intvseid,intvseii,intvseidi,intvseida,intvsnr
            ! write(20,'(8(1pe12.4))') intvaeid,intvaeii,intvaeidi,intvaeida,intvanr
            ! !
            ! !call ZDPlasKin_get_conditions(elec_power_N=powtotal,elec_power_elastic_N=powelas,elec_power_inelastic_N=powinelas)
            ! !write(*,*) powtotal,powinelas,powelas,powtotal-powelas-powinelas
            
            ! ! ---------------

        ! endif
        
        ! ! -----------------
        ! !EEDF output
        ! do iflag=1,9
        !     if(time>=1d-12*10**(iflag-1).and.flag(iflag)==1) then
        !         write(30+iflag,'(2(1pe12.4))') eedf
        !         close(30+iflag)
        !         flag(iflag)=0
        !     endif
        ! enddo
        
        ! ! get the exact source terms
        ! call ZDPlasKin_get_rates(SOURCE_TERMS=species_source(:),REACTION_RATES=source_terms(:))
        ! !------------------------------
        ! ! new timestep
        ! !------------------------------
        ! ! varstep = maxval(abs(density(:)-density_old(:))/(0.5d0*(density(:)+density_old(:))+ density_min))
        ! ! dtime = dtime/(varstep/relvar + absvar)
        ! ! varstep = abs( voltage / ( resistance * J + 1.0d-99 ) - 1.0d0 ) * (eden + 1.0d-99) / ( abs(eden) + 1.0d-99 )
        ! ! dtime = min(dtime,varstep*max_field_err)
        ! ! dtime = max(dtime,dtime_min)
        ! ! !------------------------------

        ! ! if(time<3.0d-5) dtime=min(1e-9,dtime)

        ! !if (1.0e-3<=time) dtime=1e-8
        ! !if (9.2e-3<=time) dtime=1e-4
        ! !!----------------------------------------------------

        if ((time > 5.0d-10) .and. (time < 1.0d-8) ) then 
            dtime = 1.0d-11
        else if ((time > 1.0d-8) .and. (time < 1.0d-7) ) then 
            dtime = 2.0d-11
        else if ( (time > 1.0d-7) .and. (time < 1.0d-5) ) then 
            dtime = 1.0d-10
        else if ( (time > 1.0d-5) .and. (time < 1.0d-3)) then 
            dtime = 1.0d-9
        end if

        ! write(1,'(99(1pe12.4))') time, density(species_electrons), Te, Th, 3377/log10(density(1)/density(6)), &
        !      sei, conversion_CO2, conversion_CO, conversion_O2
        ! write(99,'(99(1pe14.6))',advance='no') time, dtime, EN
        ! write(99,'(99(1pe14.6))') V, J, J/area, V*J/(length*area*density(species_electrons))
        ! ! ------------------------  

        icount = icount + 1
    enddo
    
    !
    close(1)
    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    close(18)
    close(19)
    close(20)
    ! end
    !
    call CPU_TIME(t2_CPU)
    tc_CPU = (t2_CPU - t1_CPU) / 60.
    if(tc_CPU .lt. 1.5d0) then
        tc_CPU = t2_CPU - t1_CPU
        write(*,'(A,1pe12.5,A)') 'Time of calculations:', tc_CPU,' (s)'
    else
        tc_CPU = floor( (t2_CPU-t1_CPU) / 60. )
        write(*,'(A,1pe12.5,A)') 'Time of calculations:', tc_CPU,' (min)'
    endif

    write(*,'(/,A,$)') 'PRESS ENTER TO EXIT ...'
    read(*,*)

end program co_2reac
