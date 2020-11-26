module input
    implicit none
    double precision::timeend,radius,length,pressure,igastemp,voltage,resistance,gasflow,efreq,ieden
    namelist /INPUTDATA/ timeend,radius,length,pressure,igastemp,voltage,resistance,gasflow,efreq,ieden
        
    contains
    
    subroutine init()
        timeend=1.0d-3 !s
        !
        radius=0.02  !cm
        length=0.017
        !
        ! pressure=1.01325d5*0.995 * 10 / 760 !1.33d4        
        igastemp=300.0 !gas temperature,K
        !
        !voltage and resistance of electronic circuit 
        ! voltage=4500  !V
        ! resistance=1e6 !Ohm
        !
        gasflow=1.0d2  !ml/min
        !
        !frequency of voltage, is 2.54d9 in MW.
        efreq=0 !2.54d9
        !
        !initial value of electron density
        ieden=1.0d9 !cm^-3
        !  


        ! time  = 0.0d0, time_end = 2.0d-2, dtime = 1.0d-15, & ! times, s
        ! powertime = 1.0d-2, cooltime = 2.0d-2, ptime = 0.0d0, &
        ! reduced_field = 0d0,  & ! reduced electric field, Td

        ! e_dens = 0.0d0, & ! electron density, cm-3
        ! CO2_dens = 0.0d0, & ! cm-3
        ! O2_dens = 0.0d0, & ! cm-3
        ! CO_dens = 0.0d0, & ! cm-3
        ! CO2va_dens = 0.0d0, & ! cm-3
        ! CO2vb_dens = 0.0d0, & ! cm-3
        ! CO2v1_dens = 0.0d0, & ! cm-3
        ! CO2vd_dens = 0.0d0, & ! cm-3
        ! O_dens = 0.0d0, & ! cm-3
        ! t = 0.0d0, &

        ! pressure = 10, & ! Torr
        ! torr2pa = 101325/760, &
        ! ! kb = 1.38d-23, &
        ! density_ini_CO2 = 0d0,  & ! initial Ar density, cm-3
        ! density_ini_elec = 1.0d9, &

        ! Va = 200d0, &
        ! Prf = 0d0, &
        ! Ppeak = 1d2, &
        ! frequency = 13.56d6, &
        ! gap = 1d-3, & ! m 
        ! omega = 0d0, &
        
        ! ele_charge = 1.602176487d-19, &
        ! m_ele = 9.10938215d-31, &  
        ! epsilon_CO2 = 1.0d0, &
        ! epsilon_quartz = 4d0, &
        ! ! v_m = 3d10, &

        ! v_m = 7.64d9, & ! s^-1
        ! v_dr = 0d0, &
        ! ! n0 = 2.4144d21, &
        ! ! n0 = 3.22d17, & ! cm^-3
        ! sigma = 0d0, &
        ! fpart = 0d0, &
        ! all_neutral = 0d0, &
        ! ele_mobility_n = 0d0, &
        ! current_dens = 0d0, &
        ! ele_temperature = 400.0d0, & 
        ! gas_temperature = 400.0d0, &
        ! dens = 0.0d0, &
    end subroutine
    
end
    