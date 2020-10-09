module option
    use ZDPlasKin
    implicit none
    double precision, parameter :: kb=1.3806488d-23, pi=3.14159265d0
    double precision, parameter :: ec=1.602176585d-19, me=9.10938215d-31, mc=1.667d-27
    double precision, parameter :: maxen=200.0

    integer,parameter::maxlev=70
    integer,parameter::maxreac=3000 ! (2x21x63 for CO2-CO VV)
    double precision,dimension(reactions_max)::intrates

    integer::nCO2,nCO,nO2,nevx,neva,nevs,neex,ned,nei,nedi,neda,&
             nvdm,nvdo,nvt,nvtx,nvvt,nvv,nvvi,nol,nevag
    integer,dimension(maxlev)::indCO2,indCO,indO2
    !double precision::mass(species_max),shratio(species_max)
    double precision,dimension(maxlev)::wCO2,wCO,wO2
    integer,dimension(maxreac)::ievx,ieva,ievs,ieex,ied,iei,iedi,ieda,&
                                    ivdm,ivdo,ivt,ivtx,ivvt,ivv,ivvi,iol,ievag
    double precision,dimension(maxreac)::wevx,weva,wevs,weex,wed,wei,wedi,weda,wved,wvei,wvedi,wveda,&
                                            wvdm,wvdo,wvt,wvtx,wvvt,wvv,wvvi1,wvvit
    double precision::e_evx,e_eva,e_evs,e_eex,e_ed,e_ei,e_edi,e_eda,ve_ed,ve_ei,ve_edi,ve_eda,&
                        v_vdm,z_vdm,v_vdo,z_vdo,v_vt,v_vtx,v_vvt,v_vv,v_vvi1,v_vvit,z_ol,nu_evag
    
    integer::nCO_evx,nCO_eex,nCO_ed,nCO_ei,nCO_edi
    integer,dimension(maxreac)::iCO_evx,iCO_eex,iCO_ed,iCO_ei,iCO_edi
    double precision,dimension(maxreac)::wCO_evx,wCO_eex,wCO_ed,wCO_ei,wCO_edi
    double precision::CO_evx,CO_eex,CO_ed,CO_ei,CO_edi
    
    double precision, dimension(28), parameter::CO2evib=(/ 0.0, 0.08276, 0.16570, 0.24855, 0.33140,&
                    0.29126, 0.57943, 0.86451, 1.14650, 1.42539, 1.70119, 1.97390, 2.24352, 2.51005, 2.77348,&
3.03382, 3.29107, 3.54523, 3.79629, 4.04426, 4.28914, 4.53093, 4.76962, 5.00523, 5.23774, 5.46716, 0.0, 0.0/)
    double precision, dimension(28), parameter::CO2deg=(/ 1.0, 1.0, 3.0, 3.0, 6.0,&
                          1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, &
                            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)

    double precision, dimension(2), parameter::CO2eex=(/ 6.23, 11.1 /)

    double precision, dimension(15), parameter::COevib=(/ 0.0, 0.26575, 0.52819, 0.78735, 1.04326,&
                                    1.29593, 1.54539, 1.79166, 2.03476, 2.27472, 2.51157, 0.0, 0.0, 0.0, 0.0/)  

    logical::autosave,gasheat

    integer :: CO2_ind, CO_ind, O2_ind, O3_ind, O_ind, C_ind, C2O_ind

    contains

    ! get indices of CO2 densities for summation
    subroutine getind()
        implicit none
        integer:: i  

        call ZDPlasKin_get_species_index("CO2",CO2_ind)
        call ZDPlasKin_get_species_index("CO",CO_ind)
        call ZDPlasKin_get_species_index("O2",O2_ind)
        call ZDPlasKin_get_species_index("O3",O3_ind)
        call ZDPlasKin_get_species_index("O",O_ind)
        call ZDPlasKin_get_species_index("C",C_ind)
        call ZDPlasKin_get_species_index("C2O",C2O_ind)

        nCO2=0
        nCO=0
        nO2=0

        do i=1,species_max
            if (species_name(i)=='CO2' .or. species_name(i)(1:4)=='CO2(') then
                nCO2=nCO2+1
                indCO2(nCO2)=i
            endif
            if (species_name(i)=='CO' .or. species_name(i)(1:3)=='CO(') then
                nCO=nCO+1
                indCO(nCO)=i
            endif
            if (species_name(i)=='O2' .or. species_name(i)(1:3)=='O2(') then
                nO2=nO2+1
                indO2(nO2)=i
            endif
        enddo
    end subroutine

    subroutine getmass(species,mass)
        character(species_length)::species
        double precision,intent(out)::mass
        integer::i,ic
        integer::ibeg(0:20),iend(0:20),mult(0:20)

        i=0            !position of chemical equation
        ic=0           !number of letters(element) in equation
        ibeg=0         !index of begin position
        iend=0         !index of end position
        mult=0         !number of each letter
        do while (i.lt.species_length)
            i=i+1
            if ((ichar(species(i:i)).ge.65) .and. (ichar(species(i:i)).le.90)) then ! is uppercase letter
                if (ic>0 .and. iend(ic)==0) iend(ic)=i-1
                ic=ic+1
                ibeg(ic)=i
            endif
            if ((ichar(species(i:i)).ge.48) .and. (ichar(species(i:i)).le.57)) then ! is number 
                if (iend(ic)==0) iend(ic)=i-1
                mult(ic)=mult(ic)*10+ichar(species(i:i))-48  ! decode stoichiometry number
            endif
            if ((species(i:i)==' ') .or. (species(i:i)=='(') .or. (species(i:i)=='^')) then
                if (iend(ic)==0) iend(ic)=i-1
                exit
            endif
        enddo
        if (iend(ic)==0) iend(ic)=i

        mass=0.0
        do i=1,ic
            if (mult(i)==0.0) mult(i)=1.0
            mass = mass + mult(i)*mass_element(species(ibeg(i):iend(i)))
        enddo

    end subroutine


    function mass_element(element)
        character(*)::element
        double precision::mass_element

        select case(element)
        case('H')
            mass_element=1.008
        case('C')
            mass_element=12.011
        case('N')
            mass_element=14.007
        case('O')
            mass_element=15.999
            case default 
            mass_element=0.0
        end select
    end function


    Subroutine TrimStringNumber(c)
    !// Trim real String program
    !// Fortran Coder : Gao [ gao@fcode.cn ]
    !// http://www.fcode.cn
    
        character( Len = * ) :: c
        integer :: n
        integer :: i , j , k
        c = AdjustL( c )
        n = Len_Trim( c )
        i = index( c , "E" )
        j = index( c , "e" )
        k = max( i , j )

        if( c((k+1):(k+1))=='+' ) then !//È¥µôÖ¸ÊýÀïµÄ+
            c(k+1:)=c(k+2:)
            c(n:n)=' '
        endif

        if ( k > 1 ) then !// Ö¸Êý
            j = index( c(:k-1) , "." , back = .true. )
            if ( j < 1 ) return
            Do i = k-1 , j+1 , -1
                if( c(i:i) /= '0' ) exit
            End Do    
            if ( i /= j ) i = i + 1
            c(i:n) = c(k:n)

        else !// ³£¹æÐ¡Êý
            j = index( c , "." , back = .true. )
            if ( j < 1 ) return
            Do i = n , j+1 , -1
                if( c(i:i) /= '0' ) exit
            End Do
            if ( i /= j ) i = i + 1
            c(i:n) = ' '
        end if
    End Subroutine TrimStringNumber



    ! calculate and set summed densities
    subroutine sum_densities(CO2sum,COsum,O2sum,CO2evavg,COevavg,CO2Tva,CO2Tv1,COTv)
        implicit none
        integer,parameter::rmax=10 ! number of vibrational levels used to calculate vibrational temperature
        double precision::CO2sum,COsum,O2sum,CO2evavg,COevavg,CO2Tva,CO2Tv1,COTv
        integer::i,iva,iv1

        CO2sum=0d0
        CO2evavg=0d0
        do i=1,nCO2
            if (wCO2(i)==CO2evib(2)) iva=i
            if (wCO2(i)==CO2evib(6)) iv1=i            
            CO2sum=CO2sum+density(indCO2(i))
            CO2evavg=CO2evavg+wCO2(i)*density(indCO2(i))
        enddo    
        CO2evavg=CO2evavg/CO2sum*ec/kb
        CO2Tva=ec/kb*wCO2(iva)/dlog(density(indCO2(1))/density(indCO2(iva)))
        CO2Tv1=ec/kb*wCO2(iv1)/dlog(density(indCO2(1))/density(indCO2(iv1)))
        COsum=0d0
        COevavg=0d0
        do i=1,nCO
            COsum=COsum+density(indCO(i))
            COevavg=COevavg+wCO(i)*density(indCO(i))
        enddo    
        COevavg=COevavg/COsum*ec/kb
        COTv=ec/kb*wCO(2)/dlog(density(indCO(1))/density(indCO(2)))
        O2sum=0d0
        do i=1,nO2
            O2sum=O2sum+density(indO2(i))
        enddo    
    end subroutine

    subroutine update_energy
    !    double precision::dV
        integer::i

        nu_evag=0d0 ! electron excitation frequency
        do i=1,nevag
            nu_evag=nu_evag+intrates(ievag(i))
        enddo

        e_evx=0d0 ! electron energy lost by vibrational excitation
        do i=1,nevx
            e_evx=e_evx+intrates(ievx(i))*wevx(i)
        enddo
        e_eva=0d0
        do i=1,neva
            e_eva=e_eva+intrates(ieva(i))*weva(i)
        enddo
        e_evs=0d0
        do i=1,nevs
            e_evs=e_evs+intrates(ievs(i))*wevs(i)
        enddo

        e_eex=0d0 ! electron energy lost by electronic excitation
        do i=1,neex
            e_eex=e_eex+intrates(ieex(i))*weex(i)
        enddo

        e_ed=0d0 ! electron energy lost by dissociation
        ve_ed=0d0
        do i=1,ned
            e_ed=e_ed+intrates(ied(i))*wed(i)
            ve_ed=ve_ed+intrates(ied(i))*wved(i)
        enddo
        e_ei=0d0
        ve_ei=0d0
        do i=1,nei
            e_ei=e_ei+intrates(iei(i))*wei(i)
            ve_ei=ve_ei+intrates(iei(i))*wvei(i)
        enddo
        e_edi=0d0 ! electron energy lost by dissociation
        ve_edi=0d0
        do i=1,nedi
            e_edi=e_edi+intrates(iedi(i))*wedi(i)
            ve_edi=ve_edi+intrates(iedi(i))*wvedi(i)
        enddo
        e_eda=0d0
        ve_eda=0d0
        do i=1,neda
            e_eda=e_eda+intrates(ieda(i))*weda(i)
            ve_eda=ve_eda+intrates(ieda(i))*wveda(i)
        enddo    

        v_vdm=0d0 ! vibrational energy lost by dissociation by any neutral
        z_vdm=0d0 ! total energy spent for dissociation by any neutral
        do i=1,nvdm
            v_vdm=v_vdm+intrates(ivdm(i))*wvdm(i)
            z_vdm=z_vdm+intrates(ivdm(i))  
        enddo
        !  z_vdm=z_vdm*5.5  !!!
        v_vdo=0d0 ! vibrational energy lost by dissociation by O
        z_vdo=0d0 ! total energy spent for dissociation by O
        do i=1,nvdo
            v_vdo=v_vdo+intrates(ivdo(i))*wvdo(i)
            z_vdo=z_vdo+intrates(ivdo(i))  
        enddo
        !  z_vdo=z_vdo*0.36  !!! 
        v_vt=0d0 ! vibrational energy lost by VT relaxation
        do i=1,nvt
            v_vt=v_vt+intrates(ivt(i))*wvt(i)
        enddo
        v_vtx=0d0 ! vibrational energy lost by VT relaxation due to a priori assumptions
        do i=1,nvtx
            v_vtx=v_vtx+intrates(ivtx(i))*wvtx(i)
        enddo
        v_vvt=0d0 ! vibrational energy lost by VV' relaxation
        do i=1,nvvt
            v_vvt=v_vvt+intrates(ivvt(i))*wvvt(i)
        enddo
        v_vv=0d0 ! vibrational energy lost by VV relaxation among CO2
        do i=1,nvv
            v_vv=v_vv+intrates(ivv(i))*wvv(i)
        enddo
        v_vvi1=0d0 ! CO2 vibrational energy lost by VV relaxation with CO
        v_vvit=0d0 ! total vibrational energy lost by VV relaxation with CO
        do i=1,nvvi
            v_vvi1=v_vvi1+intrates(ivvi(i))*wvvi1(i)
            v_vvit=v_vvit+intrates(ivvi(i))*wvvit(i)
        enddo
        z_ol=0d0
        do i=1,nol
            z_ol=z_ol+intrates(iol(i))
        enddo
        
        CO_evx=0d0
        do i=1,nCO_evx
            CO_evx=CO_evx+intrates(iCO_evx(i))*wCO_evx(i)
        enddo
        CO_ed=0d0
        do i=1,nCO_ed
            CO_ed=CO_ed+intrates(iCO_ed(i))*wCO_ed(i)
        enddo
        CO_ei=0d0
        do i=1,nCO_ei
            CO_ei=CO_ei+intrates(iCO_ei(i))*wCO_ei(i)
        enddo

    end subroutine
    ! finds and saves indexes to interesting reactions that are used for analysis
    subroutine reactionAnalysis
        implicit none
        integer::i
        character(reactions_length+1) :: rsign
        character(11),parameter::anypar='ANY_NEUTRAL'
        character(species_length)::species(6),anyneutral
        integer::nspec,firstprod,rtype,i1,i2,i3,i4

        anyneutral=anypar(1:species_length)

        ! get indices of densities for summation
        nCO2=0
        nCO=0
        nO2=0
        do i=1,species_max
            i1=getCO2ind(species_name(i))
            if (i1>0 .and. i1.le.28) then
                nCO2=nCO2+1
                indCO2(nCO2)=i
                wCO2(nCO2)=CO2evib(i1)
            endif
            i2=getCOind(species_name(i))
            if (i2>0 .and. i2.le.15) then
                nCO=nCO+1
                indCO(nCO)=i
                wCO(nCO)=COevib(i2)
            endif
            if (species_name(i)=='O2' .or. species_name(i)(1:3)=='O2(') then
                nO2=nO2+1
                indO2(nO2)=i
            endif
        enddo

        nevx=0; neex=0; ned=0; nvdm=0; nvdo=0; nvt=0; nvtx=0; nvvt=0; nvv=0; nvvi=0; 
        do i=1,reactions_max ! loop through all reactions
            rsign=reaction_sign(i)

            call readReaction(rsign,rtype,nspec,species,firstprod)
            !intgeid_k,intgeii_k,intgeidi_k,intgeida_k,intgnr_k
            if(i==1477) then
                print*,rsign,rtype,nspec,species
                continue
            endif
            
            
            if (rtype==1 .and. nspec==2) then
                i1=getCO2ind(species(1))
                i2=getCO2ind(species(2))
                ! electron impact excitation of vibrational levels
                if (i1>0 .and. i2>0 .and. i1.le.26 .and. i2.le.26) then 
                    nevx=nevx+1
                    ievx(nevx)=i
                    wevx(nevx)=CO2evib(i2)-CO2evib(i1)
                    if (i2.ge.2 .and. i2.le.5) then
                        nevs=nevs+1	
                        ievs(nevs)=i
                        wevs(nevs)=CO2evib(i2)-CO2evib(i1)
                    else
                        neva=neva+1	
                        ieva(neva)=i
                        weva(neva)=CO2evib(i2)-CO2evib(i1)
                    endif
                    if (i1>1 .and. i1.le.5) then ! excitation from symmetric levels
                        wevx(nevx)=CO2evib(i2)
                        weva(neva)=CO2evib(i2)
                        nvtx=nvtx+1
                        ivtx(nvtx)=i
                        wvtx(nvtx)=CO2evib(i1)  ! vibrational energy lost due to a priori assumptions
                    endif
                    if (i1==1 .and. i2>5) then ! only from ground level to asymmetric levels
                        nevag=nevag+1
                        ievag(nevag)=i
                    endif
                endif
                ! electron impact excitation of electronic levels
                if (i1>0 .and. i2.ge.27 .and. i1.le.26) then
                    neex=neex+1
                    ieex(neex)=i
                    weex(neex)=CO2eex(i2-26)
                endif
                
            endif
            ! electron impact dissociation
            if (rtype==1 .and. nspec==3) then
                i1=getCO2ind(species(1))
                if (i1>0 .and. species(2)=='CO' .and. species(3)=='O') then
                    ned=ned+1
                    ied(ned)=i
                    wed(ned)=11.9-CO2evib(i1)
                    wved(ned)=CO2evib(i1)
                endif
            endif
            ! electron impact ionization
            if (rtype==1 .and. nspec==2) then
                i1=getCO2ind(species(1))
                if (i1>0 .and. species(2)=='CO2^') then
                    nei=nei+1
                    iei(nei)=i
                    wei(nei)=13.8-CO2evib(i1)
                    wvei(nei)=CO2evib(i1)
                endif
            endif
            ! electron impact dissociative ionization
            if (rtype==1 .and. nspec==2) then
                i1=getCO2ind(species(1))
                if (i1>0 .and. species(2)=='CO^') then
                    nedi=nedi+1
                    iedi(nedi)=i
                    wedi(nedi)=19.5-CO2evib(i1)
                    wvedi(nedi)=CO2evib(i1)
                endif
                if (i1>0 .and. species(2)=='C^') then
                    nedi=nedi+1
                    iedi(nedi)=i
                    wedi(nedi)=27.8-CO2evib(i1)
                    wvedi(nedi)=CO2evib(i1)
                endif
                if (i1>0 .and. species(2)=='O^') then
                    nedi=nedi+1
                    iedi(nedi)=i
                    wedi(nedi)=19.1-CO2evib(i1)
                    wvedi(nedi)=CO2evib(i1)
                endif
            endif
            ! electron impact dissociative attachment
            if (rtype==1 .and. nspec==3) then
                i1=getCO2ind(species(1))
                if (i1>0 .and. species(2)=='CO' .and. species(3)=='O^-') then
                    neda=neda+1
                    ieda(neda)=i
                    weda(neda)=3.85d0-CO2evib(i1)
                    wveda(neda)=CO2evib(i1)
                endif
            endif
            
            ! CO2 dissociation by any neutral
            if (rtype==0 .and. nspec>=5 .and. species(firstprod)=='CO' .and. species(firstprod+1)=='O') then
                i1=getCO2ind(species(1))
                if (i1>0 .and. i1.le.26) then
                    nvdm=nvdm+1
                    ivdm(nvdm)=i
                    wvdm(nvdm)=CO2evib(i1)
                endif
            endif
            ! CO2 dissociation by O
            if (rtype==0 .and. nspec>=4 .and. species(firstprod)=='CO' .and. species(firstprod+1)=='O2') then
                i1=getCO2ind(species(1))
                if (i1>0 .and. i1.le.26) then
                    nvdo=nvdo+1
                    ivdo(nvdo)=i
                    wvdo(nvdo)=CO2evib(i1)
                endif
            endif
            ! CO2 VT relaxation
            if (rtype==0 .and. nspec>=4 .and. &
            (species(2)=="ANY_NEUT".or.species(2)=="CO2".or.species(2)=="CO".or.species(2)=="O2") .and. &
            (species(firstprod+1)=="ANY_NEUT".or.species(4)=="CO2".or.species(4)=="CO".or.species(4)=="O2") )then
                i1=getCO2ind(species(1))
                i2=getCO2ind(species(firstprod))
                if (i1>0 .and. i2>0 .and. i1.le.26 .and. i2.le.26) then 
                    nvt=nvt+1
                    ivt(nvt)=i
                    wvt(nvt)=CO2evib(i1)-CO2evib(i2)
                endif
            endif
            
            if (rtype==0 .and. nspec>=4) then
                i1=getCO2ind(species(1))
                i2=getCO2ind(species(2))
                i3=getCO2ind(species(firstprod))
                i4=getCO2ind(species(firstprod+1))
                ! VV' relaxation between asymmetric and symmetric mode of CO2
                if (i1>0 .and. i1.le.26 .and. i2>0 .and. i2.le.5 .and. i4>0 .and. i4.le.5 .and. i2/=i4) then 
                    nvvt=nvvt+1
                    ivvt(nvvt)=i
                    wvvt(nvvt)=CO2evib(i1)+CO2evib(i2)-CO2evib(i3)-CO2evib(i4)
                endif
                ! VV relaxation among CO2
                if   ((i1==1 .or. i1.ge.6) .and. i1.le.26 &
                .and.(i2==1 .or. i2.ge.6) .and. i2.le.26 &
                .and.(i3==1 .or. i3.ge.6) .and. i3.le.26 &
                .and.(i4==1 .or. i4.ge.6) .and. i4.le.26 &
                .and. i2/=i4 ) then 
                    nvv=nvv+1
                    ivv(nvv)=i
                    wvv(nvv)=CO2evib(i1)+CO2evib(i2)-CO2evib(i3)-CO2evib(i4)
                endif
                i2=getCOind(species(2))
                i4=getCOind(species(firstprod+1))
                ! VV relaxation between CO2 and CO
                if   ((i1==1 .or. i1.ge.6) .and. i1.le.26 &
                .and. i2>0 .and. i2.le.11 &
                .and.(i3==1 .or. i3.ge.6) .and. i3.le.26 &
                .and. i4>0 .and. i4.le.11) then 
                    nvvi=nvvi+1
                    ivvi(nvvi)=i
                    wvvi1(nvvi)=CO2evib(i1)-CO2evib(i3) ! lost vibrational energy of CO2
                    wvvit(nvvi)=CO2evib(i1)+COevib(i2)-CO2evib(i3)-COevib(i4) ! energy lost to translation
                endif
            endif
            ! O radical losses
            if (rtype==0 .and. nspec>=4 .and. (species(1)=='O' .or. species(2)=='O')) then
                i1=getCO2ind(species(1))
                if (i1<0) then
                    nol=nol+1
                    iol(nol)=i
                endif
            endif
            
            !CO reactions
            if(rtype==1) then
                if(species(1)(1:2)=='CO' .and. species(2)(1:4)=='CO(V') then
                    i1=getCOind(species(1))
                    i2=getCOind(species(2))
                    nCO_evx=nCO_evx+1
                    iCO_evx(nCO_evx)=i
                    wCO_evx(nCO_evx)=COevib(i2)-COevib(i1)
                endif
                if(species(1)(1:2)=='CO' .and. species(2)=='C' .and. species(3)=='O') then
                    i1=getCOind(species(1))
                    nCO_ed=nCO_ed+1
                    iCO_ed(nCO_ed)=i
                    wCO_ed(nCO_ed)=13.5-COevib(i1)
                endif
                if(species(1)(1:2)=='CO' .and. species(1)(3:3)/='2' .and. species(2)=='CO^') then
                    i1=getCOind(species(1))
                    nCO_ei=nCO_ei+1
                    iCO_ei(nCO_ei)=i
                    wCO_ei(nCO_ei)=17.5-COevib(i1)
                endif
            endif
            

        enddo ! end of reaction loop

        ! write all reactions to file    
        open(111,file='data_reactionAnalysis.txt')
        write(111,'(A)') 'CO2 species'
        do i=1,nCO2
            write(111,'(I4,A,A,F8.4)') i,'  ',species_name(indCO2(i)),wCO2(i)
        enddo
        write(111,'(/,A)') 'CO species'
        do i=1,nCO
            write(111,'(I4,A,A,F8.4)') i,'  ',species_name(indCO(i)),wCO(i)
        enddo
        write(111,'(/,A)') 'O2 species'
        do i=1,nO2
            write(111,'(I4,A,A)') i,'  ',species_name(indO2(i))
        enddo

        ! write all reactions to file    

        write(111,'(/,A)') 'Electron impact excitation of vibrational levels (CO_evx)'
        do i=1,nCO_evx
            write(111,'(I4,A,A,F8.4)') i,'  ',reaction_sign(iCO_evx(i)),wCO_evx(i)
        enddo
        write(111,'(/,A)') 'Electron impact dissociation (CO_ed)'
        do i=1,nCO_ed
            write(111,'(I4,A,A,F8.4)') i,'  ',reaction_sign(iCO_ed(i)),wCO_ed(i)
        enddo
        write(111,'(/,A)') 'Electron impact ionization (CO_ei)'
        do i=1,nCO_ei
            write(111,'(I4,A,A,F8.4)') i,'  ',reaction_sign(iCO_ei(i)),wCO_ei(i)
        enddo
        
        write(111,'(/,A)') 'Electron impact excitation of vibrational levels (evx)'
        do i=1,nevx
            write(111,'(I4,A,A,F8.4)') i,'  ',reaction_sign(ievx(i)),wevx(i)
        enddo
        write(111,'(/,A)') 'Electron impact excitation of asymmetric vibrational levels (eva)'
        do i=1,neva
            write(111,'(I4,A,A,F8.4)') i,'  ',reaction_sign(ieva(i)),weva(i)
        enddo
        write(111,'(/,A)') 'Electron impact excitation of symmetric vibrational levels (evs)'
        do i=1,nevs
            write(111,'(I4,A,A,F8.4)') i,'  ',reaction_sign(ievs(i)),wevs(i)
        enddo
        write(111,'(/,A)') 'Electron impact excitation of electronic levels (eex)'
        do i=1,neex
            write(111,'(I4,A,A,F8.4)') i,'  ',reaction_sign(ieex(i)),weex(i)
        enddo
        write(111,'(/,A)') 'Electron impact excitation from ground to asymmetric levels (evag)'
        do i=1,nevag
            write(111,'(I4,A,A)') i,'  ',reaction_sign(ievag(i))
        enddo
        write(111,'(/,A)') 'Electron impact dissociation (ed)'
        do i=1,ned
            write(111,'(I4,A,A,F8.4)') i,'  ',reaction_sign(ied(i)),wed(i)
        enddo
        write(111,'(/,A)') 'Electron impact ionization (ei)'
        do i=1,nei
            write(111,'(I4,A,A,F8.4)') i,'  ',reaction_sign(iei(i)),wei(i)
        enddo
        write(111,'(/,A)') 'Electron impact dissociative ionization (edi)'
        do i=1,nedi
            write(111,'(I4,A,A,F8.4)') i,'  ',reaction_sign(iedi(i)),wedi(i)
        enddo
        write(111,'(/,A)') 'Electron impact dissociative attachment (eda)'
        do i=1,neda
            write(111,'(I4,A,A,F8.4)') i,'  ',reaction_sign(ieda(i)),weda(i)
        enddo
        write(111,'(/,A)') 'Dissociation by any neutral (vdm)'
        do i=1,nvdm
            write(111,'(I4,A,A,F8.4)') i,'  ',reaction_sign(ivdm(i)),wvdm(i)
        enddo
        write(111,'(/,A)') 'Dissociation by O (vdo)'
        do i=1,nvdo
            write(111,'(I4,A,A,F8.4)') i,'  ',reaction_sign(ivdo(i)),wvdo(i)
        enddo
        write(111,'(/,A)') 'VT relaxation (vt)'
        do i=1,nvt
            write(111,'(I4,A,A,F8.4)') i,'  ',reaction_sign(ivt(i)),wvt(i)
        enddo
        write(111,'(/,A)') 'VT relaxation due to a priori assumptions (vtx)'
        do i=1,nvtx
            write(111,'(I4,A,A,F8.4)') i,'  ',reaction_sign(ivtx(i)),wvtx(i)
        enddo
        write(111,'(/,A)') 'VV'' relaxation between asymmetric and symmetric mode (vvt)'
        do i=1,nvvt
            write(111,'(I4,A,A,F8.4)') i,'  ',reaction_sign(ivvt(i)),wvvt(i)
        enddo
        write(111,'(/,A)') 'VV relaxation among CO2 (vv)'
        do i=1,nvv
            write(111,'(I4,A,A,F8.4)') i,'  ',reaction_sign(ivv(i)),wvv(i)
        enddo
        write(111,'(/,A)') 'VV relaxation with CO (vvi)'
        do i=1,nvvi
            write(111,'(I4,A,A,F8.4,F8.4)') i,'  ',reaction_sign(ivvi(i)),wvvi1(i),wvvit(i)
        enddo
        write(111,'(/,A)') 'O loss reactions (ol)'
        do i=1,nol
            write(111,'(I4,A,A)') i,'  ',reaction_sign(iol(i))
        enddo

        close(111)

    end subroutine

    ! rtype = 0 - standard reaction, 1 - bolsig reaction, 2 - bolsig elastic reaction
    subroutine readReaction(rsign,rtype,nspec,species,firstprod)
        implicit none
        character(reactions_length+1),intent(in) :: rsign
        character(species_length),intent(out) :: species(6)
        integer,intent(out)::nspec,firstprod,rtype
        integer::c,i

        if (rsign(1:14)=='bolsig:ELASTIC') then
            rtype=2
            nspec=1
            firstprod=0
            c=16
        elseif (rsign(1:7)=='bolsig:') then
            rtype=1
            c=8
        else
            rtype=0
            c=1
        endif
        i=1
        nspec=1
        species(1)=''
        do while (c.le.reactions_length .and. rsign(c:c).ne.' ')
            if (rsign(c:c)=='+') then
                if (rsign(c-1:c-1).ne.'^') then
                    nspec=nspec+1
                    species(nspec)=''
                    i=1
                endif
            elseif (rsign(c:c+1)=='->' .or. rsign(c:c+1)=='=>') then
                nspec=nspec+1
                species(nspec)=''
                i=1
                firstprod=nspec
                c=c+1
            else
                if (i.le.species_length) then
                    species(nspec)(i:i)=rsign(c:c)
                    i=i+1
                endif
            endif
            c=c+1
        enddo
    end subroutine


    function getCO2ind(spname)
        implicit none
        character(species_length)::spname
        integer::c,i,getCO2ind,b1,b2,lint,num
        if (spname(1:3)=='CO2' .and. spname(4:4)/='^') then
            if (spname(4:4)==' ') then
                getCO2ind=1
                return
            endif
            do i=1,species_length
                if (spname(i:i)=='(') b1=i+1
                if (spname(i:i)==')') b2=i-1
            enddo
            if (spname(b1:b1)=='E') then
                lint=b2-b1
                if (lint==1) then
                    read(spname(b1+1:b2),'(I1)') num
                else
                    read(spname(b1+1:b2),'(I2)') num
                endif
                getCO2ind=26+num
                return
            endif

            if (spname(b1:b1)=='V') then
                if ((ichar(spname(b1+1:b1+1)).ge.65) .and. (ichar(spname(b1+1:b1+1)).le.90)) then
                    getCO2ind=1+ichar(spname(b1+1:b1+1))-64
                    return
                endif
                lint=b2-b1
                if (lint==1) then
                    read(spname(b1+1:b2),'(I1)') num
                else
                    read(spname(b1+1:b2),'(I2)') num
                endif
                getCO2ind=5+num
                return
            endif
        endif
        getCO2ind=-1  
    end function


    function getCOind(spname)
        implicit none
        character(species_length)::spname
        integer::c,i,getCOind,b1,b2,lint,num

        if (spname(1:2)=='CO') then 
            if (spname(3:3)==' ') then
                getCOind=1
                return
            endif
            if (spname(3:3)=='(') then
                do i=1,species_length
                    if (spname(i:i)=='(') b1=i+1
                    if (spname(i:i)==')') b2=i-1
                enddo
                lint=b2-b1
                if (lint==1) then
                    read(spname(b1+1:b2),'(I1)') num
                else
                    read(spname(b1+1:b2),'(I2)') num
                endif
                if (spname(b1:b1)=='E') then
                    getCOind=11+num
                    return
                endif
                if (spname(b1:b1)=='V') then
                    getCOind=1+num
                    return
                endif
            endif
        endif
        getCOind=-1
    end function

end module
