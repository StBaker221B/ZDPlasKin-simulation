module input
    implicit none
    double precision::timeend,radius,length,pressure,igastemp,voltage,resistance,gasflow,efreq,ieden
    namelist /INPUTDATA/ timeend,radius,length,pressure,igastemp,voltage,resistance,gasflow,efreq,ieden
        
    contains
    
    subroutine set_default
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
    end subroutine
    
end
    