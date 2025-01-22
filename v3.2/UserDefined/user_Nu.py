# -*- coding: utf-8 -*-

def Nu_udc(L,D,Re,Pr,hflag):
    """
      User defined correlation for Nusselt number for heat transfer between 
    working fluid and pipe inside surface

      Data given in terms of arguments: 
        L = length of section (h,hl,c,cl)
        D = inner diameter of section
        Re = Reynolds number of flow in section
        Pr = Prandtl number of working fluid
        hflag = 1 (heater), 2(cooler), 0(piping)

      function returns Nui = internal Nusselt number
    
    """
    Nui = 3
    # Add your correlation here -------------------------------------------
    

    # --------- User correlation ends -------------------------------------
    return Nui

