# -*- coding: utf-8 -*-

def ThermoPhysicalProperties(T, flag):
    """
    User can enter the thermophysical property correlations for fluid and wall materials
    
    Data given in terms of arguments :
    T = Temperature of <flag> in degC
    flag = 1 (fluid), 2(wall), 3 (both)
    
    rho = density (kg/m3), 
    cp=specific heat capacity (J/kg.K),
    k=thermal conductivity (W/mK),
    mu=dynamic viscosity (Pa.s), 
    beta = volumetric thermal expansion coefficient (1/K)
    
    'flag' is introduced to reduce computational cost when property dependent
    calculations are performed.Please type in your correlations at suitable 
    location such that the philosophy of introducing 'flag' will not be lost
    
    function returns an array of properties TPP_OUT and
    an output string to be displayed in Gui 'OString2Gui' [optional]
    
    """
    # initializing property values (required)
    rhof = 0; cpf = 0; kf = 0; muf = 0; betaf = 0; rhow = 0; cpw = 0; kw = 0;
    OString2Gui = ''
    
    # ---------------User editable part starts ------------------------------% 
    
    if(flag == 1 or flag == 3):
        # write the program for calculation of working fluid properties here.
        if(T < 260):
            # conditional statement may be used for limiting temperatur range
            OString2Gui = 'Temperature is limited to lower bound 260 degC'
            # Optional messages can be written in "OString2Gui" for displaying
            # in status box
            T = 260
        elif(T > 600):
            OString2Gui = 'Temperature is limited to upper bound 600 degC'
            T = 600
        rhof = -0.64029*T+2090.1809
        if(T<460):
            cpf = 0.15072*T+1455.49915
        else:
            cpf = 0.15072*T+1459.68595
        muf = 0.02271-0.000119999*T+2.28097*0.0000001*T**2-1.47397*0.0000000001*T**3
        kf = 0.44141+1.95194*0.0001*T
        betaf = 0.64029/rhof

    if(flag == 2 or flag == 3):
        # write the program for calculation of wall properties here.
        rhow = 8700
        cpw = 500
        kw = 16
# ------------- user editable part over ----------------------------------%
    if OString2Gui == '':
        OString2Gui='Successful'
    TPP_OUT = (rhof,cpf,kf,muf,betaf,rhow,cpw,kw)
    #TPP_OUT = tuple([round(x,2) for x in TPP_OUT])
    TPP_OUT = tuple(['{:.2e}'.format(x) for x in TPP_OUT])
    return OString2Gui, TPP_OUT

