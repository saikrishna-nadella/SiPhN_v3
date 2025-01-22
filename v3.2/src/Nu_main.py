# -*- coding: utf-8 -*-

import numpy as np

# hflag = 1 (heater), 2(cooler), 0(adiabatic pipe)
def Nu_main_func(L,D,Re,Pr,hflag):
    #Hausen correlation
    #nuh  = 3.66+(0.0668*(D/L)*Re *Pr)/(1+0.04*((D/L)*Re *Pr)^0.67);
    #Skupinski correlation
    #nusk = 4.82+0.0185*(Re*Pr)^(0.827); % q = const. Pr(.003-0.05) Re(3600-9e+5) Pe(100-1e+4)
    Pe = Re*Pr;
    #Cheng and Tak correlation
    if Pe<1000:
        nuct = 4.5+0.018*(Pe)**0.8
    elif Pe>1000 and Pe<2000:
        nuct = 5.4-0.0009*Pe+0.018*Pe**0.8
    else:
        nuct = 3.6+0.018*Pe**0.8;
    
    #Shah and London correlation
    if hflag ==1 or hflag == 2:
        #Resl = L/(D*0.03*Pr);
        xstar = L/(D*Re*Pr)
        if xstar<=0.025:
            nusl = 1.953*xstar**(-1/3)
        elif xstar >=0.035:
            nusl = 4.364+0.0722/xstar
        else:
            nusl=(1.953*xstar**(-1/3)*(0.035-xstar)/0.01+(1/0.01)*(xstar-0.025)*(4.364+0.0722/xstar))
    else:
        nusl = 4.364
    #Gnielinski correlation
    f = (0.79*np.log(Re)-1.64)**(-2)
    nug  = ((f/8)*(Re -1000)*Pr)/(1+12.7*np.sqrt(f/8)*(Pr**0.67-1))
    if nug<=4.36:
        nug=4.36
    #Disttus - Boelter correlation
    if hflag == 1:
        n=0.4
    else:
        n=0.3
    
    nud  = 0.023*(Re)**0.8*Pr**n
    #Shah and London correlation -> Gnielinski correlation connector
    si1  = 1/(1+np.exp((Re -2500)/500))
    #(SL & Gnielinski) -> DB correlation connector
    si2  = 1/(1+np.exp((Re -1e+4)/100))
    #(Pr < 0.6 to Pr>0.6 connector
    si3 = 1/(1+np.exp((Pr-0.6)/0.00001))
    #Final relation for internal heat transfer coefficient
    Nui1  = ((nusl)**si1*(nug)**(1-si1))**si2*(nud)**(1-si2)
    Nui = (nuct**si3)*(Nui1**(1-si3))
    return Nui