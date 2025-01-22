# -*- coding: utf-8 -*-
"""
Created on Sun May 31 14:23:19 2020

@author: saikr
"""
import numpy as np
from Nu_main import Nu_main_func
import os
import sys
thispath = os.getcwd()
UDpath = thispath+'\\UserDefined'
sys.path.append(UDpath)
from user_f import f_udc
from user_Nu import Nu_udc


def SSA_SHF_SS_func(inpdata,Grm,hoc,uiHandle,DIM):
    
    (L1,L2,L3,L4,Lc,Lh,X1,X2,X3,X4,X5,X6,dummy,dummy)=inpdata[0]
    (H,D,tw,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alphah,alphac,dummy,dummy,dummy)=inpdata[1]
    (hoh,hohl,hocl,K,Ts,Tamb,T_UR,Re_UR,Fsp,TempVal,Tav_g,Ress_g,ModFlag,dummy)=inpdata[2]
    (rhow,cpw,kw,Prf,ReTol,TavTol,HXEffect,theta_s,kf,betaf,rhof,cpf,dummy,dummy)=inpdata[3]
    #constants
    pi = np.pi
    try:
        g = float(uiHandle.LE_g.text())
    except (ArithmeticError, TypeError, ValueError):
        g = 9.81
    # angle conversion to radians
    alpha1 = np.radians(alpha1)
    alpha2 = np.radians(alpha2)
    alpha3 = np.radians(alpha3)
    alpha4 = np.radians(alpha4)
    alpha5 = np.radians(alpha5)
    alpha6 = np.radians(alpha6)
    alphah = np.radians(alphah)
    alphac = np.radians(alphac)

    print('Prandtl number: {:g}'.format(Prf))
    
    #model selection being printed to console
    if ModFlag==1:
        if HXEffect==0:
            print('Warning!: Usage of unverified combination of model parameters! \n \t You are using BM+WTI+HX+HL model without HXEffect.')
        else:
            print('Model: BM+WTI+HX+HL')

        if (hoh < 1e-3) or (hohl < 1e-3) or (hocl <1e-3):
            print('Warning!: Usage of unverified combination of model parameters! \n \t You are using BM+WTI+HX+HL model with insignificant heat loss coefficients. \n It is recommended to use BM+WTI+HX model instead.')
        else:
            print('Model: BM+WTI+HX+HL')
    elif ModFlag==2:
        if HXEffect==0:
            print('Model: BM+WTI')
        else:
            print('Model: BM+WTI+HX')
    elif ModFlag==3:
        print('Model: BM*')
    elif ModFlag==4:
        if K==0:
            print('Model: EM')
        else:
            print('Model: BM')
        print('It should be noted that hlc input will be taken as Ui in current model!')
    else:
        print('Select proper model')
        
    #Geometric parameters

  #  D1 = D+tw
    D2 = D+2*tw
    Lhl = L1+X1+X2+X3+L2
    Lcl = L3+X4+X5+X6+L4
    Lt = Lh+Lhl+Lc+Lcl
    phi = Lt/H
    A=pi*0.25*D**2
  #  A1 = 0.25*pi*(D1**2-D**2)
  #  A2 = pi*0.25*(D2**2-D1**2)
    
    muf = Prf*kf/cpf
    if DIM==1:
        power = Grm
        Grm = (rhof**2*g*betaf*power*H*D**3)/(A*cpf*muf**3)
        theta_s_0 = theta_s
   # ndro = rhow/rhof
   # ndcp = cpw/cpf
    ndk = kw/kf
    
    ##SS code
    if ModFlag ==1 or ModFlag ==2:
        ndRw = np.log(1+2*tw/D)/2
   # Tav_o = Tav_g
   # Tav_n = Tav_o
    Ress_o = 0
    Ress_n=Ress_g
    po = 0.316
    pn=64
    b=1
    count_iter = 1

    NuCorr = uiHandle.PUM_HT.currentIndex()
    fCorr =  uiHandle.PUM_friction.currentIndex()
    
    while abs(Ress_o-Ress_n) > ReTol or abs(po-pn)>0.00001:
        Ress_o = Ress_n*Re_UR+Ress_o*(1-Re_UR)
        po = pn
        if ModFlag ==1 or ModFlag ==2:
            Xi = 4*(Lt/D)*ndk/(Ress_o*Prf*ndRw)
        else:
            Xi = np.nan
        if ModFlag ==1 or ModFlag ==2 or ModFlag ==3:
            # Nui
            if NuCorr == 0:
                Nuih = Nu_main_func(Lh,D,Ress_o,Prf,1)
                Nuihl = Nu_main_func(Lhl,D,Ress_o,Prf,0)
                Nuic = Nu_main_func(Lc,D,Ress_o,Prf,2)
                Nuicl = Nu_main_func(Lcl,D,Ress_o,Prf,0)
            elif NuCorr == 1:
                Nuih = Nu_udc(Lh,D,Ress_o,Prf,1)
                Nuihl = Nu_udc(Lhl,D,Ress_o,Prf,0)
                Nuic = Nu_udc(Lc,D,Ress_o,Prf,2)
                Nuicl = Nu_udc(Lcl,D,Ress_o,Prf,0)
            #Stmi
            Stmih = 4*(Lt/D)*Nuih/(Ress_o*Prf)
            Stmihl =4*(Lt/D)*Nuihl/(Ress_o*Prf)
            Stmic = 4*(Lt/D)*Nuic/(Ress_o*Prf)
            Stmicl = 4*(Lt/D)*Nuicl/(Ress_o*Prf)
        else:
            # Nui
            Nuih = np.nan
            Nuihl = np.nan
            Nuic = np.nan
            Nuicl = np.nan
            # Stmi
            Stmih = np.nan;
            Stmihl =np.nan;
            Stmic = np.nan;
            Stmicl = np.nan;

        # overall Stm
        if ModFlag ==1:
            if DIM==1:
                Stmoh = 4*(Lt/D)*(hoh*D2/kf)/(Ress_o*Prf)
                Stmohl = 4*(Lt/D)*(hohl*D2/kf)/(Ress_o*Prf)
                Stmocl = 4*(Lt/D)*(hocl*D2/kf)/(Ress_o*Prf)
            else:
                Stmoh = hoh
                Stmohl = hohl
                Stmocl = hocl

            Chi_h = 1/(1/Stmih+1/Xi+1/Stmoh)
            Chi_hl = 1/(1/Stmihl+1/Xi + 1/Stmohl)
            Chi_cl = 1/(1/Stmicl+1/Xi+1/Stmocl)
        else:
            Stmoh = 0
            Stmohl = 0
            Stmocl = 0
            Chi_h = 0
            Chi_hl = 0
            Chi_cl = 0


        if DIM==1:
            if ModFlag == 4 or ModFlag == 3:
                Stmoc = 4*(Lt/D)*(hoc*D/kf)/(Ress_o*Prf)
            else:
                Stmoc = 4*(Lt/D)*(hoc*D2/kf)/(Ress_o*Prf)
        else:
            Stmoc = hoc

        if ModFlag ==1 or ModFlag ==2:
            Chi_c = 1/(1/Stmic+1/Xi+1/Stmoc)
        elif ModFlag == 3:
            Chi_c = 1/(1/Stmic+1/Stmoc)
        else:
            Chi_c = Stmoc #This should be noted

        if DIM==1:
            w = Ress_o*pi*D*muf/4
            dTh = power/(w*cpf)
            theta_s = theta_s_0/dTh
        #Temp distribution and Iss calculation
        if ModFlag ==1:
            theta_clf_ss = (theta_s*np.exp(-Chi_cl*Lcl/Lt)*(1-np.exp(-Chi_c*Lc/Lt))+(Lt/(Lh*Stmoh))*(1-np.exp(-Chi_h*Lh/Lt))*np.exp(-1*(Chi_hl*Lhl+Chi_c*Lc+Chi_cl*Lcl)/Lt))/(1-np.exp(-1*(Chi_hl*Lhl+Chi_c*Lc+Chi_cl*Lcl+Chi_h*Lh)/Lt))
            theta_hf_ss = (Lt/(Lh*Stmoh))*(1-np.exp(-Chi_h*Lh/Lt))+theta_clf_ss*np.exp(-Chi_h*Lh/Lt)
            theta_hlf_ss = theta_hf_ss*np.exp(-Chi_hl*Lhl/Lt)
            theta_cf_ss = theta_s*(1-np.exp(-Chi_c*Lc/Lt))+theta_hlf_ss*np.exp(-Chi_c*Lc/Lt)
            
            theta_clf_ss_term = theta_clf_ss*((phi/Chi_h)*(1-np.exp(-Chi_h*Lh/Lt))*np.sin(alphah))
            theta_hf_ss_term = theta_hf_ss*(phi/Chi_hl)*((1-np.exp(-Chi_hl*L1/Lt))*np.sin(alphah)+np.sin(alpha1)*(np.exp(-Chi_hl*L1/Lt)-np.exp(-Chi_hl*(L1+X1)/Lt))+np.sin(alpha2)*(np.exp(-Chi_hl*(L1+X1)/Lt)-np.exp(-Chi_hl*(L1+X1+X2)/Lt))+np.sin(alpha3)*(np.exp(-Chi_hl*(X1+X2+L1)/Lt)-np.exp(-Chi_hl*(Lhl-L2)/Lt))+(np.exp(-Chi_hl*(Lhl-L2)/Lt)-np.exp(-Chi_hl*Lhl/Lt))*np.sin(alphac))
            theta_hlf_ss_term = theta_hlf_ss*(phi/Chi_c)*(1-np.exp(-Chi_c*Lc/Lt))*np.sin(alphac)
            theta_cf_ss_term = theta_cf_ss*(phi/Chi_cl)*((1-np.exp(-Chi_cl*L3/Lt))*np.sin(alphac)+np.sin(alpha4)*(np.exp(-Chi_cl*L3/Lt)-np.exp(-Chi_cl*(L3+X4)/Lt))+np.sin(alpha5)*(np.exp(-Chi_cl*(X4+L3)/Lt)-np.exp(-Chi_cl*(X4+X5+L3)/Lt))+np.sin(alpha6)*(np.exp(-Chi_cl*(X4+X5+L3)/Lt)-np.exp(-Chi_cl*(Lcl-L4)/Lt))+(np.exp(-Chi_cl*(Lcl-L4)/Lt)-np.exp(-Chi_cl*Lcl/Lt))*np.sin(alphah))
            ex_ss_term = (Lt/Lh)*((Lh/H+(phi/Chi_h)*(np.exp(-Chi_h*Lh/Lt)-1))/(Stmoh))*np.sin(alphah)
            theta_s_term = theta_s*np.sin(alphac)*(Lc/H+(phi/Chi_c)*(np.exp(-Chi_c*Lc/Lt)-1))
            Iss = theta_clf_ss_term+theta_hf_ss_term+theta_hlf_ss_term+theta_cf_ss_term+ex_ss_term+theta_s_term
        else:
            theta_clf_ss = theta_s+(np.exp(-Chi_c*Lc/Lt))/(1-np.exp(-Chi_c*Lc/Lt))
            theta_hlf_ss = theta_clf_ss+1
            theta_hf_ss = theta_hlf_ss
            theta_cf_ss = theta_clf_ss

            theta_clf_ss_term = theta_clf_ss*((Lh+L4)*np.sin(alphah)+L3*np.sin(alphac)+X4*np.sin(alpha4)+X5*np.sin(alpha5)+X6*np.sin(alpha6))/H
            theta_hlf_ss_term = theta_hlf_ss*(L1*np.sin(alphah)+X1*np.sin(alpha1)+X2*np.sin(alpha2)+X3*np.sin(alpha3)+(L2+(Lt/Chi_c)*(1-np.exp(-Chi_c*Lc/Lt)))*np.sin(alphac))/H
            ex_ss_term = Lh*np.sin(alphah)/(2*H)
            theta_s_term = theta_s*np.sin(alphac)*(Lc/H+(phi/Chi_c)*(np.exp(-Chi_c*Lc/Lt)-1))
            Iss = theta_clf_ss_term+theta_hlf_ss_term+ex_ss_term+theta_s_term

        # Re_ss prediction
        Ress_n = np.power((D*(2*Grm*Iss-K*np.power(Ress_o,3))/(po*Lt)),(1/(3-b)))
        if fCorr == 0:
            psi = 1/(1+np.exp((Ress_n-2530)/120))
            pn = (np.power(64,psi))*(np.power(0.316,(1-psi)))
            b = psi+(1-psi)*0.25
        elif fCorr == 1:
            (pn, b) = f_udc(Ress_n)
            
        count_iter = count_iter+1
        if count_iter >10000:
            if DIM==1:
                uiHandle.TE_DIMout.append('Maximum iterations over')
            else:
                uiHandle.TE_NDout.append('Maximum iterations over')
            return
    
    re_ss = Ress_n
    #p=po
    w_ss = re_ss*pi*D*muf/4
    Q = Grm*A*cpf*muf**3/(rhof**2*betaf*g*D**3*H)  #heat input in Watts
    dTh_ss = Q/(w_ss*cpf)
    Th_ss = theta_hf_ss*dTh_ss+Tamb
    Thl_ss = theta_hlf_ss*dTh_ss+Tamb
    Tc_ss = theta_cf_ss*dTh_ss+Tamb
    Tcl_ss = theta_clf_ss*dTh_ss+Tamb
    dTc_ss = Thl_ss-Tc_ss
    hih = Nuih*kf/D
    hihl = Nuihl*kf/D
    hic = Nuic*kf/D
    hicl = Nuicl*kf/D
    hlprcnt = 100*(dTh_ss-dTc_ss)/dTh_ss
    hloss = w_ss*cpf*(dTh_ss-dTc_ss)
    hoh = Stmoh*re_ss*Prf*D*kf/(4*Lt*D2)
    hohl = Stmohl*re_ss*Prf*D*kf/(4*Lt*D2)
    if (ModFlag == 4 or ModFlag == 3):
        hoc = Stmoc*re_ss*Prf*D*kf/(4*Lt*D)
    else:
        hoc = Stmoc*re_ss*Prf*D*kf/(4*Lt*D2)

    hocl = Stmocl*re_ss*Prf*D*kf/(4*Lt*D2)
    #outputs
    OString_DIM = ("""STEADY STATE RESULTS:
            Note - Tamb is taken as {:g} degC

            mass flow rate   {:g} kg/s
            Power            {:g} W
            Thi              {:g} degC
            Tho              {:g} degC
            Tci              {:g} degC
            Tco              {:g} degC
            Heat loss        {:g} ({:g} %)
            hih              {:g} W/m2/K
            hihl             {:g} W/m2/K
            hic              {:g} W/m2/K
            hicl             {:g} W/m2/K
            hoh              {:g} W/m2/K
            hohl             {:g} W/m2/K
            hoc              {:g} W/m2/K
            hocl             {:g}  W/m2/K
                    """.format(Tamb,w_ss,Q,Tcl_ss,Th_ss,Thl_ss,Tc_ss,hloss,hlprcnt,hih,hihl,hic,hicl,hoh,hohl,hoc,hocl))
    OString_ND = (""" STEADY STATE RESULTS:
            Reynolds number   {:g}
            Grm               {:g}
            Stmih             {:g}
            Stmihl            {:g}
            Stmic             {:g}
            Stmicl            {:g}
            Xi                {:g}
            Nuih              {:g}
            Nuihl             {:g}
            Nuic              {:g}
            Nuicl             {:g}
            Stmoh             {:g}
            Stmohl            {:g}
            Stmoc             {:g}
            Stmocl            {:g}
            theta_s           {:g}
                  """.format(re_ss,Grm,Stmih,Stmihl,Stmic,Stmicl,Xi,Nuih,Nuihl,Nuic,Nuicl,Stmoh,Stmohl,Stmoc,Stmocl,theta_s))
    if DIM==1:
        uiHandle.TE_DIMout.append(OString_DIM)
        uiHandle.TE_NDout.setPlainText(OString_ND)
    else:
        uiHandle.TE_DIMout.setPlainText(OString_DIM)
        uiHandle.TE_NDout.append(OString_ND)
    return None
    
