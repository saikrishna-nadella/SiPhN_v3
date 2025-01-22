# -*- coding: utf-8 -*-
"""
Created on Sun June 07 16:45:34 2020

@author: saikrishna Nadella
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

def LSA_SHF_SS_func(inpdata,Grm,Stmoc,uiHandle):
    ErrSig = 1
    
    (L1,L2,L3,L4,Lc,Lh,X1,X2,X3,X4,X5,X6,dummy,dummy)=inpdata[0]
    (H,D,tw,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alphah,alphac,dummy,dummy,dummy)=inpdata[1]
    (hoh,hohl,hocl,K,Tsin,Tamb,T_UR,Re_UR,Fsp,TempVal,Tav_g,Ress_g,ModFlag,dummy)=inpdata[2]
    (ndro,ndcp,ndk,Prf,ReTol,TavTol,HXEffect,theta_s,kf,dummy,dummy,dummy,dummy,dummy)=inpdata[3]

    #constants
    pi = np.pi

    # angle conversion to radians
    alpha1 = np.radians(alpha1)
    alpha2 = np.radians(alpha2)
    alpha3 = np.radians(alpha3)
    alpha4 = np.radians(alpha4)
    alpha5 = np.radians(alpha5)
    alpha6 = np.radians(alpha6)
    alphah = np.radians(alphah)
    alphac = np.radians(alphac)

    #geometric parameters
    D1 = D+tw
    D2 = D+2*tw
    Lhl = L1+X1+X2+X3+L2
    Lcl = L3+X4+X5+X6+L4
    Lt = Lh+Lhl+Lc+Lcl
    phi = Lt/H
    A = pi*0.25*D**2
    A1 = 0.25*pi*(D1**2-D**2)
    A2 = pi*0.25*(D2**2-D1**2)

    # SS code
    count_iter = 1
    if ModFlag ==1 or ModFlag ==2:
        ndRw = np.log(1+2*tw/D)/2
    
    Ress_o = 0
    Ress_n=Ress_g
    po = 0.316
    pn=64
    b=1

    NuCorr = uiHandle.PUM_HT.currentIndex()
    fCorr = uiHandle.PUM_friction.currentIndex()
    
    while abs(Ress_o-Ress_n) > ReTol or abs(po-pn)>0.00001:
        Ress_o = Ress_n*Re_UR+Ress_o*(1-Re_UR)
        po = pn
        if ModFlag ==1 or ModFlag ==2:
            Xi = 4*(Lt/D)*ndk/(Ress_o*Prf*ndRw)
        else:
            Xi = 0

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
            Nuih = 0
            Nuihl = 0
            Nuic = 0
            Nuicl = 0
            Stmih = 0
            Stmihl = 0
            Stmic = 0
            Stmicl = 0
        # overall Stm
        if ModFlag ==1:
            Stmoh = 4*(Lt/D)*(hoh*D2/kf)/(Ress_o*Prf)
            Stmohl = 4*(Lt/D)*(hohl*D2/kf)/(Ress_o*Prf)
            Stmocl = 4*(Lt/D)*(hocl*D2/kf)/(Ress_o*Prf)
            Chi_h = 1./(1/Stmih+1/Xi+1./Stmoh)
            Chi_hl = 1./(1/Stmihl+1/Xi + 1./Stmohl)
            Chi_cl = 1./(1/Stmicl+1/Xi+1./Stmocl)
        else:
            Stmoh = 0
            Stmohl = 0
            Stmocl = 0
            Chi_h = 0
            Chi_hl = 0
            Chi_cl = 0

        if ModFlag ==1 or ModFlag ==2:
            Chi_c = 1/(1/Stmic+1/Xi+1/Stmoc)
        elif ModFlag == 3:
            Chi_c = 1/(1/Stmic+1/Stmoc)
        else:
            Chi_c = Stmoc #This should be noted
                
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

        #Re_ss prediction
        Ress_n = np.power((D*(2*Grm*Iss-K*np.power(Ress_o,3))/(po*Lt)),(1/(3-b)))
        #friction parameter correction
        if fCorr == 0:
            psi = 1/(1+np.exp((Ress_n-2530)/120))
            pn = (np.power(64,psi))*(np.power(0.316,(1-psi)))
            b = psi+(1-psi)*0.25
        elif fCorr == 1:
            (pn, b) = f_udc(Ress_n)
        
        count_iter = count_iter+1

    re_ss = Ress_n
    p=po
    if np.isnan(re_ss):
        uiHandle.TE_TEXTOUT.setPlainText('Error: Reynolds number became unrealistic (probably, complex number)')
        ErrSig = 0
        return ErrSig,[],[]
    #preparding SS2CEdata 
    SS2CEdata = ((Lhl, Lcl,Lt,A,A1,A2,Xi,Stmoh,Stmohl,Stmocl,Stmih,Stmihl,Stmic,Stmicl,Nuih, Nuihl, Nuic, Nuicl, re_ss),(Chi_c,Chi_cl,Chi_h,Chi_hl,theta_hf_ss,theta_hlf_ss,theta_cf_ss,theta_clf_ss,p,b,NuCorr))
    return ErrSig,SS2CEdata,[]
