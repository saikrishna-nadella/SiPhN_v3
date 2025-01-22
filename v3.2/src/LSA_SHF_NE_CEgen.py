# -*- coding: utf-8 -*-
"""
Created on Sun June 07 17:27:58 2020

@author: saikrishna Nadella
"""
import numpy as np
from Nu_main import Nu_main_func
import os
import sys
thispath = os.getcwd()
UDpath = thispath+'\\UserDefined'
sys.path.append(UDpath)
from user_Nu import Nu_udc

def LSA_SHF_NE_CEgen_func(inpdata,SS2CEdata,nI,Stmoc,Grm):
    n = nI*np.complex(0,1)
    # input data assignment
    (L1,L2,L3,L4,Lc,Lh,X1,X2,X3,X4,X5,X6,dummy,dummy)=inpdata[0]
    (H,D,tw,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alphah,alphac,dummy,dummy,dummy)=inpdata[1]
    (hoh,hohl,hocl,K,Tsin,Tamb,T_UR,Re_UR,Fsp,TempVal,Tav_g,Ress_g,ModFlag,dummy)=inpdata[2]
    (ndro,ndcp,ndk,Prf,ReTol,TavTol,HXEffect,theta_s,kf,dummy,dummy,dummy,dummy,dummy)=inpdata[3]

    # SS data assignment
    (Lhl, Lcl,Lt,A,A1,A2,Xi,Stmoh,Stmohl,Stmocl,Stmih,Stmihl,Stmic,Stmicl,Nuih, Nuihl, Nuic, Nuicl, re_ss)=SS2CEdata[0]
    (Chi_c,Chi_cl,Chi_h,Chi_hl,theta_hf_ss,theta_hlf_ss,theta_cf_ss,theta_clf_ss,p,b,NuCorr)=SS2CEdata[1]

    # angle conversion to radians
    alpha1 = np.radians(alpha1)
    alpha2 = np.radians(alpha2)
    alpha3 = np.radians(alpha3)
    alpha4 = np.radians(alpha4)
    alpha5 = np.radians(alpha5)
    alpha6 = np.radians(alpha6)
    alphah = np.radians(alphah)
    alphac = np.radians(alphac)
    
    phi = Lt/H
    
    if ModFlag == 1 or ModFlag ==2:
        Sig1 = (A/A1)*(1/(ndro*ndcp))
        Sig2 = (A/A2)*(1/(ndro*ndcp))

    if ModFlag ==1:
        Mlh = n+Xi*Sig1*((n+Stmoh*Sig2)/(n+Sig2*(Xi+Stmoh)))
        Mlhl = n+Xi*Sig1*((n+Stmohl*Sig2)/(n+Sig2*(Xi+Stmohl)))
        Mlc = n+Xi*Sig1*((n+Stmoc*Sig2)/(n+Sig2*(Xi+Stmoc)))
        Mlcl = n+Xi*Sig1*((n+Stmocl*Sig2)/(n+Sig2*(Xi+Stmocl)))
    elif ModFlag ==2:
        Mlh = n+Xi*Sig1*(n/(n+Sig2*Xi))
        Mlhl = n+Xi*Sig1*(n/(n+Sig2*Xi))
        Mlc = n+Xi*Sig1*((n+Stmoc*Sig2)/(n+Sig2*(Xi+Stmoc)))
        Mlcl = n+Xi*Sig1*(n/(n+Sig2*Xi))
    
    # lam1j calculation
    if ModFlag ==1 or ModFlag == 2:
        lam1h = n+(Mlh*Stmih)/(Mlh+Stmih*Sig1)
        lam1hl = n+(Mlhl*Stmihl)/(Mlhl+Stmihl*Sig1)
        lam1c = n+(Mlc*Stmic)/(Mlc+Stmic*Sig1)
        lam1cl = n+(Mlcl*Stmicl)/(Mlcl+Stmicl*Sig1)
    elif ModFlag == 3 or ModFlag ==4:
        lam1h = n
        lam1hl = n
        lam1c = n+Chi_c
        lam1cl = n
    
    # B calculation
    if HXEffect ==1:
        if re_ss < 500:
            eps=0.01*re_ss
        else:
            eps = 5
        if NuCorr == 0:
            Bh = (re_ss/Nuih)*(Nu_main_func(Lh,D,re_ss+eps,Prf,1)-Nu_main_func(Lh,D,re_ss-eps,Prf,1))/(2*eps)
            Bhl = (re_ss/Nuihl)*(Nu_main_func(Lhl,D,re_ss+eps,Prf,0)-Nu_main_func(Lhl,D,re_ss-eps,Prf,0))/(2*eps);
            Bc = (re_ss/Nuic)*(Nu_main_func(Lc,D,re_ss+eps,Prf,2)-Nu_main_func(Lc,D,re_ss-eps,Prf,2))/(2*eps);
            Bcl = (re_ss/Nuicl)*(Nu_main_func(Lcl,D,re_ss+eps,Prf,0)-Nu_main_func(Lcl,D,re_ss-eps,Prf,0))/(2*eps);
        elif NuCorr == 1:
            Bh = (re_ss/Nuih)*(Nu_udc(Lh,D,re_ss+eps,Prf,1)-Nu_udc(Lh,D,re_ss-eps,Prf,1))/(2*eps)
            Bhl = (re_ss/Nuihl)*(Nu_udc(Lhl,D,re_ss+eps,Prf,0)-Nu_udc(Lhl,D,re_ss-eps,Prf,0))/(2*eps);
            Bc = (re_ss/Nuic)*(Nu_udc(Lc,D,re_ss+eps,Prf,2)-Nu_udc(Lc,D,re_ss-eps,Prf,2))/(2*eps);
            Bcl = (re_ss/Nuicl)*(Nu_udc(Lcl,D,re_ss+eps,Prf,0)-Nu_udc(Lcl,D,re_ss-eps,Prf,0))/(2*eps);
    else:
        Bh,Bhl,Bc,Bcl=0,0,0,0
    
    # lam2j calculation
    if ModFlag ==1:
        lam2h = Chi_h*(theta_clf_ss-Lt/(Lh*Stmoh))*(1-Bh+(Bh*Sig1*Stmih)/(Mlh+Stmih*Sig1))
        lam2hl = Chi_hl*theta_hf_ss*(1-Bhl+(Bhl*Sig1*Stmihl)/(Mlhl+Stmihl*Sig1))
        lam2c = Chi_c*(theta_hlf_ss-theta_s)*(1-Bc+(Bc*Sig1*Stmic)/(Mlc+Stmic*Sig1))
        lam2cl = Chi_cl*theta_cf_ss*(1-Bcl+(Bcl*Sig1*Stmicl)/(Mlcl+Stmicl*Sig1))
    elif ModFlag ==2 and HXEffect ==1:
        lam2h = -1*(Lt/Lh)*(1-Bh+(Bh*Sig1*Stmih)/(Mlh+Stmih*Sig1))
        lam2hl = 0
        lam2c = Chi_c*(theta_hlf_ss-theta_s)*(1-Bc+(Bc*Sig1*Stmic)/(Mlc+Stmic*Sig1))
        lam2cl = 0
    elif (ModFlag ==2 and HXEffect == 0) or ModFlag ==3 or ModFlag == 4:
        lam2h = -1*(Lt/Lh)
        lam2hl = 0
        lam2c = Chi_c*(theta_hlf_ss-theta_s)
        lam2cl = 0
         
    # aj calculation
    a_cl = lam2cl/(lam1cl-Chi_cl)
    a_c = lam2c/(lam1c-Chi_c)
    a_hl = lam2hl/(lam1hl - Chi_hl)
    a_h = lam2h/(lam1h-Chi_h)
    
    # point perturbed temperatures
    if ModFlag ==1:
        theta_clf_term = (a_cl*(np.exp(-Chi_cl*Lcl/Lt)-np.exp(-lam1cl*Lcl/Lt))+a_c*np.exp(-lam1cl*Lcl/Lt)*(np.exp(-Chi_c*Lc/Lt)-np.exp(-lam1c*Lc/Lt))+a_hl*np.exp(-1*(lam1cl*Lcl+lam1c*Lc)/Lt)*(np.exp(-Chi_hl*Lhl/Lt)-np.exp(-lam1hl*Lhl/Lt))+a_h*np.exp(-1*(lam1cl*Lcl+lam1c*Lc+lam1hl*Lhl)/Lt)*(np.exp(-Chi_h*Lh/Lt)-np.exp(-lam1h*Lh/Lt)))/(1-np.exp(-1*(lam1cl*Lcl+lam1c*Lc+lam1h*Lh+lam1hl*Lhl)/Lt))
        theta_hf_term= theta_clf_term*np.exp(-lam1h*Lh/Lt)+a_h*(np.exp(-Chi_h*Lh/Lt)-np.exp(-lam1h*Lh/Lt))
        theta_hlf_term =theta_hf_term*np.exp(-lam1hl*Lhl/Lt)+a_hl*(np.exp(-Chi_hl*Lhl/Lt)-np.exp(-lam1hl*Lhl/Lt))
        theta_cf_term = theta_hlf_term*np.exp(-lam1c*Lc/Lt)+a_c*(np.exp(-Chi_c*Lc/Lt)-np.exp(-lam1c*Lc/Lt))
    elif ModFlag ==2 or ModFlag ==3 or ModFlag ==4:
        theta_clf_term = (a_c*np.exp(-lam1cl*Lcl/Lt)*(np.exp(-Chi_c*Lc/Lt)-np.exp(-lam1c*Lc/Lt))+a_h*np.exp(-1*(lam1cl*Lcl+lam1c*Lc+lam1hl*Lhl)/Lt)*(1-np.exp(-lam1h*Lh/Lt)))/(1-np.exp(-1*(lam1cl*Lcl+lam1c*Lc+lam1h*Lh+lam1hl*Lhl)/Lt))
        theta_hf_term= theta_clf_term*np.exp(-lam1h*Lh/Lt)+a_h*(1-np.exp(-lam1h*Lh/Lt))
        theta_hlf_term =theta_hf_term*np.exp(-lam1hl*Lhl/Lt)
        theta_cf_term = theta_hlf_term*np.exp(-lam1c*Lc/Lt)+a_c*(np.exp(-Chi_c*Lc/Lt)-np.exp(-lam1c*Lc/Lt))

    # perturbed temperature integral
    if ModFlag ==1:
        # old implementation
        '''
        I_w_term1 = theta_clf_term*(phi/lam1h)*(1-np.exp(-lam1h*Lh/Lt))*np.sin(alphah)+a_h*np.sin(alphah)*((phi/Chi_h)*(1-np.exp(-1*Chi_h*Lh/Lt))+(phi/lam1h)*(np.exp(-lam1h*Lh/Lt)-1))
        I_w_term2a = theta_hf_term*(phi/lam1hl)*(1-np.exp(-lam1hl*L1/Lt))*np.sin(alphah)+a_hl*np.sin(alphah)*((phi/Chi_hl)*(1-np.exp(-Chi_hl*L1/Lt))+(phi/lam1hl)*(np.exp(-lam1hl*L1/Lt)-1))
        I_w_term2ba = theta_hf_term*(phi/lam1hl)*(np.exp(-lam1hl*L1/Lt)-np.exp(-lam1hl*(L1+X1)/Lt))*np.sin(alpha1)+a_hl*np.sin(alpha1)*((phi/Chi_hl)*(np.exp(-Chi_hl*L1/Lt)-np.exp(-Chi_hl*(L1+X1)/Lt))+(phi/lam1hl)*(np.exp(-lam1hl*(L1+X1)/Lt)-np.exp(-lam1hl*L1/Lt)))
        I_w_term2bb = theta_hf_term*(phi/lam1hl)*(np.exp(-lam1hl*(L1+X1)/Lt)-np.exp(-lam1hl*(L1+X1+X2)/Lt))*np.sin(alpha2)+a_hl*np.sin(alpha2)*((phi/Chi_hl)*(np.exp(-Chi_hl*(L1+X1)/Lt)-np.exp(-Chi_hl*(L1+X1+X2)/Lt))+(phi/lam1hl)*(np.exp(-lam1hl*(L1+X1+X2)/Lt)-np.exp(-lam1hl*(L1+X1)/Lt)))
        I_w_term2bc = theta_hf_term*(phi/lam1hl)*(np.exp(-lam1hl*(L1+X1+X2)/Lt)-np.exp(-lam1hl*(Lhl-L2)/Lt))*np.sin(alpha3)+a_hl*np.sin(alpha3)*((phi/Chi_hl)*(np.exp(-Chi_hl*(L1+X1+X2)/Lt)-np.exp(-Chi_hl*(Lhl-L2)/Lt))+(phi/lam1hl)*(np.exp(-lam1hl*(Lhl-L2)/Lt)-np.exp(-lam1hl*(L1+X1+X2)/Lt)))
        I_w_term2c = theta_hf_term*(phi/lam1hl)*(np.exp(-lam1hl*(Lhl-L2)/Lt)-np.exp(-lam1hl*Lhl/Lt))*np.sin(alphac)+a_hl*((phi/Chi_hl)*(np.exp(-Chi_hl*(Lhl-L2)/Lt)-np.exp(-Chi_hl*Lhl/Lt))+(phi/lam1hl)*(np.exp(-lam1hl*Lhl/Lt)-np.exp(-lam1hl*(Lhl-L2)/Lt)))*np.sin(alphac)
        I_w_term3 = theta_hlf_term*np.sin(alphac)*(phi/lam1c)*(1-np.exp(-lam1c*Lc/Lt))+a_c*np.sin(alphac)*((phi/Chi_c)*(1-np.exp(-Chi_c*Lc/Lt))+(phi/lam1c)*(np.exp(-lam1c*Lc/Lt)-1))
        I_w_term4a = theta_cf_term*(phi/lam1cl)*(1-np.exp(-lam1cl*L3/Lt))*np.sin(alphac)+a_cl*np.sin(alphac)*((phi/Chi_cl)*(1-np.exp(-Chi_cl*L3/Lt))+(phi/lam1cl)*(np.exp(-lam1cl*L3/Lt)-1))
        I_w_term4ba = theta_cf_term*(phi/lam1cl)*(np.exp(-lam1cl*L3/Lt)-np.exp(-lam1cl*(L3+X4)/Lt))*np.sin(alpha4)-a_cl*np.sin(alpha4)*((phi/Chi_cl)*(np.exp(-Chi_cl*(L3+X4)/Lt)-np.exp(-Chi_cl*L3/Lt))+(phi/lam1cl)*(np.exp(-lam1cl*L3/Lt)-np.exp(-lam1cl*(L3+X4)/Lt)))
        I_w_term4bb = theta_cf_term*(phi/lam1cl)*(np.exp(-lam1cl*(L3+X4)/Lt)-np.exp(-lam1cl*(L3+X4+X5)/Lt))*np.sin(alpha5)-a_cl*np.sin(alpha4)*((phi/Chi_cl)*(np.exp(-Chi_cl*(L3+X4+X5)/Lt)-np.exp(-Chi_cl*(L3+X4)/Lt))+(phi/lam1cl)*(np.exp(-lam1cl*(L3+X4)/Lt)-np.exp(-lam1cl*(L3+X4+X5)/Lt)))
        I_w_term4bc = theta_cf_term*(phi/lam1cl)*(np.exp(-lam1cl*(L3+X4+X5)/Lt)-np.exp(-lam1cl*(Lcl-L4)/Lt))*np.sin(alpha6)-a_cl*np.sin(alpha4)*((phi/Chi_cl)*(np.exp(-Chi_cl*(Lcl-L4)/Lt)-np.exp(-Chi_cl*(L3+X4+X5)/Lt))+(phi/lam1cl)*(np.exp(-lam1cl*(L3+X4+X5)/Lt)-np.exp(-lam1cl*(Lcl-L4)/Lt)))
        I_w_term4c = theta_cf_term*np.sin(alphah)*(phi/lam1cl)*(np.exp(-lam1cl*(Lcl-L4)/Lt)-np.exp(-lam1cl*Lcl/Lt))+a_cl*np.sin(alphah)*((phi/Chi_cl)*(np.exp(-Chi_cl*(Lcl-L4)/Lt)-np.exp(-Chi_cl*Lcl/Lt))+(phi/lam1cl)*(np.exp(-lam1cl*Lcl/Lt)-np.exp(-lam1cl*(Lcl-L4)/Lt)))
        '''
        # new optimized implementation
        I_w_term1 = np.sin(alphah)*((phi/lam1h)*(theta_clf_term-a_h)*(1-np.exp(-lam1h*Lh/Lt))+a_h*(phi/Chi_h)*(1-np.exp(-1*Chi_h*Lh/Lt)))
        I_w_term2a = (theta_hf_term-a_hl)*(phi/lam1hl)*(1-np.exp(-lam1hl*L1/Lt))*np.sin(alphah)+a_hl*np.sin(alphah)*(phi/Chi_hl)*(1-np.exp(-Chi_hl*L1/Lt))
        I_w_term2ba = (theta_hf_term-a_hl)*(phi/lam1hl)*(np.exp(-lam1hl*L1/Lt)-np.exp(-lam1hl*(L1+X1)/Lt))*np.sin(alpha1)+a_hl*np.sin(alpha1)*(phi/Chi_hl)*(np.exp(-Chi_hl*L1/Lt)-np.exp(-Chi_hl*(L1+X1)/Lt))
        I_w_term2bb = (theta_hf_term-a_hl)*(phi/lam1hl)*(np.exp(-lam1hl*(L1+X1)/Lt)-np.exp(-lam1hl*(L1+X1+X2)/Lt))*np.sin(alpha2)+a_hl*np.sin(alpha2)*(phi/Chi_hl)*(np.exp(-Chi_hl*(L1+X1)/Lt)-np.exp(-Chi_hl*(L1+X1+X2)/Lt))
        I_w_term2bc = (theta_hf_term-a_hl)*(phi/lam1hl)*(np.exp(-lam1hl*(L1+X1+X2)/Lt)-np.exp(-lam1hl*(Lhl-L2)/Lt))*np.sin(alpha3)+a_hl*np.sin(alpha3)*(phi/Chi_hl)*(np.exp(-Chi_hl*(L1+X1+X2)/Lt)-np.exp(-Chi_hl*(Lhl-L2)/Lt))
        I_w_term2c = (theta_hf_term-a_hl)*(phi/lam1hl)*(np.exp(-lam1hl*(Lhl-L2)/Lt)-np.exp(-lam1hl*Lhl/Lt))*np.sin(alphac)+a_hl*(phi/Chi_hl)*(np.exp(-Chi_hl*(Lhl-L2)/Lt)-np.exp(-Chi_hl*Lhl/Lt))
        I_w_term3 = (theta_hlf_term-a_c)*np.sin(alphac)*(phi/lam1c)*(1-np.exp(-lam1c*Lc/Lt))+a_c*np.sin(alphac)*(phi/Chi_c)*(1-np.exp(-Chi_c*Lc/Lt))
        I_w_term4a = (theta_cf_term-a_cl)*(phi/lam1cl)*(1-np.exp(-lam1cl*L3/Lt))*np.sin(alphac)+a_cl*np.sin(alphac)*(phi/Chi_cl)*(1-np.exp(-Chi_cl*L3/Lt))
        I_w_term4ba = (theta_cf_term-a_cl)*(phi/lam1cl)*(np.exp(-lam1cl*L3/Lt)-np.exp(-lam1cl*(L3+X4)/Lt))*np.sin(alpha4)-a_cl*np.sin(alpha4)*(phi/Chi_cl)*(np.exp(-Chi_cl*(L3+X4)/Lt)-np.exp(-Chi_cl*L3/Lt))
        I_w_term4bb = (theta_cf_term-a_cl)*(phi/lam1cl)*(np.exp(-lam1cl*(L3+X4)/Lt)-np.exp(-lam1cl*(L3+X4+X5)/Lt))*np.sin(alpha5)-a_cl*np.sin(alpha4)*(phi/Chi_cl)*(np.exp(-Chi_cl*(L3+X4+X5)/Lt)-np.exp(-Chi_cl*(L3+X4)/Lt))
        I_w_term4bc = (theta_cf_term-a_cl)*(phi/lam1cl)*(np.exp(-lam1cl*(L3+X4+X5)/Lt)-np.exp(-lam1cl*(Lcl-L4)/Lt))*np.sin(alpha6)-a_cl*np.sin(alpha4)*(phi/Chi_cl)*(np.exp(-Chi_cl*(Lcl-L4)/Lt)-np.exp(-Chi_cl*(L3+X4+X5)/Lt))
        I_w_term4c = (theta_cf_term-a_cl)*np.sin(alphah)*(phi/lam1cl)*(np.exp(-lam1cl*(Lcl-L4)/Lt)-np.exp(-lam1cl*Lcl/Lt))+a_cl*np.sin(alphah)*(phi/Chi_cl)*(np.exp(-Chi_cl*(Lcl-L4)/Lt)-np.exp(-Chi_cl*Lcl/Lt))
        
        I_w_term =I_w_term1+I_w_term2a+I_w_term2ba+I_w_term2bb+I_w_term2bc+I_w_term2c+I_w_term3+I_w_term4a+I_w_term4ba+I_w_term4bb+I_w_term4bc+I_w_term4c
    elif ModFlag ==2 or ModFlag ==3 or ModFlag ==4:
        I_w_term1 = np.sin(alphah)*((phi/lam1h)*(theta_clf_term-a_h)*(1-np.exp(-lam1h*Lh/Lt))+a_h*Lh/H)
        I_w_term2a = theta_hf_term*(phi/lam1hl)*(1-np.exp(-lam1hl*L1/Lt))*np.sin(alphah)
        I_w_term2ba = theta_hf_term*(phi/lam1hl)*(np.exp(-lam1hl*L1/Lt)-np.exp(-lam1hl*(L1+X1)/Lt))*np.sin(alpha1)
        I_w_term2bb = theta_hf_term*(phi/lam1hl)*(np.exp(-lam1hl*(L1+X1)/Lt)-np.exp(-lam1hl*(L1+X1+X2)/Lt))*np.sin(alpha2)
        I_w_term2bc = theta_hf_term*(phi/lam1hl)*(np.exp(-lam1hl*(L1+X1+X2)/Lt)-np.exp(-lam1hl*(Lhl-L2)/Lt))*np.sin(alpha3)
        I_w_term2c = theta_hf_term*(phi/lam1hl)*(np.exp(-lam1hl*(Lhl-L2)/Lt)-np.exp(-lam1hl*Lhl/Lt))*np.sin(alphac)
        I_w_term3 = (theta_hlf_term-a_c)*np.sin(alphac)*(phi/lam1c)*(1-np.exp(-lam1c*Lc/Lt))+a_c*np.sin(alphac)*(phi/Chi_c)*(1-np.exp(-Chi_c*Lc/Lt))
        I_w_term4a = theta_cf_term*(phi/lam1cl)*(1-np.exp(-lam1cl*L3/Lt))*np.sin(alphac)
        I_w_term4ba = theta_cf_term*(phi/lam1cl)*(np.exp(-lam1cl*L3/Lt)-np.exp(-lam1cl*(L3+X4)/Lt))*np.sin(alpha4)
        I_w_term4bb = theta_cf_term*(phi/lam1cl)*(np.exp(-lam1cl*(L3+X4)/Lt)-np.exp(-lam1cl*(L3+X4+X5)/Lt))*np.sin(alpha5)
        I_w_term4bc = theta_cf_term*(phi/lam1cl)*(np.exp(-lam1cl*(L3+X4+X5)/Lt)-np.exp(-lam1cl*(Lcl-L4)/Lt))*np.sin(alpha6)
        I_w_term4c = theta_cf_term*np.sin(alphah)*(phi/lam1cl)*(np.exp(-lam1cl*(Lcl-L4)/Lt)-np.exp(-lam1cl*Lcl/Lt))
        I_w_term =I_w_term1+I_w_term2a+I_w_term2ba+I_w_term2bb+I_w_term2bc+I_w_term2c+I_w_term3+I_w_term4a+I_w_term4ba+I_w_term4bb+I_w_term4bc+I_w_term4c

    cheq = n-Grm*I_w_term/np.power(re_ss,3)+(p*Lt*(2-b))/(2*D*np.power(re_ss,b))+K

    Re = np.real(cheq)
    Im = np.imag(cheq)

    return Re, Im, cheq
