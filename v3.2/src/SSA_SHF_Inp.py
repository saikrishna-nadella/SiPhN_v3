# -*- coding: utf-8 -*-
"""
Created on Sun May 31 12:56:09 2020

@author: saikr
"""

def SSA_SHF_Inp_func(uiHandle,DIM):
    State=1
    #gathering inputs from ui
    L1 = float(uiHandle.LE_L1.text())
    L2 = float(uiHandle.LE_L2.text())
    L3 = float(uiHandle.LE_L3.text())
    L4 = float(uiHandle.LE_L4.text())
    Lh = float(uiHandle.LE_Lh.text())
    Lc = float(uiHandle.LE_Lc.text())
    H = float(uiHandle.LE_H.text())
    D = float(uiHandle.LE_D.text())
    tw = float(uiHandle.LE_tw.text())
    X1 = float(uiHandle.LE_X1.text())
    X2 = float(uiHandle.LE_X2.text())
    X3 = float(uiHandle.LE_X3.text())
    X4 = float(uiHandle.LE_X4.text())
    X5 = float(uiHandle.LE_X5.text())
    X6 = float(uiHandle.LE_X6.text())
    alpha1 = float(uiHandle.LE_alpha1.text())
    alpha2 = float(uiHandle.LE_alpha2.text())
    alpha3 = float(uiHandle.LE_alpha3.text())
    alpha4 = float(uiHandle.LE_alpha4.text())
    alpha5 = float(uiHandle.LE_alpha5.text())
    alpha6 = float(uiHandle.LE_alpha6.text())
    alphah = float(uiHandle.LE_alphah.text())
    alphac = float(uiHandle.LE_alphac.text())
    

    try:
        K = float(uiHandle.LE_K.text())
    
        #property data input
        kf = float(uiHandle.LE_k.text())
        rhof = float(uiHandle.LE_rho.text())
        cpf = float(uiHandle.LE_Cp.text())
        muf = float(uiHandle.LE_mu.text())
        betaf = float(uiHandle.LE_betaf.text())
        kw = float(uiHandle.LE_kw.text())
        rhow = float(uiHandle.LE_rhow.text())
        cpw = float(uiHandle.LE_Cpw.text())
        
        Prf = muf*cpf/kf
        #solution parameters
        ReTol = float(uiHandle.LE_Retol.text())
        Re_UR = float(uiHandle.LE_ReUR.text())
        Re_g = float(uiHandle.LE_Reinit.text())
        if DIM==1:
            if uiHandle.PUM_Model.currentIndex()==5:
                hoh = float(uiHandle.LE_hoh.text())
                hohl = float(uiHandle.LE_hohl.text())
                hocl = float(uiHandle.LE_hocl.text())
            else:
                hoh = 0
                hohl = 0
                hocl = 0
                    
            Tamb = float(uiHandle.LE_Tamb.text())
            Ts = float(uiHandle.LE_Ts.text())
            theta_s = Ts-Tamb
        else:
            if uiHandle.PUM_Model.currentIndex()==5:
                Stmoh = float(uiHandle.LE_Stmoh.text())
                Stmohl = float(uiHandle.LE_Stmohl.text())
                Stmocl = float(uiHandle.LE_Stmocl.text())
            else:
                Stmoh = 0
                Stmohl = 0
                Stmocl = 0
            theta_s = float(uiHandle.LE_thetas.text())
            try:
                Tamb = float(uiHandle.LE_Tamb.text())
                #Ts = float(uiHandle.LE_Ts.text())
            except(ArithmeticError, TypeError, ValueError):
                Tamb = 0
                #Ts = 0
    except (ArithmeticError,TypeError,ValueError):
        uiHandle.TE_DIMout.setPlainText('Input Error: Please check inputs!')
        uiHandle.TE_NDout.setPlainText('Input Error: Please check inputs!')
        return None,0
    #Modflag determination
    Model = uiHandle.PUM_Model.currentIndex()
    if Model==0:
        uiHandle.TE_DIMout.setPlainText('Please select a model to proceed')
        uiHandle.TE_NDout.setPlainText('Please select a model to proceed')
        State=0
        return None,State
    elif Model==1:
        ModFlag = 4
        HXEffect=0
    elif Model==2:
         ModFlag=3
         HXEffect=0
    elif Model==3:
         ModFlag=2
         HXEffect=0
    elif Model==4:
         ModFlag=2
         HXEffect=1
    else:
         ModFlag=1
         HXEffect=1

    #avoiding zero heat loss in FullModel
    if DIM==1:
        if Model==5:
            print('checking compatibility of heat loss input to selected model')
            if max(hoh, hohl, hocl) < 0.01:
                uiHandle.TE_DIMout.append('Error: Too Low Heat Loss Coefficients Specified. Please select BM+WTI+HX model and run again')
                State=0
                return None,State
            else:
                if hoh < 0.01:
                    uiHandle.TE_DIMout.append('Warning: hoh is too low. It is made equal to 0.01 for model compatibility')
                    hoh = 0.01
                if hohl < 0.01:
                    uiHandle.TE_DIMout.append('Warning: hohl is too low. It is made equal to 0.01 for model compatibility')
                    hohl = 0.01
                if hocl < 0.01:
                    uiHandle.TE_DIMout.append('Warning: hocl is too low. It is made equal to 0.01 for model compatibility')
                    hocl = 0.01
    else:
        if Model==5:
            print('checking compatibility of heat loss input to selected model')
            if max(Stmoh, Stmohl, Stmocl) < 0.0001:
                uiHandle.TE_NDout.append('Error: Too Low Heat Loss Coefficients Specified. Please select BM+WTI+HX model and run again')
                State=0
                return None,State
            else:
                if Stmoh < 0.0001:
                    uiHandle.TE_NDout.append('Warning: Stmoh is too low. It is made equal to 0.0001 for model compatibility')
                    Stmoh = 0.0001
                if Stmohl < 0.0001:
                    uiHandle.TE_NDout.append('Warning: Stmohl is too low. It is made equal to 0.0001 for model compatibility')
                    hohl = 0.0001
                if Stmocl < 0.0001:
                    uiHandle.TE_NDout.append('Warning: Stmocl is too low. It is made equal to 0.0001 for model compatibility')
                    Stmocl = 0.0001
        
    #data assembly
    inpdata1 = (L1,L2,L3,L4,Lc,Lh,X1,X2,X3,X4,X5,X6,0,0)
    inpdata2 = (H,D,tw,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alphah,alphac,0,0,0)
    if DIM==1:
        inpdata3 = (hoh,hohl,hocl,K,Ts,Tamb,0,Re_UR,0,0,0,Re_g,ModFlag,0)
    else:
        inpdata3 = (Stmoh,Stmohl,Stmocl,K,0,Tamb,0,Re_UR,0,0,0,Re_g,ModFlag,0)
        
    inpdata4 = (rhow,cpw,kw,Prf,ReTol,0,HXEffect,theta_s,kf,betaf,rhof,cpf,0,0)
    return (inpdata1,inpdata2,inpdata3,inpdata4),State
