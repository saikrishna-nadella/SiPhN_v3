# -*- coding: utf-8 -*-
"""
Created on Sun June 07 14:41:55 2020

@author: saikrishna nadella
"""
def LSA_SHF_Inp_func(uiHandle):
    State=1
    try:
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
        
        Model = uiHandle.PUM_Model.currentIndex()
        
        #property data input
        ndrho = float(uiHandle.LE_NDrho.text())
        ndcp = float(uiHandle.LE_NDCp.text())
        ndk = float(uiHandle.LE_NDk.text())
        Prf = float(uiHandle.LE_Pr.text())
        
        #friction
        K = float(uiHandle.LE_K.text())
        
        #heat loss
        if Model==5:
            hoh = float(uiHandle.LE_hoh.text())
            hohl = float(uiHandle.LE_hohl.text())
            hocl = float(uiHandle.LE_hocl.text())
            kf = float(uiHandle.LE_kf.text())
        else:
            hoh = 0
            hohl = 0
            hocl = 0
            kf = 0
        
        theta_s = float(uiHandle.LE_thetas.text())
        
        #solution parameters
        ReTol = float(uiHandle.LE_Retol.text())
        Re_UR = float(uiHandle.LE_ReUR.text())
        Re_g = float(uiHandle.LE_Reinit.text())
     
        #Modflag determination
        if Model==0:
            uiHandle.TE_TEXTOUT.append('Please select a model to proceed')
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
    except(ArithmeticError,ValueError,TypeError):
        return None, 0
    #avoiding zero heat loss in FullModel
    if Model==5:
        if max(hoh, hohl, hocl) < 0.01:
            uiHandle.TE_TEXTOUT.append('Too Low Heat Loss Coefficients Specified. Please select BM+WTI+HX model and run again')
            State=0
            return None,State
        else:
            if hoh < 0.01:
                uiHandle.TE_TEXTOUT.append('hoh is too low. It is made equal to 0.01 for model compatibility')
                hoh = 0.01
            if hohl < 0.01:
                uiHandle.TE_TEXTOUT.append('hohl is too low. It is made equal to 0.01 for model compatibility')
                hohl = 0.01
            if hocl < 0.01:
                uiHandle.TE_TEXTOUT.append('hocl is too low. It is made equal to 0.01 for model compatibility')
                hocl = 0.01
    
    #data assembly
    inpdata1 = (L1,L2,L3,L4,Lc,Lh,X1,X2,X3,X4,X5,X6,0,0)
    inpdata2 = (H,D,tw,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alphah,alphac,0,0,0)
    inpdata3 = (hoh,hohl,hocl,K,0,0,0,Re_UR,0,0,0,Re_g,ModFlag,0)
    inpdata4 = (ndrho,ndcp,ndk,Prf,ReTol,0,HXEffect,theta_s,kf,0,0,0,0,0)
    return (inpdata1,inpdata2,inpdata3,inpdata4),State


