# -*- coding: utf-8 -*-
"""
Created on Sun Jul 26 20:06 2020

@author: saikrishna nadella
"""
import sys
import os
thispath = os.getcwd()
UDpath = thispath+'\\UserDefined'
sys.path.append(UDpath)

from userTPPC import ThermoPhysicalProperties

def TA_Inp_func(uiHandle):
    State = 1
    try:
        #Gathering inputs corresponding to earlier dat file
        H = float(uiHandle.LE_H.text())
        D = float(uiHandle.LE_D.text())
        tw = float(uiHandle.LE_tw.text())
        dZc = float(uiHandle.LE_dZc.text())
        TDP = uiHandle.PB_TempDep.isChecked()
        SUF = uiHandle.PB_SSINIT.isChecked()
        K = float(uiHandle.LE_K.text())
        #solution parameters
        dt_UR = float(uiHandle.LE_dtUR.text())
        if SUF:
            w_init = 0.00001
            tforced = float(uiHandle.LE_tf.text())
        else:
            w_init = float(uiHandle.LE_wInit.text())
            tforced = 0
        wTol = float(uiHandle.LE_wTol.text())
        Tw1i = float(uiHandle.LE_TwiInit.text())
        Tw2i = float(uiHandle.LE_TwoInit.text())
        Tfi = float(uiHandle.LE_Tinit.text())
        tMax = float(uiHandle.LE_tMax.text())
        dtsave = float(uiHandle.LE_dtSave.text())
        #property data input
        if TDP:
            (a,b) = ThermoPhysicalProperties(Tfi,3)
            b = [float(x) for x in b]
            (rhof, cpf, kf, muf, betaf, rhow, cpw, kw) = b
        else:
            kf = float(uiHandle.LE_k.text())
            rhof = float(uiHandle.LE_rho.text())
            cpf = float(uiHandle.LE_Cp.text())
            muf = float(uiHandle.LE_mu.text())
            betaf = float(uiHandle.LE_betaf.text())
            rhow = float(uiHandle.LE_rhow.text())
            cpw = float(uiHandle.LE_Cpw.text())
            kw = float(uiHandle.LE_kw.text())
        if uiHandle.PB_Monitor.isChecked():
            dtUpdate = float(uiHandle.LE_dtUpdate.text())
        else:
            dtUpdate = 1e6
        dtMax = float(uiHandle.LE_dtMax.text())
        tShift =float(uiHandle.LE_tShift.text())
    except(ArithmeticError,ValueError,TypeError):
        return None, 0
    
    inpdata1 = (rhof, cpf, muf, kf, betaf, kw, cpw, rhow, K, TDP, dZc, SUF)
    inpdata2 = (D, H, tw, dt_UR, w_init, wTol, Tw1i, Tw2i, Tfi, tMax, dtsave, tforced,dtUpdate,dtMax,tShift)
    
    return (inpdata1,inpdata2), State
