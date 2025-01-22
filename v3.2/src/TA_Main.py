# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 23:00 2020

@author: saikrishna Nadella
"""
import os
import sys
import numpy as np
from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)
thispath = os.getcwd()
UDpath = thispath+'\\UserDefined'
sys.path.append(UDpath)
from user_f import f_udc
from userTPPC import ThermoPhysicalProperties
from user_Nu import Nu_udc
from Nu_main import Nu_main_func
import time


def TA_Main_func(ThreadHandle,uiHandle,inpdata):
    starttime = time.time()
    ThreadHandle.sigTOut.emit('Reading Inputs ...')
    thisPath = os.getcwd()
    caseLabel = uiHandle.LE_CaseLabel.text()
    if len(caseLabel)==0:
        caseLabel = 'Default'
    outPath = thisPath+'\\outputs\\TA\\'+caseLabel
    os.makedirs(outPath,exist_ok=True)

    try:
        g = float(uiHandle.LE_g.text())
    except(ValueError,TypeError):
        g = 9.81
        ThreadHandle.sigTOut.emit('''\n Gravitational Acceleration is taken as 9.81 m/s2''')
    
    (rhof, cpf, muf, kf, betaf, kw, cpw, rhow, K, TDP, dZc, SUF)=inpdata[0]
    (D, H, tw, dt_UR, w_init, wTol, Tw1i, Tw2i, Tfi, tMax, dtsave, tforced,dtUpdate,dtMax,tShift)=inpdata[1]
    #some checks
    if dt_UR>1:
        ThreadHandle.sigTOut.emit('''\n Warning: Time Step Under-Relaxation specified > 1. 0.9 is taken as Default''')
        dt_UR = 0.9

    #reading Geometry, Mesh & BCs
    secNs = uiHandle.UIT.rowCount()
    inpMesh = np.zeros((secNs,11))
    for sec in range(secNs):
        for col in range(11):
            inpMesh[sec,col]=float(uiHandle.UIT.item(sec,col).text())
    secType = inpMesh[:,0].astype('int32')
    nSec = inpMesh[:,2].astype('int32')
    if max(nSec%2)!=0:
        nSec = nSec%2+nSec
        ThreadHandle.sigTOut.emit('''\nOdd number of nodes found in one or more sections.\n The node count is raised by one to make even.''')        
    LSec = inpMesh[:,1]
    angSec = inpMesh[:,3]*np.pi/180
    hlSec = inpMesh[:,8]
    QSec = inpMesh[:,4]
    Tinf = inpMesh[:,10]
    QintSec = inpMesh[:,6]
    nTotal = sum(nSec)
    Lt = sum(LSec)

    # Parameter estimation
    A = np.pi*0.25*D**2
    D1 = D+tw
    D2 = D+2*tw
    A1 = np.pi*0.25*(D1**2-D**2)
    A2 = np.pi*0.25*(D2**2-D1**2)
    Cf = rhof*cpf
    Cw = rhow*cpw
    Prf = muf*cpf/kf
    alphaf = kf/Cf
    alphaw = kf/Cw
    Rw = np.log(D2/D)/(2*np.pi*kw)
    P2 = np.pi*D2
    P = np.pi*D
    qSec = QSec/(LSec*P2)
    qintSec = QintSec/(LSec*A)

    #check closedness of loop
    ThreadHandle.sigTOut.emit('''\n Reading inputs is OVER \n Validating the geometry ...''')

    ZSec = LSec*np.sin(angSec)
    HSec = LSec*np.cos(angSec)
    if sum(ZSec+HSec) > 0.001:
        ThreadHandle.sigTOut.emit('''\n Error: Loop is not closed.''')
        return 1

    #creating files for output
    ThreadHandle.sigTOut.emit('''\n Creating Files for data output''')
    a = outPath + '\\input_sec.csv'
    np.savetxt(a,inpMesh,delimiter=',',header='',comments='')
    
    if uiHandle.CB_w.isChecked():
        a = outPath+'\\w.dat'
        fidw = open(a,mode='wt')
        fidw.write("time \t mflow \t intT \n")
    if uiHandle.CB_T.isChecked():
        a = outPath +'\\T.dat'
        fidT = open(a,mode='wt')
        fidT.write('time \t T at all nodes ')
    if uiHandle.CB_Twi.isChecked():
        a = outPath +'\\Twi.dat'
        fidTwi = open(a,mode='wt')
        fidTwi.write('time \t Twi at all nodes ')
    if uiHandle.CB_Two.isChecked():
        a = outPath +'\\Two.dat'
        fidTwo = open(a,mode='wt')
        fidTwo.write('time \t Two at all nodes ')
    if uiHandle.CB_hi.isChecked():
        a = outPath +'\\hi.dat'
        fidh = open(a,mode='wt')
        fidh.write('time \t IHTC at all sections ')
    a = outPath+'\\inputData.dat'
    with open(a,'wt') as fidInp:
        fidInp.write('rhof, cpf, muf, kf, betaf, kw, cpw, rhow, K, TDP, dZc, SUF \n')
        fidInp.write(str(inpdata[0]))
        fidInp.write('\n D, H, tw, dt_UR, w_init, wTol, Tw1i, Tw2i, Tfi, tMax, dtsave, tforced,dtUpdate,dtMax,tShift \n')
        fidInp.write(str(inpdata[1]))
        fidInp.write('\n g = {:g}'.format(g))

    # Meshing the geometry
    ThreadHandle.sigTOut.emit('''Done. \n Meshing the loop''')
    ds = np.zeros(secNs)
    endNNo = np.zeros(secNs).astype('int32')
    begNNo = endNNo.copy()
    c = 0
    d = 1
    for i in range(secNs):
        begNNo[i]=d
        d = d+nSec[i]
        endNNo[i] = c+nSec[i]
        c = endNNo[i]
        ds[i] = LSec[i]/nSec[i]
               
    #Initialization
    ThreadHandle.sigTOut.emit('''Done. \n Initializing the solver''')
    Tn = np.ones(nTotal)
    T1n = Tn.copy()
    T2n = Tn.copy()
    hi = np.ones(secNs)

    To = Tfi*np.ones(nTotal)
    T1o = Tw1i*np.ones(nTotal)
    T2o = Tw2i*np.ones(nTotal)

    if SUF:
        # forced flow for tforced seconds at approximate flow rate as per Vijayan's relation
        for i in range(secNs):
            if(secType[i] == 1):
                Q = QSec[i]+QintSec[i]
        Grm_H = rhof**2*g*betaf*D**3*Q*H/(muf**3*A*cpf)
        Grm_dZc = Grm_H*dZc/H
        Ress = ((2/64)*Grm_dZc*D/Lt)**(1/(3-1))
        if Ress > 2000:
            Ress = ((2/0.316)*Grm_dZc*D/Lt)**(1/(3-0.25))
        w_init = Ress*np.pi*D*muf/4

    wo = w_init
    t = 0
    tsave = dtsave
    tUpdate = dtUpdate
    intT=0

    # Solver starts here
    ThreadHandle.sigTOut.emit('''Done. \n Solving...''')

    # monitor data preparation
    
    wup = [w_init]
    tup = [t]
    
    #delaying progress bar update
    lastProgUpdate = time.time()
    
    while t <= tMax:
        if uiHandle.toBeStopped:
            ThreadHandle.sigTOut.emit('Interruption Signal Received')
            return 1
        if t > tShift:
            hlSec = inpMesh[:,9]
            QSec = inpMesh[:,5]
            qSec = QSec/(LSec*P2)
            QintSec = inpMesh[:,7]
            qintSec = QintSec/(LSec*A)
    
        #writing data to files
        if uiHandle.CB_w.isChecked():
            fidw.write('{:g} \t {:g} \t {:g} \n'.format(t,wo,intT))
        if t > tsave:
            if uiHandle.CB_T.isChecked():
                fidT.write('\n {:g}\t'.format(t))
                np.savetxt(fidT,To,newline=',',fmt='%g')
            if uiHandle.CB_Twi.isChecked():
                fidTwi.write('\n {:g}\t'.format(t))
                np.savetxt(fidTwi,T1o,newline=',',fmt='%g')
            if uiHandle.CB_Two.isChecked():
                fidTwo.write('\n{:g}\t'.format(t))
                np.savetxt(fidTwo,T2o,newline=',',fmt='%g')
        # monitor update
        isMonitor = uiHandle.PB_Monitor.isChecked()
        if t > tUpdate and isMonitor:
            wup.append(wo)
            tup.append(t)
            ThreadHandle.sigMonitor.emit(isMonitor,np.array([tup,wup]))
            tUpdate = tUpdate + dtUpdate
        #property update
        Tavf = np.sum(To)/nTotal  #lengt weighted avg can be done in future
        if TDP:
            (a,TPPf) = ThermoPhysicalProperties(Tavf,3)
            TPPf = [float(x) for x in TPPf]
            (rhof, cpf, kf, muf, betaf)=TPPf[0:5]
            (rhow, cpw, kw)=TPPf[5:8]
                    # OString can be added here to warn user, if reqd
            # derived values based on properties
            Cf = rhof*cpf
            Cw = rhow*cpw
            Prf = muf*cpf/kf
            alphaf = kf/Cf
            alphaw = kw/Cw
            Rw = np.log(D2/D)/(2*np.pi*kw)
        #calculation of friction factor
        Re = abs(wo)*D/(A*muf)
        fCorr = uiHandle.PUM_friction.currentIndex()
        if fCorr == 0:
            psi = 1/(1+np.exp((Re-2530)/120))
            p = (64**psi)*(0.316**(1-psi))
            b = psi + (1-psi)*0.25
        elif fCorr == 1:
            (p, b) = f_udc(Re)
        #calculation of IHTC
        NuCorr = uiHandle.PUM_HT.currentIndex()
        if NuCorr == 0:
            for i in range(secNs):
                hi[i]=(kf/D)*Nu_main_func(LSec[i],D,abs(Re),Prf,secType[i])
        elif NuCorr == 1:
            for i in range(secNs):
                hi[i]=(kf/D)*Nu_udc(LSec[i],D,abs(Re),Prf,secType[i])
        if t > tsave and uiHandle.CB_hi.isChecked():
            fidh.write('\n {:g}\t'.format(t))
            np.savetxt(fidh,hi,newline=',',fmt='%g')

        #determination of time step size
        dsm = np.min(ds)
        hlm = np.max(hlSec)
        him = np.max(hi)
        dtwo = 2*alphaw/dsm**2 + (1/(Cw*A2))*(1/Rw+hlm*P2)
        dtwo = 1/dtwo
        dtwi = 2*alphaw/dsm**2 + (1/(Cw*A1))*(1/Rw+him*P)
        dtwi = 1/dtwi
        dtf = abs(wo)/(rhof*A*dsm) + 2*alphaf/dsm**2 + him*P/(Cf*A)
        dtf = 1/dtf
        dt = dt_UR*min([dtwo,dtwi,dtf,dtMax])

        # determination of wall outer shell temperature
          # for first section
        i=0
        a2w = alphaw*dt/ds[i]**2
        a2e = a2w
        a2p = 1-dt*(2*alphaw/ds[i]**2 + (1/(Cw*A2))*(1/Rw+hlSec[i]*P2))
        if min([a2w,a2e,a2p]) < 0:
            ThreadHandle.sigTOut.emit('''\n Error: Coefficients of energy equation of wall outer shell became negative''')
            return 1

        for j in range(nSec[i]-1): 
            k = begNNo[i]+j-1
            if secType[i] == 1: #heater
                b2p = (dt/(Cw*A2))*(T1o[k]/Rw+(hlSec[i]*Tinf[i]+qSec[i])*P2)
            else:
                b2p = (dt/(Cw*A2))*(T1o[k]/Rw+(hlSec[i]*Tinf[i])*P2)

            if k==0:
                T2n[k] = a2w*T2o[nTotal-1]+a2e*T2o[k+1]+a2p*T2o[k]+b2p  
            else:
                T2n[k] = a2w*T2o[k-1]+a2e*T2o[k+1]+a2p*T2o[k]+b2p
        k = endNNo[i]-1 
        dsav = (ds[i]+ds[i+1])/2
        a2w = alphaw*dt/(dsav*ds[i])
        a2e = alphaw*dt/(dsav*ds[i+1])
        a2p = 1-dt*(2*alphaw/(ds[i]*ds[i+1]) + (1/(Cw*A2))*(1/Rw+hlSec[i]*P2))
        if min([a2w,a2e,a2p]) < 0:
            ThreadHandle.sigTOut.emit('''\n Error: Coefficients of energy equation of wall outer shell became negative''')
            return 1

        if secType[i] == 1: #heater
            b2p = (dt/(Cw*A2))*(T1o[k]/Rw+(hlSec[i]*Tinf[i]+qSec[i])*P2)
        else:
            b2p = (dt/(Cw*A2))*(T1o[k]/Rw+(hlSec[i]*Tinf[i])*P2)
        T2n[k] = a2w*T2o[k-1]+a2e*T2o[k+1]+a2p*T2o[k]+b2p

        # calculation for second to secNs-1 sections
        for i in range(1,secNs-1): 
            a2w = alphaw*dt/ds[i]**2
            a2e = a2w
            a2p = 1-dt*(2*alphaw/ds[i]**2 + (1/(Cw*A2))*(1/Rw+hlSec[i]*P2))
            if min([a2w,a2e,a2p]) < 0:
                ThreadHandle.sigTOut.emit('''\n Error: Coefficients of energy equation of wall outer shell became negative''')
                return 1
            for j in range(nSec[i]-1):
                k = begNNo[i]+j-1
                if secType[i] == 1: #heater
                    b2p = (dt/(Cw*A2))*(T1o[k]/Rw+(hlSec[i]*Tinf[i]+qSec[i])*P2)
                else:
                    b2p = (dt/(Cw*A2))*(T1o[k]/Rw+(hlSec[i]*Tinf[i])*P2)
                
                T2n[k] = a2w*T2o[k-1]+a2e*T2o[k+1]+a2p*T2o[k]+b2p

            k = endNNo[i]-1
            dsav = (ds[i]+ds[i+1])/2
            a2w = alphaw*dt/(dsav*ds[i])
            a2e = alphaw*dt/(dsav*ds[i+1])
            a2p = 1-dt*(2*alphaw/(ds[i]*ds[i+1]) + (1/(Cw*A2))*(1/Rw+hlSec[i]*P2))
            if min([a2w,a2e,a2p]) < 0:
                ThreadHandle.sigTOut.emit('''\n Error: Coefficients of energy equation of wall outer shell became negative''')
                return 1

            if secType[i] == 1: #heater
                b2p = (dt/(Cw*A2))*(T1o[k]/Rw+(hlSec[i]*Tinf[i]+qSec[i])*P2)
            else:
                b2p = (dt/(Cw*A2))*(T1o[k]/Rw+(hlSec[i]*Tinf[i])*P2)
            T2n[k] = a2w*T2o[k-1]+a2e*T2o[k+1]+a2p*T2o[k]+b2p



         #calculation for last section
        i=secNs-1
        a2w = alphaw*dt/ds[i]**2
        a2e = a2w
        a2p = 1-dt*(2*alphaw/ds[i]**2 + (1/(Cw*A2))*(1/Rw+hlSec[i]*P2))
        if min([a2w,a2e,a2p]) < 0:
            ThreadHandle.sigTOut.emit('''\n Error: Coefficients of energy equation of wall outer shell became negative''')
            return 1
        for j in range(nSec[i]-1):
            k = begNNo[i]+j-1
            if secType[i] == 1: #heater
                b2p = (dt/(Cw*A2))*(T1o[k]/Rw+(hlSec[i]*Tinf[i]+qSec[i])*P2)
            else:
                b2p = (dt/(Cw*A2))*(T1o[k]/Rw+(hlSec[i]*Tinf[i])*P2)
            T2n[k] = a2w*T2o[k-1]+a2e*T2o[k+1]+a2p*T2o[k]+b2p

        k = endNNo[i]-1
        dsav = (ds[i]+ds[0])/2
        a2w = alphaw*dt/(dsav*ds[i])
        a2e = alphaw*dt/(dsav*ds[0])
        a2p = 1-dt*(2*alphaw/(ds[i]*ds[0]) + (1/(Cw*A2))*(1/Rw+hlSec[i]*P2))
        if min([a2w,a2e,a2p]) < 0:
            ThreadHandle.sigTOut.emit('''\n Error: Coefficients of energy equation of wall outer shell became negative''')
            return 1
        if secType[i] == 1: #heater
            b2p = (dt/(Cw*A2))*(T1o[k]/Rw+(hlSec[i]*Tinf[i]+qSec[i])*P2)
        else:
            b2p = (dt/(Cw*A2))*(T1o[k]/Rw+(hlSec[i]*Tinf[i])*P2)
        T2n[k] = a2w*T2o[k-1]+a2e*T2o[0]+a2p*T2o[k]+b2p

        # determination of wall inner shell termperature (new)
            # for first section
        i=0
        a1w = alphaw*dt/ds[i]**2
        a1e = a1w
        a1p = 1-dt*(2*alphaw/ds[i]**2 + (1/(Cw*A1))*(1/Rw+hi[i]*P))
        if min([a1w,a1e,a1p]) < 0:
            ThreadHandle.sigTOut.emit('''\n Error: Coefficients of energy equation of wall inner shell became negative''')
            return 1
        for j in range(nSec[i]-1):
            k = begNNo[i]+j-1
            b1p = (dt/(Cw*A1))*(T2n[k]/Rw+hi[i]*To[k]*P)
            if k==0:
                T1n[k] = a1w*T1o[nTotal-1]+a1e*T1o[k+1]+a1p*T1o[k]+b1p
            else:
                T1n[k] = a1w*T1o[k-1]+a1e*T1o[k+1]+a1p*T1o[k]+b1p
        k = endNNo[i]-1
        dsav = (ds[i]+ds[i+1])/2
        a1w = alphaw*dt/(dsav*ds[i])
        a1e = alphaw*dt/(dsav*ds[i+1])
        a1p = 1-dt*(2*alphaw/(ds[i]*ds[i+1]) + (1/(Cw*A1))*(1/Rw+hi[i]*P))
        if min([a1w,a1e,a1p]) < 0:
            ThreadHandle.sigTOut.emit('''\n Error: Coefficients of energy equation of wall inner shell became negative''')
            return 1
        b1p = (dt/(Cw*A1))*(T2n[k]/Rw+hi[i]*To[k]*P)
        T1n[k] = a1w*T1o[k-1]+a1e*T1o[k+1]+a1p*T1o[k]+b1p
            # calculation for second to secNs-1 sections
        for i in range(1,secNs-1):
            a1w = alphaw*dt/ds[i]**2
            a1e = a1w
            a1p = 1-dt*(2*alphaw/ds[i]**2 + (1/(Cw*A1))*(1/Rw+hi[i]*P))
            if min([a1w,a1e,a1p]) < 0:
                ThreadHandle.sigTOut.emit('''\n Error: Coefficients of energy equation of wall inner shell became negative''')
                return 1
            for j in range(nSec[i]-1):
                k = begNNo[i]+j-1
                b1p = (dt/(Cw*A1))*(T2n[k]/Rw+hi[i]*To[k]*P)
                T1n[k] = a1w*T1o[k-1]+a1e*T1o[k+1]+a1p*T1o[k]+b1p

            k = endNNo[i]-1
            dsav = (ds[i]+ds[i+1])/2
            a1w = alphaw*dt/(dsav*ds[i])
            a1e = alphaw*dt/(dsav*ds[i+1])
            a1p = 1-dt*(2*alphaw/(ds[i]*ds[i+1]) + (1/(Cw*A1))*(1/Rw+hi[i]*P))
            if min([a1w,a1e,a1p]) < 0:
                ThreadHandle.sigTOut.emit('''\n Error: Coefficients of energy equation of wall inner shell became negative''')
                return 1
            b1p = (dt/(Cw*A1))*(T2n[k]/Rw+hi[i]*To[k]*P)
            T1n[k] = a1w*T1o[k-1]+a1e*T1o[k+1]+a1p*T1o[k]+b1p

            #calculation for last section
        i=secNs-1
        a1w = alphaw*dt/ds[i]**2
        a1e = a1w
        a1p = 1-dt*(2*alphaw/ds[i]**2 + (1/(Cw*A1))*(1/Rw+hi[i]*P))
        if min([a1w,a1e,a1p]) < 0:
            ThreadHandle.sigTOut.emit('''\n Error: Coefficients of energy equation of wall inner shell became negative''')
            return 1
        for j in range(nSec[i]-1):
            k = begNNo[i]+j-1
            b1p = (dt/(Cw*A1))*(T2n[k]/Rw+hi[i]*To[k]*P)
            T1n[k] = a1w*T1o[k-1]+a1e*T1o[k+1]+a1p*T1o[k]+b1p

        k = endNNo[i]-1
        dsav = (ds[i]+ds[0])/2
        a1w = alphaw*dt/(dsav*ds[i])
        a1e = alphaw*dt/(dsav*ds[0])
        a1p = 1-dt*(2*alphaw/(ds[i]*ds[0]) + (1/(Cw*A1))*(1/Rw+hi[i]*P))
        if min([a1w,a1e,a1p]) < 0:
            ThreadHandle.sigTOut.emit('''\n Error: Coefficients of energy equation of wall inner shell became negative''')
            return 1
        b1p = (dt/(Cw*A1))*(T2n[k]/Rw+hi[i]*To[k]*P)
        T1n[k] = a1w*T1o[k-1]+a1e*T1o[0]+a1p*T1o[k]+b1p

        # determination of fluid temperature (new)
        i=0
        if wo >= 0:
            aw = alphaf*dt/ds[i]**2+abs(wo)*dt/(rhof*A*ds[i])
            ae = alphaf*dt/ds[i]**2
            ap = 1-dt*(abs(wo)/(rhof*A*ds[i])+2*alphaf/ds[i]**2+hi[i]*P/(Cf*A))
        else:
            aw = alphaf*dt/ds[i]**2
            ae = alphaf*dt/ds[i]**2+abs(wo)*dt/(rhof*A*ds[i+1])
            ap = 1-dt*(abs(wo)/(rhof*A*ds[i])+2*alphaf/ds[i]**2+hi[i]*P/(Cf*A))
        
        if min([aw,ae,ap]) < 0:
            ThreadHandle.sigTOut.emit('''\n Error: Coefficients of energy equation of fluid became negative''')
            return 1
        for j in range(nSec[i]-1): 
            k = begNNo[i]+j-1
            bp = (dt/(Cf*A))*(hi[i]*T1n[k]*P+qintSec[i]*A)
            if k==0:
                Tn[k] = aw*To[nTotal-1]+ae*To[k+1]+ap*To[k]+bp
            else:
                Tn[k] = aw*To[k-1]+ae*To[k+1]+ap*To[k]+bp
        k = endNNo[i]-1
        dsav = (ds[i]+ds[i+1])/2
        if wo>=0:
            aw = alphaf*dt/(dsav*ds[i])+abs(wo)*dt/(rhof*A*ds[i])
            ae = alphaf*dt/(dsav*ds[i+1]) #corrected
            ap = 1-dt*(2*alphaf/(ds[i]*ds[i+1])+abs(wo)/(rhof*A*ds[i])+ hi[i]*P/(Cf*A))
        else:
            aw = alphaf*dt/(dsav*ds[i])
            ae = alphaf*dt/(dsav*ds[i+1])+abs(wo)*dt/(rhof*A*ds[i+1])
            ap = 1-dt*(2*alphaf/(ds[i]*ds[i+1]) +abs(wo)/(rhof*A*ds[i+1])+ hi[i]*P/(Cf*A))

        if min([aw,ae,ap]) < 0:
            ThreadHandle.sigTOut.emit('''\n Error: Coefficients of energy equation of fluid became negative''')
            return 1
        bp = (dt/(Cf*A))*(hi[i]*T1n[k]*P+qintSec[i]*A)
        Tn[k] = aw*To[k-1]+ae*To[k+1]+ap*To[k]+bp
            # calculation for second to secNs-1 sections
        for i in range(1,secNs-1): 
            if wo >= 0:
                aw = alphaf*dt/ds[i]**2+abs(wo)*dt/(rhof*A*ds[i])
                ae = alphaf*dt/ds[i]**2
                ap = 1-dt*(abs(wo)/(rhof*A*ds[i])+2*alphaf/ds[i]**2+hi[i]*P/(Cf*A))
            else:
                aw = alphaf*dt/ds[i]**2
                ae = alphaf*dt/ds[i]**2+abs(wo)*dt/(rhof*A*ds[i+1])
                ap = 1-dt*(abs(wo)/(rhof*A*ds[i])+2*alphaf/ds[i]**2+hi[i]*P/(Cf*A))
            if min([aw,ae,ap]) < 0:
                ThreadHandle.sigTOut.emit('''\n Error: Coefficients of energy equation of fluid became negative''')
                return 1
            for j in range(nSec[i]-1): 
                k = begNNo[i]+j-1
                bp = (dt/(Cf*A))*(hi[i]*T1n[k]*P+qintSec[i]*A)
                Tn[k] = aw*To[k-1]+ae*To[k+1]+ap*To[k]+bp
            
            k = endNNo[i]-1
            dsav = (ds[i]+ds[i+1])/2
            if wo >= 0:
                aw = alphaf*dt/(dsav*ds[i])+abs(wo)*dt/(rhof*A*ds[i])
                ae = alphaf*dt/(dsav*ds[i+1])
                ap = 1-dt*(2*alphaf/(ds[i]*ds[i+1]) +abs(wo)/(rhof*A*ds[i])+ hi[i]*P/(Cf*A))
            else:
                aw = alphaf*dt/(dsav*ds[i])
                ae = alphaf*dt/(dsav*ds[i+1])+abs(wo)*dt/(rhof*A*ds[i+1])
                ap = 1-dt*(2*alphaf/(ds[i]*ds[i+1]) +abs(wo)/(rhof*A*ds[i+1])+ hi[i]*P/(Cf*A))

            if min([aw,ae,ap]) < 0:
                ThreadHandle.sigTOut.emit('''\n Error: Coefficients of energy equation of fluid became negative''')
                return 1
            bp = (dt/(Cf*A))*(hi[i]*T1n[k]*P+qintSec[i]*A)
            Tn[k] = aw*To[k-1]+ae*To[k+1]+ap*To[k]+bp
        
            #calculation for last section
        i=secNs-1
        if wo >= 0:
            aw = alphaf*dt/ds[i]**2+abs(wo)*dt/(rhof*A*ds[i])
            ae = alphaf*dt/ds[i]**2
            ap = 1-dt*(abs(wo)/(rhof*A*ds[i])+2*alphaf/ds[i]**2+hi[i]*P/(Cf*A))
        else:
            aw = alphaf*dt/ds[i]**2
            ae = alphaf*dt/ds[i]**2+abs(wo)*dt/(rhof*A*ds[0])
            ap = 1-dt*(abs(wo)/(rhof*A*ds[0])+2*alphaf/ds[i]**2+hi[i]*P/(Cf*A))
      
        if min([aw,ae,ap]) < 0:
            ThreadHandle.sigTOut.emit('''\n Error: Coefficients of energy equation of fluid became negative''')
            return 1
        for j in range(nSec[i]-1):
            k = begNNo[i]+j-1
            bp = (dt/(Cf*A))*(hi[i]*T1n[k]*P+qintSec[i]*A)
            Tn[k] = aw*To[k-1]+ae*To[k+1]+ap*To[k]+bp

        k = endNNo[i]-1
        dsav = (ds[i]+ds[0])/2
        if wo >= 0:
            aw = alphaf*dt/(dsav*ds[i])+abs(wo)*dt/(rhof*A*ds[i])
            ae = alphaf*dt/(dsav*ds[0])
            ap = 1-dt*(2*alphaf/(ds[i]*ds[0]) +abs(wo)/(rhof*A*ds[i])+ hi[i]*P/(Cf*A))
        else:
            aw = alphaf*dt/(dsav*ds[i])
            ae = alphaf*dt/(dsav*ds[0])+abs(wo)*dt/(rhof*A*ds[0])
            ap = 1-dt*(2*alphaf/(ds[i]*ds[0]) +abs(wo)/(rhof*A*ds[0])+ hi[i]*P/(Cf*A))
        if min([aw,ae,ap]) < 0:
            ThreadHandle.sigTOut.emit('''\n Error: Coefficients of energy equation of fluid became negative''')
            return 1
        bp = (dt/(Cf*A))*(hi[i]*T1n[k]*P+qintSec[i]*A)
        Tn[k] = aw*To[k-1]+ae*To[0]+ap*To[k]+bp

        # integration of temperature using simpson 1/3rd rule
          # we can use scipy.integrate.simps function also.
        intT = 0
        #k=0  #old implementation
        k=-1
        for i in range(secNs):
            intTSec = Tn[k] + Tn[k+nSec[i]]
            for j in range(1,int(nSec[i]/2)):
                intTSec = intTSec + 2*Tn[k+2*j]

            for j in range(int(nSec[i]/2)):
                intTSec = intTSec + 4*Tn[k+2*j+1]

            intTSec = intTSec*ds[i]*np.sin(angSec[i])/3
            intT = intT + intTSec
            k = k + nSec[i]
        
        # determination of mass flow rate (new)
        if t > tforced or SUF != 1: 
            wg = wo
            x = K*dt/(2*rhof*A*Lt)
            y = p*muf**b*dt/(2*D**(1+b)*rhof*A**(1-b))
            z = -wo-g*rhof*betaf*A*dt*intT/Lt
            count =1

            while True:
                fwg = x*wg*abs(wg)+y*wg*abs(wg)**(1-b)+wg+z
                fdwg = 2*x*abs(wg)+y*(2-b)*abs(wg)**(1-b)+1
                wn = 1*(wg - fwg/fdwg)+0*wg   #potential to put under relaxation
                if abs((wn-wg)/wg) < wTol and fwg <= 1e-16:
                    break

                count=count+1
                if count == 200:
                    ThreadHandle.sigTOut.emit('''Error: Implicit flow rate solver did not converge''')
                    return 1
                wg=wn
        else:
            wn=wo

        #transferring data from new storage to old storage
        To = Tn
        T1o = T1n
        T2o = T2n

        wo = wn
        if wo.imag != 0:
            ThreadHandle.sigTOut.emit('''Error: Flow rate bacame complex!''')
            return 1

        if t>tsave:
            tsave=tsave+dtsave

        t=t+dt
        
        if(time.time() - lastProgUpdate)>2.0:            
            ThreadHandle.sigProgBar.emit(int(np.floor(100*t/tMax)))
            lastProgUpdate = time.time()
            
    # ------------ time loop over
    # writing End state to file
    ThreadHandle.sigTOut.emit('''Calculations over :)\n
                                Writing the last time step data to EndData.dat''')
    ThreadHandle.sigProgBar.emit(100)
    for i in range(secNs):
        if secType[i] == 1:
            Q = QSec[i]+QintSec[i]
            dTh = To[endNNo[i]-1]-To[begNNo[i]-2]
        elif secType[i] == 2:
            hic = hi[i]
            hlc = hlSec[i]
            dTc = To[begNNo[i]-2]-To[endNNo[i]-1]
    Grm_H = rhof**2*g*betaf*D**3*Q*H/(muf**3*A*cpf)
    Grm_dZc = Grm_H*dZc/H
    Iss = intT/(H*dTh)
    Stmic = (4*Lt/D)*(hic*D/kf)/(Re*Prf)
    Stmlc = (4*Lt/D)*(hlc*D2/kf)/(Re*Prf)
    Xi = (4*Lt/D)*(1/(kf*Rw*np.pi*Re*Prf))
    Stmc = 1/(1/Stmic+1/Xi+1/Stmlc)

    a = outPath+'\\EndData.dat'
    with open(a,'wt') as fidEnd:
        fidEnd.write('Check the Grm for validity in case of using both SHF and IHG at heater\n')
        fidEnd.write('Grm_H\tGrm_dZc\tIss\tRe \tStmic\tStmoc\tXi\tStmc')
        fidEnd.write('\n{:g}\t{:g}\t{:g}\t{:g}\t{:g}\t{:g}\t{:g}\t{:g}'.format(Grm_H,Grm_dZc,Iss,Re,Stmic,Stmlc,Xi,Stmc))
        fidEnd.write('\ndTh\tdTc\tTavgf\n{:g}\t{:g}\t{:g}'.format(dTh,dTc,Tavf))
        fidEnd.write('\np\tb\n{:g}\t{:g}'.format(p,b))
        fidEnd.write('\nhi\n')
        fidEnd.write(str(hi))
        fidEnd.write('\n ho\n')
        fidEnd.write(str(hlSec))
        fidEnd.write('\n ndk ndcp ndro Prf \n {:g}\t{:g}\t{:g}\t{:g}\n'.format(kw/kf,cpw/cpf,rhow/rhof,Prf))
        fidEnd.write('Mass flow rate is {:g} \n Velocity is {:g} \n\n\n Temperature profile is \n\n Tfluid\tTwi\tTwo\n'.format(wo,wo/(rhof*A)))
        for i in range(nTotal):
            fidEnd.write('{:d} \t {:g} \t {:g} \t {:g} \n'.format(i,To[i],T1o[i],T2o[i]))
    # closing open files.  we can have some wrapper class to automatically close all files at end with single command.
    if uiHandle.CB_w.isChecked():
        fidw.close()
    if uiHandle.CB_T.isChecked():
        fidT.close()
    if uiHandle.CB_Twi.isChecked():
        fidTwi.close()
    if uiHandle.CB_Two.isChecked():
        fidTwo.close()
    if uiHandle.CB_hi.isChecked():
        fidh.close()

    ThreadHandle.sigTOut.emit('''Solution completed successfully \n Check output\\TA folder for output data\n''')

    runTime = time.time() - starttime
    ThreadHandle.sigTOut.emit(r'Simulation Time = {:g} seconds'.format(runTime))
    return 1
