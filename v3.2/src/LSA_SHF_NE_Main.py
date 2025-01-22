# -*- coding: utf-8 -*-
"""
Created on Sun June 07 15:04:08 2020

@author: saikrishna Nadella
"""
import numpy as np
from LSA_SHF_SS import LSA_SHF_SS_func
from LSA_SHF_NE_CEgen import LSA_SHF_NE_CEgen_func


def LSA_SHF_NE_Main_func(inpdata, Grm, nImin, nImax, nIstep, StmocMin, StmocMax, StmocStep, uiHandle):

    # input data assignment
    ModFlag = inpdata[2][12]
    HXEffect = inpdata[3][6]
    K = inpdata[2][3]
    hoh = inpdata[2][0]
    hohl = inpdata[2][1]
    hocl = inpdata[2][2]
    Prf = inpdata[3][3]
    
    #model selection being printed to console
    if ModFlag==1:
        if HXEffect==0:
            print('Warning!: Usage of unverified combination of model parameters! \n \t You are using BM+WTI+HX+HL model without HXEffect.')
        else:
            print('Model: BM+WTI+HX+HL')

        if (hoh < 1e-2) or (hohl < 1e-2) or (hocl <1e-2):
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
        
    if ModFlag != 4:
        print('Prandtl number is {:g}'.format(Prf))
    
    size_nI = int(np.floor((nImax-nImin)/nIstep+0.001)+1) #0.001 is added to take care of truncation error)
    size_Stmoc = int(np.floor((StmocMax-StmocMin)/StmocStep+0.001)+1)

    #plot data carriers initialization
    Re = np.ones((size_Stmoc,size_nI))
    Im = np.ones((size_Stmoc,size_nI)) 
    Stmoc = StmocMin;
    cSt = 0
    while Stmoc <= StmocMax:
        # steady state calculations
        (Errsignal,SS2CEdata,SSDatai) = LSA_SHF_SS_func(inpdata,Grm,Stmoc,uiHandle) # SSDatai will be used in fout case. currently no use
        if Errsignal == 0:
            uiHandle.TE_TEXTOUT.setPlainText('Due to error signal from SS_code, Graphing code is terminated')
            return 
        nI = nImin
        cnI = 0
        while nI <= nImax:
            (Re[cSt,cnI], Im[cSt,cnI], cheq) = LSA_SHF_NE_CEgen_func(inpdata,SS2CEdata,nI,Stmoc,Grm)
            if np.isnan(Re[cSt,cnI]) or np.isnan(Im[cSt,cnI]):
                uiHandle.TE_TEXTOUT.append('''\n NaN is found in Re Im data carriers''')
            nI = nI+nIstep
            cnI=cnI+1

        Stmoc = Stmoc+StmocStep
        cSt = cSt+1
    #filling is over
    
    (Y, Z) = np.meshgrid(np.arange(nImin,nImax+nIstep/10,nIstep),np.arange(StmocMin,StmocMax+StmocStep/10,StmocStep)) #taking care of Max value to be included if comes by stepping
                
    #plotting the contours
    uiHandle.TW_plot.setCurrentIndex(0)        
    uiHandle.MplWidget_NE.canvas.axes.clear()
    '''fig.subplots_adjust(left, bottom, right, top, wspace, hspace)
        top=0.933,
        bottom=0.142,
        left=0.109,
        right=0.956,
        hspace=0.165,
        wspace=0.2'''
    cr = uiHandle.MplWidget_NE.canvas.axes.contour(Y,Z,Re,colors='black',linestyles=['-'],linewidths=1,levels = [0])
    ci = uiHandle.MplWidget_NE.canvas.axes.contour(Y,Z,Im,colors='red',linestyles=['-'],linewidths=1,levels = [0])
    uiHandle.MplWidget_NE.canvas.axes.set_xlabel('nI')
    if ModFlag == 4:
        uiHandle.MplWidget_NE.canvas.axes.set_ylabel('Stm')
    else:
        uiHandle.MplWidget_NE.canvas.axes.set_ylabel('Stmoc')
    #artists,labels=cr.legend_elements()
    #uiHandle.MplWidget_NE.canvas.axes.legend(handles=[cr,ci])
    uiHandle.MplWidget_NE.canvas.axes.set_title(r'Grm = {:g}'.format(Grm))
    uiHandle.MplWidget_NE.canvas.figure.tight_layout()
    uiHandle.MplWidget_NE.canvas.draw()

    uiHandle.TE_TEXTOUT.append('''\n Successful''')
    return None

