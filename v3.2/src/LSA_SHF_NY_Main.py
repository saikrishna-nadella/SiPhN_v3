# -*- coding: utf-8 -*-
"""
Created on Sun May 31 16:40:24 2020

@author: saikr
"""

import numpy as np
from LSA_SHF_SS import LSA_SHF_SS_func
from LSA_SHF_NE_CEgen import LSA_SHF_NE_CEgen_func


def LSA_SHF_NY_Main_func(inpdata,Grm,Stmoc,fMin,fMax,fStep,uiHandle):
    
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
        
    if ModFlag != 4:
        print('Prandtl number is {:g}'.format(Prf))

    (Errsignal,SS2CEdata,SSDatai) = LSA_SHF_SS_func(inpdata,Grm,Stmoc,uiHandle)
    kRange = int(np.floor((fMax-fMin)/fStep))+1
    store = np.ones((kRange,2)) #can be made 3 if f needs to be stored
    for k in range(kRange):
        nI = (fMin+fStep*k)
        (store[k][0], store[k][1], cheq) = LSA_SHF_NE_CEgen_func(inpdata,SS2CEdata,nI,Stmoc,Grm)

    uiHandle.TW_plot.setCurrentIndex(1)
    uiHandle.MplWidget_NY.canvas.axes.clear()
    uiHandle.MplWidget_NY.canvas.axes.plot(store[:,0],store[:,1],color='blue',linestyle='-',linewidth=1)
    uiHandle.MplWidget_NY.canvas.axes.plot([0],[0],color='black',marker='*')
    uiHandle.MplWidget_NY.canvas.axes.set_ylabel('F-Imag')
    uiHandle.MplWidget_NY.canvas.axes.set_xlabel('F-Real')
    uiHandle.MplWidget_NY.canvas.axes.set_title(r'Grm = {:g}, Stmc = {:g}'.format(Grm,Stmoc))
    uiHandle.MplWidget_NY.canvas.figure.tight_layout()
    uiHandle.MplWidget_NY.canvas.draw()

    uiHandle.TE_TEXTOUT.append('''\n Successful''')
    return None
