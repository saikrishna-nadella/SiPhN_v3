# -*- coding: utf-8 -*-
"""
Created on Fri Jun 05 21:19:20 2020

@author: saikrishna Nadella
"""


from PyQt5 import QtWidgets, uic
import sys
import os
thispath = os.getcwd()
srcpath = thispath+'\\src'
sys.path.append(srcpath)
from LSA_SHF_Inp import LSA_SHF_Inp_func
from LSA_SHF_NE_Main import LSA_SHF_NE_Main_func
from LSA_SHF_NY_Main import LSA_SHF_NY_Main_func
from Geom_Val_Plot import Geom_Val_Plot_func
from triggerGeomInputs import triggerGeomInputs_func
from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)


class clLSA_SHF(QtWidgets.QMainWindow):
    def __init__(self):
        super(clLSA_SHF,self).__init__()
        uic.loadUi('src\LSA_SHF.ui',self)
        #adding toolbars for plot manipulations
        self.widget_NE.setLayout(QtWidgets.QVBoxLayout())
        self.widget_NE.layout().addWidget(NavigationToolbar(self.MplWidget_NE.canvas, self))
        self.widget_NY.setLayout(QtWidgets.QVBoxLayout())
        self.widget_NY.layout().addWidget(NavigationToolbar(self.MplWidget_NY.canvas, self))
        self.PB_NY.setEnabled(False)
        self.PB_NE.setEnabled(False)
        #
        triggerGeomInputs_func(self)
        self.PB_NY.clicked.connect(self._forwardNY)
        self.PB_NE.clicked.connect(self._forwardNE)
        self.PB_Validate.clicked.connect(self._validateGeom)
        #changing any data disables the calculate buttons.
        self.LE_L1.textEdited.connect(self._disableCalc)
        self.LE_L2.textEdited.connect(self._disableCalc)
        self.LE_L3.textEdited.connect(self._disableCalc)
        self.LE_L4.textEdited.connect(self._disableCalc)
        self.LE_Lh.textEdited.connect(self._disableCalc)
        self.LE_Lc.textEdited.connect(self._disableCalc)
        self.LE_X3.textEdited.connect(self._disableCalc)
        self.LE_X6.textEdited.connect(self._disableCalc)
        self.LE_D.textEdited.connect(self._disableCalc)
        self.LE_H.textEdited.connect(self._disableCalc)
        self.LE_tw.textEdited.connect(self._disableCalc)
        self.LE_X1.textEdited.connect(self._disableCalc)
        self.LE_X2.textEdited.connect(self._disableCalc)
        self.LE_X4.textEdited.connect(self._disableCalc)
        self.LE_X5.textEdited.connect(self._disableCalc)
        self.LE_alphah.textEdited.connect(self._disableCalc)
        self.LE_alpha1.textEdited.connect(self._disableCalc)
        self.LE_alpha2.textEdited.connect(self._disableCalc)
        self.LE_alpha3.textEdited.connect(self._disableCalc)
        self.LE_alphac.textEdited.connect(self._disableCalc)
        self.LE_alpha4.textEdited.connect(self._disableCalc)
        self.LE_alpha5.textEdited.connect(self._disableCalc)
        self.LE_alpha6.textEdited.connect(self._disableCalc)
        self.PUM_Orientation.currentIndexChanged.connect(self._disableCalc)
        self.PUM_Orientation.currentIndexChanged.connect(self._triggerGeomInputs)
    def _triggerGeomInputs(self):
        triggerGeomInputs_func(self)
            
    def _disableCalc(self):
        self.PB_NE.setEnabled(False)
        self.PB_NY.setEnabled(False)        
        
    def _validateGeom(self):
        if self.PUM_Orientation.currentIndex()==0:
            return None
        (status, result) = Geom_Val_Plot_func(self)
        self.TW_geom.setCurrentIndex(1)
        self.TE_validation.setPlainText(result)
        if status:
            self.PB_NE.setEnabled(True)
            self.PB_NY.setEnabled(True)
        return None
    def _forwardNY(self):
        self.TE_TEXTOUT.setPlainText('Processing...')
        self.PB_NY.setEnabled(False)  #not responding.
        (inpdata,state)=LSA_SHF_Inp_func(self)
        #passing data to model for getting output printed to UI
        if state==1:
            Grm = float(self.LE_Grm2.text())
            Stmoc = float(self.LE_Stmoc.text())
            fMin = float(self.LE_nImin2.text())
            if fMin < 1e-5:
                fMin = 1e-5
                self.TE_TEXTOUT.append('Too Low nIMin detected in Nyquist plot Input. nImin is made 1e-5')
            fMax = float(self.LE_nImax2.text())
            fStep = float(self.LE_nIstep2.text())
            LSA_SHF_NY_Main_func(inpdata,Grm,Stmoc,fMin,fMax,fStep,self)
        else:
            self.TE_TEXTOUT.append('Input Error: Please check the inputs for consistency')
        self.PB_NY.setEnabled(True)
    def _forwardNE(self):
        self.TE_TEXTOUT.setPlainText('Processing...')
        self.PB_NE.setEnabled(False)  #not responding.
        (inpdata,state)=LSA_SHF_Inp_func(self)
        #passing data to model for getting output printed to UI
        if state==1:
            Grm = float(self.LE_Grm.text())
            nIMin = float(self.LE_nImin.text())
            nIMax = float(self.LE_nImax.text())
            nIStep = float(self.LE_nIstep.text())            
            StmocMin = float(self.LE_StmocMin.text())
            StmocMax = float(self.LE_StmocMax.text())
            StmocStep = float(self.LE_StmocStep.text())
            LSA_SHF_NE_Main_func(inpdata,Grm,nIMin,nIMax,nIStep,StmocMin,StmocMax,StmocStep,self)
        else:
            self.TE_TEXTOUT.append('Input Error: Please check the inputs for consistency')
        self.PB_NE.setEnabled(True)
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    uiSSA_SHF = clSSA_SHF()
    uiSSA_SHF.show()
    sys.exit(app.exec_())
    
