# -*- coding: utf-8 -*-
"""
Created on Sun May 31 12:56:41 2020

@author: saikr
"""


from PyQt5 import QtWidgets, uic
import sys
import os
thispath = os.getcwd()
srcpath = thispath+'\\src'
sys.path.append(srcpath)
from SSA_SHF_Inp import SSA_SHF_Inp_func
from SSA_SHF_SS import SSA_SHF_SS_func
from Geom_Val_Plot import Geom_Val_Plot_func
from triggerGeomInputs import triggerGeomInputs_func

class clSSA_SHF(QtWidgets.QMainWindow):
    def __init__(self):
        super(clSSA_SHF,self).__init__()
        uic.loadUi('src\SSA_SHF2.ui',self)
        
        triggerGeomInputs_func(self)
        self.PB_DIM.clicked.connect(self._forwardDIM)
        self.PB_ND.clicked.connect(self._forwardND)
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
        self.PB_DIM.setEnabled(False)
        self.PB_ND.setEnabled(False)        
        
    def _validateGeom(self):
        if self.PUM_Orientation.currentIndex()==0:
            return None
        (status, result) = Geom_Val_Plot_func(self)
        self.TW_geom.setCurrentIndex(1)
        self.TE_validation.setPlainText(result)
        if status:
            self.PB_DIM.setEnabled(True)
            self.PB_ND.setEnabled(True)
        return None
    
    def _forwardDIM(self):
        self.TE_DIMout.setPlainText('Processing...')
        #clearning ND side inputs to avoid confusion
        self.LE_Stmoh.setText('')
        self.LE_Stmohl.setText('')
        self.LE_Stmoc.setText('')
        self.LE_Stmocl.setText('')
        self.LE_Grm.setText('')
        self.LE_thetas.setText('')
        #getting input data structured 
        DIM=1
        (inpdata,state)=SSA_SHF_Inp_func(self,DIM)
        #print(inpdata)
        #print(state)
        #passing data to model for getting output printed to UI
        if state==1:
            power=float(self.LE_Q.text())
            hoc = float(self.LE_hoc.text())
            SSA_SHF_SS_func(inpdata,power,hoc,self,DIM)
    def _forwardND(self):
        self.TE_NDout.setPlainText('Processing...')
        #clearing DIM side inputs to avoid confusion
        self.LE_hoh.setText('')
        self.LE_hohl.setText('')
        self.LE_hoc.setText('')
        self.LE_hocl.setText('')
        self.LE_Q.setText('')
        #self.LE_Tamb.setText('')
        self.LE_Ts.setText('')
        #getting input data structured
        DIM=0
        (inpdata,state)=SSA_SHF_Inp_func(self,DIM)
        #passing data to model for getting output printed to UI
        if state==1:
            Grm = float(self.LE_Grm.text())
            Stmoc = float(self.LE_Stmoc.text())
            SSA_SHF_SS_func(inpdata,Grm,Stmoc,self,DIM)

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    uiSSA_SHF = clSSA_SHF()
    uiSSA_SHF.show()
    sys.exit(app.exec_())
    
