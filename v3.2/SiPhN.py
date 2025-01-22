# -*- coding: utf-8 -*-
"""
Created on Sat May 30 11:27:14 2020

@author: Saikrishna Nadella
"""

"""
SiPhN ::: Single Phase Natural circulation ::: (pronunced same as
'siphon')
It performs 
1. Steady state analysis (semi-analytical method)
2. Transient analysis (numerical method)
3. Stability analysis (semi-analytical linear model)
for single phase natural circulation loops.

The extent of modelling, the assumption involved, nomenclature etc can be
found in the associated M.Tech thesis in various chapters.

Developed by 
Saikrishna Nadella
as part of Masters Project in Homi Bhabha National Institute, Mumbai

"""

from PyQt5 import QtWidgets, uic, QtGui
import sys
import os
import src.TPPC
import src.SSA_SHF
import src.LSA_SHF
import src.TA

thispath = os.getcwd()
srcpath = thispath+'\\src'
UDpath = thispath+'\\UserDefined'

sys.path.append(srcpath)
sys.path.append(UDpath)

class clSIPHN(QtWidgets.QMainWindow):
    def __init__(self):
        super(clSIPHN,self).__init__()
        uic.loadUi('src\SiPhN.ui',self)
        self.PB_Continue.clicked.connect(self._forward)
    def _forward(self):
        AnType=str(self.PUM_Analysis_Type.currentText())
        if AnType == 'LSA_SHF':
            self.child2=src.LSA_SHF.clLSA_SHF()
        elif AnType == 'SSA_SHF':
            self.child2=src.SSA_SHF.clSSA_SHF()
        elif AnType == 'TA':
            self.child2=src.TA.clTA()
        self.child2.show()
        if self.CB_TPPC.isChecked():
            self.child1=src.TPPC.clTPPC()
            self.child1.show()

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    app.setWindowIcon(QtGui.QIcon('src/icon.jpg'))
    uiSIPHN = clSIPHN()
    uiSIPHN.show()
    sys.exit(app.exec_())
    
