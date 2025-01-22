# -*- coding: utf-8 -*-
"""
Created on Sat May 30 11:27:14 2020

@author: Saikrishna Nadella
"""

from PyQt5 import QtWidgets, uic
import sys
import os
thispath = os.getcwd()
UDpath = thispath+'\\UserDefined'
sys.path.append(UDpath)

from userTPPC import ThermoPhysicalProperties

# class for UI
class clTPPC(QtWidgets.QMainWindow):
    def __init__(self):
        super(clTPPC,self).__init__()
        uic.loadUi('src/TPPC.ui',self)
        self.PB_Calc.clicked.connect(self.TPPCcalc)
    def TPPCcalc(self):
        Temp = float(self.ET_temp.toPlainText())
        if self.RB_F.isChecked():
            flag=1
        elif self.RB_W.isChecked():
            flag=2
        else:
            flag=3
            
        (a,b)=ThermoPhysicalProperties(Temp, flag)
        self.ET_out.setPlainText(str(a))
        if flag==1 or flag==3:
            self.ET_rho.setPlainText(str(b[0]))
            self.ET_Cp.setPlainText(str(b[1]))
            self.ET_k.setPlainText(str(b[2]))
            self.ET_mu.setPlainText(str(b[3]))
            self.ET_betaf.setPlainText(str(b[4]))
            if flag ==1:
                self.ET_rhow.setPlainText('')
                self.ET_Cpw.setPlainText('')
                self.ET_kw.setPlainText('')
        if flag==2 or flag==3:
            self.ET_rhow.setPlainText(str(b[5]))
            self.ET_Cpw.setPlainText(str(b[6]))
            self.ET_kw.setPlainText(str(b[7]))
            if flag==2:
                self.ET_rho.setPlainText('')
                self.ET_Cp.setPlainText('')
                self.ET_k.setPlainText('')
                self.ET_mu.setPlainText('')
                self.ET_betaf.setPlainText('')               
# main function to create instance of application    
def main():
    # create instance of application
    app = QtWidgets.QApplication(sys.argv)
    # show GUI of application
    uiTPPC = clTPPC()
    uiTPPC.show()
    # execute the application's main loop
    sys.exit(app.exec_())
  
if __name__ == '__main__':
    main()
    
    
