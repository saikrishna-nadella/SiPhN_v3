# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 17:33:33 2020

@author: saikr
"""


def triggerGeomInputs_func(self):
    #must needed inputs for any orientation 
    self.LE_L1.setEnabled(True)
    self.LE_L2.setEnabled(True)
    self.LE_L3.setEnabled(True)
    self.LE_L4.setEnabled(True)
    self.LE_Lh.setEnabled(True)
    self.LE_Lc.setEnabled(True)
    self.LE_D.setEnabled(True)
    self.LE_H.setEnabled(True)
    self.LE_tw.setEnabled(True)
    self.LE_alphah.setEnabled(True)
    self.LE_alphac.setEnabled(True)
    # not needed except for FreeForm, VHVC2
    self.LE_X3.setEnabled(False)
    self.LE_X6.setEnabled(False)
    self.LE_alpha3.setEnabled(False)
    self.LE_alpha6.setEnabled(False)
    i = self.PUM_Orientation.currentIndex()
    if i==0:
        self.LE_L1.setEnabled(False)
        self.LE_L2.setEnabled(False)
        self.LE_L3.setEnabled(False)
        self.LE_L4.setEnabled(False)
        self.LE_Lh.setEnabled(False)
        self.LE_Lc.setEnabled(False)
        self.LE_D.setEnabled(False)
        self.LE_H.setEnabled(False)
        self.LE_tw.setEnabled(False)
        self.LE_X1.setEnabled(False)
        self.LE_X2.setEnabled(False)
        self.LE_X4.setEnabled(False)
        self.LE_X5.setEnabled(False)
        self.LE_alphah.setEnabled(False)
        self.LE_alpha1.setEnabled(False)
        self.LE_alpha2.setEnabled(False)
        self.LE_alphac.setEnabled(False)
        self.LE_alpha4.setEnabled(False)
        self.LE_alpha5.setEnabled(False)
    elif i==1 or i==2 or i==7 or i==8:
        self.LE_X1.setEnabled(True)
        self.LE_X2.setEnabled(False)
        self.LE_X4.setEnabled(True)
        self.LE_X5.setEnabled(False)
        self.LE_alpha1.setEnabled(True)
        self.LE_alpha2.setEnabled(False)
        self.LE_alpha4.setEnabled(True)
        self.LE_alpha5.setEnabled(False)
    elif i==3 or i==6:
        self.LE_X1.setEnabled(False)
        self.LE_X2.setEnabled(False)
        self.LE_X4.setEnabled(True)
        self.LE_X5.setEnabled(True)
        self.LE_alpha1.setEnabled(False)
        self.LE_alpha2.setEnabled(False)
        self.LE_alpha4.setEnabled(True)
        self.LE_alpha5.setEnabled(True)
    elif i==4 or i==5:
        self.LE_X1.setEnabled(True)
        self.LE_X2.setEnabled(True)
        self.LE_X4.setEnabled(False)
        self.LE_X5.setEnabled(False)
        self.LE_alpha1.setEnabled(True)
        self.LE_alpha2.setEnabled(True)
        self.LE_alpha4.setEnabled(False)
        self.LE_alpha5.setEnabled(False)
    elif i==9:
        self.LE_X1.setEnabled(True)
        self.LE_X2.setEnabled(True)
        self.LE_X3.setEnabled(True)
        self.LE_X4.setEnabled(False)
        self.LE_X5.setEnabled(False)
        self.LE_alpha1.setEnabled(True)
        self.LE_alpha2.setEnabled(True)
        self.LE_alpha3.setEnabled(True)
        self.LE_alpha4.setEnabled(False)
        self.LE_alpha5.setEnabled(False)
    elif i==10:
        self.LE_X1.setEnabled(False)
        self.LE_X2.setEnabled(False)
        self.LE_X4.setEnabled(True)
        self.LE_X5.setEnabled(True)
        self.LE_X6.setEnabled(True)
        self.LE_alpha1.setEnabled(False)
        self.LE_alpha2.setEnabled(False)
        self.LE_alpha4.setEnabled(True)
        self.LE_alpha5.setEnabled(True)
        self.LE_alpha6.setEnabled(True)
    elif i==11:
        self.LE_X1.setEnabled(True)
        self.LE_X2.setEnabled(True)
        self.LE_X3.setEnabled(True)
        self.LE_X4.setEnabled(True)
        self.LE_X5.setEnabled(True)
        self.LE_X6.setEnabled(True)
        self.LE_alpha1.setEnabled(True)
        self.LE_alpha2.setEnabled(True)
        self.LE_alpha3.setEnabled(True)
        self.LE_alpha4.setEnabled(True)
        self.LE_alpha5.setEnabled(True)
        self.LE_alpha6.setEnabled(True)
    if not self.LE_X1.isEnabled():
        self.LE_X1.setText('0')
    if not self.LE_X2.isEnabled():
        self.LE_X2.setText('0')
    if not self.LE_X3.isEnabled():
        self.LE_X3.setText('0')
    if not self.LE_X4.isEnabled():
        self.LE_X4.setText('0')
    if not self.LE_X5.isEnabled():
        self.LE_X5.setText('0')
    if not self.LE_X6.isEnabled():
        self.LE_X6.setText('0')
    if not self.LE_alpha1.isEnabled():
        self.LE_alpha1.setText('0')
    if not self.LE_alpha2.isEnabled():
        self.LE_alpha2.setText('0')
    if not self.LE_alpha3.isEnabled():
        self.LE_alpha3.setText('0')
    if not self.LE_alpha4.isEnabled():
        self.LE_alpha4.setText('0')
    if not self.LE_alpha5.isEnabled():
        self.LE_alpha5.setText('0')
    if not self.LE_alpha6.isEnabled():
        self.LE_alpha6.setText('0')