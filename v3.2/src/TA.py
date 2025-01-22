# -*- coding: utf-8 -*-
"""
Created on Sun Jul 26 19:52 2020

@author: saikrishna Nadella
"""

import pandas
from PyQt5 import QtCore,QtWidgets, uic
import sys
import os
import numpy as np

thispath = os.getcwd()
srcpath = thispath+'\\src'
sys.path.append(srcpath)

from TA_Inp import TA_Inp_func
from TA_Main import TA_Main_func

class  clTA(QtWidgets.QMainWindow):
    def __init__(self):
        super(clTA,self).__init__()
        uic.loadUi('src\TA.ui',self)
        #variables
        self.toBeStopped = False

        #initial states
        self.LE_wInit.setEnabled(True)
        self.LE_tf.setEnabled(False)
        
        # connections
        self.PB_addRow.clicked.connect(self._addUIT)
        self.PB_delRow.clicked.connect(self._delUIT)
        self.PB_Start.clicked.connect(self._simulate)
        self.PB_Stop.clicked.connect(self._stop)
        self.PB_SSINIT.toggled.connect(self._INITtoggle)
        self.PB_TempDep.toggled.connect(self._TDPtoggle)
        self.PB_LoadTable.clicked.connect(self._loadUIT)
        self.PB_showGrid.clicked.connect(self._showGrid)
    def _loadUIT(self):
        csvFile = QtWidgets.QFileDialog.getOpenFileName(self,'Open File',thispath,"Data file (*.csv)")
        if not csvFile[0]:
            return None
        filePath = csvFile[0]
        # using Numpy
        df = np.genfromtxt(filePath, delimiter=',')
        nRows = df.shape[0]
        nCols = df.shape[1]
        if nCols != 11:
            self.TE_TEXT_OUT.setPlainText('Incorrect Format of CSV file')
        else:
            self.UIT.setRowCount(nRows)
            for i in range(nRows):
                for j in range(11):
                    item = QtWidgets.QTableWidgetItem(str(df[i,j]))
                    self.UIT.setItem(i,j,item)
    def _addUIT(self):
        rowIndex = self.UIT.currentRow()
        self.UIT.insertRow(rowIndex+1)
    def _delUIT(self):
        rowIndex = self.UIT.currentRow()
        self.UIT.removeRow(rowIndex)
    def _simulate(self):
        self.toBeStopped = False
        self.PB_Start.setEnabled(False)
        self.TE_TEXT_OUT.setPlainText('Reading inputs...')
        (inpdata,state) = TA_Inp_func(self)
        if state:
            self.thread1 = QThread1()
            self.thread1.dataIn(self,inpdata)
            self.thread1.start()
            self.thread1.sigProgBar.connect(self._updateProgBar)
            self.thread1.sigTOut.connect(self._updateTOut)
            self.thread1.sigMonitor.connect(self._updateMonitor)
        else:
            self.TE_TEXT_OUT.setPlainText('Input Error: Please check the inputs for consistency')
    def _INITtoggle(self):
        if self.PB_SSINIT.isChecked():
            self.LE_wInit.setEnabled(False)
            self.LE_tf.setEnabled(True)
        else:
            self.LE_wInit.setEnabled(True)
            self.LE_tf.setEnabled(False)
    def _TDPtoggle(self):
        if self.PB_TempDep.isChecked():
            self.LE_rho.setEnabled(False)
            self.LE_Cp.setEnabled(False)
            self.LE_k.setEnabled(False)
            self.LE_mu.setEnabled(False)
            self.LE_betaf.setEnabled(False)
            self.LE_rhow.setEnabled(False)
            self.LE_Cpw.setEnabled(False)
            self.LE_kw.setEnabled(False)
            self.LE_TempDep.setText('ON')
        else:
            self.LE_rho.setEnabled(True)
            self.LE_Cp.setEnabled(True)
            self.LE_k.setEnabled(True)
            self.LE_mu.setEnabled(True)
            self.LE_betaf.setEnabled(True)
            self.LE_rhow.setEnabled(True)
            self.LE_Cpw.setEnabled(True)
            self.LE_kw.setEnabled(True)
            self.LE_TempDep.setText('OFF')
    def _updateProgBar(self,progress):
        self.progBar.setValue(progress)
    def _updateTOut(self,info):
        self.TE_TEXT_OUT.append(str(info))
    def _updateMonitor(self,flag,data):
        #define
        if flag:
        #    if axs not created: 
            #uiHandle.addToolBar(NavigationToolbar(uiHandle.MplWidget_monitor.canvas, uiHandle))
            self.MplWidget_monitor.canvas.axes.clear()
            Ph = self.MplWidget_monitor.canvas.axes.plot(data[0],data[1],color='blue',marker='',linestyle='-',linewidth=2)[0]#plot handle
            self.MplWidget_monitor.canvas.axes.grid()
            #self.MplWidget_monitor.canvas.axes.set_title('Monitor-1')
            self.MplWidget_monitor.canvas.axes.set_xlabel('Time (s)')
            self.MplWidget_monitor.canvas.axes.set_ylabel('Mass flow rate (kg/s)')
            self.MplWidget_monitor.canvas.draw()
            '''else:
                Ph.set_xdata(np.arange(0,t,dtUpdate))
                Ph.set_ydata(wup)
                uiHandle.MplWidget_monitor.canvas.draw()'''
    def _stop(self):
        self.toBeStopped = True

    def _showGrid(self):
        #reading geometry & grid parameters
        secNs = self.UIT.rowCount()
        inpMesh = np.zeros((secNs,4))
        for sec in range(secNs):
            for col in range(4):
                inpMesh[sec,col]=float(self.UIT.item(sec,col).text())
        secType = inpMesh[:,0].astype('int32')
        nSec = inpMesh[:,2].astype('int32')
        if max(nSec%2)!=0:
            nSec = nSec%2+nSec
            self.TE_TEXT_OUT.append('''\nOdd number of nodes found in one or more sections.\n The node count is raised by one to make even.''')
        LSec = inpMesh[:,1]
        angSec = inpMesh[:,3]*np.pi/180

        #check loop closedness
        ZSec = LSec*np.sin(angSec)
        HSec = LSec*np.cos(angSec)
        if sum(ZSec+HSec) > 0.001:
            self.TE_TEXT_OUT.append('''Error: Loop is not closed.''')
            return None
        #generating points for plot
        polygon = [(0,0)]
        xdata,ydata,xhdata,yhdata,xcdata,ycdata = [],[],[],[],[],[]
        xdata.append(0), ydata.append(0)
        j = 0
        k = 0
        for LSeci,angSeci,nSeci in zip(LSec,angSec,nSec):
            dl = LSeci/nSeci
            for i in range(nSeci):
                dx = dl*np.cos(angSeci)
                dy = dl*np.sin(angSeci)
                nextPoint = (polygon[k][0]+dx,polygon[k][1]+dy)
                if secType[j]==1:
                    if i==0:
                        xhdata.append(polygon[k][0])
                        yhdata.append(polygon[k][1])
                    xhdata.append(nextPoint[0])
                    yhdata.append(nextPoint[1])
                elif secType[j]==2:
                    if i==0:
                        xcdata.append(polygon[k][0])
                        ycdata.append(polygon[k][1])
                    xcdata.append(nextPoint[0])
                    ycdata.append(nextPoint[1])
                xdata.append(nextPoint[0])
                ydata.append(nextPoint[1])
                polygon.append(nextPoint)
                k = k + 1
            j=j+1
            
        #plotting grid   
        self.MplWidget_Grid.canvas.axes.clear()
        self.MplWidget_Grid.canvas.axes.plot(xdata,ydata,color='black',marker='o',linestyle='-',linewidth=2,label='piping')
        self.MplWidget_Grid.canvas.axes.plot(xhdata,yhdata,color='red',marker='o',linestyle='-',linewidth=5,label='heater')
        self.MplWidget_Grid.canvas.axes.plot(xcdata,ycdata,color='blue',marker='o',linestyle='-',linewidth=5,label='cooler')
        self.MplWidget_Grid.canvas.axes.legend()
        #self.MplWidget_Grid.canvas.axes.set_title('Given NCL geometry')
        self.MplWidget_Grid.canvas.draw()
        self.TW2.setCurrentIndex(0)
class QThread1(QtCore.QThread):
    sigProgBar = QtCore.pyqtSignal(int)
    sigTOut = QtCore.pyqtSignal(str)
    sigMonitor = QtCore.pyqtSignal(int,np.ndarray)
    def __init__(self, parent=None):
        QtCore.QThread.__init__(self,parent)

    def dataIn(self,uiHandle,data):
        self.uiHandle = uiHandle
        self.inpdata = data
        
    def run(self):
        status = TA_Main_func(self,self.uiHandle,self.inpdata)
        self.uiHandle.PB_Start.setEnabled(True)
        return None

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    uiTA = clTA()
    uiTA.show()
    sys.exit(app.exec_())
