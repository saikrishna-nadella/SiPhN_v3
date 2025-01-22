# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 16:57:23 2020

@author: saikr
"""
from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)
#from isConvex import is_convex_polygon

import numpy as np

def Geom_Val_Plot_func(uiHandle):
    status = 0
    result = ''
    #validation calculations
    try:
        #gathering inputs from ui
        L1 = float(uiHandle.LE_L1.text())
        L2 = float(uiHandle.LE_L2.text())
        L3 = float(uiHandle.LE_L3.text())
        L4 = float(uiHandle.LE_L4.text())
        Lh = float(uiHandle.LE_Lh.text())
        Lc = float(uiHandle.LE_Lc.text())
        H = float(uiHandle.LE_H.text())
        X1 = float(uiHandle.LE_X1.text())
        X2 = float(uiHandle.LE_X2.text())
        X3 = float(uiHandle.LE_X3.text())
        X4 = float(uiHandle.LE_X4.text())
        X5 = float(uiHandle.LE_X5.text())
        X6 = float(uiHandle.LE_X6.text())
        alpha1 = float(uiHandle.LE_alpha1.text())
        alpha2 = float(uiHandle.LE_alpha2.text())
        alpha3 = float(uiHandle.LE_alpha3.text())
        alpha4 = float(uiHandle.LE_alpha4.text())
        alpha5 = float(uiHandle.LE_alpha5.text())
        alpha6 = float(uiHandle.LE_alpha6.text())
        alphah = float(uiHandle.LE_alphah.text())
        alphac = float(uiHandle.LE_alphac.text())
    except (ArithmeticError, TypeError, ValueError):
        return 0, 'input error. check whether all inputs are given or not. Also check that there are no special characters in the input'
    
    #check for closedness
    LSec = np.array([L4,Lh,L1,X1,X2,X3,L2,Lc,L3,X4,X5,X6])
    angSec = np.array([alphah,alphah,alphah,alpha1,alpha2,alpha3,alphac,alphac,alphac,alpha4,alpha5,alpha6])
    angSec = np.radians(angSec)
    ZSec = LSec*np.sin(angSec)
    HSec = LSec*np.cos(angSec)
    if abs(np.sum(ZSec+HSec)) > 1e-4:
        result = 'Error: Loop is not closed'
    else:
        status=1
        result = 'Loop is closed. Analysis can be proceeded'
    #generating end points of each section
    polygon = [(0,0)]
    xdata,ydata,xhdata,yhdata,xcdata,ycdata = [],[],[],[],[],[]
    xdata.append(0), ydata.append(0)
    
    i = 0
    for LSeci,angSeci in zip(LSec,angSec):
        dx = LSeci*np.cos(angSeci)
        dy = LSeci*np.sin(angSeci)
        nextPoint = (polygon[i][0]+dx,polygon[i][1]+dy)
        if i==0 or i==1:
            xhdata.append(nextPoint[0])
            yhdata.append(nextPoint[1])
        elif i==6 or i==7:
            xcdata.append(nextPoint[0])
            ycdata.append(nextPoint[1])
        xdata.append(nextPoint[0])
        ydata.append(nextPoint[1])
        polygon.append(nextPoint)
        i=i+1

    if True:

        #plotting the geometry
 #       uiHandle.addToolBar(NavigationToolbar(uiHandle.MplWidget_geometry.canvas, uiHandle))
        
        uiHandle.MplWidget_geometry.canvas.axes.clear()
        uiHandle.MplWidget_geometry.canvas.axes.plot(xdata,ydata,color='black',marker='o',linestyle='-',linewidth=2,label='piping')
        uiHandle.MplWidget_geometry.canvas.axes.plot(xhdata,yhdata,color='red',linestyle='-',linewidth=5,label='heater')
        uiHandle.MplWidget_geometry.canvas.axes.plot(xcdata,ycdata,color='blue',linestyle='-',linewidth=5,label='cooler')
        uiHandle.MplWidget_geometry.canvas.axes.legend()
        uiHandle.MplWidget_geometry.canvas.axes.set_title('Given NCL geometry')
        uiHandle.MplWidget_geometry.canvas.draw()
        
    
    return status, result
