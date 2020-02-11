# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 21:19:43 2015

@author: Ultracold
"""

from PyQt5.uic import loadUiType
import numpy as np
import sys
import os

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
from PyQt5 import QtGui, QtCore  
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
import time
from scipy.optimize import curve_fit
from lmfit import Model
from lmfit.models import GaussianModel
    
ROOT_PATH = os.getcwd()

Ui_MainWindow, QMainWindow = loadUiType(os.path.join(ROOT_PATH, 'AnalysisWindow6.ui'))

# *********************************************************   
def funcLinear2(x,m,b):
# the independent variable must be sent first for the fit function
     return m*x+b    
# *********************************************************   
  
# *********************************************************   
def funcGaussian2(x,A,x0,sigma,y0):
# the independent variable must be sent first for the fit function
       return A*np.exp((-((x-x0)/sigma)**2)/2)+y0
# *********************************************************   

# *********************************************************   
def funcExponential2(x,A,x0,tau,y0):
# the independent variable must be sent first for the fit function
       return A*np.exp(-(x-x0)/tau)+y0    
# *********************************************************   

# *********************************************************   
def funcTOF2(x,sig0,v):
# the independent variable must be sent first for the fit function
       return np.sqrt(sig0**2 + (v*x)**2)    
# *********************************************************   


class Main(QMainWindow, Ui_MainWindow):
    def __init__(self, ):
        super(Main, self).__init__()
        self.setupUi(self)
       
        # Make the figure canvas, including the toolbar (top)
        self.fig1 = Figure()
        self.ax1f1 = self.fig1.add_subplot(111)
        self.canvas = FigureCanvas(self.fig1)
        self.mplvl.addWidget(self.canvas)
        self.canvas.draw()
        self.toolbar = NavigationToolbar(self.canvas, self.mplwindow, coordinates=True)
        self.toolbar.setFixedHeight(24)
        self.mplvl.addWidget(self.toolbar)
        
        # Make the figure canvas, including the toolbar (bottom)
        self.fig2 = Figure()
        self.ax1f2 = self.fig2.add_subplot(111)
        self.canvas2 = FigureCanvas(self.fig2)
        self.mplvl2.addWidget(self.canvas2)
        self.canvas2.draw()
        self.toolbar2 = NavigationToolbar(self.canvas2, self.mplwindow2, coordinates=True)
        self.toolbar2.setFixedHeight(24)
        self.mplvl2.addWidget(self.toolbar2)
       
        # Set up buttons
        self.ClrData.clicked.connect(self.clear_data)
        self.ExportData.clicked.connect(self.export_data)
        self.ExportAllData.clicked.connect(self.export_alldata)
        self.searchDataFile.clicked.connect(self.load_Datapath)
#        self.searchMmtFile.clicked.connect(self.load_Mmtpath)
        self.CollectData.setChecked(True)
        self.settingsBrowse.clicked.connect(self.load_settings)
        self.settingsSave.clicked.connect(self.save_settings)
        self.settingsApply.clicked.connect(self.apply_settings)
        self.AutoScale1.setChecked(True)
        self.AutoScale2.setChecked(True)
        self.FitOnButton.toggled.connect(self.update_figure)
        self.FitOnButton2.toggled.connect(self.update_figure2)
        
#        # default files        
        self.DataFileInput.setText("Z:\\_tof-analysis\\AllData.txt")
        self.SettingsFileLE.setText(("Z:\\_tof-analysis\\testSettings"))
        # set original data array to zero
        self.xvals=[]
        self.yvals=[]      
        self.yvals2=[]      
        self.CtrlData = []
        self.MmtData = []
#        self.initial_data()
        # and fit parameters
        self.fitParam1.setText(str(1.0))
        self.fitParam2.setText(str(1.0))
        self.fitParam3.setText(str(1.0))
        self.fitParam4.setText(str(1.0))
        self.meanLE.setText("n/a")
        self.stdLE.setText("n/a")
        self.fit2Param1.setText(str(1.0))
        self.fit2Param2.setText(str(1.0))
        self.fit2Param3.setText(str(1.0))
        self.fit2Param4.setText(str(1.0))
        self.mean2LE.setText("n/a")
        self.std2LE.setText("n/a")
        
        # Get the info for the dropdown combo boxes
#        listitems = QtCore.QString('hi hello')
#        self.get_combo()
#        self.get_comboMmt() 
        self.get_variableNames()
        self.comboControl.addItems(self.ctrlList)
        self.comboMmt.addItems(self.MmtList)
        self.comboMmt2.addItems(self.Mmt2List)
        self.updatelists.clicked.connect(self.funcSetLists)
        
       # Setup combo boxes to update things when the values change
        self.comboMmt.currentIndexChanged.connect(self.comboMmt_update)
        self.comboMmt2.currentIndexChanged.connect(self.comboMmt2_update)
        self.comboControl.currentIndexChanged.connect(self.comboCtrl_update)
        
        # Make the data table
        self.mpldata.setHorizontalHeaderLabels([self.comboControl.currentText(), self.comboMmt.currentText(), self.comboMmt2.currentText()])
        self.ax1f1.set_xlabel(self.comboControl.currentText())
        self.ax1f1.set_ylabel(self.comboMmt.currentText())
        # add a context menu for table so we can delete elements
        self.mpldata.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.mpldata.customContextMenuRequested.connect(self.dataTableMenu)

        # Set the info for the fitter 
        self.comboFitType.addItems(['Linear','Exponential','Gaussian','TOF'])
        self.comboFit2Type.addItems(['Linear','Exponential','Gaussian','TOF'])
        # Setup combo box to update things when the values change
        self.comboFitType.currentIndexChanged.connect(self.comboFitType_update)
        self.comboFit2Type.currentIndexChanged.connect(self.comboFit2Type_update)
        self.equationLabel.setText('y = mx+b')  
        self.lblFitParam1.setText('Slope, m')
        self.lblFitParam2.setText('Intcpt, b')
        self.lblFitParam3.setText('')
        self.lblFitParam4.setText('')
        self.equationLabel2.setText('y = mx+b')  
        self.lblFit2Param1.setText('Slope, m')
        self.lblFit2Param2.setText('Intcpt, b')
        self.lblFit2Param3.setText('')
        self.lblFit2Param4.setText('')
      # Set up Update Fit button
        self.UpdateFit.clicked.connect(self.update_fit_callback)
        self.UpdateFit2.clicked.connect(self.update_fit2_callback)
    
    # instead of constant queries, wait for seeing image file update
        # trigger when the image data updates         
        self.filename = [self.DataFileInput.text()]
        self.fileinfo = QtCore.QFileInfo(self.filename[0])
        #self.fileChange = self.fileinfo.lastModified
        
#        self.fs_watcher = QtCore.QFileSystemWatcher(filename)
        self.last_change = time.time()
        
        self.fs_watcher = QtCore.QFileSystemWatcher()
        self.fs_watcher.addPaths(self.filename)
        
        self.fs_watcher.fileChanged.connect(self.file_changed)


    # *********************************************************   
    def file_changed(self):
        
        if time.time() - self.last_change < 0.5:
            return
            
        self.last_change = time.time()

        if self.CollectData.isChecked():
            self.update_data()
            self.update_figure()
            self.update_figure2()
            #self.fileChange = self.fileinfo.lastModified
            if self.AutoFitOn.isChecked(): self.update_fitNew(1)
            if self.AutoFitOn2.isChecked(): self.update_fitNew(2)
            
            
    def dataTableMenu(self, pos):
        globalPos = self.mapToGlobal(pos)
        menu = QMenu()
        delRow = menu.addAction("Delete row")
        cellZero = menu.addAction("Set cell to zero")
#        selectedItem = menu.exec_(globalPos)
        selectedItem = menu.exec_(QCursor.pos())
        if selectedItem:
            row = self.mpldata.currentRow()
            col = self.mpldata.currentColumn()
            if selectedItem == delRow:
                row = self.mpldata.currentRow()
                print("selected: deleted", row,col)
                self.tableDeleteRow(row)
            if selectedItem == cellZero:
                print("selected: set zero",row,col)
                self.tableCellZero(row,col)
    # ****** end dataTableMenu ********************************   

    def tableDeleteRow(self, row):
        if row <= self.xvals.__len__():
            del self.xvals[row]
            del self.yvals[row]
            del self.yvals2[row]
            del self.CtrlData[row]
            del self.MmtData[row]
    
        # put current data in table, set last row to zero
        for i in range(len(self.xvals)):
            item = QTableWidgetItem(str(self.xvals[i]))
            self.mpldata.setItem(i,0,item)
            item = QTableWidgetItem(str(self.yvals[i]))
            self.mpldata.setItem(i,1,item)
            item = QTableWidgetItem(str(self.yvals2[i]))
            self.mpldata.setItem(i,2,item)
        item = QTableWidgetItem('')
        self.mpldata.setItem(i+1,0,item)
        item = QTableWidgetItem('')
        self.mpldata.setItem(i+1,1,item) 
        item = QTableWidgetItem('')
        self.mpldata.setItem(i+1,2,item) 
        
        self.update_figure()
        self.update_figure2()
        if self.AutoFitOn.isChecked(): self.update_fitNew(1)
    
    # ****** end tableDeleteRow ********************************

    def tableCellZero(self, row,col):
        if row <= self.xvals.__len__():
            if col == 0:
                self.xvals[row] = 0
                item = QTableWidgetItem(str(self.xvals[row]))
                self.mpldata.setItem(row,0,item)
                self.CtrlData[row][self.comboControl.currentIndex()]=0
            elif col == 1:
                self.yvals[row] = 0
                item = QTableWidgetItem(str(self.yvals[row]))
                self.mpldata.setItem(row,1,item)
                self.MmtData[row][self.comboMmt.currentIndex()] =0 
         
        self.update_figure()    
        if self.AutoFitOn.isChecked(): self.update_fitNew(1)
    
    # ****** end tableDeleteRow ********************************   

    def update_figure(self):
        # a routine to update the figure with new data
        # get plot settings before update
        xaxisdata = self.ax1f1.get_xlim()        
        yaxisdata = self.ax1f1.get_ylim()
        xlabeldata = self.ax1f1.get_xlabel()
        ylabeldata = self.ax1f1.get_ylabel()
        titledata = self.ax1f1.get_title()
        self.xvals_f = np.array(self.xvals).astype(float)
        self.yvals_f = np.array(self.yvals).astype(float)
        self.ax1f1.clear()
        self.ax1f1.plot(self.xvals_f,self.yvals_f,'ko',ms=10,markerfacecolor='g')
               
        if (self.FitOnButton.isChecked()):
            fitx = np.linspace(xaxisdata[0],xaxisdata[1], 100)
            fity = np.ones_like(fitx)
          
            if (str(self.comboFitType.currentText()) == 'Linear'):
                fity = self.funcLinear(fitx,float(self.fitParam1.text()),float(self.fitParam2.text()))
            
            elif (str(self.comboFitType.currentText()) == 'Exponential'):
                fity = self.funcExponential(fitx,float(self.fitParam1.text()),float(self.fitParam2.text()),float(self.fitParam3.text()),float(self.fitParam4.text()))
            
            elif (str(self.comboFitType.currentText()) == 'Gaussian'):
                fity = self.funcGaussian(fitx,float(self.fitParam1.text()),float(self.fitParam2.text()),float(self.fitParam3.text()),float(self.fitParam4.text()))
           
            elif (str(self.comboFitType.currentText()) == 'TOF'):
                fity = self.funcTOF(fitx,float(self.fitParam1.text()),float(self.fitParam2.text()))
            self.ax1f1.plot(fitx,fity,'k--',lw=2)
        
        # reset some of the axis properties so they don't update with new data
        if self.AutoScale1.isChecked() and len(self.xvals) > 0:
            xmin, xmax = np.nanmin(self.xvals_f), np.nanmax(self.xvals_f)
            ymin, ymax = np.nanmin(self.yvals_f), np.nanmax(self.yvals_f)
            if len(self.xvals) == 1:
                min_width, min_height = xaxisdata[1]-xaxisdata[0], yaxisdata[1]-yaxisdata[0]
            else:
                min_width, min_height = xmax - xmin, ymax - ymin
            self.ax1f1.set_xlim(xmin - 0.1*min_width, xmax + 0.1*min_width)
            self.ax1f1.set_ylim(ymin - 0.1*min_height, ymax + 0.1*min_height)
        else:
            self.ax1f1.set_xlim(xaxisdata[0],xaxisdata[1])
            self.ax1f1.set_ylim(yaxisdata[0],yaxisdata[1])
        self.ax1f1.set_xlabel(xlabeldata)
        self.ax1f1.set_ylabel(ylabeldata)
        self.ax1f1.set_title(titledata)
       
      
#        self.ax1f1.sca = self.style
        self.canvas.draw()
    # ***end update_figure **************************************
        
    def update_figure2(self):
        # a routine to update the figure with new data
        # get plot settings before update
        xaxisdata = self.ax1f2.get_xlim()        
        yaxisdata = self.ax1f2.get_ylim()
        xlabeldata = self.ax1f2.get_xlabel()
        ylabeldata = self.ax1f2.get_ylabel()
        titledata = self.ax1f2.get_title()    
        self.xvals_f = np.array(self.xvals).astype(float)
        self.yvals2_f = np.array(self.yvals2).astype(float)
        self.ax1f2.clear()
        self.ax1f2.plot(self.xvals_f,self.yvals2_f,'ko',ms=10,markerfacecolor='b')
               
        if (self.FitOnButton2.isChecked()):
            fitx = np.linspace(xaxisdata[0],xaxisdata[1])
            fity = np.ones_like(fitx)
          
            if (str(self.comboFit2Type.currentText()) == 'Linear'):
                fity = self.funcLinear(fitx,float(self.fit2Param1.text()),float(self.fit2Param2.text()))
            
            elif (str(self.comboFit2Type.currentText()) == 'Exponential'):
                fity = self.funcExponential(fitx,float(self.fit2Param1.text()),float(self.fit2Param2.text()),float(self.fit2Param3.text()),float(self.fit2Param4.text()))
            
            elif (str(self.comboFit2Type.currentText()) == 'Gaussian'):
                fity = self.funcGaussian(fitx,float(self.fit2Param1.text()),float(self.fit2Param2.text()),float(self.fit2Param3.text()),float(self.fit2Param4.text()))
           
            elif (str(self.comboFit2Type.currentText()) == 'TOF'):
                fity = self.funcTOF(fitx,float(self.fit2Param1.text()),float(self.fit2Param2.text()))
           
            self.ax1f2.plot(fitx,fity,'k--',lw=2)
        
        # reset some of the axis properties so they don't update with new data
        if self.AutoScale2.isChecked() and len(self.xvals) > 0:
            xmin, xmax = np.nanmin(self.xvals_f), np.nanmax(self.xvals_f)
            ymin, ymax = np.nanmin(self.yvals2_f), np.nanmax(self.yvals2_f)
            if len(self.xvals) == 1:
                min_width, min_height = xaxisdata[1]-xaxisdata[0], yaxisdata[1]-yaxisdata[0]
            else:
                min_width, min_height = xmax - xmin, ymax - ymin
            self.ax1f2.set_xlim(xmin - 0.1*min_width, xmax + 0.1*min_width)
            self.ax1f2.set_ylim(ymin - 0.1*min_height, ymax + 0.1*min_height)
        else:
            self.ax1f2.set_xlim(xaxisdata[0],xaxisdata[1])
            self.ax1f2.set_ylim(yaxisdata[0],yaxisdata[1])
            
        self.ax1f2.set_xlabel(xlabeldata)
        self.ax1f2.set_ylabel(ylabeldata)
        self.ax1f1.set_title(titledata)
           
#        self.ax1f1.sca = self.style
        self.canvas2.draw()
    # ***end update_figure **************************************
        
    def initial_data(self):
        # a routine to gather data from various files and compile it nicely
        
        self.mpldata.setHorizontalHeaderLabels([self.comboControl.currentText(), self.comboMmt.currentText(),self.comboMmt2.currentText()])

        allnewx = self.getData()
        self.xvals.append(allnewx[self.comboControl.currentIndex()])
        self.CtrlData.append(allnewx)
        self.yvals.append(allnewx[self.comboMmt.currentIndex()])
        self.yvals2.append(allnewx[self.comboMmt2.currentIndex()])
        self.MmtData.append(allnewx)
      
        item = QTableWidgetItem(str(self.xvals[0]))
        self.mpldata.setItem(0,0,item)
        item = QTableWidgetItem(str(self.yvals[0]))
        self.mpldata.setItem(0,1,item)
        item = QTableWidgetItem(str(self.yvals2[0]))
        self.mpldata.setItem(0,2,item)           
    # ***end intial_data **************************************
            
    def update_data(self):
        # a routine to gather data from various files and compile it nicely
        self.mpldata.setHorizontalHeaderLabels([self.comboControl.currentText(), self.comboMmt.currentText(), self.comboMmt2.currentText()])
#        allnewx= self.getCtrl()  
        allnewx = self.getData()
        self.xvals.append(allnewx[self.comboControl.currentIndex()])
        self.CtrlData.append(allnewx)
        
#        allnewy= self.getMmt()   
        self.yvals.append(allnewx[self.comboMmt.currentIndex()])
        self.MmtData.append(allnewx)
        self.yvals2.append(allnewx[self.comboMmt2.currentIndex()])
        
        # put current data in table
        ind = self.xvals.__len__()
        for i in range(len(self.xvals)):
            item = QTableWidgetItem(str(self.xvals[i]))
            self.mpldata.setItem(i,0,item)
            item = QTableWidgetItem(str(self.yvals[i]))
            self.mpldata.setItem(i,1,item)
            item = QTableWidgetItem(str(self.yvals2[i]))
            self.mpldata.setItem(i,2,item)
        
        # display mean and standard deviation of y-values
        yv = [float(i) for i in self.yvals] 
        self.meanLE.setText(str(np.mean(yv)))
        self.stdLE.setText(str(np.std(yv)))
        yv2 = [float(i) for i in self.yvals2] 
        self.mean2LE.setText(str(np.mean(yv2)))
        self.std2LE.setText(str(np.std(yv2)))
                   
    # ***end update_data **************************************
        
    def update_fit_callback(self):
        self.update_fitNew(1)
    
    def update_fit2_callback(self):
        self.update_fitNew(2)
    
    # *********************************************************   
    def update_fitNew(self,plot_num):
        if plot_num == 1:
            plot_check = self.FitOnButton.isChecked()
        else:
            plot_check = self.FitOnButton2.isChecked()
        if plot_check and (self.xvals.__len__() > 4):
                       
            xv = [float(i) for i in self.xvals]
            
            if plot_num == 1:
                yv = [float(i) for i in self.yvals]
                fit_params = [    self.fitParam1.text(), self.fitParam2.text(), 
                                self.fitParam3.text(), self.fitParam4.text()]
                                
                combo_fit_selection = str(self.comboFitType.currentText())
                
                fixed_params = [self.fixFitParam1.isChecked(), self.fixFitParam2.isChecked(),
                                self.fixFitParam3.isChecked(), self.fixFitParam4.isChecked()]
            else:
                yv = [float(i) for i in self.yvals2]
                fit_params = [    self.fit2Param1.text(), self.fit2Param2.text(), 
                                self.fit2Param3.text(), self.fit2Param4.text()]
                                
                combo_fit_selection = str(self.comboFit2Type.currentText())
                
                fixed_params = [self.fixFit2Param1.isChecked(), self.fixFit2Param2.isChecked(),
                                self.fixFit2Param3.isChecked(), self.fixFit2Param4.isChecked()]
            
            
            if combo_fit_selection == 'Linear':
                gmod = Model(funcLinear2)
                minit = float(fit_params[0])
                binit = float(fit_params[1])
                params = gmod.make_params(m=minit,b=binit)
                
                params['m'].vary = not fixed_params[0]

                params['b'].vary = not fixed_params[1]
#                 
                FitResults = gmod.fit(yv,x=xv,params=params)
                
                text1 = str(FitResults.best_values.get('m'))
                text2 = str(FitResults.best_values.get('b'))
                text3 = str(1.0)
                text4 = str(1.0)
        
            elif combo_fit_selection == 'Gaussian': #x,A,x0,sigma,y0
                
                gmod = Model(funcGaussian2)

                params = gmod.make_params(    A = float(fit_params[0]), x0=float(fit_params[1]),
                                            sigma=float(fit_params[2]),y0=float(fit_params[3]))
                
                params['A'].vary = not fixed_params[0]

                params['x0'].vary = not fixed_params[1]
#               
                params['sigma'].vary = not fixed_params[2]
                
                params['y0'].vary = not fixed_params[3]
#               
                FitResults = gmod.fit(yv,x=xv,params=params)
##                 
                text1 = str(FitResults.best_values.get('A'))
                text2 = str(FitResults.best_values.get('x0'))
                
                #############################################################################################
                text3 = str(FitResults.best_values.get('sigma')) # Needs to be fixed; always returns as 1.0 #
                #############################################################################################
                
                text4 = str(FitResults.best_values.get('y0')) 
                
            
            elif combo_fit_selection == 'Exponential':#x,A,x0,tau,y0
                gmod = Model(funcExponential2)

                params = gmod.make_params(    A = float(fit_params[0]), x0=float(fit_params[1]),
                                            tau=float(fit_params[2]),y0=float(fit_params[3]))

                params['A'].vary = not fixed_params[0]
                
                params['x0'].vary = not fixed_params[1]
    #               
                params['tau'].vary = not fixed_params[2]
                
                params['y0'].vary = not fixed_params[3]
    #                
                try:
                    FitResults = gmod.fit(yv,x=xv,params=params)
                except:
                    print('Fit did not work.')
    ##             
                text1 = str(FitResults.best_values.get('A'))
                text2 = str(FitResults.best_values.get('x0'))
                text3 = str(FitResults.best_values.get('tau'))
                text4 = str(FitResults.best_values.get('y0'))
            
            if combo_fit_selection == 'TOF':
                gmod = Model(funcTOF2)
                siginit = float(fit_params[0])
                vinit = float(fit_params[1])
                params = gmod.make_params(sig0=siginit,v=vinit)
                
                params['sig0'].vary = not fixed_params[0]
                
                params['v'].vary = not fixed_params[1]
#                 
                FitResults = gmod.fit(yv,x=xv,params=params)
##             
                text1 = str(FitResults.best_values.get('sig0'))
                text2 = str(FitResults.best_values.get('v'))
                text3 = str(1.0)
                text4 = str(1.0)
                
            if plot_num == 1:
                self.fitParam1.setText(text1)
                self.fitParam2.setText(text2)
                self.fitParam3.setText(text3)
                self.fitParam4.setText(text4)
                self.update_figure()
            else:
                self.fit2Param1.setText(text1)
                self.fit2Param2.setText(text2)
                self.fit2Param3.setText(text3)
                self.fit2Param4.setText(text4)
                self.update_figure2()
            
            
    # ***end update_fit **************************************    
                        
    
    # *********************************************************   
    def funcLinear(self,x,m,b):
    # the independent variable must be sent first for the fit function
           return m*x+b    
    # *********************************************************   
    
    # *********************************************************   
    def funcExponential(self,x,A,x0,tau,y0):
    # the independent variable must be sent first for the fit function
           return A*np.exp(-(x-x0)/tau)+y0    
    # *********************************************************   

    # *********************************************************   
    def funcGaussian(self,x,A,x0,sigma,y0):
    # the independent variable must be sent first for the fit function
           return A*np.exp(-((x-x0)/sigma)**2/2)+y0
    # *********************************************************   

    # *********************************************************   
    def funcTOF(self,x,sig0,v):
    # the independent variable must be sent first for the fit function
            return np.sqrt(sig0**2 + (v*x)**2)    
    # *********************************************************   

    
    # *********************************************************   
    def funcSetLists(self):
#        self.get_combo()
#        self.get_comboMmt() 
        # clear data first
        print('start clearing')
        self.ctrlList = []
        self.ctrlVals = []
        self.MmtList = []
        self.MmtVals = []
        self.MmtData = []
        self.Mmt2Vals = []
        self.Mmt2Data = []
        self.CtrlData = []
        self.xvals = []
        self.yvals = []
        self.yvals2 = []
        self.comboControl.clear()
        self.comboMmt.clear()
        self.comboMmt2.clear()
        self.mpldata.clear()
        print('end clearing')
        print(self.ctrlList)
        self.get_variableNames()
        self.comboControl.addItems(self.ctrlList)
        self.comboMmt.addItems(self.MmtList)
        self.comboMmt2.addItems(self.Mmt2List)
    # *********************************************************   
    
    # *********************************************************   
    def get_variableNames(self):
    # get the values that we select to plot from various files, and create a data file with them        
        filename = self.DataFileInput.text()
        print( "Filename: ", filename)
        fp = open(filename,'r')
        lines = fp.read().split('\n')
        fp.close()
        filenum = lines[0].split('=')[1]
        print filenum
        while lines[0] != '# ':
            del lines[0]
        del lines[0]
        endofdata = lines.index("# Image Data")
        for i in range(len(lines)):
            if lines[i][0] != '#':
                pass
            else:
                del lines[i:endofdata+1]
                break
        self.ctrlList = []
        self.ctrlVals = []
        self.MmtList = []
        self.MmtVals = []
        self.Mmt2List = []
        self.Mmt2Vals = []
        for i in lines:
            #print i
            try:
                self.ctrlList.append(i.split('=')[0])
                self.ctrlVals.append(i.split('=')[1])
            except:
                print(i)
                raise(ValueError)
                
            self.MmtList.append(i.split('=')[0])
            self.MmtVals.append(i.split('=')[1])
            self.Mmt2List.append(i.split('=')[0])
            self.Mmt2Vals.append(i.split('=')[1])
        self.ctrlList.append('Filenum')
        self.ctrlVals.append(filenum)
        self.MmtList.append('Filenum')
        self.MmtVals.append(filenum)
        self.Mmt2List.append('Filenum')
        self.Mmt2Vals.append(filenum)        
    # *********************************************************   
        
        
    # *********************************************************   
    def comboMmt_update(self):
        # reset labels on plot
        self.ax1f1.set_ylabel(self.comboMmt.currentText())
        #reset data label on table
        self.mpldata.setHorizontalHeaderLabels([self.comboControl.currentText(), self.comboMmt.currentText(), self.comboMmt2.currentText()])

       # put current data in table
        print 
        for i in range(len(self.yvals)):
            temp = self.MmtData[i]
            self.yvals[i] = temp[self.comboMmt.currentIndex()]
            item = QTableWidgetItem(str(self.yvals[i]))
            self.mpldata.setItem(i,1,item)
        self.update_figure()
    # *********************************************************   
   
   # *********************************************************   
    def comboMmt2_update(self):
        # reset labels on plot
        self.ax1f2.set_ylabel(self.comboMmt2.currentText())
        #reset data label on table
        self.mpldata.setHorizontalHeaderLabels([self.comboControl.currentText(), self.comboMmt.currentText(),self.comboMmt2.currentText()])

       # put current data in table
        print 
        for i in range(len(self.yvals)):
            temp = self.MmtData[i]
            self.yvals2[i] = temp[self.comboMmt2.currentIndex()]
            item = QTableWidgetItem(str(self.yvals2[i]))
            self.mpldata.setItem(i,2,item)
        self.update_figure2()
    # *********************************************************   
   
    # *********************************************************   
    def comboCtrl_update(self):
       # reset labels on plot
       self.ax1f1.set_xlabel(self.comboControl.currentText())
       self.ax1f2.set_xlabel(self.comboControl.currentText())
       #reset data label on table
       self.mpldata.setHorizontalHeaderLabels([self.comboControl.currentText(), self.comboMmt.currentText()])
       # put current data in table
       for i in range(len(self.yvals)):
            temp = self.CtrlData[i]
            self.xvals[i] = temp[self.comboControl.currentIndex()]
            item = QTableWidgetItem(str(self.xvals[i]))
            self.mpldata.setItem(i,0,item)
            
       self.update_figure()
       self.update_figure2()
    # *********************************************************   
      
    # *********************************************************   
    def comboFitType_update(self):
       # set up parameter labels
       parameterLabels1 = {'Linear':'Slope, m', 'Gaussian':'Ampl, A','Exponential':'Ampl,: A','TOF':'sig0'}            
       parameterLabels2 = {'Linear':'Int, b', 'Gaussian':'Centre, x_0','Exponential':'x Offset, x0','TOF':'v'}            
       parameterLabels3 = {'Linear':'', 'Gaussian':'Width, sigma','Exponential':'tau','TOF':''}            
       parameterLabels4 = {'Linear':'', 'Gaussian':'Offset, y0','Exponential':'y Offset, y0','TOF':''}            
       eqnlabel = {'Linear':'y = mx+b', 'Gaussian':'y = A exp[-(x-x0)^2/(2 sigma^2)]+y0','Exponential':'y = A exp(-(x-x0)/tau) + y0','TOF':'sqrt(sig0^2 + (vt)^2)'}            

       self.lblFitParam1.setText(parameterLabels1.get(str(self.comboFitType.currentText())))
       self.lblFitParam2.setText(parameterLabels2.get(str(self.comboFitType.currentText())))
       self.lblFitParam3.setText(parameterLabels3.get(str(self.comboFitType.currentText())))
       self.lblFitParam4.setText(parameterLabels4.get(str(self.comboFitType.currentText())))
       self.equationLabel.setText(eqnlabel.get(str(self.comboFitType.currentText())))  
   # *********************************************************       
    
   # *********************************************************   
    def comboFit2Type_update(self):
       # set up parameter labels
       parameterLabels1 = {'Linear':'Slope, m', 'Gaussian':'Ampl, A','Exponential':'Ampl,: A','TOF':'sig0'}            
       parameterLabels2 = {'Linear':'Int, b', 'Gaussian':'Centre, x_0','Exponential':'x Offset, x0','TOF':'v'}            
       parameterLabels3 = {'Linear':'', 'Gaussian':'Width, sigma','Exponential':'tau','TOF':''}            
       parameterLabels4 = {'Linear':'', 'Gaussian':'Offset, y0','Exponential':'y Offset, y0','TOF':''}            
       eqnlabel = {'Linear':'y = mx+b', 'Gaussian':'y = A exp[-(x-x0)^2/(2 sigma^2)]+y0','Exponential':'y = A exp(-(x-x0)/tau) + y0','TOF':'sqrt(sig0^2 + (vt)^2)'}            

       self.lblFit2Param1.setText(parameterLabels1.get(str(self.comboFit2Type.currentText())))
       self.lblFit2Param2.setText(parameterLabels2.get(str(self.comboFit2Type.currentText())))
       self.lblFit2Param3.setText(parameterLabels3.get(str(self.comboFit2Type.currentText())))
       self.lblFit2Param4.setText(parameterLabels4.get(str(self.comboFit2Type.currentText())))
       self.equationLabel2.setText(eqnlabel.get(str(self.comboFit2Type.currentText())))  
   # *********************************************************       
    

    # *********************************************************   
    def getData(self):
       
        filename = self.DataFileInput.text()
        fp = open(filename,'r')
        lines = fp.read().split('\n')
        fp.close()
       
        filenum = lines[0].split('=')[1]
        while lines[0] != '# ':
            del lines[0]
        del lines[0]
        endofdata = lines.index("# Image Data")
        
        for i in range(len(lines)):
            if lines[i][0] != '#':
                pass
            else:
                del lines[i:endofdata+1]
                break
        allvals = []
        for j in range(len(lines)):
            allvals.append(lines[j].split('=')[1])
        allvals.append(filenum)
        return  allvals
#     ***end getMmt **************************************
     
    # *********************************************************   
    def clear_data(self): 
#        print(self.xvals)
        box = QMessageBox()
        ret = box.question(self,'Clear', 'Are you sure you want to clear the data?', box.Yes | box.No)
        if ret == box.Yes:
            print('Clearing data...')
            self.xvals=[]
            self.yvals=[]
            self.yvals2=[]
            self.ctrlList = []
            self.ctrlVals = []
            self.MmtList = []
            self.MmtVals = []
            self.Mmt2List = []
            self.Mmt2Vals = []
            self.CtrlData = []
            self.mpldata.clear()
#            print(self.xvals)
            self.update_figure()
            self.update_figure2()
            self.meanLE.setText("n/a")
            self.stdLE.setText("n/a")
        else:
            print('Phew...Data is safe!')
       
    # ***** end clear_data *******************************   
     
    # *********************************************************   
    def export_data(self): 
        fname, filter = QFileDialog.getSaveFileName(self, 'Save Graph as CSV',filter='*.csv')
        headerstring = self.comboControl.currentText() +','+self.comboMmt.currentText() + ','+self.comboMmt2.currentText() + '\n'
        # make data into an array
        dataout  = np.vstack((self.xvals,self.yvals,self.yvals2))
        dataout = np.transpose(dataout)
        if fname != '':
            try:
                np.savetxt(fname, dataout, fmt = '%s',header = headerstring, delimiter=',',comments='')
            except:
                print('There was a problem saving the file.')
    # ***** end clear_data *******************************   
 
    # *********************************************************   
    def export_alldata(self): 
        fname, filter = QFileDialog.getSaveFileName(self, 'Save All Data as CSV', filter='*.csv')
        headerstring = ','.join(self.ctrlList) +','+','.join(self.MmtList) + '\n'
        print("Header: ", headerstring)
#         make data into an array
        dataout  = []
        for i in range(len(self.xvals)):
#            print(self.CtrlData[i] + self.MmtData[i])
#            dataout.append(self.CtrlData[i] + self.MmtData[i])
            #print(self.CtrlData[i])
            dataout.append(self.CtrlData[i])
#        dataout = np.transpose(dataout)
        if fname != '':
            try:
                np.savetxt(fname, dataout, fmt = '%s',header = headerstring, delimiter=',',comments='')
            except:
                print('There was a problem saving the file.')
    # ***** end clear_data *******************************   
  
 # *********************************************************   
    def load_Datapath(self): 
        self.Datapath = QFileDialog.getOpenFileName(self,'Select Data File (AllData)')
        self.DataFileInput.setText(self.Datapath)
        
        self.filename = [self.DataFileInput.text()]
        self.fileinfo = QtCore.QFileInfo(self.filename[0])

        self.file_deleted = False
        
        self.fs_watcher = QtCore.QFileSystemWatcher()
        self.fs_watcher.addPaths(self.filename)
        
        self.fs_watcher.fileChanged.connect(self.file_changed)
    # ***** end load_Datapath *******************************   
  
  # *********************************************************   
    def load_Mmtpath(self): 
        self.Mmtpath = QFileDialog.getOpenFileName(self,'Select Measurements File (Image Data)')
        self.MmtFileInput.setText(self.Mmtpath)
    # ***** end load_Mmtpath *******************************   
    
    # *********************************************************   
    def load_settings(self): 
        defaultdir =  'C:/Users/Anind/Google Drive/LindsayData/FigSettings/'
        self.Settingspath = QFileDialog.getOpenFileName(self,'Select Settings File',defaultdir)
        self.SettingsFileLE.setText(self.Settingspath)
    # ***** end load_settings ******************************* 

    # *********************************************************   
    def save_settings(self): 
        xaxisdata = self.ax1f1.get_xlim()        
        yaxisdata = self.ax1f1.get_ylim()
        xlabeldata = self.ax1f1.get_xlabel()
        ylabeldata = self.ax1f1.get_ylabel()
        titledata = self.ax1f1.get_title() 
        
        xaxisdata2 = self.ax1f2.get_xlim()        
        yaxisdata2 = self.ax1f2.get_ylim()
        xlabeldata2 = self.ax1f2.get_xlabel()
        ylabeldata2 = self.ax1f2.get_ylabel()
        titledata2 = self.ax1f2.get_title()  
        
        print('in save settings', xaxisdata)
        filename = self.SettingsFileLE.text()
        print(filename)
        np.savez(filename, xaxisdata=xaxisdata,yaxisdata=yaxisdata,xlabeldata=xlabeldata,ylabeldata=ylabeldata,titledata=titledata, xaxisdata2=xaxisdata2,yaxisdata2=yaxisdata2,xlabeldata2=xlabeldata2,ylabeldata2=ylabeldata2,titledata2=titledata2)
        print('end save settings', titledata)
       
    # ***** end save_settings *******************************      
                           
    # *********************************************************   
    def apply_settings(self): 
        filename = self.SettingsFileLE.text().split('.')[0] + ".npz"
        print(filename)
        data = np.load(filename)
        
        xaxisdata = data['xaxisdata']        
        yaxisdata = data['yaxisdata']
        xlabeldata = data['xlabeldata']
        ylabeldata = data['ylabeldata']
        titledata = data['titledata']
        print(xaxisdata, yaxisdata, xlabeldata)
        
        self.ax1f1.set_xlim(xaxisdata[0],xaxisdata[1])
        self.ax1f1.set_ylim(yaxisdata[0],yaxisdata[1])
        self.ax1f1.set_xlabel(xlabeldata)
        self.ax1f1.set_ylabel(ylabeldata)
        self.ax1f1.set_title(titledata)
        self.canvas.draw()
        
        if 'xaxisdata2' in data.files:
            xaxisdata2 = data['xaxisdata2']
            yaxisdata2 = data['yaxisdata2']
            xlabeldata2 = data['xlabeldata2']
            ylabeldata2 = data['ylabeldata2']
            titledata2 = data['titledata2']   
        else:
            xaxisdata2 = data['xaxisdata']    
            yaxisdata2 = data['yaxisdata']
            xlabeldata2 = data['xlabeldata']
            ylabeldata2 = data['ylabeldata']
            titledata2 = data['titledata']
       
        self.ax1f2.set_xlim(xaxisdata2[0],xaxisdata2[1])
        self.ax1f2.set_ylim(yaxisdata2[0],yaxisdata2[1])
        self.ax1f2.set_xlabel(xlabeldata2)
        self.ax1f2.set_ylabel(ylabeldata2)
        self.ax1f2.set_title(titledata2)
        self.canvas2.draw()
      
        print('figure updated')
    # ***** end load_Mmtpath *******************************       
        
        
# *********************************************************   
# *********************************************************   
if __name__ == '__main__':
    import sys
    from PyQt5 import QtGui
  
    app = QApplication(sys.argv)
    main = Main()
    
   
    main.show()
    sys.exit(app.exec_())