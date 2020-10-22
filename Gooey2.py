#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun August 30th 20:26:06 2020

@author: addy
"""

# any of these scripts with the suffix 2 are the first implementation of adding the black hole and dark matter options to GUI and cjam calculations 


from PyQt5 import QtCore, QtGui
from PyQt5.QtGui import QIcon, QIntValidator
from PyQt5.QtWidgets import QApplication, QWidget, QInputDialog, QLineEdit, QFileDialog, QLabel, QVBoxLayout, \
    QPushButton, QComboBox, QRadioButton, QCheckBox
import sys
import numpy as np
import ClusterDataCleanUp
import pandas as pd
import CJAM_nestsamp2 as CJAMns

class App(QWidget):
    # Constructor for object App
    def __init__(self):
        super().__init__()
        self.title = ''
        self.left = 10
        self.top = 10
        self.width = 1000
        self.height = 700

        # App data, accessible anywhere within this class
        # variable to store the file name of the database the user has selected, default is Jan2020 V2
        self.databaseFileName = "Baumgardt_GC_database.pkl"
        self.objectToDatasetMapping = {}
        self.outputDirectoryPath = ""
        self.MGEfileName = ""
        self.dataIsRadial = True
        self.NofLivePoints = "100"
        self.Interp = True
        self.addBlackHole = False
        self.addDarkMatter = True
        self.NofDMgaussians = "1"
        self.objectName = ""
        self.datasetName = ""
        self.NofInterpPts = "10"

                
        # Window Member UI Elements
        # box where db file path is displayed + drop down menu
        self.databaseFileNameBox = QLineEdit(self.databaseFileName, self)
        self.datasetDropDown = QComboBox(self)
        self.ClusterDropDown = QComboBox(self)
        self.outputDirectoryPathBox = QLineEdit(self.outputDirectoryPath, self)
        self.MGEfileNameBox = QLineEdit(self.MGEfileName, self)
        self.RadialButton = QRadioButton("Radial", self)
        self.PMButton = QRadioButton("Proper Motion", self)
        self.NofLivePointsBox = QLineEdit(self.NofLivePoints, self)
        self.NofInterpPtsBox = QLineEdit(self.NofInterpPts, self)
        self.InterpCheckBox = QCheckBox("Interpolate CJAM Results", self)
        self.addBlackHoleCheckBox = QCheckBox("Add Black Hole",self)
        self.addDarkMatterCheckBox = QCheckBox("Add Dark Matter",self)
        self.RunButton = QPushButton("Generate PyMultiNest Results", self)
        self.NofDMgaussiansBox = QLineEdit(self.NofDMgaussians, self)


        # UI member config
        # set file name box to read only
        self.databaseFileNameBox.setReadOnly(True)
        self.outputDirectoryPathBox.setReadOnly(True)
        self.MGEfileNameBox.setReadOnly(True)
        # calls following method
        self.initUI()
    
    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        self.setMouseTracking(True)
        # create labels
        labels = []
        rightLabels = []

        labels.append(QLabel("Database", self))
        labels.append(QLabel("Object and Dataset", self))
        labels.append(QLabel("Output Directory", self))
        labels.append(QLabel("MGE File", self))
        labels.append(QLabel("Velocity Dispersion Type", self))
        rightLabels.append(QLabel("Number of Live Points for PyMultiNest Sampling", self))
        # set positions of labels & style
        ypos = 20
        for label in labels:
            label.setFont(QtGui.QFont("", 24, QtGui.QFont.Normal))
            label.move(50, ypos)
            ypos = ypos + 120
        ypos = 20
        for label in rightLabels:
            label.setFont(QtGui.QFont("", 20, QtGui.QFont.Normal))
            label.move(500, ypos)
            ypos = ypos + 120

            
        # creates button to open finder window to select database file   
        databaseSelectButton = QPushButton(self)
        databaseSelectButton.setToolTip('Select a file :)')
        databaseSelectButton.setGeometry(5, 70, 40, 40)
        databaseSelectButton.setIcon(QIcon('folder-icon.png'))
        databaseSelectButton.setIconSize(QtCore.QSize(35, 35))
        # calls getDatabaseFileName when databaseSelectButton is clicked
        databaseSelectButton.clicked.connect(self.getDatabaseFileName)
        self.databaseFileNameBox.setGeometry(50, 70, 400, 40)
        
        # creates button to open finder window to select output directory
        outputDirectorySelectButton = QPushButton(self)
        outputDirectorySelectButton.setToolTip('Select a directory for pymultinest results :)')
        outputDirectorySelectButton.setGeometry(5, 300, 40, 40)
        outputDirectorySelectButton.setIcon(QIcon('folder-icon.png'))
        outputDirectorySelectButton.setIconSize(QtCore.QSize(35, 35))
        outputDirectorySelectButton.clicked.connect(self.getOutputDirectoryName)
        # box for output directory text
        self.outputDirectoryPathBox.setGeometry(50,300,400,40)

        # creates button to open finder window to select R file
        MGEfileSelectButton = QPushButton(self)
        MGEfileSelectButton.setToolTip('Select an MGE file :)')
        MGEfileSelectButton.setGeometry(5, 410, 40, 40)
        MGEfileSelectButton.setIcon(QIcon('folder-icon.png'))
        MGEfileSelectButton.setIconSize(QtCore.QSize(35, 35))
        MGEfileSelectButton.clicked.connect(self.getMGEfileName)
        self.MGEfileNameBox.setGeometry(50, 410, 400, 40)

        # velocity type radio button stuff
        
        self.RadialButton.move(50, 540)
        self.RadialButton.setChecked(True)
        self.PMButton.move(150, 540)
        self.RadialButton.toggled.connect(self.onRadioButtonClicked)
        self.PMButton.toggled.connect(self.onRadioButtonClicked)

        #interpolation checkbox stuff

        self.InterpCheckBox.setChecked(True)
        self.InterpCheckBox.move(500, 150)
        self.InterpCheckBox.clicked.connect(self.InterpCheckBoxClicked)
        self.NofInterpPtsBox.setGeometry(500,225,100,40)
        self.NofInterpPtsBoxLabel = QLabel("Number of Points to Interpolate",self)
        self.NofInterpPtsBoxLabel.setFont(QtGui.QFont("", 12, QtGui.QFont.Normal))
        self.NofInterpPtsBoxLabel.move(500,200)
        self.NofInterpPtsBox.textChanged.connect(self.onNofInterpPtsChanged)
        
        #non-luminous matter checkbox stuff

        self.addBlackHoleCheckBox.setChecked(False)
        self.addBlackHoleCheckBox.move(650, 300)
        self.addBlackHoleCheckBox.clicked.connect(self.addBlackHoleCheckBoxClicked)
        self.addDarkMatterCheckBox.setChecked(True)
        self.addDarkMatterCheckBox.move(500, 300)
        self.addDarkMatterCheckBox.clicked.connect(self.addDarkMatterCheckBoxClicked)
        self.NofDMgaussiansBox.setGeometry(500,375,45,20)
        self.NofDMgaussiansBoxLabel = QLabel("Number of Gaussians to Parameterize Dark Matter Halo",self)
        self.NofDMgaussiansBoxLabel.setFont(QtGui.QFont("", 12, QtGui.QFont.Normal))
        self.NofDMgaussiansBoxLabel.move(500,350)
        self.NofDMgaussiansBox.textChanged.connect(self.onNofDMgaussiansChanged)

        # set limit on number of DM gaussian components
        
        DMgaussianValidator = QIntValidator(0, 3, self)
        self.NofDMgaussiansBox.setValidator(DMgaussianValidator)
      
        # drop down menu for datasets + number of live points text box + run button
        self.ClusterDropDown.setGeometry(45, 170, 200, 40)
        self.datasetDropDown.setGeometry(250, 170, 200, 40)
        self.NofLivePointsBox.setGeometry(500,70,400,40)
        
        # creates button to update Updated_Baumgardt_GC_DB.pkl with most current Baumgardt website information
        updateDBbutton = QPushButton(self)
        updateDBbutton.setToolTip('Fetch most recent data from Baumgardt GC database (Default file contains GC DB V2 Jan 2020).')
        updateDBbutton.setGeometry(5, 20, 40, 40)
        updateDBbutton.setIcon(QIcon('download-butt.png'))
        updateDBbutton.setIconSize(QtCore.QSize(35, 35))
        # connect button to method updateButtonClicked
        updateDBbutton.clicked.connect(self.updateButtonClicked)
        self.ClusterDropDown.currentTextChanged.connect(self.onClusterDropDownChanged)
        self.datasetDropDown.currentTextChanged.connect(self.onDatasetDropDownChanged)
        self.RunButton.move(625,450)
        self.RunButton.clicked.connect(self.RunPyMultiNest)
        self.NofLivePointsBox.textChanged.connect(self.onNofLivePointsChanged)

        validator = QIntValidator(0, 99999, self)
        self.NofLivePointsBox.setValidator(validator)

        self.show()
        self.populateClusterOptions()
        self.updateUI()
    # define method for opening/using finder + selecting a file
    def openFileNameDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;Python Files (*.py)", options=options)
        if fileName:
            return fileName

    def onRadioButtonClicked(self):
        toggledRadioButton = self.sender()
        if toggledRadioButton.isChecked():
            self.dataIsRadial = toggledRadioButton.text() == "Radial"

    def InterpCheckBoxClicked(self):
         self.Interp = self.InterpCheckBox.isChecked()
         
    def addBlackHoleCheckBoxClicked(self):
         self.addBlackHole = self.addBlackHoleCheckBox.isChecked()
         
    def addDarkMatterCheckBoxClicked(self):
         self.addDarkMatter = self.addDarkMatterCheckBox.isChecked()
        
    # method for selecting a folder in finder
    def openFolderNameDialog(self):
        fileName = QFileDialog.getExistingDirectory(self,"Select Directory", "", QFileDialog.ShowDirsOnly)
        if fileName:
            return fileName
        
    # define method to get name of database selected w/ finder
    def getDatabaseFileName(self):
        fileName = self.openFileNameDialog()
        if fileName:
            self.databaseFileName = fileName
            self.updateUI()
            self.populateClusterOptions()
            
    def getOutputDirectoryName(self):
        fileName = self.openFolderNameDialog()
        if fileName:
            self.outputDirectoryPath = fileName
            self.updateUI()
            
    def getMGEfileName(self):
        fileName = self.openFileNameDialog()
        if fileName:
            self.MGEfileName = fileName
            self.updateUI()
            
    # define method to search through selected database and append each instance of a unique dataset to dataset drop down menu
    def populateDatasetOptions(self):
        self.datasetDropDown.clear()
        currentCluster = self.ClusterDropDown.currentText()
        if currentCluster != "":
            for dataset in self.objectToDatasetMapping[currentCluster]:
                self.datasetDropDown.addItem(dataset)
        self.datasetDropDown.model().sort(0)
        self.datasetDropDown.setCurrentIndex(0)
            
    def populateClusterOptions(self):
        self.ClusterDropDown.clear()
        ClusterSet = set()
        data = pd.read_pickle(self.databaseFileName)
        ClusterColumn = data["Cluster"]
        TypeColumn = data["Type"]
        for i in range(len(ClusterColumn)):
            cluster = ClusterColumn[i]
            ClusterSet.add(cluster)
            if cluster not in self.objectToDatasetMapping.keys():
                self.objectToDatasetMapping[cluster] = set()
            self.objectToDatasetMapping[cluster].add(TypeColumn[i])
        # add to drop down menu
        for cluster in ClusterSet:
            self.ClusterDropDown.addItem(cluster)
        
        self.populateDatasetOptions()
        self.ClusterDropDown.model().sort(0)
        self.ClusterDropDown.setCurrentIndex(0)
        
    def onClusterDropDownChanged(self):
        self.populateDatasetOptions()
        self.objectName = self.ClusterDropDown.currentText()

    def onDatasetDropDownChanged(self):
        self.datasetName = self.datasetDropDown.currentText()

    def onNofLivePointsChanged(self):
        self.NofLivePoints = self.NofLivePointsBox.text()
        
    def onNofInterpPtsChanged(self):
        self.NofInterpPts = self.NofInterpPtsBox.text()
        
    def onNofDMgaussiansChanged(self):
        self.NofDMgaussians = self.NofDMgaussiansBox.text()
        
    def updateUI(self):
        if self.databaseFileName != "":
            self.databaseFileNameBox.setText(self.databaseFileName)
            self.outputDirectoryPathBox.setText(self.outputDirectoryPath)
            self.MGEfileNameBox.setText(self.MGEfileName)
    
    def updateButtonClicked(self):
        ClusterDataCleanUp.UpdateBaumgardtDB()
        self.databaseFileName = "Updated_Baumgardt_GC_DB.pkl"
        self.updateUI()

    def RunPyMultiNest(self):
        CJAMns.GeneratePyMultiNestResults(self.databaseFileName, self.outputDirectoryPath, self.objectName, self.datasetName, int(self.NofLivePoints), self.Interp, self.dataIsRadial, self.MGEfileName, int(self.NofInterpPts),self.addBlackHole,self.addDarkMatter,int(self.NofDMgaussians))





if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())
