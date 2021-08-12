"""
 This script is calls the main sphere-overburden routine from sphereresponse.py with the user determined values that are input
 in the widgets that are defined and created in options_menu.py

 This python file is compiled into an exe by calling Pyinstaller in a directory with sphereresponse.py and options_menu.py
 the size of the executable can be decreased by omitting unecessary python libraries when compiling
"""

# Python modules
import sys
import qdarkstyle
import matplotlib.pyplot as plt
import PyQt4.QtCore as QtCore
# 3rd party modules
import matplotlib
import matplotlib.backends.backend_qt4agg as backend_qt4agg
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import PyQt4.QtGui as QtGui
import os, sys
import csv
from os.path import expanduser
import threading

# Local application modules
from sphereresponse import sphereresponse
from options_menu import OptionsMenu
import resources
import numpy as np

APP_NAME = 'EM sphere-overburden response'
AUTHOR = 'Anthony Zamperoni'


class AppForm(QtGui.QMainWindow):
    def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self, parent)

        # Set the window title
        self.setWindowTitle(APP_NAME)
        self.imported = False

        # Create the options menu in a dock widget
        self.options_menu = OptionsMenu()
        dock = QtGui.QDockWidget('Options', self)
        dock.setFeatures(
            QtGui.QDockWidget.NoDockWidgetFeatures |
            QtGui.QDockWidget.DockWidgetMovable |
            QtGui.QDockWidget.DockWidgetFloatable
        )
        dock.setAllowedAreas(
            QtCore.Qt.LeftDockWidgetArea | QtCore.Qt.RightDockWidgetArea,
        )
        dock.setWidget(self.options_menu)
        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, dock)

        # Connect the signals from the options menu
        self.connect(self.options_menu.update_btn, QtCore.SIGNAL(
            'clicked()'),
            self.launch_selenium_Thread,
        )
        self.connect(self.options_menu.clear_graph_btn, QtCore.SIGNAL(
            'clicked()'),
            self.clear_graph,
        )
        self.connect(self.options_menu.legend_cb, QtCore.SIGNAL(
            'stateChanged(int)'),
            self.redraw_graph,
        )
        self.connect(self.options_menu.grid_cb, QtCore.SIGNAL(
            'stateChanged(int)'),
            self.redraw_graph,
        )
        self.connect(self.options_menu.read_data_btn, QtCore.SIGNAL(
            'clicked()'), self.read_csv
        )

        if self.options_menu.grid_cb.isChecked() == True and self.options_menu.legend_cb.isChecked() == True:

            self.fig = Figure()
            self.canvas = backend_qt4agg.FigureCanvasQTAgg(self.fig)
            self.canvas.setParent(self)    
            self.ax1 = self.axes = self.fig.add_subplot(111)
            self.ax2 = self.axes = self.fig.add_subplot(211)

        if self.options_menu.grid_cb.isChecked() == False or self.options_menu.legend_cb.isChecked() == False:
            self.fig = Figure()
            self.canvas = backend_qt4agg.FigureCanvasQTAgg(self.fig)
            self.canvas.setParent(self)


        self.status_text = QtGui.QLabel("Set paramters and select response components to be plotted")
        self.statusBar().addWidget(self.status_text, 0)
        self.statusBar().setFont(QtGui.QFont("Times", 20, QtGui.QFont.Bold))
        self.progressBar = QtGui.QProgressBar(self)
        self.statusBar().addPermanentWidget(self.progressBar, 1)
        self.statusBar().addWidget(self.progressBar)




        # Initialize the graph
        self.clear_graph()

        # Set the graph as the main window widget
        self.setCentralWidget(self.canvas)

        # Create the exit application function in the menubar

        file_exit_action = QtGui.QAction('E&xit', self)
        file_exit_action.setToolTip('Exit')
        file_exit_action.setIcon(QtGui.QIcon(':/resources/door_open.png'))
        self.connect(
            file_exit_action,
            QtCore.SIGNAL('triggered()'),
            self.close,
        )

        about_action = QtGui.QAction('&About', self)
        about_action.setToolTip('About')
        about_action.setIcon(QtGui.QIcon(':/resources/icon_info.gif'))
        self.connect(
            about_action,
            QtCore.SIGNAL('triggered()'),
            self.show_about,
        )

        # Create the menubar add further functionality at later date

        file_menu = self.menuBar().addMenu('&File')
        #file_menu.addAction(file_preview_waveform)
        file_menu.addAction(file_exit_action)
        file_menu = self.menuBar().addMenu('&Edit')
        file_menu = self.menuBar().addMenu('&View')
        help_menu = self.menuBar().addMenu('&Help')
        help_menu.addAction(about_action)

    def read_csv(self):
        """
                Reads in waveform from csv to convolve
        """
        self.filePath = QtGui.QFileDialog.getOpenFileName(self, 'Select CSV File:', '', "CSV data files (*.csv)")
        if self.filePath == "":
            self.status_text.setText("No waveform data selected")
            return
        else:
            with open(self.filePath) as input_file:
                self.imported = True
                self.waveformdata = np.genfromtxt(input_file,delimiter = ',')
                self.status_text.setText("Successfully loaded and updated waveform parameters")
           
    
                if self.waveformdata.shape[1] >= 2:
                    self.windows = self.waveformdata.T[0]
                    self.wave = self.waveformdata.T[1]
                else:
                    self.windows = self.waveformdata.T[0]


    def calculate_data(self):

        """
            A function call to the main function that calculates the response using the user inputted values being
            called from the widgets defined in options_menu.py
        """

        sphere = sphereresponse()
        sphere.a = self.options_menu.a_sb.value()
        sphere.rsp = np.array([int(n) for n in self.options_menu.rspop.text().split(',')],dtype=np.int64)
        sphere.offset_tx_rx = np.array([int(n) for n in self.options_menu.txrx.text().split(',')], dtype = np.int64)
        sphere.rtx = np.array([int(n) for n in self.options_menu.tx.text().split(',')], dtype = np.int64)
        sphere.thick_ob = self.options_menu.thick_ob_sb.value()
        sphere.sigma_sp = self.options_menu.sigma_sp_sb.value()
        sphere.sigma_ob = self.options_menu.sigma_ob_sb.value()
        sphere.strike = self.options_menu.strike.value()
        sphere.dip = self.options_menu.dip.value()
        sphere.P = self.options_menu.pulse.value()
        sphere.T = self.options_menu.timedelta_sb.value()
        if self.imported == True:
            sphere.wc = self.windows
            sphere.wave = self.wave

        # Checks if the sphere is dipping or not passed as value to main routine

        if self.options_menu.dip.value() == 0:
            sphere.apply_dip = 0
        else: sphere.apply_dip = 1

        if sphere.sigma_sp == 0:
            sphere.sigma_sp = 0.00000000000001
        if sphere.sigma_ob == 0:
            sphere.sigma_ob = 0.00000000000001



        results = sphere.calculate()


        """
            The following is an if-then loop for plotting the different components of the response given which boxes are checked
            This will be rewritten more efficiently 
        """


        if (self.options_menu.grid_cb.isChecked() and self.options_menu.legend_cb.isChecked()) and self.options_menu.gridy_cb.isChecked():
            

            self.axes = self.fig.add_subplot(3,1,1)
            #self.subplot(2, 1, 1)
            x=sphere.H_tot_x
            i=0
            while i < len(sphere.H_tot_x):
                self.axes.plot(np.linspace(sphere.profile[0][0], sphere.profile[0][100], 101), sphere.H_tot_x[i],color= '0.4')  # will have to change x axis for changing param
                i += 1

            self.axes.set_xlabel('Profile (m)')
            self.axes.set_ylabel('x-component (A/m)')
            self.axes.grid(True, which = 'both', ls = '-')
       # the first subplot in the first figure

            self.axes = self.fig.add_subplot(3,1,2)
            z=sphere.H_tot_z
            k=0
            while k < len(sphere.H_tot_z):
                self.axes.plot(np.linspace(sphere.profile[0][0], sphere.profile[0][100], 101), sphere.H_tot_z[k],color= '0.4')  # will have to change x axis for changing param
                k += 1

            self.axes.set_xlabel('Profile (m)')
            self.axes.set_ylabel('z-component (A/m)')
            self.axes.grid(True, which = 'both', ls = '-')
            
            self.axes = self.fig.add_subplot(3,1,3)
            #self.subplot(2, 1, 1)
            y=sphere.H_tot_y
            i=0
            while i < len(sphere.H_tot_y):
                self.axes.plot(np.linspace(sphere.profile[0][0], sphere.profile[0][100], 101), sphere.H_tot_y[i],color= '0.4')  # will have to change x axis for changing param
                i += 1

            self.axes.set_xlabel('Profile (m)')
            self.axes.set_ylabel('y-component (A/m)')
            self.axes.grid(True, which = 'both', ls = '-')            
            #self.canvas.addWidget(self.navi_toolbar)
            
            self.canvas.draw()
            
            
            
        elif self.options_menu.grid_cb.isChecked() and self.options_menu.legend_cb.isChecked() and self.options_menu.gridy_cb.isChecked() == False:
                
    
            self.axes = self.fig.add_subplot(2,1,1)
                #self.subplot(2, 1, 1)
            x=sphere.H_tot_x
            i=0
            while i < len(sphere.H_tot_x):
                self.axes.plot(np.linspace(sphere.profile[0][0], sphere.profile[0][100], 101), sphere.H_tot_x[i],color= '0.4')  # will have to change x axis for changing param
                i += 1
    
            self.axes.set_xlabel('Profile (m)')
            self.axes.set_ylabel('x-component (A/m)')
            self.axes.grid(True, which = 'both', ls = '-')
           # the first subplot in the first figure
    
            self.axes = self.fig.add_subplot(2,1,2)
            z=sphere.H_tot_z
            k=0
            while k < len(sphere.H_tot_z):
                self.axes.plot(np.linspace(sphere.profile[0][0], sphere.profile[0][100], 101), sphere.H_tot_z[k],color= '0.4')  # will have to change x axis for changing param
                k += 1
    
            self.axes.set_xlabel('Profile (m)')
            self.axes.set_ylabel('z-component (A/m)')
            self.axes.grid(True, which = 'both', ls = '-')            
                #self.canvas.addWidget(self.navi_toolbar)
                
            self.canvas.draw()
            
        elif self.options_menu.grid_cb.isChecked() == False and self.options_menu.legend_cb.isChecked() and self.options_menu.gridy_cb.isChecked():
    
    
            self.axes = self.fig.add_subplot(2,1,1)
                #self.subplot(2, 1, 1)
            x=sphere.H_tot_x
            i=0
            while i < len(sphere.H_tot_x):
                self.axes.plot(np.linspace(sphere.profile[0][0], sphere.profile[0][100], 101), sphere.H_tot_x[i],color= '0.4')  # will have to change x axis for changing param
                i += 1
    
            self.axes.set_xlabel('Profile (m)')
            self.axes.set_ylabel('x-component (A/m)')
            self.axes.grid(True, which = 'both', ls = '-')
        # the first subplot in the first figure
    
            self.axes = self.fig.add_subplot(2,1,2)
            y=sphere.H_tot_y
            k=0
            while k < len(sphere.H_tot_z):
                self.axes.plot(np.linspace(sphere.profile[0][0], sphere.profile[0][100], 101), sphere.H_tot_y[k],color= '0.4')  # will have to change x axis for changing param
                k += 1
    
            self.axes.set_xlabel('Profile (m)')
            self.axes.set_ylabel('y-component (A/m)')
            self.axes.grid(True, which = 'both', ls = '-')            
                #self.canvas.addWidget(self.navi_toolbar)
    
            self.canvas.draw()
        elif self.options_menu.grid_cb.isChecked() and self.options_menu.legend_cb.isChecked() == False and self.options_menu.gridy_cb.isChecked():
    
    
            self.axes = self.fig.add_subplot(2,1,1)
                #self.subplot(2, 1, 1)
            z=sphere.H_tot_z
            i=0
            while i < len(sphere.H_tot_x):
                self.axes.plot(np.linspace(sphere.profile[0][0], sphere.profile[0][100], 101), sphere.H_tot_z[i],color= '0.4')  # will have to change x axis for changing param
                i += 1
    
            self.axes.set_xlabel('Profile (m)')
            self.axes.set_ylabel('z-component (A/m)')
            self.axes.grid(True, which = 'both', ls = '-')
        # the first subplot in the first figure
    
            self.axes = self.fig.add_subplot(2,1,2)
            y=sphere.H_tot_y
            k=0
            while k < len(sphere.H_tot_z):
                self.axes.plot(np.linspace(sphere.profile[0][0], sphere.profile[0][100], 101), sphere.H_tot_y[k],color= '0.4')  # will have to change x axis for changing param
                k += 1
    
            self.axes.set_xlabel('Profile (m)')
            self.axes.set_ylabel('y-component (A/m)')
            self.axes.grid(True, which = 'both', ls = '-')            
                #self.canvas.addWidget(self.navi_toolbar)
    
            self.canvas.draw()    


        elif self.options_menu.legend_cb.isChecked() and self.options_menu.grid_cb.isChecked() == False and self.options_menu.gridy_cb.isChecked() == False:

            self.fig.clf()
            self.axes = self.fig.add_subplot(111)
            x=sphere.H_tot_x
            i=0
            while i < len(sphere.H_tot_x):
                self.axes.plot(np.linspace(sphere.profile[0][0], sphere.profile[0][100], 101), sphere.H_tot_x[i],color= '0.4')  # will have to change x axis for changing param
                i += 1


            self.axes.set_xlabel('Profile (m)')
            self.axes.set_ylabel('x-component (A/m)')
            self.axes.grid(True, which = 'both', ls = '-')
            #self.axes.plot(x)
            self.canvas.draw()


        elif self.options_menu.grid_cb.isChecked() and self.options_menu.legend_cb.isChecked() == False and self.options_menu.gridy_cb.isChecked() == False:

            self.fig.clf()
            self.axes = self.fig.add_subplot(111)
            z=sphere.H_tot_z
            k=0
            while k < len(sphere.H_tot_z):
                self.axes.plot(np.linspace(sphere.profile[0][0], sphere.profile[0][100], 101), sphere.H_tot_z[k],color= '0.4')  # will have to change x axis for changing param
                k += 1

            self.axes.set_xlabel('Profile (m)')
            self.axes.set_ylabel('z-component (A/m)')
            self.axes.grid(True, which = 'both', ls = '-')
            #self.axes.plot(x)
            self.canvas.draw()
        elif self.options_menu.grid_cb.isChecked() == False and self.options_menu.legend_cb.isChecked() == False and self.options_menu.gridy_cb.isChecked():

            self.fig.clf()
            self.axes = self.fig.add_subplot(111)
            y=sphere.H_tot_y
            k=0
            while k < len(sphere.H_tot_z):
                self.axes.plot(np.linspace(sphere.profile[0][0], sphere.profile[0][100], 101), sphere.H_tot_y[k],color= '0.4')  # will have to change x axis for changing param
                k += 1

            self.axes.set_xlabel('Profile (m)')
            self.axes.set_ylabel('y-component (A/m)')
            self.axes.grid(True, which = 'both', ls = '-')
            #self.axes.plot(x)
            self.canvas.draw()
            
            
            
            
        self.progressBar.setRange(0,1)
        self.status_text.setText("Finished")
        self.statusBar().setFont(QtGui.QFont("Times", 20, QtGui.QFont.Bold))




    def clear_graph(self):

        self.redraw_graph()

    def redraw_graph(self):

        self.fig.clf()

        self.canvas.draw()

    def launch_selenium_Thread(self):

        """
            A function to prevent the program from becoming unresponsive while the response is being calculated/plotted
        """
        t = threading.Thread(target=self.calculate_data)
        self.status_text.setText("Generating response")

        # Create updating progress bar
        self.statusBar().setFont(QtGui.QFont("Times", 20, QtGui.QFont.Bold))
        self.progressBar.setRange(0,0)
        t.start()    

   


    def show_about(self):
        """
        Display the "about" dialog box.
        """
        message = '''<font size="+2">%s</font>
            <p>A sphere - overburden response plotter written in Python.
            <p>Written by %s,
            <a href="http://opensource.org/licenses/MIT">MIT Licensed</a>
            <p>Icons from <a href="http://www.famfamfam.com/">famfamfam</a> and
            <a href="http://commons.wikimedia.org/">Wikimedia
            Commons</a>.''' % (APP_NAME, AUTHOR)

        QtGui.QMessageBox.about(self, 'About ' + APP_NAME, message)

if __name__ == "__main__":
    # Dark UI theme
    dark_stylesheet = qdarkstyle.load_stylesheet_pyqt()
    app = QtGui.QApplication(sys.argv)
    app.setWindowIcon(QtGui.QIcon(':/resources/icon.svg'))
    app.setStyleSheet(dark_stylesheet)
    form = AppForm()
    form.show()
    app.exec_()
