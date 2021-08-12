# options_menu
# PyQt4 framework used to generate and modify the GUI used in the sphere-overburden program
# Generates and sets initial values for the input box widgets on the left side of the user interface
# Creates checkbox and button widgets to allow for the plotting of specific components of the sphere-overburden response



# 3rd party modules import
import PyQt4.QtCore as QtCore
import PyQt4.QtGui as QtGui



class OptionsMenu(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)

        # Create the spinbox widgets used for inputting sphere-overburden parameters (thickness, conductivity)

        self.thick_ob_sb = QtGui.QDoubleSpinBox()
        self.sigma_sp_sb = QtGui.QDoubleSpinBox()
        self.sigma_ob_sb = QtGui.QDoubleSpinBox()

        # Create text input widget for transmitter-receiver geometry
        # Set position / offset values using .setText

        self.tx = QtGui.QLineEdit()
        self.tx.setText('0,0,1')



        self.rx = QtGui.QLineEdit()
        self.tx.setText('0,0,1')

        self.txrx = QtGui.QLineEdit()
        self.txrx.setText('125,0,56')

        # Setting labels and positions of variables in widgets
        # Defining limits for integers values

        for widget in (self.thick_ob_sb, self.sigma_sp_sb, self.sigma_ob_sb):
            widget.setRange(0, 1000)
            widget.setSingleStep(1)

        coeff_grid = QtGui.QGridLayout()
        coeff_grid.addWidget(QtGui.QLabel('Receiver Position(m)'), 0, 0)
        coeff_grid.addWidget(self.tx, 0, 1)
        coeff_grid.addWidget(QtGui.QLabel('Transmitter Position(m)'), 1, 0)
        coeff_grid.addWidget(self.txrx, 1, 1)
        coeff_grid.addWidget(QtGui.QLabel('Overburden Thickness (m)'), 2, 0)
        coeff_grid.addWidget(self.thick_ob_sb, 2, 1)
        coeff_grid.addWidget(QtGui.QLabel('Sphere Conductivity (S/m)'), 3, 0)
        coeff_grid.addWidget(self.sigma_sp_sb, 3, 1)
        coeff_grid.addWidget(QtGui.QLabel('Overburden Conductivity (S/m)'), 4, 0)
        coeff_grid.addWidget(self.sigma_ob_sb, 4, 1)
        coeff_gb = QtGui.QGroupBox('')
        coeff_gb.setLayout(coeff_grid)

        # Create the second parameters widget (dip, strike, dipole moment, radius etc)

        self.strike = QtGui.QDoubleSpinBox()
        self.strike.setRange(0, 360)
        self.strike.setSingleStep(1)

        self.dip = QtGui.QDoubleSpinBox()
        self.dip.setRange(0, 360)
        self.dip.setSingleStep(1)

        self.dipole = QtGui.QDoubleSpinBox()
        self.dipole.setRange(0, 1000000000)
        self.dipole.setDecimals(2)


        self.pulse = QtGui.QDoubleSpinBox()
        self.pulse.setRange(0, 100000000)
        self.pulse.setDecimals(8)
        

        self.rspop = QtGui.QLineEdit()
        self.rspop.setText('0,0,-200')


        self.a_sb = QtGui.QDoubleSpinBox()
        self.a_sb.setRange(0,1000)
        self.a_sb.setSingleStep(1)

        self.timedelta_sb = QtGui.QDoubleSpinBox()
        self.timedelta_sb.setRange(0, 10000)
        self.timedelta_sb.setDecimals(4)

        # Setting labels and positions of variables in widgets
        # Defining limits for integers values

        other_grid = QtGui.QGridLayout()
        other_grid.addWidget(QtGui.QLabel('Strike'), 2, 0)
        other_grid.addWidget(self.strike, 2, 1)
        other_grid.addWidget(QtGui.QLabel('Dip'), 3, 0)
        other_grid.addWidget(self.dip, 3, 1)
        other_grid.addWidget(QtGui.QLabel('Pulse Length'), 4, 0)
        other_grid.addWidget(self.pulse, 4, 1)
        other_grid.addWidget(QtGui.QLabel('Period'), 5, 0)
        other_grid.addWidget(self.timedelta_sb, 5, 1)
        other_grid.addWidget(QtGui.QLabel('Dipole Moment'), 6, 0)
        other_grid.addWidget(self.dipole, 6, 1)


        other_grid.addWidget(QtGui.QLabel('Sphere position (m)'), 0, 0)
        other_grid.addWidget(self.rspop, 0, 1)
        other_grid.addWidget(QtGui.QLabel('Sphere radius (m)'), 1, 0)
        other_grid.addWidget(self.a_sb, 1, 1)
        other_gb = QtGui.QGroupBox('')
        other_gb.setLayout(other_grid)



        # Create the "Graph Options" widgets
        # Create checkbox for user to choose which components of response to plot

        self.legend_cb = QtGui.QCheckBox('x-component')
        self.legend_cb.setChecked(False)
        self.connect(self.legend_cb, QtCore.SIGNAL(
            'stateChanged(int)'),
            self.legend_change,
        )
        self.grid_cb = QtGui.QCheckBox('z-component')
        self.grid_cb.setChecked(False)

        self.gridy_cb = QtGui.QCheckBox('y-component')
        self.gridy_cb.setChecked(False)

        # Self.legend_loc_lbl = QtGui.QLabel('waveform data')
        # Self.legend_loc_cb = QtGui.QPushButton('read in from csv',self)

        cb_box = QtGui.QHBoxLayout()

        # Create plot area

        cb_box.addWidget(self.legend_cb)
        cb_box.addWidget(self.grid_cb)
        cb_box.addWidget(self.gridy_cb)

        legend_box = QtGui.QHBoxLayout()
        #legend_box.addWidget(self.legend_loc_cb)
        legend_box.addStretch()

        graph_box = QtGui.QVBoxLayout()
        graph_box.addLayout(cb_box)
        #graph_box.addWidget(self.legend_loc_lbl)
        graph_box.addLayout(legend_box)

        graph_gb = QtGui.QGroupBox('')
        graph_gb.setLayout(graph_box)

        # Create the update/reset plot buttons

        self.update_btn = QtGui.QPushButton(
            QtGui.QIcon(':/resources/calculator.png'),
            'Plot Response',
        )
        self.reset_values_btn = QtGui.QPushButton(
            QtGui.QIcon(':/resources/arrow_undo.png'),
            'Reset Values',
        )
        self.clear_graph_btn = QtGui.QPushButton(
            QtGui.QIcon(':/resources/chart_line_delete.png'),
            'Clear Plot',
        )
        self.connect(self.reset_values_btn, QtCore.SIGNAL(
            'clicked()'),
            self.reset_values,
        )

        self.read_data_btn = QtGui.QPushButton(
            QtGui.QIcon(':/resources/chart_line_delete.png'),'Import Waveform Data',
        )

        # Create the main layout with widgets and plotting area
        container = QtGui.QVBoxLayout()
        container.addWidget(coeff_gb)
        container.addWidget(other_gb)
        container.addWidget(graph_gb)
        container.addWidget(self.update_btn)
        container.addStretch()
        container.addWidget(self.reset_values_btn)
        container.addWidget(self.clear_graph_btn)
        container.addWidget(self.read_data_btn)
        self.setLayout(container)

        # Populate and reset the widgets with values
        self.reset_values()

    def reset_values(self):
        """
        Sets the default values of the option widgets.
        """
        self.a_sb.setValue(100)
        self.rspop.setText('0,0,-200')
        self.thick_ob_sb.setValue(4)
        self.sigma_sp_sb.setValue(0.5)
        self.sigma_ob_sb.setValue(0.03)
        self.strike.setValue(90)
        self.dip.setValue(0)
        self.dipole.setValue(1847300)
        self.pulse.setValue(0.00398)
        self.timedelta_sb.setValue(0.03)
        self.txrx.setText('125,0,56')
        self.tx.setText('0,0,120')

    def legend_change(self):
        return
        #self.legend_loc_cb.setEnabled(self.legend_cb.isChecked())
        #self.legend_loc_lbl.setEnabled(self.legend_cb.isChecked())
