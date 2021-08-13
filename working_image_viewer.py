"""
teophile lemay, summer 2021

This file contains code for the roi-analysis/image viewer gui
"""


from PyQt5 import QtCore, QtGui, QtWidgets
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import fluorescence_functions as func
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from qtrangeslider import QRangeSlider
import json


class Canvas(FigureCanvas):
    def __init__(self, parent):
        self.fig = plt.figure()
        self.ax = plt.axes([0., 0., 1., 1.])
        self.fig.add_axes(self.ax)
        super().__init__(self.fig)
        self.setParent(parent)

    def no_ticks(self):
        #this function makes sure no ticks are shown
        self.ax.set_axis_off()

    def wipe(self):
        #this function clears what was proviously drawn on the ax.
        self.ax.clear()
    
    def define_axis(self, left=0.15, bot=0.12, width=0.8, height=0.75):
        #this function lets the user define a custom position for the ax inside the figure (to show ticks, labels, etc.)
        self.ax = plt.axes([left, bot, width, height])


class Display_Plots(QtWidgets.QDialog):
    #this class makes a popup window to show decay plots, spectral intensity plots and all that
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi()

    def setupUi(self):
        self.resize(1260, 790)
        self.donorPlotLabel = QtWidgets.QLabel(self)
        self.donorPlotLabel.setGeometry(10, 10, 200, 20)
        self.donorDecayPlot = Canvas(self)
        self.donorDecayPlot.define_axis()
        self.donorDecayPlot.setGeometry(10, 35, 600, 500)
        self.donorDecayToolbar = NavigationToolbar(self.donorDecayPlot, self)
        self.donorDecayToolbar.setGeometry(10, 540, 600, 50)
        self.secondaryPlotLabel = QtWidgets.QLabel(self)
        self.secondaryPlotLabel.setGeometry(630, 10, 200, 20)
        self.secondaryPlot = Canvas(self)
        self.secondaryPlot.define_axis()
        self.secondaryPlot.setGeometry(630, 35, 600, 500)
        self.secondaryPlotToolbar = NavigationToolbar(self.secondaryPlot, self)
        self.secondaryPlotToolbar.setGeometry(630, 540, 600, 50)
        self.donorT0Label = QtWidgets.QLabel(self)
        self.donorT0Label.setGeometry(30, 605, 250, 20)
        self.donorT0Label.setText('Current donor FLIM t0 bounds:')
        self.donorT0Slider = QRangeSlider(self)
        self.donorT0Slider.setOrientation(QtCore.Qt.Horizontal)
        self.donorT0Slider.setGeometry(30, 630, 550, 30)
        self.donorLowerLabel = QtWidgets.QLabel(self)
        self.donorLowerLabel.setGeometry(40, 665, 180, 20)
        self.donorUpperLabel = QtWidgets.QLabel(self)
        self.donorUpperLabel.setGeometry(400, 665, 180, 20)
        self.donorWidthSpinBox = QtWidgets.QSpinBox(self)
        self.donorWidthSpinBox.setGeometry(60, 725, 45, 25)
        self.donorWidthSpinBox.setMinimum(0)
        self.donorWidthSpinBox.setMaximum(50)
        self.donorSpinBoxLabel = QtWidgets.QLabel(self)
        self.donorSpinBoxLabel.setGeometry(40, 695, 150, 25)
        self.donorSpinBoxLabel.setText('Width from peak:')
        self.donorSetButton = QtWidgets.QPushButton(self)
        self.donorSetButton.setGeometry(115, 725, 50, 25)
        self.donorSetButton.setText('Set')
        self.showDonorCheckBox = QtWidgets.QCheckBox(self)
        self.showDonorCheckBox.setGeometry(400, 695, 150, 20)
        self.showDonorCheckBox.setLayoutDirection(QtCore.Qt.LayoutDirection(QtCore.Qt.RightToLeft))
        self.showDonorCheckBox.setText('Show T0 Bounds')
        self.acceptorT0Label = QtWidgets.QLabel(self)
        self.acceptorT0Label.setGeometry(650, 605, 250, 20)
        self.acceptorT0Label.setText('Current acceptor FLIM t0 bounds:')
        self.acceptorT0Slider = QRangeSlider(self)
        self.acceptorT0Slider.setOrientation(QtCore.Qt.Horizontal)
        self.acceptorT0Slider.setGeometry(650, 630, 550, 30)
        self.acceptorLowerLabel = QtWidgets.QLabel(self)
        self.acceptorLowerLabel.setGeometry(660, 665, 180, 20)
        self.acceptorUpperLabel = QtWidgets.QLabel(self)
        self.acceptorUpperLabel.setGeometry(1020, 665, 190, 20)
        self.acceptorWidthSpinBox = QtWidgets.QSpinBox(self)
        self.acceptorWidthSpinBox.setGeometry(680, 725, 45, 25)
        self.acceptorWidthSpinBox.setMinimum(0)
        self.acceptorWidthSpinBox.setMaximum(50)
        self.acceptorSpinBoxLabel = QtWidgets.QLabel(self)
        self.acceptorSpinBoxLabel.setGeometry(660, 695, 150, 25)
        self.acceptorSpinBoxLabel.setText('Width from peak:')
        self.acceptorSetButton = QtWidgets.QPushButton(self)
        self.acceptorSetButton.setGeometry(735, 725, 50, 25)
        self.acceptorSetButton.setText('Set')
        self.showAcceptorCheckBox = QtWidgets.QCheckBox(self)
        self.showAcceptorCheckBox.setGeometry(1040, 695, 150, 20)
        self.showAcceptorCheckBox.setLayoutDirection(QtCore.Qt.LayoutDirection(QtCore.Qt.RightToLeft))
        self.showAcceptorCheckBox.setText('Show T0 Bounds')
        self.okButton = QtWidgets.QPushButton(self)
        self.okButton.setGeometry(1000, 735, 125, 35)
        self.okButton.setText('Save T0 Bounds')
        self.cancelButton = QtWidgets.QPushButton(self)
        self.cancelButton.setGeometry(1135, 735, 115, 35)
        self.cancelButton.setText('Cancel')
        #signal connections
        self.donorT0Slider.valueChanged.connect(self.donor_slider_changed)
        self.acceptorT0Slider.valueChanged.connect(self.acceptor_slider_changed)
        self.showDonorCheckBox.clicked.connect(self.show_donor_bounds)
        self.showAcceptorCheckBox.clicked.connect(self.show_acceptor_bounds)
        self.donorSetButton.clicked.connect(self.donor_width)
        self.acceptorSetButton.clicked.connect(self.acceptor_width)
        self.okButton.clicked.connect(self.accept)
        self.cancelButton.clicked.connect(self.reject)
        #disable acceptor T0 functions in case of INO file
        self.acceptorT0Label.setEnabled(False)
        self.acceptorT0Slider.setEnabled(False)
        self.acceptorLowerLabel.setEnabled(False)
        self.acceptorUpperLabel.setEnabled(False)
        self.showAcceptorCheckBox.setEnabled(False)
        self.acceptorWidthSpinBox.setEnabled(False)   
        self.acceptorSpinBoxLabel.setEnabled(False)  
        self.acceptorSetButton.setEnabled(False)

    def donor_width(self):
        #this function sets the bounds to around the peak by a user defined amount
        width = self.donorWidthSpinBox.value()
        if width > self.donortop:
            width = self.donortop
            self.donorWidthSpinBox.setValue(width)
        self.donorT0Lims = (self.donortop - width, self.donortop + width)
        self.donorT0Slider.setValue(self.donorT0Lims)

    def acceptor_width(self):
        #this function sets the bounds to around the peak by a user defined amount
        width = self.acceptorWidthSpinBox.value()
        if width > self.acceptortop:
            width = self.acceptortop
            self.acceptorWidthSpinBox.setValue(width)
        self.acceptorT0Lims = (self.acceptortop - width, self.acceptortop + width)
        self.acceptorT0Slider.setValue(self.acceptorT0Lims)

    def accept(self):
        #this function makes the dialog accept
        QtWidgets.QDialog.accept(self)

    def reject(self):
        #this function makes the dialog reject
        QtWidgets.QDialog.reject(self)

    def donor_slider_changed(self):
        #this function changes the T0 limits when the donor slider is changed
        self.donorT0Lims = self.donorT0Slider.value() #{(self.taumin):.2e}
        self.donorLowerLabel.setText(f'Lower bound: {(self.donorTime[self.donorT0Lims[0]]):.3e}')
        self.donorUpperLabel.setText(f'Upper bound: {(self.donorTime[self.donorT0Lims[1]]):.3e}')

    def acceptor_slider_changed(self):
        #this function changes the T0 limits when the acceptor slider is changed
        self.acceptorT0Lims = self.acceptorT0Slider.value()
        self.acceptorLowerLabel.setText(f'Lower bound: {(self.acceptorTime[self.acceptorT0Lims[0]]):.3e}')
        self.acceptorUpperLabel.setText(f'Upper Bound: {(self.acceptorTime[self.acceptorT0Lims[1]]):.3e}')

    def show_donor_bounds(self):
        #this function will set vlines at the position of the T0 bounds on the donor FLIM plot
        if self.showDonorCheckBox.isChecked():
            self.donorDecayPlot.wipe()
            self.donorDecayPlot.ax.plot(self.donorTime, self.donorDecay)
            self.donorDecayPlot.ax.set(xlabel='Time (s)', ylabel='Photon Count')
            self.donorDecayPlot.ax.axvline(x=self.donorTime[self.donorT0Lims[0]], color='red')
            self.donorDecayPlot.ax.axvline(x=self.donorTime[self.donorT0Lims[1]], color='red')
            self.donorDecayPlot.fig.canvas.draw_idle()
        else:
            self.donorDecayPlot.wipe()
            self.donorDecayPlot.ax.plot(self.donorTime, self.donorDecay)
            self.donorDecayPlot.ax.set(xlabel='Time (s)', ylabel='Photon Count')
            self.donorDecayPlot.fig.canvas.draw_idle()
    
    def show_acceptor_bounds(self):
        #this function will set Vlines at the position of the T0 bounds on the acceptor FLIM Plot
        if self.showAcceptorCheckBox.isChecked():
            self.secondaryPlot.wipe()
            self.secondaryPlot.ax.plot(self.acceptorTime, self.acceptorDecay)
            self.secondaryPlot.ax.set(xlabel='Time (s)', ylabel='Photon Count')
            self.secondaryPlot.ax.axvline(x=self.acceptorTime[self.acceptorT0Lims[0]], color='red')
            self.secondaryPlot.ax.axvline(x=self.acceptorTime[self.acceptorT0Lims[1]], color='red')
            self.secondaryPlot.fig.canvas.draw_idle()
        else:
            self.secondaryPlot.wipe()
            self.secondaryPlot.ax.plot(self.acceptorTime, self.acceptorDecay)
            self.secondaryPlot.ax.set(xlabel='Time (s)', ylabel='Photon Count')
            self.secondaryPlot.fig.canvas.draw_idle()

    def show_ino_plots(self):
        #this function is used to show the plots for an ISS file.
        self.donorPlotLabel.setText('Donor FLIM decay:')
        self.donorDecay = np.sum(np.sum(self.donor_cube, axis=0), axis=0)
        self.donorTime = np.linspace(0, len(self.donorDecay)*self.FLIM_time_res, len(self.donorDecay))
        self.donorDecayPlot.ax.plot(self.donorTime, self.donorDecay)
        self.donorDecayPlot.ax.set(xlabel='Time (s)', ylabel='Photon Count')
        self.donorDecayPlot.fig.canvas.draw_idle()
        self.secondaryPlotLabel.setText('Spectral intensity:')
        spectral = self.spectral_cube[:, :, self.spectralRange0[0]:self.spectralRange0[1]]
        self.secondaryPlot.ax.plot(np.sum(np.sum(spectral, axis=0), axis=0))
        steps = min(7, self.spectralRange0[1]-self.spectralRange0[0])
        label_ticks = [int(i) for i in np.linspace(self.spectralRange0[0], self.spectralRange0[1], num=steps)]
        tick_positions = [int(i)-self.spectralRange0[0] for i in np.linspace(self.spectralRange0[0], self.spectralRange0[1], num=steps)]
        self.secondaryPlot.ax.set_xticks(tick_positions)
        self.secondaryPlot.ax.set_xticklabels([f'{round(self.INO_spectral_array[i], 1)}' for i in label_ticks])
        self.secondaryPlot.ax.set(xlabel='Wavelength (nm)', ylabel='Photon Count')
        self.secondaryPlot.ax.set_ylim(bottom=0, top=max(np.sum(np.sum(spectral, axis=0), axis=0))+10)
        self.secondaryPlot.fig.canvas.draw_idle()
        self.donortop = np.argmax(self.donorDecay) 
        #slider things
        self.donorT0Slider.setMinimum(0)
        self.donorT0Slider.setMaximum(len(self.donorDecay))
        self.donorT0Slider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.donorT0Slider.setTickInterval(1)
        self.donorT0Slider.setValue(self.donorT0Lims)
        #spin box things
        self.donorWidthSpinBox.setValue(self.donortop - self.donorT0Lims[0])

    def show_iss_plots(self):
        #this function is used to show the plots for an INO file.
        self.donorDecay = np.sum(np.sum(self.donor_cube, axis=0), axis=0)
        self.donorTime = np.linspace(0, len(self.donorDecay)*self.FLIM_time_res, len(self.donorDecay))
        self.donorPlotLabel.setText('Donor FLIM decay:')
        self.donorDecayPlot.ax.plot(self.donorTime, self.donorDecay)
        self.donorDecayPlot.ax.set(xlabel='Time (s)', ylabel='Photon Count')
        self.donorDecayPlot.fig.canvas.draw_idle()
        self.secondaryPlotLabel.setText('Acceptor FLIM decay:')
        self.acceptorDecay = np.sum(np.sum(self.acceptor_cube, axis=0), axis=0)
        self.acceptorTime = np.linspace(0, len(self.acceptorDecay)*self.FLIM_time_res, len(self.acceptorDecay))
        self.secondaryPlot.ax.plot(self.acceptorTime, self.acceptorDecay)
        self.secondaryPlot.ax.set(xlabel='Time (s)', ylabel='Photon Count')
        self.secondaryPlot.fig.canvas.draw_idle()
        self.acceptorT0Label.setEnabled(True)
        self.acceptorT0Slider.setEnabled(True)
        self.acceptorLowerLabel.setEnabled(True)
        self.acceptorUpperLabel.setEnabled(True)
        self.showAcceptorCheckBox.setEnabled(True)
        self.acceptorWidthSpinBox.setEnabled(True)   
        self.acceptorSpinBoxLabel.setEnabled(True)  
        self.acceptorSetButton.setEnabled(True)   
        self.donortop = np.argmax(self.donorDecay)
        self.acceptortop = np.argmax(self.acceptorDecay)
        #dont forget to set the sliders properly
        self.donorT0Slider.setMinimum(0)
        self.donorT0Slider.setMaximum(len(self.donorDecay))
        self.donorT0Slider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.donorT0Slider.setTickInterval(1)
        self.donorT0Slider.setValue(self.donorT0Lims)
        self.acceptorT0Slider.setMinimum(0)
        self.acceptorT0Slider.setMaximum(len(self.acceptorDecay))
        self.acceptorT0Slider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.acceptorT0Slider.setTickInterval(1)
        self.acceptorT0Slider.setValue(self.acceptorT0Lims)
        #set spin boxes
        self.donorWidthSpinBox.setValue(self.donortop - self.donorT0Lims[0])
        self.acceptorWidthSpinBox.setValue(self.acceptortop - self.acceptorT0Lims[0])
        

class Lifetime_Map_Window(QtWidgets.QMainWindow):
    #this class makes a popup window to show the lifetime map for  plots, spectral intensity plots and all that
    def __init__(self, parent=None):
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        super().__init__(parent)
        self.setupUi()

    def setupUi(self):
        self.setObjectName("LifetimeMapWindow")
        self.resize(1200, 870)
        self.centralwidget = QtWidgets.QWidget()
        self.setCentralWidget(self.centralwidget)
        self.lifetimeMapLabel = QtWidgets.QLabel(self.centralwidget)
        self.lifetimeMapLabel.setGeometry(440, 20, 200, 20)
        self.lifetimeMapLabel.setText('ROI Lifetime Map:')
        self.lifetimeMap = Canvas(self.centralwidget)
        self.lifetimeMap.define_axis(left=0.01, bot=0.01, width=0.94, height=0.94)
        self.lifetimeMap.ax.set(xticks=[], yticks=[])
        self.lifetimeMap.setGeometry(440, 40, 750, 750)
        self.lifetime_num = plt.gcf().number
        self.lifetimeMapToolbar = NavigationToolbar(self.lifetimeMap, self.centralwidget)
        self.lifetimeMapToolbar.setGeometry(440, 800, 750, 50)
        self.lifetimeSliderLabel = QtWidgets.QLabel(self.centralwidget)
        self.lifetimeSliderLabel.setGeometry(10, 30, 250, 20)
        self.lifetimeSliderLabel.setText('Custom Lifetime Range (s)')
        self.lifetimeSlider = QRangeSlider(self.centralwidget)
        self.lifetimeSlider.setOrientation(QtCore.Qt.Horizontal)
        self.lifetimeSlider.setGeometry(10, 60, 420, 20)
        self.lifetimeMinLabel = QtWidgets.QLabel(self.centralwidget)
        self.lifetimeMinLabel.setGeometry(15, 90, 120, 20)
        self.lifetimeMaxLabel = QtWidgets.QLabel(self.centralwidget)
        self.lifetimeMaxLabel.setGeometry(320, 90, 120, 20)
        self.changeLifetimeLimsButton = QtWidgets.QPushButton(self.centralwidget)
        self.changeLifetimeLimsButton.setGeometry(20, 125, 55, 35)
        self.changeLifetimeLimsButton.setText('Apply')
        self.donorCheckBox = QtWidgets.QCheckBox(self.centralwidget)
        self.donorCheckBox.setGeometry(10, 175, 140, 25)
        self.donorCheckBox.setLayoutDirection(QtCore.Qt.LayoutDirection(QtCore.Qt.RightToLeft))
        self.donorCheckBox.setText('Donor Lifetimes')
        self.donorCheckBox.setEnabled(False)
        self.acceptorCheckBox = QtWidgets.QCheckBox(self.centralwidget)
        self.acceptorCheckBox.setGeometry(160, 175, 170, 25)
        self.acceptorCheckBox.setLayoutDirection(QtCore.Qt.LayoutDirection(QtCore.Qt.RightToLeft))
        self.acceptorCheckBox.setText('Acceptor Lifetimes')
        self.acceptorCheckBox.setEnabled(False)
        self.histLabel = QtWidgets.QLabel(self.centralwidget)
        self.histLabel.setGeometry(10, 210, 120, 20)
        self.histCanvas = Canvas(self.centralwidget)                
        self.histCanvas.setGeometry(10, 235, 420, 300)
        self.histCanvas.define_axis(left=0.15, bot=0.16, width=0.8, height=0.8)
        self.hist_num = plt.gcf().number
        self.histToolbar = NavigationToolbar(self.histCanvas, self.centralwidget)
        self.histToolbar.setGeometry(10, 540, 420, 30)
        self.binSpinBoxLabel = QtWidgets.QLabel(self.centralwidget)
        self.binSpinBoxLabel.setGeometry(20, 595, 120, 20)
        self.binSpinBoxLabel.setText('Histogram bins:')
        self.binSpinBox = QtWidgets.QSpinBox(self.centralwidget)
        self.binSpinBox.setGeometry(140, 595, 55, 25)
        self.binSpinBox.setMinimum(2)
        self.binSpinBox.setMaximum(9999)
        self.setBinsButton = QtWidgets.QPushButton(self.centralwidget)
        self.setBinsButton.setGeometry(200, 595, 50, 25)
        self.setBinsButton.setText('Set')
        #signal calls
        self.lifetimeSlider.valueChanged.connect(self.lifetime_lims_changed)
        self.changeLifetimeLimsButton.clicked.connect(self.apply_changes)
        self.donorCheckBox.clicked.connect(self.donor_clicked)
        self.acceptorCheckBox.clicked.connect(self.acceptor_clicked)
        self.setBinsButton.clicked.connect(self.set_bins)

    def set_bins(self):
        #this function sets the histogram to have a user defined number of bins!
        plt.figure(num=self.hist_num)
        self.histCanvas.wipe()
        self.histCanvas.ax.hist(self.lifetime_vals, bins=self.binSpinBox.value())
        self.histCanvas.ax.set(xlabel='Lifetimes (s)', ylabel='Count')
        self.histCanvas.fig.canvas.draw_idle()

    def donor_clicked(self):
        #this function makes sure that only one channel is checked at a time
        if self.donorCheckBox.isChecked():
            self.acceptorCheckBox.setChecked(False)
            self.flimCube = self.donorCube
            plt.figure(num=self.lifetime_num)
            self.lifetimeMap.fig.clf()
            self.lifetimeMap.define_axis(left=0.01, bot=0.01, width=0.92, height=0.92)
            self.lifetimeMap.ax.set(xticks=[], yticks=[])
            self.calculate_lifetimes()
            self.update_slider()
        else:
            self.acceptorCheckBox.setChecked(True)

    def acceptor_clicked(self):
        #this function makes sure that only one channel is checked at a time
        if self.acceptorCheckBox.isChecked():
            self.donorCheckBox.setChecked(False)
            self.flimCube = self.acceptorCube
            plt.figure(num=self.lifetime_num)
            self.lifetimeMap.fig.clf()
            self.lifetimeMap.define_axis(left=0.01, bot=0.01, width=0.92, height=0.92)
            self.lifetimeMap.ax.set(xticks=[], yticks=[])
            self.calculate_lifetimes()
            self.update_slider()
        else:
            self.donorCheckBox.setChecked(True)

    def lifetime_lims_changed(self):
        #this function updates the limits on the lifetimes shown
        botlim, toplim = self.lifetimeSlider.value()
        self.taumin = self.unique_lifetimes[botlim]
        self.taumax = self.unique_lifetimes[toplim]
        if self.taumin >= self.taumax:
            self.taumin, self.taumax = self.taumax, self.taumin
        self.lifetimeMinLabel.setText(f'Min: {(self.taumin):.2e}')
        self.lifetimeMaxLabel.setText(f'Max: {(self.taumax):.2e}')

    def apply_changes(self):
        #this function applies the changes made to the lifetime limits and shows the map with the new limits
        plt.figure(num=self.lifetime_num)  #because the axes need to be removed (for the colorbar), and new axes specified, the currrent figure also must be set before creating new axes.
        self.lifetimeMap.fig.clf()
        self.lifetimeMap.define_axis(left=0.01, bot=0.01, width=0.92, height=0.92)
        self.lifetimeMap.ax.set(xticks=[], yticks=[])
        self.lifetimeMap.ax.imshow(self.FLIMimg, 'gray')
        self.lifetimeMap.ax.imshow(self.lifetime_mask, vmin=self.taumin, vmax=self.taumax, cmap='seismic_r', interpolation='nearest')
        divider = make_axes_locatable(self.lifetimeMap.ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        self.lifetimeMap.fig.colorbar(self.lifetimeMap.ax.imshow(self.lifetime_mask, vmin=self.taumin, vmax=self.taumax, cmap='seismic_r', interpolation='nearest'), cax=cax)
        self.lifetimeMap.fig.canvas.draw_idle()

    def calculate_lifetimes(self):
        #this function calculates a lifetime map for the given FLIM cube over the provided ROI map
        decay_profile = np.sum(np.sum(self.flimCube, axis=0), axis=0)
        decaymin = np.argmax(decay_profile)
        decaymax = np.argmin(decay_profile[decaymin:]) + decaymin
        decaymin += 10
        decaymax -= 10
        self.FLIMimg = np.sum(self.flimCube, axis=2)
        self.lifetime_map, self.lifetime_vals = func.ROI_phasor_lifetime(self.segmented_img, self.labels, self.flimCube[:, :, decaymin:decaymax], freq=self.laser_freq, delta_t=self.FLIM_time_res, list_taus=True)
        self.lifetime_mask = np.ma.masked_where(self.segmented_img == 0, self.lifetime_map)
        self.unique_lifetimes = np.unique(self.lifetime_vals)
        self.taumin = np.min(self.unique_lifetimes)
        self.taumax = np.max(self.unique_lifetimes)
        plt.figure(num=self.lifetime_num)
        self.lifetimeMap.ax.imshow(self.FLIMimg, 'gray')
        self.lifetimeMap.ax.imshow(self.lifetime_mask, vmin=self.taumin, vmax=self.taumax, cmap='seismic_r', interpolation='nearest')
        divider = make_axes_locatable(self.lifetimeMap.ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        self.lifetimeMap.fig.colorbar(mappable=self.lifetimeMap.ax.imshow(self.lifetime_mask, vmin=self.taumin, vmax=self.taumax, cmap='seismic_r', interpolation='nearest'), cax=cax)
        self.lifetimeMap.fig.canvas.draw_idle()

    def update_slider(self):
        #this function updates the slider limits and labels to take on current values (called by update_display)
        self.lifetimeSlider.setMinimum(0)
        self.lifetimeSlider.setMaximum(len(self.unique_lifetimes)-1)
        self.lifetimeSlider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.lifetimeSlider.setTickInterval(1)
        self.lifetimeSlider.setValue((0, (len(self.unique_lifetimes)-1)))
        self.lifetimeMinLabel.setText(f'Min: {(self.taumin):.2e}')
        self.lifetimeMaxLabel.setText(f'Max: {(self.taumax):.2e}')

    def update_display(self):
        #this function ensures that labels and sliders are set to what they should be
        self.update_slider()
        if self.file_type == 'iss-tdflim':
            self.donorCheckBox.setEnabled(True)
            self.acceptorCheckBox.setEnabled(True)
            if self.channel == 'donor':
                self.donorCheckBox.setChecked(True)
            elif self.channel == 'secondary':
                self.acceptorCheckBox.setChecked(True)              
        plt.figure(num=self.hist_num)
        self.histLabel.setText('ROI Lifetime Histogram')
        n, _bins, _patches = self.histCanvas.ax.hist(self.lifetime_vals, 'auto')
        self.histCanvas.ax.set(xlabel='Lifetimes (s)', ylabel='Count')
        self.histCanvas.fig.canvas.draw_idle()
        self.binSpinBox.setValue(len(n))


class Multiple_Fluorophore_Window(QtWidgets.QDialog):
    #this class makes a popup window to allow the user to look at ROI intensity for a given set of wavelengths
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi()

    def setupUi(self):
        self.setObjectName("FluorophoreWindow")
        self.resize(1200, 770)
        # self.centralwidget = QtWidgets.QWidget()
        # self.setCentralWidget(self.centralwidget)
        self.spectralLabel = QtWidgets.QLabel(self) #.centralwidget)
        self.spectralLabel.setGeometry(540, 20, 200, 20)
        self.spectralLabel.setText('Spectral Image:')
        self.spectralMap = Canvas(self) #.centralwidget)
        self.spectralMap.setGeometry(540, 40, 650, 650)
        self.spectralMap.no_ticks()
        self.spectralMapToolbar = NavigationToolbar(self.spectralMap, self) #.centralwidget)
        self.spectralMapToolbar.setGeometry(540, 700, 650, 50)
        self.slider1Label = QtWidgets.QLabel(self) #.centralwidget)
        self.slider1Label.setGeometry(10, 40, 150, 20)
        self.slider1Label.setText('Spectral Range 1:')
        self.slider1 = QRangeSlider(self) #.centralwidget)
        self.slider1.setOrientation(QtCore.Qt.Horizontal)
        self.slider1.setGeometry(10, 65, 370, 20)
        self.slider1MinLabel = QtWidgets.QLabel(self) #.centralwidget)
        self.slider1MinLabel.setGeometry(390, 50, 40, 20)
        self.slider1MinLabel.setText('Min:')
        self.slider1MinBox = QtWidgets.QDoubleSpinBox(self) #.centralwidget)
        self.slider1MinBox.setGeometry(430, 50, 80, 25)
        self.slider1MaxLabel = QtWidgets.QLabel(self) #.centralwidget)
        self.slider1MaxLabel.setGeometry(390, 75, 40, 20)
        self.slider1MaxLabel.setText('Max:')
        self.slider1MaxBox = QtWidgets.QDoubleSpinBox(self) #.centralwidget)
        self.slider1MaxBox.setGeometry(430, 75, 80, 25)
        self.slider2Label = QtWidgets.QLabel(self) #.centralwidget)
        self.slider2Label.setGeometry(10, 125, 150, 20)
        self.slider2Label.setText('Spectral Range 2:')
        self.slider2 = QRangeSlider(self) #.centralwidget)
        self.slider2.setOrientation(QtCore.Qt.Horizontal)
        self.slider2.setGeometry(10, 150, 370, 20)
        self.slider2MinLabel = QtWidgets.QLabel(self) #.centralwidget)
        self.slider2MinLabel.setGeometry(390, 135, 40, 20)
        self.slider2MinLabel.setText('Min:')
        self.slider2MinBox = QtWidgets.QDoubleSpinBox(self) #.centralwidget)
        self.slider2MinBox.setGeometry(430, 135, 80, 25)
        self.slider2MaxLabel = QtWidgets.QLabel(self) #.centralwidget)
        self.slider2MaxLabel.setGeometry(390, 160, 40, 20)
        self.slider2MaxLabel.setText('Max:')
        self.slider2MaxBox = QtWidgets.QDoubleSpinBox(self) #.centralwidget)
        self.slider2MaxBox.setGeometry(430, 160, 80, 25)
        self.slider3Label = QtWidgets.QLabel(self) #.centralwidget)
        self.slider3Label.setGeometry(10, 210, 150, 20)
        self.slider3Label.setText('Spectral Range 3:')
        self.slider3 = QRangeSlider(self) #.centralwidget)
        self.slider3.setOrientation(QtCore.Qt.Horizontal)
        self.slider3.setGeometry(10, 235, 370, 20)
        self.slider3MinLabel = QtWidgets.QLabel(self) #.centralwidget)
        self.slider3MinLabel.setGeometry(390, 220, 40, 20)
        self.slider3MinLabel.setText('Min:')
        self.slider3MinBox = QtWidgets.QDoubleSpinBox(self) #.centralwidget)
        self.slider3MinBox.setGeometry(430, 220, 80, 25)
        self.slider3MaxLabel = QtWidgets.QLabel(self) #.centralwidget)
        self.slider3MaxLabel.setGeometry(390, 245, 40, 20)
        self.slider3MaxLabel.setText('Max:')
        self.slider3MaxBox = QtWidgets.QDoubleSpinBox(self) #.centralwidget)
        self.slider3MaxBox.setGeometry(430, 245, 80, 25)
        self.applyRangesButton = QtWidgets.QPushButton(self) #.centralwidget)
        self.applyRangesButton.setGeometry(40, 270, 130, 30)
        self.applyRangesButton.setText('Apply Ranges')
        self.intensityTitle = QtWidgets.QLabel(self) #.centralwidget)
        self.intensityTitle.setGeometry(20, 340, 400, 20)
        self.intensityTitle.setText('Mean Spectral Range Intensities:')
        self.slider1Intensity = QtWidgets.QLabel(self) #.centralwidget)
        self.slider1Intensity.setGeometry(30, 370, 250, 20)
        self.slider1Intensity.setText('Click an ROI to view spectral range intensity values')
        self.slider2Intensity = QtWidgets.QLabel(self) #.centralwidget)
        self.slider2Intensity.setGeometry(30, 400, 200, 20)
        self.slider3Intensity = QtWidgets.QLabel(self) #.centralwidget)
        self.slider3Intensity.setGeometry(30, 430, 200, 20)
        self.saveButton = QtWidgets.QPushButton(self)
        self.saveButton.setGeometry(200, 700, 150, 35)
        self.saveButton.setText('Save Spectral Ranges')
        self.cancelButton = QtWidgets.QPushButton(self)
        self.cancelButton.setGeometry(370, 700, 150, 35)
        self.cancelButton.setText('Cancel')
        self.setTabOrder(self.saveButton, self.applyRangesButton)
        #signal calls
        self.slider1.valueChanged.connect(self.slider1_changed)
        self.slider2.valueChanged.connect(self.slider2_changed)
        self.slider3.valueChanged.connect(self.slider3_changed)
        self.applyRangesButton.clicked.connect(self.apply_slider_lims)
        self.spectralMap.mpl_connect('button_press_event', self.image_click)
        self.saveButton.clicked.connect(self.accept)
        self.cancelButton.clicked.connect(self.reject)

    def accept(self):
        self.spectralRange1 = self.slider1Lims
        self.spectralRange2 = self.slider2Lims
        self.spectralRange3 = self.slider3Lims
        QtWidgets.QDialog.accept(self)

    def reject(self):
        QtWidgets.QDialog.reject(self)

    def image_click(self, event):
        col = int(np.floor(event.xdata))
        row = int(np.floor(event.ydata))
        roi = self.segmented_img[row, col]
        if roi != 0:
            self.intensityTitle.setText(f'Mean Spectral Range Intensities:  ROI {roi}')
            self.slider1Intensity.setText(f'Spectral range 1: {round(np.mean(self.spectral_img1[self.segmented_img == roi]), 3)}')
            self.slider2Intensity.setText(f'Spectral range 2: {round(np.mean(self.spectral_img2[self.segmented_img == roi]), 3)}')
            self.slider3Intensity.setText(f'spectral range 3: {round(np.mean(self.spectral_img3[self.segmented_img == roi]), 3)}')
        else:
            self.intensityTitle.setText(f'Mean Spectral Range Intensities:  NO ROI SELECTED')

    def slider1_changed(self):
        #this function updates the spectral limits from slider1
        self.slider1Lims = self.slider1.value()
        self.slider1MinBox.setValue(self.spectral_array[self.slider1Lims[0]])
        self.slider1MaxBox.setValue(self.spectral_array[self.slider1Lims[1]])

    def slider2_changed(self):
        #this function updates the spectral limits from slider2
        self.slider2Lims = self.slider2.value()
        self.slider2MinBox.setValue(self.spectral_array[self.slider2Lims[0]])
        self.slider2MaxBox.setValue(self.spectral_array[self.slider2Lims[1]])

    def slider3_changed(self):
        #this function updates the spectral limits from slider3
        self.slider3Lims = self.slider3.value()
        self.slider3MinBox.setValue(self.spectral_array[self.slider3Lims[0]])
        self.slider3MaxBox.setValue(self.spectral_array[self.slider3Lims[1]])
    
    def apply_slider_lims(self):
        #this function updates the spectral limits to those shown on the spin boxes.
        min1 = self.slider1MinBox.value()
        max1 = self.slider1MaxBox.value()
        if max1 < min1:
            min1, max1 = max1, min1
        if min1 < self.spectral_array[0]:
            min1 = self.spectral_array[0]
        if max1 > self.spectral_array[-1]:
            max1 = self.spectral_array[-1]
        min_index = func.closest_index(min1, self.spectral_array)
        max_index = func.closest_index(max1, self.spectral_array)
        if max_index == min_index:
            max_index += 1
        self.slider1Lims = (min_index, max_index)
        self.slider1.setValue(self.slider1Lims)   #changing the slider here forces the spin boxes to go to the correct spectral array value
        self.spectral_img1 = np.sum(self.spectralCube[:, :, self.slider1Lims[0]:self.slider1Lims[1]], axis=2)
        min2 = self.slider2MinBox.value()
        max2 = self.slider2MaxBox.value()
        if max2 < min2:
            min2, max2 = max2, min2
        if min2 < self.spectral_array[0]:
            min2 = self.spectral_array[0]
        if max2 > self.spectral_array[-1]:
            max2 = self.spectral_array[-1]
        min_index = func.closest_index(min2, self.spectral_array)
        max_index = func.closest_index(max2, self.spectral_array)
        if max_index == min_index:
            max_index += 1
        self.slider2Lims = (min_index, max_index)
        self.slider2.setValue(self.slider2Lims)
        self.spectral_img2 = np.sum(self.spectralCube[:, :, self.slider2Lims[0]:self.slider2Lims[1]], axis=2)
        min3 = self.slider3MinBox.value()
        max3 = self.slider3MaxBox.value()
        if max3 < min3:
            min3, max3 = max3, min3
        if min3 < self.spectral_array[0]:
            min3 = self.spectral_array[0]
        if max3 > self.spectral_array[-1]:
            max3 = self.spectral_array[-1]
        min_index = func.closest_index(min3, self.spectral_array)
        max_index = func.closest_index(max3, self.spectral_array)
        if max_index == min_index:
            max_index += 1
        self.slider3Lims = (min_index, max_index)
        self.slider3.setValue(self.slider3Lims)
        self.spectral_img3 = np.sum(self.spectralCube[:, :, self.slider3Lims[0]:self.slider3Lims[1]], axis=2)

    def updateUi(self):
        #image/map
        self.spectralMap.ax.imshow(np.sum(self.spectralCube, axis=2), 'gray')
        self.segmented_mask = np.ma.masked_where(self.segmented_img == 0, self.segmented_img, copy=True)
        self.spectralMap.ax.imshow(self.segmented_mask, cmap=func.random_cmap(len(self.labels)), interpolation='none')
        self.spectralMap.fig.canvas.draw_idle()
        #sliders/spectral array
        self.slider1MinBox.setMinimum(0)
        self.slider1MinBox.setMaximum(99999)
        self.slider1MaxBox.setMinimum(0)
        self.slider1MaxBox.setMaximum(99999)
        self.slider2MinBox.setMinimum(0)
        self.slider2MinBox.setMaximum(99999)
        self.slider2MaxBox.setMinimum(0)
        self.slider2MaxBox.setMaximum(99999)
        self.slider3MinBox.setMinimum(0)
        self.slider3MinBox.setMaximum(99999)
        self.slider3MaxBox.setMinimum(0)
        self.slider3MaxBox.setMaximum(99999)
        sliderLen = len(self.spectral_array)-1
        self.slider1.setMinimum(0)
        self.slider1.setMaximum(sliderLen)
        self.slider1.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.slider1.setTickInterval(1)
        self.slider1Lims = (0, sliderLen//3 )
        self.slider1.setValue(self.slider1Lims)
        self.spectral_img1 = np.sum(self.spectralCube[:, :, self.slider1Lims[0]:self.slider1Lims[1]], axis=2)
        self.slider2.setMinimum(0)
        self.slider2.setMaximum(sliderLen)
        self.slider2.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.slider2.setTickInterval(1)
        self.slider2Lims = (sliderLen//3, 2*sliderLen//3)
        self.slider2.setValue(self.slider2Lims)
        self.spectral_img2 = np.sum(self.spectralCube[:, :, self.slider2Lims[0]:self.slider2Lims[1]], axis=2)
        self.slider3.setMinimum(0)
        self.slider3.setMaximum(sliderLen)
        self.slider3.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.slider3.setTickInterval(1)
        self.slider3Lims = (2*sliderLen//3, sliderLen)
        self.slider3.setValue(self.slider3Lims)
        self.spectral_img3 = np.sum(self.spectralCube[:, :, self.slider3Lims[0]:self.slider3Lims[1]], axis=2)


class Segmentation_Dialog(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)

    def setupUi(self):
        self.resize(460, 470)
        self.algoLabel = QtWidgets.QLabel(self)  #label
        self.algoLabel.setGeometry(10, 10, 200, 20)
        self.watershedBox = QtWidgets.QCheckBox(self)      #checkbox
        self.watershedBox.setGeometry(20, 35, 180, 20)
        self.otsuBox = QtWidgets.QCheckBox(self)      #checkbox
        self.otsuBox.setGeometry(240, 35, 120, 20)
        self.kernelLabel = QtWidgets.QLabel(self)
        self.kernelLabel.setGeometry(10, 70, 100, 20)
        self.kernelSpinBox = QtWidgets.QSpinBox(self)
        self.kernelSpinBox.setGeometry(30, 95, 50, 25)
        self.sigmaLabel = QtWidgets.QLabel(self)
        self.sigmaLabel.setGeometry(140, 70, 200, 20)
        self.sigmaSpinBox = QtWidgets.QDoubleSpinBox(self)
        self.sigmaSpinBox.setGeometry(160, 95, 55, 25)
        self.noiseLabel = QtWidgets.QLabel(self)   #label
        self.noiseLabel.setGeometry(10, 140, 90, 20)
        self.noiseSpinBox = QtWidgets.QSpinBox(self) #spin box
        self.noiseSpinBox.setGeometry(30, 165, 50, 25)
        self.smallROILabel = QtWidgets.QLabel(self) #label
        self.smallROILabel.setGeometry(120, 140, 110, 20)
        self.smallROISpinBox = QtWidgets.QSpinBox(self) #spin box
        self.smallROISpinBox.setGeometry(140, 165, 50, 25)
        self.ROIsizeLabel = QtWidgets.QLabel(self)
        self.ROIsizeLabel.setGeometry(260, 140, 180, 20)
        self.ROIsizeSpinBox = QtWidgets.QSpinBox(self)
        self.ROIsizeSpinBox.setGeometry(280, 165, 50, 25)
        self.logBinaryBox = QtWidgets.QCheckBox(self)      #check Box
        self.logBinaryBox.setGeometry(10, 210, 140, 20)
        self.thresholdROIBox = QtWidgets.QCheckBox(self)#check box
        self.thresholdROIBox.setGeometry(10, 255, 170, 20)
        self.meanScalarLabel = QtWidgets.QLabel(self)  #label
        self.meanScalarLabel.setGeometry(30, 285, 200, 20)
        self.meanScalarSpinBox = QtWidgets.QDoubleSpinBox(self)# double spin box
        self.meanScalarSpinBox.setGeometry(50, 310, 55, 25)
        self.channelLabel = QtWidgets.QLabel(self)
        self.channelLabel.setGeometry(10, 350, 100, 20)
        self.donorCheckBox = QtWidgets.QCheckBox(self)
        self.donorCheckBox.setGeometry(20, 375, 120, 20)
        self.secondaryCheckBox = QtWidgets.QCheckBox(self)
        self.secondaryCheckBox.setGeometry(160, 375, 245, 20)
        self.OKButton = QtWidgets.QPushButton(self)
        self.OKButton.setGeometry(200, 430, 112, 34)
        self.CancelButton = QtWidgets.QPushButton(self)
        self.CancelButton.setGeometry(322, 430, 112, 34)

        self.retranslateUi()
        self.CancelButton.clicked.connect(self.reject)
        self.OKButton.clicked.connect(self.accept)
        self.watershedBox.clicked.connect(self.watershed_click)
        self.otsuBox.clicked.connect(self.otsu_click)
        self.thresholdROIBox.clicked.connect(self.threshold_click)
        self.donorCheckBox.clicked.connect(self.donor_click)
        self.secondaryCheckBox.clicked.connect(self.secondary_click)

    def donor_click(self):
        #this function makes sure that only one channel is checked at a time
        if self.donorCheckBox.isChecked():
            self.secondaryCheckBox.setChecked(False)
        else:
            self.secondaryCheckBox.setChecked(True)

    def secondary_click(self):
        #this function makes sure that only one channel is checked at a time
        if self.secondaryCheckBox.isChecked():
            self.donorCheckBox.setChecked(False)
        else:
            self.donorCheckBox.setChecked(True)

    def watershed_click(self):
        #this function makes sure that only one segmentation algorithm is checked at a time
        if self.watershedBox.isChecked():
            self.otsuBox.setChecked(False)
            self.ROIsizeLabel.setEnabled(True)
            self.ROIsizeSpinBox.setEnabled(True)
            self.logBinaryBox.setEnabled(True)
            self.noiseSpinBox.setEnabled(True)
            self.noiseLabel.setEnabled(True)
        else:
            self.otsuBox.setChecked(True)
            self.ROIsizeLabel.setEnabled(False)
            self.ROIsizeSpinBox.setEnabled(False)
            self.logBinaryBox.setEnabled(False)
            self.noiseSpinBox.setEnabled(False)
            self.noiseLabel.setEnabled(False)

    def otsu_click(self):
        #this function makes sure that only one segmentation algorithm is checked at a time
        if self.otsuBox.isChecked():
            self.watershedBox.setChecked(False)
            self.ROIsizeLabel.setEnabled(False)
            self.ROIsizeSpinBox.setEnabled(False)
            self.logBinaryBox.setEnabled(False)
            self.noiseSpinBox.setEnabled(False)
            self.noiseLabel.setEnabled(False)
        else:
            self.watershedBox.setChecked(True)
            self.ROIsizeLabel.setEnabled(True)
            self.ROIsizeSpinBox.setEnabled(True)
            self.logBinaryBox.setEnabled(True)
            self.noiseSpinBox.setEnabled(True)
            self.noiseLabel.setEnabled(True)

    def threshold_click(self):
        if self.thresholdROIBox.isChecked():
            self.meanScalarLabel.setEnabled(True)
            self.meanScalarSpinBox.setEnabled(True)
        else:
            self.meanScalarLabel.setEnabled(False)
            self.meanScalarSpinBox.setEnabled(False)

    def accept(self):
        if self.watershedBox.isChecked():
            self.segParams['algorithm'] = 'watershed'
        elif self.otsuBox.isChecked():
            self.segParams['algorithm'] = 'otsu'
        self.segParams['kernel_size'] = self.kernelSpinBox.value()
        self.segParams['sigma'] = self.sigmaSpinBox.value()
        self.segParams['small_roi_radius'] = self.smallROISpinBox.value()
        self.segParams['noise_level'] = self.noiseSpinBox.value()
        self.segParams['roi_size'] = self.ROIsizeSpinBox.value()
        if self.logBinaryBox.isChecked():
            self.segParams['log_binary'] = True
        else:
            self.segParams['log_binary'] = False
        if self.thresholdROIBox.isChecked():
            self.segParams['threshold'] = True
        else:
            self.segParams['threshold'] = False
        self.segParams['mean_scalar'] = self.meanScalarSpinBox.value()
        if self.donorCheckBox.isChecked():
            self.segParams['channel'] = 'donor'
        elif self.secondaryCheckBox.isChecked():
            self.segParams['channel'] = 'secondary'
        QtWidgets.QDialog.accept(self)
    
    def reject(self):
        QtWidgets.QDialog.reject(self)

    def retranslateUi(self):
        self.setWindowTitle( "DialogWindow")
        self.algoLabel.setText( 'Segmentation Algorithm:')
        self.watershedBox.setLayoutDirection(QtCore.Qt.LayoutDirection(QtCore.Qt.RightToLeft))
        self.watershedBox.setText( 'Multiseed Watershed') 
        self.watershedBox.setChecked(False) #make sure no boxes are checked ahead of time
        self.otsuBox.setLayoutDirection(QtCore.Qt.LayoutDirection(QtCore.Qt.RightToLeft))
        self.otsuBox.setText( 'Otsu (quick)')
        self.otsuBox.setChecked(False)       #same as above
        self.kernelLabel.setText('Kernel Size:')
        self.kernelLabel.setToolTip('Side length of the Laplacian of Gaussian kernel used to make binary maps.')
        self.kernelSpinBox.setMinimum(0)
        self.kernelSpinBox.setMaximum(999999)
        self.kernelSpinBox.setValue(self.segParams['kernel_size'])
        self.sigmaLabel.setText('kernel Standard Deviation:')
        self.sigmaLabel.setToolTip('Standard deviation of Gaussian for Laplacian of Gaussian kernel.')
        self.sigmaSpinBox.setMinimum(0)
        self.sigmaSpinBox.setMaximum(999999)
        self.sigmaSpinBox.setValue(self.segParams['sigma'])
        self.noiseLabel.setText( 'Noise Level:')
        self.noiseLabel.setToolTip( 'Pixels with photon count lower than "Noise Level" will be treated as background.')
        self.noiseSpinBox.setMinimum(0)
        self.noiseSpinBox.setMaximum(999999)
        self.noiseSpinBox.setValue(self.segParams['noise_level'])
        self.smallROILabel.setText( 'Min ROI Size:')
        self.smallROILabel.setToolTip( 'Bright spots with diameter smaller than "minimum ROI size" will be treated as background.')
        self.smallROISpinBox.setMinimum(0)
        self.smallROISpinBox.setMaximum(999999)
        self.smallROISpinBox.setValue(self.segParams['small_roi_radius'])
        self.ROIsizeLabel.setText('Seed separation:')
        self.ROIsizeLabel.setToolTip('Minimum distance between local maxima for multiseed watershed.')
        self.ROIsizeSpinBox.setMinimum(0)
        self.ROIsizeSpinBox.setMaximum(999999)
        self.ROIsizeSpinBox.setValue(self.segParams['roi_size'])
        self.logBinaryBox.setText( 'log Binary map')
        self.logBinaryBox.setChecked(self.segParams['log_binary'])
        self.logBinaryBox.setLayoutDirection(QtCore.Qt.LayoutDirection(QtCore.Qt.RightToLeft))
        self.logBinaryBox.setToolTip( 'Binary map made from log of image instead of raw image.')
        self.thresholdROIBox.setText( 'Threshold ROI Map')
        self.thresholdROIBox.setChecked(self.segParams['threshold'])
        self.thresholdROIBox.setLayoutDirection(QtCore.Qt.LayoutDirection(QtCore.Qt.RightToLeft))
        self.thresholdROIBox.setToolTip( 'Apply mean thresholding within each ROI.')
        self.meanScalarLabel.setText( 'ROI Threshold Mean Factor:')
        self.meanScalarLabel.setToolTip( 'Factor of mean ROI value to use as threshold.')
        self.meanScalarSpinBox.setMinimum(0)
        self.meanScalarSpinBox.setMaximum(100)
        self.meanScalarSpinBox.setSingleStep(0.01)
        self.meanScalarSpinBox.setValue(self.segParams['mean_scalar'])
        self.channelLabel.setText('Channel:')
        self.donorCheckBox.setText('Donor FLIM')
        self.donorCheckBox.setLayoutDirection(QtCore.Qt.LayoutDirection(QtCore.Qt.RightToLeft))
        self.donorCheckBox.setChecked(False)
        self.secondaryCheckBox.setText('Acceptor FLIM / Hyperspectral')
        self.secondaryCheckBox.setLayoutDirection(QtCore.Qt.LayoutDirection(QtCore.Qt.RightToLeft))
        self.secondaryCheckBox.setChecked(False)
        self.OKButton.setText('OK')
        self.CancelButton.setText('Cancel')
        if self.segParams['algorithm'] == 'watershed':
            self.watershedBox.setChecked(True)
        elif self.segParams['algorithm'] == 'otsu':
            self.otsuBox.setChecked(True)
            self.ROIsizeLabel.setEnabled(False)
            self.ROIsizeSpinBox.setEnabled(False)
            self.noiseLabel.setEnabled(False)
            self.noiseSpinBox.setEnabled(False)
            self.logBinaryBox.setEnabled(False)
        if self.segParams['channel'] == 'donor':
            self.donorCheckBox.setChecked(True)
        elif self.segParams['channel'] == 'secondary':
            self.secondaryCheckBox.setChecked(True)
        self.meanScalarLabel.setEnabled(self.segParams['threshold'])
        self.meanScalarSpinBox.setEnabled(self.segParams['threshold'])


class Config_Window(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi()

    def setupUi(self):
        self.resize(300, 150)
        self.inputLabel = QtWidgets.QLabel(self)
        self.inputLabel.setGeometry(20, 20, 150, 20)
        self.inputLabel.setText('Set file name:')
        self.inputLine = QtWidgets.QLineEdit(self)
        self.inputLine.setGeometry(20, 45, 260, 25)
        self.inputLine.setText('new_config.json')
        self.saveButton = QtWidgets.QPushButton(self)
        self.saveButton.setGeometry(70, 110, 100, 35)
        self.saveButton.setText('Save')
        self.cancelButton = QtWidgets.QPushButton(self)
        self.cancelButton.setGeometry(180, 110, 100, 35)
        self.cancelButton.setText('Cancel')
        #signals
        self.saveButton.clicked.connect(self.accept)
        self.cancelButton.clicked.connect(self.reject)

    def accept(self):
        #this function accepts the dialog and saves the input line
        self.file_name = self.inputLine.text()
        if self.file_name[-5:] != '.json':
            self.inputLine.setText('file name must end with ".json"')
        else:
            QtWidgets.QDialog.accept(self)
        
    def reject(self):
        #this function rejects the dialog
        QtWidgets.QDialog.reject(self)


class Batch_Process_Window(QtWidgets.QMainWindow):
    #this class will hold a batch processing window.
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi()
        self.line = 0
        self.old_message = None
    
    def setupUi(self):
        self.resize(800, 700)
        self.centralwidget = QtWidgets.QWidget()
        self.setCentralWidget(self.centralwidget)
        self.pathLabel = QtWidgets.QLabel(self.centralwidget)
        self.pathLabel.setGeometry(20, 20, 650, 20)
        self.jsonLabel = QtWidgets.QLabel(self.centralwidget)
        self.jsonLabel.setGeometry(20, 60, 200, 20)
        self.jsonLabel.setText('Configuration JSON file:')
        self.jsonLineEdit = QtWidgets.QLineEdit(self.centralwidget)
        self.jsonLineEdit.setGeometry(20, 85, 220, 25)
        self.csvLabel = QtWidgets.QLabel(self.centralwidget)
        self.csvLabel.setGeometry(20, 130, 180, 20)
        self.csvLabel.setText('Output CSV file name:')
        self.csvLineEdit = QtWidgets.QLineEdit(self.centralwidget)
        self.csvLineEdit.setGeometry(20, 155, 220, 25)
        self.autoT0BoxLabel = QtWidgets.QLabel(self.centralwidget)
        self.autoT0BoxLabel.setGeometry(50, 200, 180, 20)
        self.autoT0BoxLabel.setText('Automatic T0 bounds')
        self.autoT0BoxLabel.setToolTip('Use automatic T0 limits, useful when decay peaks may shift between images.')
        self.autoT0Box = QtWidgets.QCheckBox(self.centralwidget)
        self.autoT0Box.setGeometry(20, 200, 20, 20)
        self.T0BoundsLabel = QtWidgets.QLabel(self.centralwidget)
        self.T0BoundsLabel.setGeometry(40, 230, 200, 20)
        self.T0BoundsLabel.setText('Width from decay peak:')
        self.T0BoundsLabel.setEnabled(False)
        self.T0BoundSpinBox = QtWidgets.QSpinBox(self.centralwidget)
        self.T0BoundSpinBox.setGeometry(60, 255, 40, 25)
        self.T0BoundSpinBox.setValue(1)
        self.T0BoundSpinBox.setEnabled(False)
        self.noAcceptorBoxLabel = QtWidgets.QLabel(self.centralwidget)
        self.noAcceptorBoxLabel.setGeometry(50, 300, 200, 20)
        self.noAcceptorBoxLabel.setText('Custom donor lifetime')
        self.noAcceptorBoxLabel.setToolTip('Set lifetime of donor fluorophore when not in presence of acceptor.')
        self.noAcceptorCheckBox = QtWidgets.QCheckBox(self.centralwidget)
        self.noAcceptorCheckBox.setGeometry(20, 300, 20, 20)
        self.noAcceptorLabel = QtWidgets.QLabel(self.centralwidget)
        self.noAcceptorLabel.setGeometry(40, 330, 200, 20)
        self.noAcceptorLabel.setText('Donor Lifetime (seconds)')
        self.noAcceptorLabel.setEnabled(False)
        self.noAcceptorLineEdit = QtWidgets.QLineEdit(self.centralwidget)
        self.noAcceptorLineEdit.setGeometry(60, 355, 70, 25)
        self.noAcceptorLineEdit.setText('3.8e-9')
        self.noAcceptorLineEdit.setEnabled(False)
        self.processButton = QtWidgets.QPushButton(self.centralwidget)
        self.processButton.setGeometry(30, 400, 130, 35)
        self.processButton.setText('Process Stack')
        self.progressTitle = QtWidgets.QLabel(self.centralwidget)
        self.progressTitle.setGeometry(400, 100, 200, 20)
        self.progressTitle.setText('ANALYSIS PROGRESS:')
        self.progressLabel0 = QtWidgets.QLabel(self.centralwidget)
        self.progressLabel0.setGeometry(420, 130, 200, 25)
        self.progressLabel0.setText('Process not started')
        #signal calls
        self.autoT0Box.clicked.connect(self.auto_bounds)
        self.noAcceptorCheckBox.clicked.connect(self.custom_lifetime)
        self.processButton.clicked.connect(self.check_params)

    def custom_lifetime(self):
        #this function enables custom no acceptor donor lifetime
        if self.noAcceptorCheckBox.isChecked():
            self.noAcceptorLabel.setEnabled(True)
            self.noAcceptorLineEdit.setEnabled(True)
        else:
            self.noAcceptorLineEdit.setEnabled(False)
            self.noAcceptorLabel.setEnabled(True)

    def auto_bounds(self):
        #this function enables the automatic T0 bounds control
        if self.autoT0Box.isChecked():
            self.T0BoundsLabel.setEnabled(True)
            self.T0BoundSpinBox.setEnabled(True)
        else:
            self.T0BoundsLabel.setEnabled(False)
            self.T0BoundSpinBox.setEnabled(False)

    def check_params(self):
        #this function checks that all the given params are ok
        #check config file
        # all_params_good = True
        json_name = self.jsonLineEdit.text()
        if json_name[-5:] != '.json':
            self.jsonLineEdit.setText('File name must end with ".json"')
            return None
        else:
            try:
                file = open(json_name)
                self.config = json.load(file)
                file.close()
            except Exception as e:
                print('JSON loading error')
                print(e)
                self.jsonLineEdit.setText('Failed to load configuration file')
                return None
        if self.autoT0Box.isChecked():
            self.config['autoT0'] = True
            self.config['donorT0Lims'] = (self.T0BoundSpinBox.value(), self.T0BoundSpinBox.value())
            if self.config['acceptorT0Lims'] != None:
                self.config['acceptorT0Lims'] = (self.T0BoundSpinBox.value(), self.T0BoundSpinBox.value())
        if self.noAcceptorCheckBox.isChecked():
            try:
                no_acceptor_lifetime = float(self.noAcceptorLineEdit.text())
                self.config['donor_lifetime'] = no_acceptor_lifetime
            except:
                self.noAcceptorLineEdit.setText('Value Error')
                return None
        else:
            self.config['donor_lifetime'] = 3.8e-9
        csv_name = self.csvLineEdit.text()
        if csv_name[-4:] != '.csv':
            self.csvLineEdit.setText('File name must end with ".csv"')
            return None
        else:
            self.config['csv_name'] = csv_name
        self.progressLabel0.setText('Analysis started ... ')
        self.process_stack()
    
    def process_stack(self):
        #this function calls the nonstructured data functions when the params are ok.
        if self.config['file_type'] == 'iss-tdflim':
            df = func.ISS_stack_nonstructured_data_withFRET(self.stack_path, self.ISS_stack, self.config['csv_name'], self.config['algorithm'], self.config['kernel_size'], self.config['sigma'], self.config['noise_level'], self.config['small_roi_radius'], self.config['roi_size'], self.config['log_binary'], self.config['mean_scalar'], self.config['donorT0Lims'], self.config['acceptorT0Lims'], autoT0=self.config['autoT0'], ROI_thresholding=self.config['threshold'], no_acceptor_lifetime=self.config['donor_lifetime'])                                                      
        elif self.config['file_type'] == 'ino-folder':
            df = func.INO_stack_nonstructured_data_withFRET(self.stack_path, self.config['csv_name'], self.config['spectralRange0'], self.config['spectralRange1'], self.config['spectralRange2'], self.config['spectralRange3'], self.config['algorithm'], self.config['kernel_size'], self.config['sigma'], self.config['noise_level'], self.config['small_roi_radius'], self.config['roi_size'], self.config['log_binary'], self.config['mean_scalar'], self.config['donorT0Lims'], autoT0=self.config['autoT0'], ROI_thresholding=self.config['threshold'], no_acceptor_lifetime=self.config['donor_lifetime'])
        self.progressLabel0.setText('Finished! CSV file saved.')

    def retranslateUi(self):
        #this function writes correct values into line edits.F
        self.pathLabel.setText(f'Current Image Stack: {self.stack_path}')
        if self.config_name == None:
            self.jsonLineEdit.setText('ENTER configuration file')
        else:
            self.jsonLineEdit.setText(self.config_name)
        period = self.stack_path.rfind('.')
        if period == -1:
            self.csvLineEdit.setText(self.stack_path + '.csv')
        else:
            self.csvLineEdit.setText(self.stack_path[:period] + '.csv')


class Ui_MainWindow(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        super(Ui_MainWindow, self).__init__(parent)
        self.new_file = False #bool. true if new file is loaded, false if some error
        self.current_path = None #path to current image file
        self.file_type = None    #'iss-tdflim' ; 'ino-file' ; 'ino-folder'
        self.all_current_file = None #all info for current image file
        self.ISS_stack = None #iss_image stack
        self.INO_stack = None #list of INO files in folder
        self.slice_list = []
        self.stack_index = 0 #current index in image stack
        self.ISS_cubes = None #current iss cubes [acceptor flim, donor flim]
        self.INO_cubes = None #current INO cubes (donor flim, hyperspectral)
        self.spectralRange0 = [0, 64]
        self.segParams = ({'channel':'donor', 'algorithm':'watershed', 'kernel_size':50, 'sigma':0.25, 'small_roi_radius':5, 'noise_level':20, 'roi_size':10, 'log_binary':False, 'threshold':True, 'mean_scalar':1})
        self.INO_spectral_array = None
        self.spectralRange1 = None
        self.spectralRange2 = None
        self.spectralRange3 = None 
        self.acceptorT0Lims = None
        self.donorT0Lims = None
        self.jsonName = None
        self.setupUi()

    def setupUi(self):
        self.setObjectName("MainWindow")
        self.resize(1600, 900)
        self.centralwidget = QtWidgets.QWidget()
        self.setCentralWidget(self.centralwidget)
        self.stackList = QtWidgets.QListWidget(self.centralwidget)
        self.stackList.setGeometry(10, 110, 220, 400 )
        font = QtGui.QFont()
        font.setPointSize(8)
        self.stackList.setFont(font)
        item = QtWidgets.QListWidgetItem()
        self.stackList.addItem(item)
        self.fileslabel = QtWidgets.QLabel(self.centralwidget)
        self.fileslabel.setGeometry(10, 15, 200, 20)
        self.stackPath = QtWidgets.QLineEdit(self.centralwidget)
        self.stackPath.setGeometry(10, 40, 220, 20)
        self.donorFLIMimage = Canvas(self.centralwidget)
        self.donorFLIMimage.setGeometry(250, 110, 650, 650)
        self.donorFLIMimage.ax.imshow(plt.imread('none.png'))
        self.donorFLIMimage.no_ticks()
        self.donorFLIMimage.fig.canvas.draw_idle()
        self.donorFLIMToolbar = NavigationToolbar(self.donorFLIMimage, self.centralwidget)
        self.donorFLIMToolbar.setGeometry(240, 765, 650, 50)
        self.donorFLIM_name = QtWidgets.QLabel(self.centralwidget)
        self.donorFLIM_name.setEnabled(True)
        self.donorFLIM_name.setText('Donor FLIM:')
        self.donorFLIM_name.setGeometry(250, 88, 600, 21)
        self.secondaryImage = Canvas(self.centralwidget)
        self.secondaryImage.setGeometry(920, 110, 650, 650)
        self.secondaryImage.ax.imshow(plt.imread('none.png'))
        self.secondaryImage.no_ticks()
        self.secondaryImage.fig.canvas.draw_idle()
        self.secondaryImageToolbar = NavigationToolbar(self.secondaryImage, self.centralwidget)
        self.secondaryImageToolbar.setGeometry(910, 765, 650, 50)
        self.secondary_name = QtWidgets.QLabel(self.centralwidget)
        self.secondary_name.setEnabled(True)
        self.secondary_name.setText('Spectral / Acceptor FLIM:')
        self.secondary_name.setGeometry(920, 88, 600, 21)
        self.nextButton = QtWidgets.QPushButton(self.centralwidget)
        self.nextButton.setGeometry(125, 70, 105, 30)
        self.previousButton = QtWidgets.QPushButton(self.centralwidget)
        self.previousButton.setGeometry(10, 70, 105, 30)
        self.segmentButton = QtWidgets.QPushButton(self.centralwidget)
        self.segmentButton.setGeometry(360, 40, 130, 40)
        self.lifetimeButton = QtWidgets.QPushButton(self.centralwidget)
        self.lifetimeButton.setGeometry(510, 40, 130, 40)
        self.clearRoiButton = QtWidgets.QPushButton(self.centralwidget)
        self.clearRoiButton.setGeometry(660, 40, 130, 40)
        self.spectralRangeLabel = QtWidgets.QLabel(self.centralwidget)
        self.spectralRangeLabel.setGeometry(940, 20, 150, 20)
        self.spectralRangeSlider = QRangeSlider(self.centralwidget)
        self.spectralRangeSlider.setOrientation(QtCore.Qt.Horizontal)
        self.spectralRangeSlider.setGeometry(940, 45, 300, 20)
        self.spectralMinLabel = QtWidgets.QLabel(self.centralwidget)
        self.spectralMinLabel.setGeometry(1255, 30, 40, 20)
        self.spectralMinSpinBox = QtWidgets.QDoubleSpinBox(self.centralwidget)
        self.spectralMinSpinBox.setGeometry(1295, 30, 80, 25)
        self.spectralMaxLabel = QtWidgets.QLabel(self.centralwidget)
        self.spectralMaxLabel.setGeometry(1255, 55, 40, 20)
        self.spectralMaxSpinBox = QtWidgets.QDoubleSpinBox(self.centralwidget)
        self.spectralMaxSpinBox.setGeometry(1295, 55, 80, 25)
        self.applySpectralButton = QtWidgets.QPushButton(self.centralwidget)
        self.applySpectralButton.setGeometry(1390, 40, 70, 30)
        self.multiplefluorophoreButton = QtWidgets.QPushButton(self.centralwidget)
        self.multiplefluorophoreButton.setGeometry(1470, 40, 120, 30)
        self.plotButton = QtWidgets.QPushButton(self.centralwidget)
        self.plotButton.setGeometry(55, 550, 130, 40)
        self.menuBar = QtWidgets.QMenuBar()
        self.menuBar.setGeometry(0,0, 1130, 21)
        self.menuSettings = QtWidgets.QMenu(self.menuBar)
        self.setMenuBar(self.menuBar)
        self.actionSegmentationSettings = QtWidgets.QAction()
        self.actionSaveConfig = QtWidgets.QAction()
        self.actionBatchProcess = QtWidgets.QAction()
        self.menuSettings.addAction(self.actionSegmentationSettings)
        self.menuSettings.addAction(self.actionSaveConfig)
        self.menuSettings.addAction(self.actionBatchProcess)
        self.menuBar.addAction(self.menuSettings.menuAction())
        #signal calls
        self.retranslateUi()
        self.stackPath.returnPressed.connect(self.get_files)
        self.nextButton.clicked.connect(self.next_img)
        self.previousButton.clicked.connect(self.previous_img)
        self.stackList.itemClicked.connect(self.list_item_clicked)
        self.segmentButton.clicked.connect(self.show_roi_map)
        self.clearRoiButton.clicked.connect(self.remove_roi_map)
        self.lifetimeButton.clicked.connect(self.display_lifetime_map)
        self.spectralRangeSlider.valueChanged.connect(self.slider_changed)
        self.applySpectralButton.clicked.connect(self.apply_spectral_lims)
        self.multiplefluorophoreButton.clicked.connect(self.display_multiple_fluorophore)
        self.actionSegmentationSettings.triggered.connect(self.segmentation_settings_window)
        self.plotButton.clicked.connect(self.display_plots)
        self.actionSaveConfig.triggered.connect(self.save_config)
        self.actionBatchProcess.triggered.connect(self.open_batch_process)

    def open_batch_process(self):
        if self.new_file:
            batch_window = Batch_Process_Window(self)
            batch_window.stack_path = self.current_path
            batch_window.config_name = self.jsonName
            if self.file_type == 'iss-tdflim':
                batch_window.ISS_stack = self.ISS_stack
            batch_window.retranslateUi()
            batch_window.show()

    def save_config(self):
        #this function saves the current configuration in a json file
        if self.new_file:
            config_dict = self.segParams
            config_dict ['autoT0'] = False
            config_dict['file_type'] = self.file_type
            config_dict['spectralRange0'] = (int(self.spectralRange0[0]), int(self.spectralRange0[1]))
            if self.spectralRange1 == None:
                config_dict['spectralRange1'] = self.spectralRange1
            else:
                config_dict['spectralRange1'] = (int(self.spectralRange1[0]), int(self.spectralRange1[1]))
            if self.spectralRange2 == None:
                config_dict['spectralRange2'] = self.spectralRange2
            else:
                config_dict['spectralRange2'] = (int(self.spectralRange2[0]), int(self.spectralRange2[1]))
            if self.spectralRange3 == None:
                config_dict['spectralRange3'] = self.spectralRange3
            else:
                config_dict['spectralRange3'] = (int(self.spectralRange3[0]), int(self.spectralRange3[1]))
            config_dict['donorT0Lims'] = (int(self.donorT0Lims[0]), int(self.donorT0Lims[1])) #no if statement because show cubes already calculates the limits (all file types have donor T0)
            if self.acceptorT0Lims == None:
                config_dict['acceptorT0Lims'] = self.acceptorT0Lims
            else:
                config_dict['acceptorT0Lims'] = (int(self.acceptorT0Lims[0]), int(self.acceptorT0Lims[1]))
            fileNameWindow = Config_Window(self)
            result = fileNameWindow.exec()
            if result == 1:
                self.jsonName = fileNameWindow.file_name
                with open(self.jsonName, 'w') as newfile:
                    json.dump(config_dict, newfile, indent=4)
                fileNameWindow.close()
            elif result == 0:
                fileNameWindow.close()

    def display_multiple_fluorophore(self):
        fluorophoreWindow = Multiple_Fluorophore_Window(self)
        fluorophoreWindow.spectralCube = self.INO_cubes[1]
        fluorophoreWindow.segmented_img = self.segmented_img
        fluorophoreWindow.labels = self.labels
        fluorophoreWindow.spectral_array = self.INO_spectral_array
        fluorophoreWindow.updateUi()
        result = fluorophoreWindow.exec()
        if result == 1:
            self.spectralRange1 = fluorophoreWindow.spectralRange1
            self.spectralRange2 = fluorophoreWindow.spectralRange2
            self.spectralRange3 = fluorophoreWindow.spectralRange3
            fluorophoreWindow.close()
        else:
            self.spectralRange1 = None
            self.spectralRange2 = None
            self.spectralRange3 = None
            fluorophoreWindow.close()

    def segmentation_settings_window(self):
        #this function opens a dialog window, so that the user can adjust the ROI segmentation settings
        settingsWindow = Segmentation_Dialog(self)
        settingsWindow.segParams = self.segParams.copy()
        settingsWindow.setupUi()   #setupUi() is here instead of in the __init__ for the class b/c segParams needs to be given first
        result = settingsWindow.exec()
        if result == 1:
            self.segParams = settingsWindow.segParams.copy()
            settingsWindow.close()
        else:
            settingsWindow.close()

    def display_lifetime_map(self):
        #this function calls a pop up window displaying a lifetime map for the ROIs in the most recent ROI map.
        lifetimeWindow = Lifetime_Map_Window(self)
        lifetimeWindow.channel = self.segParams['channel']
        lifetimeWindow.file_type = self.file_type
        lifetimeWindow.segmented_img = self.segmented_img
        lifetimeWindow.labels = self.labels
        lifetimeWindow.laser_freq = self.laser_freq
        lifetimeWindow.FLIM_time_res = self.FLIM_time_res
        if self.segParams['channel'] == 'donor':
            if self.file_type == 'iss-tdflim':
                lifetimeWindow.donorCube = self.ISS_cubes[1]
                lifetimeWindow.acceptorCube = self.ISS_cubes[0]
                lifetimeWindow.flimCube = lifetimeWindow.donorCube
            elif (self.file_type == 'ino-file') or (self.file_type == 'ino-folder'):
                lifetimeWindow.donorCube = self.INO_cubes[0]
                lifetimeWindow.flimCube = lifetimeWindow.donorCube
        elif self.segParams['channel'] == 'secondary':
            if self.file_type == 'iss-tdflim':
                lifetimeWindow.acceptorCube = self.ISS_cubes[0]
                lifetimeWindow.donorCube = self.ISS_cubes[1]
                lifetimeWindow.flimCube = lifetimeWindow.acceptorCube
            elif (self.file_type == 'ino-file') or (self.file_type == 'ino-folder'):
                lifetimeWindow.donorCube = self.INO_cubes[0]
                lifetimeWindow.flimCube = lifetimeWindow.donorCube
        lifetimeWindow.calculate_lifetimes()
        lifetimeWindow.update_display()
        lifetimeWindow.show()

    def display_plots(self):
        #this function calls a popup window to show decay plots and spectral plots
        plot_window = Display_Plots(self)
        if self.file_type == 'iss-tdflim':
            plot_window.FLIM_time_res = self.FLIM_time_res
            plot_window.donor_cube = self.ISS_cubes[1]
            plot_window.acceptor_cube = self.ISS_cubes[0]
            plot_window.donorT0Lims = self.donorT0Lims
            plot_window.acceptorT0Lims = self.acceptorT0Lims
            plot_window.show_iss_plots()
            result = plot_window.exec()
            if result == 1:
                self.donorT0Lims = plot_window.donorT0Lims
                self.acceptorT0Lims = plot_window.acceptorT0Lims
                plot_window.close()
            elif result == 0:
                plot_window.close()
        elif (self.file_type == 'ino-file') or (self.file_type == 'ino-folder'):
            plot_window.FLIM_time_res = self.FLIM_time_res
            plot_window.donor_cube = self.INO_cubes[0]
            plot_window.spectral_cube = self.INO_cubes[1]
            plot_window.INO_spectral_array = self.INO_spectral_array
            plot_window.spectralRange0 = self.spectralRange0
            plot_window.donorT0Lims = self.donorT0Lims
            plot_window.show_ino_plots()
            result = plot_window.exec()
            if result == 1:
                self.donorT0Lims = plot_window.donorT0Lims
                plot_window.close()
            elif result == 0:
                plot_window.close()

    def apply_spectral_lims(self):
        #this function updates the spectral limits when changed on the spin boxes
        new_min = self.spectralMinSpinBox.value()
        new_max = self.spectralMaxSpinBox.value()
        if new_max < new_min:
            new_min, new_max = new_max, new_min
        if new_min < self.INO_spectral_array[0]:
            new_min = self.INO_spectral_array[0]
        if new_max > self.INO_spectral_array[-1]:
            new_max = self.INO_spectral_array[-1]
        min_index = func.closest_index(new_min, self.INO_spectral_array)
        max_index = func.closest_index(new_max, self.INO_spectral_array)
        if max_index == min_index:
            max_index += 1
        self.spectralRange0 = [min_index, max_index]
        self.spectralRangeSlider.setValue([min_index, max_index]) #changing the slider will force the spin boxes to go to the correct spectral array value
        self.show_INO_spectral_range()
        
    def slider_changed(self):
        #this function updates the spectral limits when changed using the slider
        new_min_index, new_max_index = self.spectralRangeSlider.value()
        self.spectralMinSpinBox.setValue(self.INO_spectral_array[new_min_index])
        self.spectralMaxSpinBox.setValue(self.INO_spectral_array[new_max_index])
        self.spectralRange0[0] = new_min_index
        self.spectralRange0[1] = new_max_index

    def show_INO_spectral_range(self):
        #this function shows the spectral INO cube, restricted to the spectral limits defined by the user.
        if self.file_type == 'ino-file':
            spectral = self.INO_cubes[1]
            spectral = spectral[:, :, self.spectralRange0[0]:self.spectralRange0[1]]
        elif self.file_type == 'ino-folder':
            file_name = self.INO_stack[self.stack_index]
            self.INO_cubes = func.read_multipage_tiff(self.current_path + '\\' + file_name)
            spectral = self.INO_cubes[1]
            spectral = spectral[:, :, self.spectralRange0[0]:self.spectralRange0[1]]
        self.secondaryImage.wipe()
        self.secondaryimg = np.sum(spectral, axis=2)
        self.secondaryImage.ax.imshow(self.secondaryimg, 'gray')
        self.secondaryImage.fig.canvas.draw_idle()
        self.secondary_name.setText(f'Spectral Image:   max [{np.max(self.secondaryimg)}]  min [{np.min(self.secondaryimg)}]')

    def show_roi_map(self):
        #this function produces an ROI map for the current donor image.
        #roi map is displayed overlayed onto the donor flim image
        if self.new_file:
            if self.segParams['channel'] == 'donor':
                donor_img = self.donorimg 
                self.segmented_img, self.labels = func.segmentation_choice(donor_img, self.segParams['algorithm'], kernel_size=self.segParams['kernel_size'], sigma=self.segParams['sigma'], small_roi_radius=self.segParams['small_roi_radius'], noise_level=self.segParams['noise_level'], roi_size=self.segParams['roi_size'], log_binary=self.segParams['log_binary'], threshold=self.segParams['threshold'], mean_scalar=self.segParams['mean_scalar'])
                segment_mask = np.ma.masked_where(self.segmented_img == 0, self.segmented_img)
                self.donorFLIMimage.wipe()
                self.donorFLIMimage.ax.imshow(donor_img, 'gray')
                self.donorFLIMimage.ax.imshow(segment_mask, cmap=func.random_cmap(len(self.labels)), interpolation='none')
                self.donorFLIMimage.fig.canvas.draw_idle()
            elif self.segParams['channel'] == 'secondary':
                spectral_img = self.secondaryimg 
                self.segmented_img, self.labels = func.segmentation_choice(spectral_img, self.segParams['algorithm'], kernel_size=self.segParams['kernel_size'], sigma=self.segParams['sigma'], small_roi_radius=self.segParams['small_roi_radius'], noise_level=self.segParams['noise_level'], roi_size=self.segParams['roi_size'], log_binary=self.segParams['log_binary'], threshold=self.segParams['threshold'], mean_scalar=self.segParams['mean_scalar'])
                segment_mask = np.ma.masked_where(self.segmented_img == 0, self.segmented_img)
                self.secondaryImage.wipe()
                self.secondaryImage.ax.imshow(spectral_img, 'gray')
                self.secondaryImage.ax.imshow(segment_mask, cmap=func.random_cmap(len(self.labels)), interpolation='none')
                self.secondaryImage.fig.canvas.draw_idle()
            if len(self.labels) >= 1:
                self.lifetimeButton.setEnabled(True)
            if (self.file_type == 'ino-file') or (self.file_type == 'ino-folder'):
                self.multiplefluorophoreButton.setEnabled(True)
            
    def remove_roi_map(self):
        #this function shows the images again and deletes any old roi maps
        if self.new_file:
            self.show_cubes() 
            self.lifetimeButton.setEnabled(False) 
            self.multiplefluorophoreButton.setEnabled(False)

    def clear_old_stacks(self):
        #this function gets rid of any old image stacks/cubes that might be held on to when a new file is loaded.
        if self.file_type == 'iss-tdflim':
            self.INO_stack = None 
            self.INO_cubes = None 
        elif (self.file_type == 'ino-file') or (self.file_type == 'ino-folder'):
            self.ISS_stack = None 
            self.ISS_cubes = None 
        else:
            pass

    def list_item_clicked(self):
        #this function allows the user to choose which slice is shown 
        if self.new_file:
            self.stack_index = self.stackList.currentRow()
            self.stackList.item(self.stack_index).setSelected(True)
            self.show_cubes()
            self.multiplefluorophoreButton.setEnabled(False)

    def next_img(self):
        if self.new_file:
            #this function increases the stack index by 1. and shows next image (or first image if at end of stack.)
            if self.file_type == 'ino-file':
                self.show_cubes()
            elif self.file_type == 'iss-tdflim':
                if self.stack_index == len(self.ISS_stack)-1:
                    self.stack_index = 0
                else:
                    self.stack_index += 1
                self.show_cubes()
            elif self.file_type == 'ino-folder':
                if self.stack_index == len(self.INO_stack)-1:
                    self.stack_index = 0
                else:
                    self.stack_index +=1
                self.show_cubes()
            self.stackList.item(self.stack_index).setSelected(True)
            self.multiplefluorophoreButton.setEnabled(False)

    def previous_img(self):
        #this function decreases the stack index by 1. and shows previous image (or last if at start of stack)
        if self.new_file:
            if self.file_type == 'ino-file':
                self.show_cubes()
            elif self.file_type == 'iss-tdflim':
                if self.stack_index == 0:
                    self.stack_index = len(self.ISS_stack) -1
                else:
                    self.stack_index -=1
                self.show_cubes()
            elif self.file_type == 'ino-folder':
                if self.stack_index == 0:
                    self.stack_index = len(self.INO_stack)-1
                else:
                    self.stack_index -=1
                self.show_cubes()
            self.stackList.item(self.stack_index).setSelected(True)
            self.multiplefluorophoreButton.setEnabled(False)

    def get_files(self):
        #this function gets collects the data from the file/folder given by the user.
        #defines the current path, current file if single file, file/folder type, image stack, stack index.
        self.stackList.clear()
        self.new_file = False
        try:
            self.current_path = self.stackPath.text()
            if self.current_path[-11:] == '.iss-tdflim': #if single file
                self.file_type = 'iss-tdflim'
                self.clear_old_stacks()
                self.all_current_file = func.open_file(self.current_path)
                self.FLIM_time_res = self.all_current_file[3]
                self.laser_freq = self.all_current_file[4]
                self.ISS_stack = self.all_current_file[1]
                self.stack_index = 0 #set initial iss stack index to 0
                self.new_file = True
                slice_nums = np.arange(0, len(self.ISS_stack), 1)
                self.slice_list = [f'slice {i}' for i in slice_nums]
                self.spectralRangeLabel.setVisible(False)
                self.spectralRangeSlider.setVisible(False)
                self.spectralMaxLabel.setVisible(False)
                self.spectralMinLabel.setVisible(False)
                self.spectralMaxSpinBox.setVisible(False)
                self.spectralMinSpinBox.setVisible(False)
                self.applySpectralButton.setVisible(False)
                self.multiplefluorophoreButton.setVisible(False)
            elif self.current_path[-4:]== '.tif':
                self.file_type = 'ino-file'
                self.clear_old_stacks()
                self.slice_list = [self.current_path[self.current_path.rfind('\\'):]]
                self.all_current_file = func.open_file(self.current_path)
                self.FLIM_time_res = self.all_current_file[3]
                self.INO_cubes = self.all_current_file[1]
                self.laser_freq = self.all_current_file[4]
                INO_spectral_array = func.get_spectral_map(self.current_path)
                self.INO_spectral_array = INO_spectral_array[INO_spectral_array != 0]
                self.new_file = True
                self.spectralRange0[0] = 0
                self.spectralRange0[1] = np.argmax(self.INO_spectral_array)
                self.spectralRangeSlider.setMaximum(self.spectralRange0[1])
                self.spectralRangeSlider.setValue(self.spectralRange0)    
                self.spectralMinSpinBox.setValue(self.INO_spectral_array[0])
                self.spectralMaxSpinBox.setValue(self.INO_spectral_array[np.argmax(self.INO_spectral_array)])
                self.spectralRangeLabel.setVisible(True)
                self.spectralRangeSlider.setVisible(True)
                self.spectralMinLabel.setVisible(True)
                self.spectralMaxLabel.setVisible(True)
                self.spectralMinSpinBox.setVisible(True)
                self.spectralMaxSpinBox.setVisible(True)
                self.applySpectralButton.setVisible(True)
                self.multiplefluorophoreButton.setVisible(True)
            else:
                INO_files = func.list_INO_in_folder(self.current_path)
                if INO_files != []:
                    self.file_type = 'ino-folder'
                    self.clear_old_stacks()
                    self.INO_stack = INO_files #list of ino files in folder
                    self.slice_list = self.INO_stack
                    INO_spectral_array = func.get_spectral_map(self.current_path + '\\' + self.INO_stack[0])
                    self.INO_spectral_array = INO_spectral_array[INO_spectral_array != 0]
                    self.stack_index = 0
                    self.new_file = True
                    self.spectralRange0[0] = 0
                    self.spectralRange0[1] = func.closest_index(np.max(self.INO_spectral_array), self.INO_spectral_array)
                    self.spectralRangeSlider.setMaximum(self.spectralRange0[1])
                    self.spectralRangeSlider.setValue(self.spectralRange0)    
                    self.spectralMinSpinBox.setValue(self.INO_spectral_array[0])
                    self.spectralMaxSpinBox.setValue(np.max(self.INO_spectral_array))
                    self.spectralRangeLabel.setVisible(True)
                    self.spectralRangeSlider.setVisible(True)
                    self.spectralMinLabel.setVisible(True)
                    self.spectralMaxLabel.setVisible(True)
                    self.spectralMinSpinBox.setVisible(True)
                    self.spectralMaxSpinBox.setVisible(True)
                    self.applySpectralButton.setVisible(True)
                    self.multiplefluorophoreButton.setVisible(True)
            if self.new_file:
                self.show_cubes() 
                self.plotButton.setEnabled(True)
            for i, slice in enumerate(self.slice_list):
                item = QtWidgets.QListWidgetItem()
                self.stackList.addItem(item)
                item = self.stackList.item(i)
                item.setText(slice)
            self.stackList.item(self.stack_index).setSelected(True)
        except Exception as e:
            print('Image loading error')
            print(e) 
            self.new_file = False
            self.plotButton.setEnabled(False)
            self.donorFLIMimage.wipe()
            self.donorFLIMimage.ax.imshow(plt.imread('bad_path.png'))
            self.donorFLIMimage.fig.canvas.draw_idle()
            self.secondaryImage.wipe()
            self.secondaryImage.ax.imshow(plt.imread('bad_path.png'))
            self.secondaryImage.fig.canvas.draw_idle()

    def show_cubes(self):
        #this function shows the image cubes at the given index in the loaded image stack
        #this function also defines the current ISS_cubes or INO_cubes attributes
        if self.file_type == 'iss-tdflim':
            self.ISS_cubes = self.ISS_stack[self.stack_index]
            acceptor, donor = self.ISS_cubes                   #define the cubes for plotting
            top = np.argmax(np.sum(np.sum(acceptor, axis=0), axis=0))
            self.acceptorT0Lims = (top-1, top+1)
            top = np.argmax(np.sum(np.sum(donor, axis=0), axis=0))
            self.donorT0Lims = (top-1, top+1)
            self.secondary_name.setText('Acceptor FLIM:')
            self.donorFLIMimage.wipe()
            self.donorimg = np.sum(donor, axis=2)
            self.donorFLIMimage.ax.imshow(self.donorimg, 'gray')
            self.donorFLIMimage.fig.canvas.draw_idle()
            self.secondaryImage.wipe()
            self.secondaryimg = np.sum(acceptor, axis=2)
            self.secondaryImage.ax.imshow(self.secondaryimg, 'gray')
            self.secondaryImage.fig.canvas.draw_idle()
        elif self.file_type == 'ino-file':
            donor, spectral = self.INO_cubes
            top = np.argmax(np.sum(np.sum(donor, axis=0), axis=0))
            self.donorT0Lims = (top-1, top+1)
            self.acceptorT0Lims = None
            self.secondary_name.setText('Spectral Image:')
            self.donorFLIMimage.wipe()
            self.donorimg = np.sum(donor, axis=2)
            self.donorFLIMimage.ax.imshow(self.donorimg, 'gray')
            self.donorFLIMimage.fig.canvas.draw_idle()
            self.secondaryImage.wipe()
            self.secondaryimg = np.sum(spectral[:, :, self.spectralRange0[0]:self.spectralRange0[1]], axis=2)
            self.secondaryImage.ax.imshow(self.secondaryimg, 'gray')
            self.secondaryImage.fig.canvas.draw_idle()
        elif self.file_type == 'ino-folder':
            file_name = self.INO_stack[self.stack_index]
            self.all_current_file = func.open_file(self.current_path + '\\' + file_name)
            self.INO_cubes = self.all_current_file[1]
            self.FLIM_time_res = self.all_current_file[3]
            self.laser_freq = self.all_current_file[4]
            donor, spectral = self.INO_cubes
            top = np.argmax(np.sum(np.sum(donor, axis=0), axis=0))
            self.donorT0Lims = (top-1, top+1)
            self.acceptorT0Lims = None
            self.secondary_name.setText('Spectral Image:')
            self.donorFLIMimage.wipe()
            self.donorimg = np.sum(donor, axis=2)
            self.donorFLIMimage.ax.imshow(self.donorimg, 'gray')
            self.donorFLIMimage.fig.canvas.draw_idle()
            self.secondaryImage.wipe()
            self.secondaryimg = np.sum(spectral[:, :, self.spectralRange0[0]:self.spectralRange0[1]], axis=2)
            self.secondaryImage.ax.imshow(self.secondaryimg, 'gray')
            self.secondaryImage.fig.canvas.draw_idle()
        self.lifetimeButton.setEnabled(False)
        self.donorFLIM_name.setText(f'Donor FLIM:   max [{np.max(self.donorimg)}]   min [{np.min(self.donorimg)}]')
        if self.file_type == 'iss-tdflim':
            self.secondary_name.setText(f'Acceptor FLIM:   max [{np.max(self.secondaryimg)}]   min [{np.min(self.secondaryimg)}]')
        elif (self.file_type == 'ino-file') or (self.file_type == 'ino-folder'):
            self.secondary_name.setText(f'Spectral Image:   max [{np.max(self.secondaryimg)}]   min [{np.min(self.secondaryimg)}]')
 
    def retranslateUi(self):
        self.setWindowTitle("MainWindow")
        item = self.stackList.item(0)
        item.setText("No file/folder selected")
        self.stackList.setSortingEnabled(False)
        self.fileslabel.setText("Image file / folder")
        self.stackPath.setText("ENTER file/folder path")
        self.nextButton.setText("Next")
        self.previousButton.setText("Previous")
        self.segmentButton.setText( 'Segment Image')
        self.lifetimeButton.setText('ROI Lifetimes')
        self.lifetimeButton.setToolTip('show lifetimes for most recent ROI map (any channel)')
        self.lifetimeButton.setEnabled(False)  #button off until roi map is made
        self.clearRoiButton.setText( 'Clear Overlay')
        self.menuSettings.setTitle( 'Options')
        self.actionSegmentationSettings.setText( 'segmentation settings')
        self.actionSegmentationSettings.setToolTip( 'Change ROI segmentation settings (CTRL + 1)')
        self.actionSegmentationSettings.setShortcut( 'Ctrl+1')
        self.spectralRangeLabel.setText('Spectral range')
        self.spectralRangeSlider.setMinimum(0)
        self.spectralRangeSlider.setMaximum(63)
        self.spectralRangeSlider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.spectralRangeSlider.setTickInterval(1)
        self.spectralMinLabel.setText('Min:')
        self.spectralMinSpinBox.setMinimum(0)
        self.spectralMinSpinBox.setMaximum(99999)
        self.spectralMaxSpinBox.setMinimum(0)
        self.spectralMaxSpinBox.setMaximum(99999)
        self.spectralMaxLabel.setText('Max:')
        self.spectralRangeLabel.setVisible(False)
        self.spectralRangeSlider.setVisible(False)
        self.spectralMinLabel.setVisible(False)
        self.spectralMaxLabel.setVisible(False)
        self.spectralMaxSpinBox.setVisible(False)
        self.spectralMinSpinBox.setVisible(False)
        self.applySpectralButton.setText('Apply')
        self.applySpectralButton.setToolTip('apply spectral range')
        self.applySpectralButton.setVisible(False)
        self.multiplefluorophoreButton.setText('ROI Intensities')
        self.multiplefluorophoreButton.setToolTip('View mean ROI intensities for multiple spectral ranges.')
        self.multiplefluorophoreButton.setEnabled(False)
        self.multiplefluorophoreButton.setVisible(False)
        self.plotButton.setText('Decay Plots')
        self.plotButton.setToolTip('Show FLIM decay plots/spectral intensity plots')
        self.plotButton.setEnabled(False)
        self.actionSaveConfig.setText('Save configuration')
        self.actionSaveConfig.setToolTip('Save current configuration in a json file for batch processing.')
        self.actionSaveConfig.setShortcut('Ctrl+2')
        self.actionBatchProcess.setText('Batch Process')
        self.actionBatchProcess.setShortcut('Ctrl+3')


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    main_window = Ui_MainWindow()
    main_window.show()
    main_window.activateWindow()
    main_window.raise_()
    sys.exit(app.exec_())
    
    












