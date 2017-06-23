from astropy.io import fits
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from astropy.stats import sigma_clipped_stats
from photutils import datasets
from photutils import DAOStarFinder
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils import CircularAperture

from PyQt5.QtWidgets import QApplication, QMainWindow, QMenu, QHBoxLayout, QVBoxLayout, QSizePolicy, QMessageBox, QWidget, QPushButton, QGridLayout, QFileDialog, QTextEdit
from PyQt5.QtGui import QIcon
from PyQt5.QtCore  import QTimer
from PyQt5.QtCore import pyqtSlot

from qgmap import *

import matplotlib
matplotlib.use("Qt5Agg")

import random
import time

class PlotCanvas(FigureCanvas):
    ax = 0
    def __init__(self, parent=None, width=7, height=7, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)

        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                QSizePolicy.Expanding,
                QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def plotFit(self):
        ax = self.figure.add_subplot(111)
        ax.imshow(image_data, cmap='gray')
        ax.set_title('Image')
        self.draw()

    def plotFitLog(self):
        ax = self.figure.add_subplot(111)
        ax.imshow(image_data, cmap='gray', vmin=2000, vmax=2500, norm=LogNorm())
        ax.set_title('Image')
        self.draw()

    def addSubplot(self) :
        self.ax = self.figure.add_subplot(111)

    def plotHist(self, num):
        NBINS = 1000
        self.ax.cla()
        self.ax.hist(image_data.ravel()[:num], NBINS)
        self.ax.set_title('Distribution')
        self.ax.set_xlabel('Brightness')
        self.ax.set_ylabel('Number of pixels')
        self.draw()

def goCoords() :
    def resetError() :
            coordsEdit.setStyleSheet('')
    try : latitude, longitude = coordsEdit.text().split(",")
    except ValueError :
            coordsEdit.setStyleSheet("color: red;")
            QtCore.QTimer.singleShot(500, resetError)
    else :
            gmap.centerAt(latitude, longitude)
            gmap.moveMarker("MyDragableMark", latitude, longitude)

def goAddress() :
    def resetError() :
            addressEdit.setStyleSheet('')
    coords = gmap.centerAtAddress(addressEdit.text())
    if coords is None :
            addressEdit.setStyleSheet("color: red;")
            QtCore.QTimer.singleShot(500, resetError)
            return
    gmap.moveMarker("MyDragableMark", *coords)
    coordsEdit.setText("{}, {}".format(*coords))

def onMarkerMoved(key, latitude, longitude) :
    print("Moved!!", key, latitude, longitude)
    coordsEdit.setText("{}, {}".format(latitude, longitude))
def onMarkerRClick(key) :
    print("RClick on ", key)
    gmap.setMarkerOptions(key, draggable=False)
def onMarkerLClick(key) :
    print("LClick on ", key)
def onMarkerDClick(key) :
    print("DClick on ", key)
    gmap.setMarkerOptions(key, draggable=True)

def onMapMoved(latitude, longitude) :
    print("Moved to ", latitude, longitude)
def onMapRClick(latitude, longitude) :
    print("RClick on ", latitude, longitude)
def onMapLClick(latitude, longitude) :
    print("LClick on ", latitude, longitude)
def onMapDClick(latitude, longitude) :
    print("DClick on ", latitude, longitude)

@pyqtSlot()
def analyse():
    completed = 0
    textOutput.append("Analysing...")
    while (completed < 5):
        completed += 0.1
        progress.setValue(completed)
        time.sleep(0.02)

    canvas.plotFitLog()
    while (completed < 10):
        completed += 0.1
        progress.setValue(completed)
        time.sleep(0.02)

    for index, pixel in enumerate(image_data.ravel()) :
        if (index%10000 == 0) :
            textOutput.append("Plotting Hist " + str(index))
            hist.plotHist(index)
            completed += 1
            progress.setValue(completed)

    textOutput.append(str('Min:') + str(np.min(image_data)))
    textOutput.append(str('Max:') + str(np.max(image_data)))
    textOutput.append(str('Mean:') + str(np.mean(image_data)))
    textOutput.append(str('Stdev:') + str(np.std(image_data)))
'''
    mean, median, std = sigma_clipped_stats(image_data, sigma=3.0, iters=5)
    textOutput.append(str(mean, median, std))

    daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)
    sources = daofind(image_data - median)
    textOutput.append(str(sources))

    positions = (sources['xcentroid'], sources['ycentroid'])
    apertures = CircularAperture(positions, r=10.)
    #norm = ImageNormalize(stretch=SqrtStretch())
    plt.imshow(image_data, cmap='gray', vmin=2000, vmax=2500, norm=LogNorm())
    apertures.plot(color='white', lw=1, alpha=1)
    plt.show()
'''
#-------------------------------------------------------

app = QtWidgets.QApplication([])
w = QtWidgets.QDialog()
w.setWindowTitle("Atmosphere analyser")
w.resize(1280, 720)
vertLayout = QtWidgets.QVBoxLayout(w)

g = QtWidgets.QGridLayout()
vertLayout.addLayout(g)

h = QtWidgets.QHBoxLayout()
vertLayout.addLayout(h)

fname = QFileDialog.getOpenFileName(None, 'Open file', '',"Fit files (*.fit *.fits)")

w.show()

button = QPushButton('Start analysis')
button.setToolTip('Press the button to start analysis')
button.clicked.connect(analyse)
h.addWidget(button)

progress = QtWidgets.QProgressBar()
h.addWidget(progress)

textOutput = QtWidgets.QTextEdit()
textOutput.setReadOnly(True)
g.addWidget(textOutput, 1, 1)

image_data, header = fits.getdata(fname[0], header=True)

textOutput.append(str(type(image_data)))
textOutput.append(str(image_data.shape))
textOutput.append("-----------------------------------------------------------------------------")
for field in header:
    textOutput.append(str(field) + ":  \t " + str(header[field]) + "\t //"  + str(header.comments[field]))

textOutput.append("-----------------------------------------------------------------------------")

canvas = PlotCanvas()
canvas.plotFit()
g.addWidget(canvas, 0, 0)

hist = PlotCanvas()
hist.addSubplot()
g.addWidget(hist, 1, 0)

gmap = QGoogleMap(w)
gmap.mapMoved.connect(onMapMoved)
gmap.markerMoved.connect(onMarkerMoved)
gmap.mapClicked.connect(onMapLClick)
gmap.mapDoubleClicked.connect(onMapDClick)
gmap.mapRightClicked.connect(onMapRClick)
gmap.markerClicked.connect(onMarkerLClick)
gmap.markerDoubleClicked.connect(onMarkerDClick)
gmap.markerRightClicked.connect(onMarkerRClick)
g.addWidget(gmap, 0, 1)
gmap.setSizePolicy(
        QtWidgets.QSizePolicy.Expanding,
        QtWidgets.QSizePolicy.Expanding)
gmap.waitUntilReady()
gmap.addMarker("Observatory", header["SITELAT"], header["SITELON"])
coords = gmap.centerAt(header["SITELAT"], header["SITELON"])
gmap.setZoom(8)

app.exec_()
exit()


