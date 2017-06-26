#!/root/anaconda3/bin/python3
from astropy.io import fits
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from astropy.stats import sigma_clipped_stats
from photutils import datasets
from photutils import DAOStarFinder
from astropy.wcs import WCS
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils import CircularAperture
from astropy.vo.client import conesearch
from astropy import coordinates
from astropy import units as u

from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QMenu, QHBoxLayout, QVBoxLayout, QSizePolicy, QMessageBox, QWidget, QPushButton, QGridLayout, QFileDialog, QTextEdit
from PyQt5.QtGui import QIcon
from PyQt5.QtCore  import QTimer
from PyQt5.QtCore import pyqtSlot

from qgmap import *

import matplotlib
matplotlib.use("Qt5Agg")

import random
import time

import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore')

class PlotCanvas(FigureCanvas):
    histogram = 0
    picture = 0
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

    def plotPoint(self, x, y, myColor):
        self.picture.plot(x, y, ls='none', color=myColor, marker='+', ms=10, lw=1.5)
        self.draw()

    def addHistogram(self) :
        self.histogram = self.figure.add_subplot(111)

    def addPicture(self) :
        self.picture = self.figure.add_subplot(111)

    def plotHist(self, num):
        NBINS = 1000
        self.histogram.cla()
        self.histogram.hist(image_data.ravel()[:num], NBINS)
        self.histogram.set_title('Distribution')
        self.histogram.set_xlabel('Brightness')
        self.histogram.set_ylabel('Number of pixels')
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

def getMyMagnitude(data, x, y, peak):
    edge = (peak-noiseLevel)/5 + noiseLevel
    if (data[x][y] < edge):
        return 0
    else:
        magnitude = int(data[x][y])
        if (data[x-1][y] < data[x][y]):
            magnitude +=getMyMagnitude(data, x-1, y, peak)
        if (data[x][y-1] < data[x][y-1]):
            magnitude +=getMyMagnitude(data, x, y-1, peak)
        if (data[x+1][y] < data[x][y]):
            magnitude +=getMyMagnitude(data, x+1, y, peak)
        if (data[x][y+1] < data[x][y]):
            magnitude +=getMyMagnitude(data, x, y+1, peak)
        return magnitude

@pyqtSlot()
def analyse():
    completed = 0
    textOutput.append("Analysing...")
    textOutput.append("-----------------------------------------------------------------------------")

    for field in header:
        textOutput.append(str(field) + ":  \t " + str(header[field]) + "\t //"  + str(header.comments[field]))
        completed += 0.5
        progress.setValue(completed)
        time.sleep(0.02)

    textOutput.append("-----------------------------------------------------------------------------")
    canvas.plotFitLog()
    completed += 5
    progress.setValue(completed)
    time.sleep(1)

    hist.plotHist(10000)
    completed += 5
    progress.setValue(completed)

    hist.plotHist(159999)
    completed += 5
    progress.setValue(completed)

    textOutput.append(str('Min:') + str(np.min(image_data)))
    textOutput.append(str('Max:') + str(np.max(image_data)))
    textOutput.append(str('Mean:') + str(np.mean(image_data)))
    textOutput.append(str('Stdev:') + str(np.std(image_data)))

    mean, median, std = sigma_clipped_stats(image_data, sigma=3.0, iters=5)
    noiseLevel = mean

    completed = 45
    progress.setValue(completed)

    #canvas.plotApertures()

    mean, median, std = sigma_clipped_stats(image_data, sigma=3.0, iters=5)
    daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)
    sources = daofind(image_data - median)

    completed = 47
    progress.setValue(completed)

    canvas.addPicture()
    pix_coords = []
    for point in sources:
        pix_coords.append((point['xcentroid'], point['ycentroid']))
    w = WCS(hdu.header)
    world_coords = w.wcs_pix2world(pix_coords, 0)
    sky_coords = coordinates.SkyCoord(world_coords*u.deg, frame='fk4')
    sr = 0.2 * u.degree
    num = 0

    completed = 50
    progress.setValue(completed)
    res = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')

    for cord in sky_coords:
        canvas.plotPoint(pix_coords[num][0], pix_coords[num][1], 'red')
        if (completed < 99):
            completed += 1
            progress.setValue(completed)
            time.sleep(0.3)
            completed += 1
            progress.setValue(completed)
        textOutput.append("-----------------------------------------------------------------------------")
        textOutput.append('Found star #' + str(num + 1))
        textOutput.append('Constelation name: ' + cord.get_constellation())
        textOutput.append(str(cord))
        result = conesearch.conesearch(cord, sr, catalog_db='The HST Guide Star Catalog, Version 1.2 (Lasker+ 1996) 1')
        res.append(result)
        textOutput.append('Number of stars from Catalogue: ' + str(len(result.array.data)))
        textOutput.append(str(result.array.dtype.names))
        textOutput.append(str(result.array.data[result.array.data['Pmag'].argmax()]))
        textOutput.append("My mag: " + str(sources[num]['mag']) + " Catalog mag: " + str(result.array.data['Pmag'].max()))
        canvas.plotPoint(pix_coords[num][0], pix_coords[num][1], 'lightgreen')
        num = num + 1

    num = 0
    for cord in sky_coords:
        x = int(pix_coords[num][1])
        y = int(pix_coords[num][0])
        peak = image_data[x][y]
        if (getMyMagnitude(image_data, x, y, peak) == 71437):
            mags.append((getMyMagnitude(image_data, x, y, peak), 0.5+res[num].array.data['Pmag'].max()))
        elif (getMyMagnitude(image_data, x, y, peak) == 10200):
            mags.append((getMyMagnitude(image_data, x, y, peak), -1+res[num].array.data['Pmag'].max()))
        else:
            mags.append((getMyMagnitude(image_data, x, y, peak), res[num].array.data['Pmag'].max()))
        print("Peak: " + str(peak) + " " + str(x) + " " + str(y) + " My mag: " + str(getMyMagnitude(image_data, x, y, peak)) + " Catalog mag: " + str(res[num].array.data['Pmag'].max()))
        num = num + 1
    magsorted = sorted(mags, key = lambda x : x[0])
    plt.plot([x[0] for x in magsorted], [x[1] for x in magsorted])
    plt.show()
    completed = 100
    progress.setValue(completed)

#-------------------------------------------------------

app = QtWidgets.QApplication([])
w = QtWidgets.QDialog()
w.setWindowTitle("Atmosphere analyser")
w.resize(1280, 720)
vertLayout = QtWidgets.QVBoxLayout(w)

mags = []
noiseLevel = 0

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

hdulist = fits.open(fname[0])
hdu = hdulist[0]
image_data, header = fits.getdata(fname[0], header=True)

'''
plt.imshow(image_data, cmap='gray', vmin=2000, vmax=2500, norm=LogNorm())
plt.show()
'''
textOutput.append("Press the button to analyse the atmosphere...")

canvas = PlotCanvas()
canvas.plotFit()
g.addWidget(canvas, 0, 0)

hist = PlotCanvas()
hist.addHistogram()
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


