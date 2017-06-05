#!/usr/bin/python3
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from astropy.stats import sigma_clipped_stats
from photutils import datasets
from photutils import DAOStarFinder
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils import CircularAperture

from qgmap import *

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

app = QtGui.QApplication([])
w = QtGui.QDialog()
h = QtGui.QVBoxLayout(w)
l = QtGui.QFormLayout()
h.addLayout(l)

gmap = QGoogleMap(w)
gmap.mapMoved.connect(onMapMoved)
gmap.markerMoved.connect(onMarkerMoved)
gmap.mapClicked.connect(onMapLClick)
gmap.mapDoubleClicked.connect(onMapDClick)
gmap.mapRightClicked.connect(onMapRClick)
gmap.markerClicked.connect(onMarkerLClick)
gmap.markerDoubleClicked.connect(onMarkerDClick)
gmap.markerRightClicked.connect(onMarkerRClick)
h.addWidget(gmap)
gmap.setSizePolicy(
        QtGui.QSizePolicy.MinimumExpanding,
        QtGui.QSizePolicy.MinimumExpanding)
w.show()

gmap.waitUntilReady()

image_data, header = fits.getdata("frame.fit", header=True)

print(type(image_data))
print(image_data.shape)
print("-----------------------------------------------------------------------------")
for field in header:
    print(str(field) + ":  \t " + str(header[field]) + "\t //"  + str(header.comments[field]))

print("-----------------------------------------------------------------------------")

coords = gmap.centerAt(header["SITELAT"], header["SITELON"])
gmap.setZoom(8)

gmap.addMarker("MyDragableMark", header["SITELAT"], header["SITELON"])

app.exec_()
exit()

plt.imshow(image_data, cmap='gray')
plt.colorbar()
plt.show()

print('Min:', np.min(image_data))
print('Max:', np.max(image_data))
print('Mean:', np.mean(image_data))
print('Stdev:', np.std(image_data))

NBINS = 1000
histogram = plt.hist(image_data.flat, NBINS)
plt.show()

plt.imshow(image_data, cmap='gray', vmin=2000, vmax=2500, norm=LogNorm())
plt.colorbar()
plt.show();

mean, median, std = sigma_clipped_stats(image_data, sigma=3.0, iters=5)    
print((mean, median, std))  

daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)    
sources = daofind(image_data - median)    
print(sources) 

positions = (sources['xcentroid'], sources['ycentroid'])
apertures = CircularAperture(positions, r=10.)
#norm = ImageNormalize(stretch=SqrtStretch())
plt.imshow(image_data, cmap='gray', vmin=2000, vmax=2500, norm=LogNorm())
apertures.plot(color='white', lw=1, alpha=1)
plt.show()
