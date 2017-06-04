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

hdu_list = fits.open('frame.fit')

image_data = hdu_list[0].data

print(type(image_data))
print(image_data.shape)
print(fits.getheader('frame.fit'))

hdu_list.close()

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
