'''
Created on Nov 4, 2014

@author: karmel

We have confocal microscopy images of cells stained with two colors--
CD4 and another antibody of interest. The primary goal is to determine
the membrane localization of the second molecule of interest in the CD4
cells.

So, the first step is to identify CD4 cells and throw out the rest. The
dendritic cells seem to take up some antibody, so we will want to find
blobs in the CD4 channel and then filter out the small or large ones.

Next we will use the identified CD4s as a mask, and compare spectral overlap
in the second channel with the CD4 stain, as the CD4 stain is 
membrane-localized.
'''
from math import sqrt
import os

from matplotlib import pyplot as plt
from skimage import io, exposure
from skimage.color import rgb2gray
from skimage.draw.draw import circle
from skimage.feature import blob_log

import numpy as np


class ImageAnalyzer(object):

    def import_image(self, *fileparts):
        filename = os.path.join(*fileparts)
        return io.imread(filename)

    def rescale_intensity(self, image):
        '''
        Rescale the intensity at the ends of the spectral range
        to make blobbing better.
        '''
        low, high = np.percentile(image, (2, 98))
        image = exposure.rescale_intensity(image, in_range=(low, high))

        return image

    def increase_contrast(self, image):
        '''
        Expand spectrum at frequently used intensities.
        '''
        return exposure.equalize_hist(image)

    def grayscale(self, image):
        '''
        Render as grayscale
        '''
        return rgb2gray(image)

    def get_blobs(self, image, **kwargs):
        '''
        Use Laplacian of Gaussian to find blobs in passed grayscale image.
        '''
        blobs = blob_log(image, **kwargs)
        # Compute radii in the 3rd column.
        # Expand the radius to get a circle around the whole cell.
        blobs[:, 2] = blobs[:, 2] * 1.3 * sqrt(2)
        return blobs

    def filter_blobs(self, blobs, min, max):
        '''
        Filter blobs for radius range.
        '''
        blobs = blobs[blobs[:, 2] >= min]
        blobs = blobs[blobs[:, 2] <= max]
        return blobs

    def plot_blobs(self, image, blobs, savepath=None):
        _, ax = plt.subplots(1, 1)
        ax.imshow(image, interpolation='nearest')
        for blob in blobs:
            y, x, r = blob
            c = plt.Circle((x, y), r, color='black', linewidth=2, fill=False)
            ax.add_patch(c)

        if savepath:
            plt.savefig(savepath)

        plt.show()

    def make_mask(self, image, blobs):
        '''
        Given a set of blobs, create a mask with circles appropriately
        placed, the same size as the image.
        '''

        mask_array = np.zeros(image.shape)
        for x, y, r in blobs:
            rows, cols = circle(x, y, r)
            # Make sure we keep our circles in bounds
            rows = [max(0, min(r, mask_array.shape[0] - 1)) for r in rows]
            cols = [max(0, min(c, mask_array.shape[1] - 1)) for c in cols]
            mask_array[rows, cols] = 1

        return mask_array

    def mask_image(self, image, mask):
        '''
        Zero out non-CD4 elements of the image.
        '''
        masked = np.minimum(image, mask)
        _, (ax1, ax2) = plt.subplots(2, 1)
        ax1.imshow(image)
        ax2.imshow(masked)
        plt.show()
