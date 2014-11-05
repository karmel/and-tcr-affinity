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

    def rescale_intensity(self, image, bottom=2, top=98):
        '''
        Rescale the intensity at the ends of the spectral range
        to make blobbing better.
        '''
        low, high = np.percentile(image, (bottom, top))
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

    def plot_blobs(self, image, blobs, savepath=None, show_plot=False):
        _, ax = plt.subplots(1, 1)
        ax.imshow(image, interpolation='nearest')
        for blob in blobs:
            y, x, r = blob
            c = plt.Circle((x, y), r, color='black', linewidth=2, fill=False)
            ax.add_patch(c)

        if savepath:
            plt.savefig(savepath)

        if show_plot:
            plt.show()

        plt.close()

    def make_mask(self, image, blobs):
        '''
        Given a set of blobs, create a mask with circles appropriately
        placed, the same size as the image.

        At the same time, save the coordinates of the circles so that we can
        extract the squares with contained blobs if desired.

        Return both the mask and the set of squares as
        [row_start, row_end, col_start, col_end]
        '''

        mask_array = np.zeros(image.shape)
        squares = []
        for x, y, r in blobs:
            rows, cols = circle(x, y, r)
            # Make sure we keep our circles in bounds
            rows = [max(0, min(r, mask_array.shape[0] - 1)) for r in rows]
            cols = [max(0, min(c, mask_array.shape[1] - 1)) for c in cols]
            mask_array[rows, cols] = 1
            squares.append([min(rows), max(rows), min(cols), max(cols)])
        return mask_array, squares

    def mask_image(self, image, mask, savepath=None, show_plot=False):
        '''
        Zero out non-CD4 elements of the image.
        '''
        masked = np.minimum(image, mask)
        _, (ax1, ax2) = plt.subplots(1, 2)
        ax1.imshow(image)
        ax2.imshow(masked)

        if savepath:
            plt.savefig(savepath)

        if show_plot:
            plt.show()

        plt.close()
        return masked

    def extract_squares(self, image, squares, savepath=None, show_plot=False):
        '''
        Given an image and a set of [row_start, row_end, col_start, col_end]
        coordinates, extract the squares and return the set of image segments.

        Savepath should have a formattable place for the index of the segment
        if passed-- ie, `'extracted_cell_cd4_stain_{}.png'`
        '''
        segments = []
        for i, (row_start, row_end, col_start, col_end) in enumerate(squares):
            segment = image[row_start:row_end, col_start:col_end]
            segments.append(segment)
            _, ax1 = plt.subplots(1, 1)
            ax1.imshow(segment)
            if savepath:
                plt.savefig(savepath.format(i))

            if show_plot:
                plt.show()

            plt.close()

        return segments

    def find_low_intensity(self, segments, min_intensity=.5, threshold=.9):
        '''
        Returns the indices of segments where more than `threshold` fraction
        of the total pixels are below the `min_intensity` intensity value.
        '''
        low_intensity = []
        for i, segment in enumerate(segments):
            total_pixels = segment.size
            below_min = sum(segment.flatten() < min_intensity)
            if below_min / total_pixels > threshold:
                low_intensity.append(i)
        return low_intensity
