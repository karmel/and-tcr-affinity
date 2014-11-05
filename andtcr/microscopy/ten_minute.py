'''
Created on Nov 4, 2014

@author: karmel

Ten minute timepoint images.
'''

from andtcr.microscopy.base import ImageAnalyzer
if __name__ == '__main__':

    yzer = ImageAnalyzer()

    antibody = 'pzap70'
    pathing = ['/Users/karmel/GlassLab/Notes_and_Reports',
               'AND_TCR', 'Microscopy',
               '2014-9-24 2 & 10 min timepoint', '10 min',
               'No peptide', 'pZap70']
    filename = '2014-9-25 AND CD4+ 10 min_NP_pZap70_1_{}.tif'
    orig_image = yzer.import_image(*pathing +
                                   [filename.format('cd4')])

    # Prep CD4 image
    image = yzer.rescale_intensity(orig_image)
    image = yzer.grayscale(image)
    blobs = yzer.get_blobs(
        image, min_sigma=5, max_sigma=25, num_sigma=100, threshold=.01)
    blobs = yzer.filter_blobs(blobs, 14, 22)
    yzer.plot_blobs(image, blobs)
    mask = yzer.make_mask(image, blobs)

    # Isolate CD4s
    cd4_image = yzer.grayscale(yzer.increase_contrast(orig_image))
    masked_cd4 = yzer.mask_image(cd4_image, mask)

    # Get second antibody
    second_image = yzer.import_image(*pathing +
                                     [filename.format(antibody)])
    second_image = yzer.grayscale(yzer.increase_contrast(second_image))
    masked_second = yzer.mask_image(second_image, mask)

    # Pick out CD4s from
