'''
Created on Nov 4, 2014

@author: karmel

Ten minute timepoint images.
'''

import os

from andtcr.microscopy.base import ImageAnalyzer
if __name__ == '__main__':

    yzer = ImageAnalyzer()

    antibody = 'pzap70'
    pathing = ['/Users/karmel/GlassLab/Notes_and_Reports',
               'AND_TCR', 'Microscopy',
               '2014-9-24 2 & 10 min timepoint', '10 min',
               'No peptide', 'pZap70']
    save_dir = os.path.join(*pathing + ['output'])
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
    mask, squares = yzer.make_mask(image, blobs)

    # Isolate CD4s
    cd4_image = yzer.grayscale(orig_image)
    masked_cd4 = yzer.mask_image(cd4_image, mask)

    # Get second antibody
    second_image = yzer.import_image(*pathing +
                                     [filename.format(antibody)])
    second_image = yzer.grayscale(second_image)
    masked_second = yzer.mask_image(second_image, mask)

    # Pick out CD4s from the masked image.
    savepath = os.path.join(save_dir, 'segment_cd4_stain_{}.png')
    cells_cd4 = yzer.extract_squares(masked_cd4, squares, savepath=savepath)
    savepath = os.path.join(
        save_dir, 'segment_' + antibody + '_cd4_stain_{}.png')
    cells_second = yzer.extract_squares(
        masked_second, squares, savepath=savepath)

    # Now filter out cells without a sufficient amount of cd4--
    # These are likely DCs, etc, that have taken up some antibody.
    cells_cd4 = [yzer.rescale_intensity(cell) for cell in cells_cd4]
    cells_second = [yzer.rescale_intensity(cell) for cell in cells_second]

    to_remove = yzer.find_low_intensity(cells_cd4)
    print(len(cells_cd4))
    [cells_cd4.remove(i) for i in to_remove]
    [cells_second.remove(i) for i in to_remove]
    print(len(cells_cd4))
    print(len(cells_second))

    # Pair the cells
    cd4_second_pairs = zip(cells_cd4, cells_second)
