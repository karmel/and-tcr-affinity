'''
Created on Nov 4, 2014

@author: karmel

Ten minute timepoint images.
'''

import os

from andtcr.microscopy.base import ImageAnalyzer
import numpy as np


def process_sample_set(yzer, timepoint, sample, antibody, image_num):
    pathing = ['/Users/karmel/GlassLab/Notes_and_Reports',
               'AND_TCR', 'Microscopy',
               'for_processing', timepoint,
               sample, antibody]
    filename = '{}_{}_{}_{}_{{}}.tif'.format(
        timepoint, sample, antibody, image_num)
    orig_image = yzer.import_image(*pathing +
                                   [filename.format('cd4')])

    save_dir = os.path.join(*pathing + ['output_{}'.format(image_num)])
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    # Prep CD4 image
    image = yzer.rescale_intensity(orig_image)
    image = yzer.grayscale(image)
    blobs = yzer.get_blobs(
        image, min_sigma=5, max_sigma=25, num_sigma=100, threshold=.01)
    blobs = yzer.filter_blobs(blobs, 14, 22)
    savepath = os.path.join(save_dir, 'blob_localization.png')
    yzer.plot_blobs(image, blobs, savepath=savepath)
    mask, squares = yzer.make_mask(image, blobs)

    # Isolate CD4s
    cd4_image = yzer.grayscale(orig_image)
    savepath = os.path.join(save_dir, 'masked_image_cd4.png')
    masked_cd4 = yzer.mask_image(cd4_image, mask, savepath=savepath)

    # Get second antibody
    second_image = yzer.import_image(*pathing +
                                     [filename.format(antibody)])
    second_image = yzer.grayscale(second_image)
    savepath = os.path.join(save_dir, 'masked_image_{}.png'.format(antibody))
    masked_second = yzer.mask_image(second_image, mask, savepath=savepath)

    # Pick out CD4s from the masked image.
    savepath = os.path.join(save_dir, 'segment_cd4_stain_{}.png')
    cells_cd4 = yzer.extract_squares(masked_cd4, squares, savepath=savepath)
    savepath = os.path.join(
        save_dir, 'segment_{}_cd4_stain_{{}}.png'.format(antibody))
    cells_second = yzer.extract_squares(
        masked_second, squares, savepath=savepath)

    # Now filter out cells without a sufficient amount of cd4--
    # These are likely DCs, etc, that have taken up some antibody.
    cells_cd4 = [yzer.rescale_intensity(cell) for cell in cells_cd4]
    cells_second = [yzer.rescale_intensity(cell) for cell in cells_second]

    to_remove = yzer.find_low_intensity(cells_cd4)
    print('Removing cells with indices', to_remove)
    cells_cd4 = [el for i, el in enumerate(cells_cd4) if i not in to_remove]
    cells_second = [
        el for i, el in enumerate(cells_second) if i not in to_remove]

    # Determine similarity for all cells.
    cd4_scores = yzer.determine_scores(cells_cd4)
    print('CD4 scores\t', cd4_scores)
    print('CD4 mean, std\t{}\t{}'.format(
        np.mean(cd4_scores), np.std(cd4_scores)))

    second_scores = yzer.determine_scores(cells_second)
    print(antibody + ' scores\t', second_scores)
    return second_scores

if __name__ == '__main__':

    yzer = ImageAnalyzer()

    yzer.skip_plotting = False

    timepoints = ['10min']
    samples = ['no_peptide', '10um_k99a', '100um_k99a',
               '0_001um_pcc', '10um_pcc']
    antibodies = ['perk', 'plat', 'pzap70', 'rasgrp', 'sos']

    for timepoint in timepoints:
        for sample in samples:
            for antibody in antibodies:
                scores = []
                for image_num in range(1, 5):
                    scores += process_sample_set(yzer,
                                                 timepoint,
                                                 sample,
                                                 antibody,
                                                 image_num)

                above = [s for s in scores if s > 2]
                print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                    timepoint,
                    sample,
                    antibody,
                    np.mean(scores),
                    np.std(scores),
                    len(above),
                    len(scores),
                    len(above) / len(scores)))
