import os

import numpy as np
from aicsimageio.writers import OmeTiffWriter
from matplotlib import pyplot as plt
from skimage import exposure, morphology

# from src.stackio.Channel import channel
from src.segmentation.segment_GFP import get_stack_channel


def create_histogram(plotimg, title, savepath):
    fig, ax1 = plt.subplots()
    ax1.hist(plotimg.ravel(), bins=256, histtype='step', color='black')
    plt.title(title)
    ax1.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))
    plt.xlabel('Pixel intensity')
    plt.ylabel('No. of Pixels')
    ax2 = ax1.twinx()
    img_cdf, bins = exposure.cumulative_distribution(plotimg, 256)
    ax2.plot(bins, img_cdf, 'r')
    # ax2.set_yticks([])
    plt.ylabel("Fraction of total intensity")
    plt.savefig(f"{savepath}_{title}.png")
    plt.close()
    plt.clf()


####################################ACTB####################################

maindirpath = "C:/Users/satheps/PycharmProjects/Results/2022/Apr8/bkgcorr_tests/SLC_ex/"
savedirpath = "C:/Users/satheps/PycharmProjects/Results/2022/Apr8/bkgcorr_tests/"

assert os.path.exists(maindirpath)
assert os.path.exists(savedirpath)
# usechanneldirnames  = [channel.dirnames[c] for c in usechannels]
channelname = 'slc25a17'
savesegmentation = False
plot_histograms = True

for f in os.listdir(maindirpath):
    bn, ext1, ext2 = f.split(".")
    fpath = os.path.join(maindirpath, f).replace("\\", "/")
    IMG = get_stack_channel(fpath, channelname, verbose=False)
    IMG = IMG.astype(np.uint16)
    # 2.75 * 2
    # footprint = morphology.ball(27/2)
    footprint = morphology.cube(27)
    print(footprint.shape)
    print(IMG.shape)
    wres = morphology.white_tophat(IMG.copy(), footprint)
    bres = morphology.black_tophat(IMG, footprint)
    IMG_tophat = wres  # IMG - wres to remove small objects
    IMG_bottomhat = bres
    IMG_clahe = exposure.equalize_adapthist(IMG, clip_limit=0.03)
    if savesegmentation:
        OmeTiffWriter.save(data=IMG_tophat, uri=f'{savedirpath}{bn}_wtophat.tiff', compress=6)
        OmeTiffWriter.save(data=IMG_tophat, uri=f'{savedirpath}{bn}_btophat.tiff', compress=6)
        OmeTiffWriter.save(data=IMG_clahe, uri=f'{savedirpath}{bn}_clahe.tiff', compress=6)
        OmeTiffWriter.save(data=IMG, uri=f'{savedirpath}{f}', compress=6)
    if plot_histograms:
        create_histogram(IMG_tophat, title="wtophat", savepath=f'{savedirpath}{bn}')
        create_histogram(IMG_tophat, title="btophat", savepath=f'{savedirpath}{bn}')
        create_histogram(IMG_clahe, title="clahe", savepath=f'{savedirpath}{bn}')
        create_histogram(IMG, title="orig", savepath=f'{savedirpath}{bn}')
