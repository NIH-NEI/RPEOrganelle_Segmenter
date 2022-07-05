import os

from aicsimageio import AICSImage
from matplotlib import pyplot as plt

from src.stackio import experimentalparams

# from src.stackio import stack3d

path_all_stacks = "C:/Users/satheps/PycharmProjects/Results/2022/final_stacks"
path_all_segmentations = "C:/Users/satheps/PycharmProjects/Results/2022/final_segmentations"
dirlist_stacks = os.listdir(path_all_stacks)
dirlist_segs = os.listdir(path_all_segmentations)
savepath = "C:/Users/satheps/PycharmProjects/Results/2022/Jun24/illustration_imgs/"
assert dirlist_segs == dirlist_stacks
assert os.path.exists(savepath)
roi_r = 3
roi_f = 2
roi_w = "W3"
roi_z = 9  # middle slice =13
x, y = 50, 50
roi_s = 750


# usevars = {}
def getwr(filename, gfp):
    """
    :param filename: gfp_channel_filenames
    :return: week id, replicate id, week no. replicate number, common base string
    """
    filename = filename.split(".")[0]
    if gfp:
        filename = filename.replace("s2", "")  # NOTE: temporary for lamp1
        basestring = "_".join(filename.split("_")[:3])
    else:
        basestring = "_".join(filename.split("_")[:3])
    s1, r, fov = basestring.split("_")
    # print(s1, r, fov, basestring)
    w = s1.split("-")[1]
    w_ = experimentalparams.WS.index(w)
    r_ = int(r[1:]) - 2  # r goes from 2 to 11 - change it to 0-9
    fov_ = int(fov[-1:]) - 1  # fov goes from 1 to 6 - change it to 0-5
    return w, r, w_, r_, fov, fov_, basestring


for i, (org, seg) in enumerate(zip(dirlist_stacks, dirlist_segs)):
    assert org == seg
    IMG_org = None
    IMG_seg = None
    orgstring = None
    segstring = None
    wlist_org = os.listdir(os.path.join(path_all_stacks, org))

    for wstacks in wlist_org:
        if roi_w in wstacks:
            flist = [f for f in os.listdir(os.path.join(path_all_stacks, org, wstacks)) if f.__contains__(".tif")]
            # print(flist)
            for f in flist:
                w, r, w_, r_, fov, fov_, basestring = getwr(f, gfp=False)
                if roi_w == w and roi_r == r_ and roi_f == fov_:
                    IMG_org = AICSImage(os.path.join(path_all_stacks, org, wstacks, f)).data
                    print(w, r, w_, r_, fov, fov_, basestring)
                    orgstring = basestring

    wlist_seg = os.listdir(os.path.join(path_all_segmentations, seg))
    flist_seg = [f for f in wlist_seg if f.__contains__(".tif")]
    for fs in flist_seg:
        # print(fs)
        w, r, w_, r_, fov, fov_, basestring = getwr(fs, gfp=True)
        if roi_w == w and roi_r == r_ and roi_f == fov_:
            IMG_seg = AICSImage(os.path.join(path_all_segmentations, seg, fs)).data
            print(w, r, w_, r_, fov, fov_, basestring)
            segstring = basestring
    assert orgstring == segstring
    roi_org = IMG_org[0, 0, roi_z, x:x + roi_s, y:y + roi_s]
    roi_seg = IMG_seg[0, 0, roi_z, x:x + roi_s, y:y + roi_s] * 255
    fig, axarr = plt.subplots(1, 2)
    ch = orgstring.split("_")[0].split("-")[-1]
    fig.suptitle(ch, y=0.8)
    axarr[0].imshow(roi_org, cmap='gray')
    axarr[0].axis('off')
    axarr[1].imshow(roi_seg, cmap='gray')
    axarr[1].axis('off')
    # plt.show()
    plt.savefig(savepath + f"{orgstring}", bbox_inches='tight')
    plt.close()
    plt.clf()
