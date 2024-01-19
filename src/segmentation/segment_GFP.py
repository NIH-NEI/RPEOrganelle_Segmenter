import os
import sys

import numpy as np
import tifffile
from aicsimageio import AICSImage
from aicsimageio.writers import OmeTiffWriter
from aicssegmentation.core.MO_threshold import MO
from aicssegmentation.core.pre_processing_utils import intensity_normalization, image_smoothing_gaussian_slice_by_slice, \
    image_smoothing_gaussian_3d, edge_preserving_smoothing_3d
from aicssegmentation.core.seg_dot import dot_2d_slice_by_slice_wrapper, dot_3d_wrapper
from aicssegmentation.core.utils import get_middle_frame, hole_filling, get_3dseed_from_mid_frame
from aicssegmentation.core.utils import topology_preserving_thinning
from aicssegmentation.core.vessel import filament_2d_wrapper, filament_3d_wrapper
from scipy.ndimage import distance_transform_edt
from scipy.stats import norm
from skimage.feature import peak_local_max
from skimage.measure import label
from skimage.morphology import remove_small_objects, dilation, ball, closing, disk
from skimage.segmentation import watershed


def get_stack_channel(fpath, channelname, verbose=False):
    if verbose:
        print(f"getting {channelname} stack")
    reader = AICSImage(fpath)
    print(reader.data.shape)
    IMG = reader.data.astype(np.float32)
    # channelobject = Channel.channel(channelname)
    return IMG[0, 0, :, :, :].copy()


def parseparam(params, string):
    """
    parse input parameter from list or dict

    params: parameter obtained from input
    string:  expected name of parameter
    :return: parameter value
    """
    paramvalue = None
    try:
        if isinstance(params, dict):
            paramvalue = params[string]
        if isinstance(params, list):
            paramvalue = params
        print(paramvalue)
        return paramvalue
    except Exception as e:
        print(e)
        return None


def savesegmented(data, savepath, basename, parametertype=None, postprocessinginfo=None, useparameterinfo=True):
    """
    Save segmented tif file
    
    Args:
        data:
        savepath:
        basename:
        parametertype:
        postprocessinginfo: e.g. 's3p'
    :return:
    """
    filename = basename
    if useparameterinfo:
        filename = f'{savepath}{basename}{parametertype}{postprocessinginfo}_GFP.tiff'
    OmeTiffWriter.save(data=data, uri=filename, compress=6)


def autoselect_normalization_parameters(imgstack, debug=False, percentile=99.99):
    """
    return suggested scaling parameter assuming the image is a representative example of this cell structure.
    Returns the values of selected parameters

    Original code at:
    long-url: https://github.com/AllenCell/aics-segmentation
    Args:
        imgstack: input image stack
        debug: toggle to display additional details for debugging or verbosity
        percentile: percentile of data to include in the calculation
    Returns:
        Suggested low and high normalization parameters

    """

    m, s = norm.fit(imgstack.flat)
    if debug:
        print(f" stack mean intensity: {m}\t\tstandard deviation: {s}")

    perc_intensity = np.percentile(imgstack, percentile)
    pmin = imgstack.min()
    pmax = imgstack.max()
    if debug:
        print(f"{percentile} percentile of the stack intensity is: {perc_intensity}")
    up_ratio = 0
    for up_i in np.arange(0.5, 1000, 0.5):
        if m + s * up_i > perc_intensity:
            if m + s * up_i > pmax:
                if debug:
                    print(f"suggested upper range is {up_i - 0.5}, which is {m + s * (up_i - 0.5)}")
                up_ratio = up_i - 0.5
            else:
                if debug:
                    print(f"suggested upper range is {up_i}, which is {m + s * up_i}")
                up_ratio = up_i
            break

    low_ratio = 0
    for low_i in np.arange(0.5, 1000, 0.5):
        if m - s * low_i < pmin:
            if debug:
                print(f"suggested lower range is {low_i - 0.5}, which is {m - s * (low_i - 0.5)}")
            low_ratio = low_i - 0.5
            break
    if debug:
        print(
            f"Suggested parameter for normalization is [{low_ratio}, {up_ratio}]\n To further enhance the contrast:"
            f" You may increase the first value  (may loss some dim parts), or decrease the second value (may loss some"
            f" texture in super bright regions)\n To slightly reduce the contrast: You may decrease the first value, or"
            f" increase the second value")
    return [low_ratio, up_ratio]


def segmentsec61tacks(fpath, savepath, params, channel="sec61b", minarea=None):
    """
       Segments sec61b structures in 3D image stacks.

       Args:
           fpath (str): Path to the input image stack.
           savepath (str): Path to save the segmented output.
           params (dict or list):
               - If dict: Dictionary containing segmentation parameters:
                   - "f2params": List of Frangi filter parameters (list of lists).
               - If list: List of Frangi vesselness filter parameters directly.
           channel (str, optional): Name of the channel to segment. Defaults to "sec61b".
           minarea (int, optional): Minimum object area for filtering. Defaults to None.

       Returns:
           None
       """
    struct_img0 = get_stack_channel(fpath, channel)

    intensity_scaling_param = autoselect_normalization_parameters(struct_img0)
    # [0.5, 5]  # [1, 2]  # [1.0, 9.0]
    ################################ code modified to return suggested params

    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)

    structure_img_smooth = edge_preserving_smoothing_3d(struct_img)
    ################################
    # scale = params["scale"]             # [0.75, 1]  # [0.75]
    # cutoff = params["cutoff"]            #[0.15]  # [0.15]
    # f2params = [[[scale, cutoff]] for scale in scales for cutoff in cutoffs]]
    print("params: ", params)
    try:
        if isinstance(params, dict):
            f2_param = params["f2params"]
        if isinstance(params, list):
            f2_param = params
    except Exception as e:
        print(e)
    # [10,15,20]
    bw = filament_2d_wrapper(structure_img_smooth.copy(), f2_param)

    seg = remove_small_objects(bw > 0, min_size=minarea, connectivity=1, in_place=False)
    seg = seg > 0
    out = seg.astype(np.uint8)
    out[out > 0] = 255
    #         s = "sec61"
    s = ""
    for i in range(len(f2_param)):
        for j in range(len(f2_param[i])):
            s = "_".join([s, str(f2_param[i][j])])
    name = fpath.split("/")[-1].split(".")[0]
    tifffile.imwrite(f'{savepath}{name}{s}_GFP.tiff', out, compress=6)


def segmentlaminstacks(fpath, savepath, params, channel="lmnb1", minarea=None, padsize=0):
    """
    Segments lamin structures in 3D image stacks.

    Args:
        fpath (str): Path to the input image stack.
        savepath (str): Path to save the segmented output.
        params (dict): Dictionary containing segmentation parameters:
            - "f2params": List of Frangi filter parameters (list of lists).
            - "useclosing": Boolean indicating whether to apply morphological closing (bool).
        channel (str, optional): Name of the channel to segment. Defaults to "lmnb1".
        minarea (int, optional): Minimum object area for filtering. Defaults to None.
        padsize (int, optional): Amount of padding to add to the image. Defaults to 0.

    Returns:
        None

    """
    struct_img0 = get_stack_channel(fpath, channel)
    intensity_scaling_param = autoselect_normalization_parameters(struct_img0)
    gscale = (0.5 / 0.216666)
    gaussian_smoothing_sigma = [gscale, 1, 1]
    struct_img = intensity_normalization(struct_img0.copy(), scaling_param=intensity_scaling_param)
    structure_img = struct_img.copy()
    ####################PAD####################
    if padsize:
        structure_img = np.pad(struct_img, [(padsize, padsize), (0, 0), (0, 0)])
    ###########################################
    # smoothing with gaussian filter
    structure_img_smooth = image_smoothing_gaussian_3d(structure_img.copy(),
                                                       sigma=gaussian_smoothing_sigma)

    ################################
    bw = None
    methodinfo = ""
    # f2_param = [[0.1, 0.01], [0.2, 0.01], [0.4, 0.01], [0.8, 0.01], [1.6, 0.01]]
    #     f2_param = [[0.1, 0.01], [0.25, 0.01], [0.5, 0.01], [1, 0.01], [1.5, 0.01], [2, 0.2], [3,0.5]]  # suggested for lamin [[0.5, 0.01]]
    f2_param = params["f2params"]
    closingstate = params["useclosing"]
    closingused = ""
    for f2pi in f2_param:
        methodinfo = methodinfo + f"_{f2pi[0]}_{f2pi[1]}"
    print(f"f2 params : {methodinfo}")
    #### middle frame analysis
    mid_z = get_middle_frame(structure_img_smooth, method='intensity')
    bw_mid_z = filament_2d_wrapper(structure_img_smooth[mid_z, :, :], f2_param)
    bw_dil_z = bw_mid_z.copy()
    if closingstate:
        bw_dil_z = closing(bw_mid_z.copy(), selem=disk(2))
        closingused = "_wclosing"
    ####################
    hole_max = 4000  # 40000
    hole_min = 1  # 400
    area_min = hole_min
    bw_fill_mid_z = hole_filling(bw_mid_z, hole_min, hole_max)
    seed = get_3dseed_from_mid_frame(np.logical_xor(bw_fill_mid_z, bw_dil_z), structure_img.shape,
                                     mid_z,
                                     hole_min)
    bw_filled = watershed(structure_img, seed.astype(int), watershed_line=True) > 0
    methodinfo = methodinfo + f"_hmax_{hole_max}_hmin{hole_min}{closingused}"
    # get the shell
    seg = np.logical_xor(bw_filled, dilation(bw_filled, selem=ball(1)))
    ####################
    print(structure_img_smooth.shape, seg.shape)

    #     seg = seg[padsize:-padsize, :, :]
    ###############################
    #     seg = algos.exclude_edgeobjs(seg)
    final_seg = seg > 0
    out = (final_seg * 255).astype(np.uint8)
    out[out > 0] = 255

    fname = (fpath.split("/")[-1]).split(".")[0]
    print(f"saving{fname}")
    # name = fpath.split("/")[-1].split(".")[0]
    # print(name)
    tifffile.imwrite(f'{savepath}{fname}{methodinfo}_GFP.tiff', out, compress=6)
    # writer.save(out)


def segmenttom(fpath, savepath, params, channel="tom20", minarea=None):
    """
       Segments tom20 structures in 3D image stacks.

       Args:
           fpath (str): Path to the input image stack.
           savepath (str): Path to save the segmented output.
           params (dict or list):
               - If dict: Dictionary containing segmentation parameters:
                   - "f2params": List of Frangi filter parameters (list of lists).
               - If list: List of Frangi filter parameters directly.
           channel (str, optional): Name of the channel to segment. Defaults to "tom20".
           minarea (int, optional): Minimum object area for filtering. Defaults to None.

       Returns:
           None
       """
    struct_img0 = get_stack_channel(fpath, channel)

    ################################
    intensity_scaling_param = autoselect_normalization_parameters(struct_img0)
    gaussian_smoothing_sigma = 1
    # intensity normalization
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)
    structure_img_smooth = image_smoothing_gaussian_3d(struct_img, sigma=gaussian_smoothing_sigma)
    ################################
    # f2_param = [[1.2, 0.16]]
    try:
        if isinstance(params, dict):
            f2_param = params["f2params"]
        if isinstance(params, list):
            f2_param = params
    except Exception as e:
        print(e)
    bw = filament_2d_wrapper(structure_img_smooth, f2_param)
    ################################

    seg = remove_small_objects(bw > 0, min_size=minarea, connectivity=1, in_place=False)

    seg = seg > 0
    out = seg.astype(np.uint8)
    out[out > 0] = 255
    s = ""
    for i in range(len(f2_param)):
        for j in range(len(f2_param[i])):
            s = "_".join([s, str(round(f2_param[i][j], 3))])
    # name = fpath.split("/")[-1].split(".")[0]
    ################################
    # fpath = str.replace()
    fname = (fpath.split("/")[-1]).split(".")[0]
    print(f"saving {fname}")
    if not os.path.exists(savepath):
        print(savepath)
        os.mkdir(savepath)
    tifffile.imsave(f'{savepath}{fname}_f2{s}_GFP.tiff', out, compress=6)
    # writer.save(out)


def segmentlampstacks(fpath, savepath, params, channel="lamp1", minarea=None):  # DO NOT USE FILAMENTS FOR NOW
    """
     Segments Lamp1 structures in 3D image stacks.

     Args:
         fpath (str): Path to the input image stack.
         savepath (str): Path to save the segmented output.
         params (dict or list):
             - If dict: Dictionary containing segmentation parameters:
                 - "s2params": List of spot detection parameters (list of lists).
             - If list: List of spot detection parameters directly.
         channel (str, optional): Name of the channel to segment. Defaults to "lamp1".
         minarea (int, optional): Minimum object area for filtering. Defaults to None.

     Returns:
            None
     """
    struct_img0 = get_stack_channel(fpath, channel)

    ################################
    gaussian_smoothing_sigma = 0.5
    intensity_scaling_param = autoselect_normalization_parameters(struct_img0)
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)
    # smoothing with 2d gaussian filter individual slices
    structure_img_smooth = image_smoothing_gaussian_slice_by_slice(struct_img,
                                                                   sigma=gaussian_smoothing_sigma)
    # s2_param = [[4, 0.12], [2, 0.09], [1, 0.02]]
    # f2_param = [[0.75, 0.15]]
    print("params: ", params)
    try:
        if isinstance(params, dict):
            s2_param = params["s2params"]
        if isinstance(params, list):
            s2_param = params
    except Exception as e:
        print(e)
    ################################
    bw_spot = dot_2d_slice_by_slice_wrapper(structure_img_smooth, s2_param)
    # bw_filament = filament_2d_wrapper(structure_img_smooth, f2_param)
    bw = bw_spot
    # bw = np.logical_or(bw_spot, bw_filament)
    ################################
    fill_2d = True
    fill_max_size = 400  # 1600

    bw_fill = hole_filling(bw, 0, fill_max_size, fill_2d)
    seg = remove_small_objects(bw_fill > 0, min_size=minarea, connectivity=1, in_place=False)
    seg = seg > 0
    out = seg.astype(np.uint8)
    out[out > 0] = 255
    ################################
    fname = (fpath.split("/")[-1]).split(".")[0]
    print(f"saving{fname}")
    ss = ""
    for i in range(len(s2_param)):
        for j in range(len(s2_param[i])):
            ss = "_".join([ss, str(s2_param[i][j])])
    sf = ""
    # for i in range(len(f2_param)):
    #     for j in range(len(f2_param[i])):
    #         sf = "_".join([sf, str(f2_param[i][j])])
    tifffile.imwrite(f'{savepath}{fname}s2{ss}_f2{sf}_GFP.tiff', out, compress=6)
    # writer.save(out)


def segmentstgal(fpath, savepath, params, channel="st6gal1", method="nomo", minarea=None, thin_dist=2):
    """
    Segments st6gal1 structures in 3D image stacks using spot detection, Gaussian smoothing,
    and optionally morphological opening and topology-preserving thinning.

    Args:
        fpath (str): Path to the input image stack.
        savepath (str): Path to save the segmented output.
        params (dict or list):
            - If dict: Dictionary containing segmentation parameters:
                - "s3params": List of 3D spot detection parameters (list of lists).
                - "topothin" (optional): Distance for topology-preserving thinning.
            - If list: List of 3D spot detection parameters directly.
        channel (str, optional): Name of the channel to segment. Defaults to "st6gal1".
        method (str, optional): Segmentation method, options:
            - "nomo": Spot detection only.
            - "mo": Spot detection with morphological opening.
            - "mothin": Spot detection, morphological opening, and topology-preserving thinning.
        minarea (int, optional): Minimum object area for filtering. Defaults to None.
        thin_dist (int, optional): Distance for thinning. Used only if "mothin" method is selected.

    Returns:
        None
    """

    struct_img0 = get_stack_channel(fpath, channel)

    intensity_scaling_param = autoselect_normalization_parameters(struct_img0)

    gaussian_smoothing_sigma = 1
    ################################
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)
    # intensityparam = '_'.join(map(str, intensity_scaling_param))

    # smoothing with gaussian filter
    # bw_combine = None
    bw_thin = None

    print(f"method: {method}")
    structure_img_smooth = image_smoothing_gaussian_3d(struct_img, sigma=gaussian_smoothing_sigma)
    try:
        thin_dist_preserve = params["topothin"]
    except:
        pass
    if not method.__contains__("nomo"):
        print("MO")
        bw, object_for_debug = MO(structure_img_smooth, global_thresh_method='tri', object_minArea=300,
                                  return_object=True)
        if not method.__contains__("nothin"):
            bw_thin = topology_preserving_thinning(bw > 0, thin_dist_preserve, thin_dist)
            print("thinned")
    else:
        print("no MO thresholding.")
    ################################
    ## PARAMETERS for this step ##
    print("params: ", params)
    try:
        if isinstance(params, dict):
            s3_param = params["s3params"]
        if isinstance(params, list):
            s3_param = params
    except Exception as e:
        print(e)
    ################################
    bw_extra = dot_3d_wrapper(structure_img_smooth, s3_param)
    if method.__contains__("nomo"):
        bw_combine = bw_extra
    else:
        bw_combine = np.logical_or(bw_extra > 0, bw)
        if not method.__contains__("nothin"):
            bw_combine = np.logical_or(bw_extra > 0, bw_thin)

    ################################
    ## PARAMETERS for this step ##

    ################################
    s3p = ""
    for i in range(len(s3_param)):
        for j in range(len(s3_param[i])):
            s3p = "_".join([s3p, str(s3_param[i][j])])
    seg = remove_small_objects(bw_combine > 0, min_size=minarea, connectivity=1, in_place=False)

    seg = seg > 0
    out = seg.astype(np.uint8)
    out[out > 0] = 255
    name = fpath.split("/")[-1].split(".")[0]
    if method == "nomo":
        tifffile.imwrite(f'{savepath}{name}_s3{s3p}_GFP.tif', out, compress=6)
    else:
        tifffile.imwrite(f'{savepath}{name}_thin{thin_dist_preserve}_method{method}_s3{s3p}_GFP.tif', out, compress=6)
    # writer.save(out)


def segmentfbl(fpath, savepath, params, channel="fbl", minarea=None):
    """
    Segments fbl structures in 3D image stacks using spot detection and Gaussian smoothing.

    Args:
        fpath (str): Path to the input image stack.
        savepath (str): Path to save the segmented output.
        params (dict or list):
            - If dict: Dictionary containing segmentation parameters:
                - "s2params": List of spot detection parameters (list of lists).
            - If list: List of spot detection parameters directly.
        channel (str, optional): Name of the channel to segment. Defaults to "fbl".
        minarea (int, optional): Minimum object area for filtering. Defaults to None.

    Returns:
        None
    """
    struct_img0 = get_stack_channel(fpath, channel)

    ################################
    # intensity_scaling_param = [0, 32]
    gaussian_smoothing_sigma = 1.5
    ################################

    intensity_scaling_param = autoselect_normalization_parameters(struct_img0)
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)
    # smoothing with gaussian filter
    structure_img_smooth = image_smoothing_gaussian_3d(struct_img, sigma=gaussian_smoothing_sigma)

    ################################
    # s2_param = [[1, 0.01]]
    print("params: ", params)
    try:
        if isinstance(params, dict):
            s2_param = params["s2params"]
        if isinstance(params, list):
            s2_param = params
    except Exception as e:
        print(e)
    ################################

    bw = dot_2d_slice_by_slice_wrapper(structure_img_smooth, s2_param)

    s2p = ""
    for i in range(len(s2_param)):
        for j in range(len(s2_param[i])):
            s2p = "_".join([s2p, str(s2_param[i][j])])
    name = fpath.split("/")[-1].split(".")[0]

    final_seg = remove_small_objects(bw > 0, min_size=minarea, connectivity=1, in_place=False)
    final_seg = final_seg > 0
    out = final_seg.astype(np.uint8)
    out[out > 0] = 255

    OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_s2{s2p}_GFP.tiff', compress=6)


def segmenttub(fpath, savepath, params, channel="tuba1b", minarea=None):
    """
    Segments Tubulin structures in 3D image stacks.

    Args:
        fpath (str): Path to the input image stack.
        savepath (str): Path to save the segmented output.
        params (dict): Dictionary containing segmentation parameters:
            - "f3params": List of Frangi filter parameters (list of lists).
        channel (str, optional): Name of the channel to segment. Defaults to "tuba1b".
        minarea (int, optional): Minimum object area for filtering. Defaults to None.

    Returns:
        None
    """
    struct_img0 = get_stack_channel(fpath, channel)

    intensity_scaling_param = autoselect_normalization_parameters(struct_img0)
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)

    # smoothing with edge preserving smoothing
    structure_img_smooth = edge_preserving_smoothing_3d(struct_img)
    ################################
    # f3_param = [[1, 0.01]]
    f3_param = params['f3params']
    ################################

    bw = filament_3d_wrapper(structure_img_smooth, f3_param)

    ################################
    ## PARAMETERS for this step ##

    ################################

    final_seg = remove_small_objects(bw > 0, min_size=minarea, connectivity=1, in_place=False)

    final_seg = final_seg > 0
    out = final_seg.astype(np.uint8)
    out[out > 0] = 255
    f3p = ""
    for i in range(len(f3_param)):
        for j in range(len(f3_param[i])):
            f3p = "_".join([f3p, str(f3_param[i][j])])
    name = fpath.split("/")[-1].split(".")[0]

    OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_f3{f3p}_GFP.tiff', compress=6)


def segmentrab5(fpath, savepath, params, channel="rab5", minarea=None):
    """
    Segments rab5 structures in 3D image stacks using spot detection, Gaussian smoothing, and hole filling.

    Args:
        fpath (str): Path to the input image stack.
        savepath (str): Path to save the segmented output.
        params (dict or list):
            - If dict: Dictionary containing segmentation parameters:
                - "s2params": List of spot detection parameters (list of lists).
            - If list: List of spot detection parameters directly.
        channel (str, optional): Name of the channel to segment. Defaults to "rab5".
        minarea (int, optional): Minimum object area for filtering. Defaults to None.

    Returns:
        None
    """
    struct_img0 = get_stack_channel(fpath, channel)

    ################################
    ## PARAMETERS for this step ##
    # intensity_scaling_param = [3, 19]
    gaussian_smoothing_sigma = 0.5
    ################################
    intensity_scaling_param = autoselect_normalization_parameters(struct_img0)
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)

    # smoothing with 2d gaussian filter slice by slice
    structure_img_smooth = image_smoothing_gaussian_slice_by_slice(struct_img, sigma=gaussian_smoothing_sigma)

    ################################
    # s2_param = [[4, 0.12], [2, 0.09], [1, 0.02]]
    # s2_param = params['s2params']
    print("params: ", params)
    try:
        if isinstance(params, dict):
            s2_param = params["s2params"]
        if isinstance(params, list):
            s2_param = params
    except Exception as e:
        print(e)
    ################################

    bw = dot_2d_slice_by_slice_wrapper(structure_img_smooth, s2_param)

    out1 = bw.astype(np.uint8)
    out1[out1 > 0] = 255
    ################################
    fill_2d = True
    fill_max_size = 400  # 1600

    ################################

    bw_fill = hole_filling(bw, 0, fill_max_size, fill_2d)
    # bw_fill = bw
    final_seg = remove_small_objects(bw_fill > 0, min_size=minarea, connectivity=1, in_place=False)

    final_seg = final_seg > 0
    out = final_seg.astype(np.uint8)
    out[out > 0] = 255
    s2p = ""
    for i in range(len(s2_param)):
        for j in range(len(s2_param[i])):
            s2p = "_".join([s2p, str(s2_param[i][j])])
    name = fpath.split("/")[-1].split(".")[0]

    OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_s2{s2p}_GFP.tiff', compress=6)


def segmentmyh(fpath, savepath, params, channel="myh10", minarea=None, method="2dvessellness", savesegmentation=True,
               returnsegmentation=True):
    """
    Segments myh10 structures in 3D image stacks using.

    Args:
        fpath (str): Path to the input image stack.
        savepath (str): Path to save the segmented output.
        params (dict): Dictionary containing segmentation parameters:
            - "f2params": List of 2D Frangi filter parameters (list of lists).
            - "f3params": List of 3D Frangi filter parameters (list of lists).
            - "both": List of Frangi filter parameters to be used for both 2D and 3D (list of lists).
        channel (str, optional): Name of the channel to segment. Defaults to "myh10".
        minarea (int, optional): Minimum object area for filtering. Defaults to None.
        method (str, optional): Segmentation method, options:
            - "2dvessellness": 2D Frangi filtering.
            - "3dvessellness": 3D Frangi filtering.
            - "both": Both 2D and 3D Frangi filtering combined.
        savesegmentation (bool, optional): Whether to save the segmentation output. Defaults to True.
        returnsegmentation (bool, optional): Whether to return the segmentation output. Defaults to True.

    Returns:
        numpy.ndarray: The segmented image (if returnsegmentation is True).
    """
    struct_img0 = get_stack_channel(fpath, channel)
    # intensity_scaling_param = [2.5, 17]
    intensity_scaling_param = autoselect_normalization_parameters(struct_img0)
    ################################
    # intensity normalization
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)

    # smoothing with edge preserving smoothing
    structure_img_smooth = edge_preserving_smoothing_3d(struct_img)
    ################################
    ################################
    useparam = None
    print(f"using method: {method}")
    if method == "3dvessellness":
        # f3_param = [[1, 0.1]]
        f3param = parseparam(params, 'f3params')
        bw = filament_3d_wrapper(structure_img_smooth, f3param)
        useparam = f3param
    elif method == "2dvessellness":
        f2param = parseparam(params, "f2params")
        bw = filament_2d_wrapper(structure_img_smooth, f2param)
        useparam = f2param
    if method == "both":
        both = parseparam(params, 'both')
        bw3 = filament_3d_wrapper(structure_img_smooth, both)
        bw2 = filament_2d_wrapper(structure_img_smooth, both)
        bw = bw3 | bw2
        useparam = both
    ################################
    final_seg = remove_small_objects(bw > 0, min_size=minarea, connectivity=1, in_place=False)

    final_seg = final_seg > 0
    out = final_seg.astype(np.uint8)
    out[out > 0] = 255
    saveparaminfo = ""
    for i in range(len(useparam)):
        for j in range(len(useparam[i])):
            saveparaminfo = "_".join([saveparaminfo, str(useparam[i][j])])
    name = fpath.split("/")[-1].split(".")[0]
    if savesegmentation:
        OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_method{method}_s3{saveparaminfo}_minarea{minarea}_GFP.tiff',
                           compress=6)
    if returnsegmentation:
        return out
    # OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_f3{f3p}_GFP.tiff', compress=6)


def segmentpxn(fpath, savepath, params, channel="pxn", minarea=None):
    """
    Segments pxn structures in 3D image stacks.

    Args:
        fpath (str): Path to the input image stack.
        savepath (str): Path to save the segmented output.
        params (dict): Dictionary containing segmentation parameters:
            - "f3params": List of Frangi filter parameters (list of lists).
        channel (str, optional): Name of the channel to segment. Defaults to "pxn".
        minarea (int, optional): Minimum object area for filtering. Defaults to None.

    Returns:
        None
    """
    struct_img0 = get_stack_channel(fpath, channel)

    intensity_scaling_param = autoselect_normalization_parameters(struct_img0)
    ################################
    # intensity normalization
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)
    structure_img_smooth = edge_preserving_smoothing_3d(struct_img)
    ################################[[1, 0.01]]
    f3_param = params["f3params"]
    ################################

    bw = filament_3d_wrapper(structure_img_smooth, f3_param)
    ################################

    ################################

    seg = remove_small_objects(bw > 0, min_size=minarea, connectivity=1, in_place=False)
    seg = seg > 0
    out = seg.astype(np.uint8)
    out[out > 0] = 255
    f3p = ""
    for i in range(len(f3_param)):
        for j in range(len(f3_param[i])):
            f3p = "_".join([f3p, str(f3_param[i][j])])
    name = fpath.split("/")[-1].split(".")[0]

    OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_f3{f3p}_GFP.tiff', compress=6)


def segmentdsp(fpath, savepath, params, channel="dsp", minarea=None, method="3dspot", savesegmentation=True,
               returnsegmentation=True):
    """
    Segments dsp structures in 3D image stacks.

    Args:
        fpath (str): Path to the input image stack.
        savepath (str): Path to save the segmented output.
        params (dict): Dictionary containing segmentation parameters:
            - "s2params": List of 2D spot detection parameters (list of lists).
            - "s3params": List of 3D spot detection parameters (list of lists).
            - "both": List of spot detection parameters to be used for both 2D and 3D (list of lists).
        channel (str, optional): Name of the channel to segment. Defaults to "dsp".
        minarea (int, optional): Minimum object area for filtering. Defaults to None.
        method (str, optional): Segmentation method, options:
            - "2dspot": 2D spot detection.
            - "3dspot": 3D spot detection.
            - "both": Both 2D and 3D spot detection combined.
        savesegmentation (bool, optional): Whether to save the segmentation output. Defaults to True.
        returnsegmentation (bool, optional): Whether to return the segmentation output. Defaults to True.

    Returns:
        numpy.ndarray: The segmented image (if returnsegmentation is True).
    """
    struct_img0 = get_stack_channel(fpath, channel)

    ################################
    # intensity_scaling_param = [8000]
    gaussian_smoothing_sigma = 1
    ################################
    intensity_scaling_param = autoselect_normalization_parameters(struct_img0)

    # intensity normalization
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)

    # smoothing with gaussian filter
    structure_img_smooth = image_smoothing_gaussian_slice_by_slice(struct_img,
                                                                   sigma=gaussian_smoothing_sigma)
    ################################
    # [[1, 0.04]]
    useparam = None
    print(f"using method: {method}")
    if method == "2dspot":
        s2_param = parseparam(params, "s2params")
        bw = dot_2d_slice_by_slice_wrapper(structure_img_smooth, s2_param)
        useparam = s2_param
    elif method == "3dspot":
        s3_param = parseparam(params, "s3params")
        bw = dot_3d_wrapper(structure_img_smooth, s3_param)
        useparam = s3_param
    if method == "both":
        both = parseparam(params, 'both')
        sp2 = dot_2d_slice_by_slice_wrapper(structure_img_smooth, both)
        sp3 = dot_3d_wrapper(structure_img_smooth, both)
        bw = sp2 | sp3
        useparam = both
    ################################

    ################################

    ################################
    # watershed
    maskminarea = minarea
    Mask = remove_small_objects(bw > 0, min_size=maskminarea, connectivity=1, in_place=False)
    Seed = dilation(peak_local_max(struct_img, labels=label(Mask), min_distance=2, indices=False),
                    selem=ball(1))
    Watershed_Map = -1 * distance_transform_edt(bw)
    seg = watershed(Watershed_Map, label(Seed), mask=Mask, watershed_line=True)
    seg = remove_small_objects(seg > 0, min_size=minarea, connectivity=1, in_place=False)

    seg = seg > 0
    out = seg.astype(np.uint8)
    out[out > 0] = 255
    saveparaminfo = ""
    for i in range(len(useparam)):
        for j in range(len(useparam[i])):
            saveparaminfo = "_".join([saveparaminfo, str(useparam[i][j])])
    name = fpath.split("/")[-1].split(".")[0]

    if savesegmentation:
        OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_method{method}_s3{saveparaminfo}_minarea{minarea}_GFP.tiff',
                           compress=6)
    if returnsegmentation:
        return out


def segmentslc(fpath, savepath, params, channel="slc25a17", minarea=None, method="3dspot", savesegmentation=True,
               returnsegmentation=True):
    """
    Segments slc25a17 structures in 3D image stacks.

    Args:
        fpath (str): Path to the input image stack.
        savepath (str): Path to save the segmented output.
        params (dict): Dictionary containing segmentation parameters:
            - "s2params": List of 2D spot detection parameters (list of lists).
            - "s3params": List of 3D spot detection parameters (list of lists).
            - "both": List of spot detection parameters to be used for both 2D and 3D (list of lists).
        channel (str, optional): Name of the channel to segment. Defaults to "slc25a17".
        minarea (int, optional): Minimum object area for filtering. Defaults to None.
        method (str, optional): Segmentation method, options:
            - "2dspot": 2D spot detection.
            - "3dspot": 3D spot detection.
            - "both": Both 2D and 3D spot detection combined.
        savesegmentation (bool, optional): Whether to save the segmentation output. Defaults to True.
        returnsegmentation (bool, optional): Whether to return the segmentation output. Defaults to True.

    Returns:
        numpy.ndarray: The segmented image (if returnsegmentation is True).
    """
    struct_img0 = get_stack_channel(fpath, channel)
    gaussian_smoothing_sigma = 1
    ################################
    intensity_scaling_param = autoselect_normalization_parameters(struct_img0)
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)

    # smoothing with gaussian filter
    structure_img_smooth = image_smoothing_gaussian_slice_by_slice(struct_img, sigma=gaussian_smoothing_sigma)
    ################################
    # [[1, 0.04]]
    useparam = None
    print(f"using method: {method}")
    if method == "2dspot":
        s2_param = parseparam(params, "s2params")
        bw = dot_2d_slice_by_slice_wrapper(structure_img_smooth, s2_param)
        useparam = s2_param
    elif method == "3dspot":
        s3_param = parseparam(params, "s3params")
        bw = dot_3d_wrapper(structure_img_smooth, s3_param)
        useparam = s3_param
    if method == "both":
        both = parseparam(params, 'both')
        sp2 = dot_2d_slice_by_slice_wrapper(structure_img_smooth, both)
        sp3 = dot_3d_wrapper(structure_img_smooth, both)
        bw = sp2 | sp3
        useparam = both
    ################################

    # watershed
    maskminarea = minarea
    Mask = remove_small_objects(bw > 0, min_size=maskminarea, connectivity=1, in_place=False)
    Seed = dilation(peak_local_max(struct_img, labels=label(Mask), min_distance=2, indices=False), selem=ball(1))
    Watershed_Map = -1 * distance_transform_edt(bw)
    seg = watershed(Watershed_Map, label(Seed), mask=Mask, watershed_line=True)

    ################################

    ################################

    seg = remove_small_objects(seg > 0, min_size=minarea, connectivity=1, in_place=False)

    seg = seg > 0
    out = seg.astype(np.uint8)
    out[out > 0] = 255
    saveparaminfo = ""
    for i in range(len(useparam)):
        for j in range(len(useparam[i])):
            saveparaminfo = "_".join([saveparaminfo, str(useparam[i][j])])
    name = fpath.split("/")[-1].split(".")[0]

    if savesegmentation:
        OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_method{method}_s3{saveparaminfo}_minarea{minarea}_GFP.tiff',
                           compress=6)
    if returnsegmentation:
        return out


def segmentgja(fpath, savepath, params, channel="gja1", minarea=None, method="3dspot", savesegmentation=True,
               returnsegmentation=True):
    """
    Segments gja1 structures in 3D image stacks using spot detection, Gaussian smoothing, and optional 2D or 3D methods.

    Args:
        fpath (str): Path to the input image stack.
        savepath (str): Path to save the segmented output.
        params (dict): Dictionary containing segmentation parameters:
            - "s2params": List of 2D spot detection parameters (list of lists).
            - "s3params": List of 3D spot detection parameters (list of lists).
            - "both": List of spot detection parameters to be used for both 2D and 3D (list of lists).
        channel (str, optional): Name of the channel to segment. Defaults to "gja1".
        minarea (int, optional): Minimum object area for filtering. Defaults to None.
        method (str, optional): Segmentation method, options:
            - "2dspot": 2D spot detection.
            - "3dspot": 3D spot detection.
            - "both": Both 2D and 3D spot detection combined.
        savesegmentation (bool, optional): Whether to save the segmentation output. Defaults to True.
        returnsegmentation (bool, optional): Whether to return the segmentation output. Defaults to True.

    Returns:
        numpy.ndarray: The segmented image (if returnsegmentation is True).
    """
    struct_img0 = get_stack_channel(fpath, channel)

    ################################
    # intensity_scaling_param = [1, 40]
    intensity_scaling_param = autoselect_normalization_parameters(struct_img0)
    gaussian_smoothing_sigma = 1
    ################################

    # intensity normalization
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)

    # smoothing with 2d gaussian filter slice by slice
    structure_img_smooth = image_smoothing_gaussian_slice_by_slice(struct_img,
                                                                   sigma=gaussian_smoothing_sigma)
    ################################
    ## PARAMETERS for this step ##
    # s3_param = [[1, 0.031]]
    useparam = None
    print(f"using method: {method}")
    if method == "2dspot":
        # s2_param = [[1.5, 0.01]]
        s2_param = parseparam(params, "s2params")
        bw = dot_2d_slice_by_slice_wrapper(structure_img_smooth, s2_param)
        useparam = s2_param
    elif method == "3dspot":
        s3_param = parseparam(params, "s3params")
        bw = dot_3d_wrapper(structure_img_smooth, s3_param)
        useparam = s3_param
    if method == "both":
        both = parseparam(params, 'both')
        sp2 = dot_2d_slice_by_slice_wrapper(structure_img_smooth, both)
        sp3 = dot_3d_wrapper(structure_img_smooth, both)
        bw = sp2 | sp3
        useparam = both
    ################################
    seg = remove_small_objects(bw > 0, min_size=minarea, connectivity=1, in_place=False)
    seg = seg > 0
    out = seg.astype(np.uint8)
    out[out > 0] = 255
    saveparaminfo = ""
    for i in range(len(useparam)):
        for j in range(len(useparam[i])):
            saveparaminfo = "_".join([saveparaminfo, str(useparam[i][j])])
    name = fpath.split("/")[-1].split(".")[0]

    if savesegmentation:
        OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_method{method}_s3{saveparaminfo}_minarea{minarea}_GFP.tiff',
                           compress=6)
    if returnsegmentation:
        return out


def segmentctnnb(fpath, savepath, params, channel="ctnnb1", minarea=None, method="both", savesegmentation=True,
                 returnsegmentation=True):
    """
    Segments ctnnb1 structures in 3D image stacks.

    Args:
        fpath (str): Path to the input image stack.
        savepath (str): Path to save the segmented output.
        params (dict): Dictionary containing segmentation parameters:
            - "s2params": List of 2D spot detection parameters (list of lists).
            - "f2params": List of 2D Frangi filtering parameters (list of lists).
            - "both": List of parameters to be used for both 2D spot detection and Frangi filtering (list of lists).
        channel (str, optional): Name of the channel to segment. Defaults to "ctnnb1".
        minarea (int, optional): Minimum object area for filtering. Defaults to None.
        method (str, optional): Segmentation method, options:
            - "2dspot": 2D spot detection.
            - "2dvessellness": 2D Frangi filtering.
            - "both": Both 2D spot detection and Frangi filtering combined.
        savesegmentation (bool, optional): Whether to save the segmentation output. Defaults to True.
        returnsegmentation (bool, optional): Whether to return the segmentation output. Defaults to True.

    Returns:
        numpy.ndarray: The segmented image (if returnsegmentation is True).
    """

    struct_img0 = get_stack_channel(fpath, channel)

    ################################
    ## PARAMETERS for this step ##
    # intensity_scaling_param =  [4, 27]
    gaussian_smoothing_sigma = 1
    ################################

    intensity_scaling_param = autoselect_normalization_parameters(struct_img0)
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)
    # smoothing with gaussian filter
    structure_img_smooth = image_smoothing_gaussian_3d(struct_img, sigma=gaussian_smoothing_sigma)

    ################################
    ## PARAMETERS for this step ##
    useparam = None
    print(f"using method: {method}")
    if method == "2dspot":
        # s2_param = [[1.5, 0.01]]
        s2_param = parseparam(params, "s2params")
        bw = dot_2d_slice_by_slice_wrapper(structure_img_smooth, s2_param)
        useparam = s2_param
    elif method == "2dvessellness":
        f2param = parseparam(params, "f2params")
        bw = filament_2d_wrapper(structure_img_smooth, f2param)
        useparam = f2param
    if method == "both":
        both = parseparam(params, 'both')
        sp2 = dot_2d_slice_by_slice_wrapper(structure_img_smooth, both)
        bw2 = filament_2d_wrapper(structure_img_smooth, both)
        bw = sp2 | bw2
        useparam = both
    ################################

    ################################

    ################################

    seg = remove_small_objects(bw > 0, min_size=minarea, connectivity=1, in_place=False)
    seg = seg > 0
    out = seg.astype(np.uint8)
    out[out > 0] = 255
    saveparaminfo = ""
    for i in range(len(useparam)):
        for j in range(len(useparam[i])):
            saveparaminfo = "_".join([saveparaminfo, str(useparam[i][j])])
    name = fpath.split("/")[-1].split(".")[0]

    if savesegmentation:
        OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_method{method}_s3{saveparaminfo}_minarea{minarea}_GFP.tiff',
                           compress=6)
    if returnsegmentation:
        return out


def segmentactb(fpath, savepath, params, channel="actb", minarea=0, method="2dvessellness", returnsegmentation=True,
                savesegmentation=True):
    """
    Segments actb structures in 3D image stacks.

    Args:
        fpath (str): Path to the input image stack.
        savepath (str): Path to save the segmented output.
        params (dict): Dictionary containing segmentation parameters:
            - "f2params": List of 2D Frangi filtering parameters (list of lists).
            - "f3params": List of 3D Frangi filtering parameters (list of lists).
            - "both": List of parameters to be used for both 2D and 3D Frangi filtering (list of lists).
        channel (str, optional): Name of the channel to segment. Defaults to "actb".
        minarea (int, optional): Minimum object area for filtering. Defaults to 0.
        method (str, optional): Segmentation method, options:
            - "2dvessellness": 2D Frangi filtering.
            - "3dvessellness": 3D Frangi filtering.
            - "both": Both 2D and 3D Frangi filtering combined.
        savesegmentation (bool, optional): Whether to save the segmentation output. Defaults to True.
        returnsegmentation (bool, optional): Whether to return the segmented image. Defaults to True.

    Returns:
        numpy.ndarray: The segmented image (if returnsegmentation is True).
    """
    struct_img0 = get_stack_channel(fpath, channel)

    ################################
    ## PARAMETERS for this step ##
    # intensity_scaling_param = [2.5, 17]

    intensity_scaling_param = autoselect_normalization_parameters(struct_img0)
    ################################
    # intensity normalization
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)

    # smoothing with edge preserving smoothing
    structure_img_smooth = edge_preserving_smoothing_3d(struct_img)
    ################################
    useparam = None
    print(f"using method: {method}")
    if method == "3dvessellness":
        # f3_param = [[1, 0.1]]
        f3param = parseparam(params, 'f3params')
        bw = filament_3d_wrapper(structure_img_smooth, f3param)
        useparam = f3param
    elif method == "2dvessellness":
        f2param = parseparam(params, "f2params")
        bw = filament_2d_wrapper(structure_img_smooth, f2param)
        useparam = f2param
    if method == "both":
        # f3_param = [[1, 0.1]]
        both = parseparam(params, 'both')
        bw3 = filament_3d_wrapper(structure_img_smooth, both)
        bw2 = filament_2d_wrapper(structure_img_smooth, both)
        bw = bw3 | bw2
        useparam = both
    ################################

    seg = remove_small_objects(bw > 0, min_size=minarea, connectivity=1, in_place=False)
    seg = seg > 0
    out = seg.astype(np.uint8)
    out[out > 0] = 255
    saveparaminfo = ""
    for i in range(len(useparam)):
        for j in range(len(useparam[i])):
            saveparaminfo = "_".join([saveparaminfo, str(useparam[i][j])])
    name = fpath.split("/")[-1].split(".")[0]
    if savesegmentation:
        OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_method{method}_s3{saveparaminfo}_minarea{minarea}_GFP.tiff',
                           compress=6)
    if returnsegmentation:
        return out


def segmenttjp(fpath, savepath, params, channel="tjp1", minarea=0, method="2dvessellness", returnsegmentation=True,
               savesegmentation=True):
    """
    Segments tjp1 structures in 3D image stacks.

    Args:
        fpath (str): Path to the input image stack.
        savepath (str): Path to save the segmented output.
        params (dict): Dictionary containing segmentation parameters:
            - "f2params": List of 2D Frangi filtering parameters (list of lists).
            - "f3params": List of 3D Frangi filtering parameters (list of lists).
            - "both": List of parameters to be used for both 2D and 3D Frangi filtering (list of lists).
        channel (str, optional): Name of the channel to segment. Defaults to "tjp1".
        minarea (int, optional): Minimum object area for filtering. Defaults to 0.
        method (str, optional): Segmentation method, options:
            - "2dvessellness": 2D Frangi filtering.
            - "3dvessellness": 3D Frangi filtering.
            - "both": Both 2D and 3D Frangi filtering combined.
        savesegmentation (bool, optional): Whether to save the segmentation output. Defaults to True.
        returnsegmentation (bool, optional): Whether to return the segmented image. Defaults to True.

    Returns:
        numpy.ndarray: The segmented image (if returnsegmentation is True).
    """
    struct_img0 = get_stack_channel(fpath, channel)

    ################################
    ## PARAMETERS for this step ##
    # intensity_scaling_param = [2.5, 17]

    intensity_scaling_param = autoselect_normalization_parameters(struct_img0)
    ################################
    # intensity normalization
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)

    # smoothing with edge preserving smoothing
    structure_img_smooth = edge_preserving_smoothing_3d(struct_img)
    ################################
    useparam = None
    print(f"using method: {method}")
    if method == "3dvessellness":
        # f3_param = [[1, 0.1]]
        f3param = parseparam(params, 'f3params')
        bw = filament_3d_wrapper(structure_img_smooth, f3param)
        useparam = f3param
    elif method == "2dvessellness":
        f2param = parseparam(params, "f2params")
        bw = filament_2d_wrapper(structure_img_smooth, f2param)
        useparam = f2param
    if method == "both":
        # f3_param = [[1, 0.1]]
        both = parseparam(params, 'both')
        bw3 = filament_3d_wrapper(structure_img_smooth, both)
        bw2 = filament_2d_wrapper(structure_img_smooth, both)
        bw = bw3 | bw2
        useparam = both
    ################################

    seg = remove_small_objects(bw > 0, min_size=minarea, connectivity=1, in_place=False)
    seg = seg > 0
    out = seg.astype(np.uint8)
    out[out > 0] = 255
    saveparaminfo = ""
    for i in range(len(useparam)):
        for j in range(len(useparam[i])):
            saveparaminfo = "_".join([saveparaminfo, str(useparam[i][j])])
    name = fpath.split("/")[-1].split(".")[0]
    if savesegmentation:
        OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_method{method}_s3{saveparaminfo}_minarea{minarea}_GFP.tiff',
                           compress=6)
    if returnsegmentation:
        return out


def segmentcetn2(fpath, savepath, params, channel="cetn2", minarea=None):
    """
    Segments cetn2 structures in 3D image stacks.

    Args:
        fpath (str): Path to the input image stack.
        savepath (str): Path to save the segmented output.
        params (dict or list): Dictionary containing segmentation parameters, or a list of spot detection parameters.
            - If a dictionary, it should contain the key "s3params" with a list of spot detection parameters.
            - If a list, it will be directly used as spot detection parameters.
        channel (str, optional): Name of the channel to segment. Defaults to "cetn2".
        minarea (int, optional): Minimum object area for filtering. Defaults to None.
    """
    struct_img0 = get_stack_channel(fpath, channel)
    gaussian_smoothing_sigma = 1
    intensity_scaling_param = autoselect_normalization_parameters(struct_img0)

    # intensity normalization
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)

    # smoothing with gaussian filter
    structure_img_smooth = image_smoothing_gaussian_slice_by_slice(struct_img,
                                                                   sigma=gaussian_smoothing_sigma)
    # s3_param = [[1, 0.04]]
    print("params: ", params)
    try:
        if isinstance(params, dict):
            s3_param = params['s3params']
        if isinstance(params, list):
            s3_param = params
    except Exception as e:
        print(e)

    ################################

    bw = dot_3d_wrapper(structure_img_smooth, s3_param)

    # watershed
    ################################

    maskminarea = minarea
    Mask = remove_small_objects(bw > 0, min_size=maskminarea, connectivity=1, in_place=False)
    Seed = dilation(peak_local_max(struct_img, labels=label(Mask), min_distance=2, indices=False),
                    selem=ball(1))
    Watershed_Map = -1 * distance_transform_edt(bw)
    seg = watershed(Watershed_Map, label(Seed), mask=Mask, watershed_line=True)
    seg = remove_small_objects(seg > 0, min_size=minarea, connectivity=1, in_place=False)
    seg = seg > 0
    out = seg.astype(np.uint8)
    out[out > 0] = 255
    s3p = ""
    for i in range(len(s3_param)):
        for j in range(len(s3_param[i])):
            s3p = "_".join([s3p, str(s3_param[i][j])])
    name = fpath.split("/")[-1].split(".")[0]

    OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_s3{s3p}_GFP.tiff', compress=6)


def segmentlc3b(fpath, savepath, params, channel="lc3b", minarea=None):
    """
    Segments lc3b structures in 3D image stacks.

    Args:
        fpath (str): Path to the input image stack.
        savepath (str): Path to save the segmented output.
        params (dict or list): Dictionary containing segmentation parameters, or a list of spot detection parameters.
            - If a dictionary, it should contain the key "s3params" with a list of spot detection parameters.
            - If a list, it will be directly used as spot detection parameters.
        channel (str, optional): Name of the channel to segment. Defaults to "lc3b".
        minarea (int, optional): Minimum object area for filtering. Defaults to None.
    """
    struct_img0 = get_stack_channel(fpath, channel)

    ################################
    # intensity_scaling_param = [8000]
    gaussian_smoothing_sigma = 1
    ################################
    intensity_scaling_param = autoselect_normalization_parameters(struct_img0)

    # intensity normalization
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)

    # smoothing with gaussian filter
    structure_img_smooth = image_smoothing_gaussian_slice_by_slice(struct_img,
                                                                   sigma=gaussian_smoothing_sigma)
    ################################
    # s3_param = [[1, 0.04]]
    print("params: ", params)
    try:
        if isinstance(params, dict):
            s3_param = params['s3params']
        if isinstance(params, list):
            s3_param = params
    except Exception as e:
        print(e)
    ################################

    bw = dot_3d_wrapper(structure_img_smooth, s3_param)
    ################################
    # watershed
    maskminarea = minarea
    Mask = remove_small_objects(bw > 0, min_size=maskminarea, connectivity=1, in_place=False)
    Seed = dilation(peak_local_max(struct_img, labels=label(Mask), min_distance=2, indices=False),
                    selem=ball(1))
    Watershed_Map = -1 * distance_transform_edt(bw)
    seg = watershed(Watershed_Map, label(Seed), mask=Mask, watershed_line=True)

    seg = remove_small_objects(seg > 0, min_size=minarea, connectivity=1, in_place=False)
    seg = seg > 0
    out = seg.astype(np.uint8)
    out[out > 0] = 255
    s3p = ""
    for i in range(len(s3_param)):
        for j in range(len(s3_param[i])):
            s3p = "_".join([s3p, str(s3_param[i][j])])
    name = fpath.split("/")[-1].split(".")[0]

    OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_s3{s3p}_GFP.tiff', compress=6)


def parse_segmentation_args(resp_path):
    """
    function to read arguments from a file
    Args:
        resp_path:
    :return: split argument list
    """
    global DEFAULT_DATA_DIR
    fpath = os.path.abspath(resp_path)
    with open(fpath, 'rt') as fi:
        lines = fi.read().split('\n')
    arglist = []
    # for _line in lines:
    #     line = _line.strip()
    #     if len(line) == 0 or line.startswith('#'):
    #         continue
    #     idx = line.find('=')
    #     if idx < 0:
    #         idx = line.find(' ')
    #     if idx < 0:
    #         arglist.append(line)
    #         continue
    #     arglist.append(line[:idx])
    #     arglist.append(line[idx + 1:])
    # DEFAULT_DATA_DIR = os.path.dirname(fpath)
    return arglist


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].startswith('@'):
        try:
            arglist = parse_segmentation_args(sys.argv[1][1:])
        except Exception as ex:
            print('Error parsing response file:', str(ex))
            sys.exit(1)
    else:
        arglist = sys.argv[1:]
