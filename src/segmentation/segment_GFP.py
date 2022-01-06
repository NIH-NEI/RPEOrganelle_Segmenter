import numpy as np
import tifffile
from aicsimageio import AICSImage
from aicsimageio.writers import OmeTiffWriter
from aicssegmentation.core.MO_threshold import MO
from aicssegmentation.core.pre_processing_utils import intensity_normalization, image_smoothing_gaussian_slice_by_slice, \
    suggest_normalization_param, image_smoothing_gaussian_3d, edge_preserving_smoothing_3d
from aicssegmentation.core.seg_dot import dot_2d_slice_by_slice_wrapper, dot_3d_wrapper
from aicssegmentation.core.utils import get_middle_frame, hole_filling, get_3dseed_from_mid_frame
from aicssegmentation.core.utils import topology_preserving_thinning
from aicssegmentation.core.vessel import filament_2d_wrapper, filament_3d_wrapper
from scipy.ndimage import distance_transform_edt
from skimage.feature import peak_local_max
from skimage.measure import label
from skimage.morphology import remove_small_objects, watershed, dilation, ball, closing, disk


def segmentsec61tacks(fpath, savepath, params):
    reader = AICSImage(fpath)
    IMG = reader.data.astype(np.float32)
    struct_img0 = IMG[0, 0, :, :, :].copy()

    ################################
    ## PARAMETERS for this step ##
    intensity_scaling_param = suggest_normalization_param(struct_img0)  # [0.5, 5]  # [1, 2]  # [1.0, 9.0]
    ################################ code modified to return suggested params

    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)

    structure_img_smooth = edge_preserving_smoothing_3d(struct_img)
    ################################
    # scale = params["scale"]             # [0.75, 1]  # [0.75]
    # cutoff = params["cutoff"]            #[0.15]  # [0.15]
    # f2params = [[[scale, cutoff]] for scale in scales for cutoff in cutoffs]]
    f2_param = params["f2params"]
    minareas = [50]  # [10,15,20]
    bw = filament_2d_wrapper(structure_img_smooth.copy(), f2_param)
    for minarea in minareas:
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
        tifffile.imwrite(f'{savepath}{name}{s}_minarea{minarea}.tiff', out, compress=6)
        # writer.save(out)


def segmentlaminstacks(fpath, savepath, params):
    #     padsize = 0
    #     fpath = filepath + file
    reader = AICSImage(fpath)
    IMG = reader.data.astype(np.float32)
    print("getting lamin stack:", IMG.data.shape)
    struct_img0 = IMG[0, 0, :, :, :].copy()
    ################################
    intensity_scaling_param = suggest_normalization_param(struct_img0)
    gscale = (0.5 / 0.216666)
    gaussian_smoothing_sigma = [gscale, 1, 1]
    struct_img = intensity_normalization(struct_img0.copy(), scaling_param=intensity_scaling_param)
    structure_img = struct_img.copy()
    ####################PAD
    #     structure_img = np.pad(struct_img, [(padsize, padsize), (0, 0), (0, 0)])
    ####################
    # smoothing with gaussian filter
    structure_img_smooth = image_smoothing_gaussian_3d(structure_img.copy(), sigma=gaussian_smoothing_sigma)

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
    seed = get_3dseed_from_mid_frame(np.logical_xor(bw_fill_mid_z, bw_dil_z), structure_img.shape, mid_z,
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
    tifffile.imwrite(f'{savepath}{fname}{methodinfo}.tiff', out, compress=6)
    # writer.save(out)


def segmenttomstacks(fpath, savepath, params):
    reader = AICSImage(fpath)
    IMG = reader.data.astype(np.float32)
    print("getting TOM20 stack:", IMG.data.shape)
    struct_img0 = IMG[0, 0, :, :, :].copy()

    ################################
    ## PARAMETERS for this step ##
    intensity_scaling_param = suggest_normalization_param(struct_img0)
    gaussian_smoothing_sigma = 1
    ################################
    # intensity normalization
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)
    structure_img_smooth = image_smoothing_gaussian_3d(struct_img, sigma=gaussian_smoothing_sigma)
    ################################
    # f2_param = [[1.2, 0.16]]
    f2_param = params['f2params']
    bw = filament_2d_wrapper(structure_img_smooth, f2_param)
    ################################
    minArea = 2
    seg = remove_small_objects(bw > 0, min_size=minArea, connectivity=1, in_place=False)

    seg = seg > 0
    out = seg.astype(np.uint8)
    out[out > 0] = 255
    s = ""
    for i in range(len(f2_param)):
        for j in range(len(f2_param[i])):
            s = "_".join([s, str(f2_param[i][j])])
    # name = fpath.split("/")[-1].split(".")[0]
    ################################
    fname = (fpath.split("/")[-1]).split(".")[0]
    print(f"saving{fname}")
    tifffile.imwrite(f'{savepath}{fname}_TOM20_f2{s}.tiff', out, compress=6)
    # writer.save(out)


def segmentlampstacks(fpath, savepath, params):  # DO NOT USE FILAMENTS FOR NOW
    reader = AICSImage(fpath)
    IMG = reader.data.astype(np.float32)
    print("getting LAMP1 stack:", IMG.data.shape)
    struct_img0 = IMG[0, 0, :, :, :].copy()
    ################################
    gaussian_smoothing_sigma = 0.5
    intensity_scaling_param = suggest_normalization_param(struct_img0)
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)
    # smoothing with 2d gaussian filter individual slices
    structure_img_smooth = image_smoothing_gaussian_slice_by_slice(struct_img, sigma=gaussian_smoothing_sigma)
    # s2_param = [[4, 0.12], [2, 0.09], [1, 0.02]]
    # f2_param = [[0.75, 0.15]]
    s2_param = params["s2params"]
    # f2_param = params["f2params"]
    ################################
    bw_spot = dot_2d_slice_by_slice_wrapper(structure_img_smooth, s2_param)
    # bw_filament = filament_2d_wrapper(structure_img_smooth, f2_param)
    bw = bw_spot
    # bw = np.logical_or(bw_spot, bw_filament)
    ################################
    fill_2d = True
    fill_max_size = 400  # 1600
    minArea = 3

    bw_fill = hole_filling(bw, 0, fill_max_size, fill_2d)
    seg = remove_small_objects(bw_fill > 0, min_size=minArea, connectivity=1, in_place=False)
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
    tifffile.imwrite(f'{savepath}{fname}s2{ss}_f2{sf}_LAMP1.tiff', out, compress=6)
    # writer.save(out)


def segmentstgal(fpath, savepath, params):
    reader = AICSImage(fpath)
    IMG = reader.data.astype(np.float32)
    print("getting ST6GAL1 stack:", IMG.data.shape)
    struct_img0 = IMG[0, 0, :, :, :].copy()
    intensity_scaling_param = suggest_normalization_param(struct_img0)

    gaussian_smoothing_sigma = 1
    ################################
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)
    # intensityparam = '_'.join(map(str, intensity_scaling_param))

    # smoothing with gaussian filter
    structure_img_smooth = image_smoothing_gaussian_3d(struct_img, sigma=gaussian_smoothing_sigma)
    bw, object_for_debug = MO(structure_img_smooth, global_thresh_method='tri', object_minArea=300, return_object=True)
    thin_dist_preserve = params["topothin"]
    # thin_dist_preserve=1.6
    # thin_dist_preserve = 0.8
    thin_dist = 1
    bw_thin = topology_preserving_thinning(bw > 0, thin_dist_preserve, thin_dist)

    ################################
    ## PARAMETERS for this step ##
    # s3_param = [[0.8, 0.02]]
    s3_param = params["s3params"]
    ################################
    bw_extra = dot_3d_wrapper(structure_img_smooth, s3_param)

    bw_combine = np.logical_or(bw_extra > 0, bw_thin)
    ################################
    ## PARAMETERS for this step ##
    minArea = 4
    ################################
    s3p = ""
    for i in range(len(s3_param)):
        for j in range(len(s3_param[i])):
            s3p = "_".join([s3p, str(s3_param[i][j])])
    seg = remove_small_objects(bw_combine > 0, min_size=minArea, connectivity=1, in_place=False)

    seg = seg > 0
    out = seg.astype(np.uint8)
    out[out > 0] = 255
    name = fpath.split("/")[-1].split(".")[0]

    tifffile.imwrite(f'{savepath}{name}_thin{thin_dist_preserve}_s3_{s3p}_minAr{minArea}.tif', out, compress=6)
    # writer.save(out)


def segmentfbl(fpath, savepath, params):
    reader = AICSImage(fpath)
    IMG = reader.data.astype(np.float32)
    struct_img0 = IMG[0, 0, :, :, :]
    print("FOR DEBUG: ", IMG.shape, struct_img0.shape)
    ################################
    ## PARAMETERS for this step ##
    # intensity_scaling_param = [0, 32]
    gaussian_smoothing_sigma = 1.5
    ################################

    intensity_scaling_param = suggest_normalization_param(struct_img0)
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)
    # smoothing with gaussian filter
    structure_img_smooth = image_smoothing_gaussian_3d(struct_img, sigma=gaussian_smoothing_sigma)

    ################################
    ## PARAMETERS for this step ##
    # s2_param = [[1, 0.01]]
    s2_param = params["s2params"]

    ################################

    bw = dot_2d_slice_by_slice_wrapper(structure_img_smooth, s2_param)

    ################################
    ## PARAMETERS for this step ##
    minArea = 5
    ################################
    s2p = ""
    for i in range(len(s2_param)):
        for j in range(len(s2_param[i])):
            s2p = "_".join([s2p, str(s2_param[i][j])])
    name = fpath.split("/")[-1].split(".")[0]

    final_seg = remove_small_objects(bw > 0, min_size=minArea, connectivity=1, in_place=False)
    final_seg = final_seg > 0
    out = final_seg.astype(np.uint8)
    out[out > 0] = 255

    OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_s2{s2p}_minar{minArea}.tiff', compress=6)


def segmenttub(fpath, savepath, params):
    reader = AICSImage(fpath)
    IMG = reader.data.astype(np.float32)
    struct_img0 = IMG[0, 0, :, :, :]
    intensity_scaling_param = suggest_normalization_param(struct_img0)
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)
    print("FOR DEBUG: ", IMG.shape, struct_img0.shape)

    # smoothing with edge preserving smoothing
    structure_img_smooth = edge_preserving_smoothing_3d(struct_img)
    ################################
    # f3_param = [[1, 0.01]]
    f3_param = params['f3params']
    ################################

    bw = filament_3d_wrapper(structure_img_smooth, f3_param)

    ################################
    ## PARAMETERS for this step ##
    minArea = 4
    ################################

    final_seg = remove_small_objects(bw > 0, min_size=minArea, connectivity=1, in_place=False)

    final_seg = final_seg > 0
    out = final_seg.astype(np.uint8)
    out[out > 0] = 255
    f3p = ""
    for i in range(len(f3_param)):
        for j in range(len(f3_param[i])):
            f3p = "_".join([f3p, str(f3_param[i][j])])
    name = fpath.split("/")[-1].split(".")[0]

    OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_f3{f3p}_minar{minArea}.tiff', compress=6)


def segmentrab5(fpath, savepath, params):
    reader = AICSImage(fpath)
    IMG = reader.data.astype(np.float32)
    struct_img0 = IMG[0, 0, :, :, :]
    ################################
    ## PARAMETERS for this step ##
    # intensity_scaling_param = [3, 19]
    gaussian_smoothing_sigma = 0.5
    ################################
    intensity_scaling_param = suggest_normalization_param(struct_img0)
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)

    # smoothing with 2d gaussian filter slice by slice
    structure_img_smooth = image_smoothing_gaussian_slice_by_slice(struct_img, sigma=gaussian_smoothing_sigma)

    ################################
    # s2_param = [[4, 0.12], [2, 0.09], [1, 0.02]]
    s2_param = params['s2params']
    ################################

    bw = dot_2d_slice_by_slice_wrapper(structure_img_smooth, s2_param)

    out1 = bw.astype(np.uint8)
    out1[out1 > 0] = 255
    ################################
    fill_2d = True
    fill_max_size = 400  # 1600
    minArea = 4  # 15
    ################################

    bw_fill = hole_filling(bw, 0, fill_max_size, fill_2d)
    # bw_fill = bw
    final_seg = remove_small_objects(bw_fill > 0, min_size=minArea, connectivity=1, in_place=False)

    final_seg = final_seg > 0
    out = final_seg.astype(np.uint8)
    out[out > 0] = 255
    s2p = ""
    for i in range(len(s2_param)):
        for j in range(len(s2_param[i])):
            s2p = "_".join([s2p, str(s2_param[i][j])])
    name = fpath.split("/")[-1].split(".")[0]

    OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_s2{s2p}_minar{minArea}.tiff', compress=6)


def segmentmyh(fpath, savepath, params):
    reader = AICSImage(fpath)
    IMG = reader.data.astype(np.float32)
    struct_img0 = IMG[0, 0, :, :, :]
    print("FOR DEBUG: ", IMG.shape, struct_img0.shape)

    # intensity_scaling_param = [2.5, 17]
    intensity_scaling_param = suggest_normalization_param(struct_img0)
    ################################
    # intensity normalization
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)

    # smoothing with edge preserving smoothing
    structure_img_smooth = edge_preserving_smoothing_3d(struct_img)
    ################################
    # f3_param = [[1, 0.01]]
    f3_param = params['f3params']
    ################################

    bw = filament_3d_wrapper(structure_img_smooth, f3_param)

    ################################
    minArea = 4
    ################################

    final_seg = remove_small_objects(bw > 0, min_size=minArea, connectivity=1, in_place=False)

    final_seg = final_seg > 0
    out = final_seg.astype(np.uint8)
    out[out > 0] = 255
    f3p = ""
    for i in range(len(f3_param)):
        for j in range(len(f3_param[i])):
            f3p = "_".join([f3p, str(f3_param[i][j])])
    name = fpath.split("/")[-1].split(".")[0]

    OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_f3{f3p}_minar{minArea}.tiff', compress=6)


def segmentpxn(fpath, savepath, params):
    reader = AICSImage(fpath)
    IMG = reader.data.astype(np.float32)
    struct_img0 = IMG[0, 0, :, :, :]

    intensity_scaling_param = suggest_normalization_param(struct_img0)
    ################################
    # intensity normalization
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)
    structure_img_smooth = edge_preserving_smoothing_3d(struct_img)
    ################################[[1, 0.01]]
    f3_param = params["f3params"]
    ################################

    bw = filament_3d_wrapper(structure_img_smooth, f3_param)
    ################################
    minArea = 4
    ################################

    seg = remove_small_objects(bw > 0, min_size=minArea, connectivity=1, in_place=False)
    seg = seg > 0
    out = seg.astype(np.uint8)
    out[out > 0] = 255
    f3p = ""
    for i in range(len(f3_param)):
        for j in range(len(f3_param[i])):
            f3p = "_".join([f3p, str(f3_param[i][j])])
    name = fpath.split("/")[-1].split(".")[0]

    OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_f3{f3p}_minar{minArea}.tiff', compress=6)


def segmentdsp(fpath, savepath, params):
    reader = AICSImage(fpath)
    IMG = reader.data.astype(np.float32)
    struct_img0 = IMG[0, 0, :, :, :]
    ################################
    # intensity_scaling_param = [8000]
    gaussian_smoothing_sigma = 1
    ################################
    intensity_scaling_param = suggest_normalization_param(struct_img0)

    # intensity normalization
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)

    # smoothing with gaussian filter
    structure_img_smooth = image_smoothing_gaussian_slice_by_slice(struct_img, sigma=gaussian_smoothing_sigma)
    ################################

    s3_param = params['s3params']  # [[1, 0.04]]
    ################################

    bw = dot_3d_wrapper(structure_img_smooth, s3_param)

    # watershed
    minArea = 4
    Mask = remove_small_objects(bw > 0, min_size=minArea, connectivity=1, in_place=False)
    Seed = dilation(peak_local_max(struct_img, labels=label(Mask), min_distance=2, indices=False), selem=ball(1))
    Watershed_Map = -1 * distance_transform_edt(bw)
    seg = watershed(Watershed_Map, label(Seed), mask=Mask, watershed_line=True)

    ################################
    ## PARAMETERS for this step ##
    minArea = 4
    ################################

    seg = remove_small_objects(seg > 0, min_size=minArea, connectivity=1, in_place=False)

    seg = seg > 0
    out = seg.astype(np.uint8)
    out[out > 0] = 255
    s3p = ""
    for i in range(len(s3_param)):
        for j in range(len(s3_param[i])):
            s3p = "_".join([s3p, str(s3_param[i][j])])
    name = fpath.split("/")[-1].split(".")[0]

    OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_s3{s3p}_minar{minArea}.tiff', compress=6)


def segmentslc(fpath, savepath, params):
    reader = AICSImage(fpath)
    IMG = reader.data.astype(np.float32)
    struct_img0 = IMG[0, 0, :, :, :]

    gaussian_smoothing_sigma = 1
    ################################
    intensity_scaling_param = suggest_normalization_param(struct_img0)

    # intensity normalization
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)

    # smoothing with gaussian filter
    structure_img_smooth = image_smoothing_gaussian_slice_by_slice(struct_img, sigma=gaussian_smoothing_sigma)
    ################################
    ## PARAMETERS for this step ##
    s3_param = params['s3params']  # [[1, 0.04]]
    ################################

    bw = dot_3d_wrapper(structure_img_smooth, s3_param)

    # watershed
    minArea = 4
    Mask = remove_small_objects(bw > 0, min_size=minArea, connectivity=1, in_place=False)
    Seed = dilation(peak_local_max(struct_img, labels=label(Mask), min_distance=2, indices=False), selem=ball(1))
    Watershed_Map = -1 * distance_transform_edt(bw)
    seg = watershed(Watershed_Map, label(Seed), mask=Mask, watershed_line=True)

    ################################
    ## PARAMETERS for this step ##
    minArea = 4
    ################################

    seg = remove_small_objects(seg > 0, min_size=minArea, connectivity=1, in_place=False)

    seg = seg > 0
    out = seg.astype(np.uint8)
    out[out > 0] = 255
    s3p = ""
    for i in range(len(s3_param)):
        for j in range(len(s3_param[i])):
            s3p = "_".join([s3p, str(s3_param[i][j])])
    name = fpath.split("/")[-1].split(".")[0]

    OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_s3{s3p}_minar{minArea}.tiff', compress=6)


def segmentgja(fpath, savepath, params):
    reader = AICSImage(fpath)
    IMG = reader.data.astype(np.float32)
    struct_img0 = IMG[0, 0, :, :, :]
    ################################
    ## PARAMETERS for this step ##
    # intensity_scaling_param = [1, 40]
    intensity_scaling_param = suggest_normalization_param(struct_img0)
    gaussian_smoothing_sigma = 1
    ################################

    # intensity normalization
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)

    # smoothing with 2d gaussian filter slice by slice
    structure_img_smooth = image_smoothing_gaussian_slice_by_slice(struct_img, sigma=gaussian_smoothing_sigma)
    ################################
    # s3_param = [[1, 0.031]]
    s3_param = params['s3params']

    ################################

    bw = dot_3d_wrapper(structure_img_smooth, s3_param)

    ################################
    minArea = 5
    ################################

    seg = remove_small_objects(bw > 0, min_size=minArea, connectivity=1, in_place=False)
    seg = seg > 0
    out = seg.astype(np.uint8)
    out[out > 0] = 255
    s3p = ""
    for i in range(len(s3_param)):
        for j in range(len(s3_param[i])):
            s3p = "_".join([s3p, str(s3_param[i][j])])
    name = fpath.split("/")[-1].split(".")[0]

    OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_s3{s3p}_minar{minArea}.tiff', compress=6)


def segmentctnnb(fpath, savepath, params):
    reader = AICSImage(fpath)
    IMG = reader.data.astype(np.float32)
    struct_img0 = IMG[0, 0, :, :, :]
    ################################
    ## PARAMETERS for this step ##
    # intensity_scaling_param =  [4, 27]
    gaussian_smoothing_sigma = 1
    ################################

    intensity_scaling_param = suggest_normalization_param(struct_img0)
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)
    # smoothing with gaussian filter
    structure_img_smooth = image_smoothing_gaussian_3d(struct_img, sigma=gaussian_smoothing_sigma)

    ################################
    ## PARAMETERS for this step ##
    # s2_param = [[1.5, 0.01]]
    s2_param = params['s2params']
    ################################

    bw = dot_2d_slice_by_slice_wrapper(structure_img_smooth, s2_param)

    ################################
    ## PARAMETERS for this step ##
    minArea = 5
    ################################

    seg = remove_small_objects(bw > 0, min_size=minArea, connectivity=1, in_place=False)
    seg = seg > 0
    out = seg.astype(np.uint8)
    out[out > 0] = 255
    s3p = ""
    for i in range(len(s2_param)):
        for j in range(len(s2_param[i])):
            s3p = "_".join([s3p, str(s2_param[i][j])])
    name = fpath.split("/")[-1].split(".")[0]

    OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_s3{s3p}_minar{minArea}.tiff', compress=6)


def segmentactb(fpath, savepath, params):
    reader = AICSImage(fpath)
    IMG = reader.data.astype(np.float32)
    struct_img0 = IMG[0, 0, :, :, :]
    ################################
    ## PARAMETERS for this step ##
    # intensity_scaling_param = [2.5, 17]

    intensity_scaling_param = suggest_normalization_param(struct_img0)
    ################################
    # intensity normalization
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)

    # smoothing with edge preserving smoothing
    structure_img_smooth = edge_preserving_smoothing_3d(struct_img)
    ################################
    # f3_param = [[1, 0.1]]
    f3_param = params['f3params']
    ################################

    bw = filament_3d_wrapper(structure_img_smooth, f3_param)
    ################################
    minArea = 4
    ################################

    seg = remove_small_objects(bw > 0, min_size=minArea, connectivity=1, in_place=False)
    seg = seg > 0
    out = seg.astype(np.uint8)
    out[out > 0] = 255
    s3p = ""
    for i in range(len(f3_param)):
        for j in range(len(f3_param[i])):
            s3p = "_".join([s3p, str(f3_param[i][j])])
    name = fpath.split("/")[-1].split(".")[0]

    OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_s3{s3p}_minar{minArea}.tiff', compress=6)


def segmentcetn2(fpath, savepath, params):
    reader = AICSImage(fpath)
    IMG = reader.data.astype(np.float32)
    struct_img0 = IMG[0, 0, :, :, :]
    ################################
    gaussian_smoothing_sigma = 1
    ################################
    intensity_scaling_param = suggest_normalization_param(struct_img0)

    # intensity normalization
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)

    # smoothing with gaussian filter
    structure_img_smooth = image_smoothing_gaussian_slice_by_slice(struct_img, sigma=gaussian_smoothing_sigma)
    ################################
    # s3_param = [[1, 0.04]]
    s3_param = params['s3params']
    ################################

    bw = dot_3d_wrapper(structure_img_smooth, s3_param)

    # watershed
    minArea = 4
    Mask = remove_small_objects(bw > 0, min_size=minArea, connectivity=1, in_place=False)
    Seed = dilation(peak_local_max(struct_img, labels=label(Mask), min_distance=2, indices=False), selem=ball(1))
    Watershed_Map = -1 * distance_transform_edt(bw)
    seg = watershed(Watershed_Map, label(Seed), mask=Mask, watershed_line=True)

    ################################
    minArea = 4
    ################################

    seg = remove_small_objects(seg > 0, min_size=minArea, connectivity=1, in_place=False)
    seg = seg > 0
    out = seg.astype(np.uint8)
    out[out > 0] = 255
    s3p = ""
    for i in range(len(s3_param)):
        for j in range(len(s3_param[i])):
            s3p = "_".join([s3p, str(s3_param[i][j])])
    name = fpath.split("/")[-1].split(".")[0]

    OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_s3{s3p}_minar{minArea}.tiff', compress=6)


def segmentlc3b(fpath, savepath, params):
    reader = AICSImage(fpath)
    IMG = reader.data.astype(np.float32)
    struct_img0 = IMG[0, 0, :, :, :]
    ################################
    # intensity_scaling_param = [8000]
    gaussian_smoothing_sigma = 1
    ################################
    intensity_scaling_param = suggest_normalization_param(struct_img0)

    # intensity normalization
    struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)

    # smoothing with gaussian filter
    structure_img_smooth = image_smoothing_gaussian_slice_by_slice(struct_img, sigma=gaussian_smoothing_sigma)
    ################################
    # s3_param = [[1, 0.04]]
    s3_param = params['s3params']
    ################################

    bw = dot_3d_wrapper(structure_img_smooth, s3_param)

    # watershed
    minArea = 4
    Mask = remove_small_objects(bw > 0, min_size=minArea, connectivity=1, in_place=False)
    Seed = dilation(peak_local_max(struct_img, labels=label(Mask), min_distance=2, indices=False), selem=ball(1))
    Watershed_Map = -1 * distance_transform_edt(bw)
    seg = watershed(Watershed_Map, label(Seed), mask=Mask, watershed_line=True)
    ################################
    ## PARAMETERS for this step ##
    minArea = 4
    ################################

    seg = remove_small_objects(seg > 0, min_size=minArea, connectivity=1, in_place=False)
    seg = seg > 0
    out = seg.astype(np.uint8)
    out[out > 0] = 255
    s3p = ""
    for i in range(len(s3_param)):
        for j in range(len(s3_param[i])):
            s3p = "_".join([s3p, str(s3_param[i][j])])
    name = fpath.split("/")[-1].split(".")[0]

    OmeTiffWriter.save(data=out, uri=f'{savepath}{name}_s3{s3p}_minar{minArea}.tiff', compress=6)
