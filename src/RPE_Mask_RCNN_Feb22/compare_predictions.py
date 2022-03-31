
import os, sys

import numpy as np
import tifffile

import src.RPE_Mask_RCNN_Feb22.rpefs
from src.RPE_Mask_RCNN_Feb22 import rpesegm,rpeutil

hint = '''Usage: python compare_predictions.py </path/to/primary/rpe/annotations.csv> </path/to/predictions[.csv|.tif]>

First parameter is a path to CSV file (_RPE.csv) containing RPE annotations, usually from the Manual directory.
Second parameter is either another _RPE.csv (usually from Predicted directory) or
    16-bit TIFF stack (pixel value = Slice ID) containing predictions from Cellpose.
If 2nd parameter is TIFF, a conversion is performed and a pair of files _RPE.csv/_RPE.tif is created in the
directory containing Cellpose TIFF with names matching the first parameter base name.

*IMPORTANT* To prevent the first file (_RPE.csv) from being overwritten,
place the second one (Cellpose's .tif) in a different directory.

Examples:
python compare_predictions.py C:\\rpemrcnn\\W1\\Manual\\P1-W1-ZO1_D02_F002_Actin_RPE.csv C:\\rpemrcnn\\W1\\Cellpose\\W1.tif
python compare_predictions.py C:\\rpemrcnn\W1\Manual\\P1-W1-ZO1_D02_F002_Actin_RPE.csv C:\\rpemrcnn\\W1\\Predicted\\P1-W1-ZO1_D02_F002_Actin_RPE.csv
'''
    
iou_pcts = (50, 60, 70, 75, 80, 85, 90, 95)

if __name__ == '__main__':
    
    if len(sys.argv) < 3:
        print(hint)
        sys.exit(0)
        
    base_csv = os.path.abspath(sys.argv[1])
    cp_tif = os.path.abspath(sys.argv[2])
    assert os.path.isfile(base_csv), 'CSV annotation file must exist'
    assert os.path.isfile(cp_tif), 'TIFF ID file (predictions from cellpose) must exist'
    
    base_dir, base_fn = os.path.split(base_csv)
    base_bn, base_ext = os.path.splitext(base_fn)
    
    cp_dir, cp_fn = os.path.split(cp_tif)
    cp_bn, cp_ext = os.path.splitext(cp_fn)
    if cp_ext.lower() in ('.tif', '.tiff'):
        assert cp_dir != base_dir, 'CSV annotations and TIFF ID predictions must be in different directories'
        
        cmp_csv = os.path.join(cp_dir, base_bn+'.csv')
        cmp_tif = os.path.join(cp_dir, base_bn+'.tif')
        
        print('Read:', cp_tif)
        img_data = tifffile.imread(cp_tif)
        assert len(img_data.shape) == 3 and img_data.dtype == np.uint16, 'Expected 16-bit TIFF stack'
        
        msk_data = np.empty(shape=img_data.shape, dtype=np.uint8)
        rpesegm.import_mask_id_3d(img_data, msk_data, cmp_csv);
        print('Write:', cmp_tif)
        tifffile.imwrite(cmp_tif, msk_data, compress=6, photometric='minisblack')
    else:
        cmp_csv = cp_tif
        base_tif = os.path.join(base_dir, base_bn+'.tif')
        cmp_tif = os.path.join(cp_dir, cp_bn+'.tif')
        if os.path.isfile(base_tif):
            img_data = tifffile.imread(base_tif)
        elif os.path.isfile(cmp_tif):
            img_data = tifffile.imread(cmp_tif)
        else:
            print('Neither %s nor %s exists.\nCannot perform comparison because the exact stack size is not known' % \
                  (base_tif, cmp_tif))
            sys.exit(2)
    
    d, h, w = tuple(img_data.shape)
        
    print('Compare:', base_csv)
    print('     To:', cmp_csv)
    res = rpeutil.compare_3d_annotations(w, h, d, base_csv, cmp_csv)
    print ('==2D== Primary Slices:', res.base_slices, ' - Secondary Slices:', res.cmp_slices)
    print ('False Positives:', res.f_pos, ' | False Negatives:', res.f_neg,
           ' | Fragmented:', res.fragm, ' | Fused:', res.fused)
    for idx in iou_pcts:
        cnt = 0
        for i in range(idx, 100):
            cnt += res.pct_match[i]
        print ('IOU_2D_PCT_%d : %5d / %5.1f%%' % (idx, cnt, cnt*100./res.base_slices))
    print ('==3D== Primary Cells:', res.base_cells, ' - Secondary Cells:', res.cmp_cells)
    print ('False Positives:', res.f_pos_3d, ' | False Negatives:', res.f_neg_3d,
           ' | Fragmented:', res.fragm_3d, ' | Fused:', res.fused_3d)
    for idx in iou_pcts:
        cnt = 0
        for i in range(idx, 100):
            cnt += res.pct_match_3d[i]
        print ('IOU_3D_PCT_%d : %5d / %5.1f%%' % (idx, cnt, cnt*100./res.base_cells))
    
    sys.exit(0)
    
