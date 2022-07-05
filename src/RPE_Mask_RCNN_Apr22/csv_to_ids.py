import sys, os, datetime

import numpy as np
import tifffile
from rpefs import RpeTiffStack
import rpesegm

def usage():
    txt = '''Usage: python csv_to_ids.py <csvfile>.csv [<mask>.tif|<stack>.ome.tif|<RPEMeta>.rpe.json]
Output (if successful): <csvfile>_ids.tif

If the second parameter is omitted, the script will try to find corresponding tif (either source stack or mask).
This is only needed to determine the exact shape (n_frames, height, width) of the output .tif stack.
'''
    print(txt)
    
def error(txt):
    print(txt)
    print()
    usage()
    sys.exit(1)

if __name__ == '__main__':
    
    if len(sys.argv) < 2:
        usage()
        sys.exit(0)
        
    csvpath = os.path.abspath(sys.argv[1])
    if not os.path.isfile(csvpath):
        error(f'ERROR: No such file: {csvpath}')
        
    print(f'Input segmentation CSV: {csvpath}')
    basedir, fn = os.path.split(csvpath)
    bn, ext = os.path.splitext(fn)
    tifpath = os.path.join(basedir, bn+'_ids.tif')
    num_dilations = 0
    if bn.find('_Actin') >= 0:
        # Extra dilations for Actin only
        num_dilations = 0
    
    reftif = None
    if len(sys.argv) > 2:
        reftif = os.path.abspath(sys.argv[2])
        if not os.path.isfile(reftif):
            error(f'ERROR: No such file: {reftif}')
    
    if reftif is None:
        # Try accompanying mask .tif
        reftif = os.path.join(basedir, bn+'.tif')
        if not os.path.isfile(reftif):
            reftif = None
            
    if reftif is None:
        # Try source tif or .rpe.json

        # Chop off _DNA_RPE or _Actin_RPE from the end of the basename
        for ch in ('_DNA', '_Actin'):
            idx = bn.find(ch)
            if idx >= 0:
                bn = bn[:idx]
                break
        # Go one directory up to look for the source file
        srcdir = os.path.dirname(basedir)
        
        # Try .rpe.json, then .ome.tif
        reftif = os.path.join(srcdir, bn+'.rpe.json')
        if not os.path.isfile(reftif):
            reftif = os.path.join(srcdir, bn+'.ome.tif')
            if not os.path.isfile(reftif):
                error(f'ERROR: could not find an accompanying TIFF for {csvpath}')
        
    # Determine ref.stack shape (n_frames, height, width)
    if reftif.endswith('.rpe.json') or reftif.endswith('.ome.tif'):
        stk = RpeTiffStack(reftif)
        if stk.n_frames < 2:
            error(f'ERROR: bad format {reftif}')
        shape = (stk.n_frames, stk.height, stk.width)
    else:
        stk = tifffile.TiffFile(reftif)
        n_pages = len(stk.pages)
        if n_pages < 2:
            error(f'ERROR: bad format {reftif}')
        page = stk.pages[0]
        shape = (n_pages, page.shape[0], page.shape[1])
    
    print(f'Reference Stack: {reftif}, shape: {shape}')

    mdata = np.empty(shape=shape, dtype=np.uint16)
    # The last parameter is number of 1-pix dilations of the cells before exporting
    id = rpesegm.export_mask_id_3d(mdata, csvpath, num_dilations)
    
    if id < 1:
        print(f'ERROR: No objects imported from {csvpath}. Wrong CSV file maybe?')
        print('No output generated.')
    else:
        print()
        print(f'Write output segmentation (TIFF-16) to: {tifpath}')
        tifffile.imwrite(tifpath, mdata, compress=6, photometric='minisblack')
    
    print('OK.')
    sys.exit(0)












