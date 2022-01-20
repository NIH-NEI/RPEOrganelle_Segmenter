import os

# DEBUG: Disable GPU(s)
# os.environ['CUDA_VISIBLE_DEVICES']='-1'

import sys
import datetime
import csv
import numpy as np
import scipy
import skimage
import tifffile
from src.RPE_Mask_RCNN_Jan22.savgol import SavitzkyGolay2D

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
#tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

# Root directory of the project
ROOT_DIR = os.path.abspath(os.path.dirname(__file__))

# Import Mask RCNN
if not ROOT_DIR in sys.path:
    sys.path.append(ROOT_DIR)  # To find local version of the library
from src.RPE_Mask_RCNN_Jan22.mrcnn.config import Config
from src.RPE_Mask_RCNN_Jan22.mrcnn import model as modellib

from src.RPE_Mask_RCNN_Jan22 import rpesegm, rpeutil
from src.RPE_Mask_RCNN_Jan22.util import slice_area
from src.RPE_Mask_RCNN_Jan22.rpefs import RpeTiffStack

class RpeInferenceConfig(Config):
    # Give the configuration a recognizable name
    NAME = 'TBD'
    LOGS_DIR = 'logs'
    PREDICTED_SUBDIR = 'Predicted'
    # Image tile size (width and height)
    TILESIZE = 768
    # Output raw predictions (.raw.tif)
    WRITE_RAW = False
    # Adjust borders if True (also filter out "phantom slices")
    BORDER_ADJUST = 0 # 0, rpeutil.FIX_BORDERS_DNA or rpeutil.FIX_BORDERS_ACTIN
    # Number of GPUs to use (1 if no GPUs)
    GPU_COUNT = 1
    # Adjust down if you use a smaller GPU.
    IMAGES_PER_GPU = 1
    # Number of classes (including background)
    NUM_CLASSES = 1 + 1  # Background + DNA/Actin
    # Skip detections with < nn% confidence
    DETECTION_MIN_CONFIDENCE = 0.5
    DETECTION_MAX_INSTANCES = 1000
    #
    
def iter_stacks(data_dir, rpefile, chname, recurse=False):
    def open_stack(fpath, chname):
        tifstk = RpeTiffStack(fpath)
        if not tifstk.hasChannel(chname):
            print ('%s -- not a valid RPE Meta file, skipped.' % (fpath,))
            return None
        return tifstk
    for rpath in rpefile:
        #print({'data_dir':data_dir, 'rpath':rpath})
        fpath = os.path.abspath(os.path.join(data_dir, rpath))
        if os.path.isfile(fpath):
            tifstk = open_stack(fpath, chname)
            if not tifstk is None:
                yield tifstk
        elif os.path.isdir(fpath):
            fdir = fpath
            print('Scanning directory:', fdir, '\n  with recurse option' if recurse else '')
            for fn in os.listdir(fdir):
                fpath = os.path.join(fdir, fn)
                if os.path.isdir(fpath):
                    if recurse:
                        for tifstk in iter_stacks('', [fpath], chname, False):
                            yield tifstk
                    continue;
                if not fn.endswith('.rpe.json'):
                    continue
                tifstk = open_stack(fpath, chname)
                if not tifstk is None:
                    yield tifstk
    #

def assemble_mask_stack(masks):
    fr_mask = None
    nobj = 0
    h = len(masks)
    w = None
    for i, mask in enumerate(masks):
        mask = mask.astype(np.uint8)
        mask[mask != 0] = 0xFF
        if fr_mask is None:
            w = mask.shape[0]
            nobj = mask.shape[1]
            fr_mask = np.empty(shape=(nobj, h, w), dtype=np.uint8)
        for j in range(nobj):
            fr_mask[j, i, :] = mask[:, j]
    return fr_mask
    #

class RpeDnaSegmenter(object):
    def __init__(self, weights_path, args):
        self.weights_path = weights_path
        self.args = args
        #
        class DnaInferenceConfig(RpeInferenceConfig): 
            NAME = 'DNA'
            LOGS_DIR = os.path.join(args.data_dir, 'logs')
            PREDICTED_SUBDIR = 'Predicted'
            TILESIZE = args.tilesize
            WRITE_RAW = args.raw
            BORDER_ADJUST = 0 if args.no_adjust else rpeutil.FIX_BORDERS_DNA
            GPU_COUNT = args.gpu_count
            IMAGES_PER_GPU = args.gpu_images
        if args.no_adjust:
            DnaInferenceConfig.PREDICTED_SUBDIR = 'Predicted_N'
            DnaInferenceConfig.BORDER_ADJUST = 0
            print ('Predicted object border adjustment (post-processing) disabled.')
            print ('Output directory changed to', DnaInferenceConfig.PREDICTED_SUBDIR)
        self.config = DnaInferenceConfig()
        #
        self.model = modellib.MaskRCNN(mode="inference", config=self.config, model_dir=self.config.LOGS_DIR)
        self.model.load_weights(self.weights_path, by_name=True)
        #
        sg = SavitzkyGolay2D(5, 2)
        self.smooth_kern = sg.kernel(0, 0)
        #
        self.tifstk = None
        self.sc_norm = 1.
        self.glob_otsu = None
        self.tiles = []
        self.mask_cache = {}
    #
    # called from _iter_batches()
    def _iter_tiles(self):
        chname = self.config.NAME
        tilesize = self.config.TILESIZE
        for cframe in range(self.tifstk.n_frames):
            fr_data = self.tifstk.getFrame(cframe, chname, nocache=True)
            #
            if self.config.BORDER_ADJUST:
                fr_smooth = scipy.signal.convolve2d(fr_data, self.smooth_kern, boundary='symm', mode='same')
                fr_mask = np.zeros(fr_data.shape, np.uint8)
                fr_mask[fr_smooth > self.glob_otsu] = 0xFF
                del fr_smooth
                self.mask_cache[cframe] = fr_mask
            #
            fr_data = fr_data.astype(np.float32) * self.sc_norm
            fr_data[fr_data > 65530.] = 65530.
            fr_data = fr_data.astype(np.uint16)
            for x0, x1, y0, y1 in self.tiles:
                image = skimage.color.gray2rgb(fr_data[y0:y1+1, x0:x1+1])
                yield cframe, (x0, x1, y0, y1), image
        for j in range(self.config.BATCH_SIZE):
            image = np.zeros(shape=(tilesize, tilesize, 3), dtype=np.uint8)
            yield cframe, (0, tilesize-1, 0, tilesize-1), image
    #
    # called from segment()
    def _iter_batches(self):
        cframes = []
        ctiles = []
        cimages = []
        for cframe, ctile, cimage in self._iter_tiles():
            cframes.append(cframe)
            ctiles.append(ctile)
            cimages.append(cimage)
            if len(cimages) == self.config.BATCH_SIZE:
                results = self.model.detect(cimages, verbose=0)
                for _cframe, _ctile, res in zip(cframes, ctiles, results):
                    res['frame'] = _cframe
                    res['tile'] = _ctile
                    yield res
                cframes = []
                ctiles = []
                cimages = []
    #
    def segment(self, tifstk):
        all_start_ts = datetime.datetime.now()
        self.tifstk = tifstk
        chname = self.config.NAME
        tilesize = self.config.TILESIZE
        otsu = tifstk.getChannelOtsu(chname)
        self.sc_norm = 16383./otsu
        
        if self.config.BORDER_ADJUST:
            otsu3 = tifstk.getChannelOtsu3(chname)
            self.glob_otsu = int(otsu3[0])
        else:
            self.glob_otsu = None
            
        self.tiles = slice_area((tifstk.height, tifstk.width), (tilesize, tilesize))
        self.mask_cache = {}
    
        bn = tifstk.base_name+'_'+chname+'_RPE'
        ml_dir = tifstk.cache.makedir(self.config.PREDICTED_SUBDIR)
        ml_csv = os.path.join(ml_dir, bn+'.csv')
        ml_tif = os.path.join(ml_dir, bn+'.tif')
        ml_raw_tif = os.path.join(ml_dir, bn+'.raw.tif')
        ml_score_csv = os.path.join(ml_dir, bn+'.score.csv')
        if self.config.WRITE_RAW:
            fo = open(ml_score_csv, 'w')
            wr = csv.writer(fo, dialect='excel')
            wr.writerow(['Frame', 'XC', 'YC', 'Score'])
        #
        NT = len(self.tiles)
        tileg = self._iter_batches()
        dna_frames = np.empty(shape=(tifstk.n_frames, tifstk.height, tifstk.width), dtype=np.uint8)
        #
        #print ('Segmenting %s of %s' % (chname, tifstk.base_name))
        for ifr in range(tifstk.n_frames):
            start_ts = datetime.datetime.now()
            
            scores = []
            particles = []
            for itl in range(NT):
                results = next(tileg)
                x0, x1, y0, y1 = results['tile']
                tile_mask = assemble_mask_stack(results.pop('masks'))
                if tile_mask is None: continue
                scores.extend([float(x) for x in results['scores']])
                particles.extend(rpeutil.mask_rcnn_to_particles(tile_mask, results['rois'], x0, y0))
                del tile_mask
    
            if self.config.BORDER_ADJUST:
                fr_mask = self.mask_cache.pop(ifr)
            else:
                fr_mask = np.zeros(shape=(tifstk.height, tifstk.width), dtype=np.uint8)
            flat_scores = rpeutil.recombine_flat_particles(fr_mask, particles, scores, self.config.BORDER_ADJUST)
            dna_frames[ifr] = fr_mask
            del fr_mask
            #
            if self.config.WRITE_RAW:
                n_scores = flat_scores[0]
                j = 1
                for i in range(n_scores):
                    xc = flat_scores[j]
                    yc = flat_scores[j+1]
                    sc = flat_scores[j+2]
                    j += 3
                    wr.writerow([ifr, xc, yc, sc])
            #
            elapsed = str(datetime.datetime.now() - start_ts)
            print ('Frame %d of %d segmented in %s' % ((ifr+1), tifstk.n_frames, elapsed))
        #  
        if self.config.WRITE_RAW:
            print ('Write raw predictions:', ml_raw_tif)
            tifffile.imwrite(ml_raw_tif, dna_frames, compress=6, photometric='minisblack')
            print ('Write scores:', ml_score_csv)
            fo.close()
        #
        print ('Assembling 3D %s of %s' % (chname, tifstk.base_name))
        start_ts = datetime.datetime.now()
        rpesegm.assemble_basic(dna_frames, ml_csv)
        elapsed = str(datetime.datetime.now() - start_ts)
        print('3D Assembly done in', elapsed)
        
        print ('Write:', ml_tif)
        tifffile.imwrite(ml_tif, dna_frames, compress=6, photometric='minisblack')
        
        elapsed = str(datetime.datetime.now() - all_start_ts)
        print ('Segmenting %s of %s done in: %s' % (chname, tifstk.base_name, elapsed))
    #

class RpeActinSegmenter(object):
    RGB_FRAMES = ['DNA', 'Actin', 'Membrane']
    def __init__(self, weights_path, args):
        self.weights_path = weights_path
        self.args = args
        #
        class ActinInferenceConfig(RpeInferenceConfig): 
            NAME = 'Actin'
            LOGS_DIR = os.path.join(args.data_dir, 'logs')
            PREDICTED_SUBDIR = 'Predicted'
            TILESIZE = args.tilesize
            WRITE_RAW = args.raw
            BORDER_ADJUST = 0 if args.no_adjust else rpeutil.FIX_BORDERS_ACTIN
            GPU_COUNT = args.gpu_count
            IMAGES_PER_GPU = args.gpu_images
        if args.no_adjust:
            ActinInferenceConfig.PREDICTED_SUBDIR = 'Predicted_N'
            ActinInferenceConfig.BORDER_ADJUST = 0
            print ('Predicted object border adjustment (post-processing) disabled.')
            print ('Output directory changed to', ActinInferenceConfig.PREDICTED_SUBDIR)
        self.config = ActinInferenceConfig()
        #
        self.model = modellib.MaskRCNN(mode="inference", config=self.config, model_dir=self.config.LOGS_DIR)
        self.model.load_weights(self.weights_path, by_name=True)
        #
        sg = SavitzkyGolay2D(5, 2)
        self.smooth_kern = sg.kernel(0, 0)
        #
        self.tifstk = None
        self.sc_norm = {}
        self.glob_otsu = None
        self.tiles = []
        self.mask_cache = {}
    #
    def _prepare_rgb_frame(self, iframe):
        fr_data = np.empty(shape=(self.tifstk.height, self.tifstk.width, 3), dtype=np.uint8)
        for i, chname in enumerate(self.RGB_FRAMES):
            norm = self.sc_norm[chname]
            ch_data = self.tifstk.getFrame(iframe, chname, nocache=True).astype(np.float32)
            #
            if chname == self.config.NAME and self.config.BORDER_ADJUST == rpeutil.FIX_BORDERS_DNA:
                fr_smooth = scipy.signal.convolve2d(ch_data, self.smooth_kern, boundary='symm', mode='same')
                fr_mask = np.zeros(ch_data.shape, np.uint8)
                fr_mask[fr_smooth <= self.glob_otsu] = 0xFF
                del fr_smooth
                self.mask_cache[iframe] = fr_mask
            #   
            ch_data = ch_data * norm
            ch_data[ch_data > 255.] = 255.
            fr_data[:,:,i] = ch_data.astype(np.uint8)
        del ch_data
        return fr_data
    #
    # called from _iter_batches()
    def _iter_tiles(self):
        chname = self.config.NAME
        tilesize = self.config.TILESIZE
        for cframe in range(self.tifstk.n_frames):
            fr_data = self._prepare_rgb_frame(cframe)
            for x0, x1, y0, y1 in self.tiles:
                image = fr_data[y0:y1+1, x0:x1+1]
                yield cframe, (x0, x1, y0, y1), image
        for j in range(self.config.BATCH_SIZE):
            image = np.zeros(shape=(tilesize, tilesize, 3), dtype=np.uint8)
            yield cframe, (0, tilesize-1, 0, tilesize-1), image
    #
    # called from segment()
    def _iter_batches(self):
        cframes = []
        ctiles = []
        cimages = []
        for cframe, ctile, cimage in self._iter_tiles():
            cframes.append(cframe)
            ctiles.append(ctile)
            cimages.append(cimage)
            if len(cimages) == self.config.BATCH_SIZE:
                results = self.model.detect(cimages, verbose=0)
                for _cframe, _ctile, res in zip(cframes, ctiles, results):
                    res['frame'] = _cframe
                    res['tile'] = _ctile
                    yield res
                cframes = []
                ctiles = []
                cimages = []
    #
    # --- for debug only ---
    def segment_one(self, tifstk, ifr):
        self.tifstk = tifstk
        
        for i, chname in enumerate(self.RGB_FRAMES):
            otsu = tifstk.getChannelOtsu(chname)
            self.sc_norm[chname] = 63./otsu
        
        chname = self.config.NAME
        tilesize = self.config.TILESIZE
        self.tiles = slice_area((tifstk.height, tifstk.width), (tilesize, tilesize))
        
        scores = []
        particles = []
        fr_data = self._prepare_rgb_frame(ifr)
        for x0, x1, y0, y1 in self.tiles:
            image = fr_data[y0:y1+1, x0:x1+1]
            results = self.model.detect([image], verbose=0)[0]
            tile_mask = assemble_mask_stack(results.pop('masks'))
            if tile_mask is None: continue
            scores.extend([float(x) for x in results['scores']])
            particles.extend(rpeutil.mask_rcnn_to_particles(tile_mask, results['rois'], x0, y0))
            del tile_mask

        fr_mask = np.zeros(shape=(tifstk.height, tifstk.width), dtype=np.uint8)
        flat_scores = rpeutil.recombine_flat_particles(fr_mask, particles, scores, rpeutil.FIX_BORDERS_ACTIN)
        return fr_mask
    #
    def segment(self, tifstk):
        all_start_ts = datetime.datetime.now()
        self.tifstk = tifstk
        
        for i, chname in enumerate(self.RGB_FRAMES):
            otsu = tifstk.getChannelOtsu(chname)
            self.sc_norm[chname] = 63./otsu
        
        chname = self.config.NAME
        tilesize = self.config.TILESIZE
        if self.config.BORDER_ADJUST == rpeutil.FIX_BORDERS_DNA:
            otsu3 = tifstk.getChannelOtsu3(chname)
            self.glob_otsu = int(otsu3[1])
        else:
            self.glob_otsu = None
            
        self.tiles = slice_area((tifstk.height, tifstk.width), (tilesize, tilesize))
        self.mask_cache = {}
    
        bn = tifstk.base_name+'_'+chname+'_RPE'
        ml_dir = tifstk.cache.makedir(self.config.PREDICTED_SUBDIR)
        ml_csv = os.path.join(ml_dir, bn+'.csv')
        ml_tif = os.path.join(ml_dir, bn+'.tif')
        ml_raw_tif = os.path.join(ml_dir, bn+'.raw.tif')
        ml_score_csv = os.path.join(ml_dir, bn+'.score.csv')
        if self.config.WRITE_RAW:
            fo = open(ml_score_csv, 'w')
            wr = csv.writer(fo, dialect='excel')
            wr.writerow(['Frame', 'XC', 'YC', 'Score'])
        #
        NT = len(self.tiles)
        tileg = self._iter_batches()
        actin_frames = np.empty(shape=(tifstk.n_frames, tifstk.height, tifstk.width), dtype=np.uint8)
        #
        #print ('Segmenting %s of %s' % (chname, tifstk.base_name))
        for ifr in range(tifstk.n_frames):
            start_ts = datetime.datetime.now()
            
            scores = []
            particles = []
            for itl in range(NT):
                results = next(tileg)
                x0, x1, y0, y1 = results['tile']
                tile_mask = assemble_mask_stack(results.pop('masks'))
                if tile_mask is None: continue
                scores.extend([float(x) for x in results['scores']])
                particles.extend(rpeutil.mask_rcnn_to_particles(tile_mask, results['rois'], x0, y0))
                del tile_mask
    
            if self.config.BORDER_ADJUST == rpeutil.FIX_BORDERS_DNA:
                fr_mask = self.mask_cache.pop(ifr)
            else:
                fr_mask = np.zeros(shape=(tifstk.height, tifstk.width), dtype=np.uint8)
            flat_scores = rpeutil.recombine_flat_particles(fr_mask, particles, scores, self.config.BORDER_ADJUST)
            
            actin_frames[ifr] = fr_mask
            del fr_mask
            #
            if self.config.WRITE_RAW:
                n_scores = flat_scores[0]
                j = 1
                for i in range(n_scores):
                    xc = flat_scores[j]
                    yc = flat_scores[j+1]
                    sc = flat_scores[j+2]
                    j += 3
                    wr.writerow([ifr, xc, yc, sc])
            #
            elapsed = str(datetime.datetime.now() - start_ts)
            print ('Frame %d of %d segmented in %s' % ((ifr+1), tifstk.n_frames, elapsed))
        #  
        if self.config.WRITE_RAW:
            print ('Write raw predictions:', ml_raw_tif)
            tifffile.imwrite(ml_raw_tif, actin_frames, compress=6, photometric='minisblack')
            print ('Write scores:', ml_score_csv)
            fo.close()
        #
        print ('Assembling 3D %s of %s' % (chname, tifstk.base_name))
        start_ts = datetime.datetime.now()
        rpesegm.assemble_basic(actin_frames, ml_csv)
        elapsed = str(datetime.datetime.now() - start_ts)
        print('3D Assembly done in', elapsed)
        
        print ('Write:', ml_tif)
        tifffile.imwrite(ml_tif, actin_frames, compress=6, photometric='minisblack')
        
        elapsed = str(datetime.datetime.now() - all_start_ts)
        print ('Segmenting %s of %s done in: %s' % (chname, tifstk.base_name, elapsed))
    #

def get_segmenter(channel, weights_path, args):
    segmenter_map = {
        'DNA': RpeDnaSegmenter,
        'Actin': RpeActinSegmenter,
    }
    sclass = segmenter_map.get(channel)
    if sclass is None:
        return None
    return sclass(weights_path, args)

