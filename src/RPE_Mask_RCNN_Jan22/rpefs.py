import os, sys, datetime
import json
import csv
import enum
import xml.etree.ElementTree as ET

import numpy as np
from tifffile import TiffFile, imwrite, imread
import skimage.filters as ski_filters

from src.RPE_Mask_RCNN_Jan22.rpecontour import *
from src.RPE_Mask_RCNN_Jan22.util import *

class RpeStackPage(object):
    def __init__(self, fpath, shape, dtype):
        self.fpath = fpath
        self.shape = shape
        self.dtype = dtype
        #
        self.tags = {'RpeStackPage':True}
    #
    def asarray(self):
        try:
            fr_data = imread(self.fpath)
            self.shape = fr_data.shape
            self.dtype = fr_data.dtype
            self.tags['valid'] = True
            return fr_data
        except Exception:
            pass
        self.tags['valid'] = False
        try:
            return np.zeros(self.shape, self.dtype)
        except Exception:
            return None
#

class RpeJsonStack(object):
    def __init__(self, fpath, fmeta=None):
        self.fpath = fpath
        #
        self.basedir = os.path.abspath(os.path.dirname(fpath))
        #
        self.n_frames = 0
        self.n_channels = 0
        self.dimorder = ''
        self.shape = None
        self.dtype = None
        self.pages = []
        filelist = []
        try:
            if fmeta is None:
                with open(self.fpath, 'r') as fi:
                    o = json.load(fi)
            else:
                o = fmeta
            self.n_frames = o['n_frames']
            self.n_channels = o['n_channels']
            self.dimorder = o['dimorder']
            for rpath in o['filelist']:
                if rpath:
                    fpath = os.path.abspath(os.path.join(self.basedir, rpath))
                else:
                    fpath = None
                filelist.append(fpath)
            self.labels = o['labels']
        except Exception:
            pass
        for fpath in filelist:
            try:
                fr_data = imread(fpath)
                self.shape = fr_data.shape
                self.dtype = fr_data.dtype
                break
            except Exception:
                pass
        else:
            self.nframes = self.nchannels = 0
            return
        self.pages = [RpeStackPage(fpath, self.shape, self.dtype) for fpath in filelist]
    #

def et_find_elem(root, path):
    if not path:
        return root
    nm = path[0]
    for elem in root:
        tag = elem.tag.split('}')[-1]
        if tag == nm:
            return et_find_elem(elem, path[1:])
    return None

class RpeTiffStack(object):
    def __init__(self, fpath, default_channel='DNA'):
        if isinstance(fpath, dict):
            self.fmeta = fmeta = fpath
            self.fpath = os.path.join(fmeta.pop('basedir'), fmeta.pop('basename')+'.rpe.json')
        else:
            self.fmeta = fmeta = None
            self.fpath = fpath
        self.default_channel = default_channel
        #
        self.base_dir, self.name = os.path.split(self.fpath)
        self.base_name, self.ext = os.path.splitext(self.name)
        if self.base_name.endswith('.ome') or self.base_name.endswith('.rpe'):
            self.base_name = self.base_name[:-4]
        #
        if self.ext.lower() == '.json':
            self.tif = RpeJsonStack(self.fpath, self.fmeta)
        else:
            self.tif = TiffFile(self.fpath)
        self.n_pages = len(self.tif.pages)
        self.n_frames = 0
        self.n_channels = 0
        self.dimorder = ''
        self.height = 0
        self.width = 0
        #
        self.shape = None
        self.dtype = None
        self.tags = {}
        if self.n_pages > 0:
            # self.nframes = self.n_pages / 4
            page = self.tif.pages[0]
            self.tags = page.tags
            self.shape = page.shape
            self.height = self.shape[0]
            self.width = self.shape[1]
            self.dtype = page.dtype
            self._parse_ome()
        #
        if self.n_frames > 0:
            if self.n_channels == 4:
                if hasattr(self.tif, 'labels') and len(self.tif.labels) == 4:
                    self.channels = self.tif.labels
                else:
                    self.channels = ['GFP', 'DNA', 'Actin', 'Membrane',]
            elif self.n_channels == 1:
                self.channels = [self.default_channel,]
        #
        self.framecache = {}
        self.cache = RpeCache(self.base_dir, self.base_name, self.width, self.height, self.n_frames)
    #
    def _parse_ome(self):
        try:
            tag = self.tags['RpeStackPage']
            self.dimorder = self.tif.dimorder
            self.n_channels = self.tif.n_channels
            self.n_frames = self.tif.n_frames
            return
        except Exception:
            pass
        try:
            tag = self.tags['ImageDescription']
            elem = et_find_elem(ET.fromstring(tag.value), ['Image', 'Pixels'])
            self.dimorder = elem.get('DimensionOrder')
            self.n_channels = int(elem.get('SizeC'))
            self.n_frames = int(elem.get('SizeZ'))
        except Exception:
            self.dimorder = ''
            self.n_channels = 1
            self.n_frames = self.n_pages
    #
    def saveMeta(self):
        if self.fmeta is None:
            return None
        try:
            with open(self.fpath, 'w') as fo:
                json.dump(self.fmeta, fo, indent=2)
            return self.fpath
        except Exception as ex:
            print (ex)
        return None
    #
    def channelName(self, ichan):
        try:
            chname = self._chlist[ichan]
        except Exception:
            chname = 'ch%02d' % (ichan,)
        return chname
    def channelId(self, chname, validated=False):
        if chname in self._chmap:
            return self._chmap[chname]
        if validated:
            raise ValueError('No such channel: '+chname)
        return -1
    def hasChannel(self, chname):
        return chname in self._chmap
    #
    @property
    def channels(self):
        return [self.channelName(ichan) for ichan in range(self.n_channels)]
    @channels.setter
    def channels(self, ch_list):
        self._chlist = [str(ch) for ch in ch_list]
        # self._chmap = dict([(self.channelName(ichan), ichan) for ichan in range(self.n_channels)])
        self._chmap = {}
        for ichan, ch in enumerate(self._chlist):
            self._chmap[ichan] = ichan
            self._chmap[ch] = ichan
            self._chmap[ch.replace('0', 'O')] = ichan
            self._chmap[ch.replace('O', '0')] = ichan
    #
    def isZ01(self, chname):
        ichan = self.channelId(chname)
        if ichan < 0: return False
        return ichan == self.channelId('Z01')
    #
    def getChannelData(self, chname=None):
        self.framecache.clear()
        if chname is None:
            chname = self.default_channel
        ichan = self.channelId(chname, True)
        datayz = np.empty(shape=(self.n_frames, self.height, self.width), dtype=self.dtype)
        for iframe in range(self.n_frames):
            page = self.tif.pages[iframe*self.n_channels + ichan]
            datayz[iframe] = page.asarray()
        return datayz
    #
    def getChannelOtsu(self, chname=None, cachedonly=False):
        ichan = self.channelId(chname, True)
        chname = self.channelName(ichan)
        otsu = self.cache.get(chname, 'otsu')
        if otsu is None and not cachedonly:
            ch_data = np.float32(self.getChannelData(ichan))
            otsu = float(ski_filters.threshold_otsu(ch_data, 4096))
            self.cache.put(chname, 'otsu', otsu)
            self.cache.save()
        return otsu
    #
    def getChannelOtsu3(self, chname=None, cachedonly=False):
        ichan = self.channelId(chname, True)
        chname = self.channelName(ichan)
        otsu3 = self.cache.get(chname, 'otsu3')
        if otsu3 is None and not cachedonly:
            ch_data = np.float32(self.getChannelData(ichan))
            motsu = ski_filters.threshold_multiotsu(ch_data, classes=3, nbins=512)
            otsu3 = [float(motsu[0]), float(motsu[1])]
            self.cache.put(chname, 'otsu3', otsu3)
            self.cache.save()
        return otsu3
    #
    def getZ01FlatOtsu(self, cachedonly=False):
        chname = self.cache.z01flatCh
        otsu = self.cache.get(chname, 'otsu')
        if otsu is None and not cachedonly:
            ch_data = self.cache.getZ01flat()
            if not ch_data is None:
                otsu = float(ski_filters.threshold_otsu(ch_data, 4096))
                self.cache.put(chname, 'otsu', otsu)
                self.cache.save()
        return otsu       
    #
    def getFrame(self, iframe, chname=None, nocache=False):
        if chname is None:
            chname = self.default_channel
        ichan = self.channelId(chname, True)
        fr_key = (ichan, iframe)
        if fr_key in self.framecache:
            return self.framecache[fr_key]
        page = self.tif.pages[iframe*self.n_channels + ichan]
        data = page.asarray()
        if not nocache:
            self.framecache[fr_key] = data
        return data
    #
    def needToLoadSegm(self, iframe, chname=None):
        if chname is None:
            chname = self.default_channel
        if self.isZ01(chname):
            return not self.cache.isZ01Cached()
        ichan = self.channelId(chname, True)
        chname = self.channelName(ichan)
        return not self.cache.isChannelCached(chname)
    #
    def getFrameAndSegmentation(self, iframe, chname=None):
        if chname is None:
            chname = self.default_channel
        if self.isZ01(chname):
            key = 'Z01flat'
            if key in self.framecache:
                fr_data = self.framecache[key]
            else:
                fr_data = self.cache.getZ01flat()
                if fr_data is None:
                    fr_data = self.getFrame(iframe, chname)
                else:
                    self.framecache[key] = fr_data
            return fr_data, self.cache.getZ01Contours()
        ichan = self.channelId(chname, True)
        chname = self.channelName(ichan)
        return self.getFrame(iframe, chname), self.cache.getFrameContours(chname, iframe)
    #
    def getChannelManager(self, chname=None):
        if chname is None:
            chname = self.default_channel
        if self.isZ01(chname):
            return self.cache.getZ01Manager()
        ichan = self.channelId(chname)
        if ichan < 0:
            return None
        chname = self.channelName(ichan)
        return self.cache.getManager(chname)
    def isChannelEmpty(self, chname=None):
        mgr = self.getChannelManager(chname)
        if mgr is None:
            return True
        return mgr.empty
    #
    def eraseChannel(self, chname=None):
        mgr = self.getChannelManager(chname)
        if mgr is None:
            return False
        return mgr.deleteAll()
    #
    @staticmethod
    def is_acceptable_name(fn):
        base, ext = os.path.splitext(fn.lower())
        if ext in ('.tif', '.tiff'):
            return True
        if not ext in ('.json',):
            return False
        return base.endswith('.rpe')
    #
    @staticmethod
    def channelFromFilename(fpath):
        try:
            fn, ext = os.path.splitext(os.path.basename(fpath))
            parts = fn.split('_')
            if parts[-1] == 'RPE':
                return parts[-2]
        except Exception:
            pass
        return None
#

class RpeCache(object):
    SEGMENTED = 'Segmented'
    MANUAL = 'Manual'
    RESHAPE = 'REShAPE'
    CACHE = '.cache'
    Z01_synonyms = set(['ZO1', 'Z01', 'zo1', 'z01'])
    def __init__(self, base_dir, base_name, w, h, d):
        self.base_dir = base_dir
        self.base_name = base_name
        self.w = w
        self.h = h
        self.d = d
        #
        self.cache_rel = self.CACHE
        self.cache_dir = os.path.join(self.base_dir, self.cache_rel)
        self.cache_fpath = os.path.join(self.cache_dir, self.base_name+'.cache.json')
        #
        self.segm_dir = os.path.join(self.base_dir, self.SEGMENTED)
        self.reshape_dir = os.path.join(self.base_dir, self.RESHAPE)
        self.manual_dir = os.path.join(self.base_dir, self.MANUAL)
        #
        self.z01flatCh = 'Z01flat'
        self.z01flat = None
        self.cache = {}
        self.mgrmap = {}
        self.load()
    #
    def makedir(self, relpath):
        fullpath = os.path.join(self.base_dir, relpath)
        if not os.path.isdir(fullpath):
            try:
                os.makedirs(fullpath)
            except Exception as ex:
                print ('RpeCache.makedir(', relpath, '):', str(ex))
                return None
        return fullpath
    #
    def save(self, contours=True):
        cache_dir = self.makedir(self.cache_rel)
        if cache_dir is None:
            return False
        try:
            if contours:
                for chname, mgr in self.mgrmap.items():
                    if hasattr(mgr, 'save'):
                        # print (chname, 'pickle save:', mgr.cache_fpath)
                        mgr.save()
            with open(self.cache_fpath, 'w') as fo:
                json.dump(self.cache, fo)
                return True
        except Exception:
            pass
        return False
    #
    def load(self):
        try:
            with open(self.cache_fpath, 'r') as fi:
                self.cache = json.load(fi)
                return True
        except Exception:
            pass
        self.cache = {}
        return False
    #
    def get(self, chname, key, dflt=None):
        try:
            return self.cache[chname][key]
        except Exception:
            return dflt
    #
    def put(self, chname, key, value):
        if not chname in self.cache:
            self.cache[chname] = {}
        self.cache[chname][key] = value
    #
    def getCacheName(self, chname, ext='.pkl'):
        return os.path.join(self.cache_dir, self.base_name+'_'+chname+'_RPE'+ext)
    #
    def getSegmName(self, chname, ext='.csv'):
        return os.path.join(self.segm_dir, self.base_name+'_'+chname+'_RPE'+ext)
    def getManualName(self, chname, ext='.csv'):
        return os.path.join(self.manual_dir, self.base_name+'_'+chname+'_RPE'+ext)
    def isChannelCached(self, chname):
        return chname in self.mgrmap
    def getManager(self, chname):
        if not chname in self.mgrmap:
            cache_fpath = self.getCacheName(chname)
            mgr = ContourManager3d(self.w, self.h, self.d, cache_fpath)
            self.mgrmap[chname] = mgr
        else:
            mgr = self.mgrmap[chname]
        return mgr
    def getFrameContours(self, chname, iframe):
        mgr = self.getManager(chname)
        return mgr.getContourList(iframe)
    #
    def loadSegmentation(self, chname, csv_file):
        mgr = self.getManager(chname)
        mgr.loadSegmentedCells(csv_file)
        return len(mgr.cells)
    def loadPickled(self, chname, pkl_file):
        mgr = self.getManager(chname)
        mgr.loadPickled(pkl_file)
        return len(mgr.cells)
    #
    def saveSegmentation(self, chname, csv_file, tif_file, validate=True, separate=True):
        mgr = self.getManager(chname)
        if tif_file is None:
            # If no tif given, csv_file must be actually a pickle file name
            mgr.savePickled(csv_file)
        else:
            mgr.saveSegmented(csv_file, tif_file, validate, separate)
    #
    def deleteAll(self, chname):
        mgr = self.getManager(chname)
        if mgr.empty:
            return False
        mgr.deleteAll()
        return True
    #
    def getZ01flatName(self):
        return os.path.join(self.reshape_dir, self.base_name+'_Z01_REShAPE.tif')
    def getZ01flat(self):
        if self.z01flat is None:
            fpath = self.getZ01flatName()
            if os.path.isfile(fpath):
                self.z01flat = imread(fpath)
        return self.z01flat
    #
    def getZ01SegmName(self, ext='.csv'):
        return os.path.join(self.segm_dir, self.base_name+'_Z01_RPE'+ext)
    def getZ01ManualName(self, ext='.csv'):
        return os.path.join(self.manual_dir, self.base_name+'_Z01_RPE'+ext)
    def isZ01Cached(self):
        return self.z01flatCh in self.mgrmap
    def getZ01Manager(self):
        chname = self.z01flatCh
        if not chname in self.mgrmap:
            cache_fpath = self.getCacheName(chname)
            mgr = ContourManager2d(self.w, self.h, cache_fpath)
            self.mgrmap[chname] = mgr
        else:
            mgr = self.mgrmap[chname]
        return mgr
    def getZ01Contours(self):
        mgr = self.getZ01Manager()
        return mgr.getContourList()
    #
    def loadZ01Segmentation(self, csv_file):
        mgr = self.getZ01Manager()
        mgr.loadSegmentedContours(csv_file)
        return len(mgr.contours)
    def loadZ01Pickled(self, pkl_file):
        mgr = self.getZ01Manager()
        mgr.loadPickled(pkl_file)
        return len(mgr.contours)
    #
    def saveZ01Segmentation(self, csv_file, tif_file):
        mgr = self.getZ01Manager()
        if tif_file is None:
            mgr.savePickled(csv_file)
        else:
            mgr.saveSegmented(csv_file, tif_file)
    #
    def loadReshapeSegmentation(self, csv_file):
        csv_dir, csv_fn = os.path.split(csv_file)
        tif_file = os.path.join(os.path.dirname(csv_dir), 'Segmented Images', csv_fn[:-9]+'_Outlines.tif')
        mgr = self.getZ01Manager()
        mgr.loadReshapeContours(tif_file, csv_file)
        return len(mgr.contours)
    #
    def deleteZ01All(self):
        mgr = self.getZ01Manager()
        if mgr.empty:
            return False
        mgr.deleteAll()
        return True
    #
#

class RpeDataDir(object):
    def __init__(self, fpath, default_channel='DNA'):
        self.fpath = os.path.abspath(fpath)
        self.default_channel = default_channel
        #
        if os.path.isdir(fpath):
            self.dirpath = fpath
            self.flist = []
            for fn in os.listdir(self.dirpath):
                if RpeTiffStack.is_acceptable_name(fn):
                    self.flist.append(fn)
        elif os.path.exists(fpath):
            self.dirpath, fn = os.path.split(fpath)
            if not RpeTiffStack.is_acceptable_name(fn):
                raise ValueError('Not recognized file. ' +
                    'Acceptable extensions are ".tif", ".tiff", ".rpe.json"')
            self.flist = [fn,]
        else:
            raise ValueError('No such file or directory: '+fpath)
    #
    def isEmpty(self):
        return len(self.flist) == 0
    #
    def getSubdir(self, subdir):
        spath = os.path.join(self.dirpath, subdir)
        if not os.path.isdir(spath):
            os.mkdir(spath)
        return spath
    #
    def iterstacks(self):
        for i, fn in enumerate(self.flist):
            fpath = os.path.join(self.dirpath, fn)
            item = '%d of %d' % ((i+1), len(self.flist))
            yield item, RpeTiffStack(fpath, self.default_channel)
#

@enum.unique
class CsvType(enum.IntEnum):
    Unknown = 0
    Flat = 1
    Multi = 2
    Reshape = 3

def getCsvType(csv_fpath):
    try:
        csv_dir, csv_fn = os.path.split(csv_fpath)
        with open(csv_fpath, 'r') as fi:
            rd = csv.reader(fi, dialect='excel')
            hdrs = set(next(rd))
            for hd in ['ID', 'y', 'xL', 'xR']:
                if not hd in hdrs:
                    break
            else:
                return CsvType.Multi if 'Frame' in hdrs else CsvType.Flat
            # Try REShAPE
            if not csv_fn.endswith('_Data.csv'):
                return CsvType.Unknown
            tif_fpath = os.path.join(csv_dir, '..', 'Segmented Images', csv_fn[:-9]+'_Outlines.tif')
            if not os.path.isfile(tif_fpath):
                return CsvType.Unknown
            for hd in ['XStart', 'YStart']:
                if not hd in hdrs:
                    return CsvType.Unknown
            return CsvType.Reshape
    except Exception:
        return CsvType.Unknown
#

if __name__ == '__main__':
    cdir = '.'
    if len(sys.argv) > 1:
        cdir = sys.argv[1]
    rd = RpeDataDir(cdir)
    print ('Found %d files in %s' % (len(rd.flist), rd.fpath))
    for _, stk in rd.iterstacks():
        print (stk.base_name, ':', stk.channels)
        break
        
    