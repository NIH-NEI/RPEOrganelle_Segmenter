__all__ = ('isPointInside', 'isIntersected', 'findContour', 'deleteMultiple', 'findIntersecting',
           'splitContours', 'joinContours', 'dist', 'optimizeContour', 'sameContour', 'contourChanged',
           'importContours', 'importReshapeContours',
           'contourBoundary', 'flattenContour', 'ContourManager2d', 'ContourManager3d',)

import datetime
import math
import pickle
from collections import defaultdict

import numpy as np
from tifffile import imread, imwrite

from src.RPE_Mask_RCNN_Jan22 import rpeutil
from src.RPE_Mask_RCNN_Jan22.util import *

# Winding number test for a point in a polygon
# Adapted from: http://geomalgorithms.com/a03-_inclusion.html

# isLeft(): tests if a point is Left|On|Right of an infinite line.
#    Input:  three points P0, P1, and P2
#    Return: >0 for P2 left of the line through P0 and P1
#            =0 for P2  on the line
#            <0 for P2  right of the line
#
# P1, P2, P3 are lists [x,y] or tuples (x,y)
def isLeft(P0, P1, P2):
    return (P1[0]-P0[0])*(P2[1]-P0[1]) - (P2[0]-P0[0])*(P1[1]-P0[1])

# wn_PnPoly(): winding number test for a point in a polygon
#      Input:   pt = a point, list or tuple (x,y)
#               contour = vertex points of a polygon, a collection of points
#                  [(x0,y0), (x1,y1), ...]
#      Return:  wn = the winding number (=0 only when pt is outside)
def wn_PnPoly(pt, contour):
    poly = list(contour)
    n = len(poly)
    poly.append(poly[0])    # make poly[0] == poly[n+1]
    wn = 0
    for i in range(n):
        if poly[i][1] < pt[1]:          # start y <= pt.y
            if poly[i+1][1] > pt[1]:    # an upward crossing
                if isLeft(poly[i], poly[i+1], pt) > 0:  # pt left of edge
                    wn += 1             # have a valid up intersect
        else:                           # start y > P.y (no test needed)
            if poly[i+1][1] <= pt[1]:   # a downward crossing
                if isLeft(poly[i], poly[i+1], pt) < 0:  # pt right of edge
                    wn -= 1             # have a valid down intersect
    return wn

def isPointInside(pt, contour):
    return wn_PnPoly(pt, contour) != 0

# Simplified -- looks for the edges of one polygon inside another
def isIntersected(contour1, contour2):
    for pt in contour1:
        if wn_PnPoly(pt, contour2) != 0:
            return True
#     for pt in contour2:
#         if wn_PnPoly(pt, contour1) != 0:
#             return True
    return False

# Find contour containing point pt
def findContour(pt, contours):
    if isinstance(contours, zContourList):
        return contours.find(pt)
    for i, contour in enumerate(contours):
        if wn_PnPoly(pt, contour) != 0:
            return i
    return -1

# Delete multiple indexes from contour list
def deleteMultiple(contours, idxlist):
    if len(idxlist) == 0:
        return
    if hasattr(contours, 'delMultiple'):
        contours.delMultiple(idxlist)
    else:
        upd = [cont for idx, cont in enumerate(contours) if not idx in idxlist]
        contours[:] = upd
#

# Find contours intersecting with cont, return list of their indexes
def findIntersecting(contours, cont):
    if hasattr(contours, 'findIntersecting'):
        return contours.findIntersecting(cont)
    return [idx for idx, _cont in enumerate(contours) if isIntersected(_cont, cont)]

# Find contours split by a line drawn by user
def splitContours(contours, cont):
    if not hasattr(contours, '_mgr'):
        idxlist = []
        w = h = 0
    else:
        idxlist = findIntersecting(contours, cont)
        w = contours._mgr.w
        h = contours._mgr.h
    if len(idxlist) == 0:
        return [], []
    flat_slices = []
    for idx in idxlist:
        flat_slices.extend(flattenContour(contours[idx]))
    res = rpeutil.cut_slices(w, h, flat_slices, flattenContour(cont))
    if len(res) > 0:
        dlen = res[0]
        dellist = [idxlist[i] for i in res[1:1+dlen]]
        addlist = _unflatten_contours(res[1+dlen:])
        return dellist, addlist
    return [], []

# Join contours intersecting with the contour drawn by the user
def joinContours(contours, cont):
    if not hasattr(contours, '_mgr'):
        idxlist = []
        w = h = 0
    else:
        idxlist = findIntersecting(contours, cont)
        w = contours._mgr.w
        h = contours._mgr.h
    if len(idxlist) == 0:
        return [], []
    flat_slices = []
    for idx in idxlist:
        flat_slices.extend(flattenContour(contours[idx]))
    res = rpeutil.join_slices(w, h, flat_slices, flattenContour(cont))
    if len(res) > 0:
        dlen = res[0]
        dellist = [idxlist[i] for i in res[1:1+dlen]]
        addlist = _unflatten_contours(res[1+dlen:])
        return dellist, addlist
    return [], []

# Geometric distance between two points
def dist(pt1, pt2):
    return math.sqrt((pt1[0]-pt2[0])**2 + (pt1[1]-pt2[1])**2)

# Optimize a contour by removing vertices too close to each other
def optimizeContour_0(contour, min_dist=1.5):
    n = len(contour)
    if n < 5:
        return contour
    pt1 = contour[0]
    res = [pt1]
    for i in range(1, n):
        pt2 = contour[i]
        cur_dist = dist(pt1, pt2)
        if cur_dist >= min_dist:
            res.append(pt2)
            pt1 = pt2
    return res

# Optimize a contour by removing vertices in the middle of straight lines
def optimizeContour(contour):
    res = rpeutil.simplify_contour(flattenContour(contour))
    return _unflatten_contours(res)[0]

def sameContour(cont1, cont2):
    if len(cont1) != len(cont2):
        return False
    for pt1, pt2 in zip(cont1, cont2):
        if pt1[0] != pt2[0] or pt1[1] != pt2[1]:
            return False
    return True

# Test if a contour is sufficiently changed (after an edit operation)
def contourChanged(oldc, newc, tolerance=0.005):
    if len(oldc) != len(newc):
        return True
    for pt1, pt2 in zip(oldc, newc):
        if math.fabs(pt1[0] - pt2[0]) > tolerance:
            return True
        if math.fabs(pt1[1] - pt2[1]) > tolerance:
            return True
    return False

# Calculate bounding rectangle
def contourBoundary(cont):
    xmin = xmax = ymin = ymax = None
    for x, y in cont:
        if xmin is None:
            xmin = xmax = x
            ymin = ymax = y
        else:
            if x < xmin: xmin = x
            if x > xmax: xmax = x
            if y < ymin: ymin = y
            if y > ymax: ymax = y
    return xmin, xmax, ymin, ymax
#

# Convert contour to flat int array (for passing to native code)
def flattenContour(cont):
    res = [len(cont),]
    for x, y in cont:
        res.append(int(x+0.5))
        res.append(int(y+0.5))
    return res
#

def _unflatten_contours(flat):
    res = []
    i = 0
    while i < len(flat):
        sz = flat[i];
        i += 1
        cont = []
        for n in range(sz):
            cont.append((flat[i], flat[i+1]))
            i += 2
        res.append(cont)
    return res

def importContours(w, h, csv_file):
    return _unflatten_contours(rpeutil.import_contours(w, h, csv_file))

def importReshapeContours(rs_mask, rs_csv):
    return _unflatten_contours(rpeutil.import_reshape_contours(rs_mask, rs_csv))

# 2D segmentation

class ContourList(object):
    def __init__(self, _mgr):
        self._mgr = _mgr
        #
        self._undo = _mgr._undo
    #
    @property
    def contours(self):
        return self._mgr.contours[:]
    #
    @contours.setter
    def contours(self, upd_contours):
        self._mgr.contours[:] = upd_contours
        self._mgr.dirty = True
    #
    #
    def __len__(self):
        return len(self.contours)
    #
    def __getitem__(self, idx):
        if isinstance(idx, slice):
            raise KeyError('Slicing not supported')
        return self.contours[idx]
    #
    def __setitem__(self, idx, val):
        if isinstance(idx, slice):
            raise KeyError('Slicing not supported')
        upd_contours = self.contours
        upd_contours[idx] = self._mgr.validateContour(val, idx)
        self.contours = upd_contours
    #
    def __delitem__(self, idx):
        upd_contours = self.contours
        del upd_contours[idx]
        self.contours = upd_contours
    #
    def __iter__(self):
        return self.contours.__iter__()
    #
    def append(self, val):
        idx = self._mgr.addContour(val)
        if idx >= 0:
            self._mgr.dirty = True
    #
    def delMultiple(self, idxlist):
        upd_contours = [item for idx, item in enumerate(self.contours) if not idx in idxlist]
        self.contours = upd_contours
    #
    def find(self, pt):
        for idx, _cont in enumerate(self.contours):
            if wn_PnPoly(pt, _cont) != 0:
                return idx
        return -1
    #
    def findIntersecting(self, cont):
        flat_cont = flattenContour(cont)
        flat_slices = []
        for _cont in self.contours:
            flat_slices.extend(flattenContour(_cont))
        return rpeutil.intersecting_slices(flat_slices, flat_cont)
        #return [idx for idx, (_, _, _cont) in enumerate(self._ided) if isIntersected(cont, _cont)]
    #
    #
    def indexOf(self, cont):
        for idx, _cont in enumerate(self._mgr.contours):
            if sameContour(cont, _cont):
                return idx
        return -1
    #
    def undo_empty(self):
        return self._undo.is_empty()
    #
    def push_undo(self, to_remove, to_add):
        return self._undo.push_undo(to_remove, to_add)
    #
    def undo(self):
        to_remove, to_add = self._undo.pop_undo()
        rc = False
        if not to_remove is None:
            idxlist = []
            for cont in to_remove:
                idx = self.indexOf(cont)
                if idx >= 0:
                    idxlist.append(idx)
            if len(idxlist) > 0:
                self.delMultiple(idxlist)
                rc = True
        if not to_add is None:
            for cont in to_add:
                self.append(cont)
            rc = True
        return rc
    #


class ContourManager2d(object):
    def __init__(self, w, h, cache_fpath):
        self.w = w
        self.h = h
        self.cache_fpath = cache_fpath
        #
        self.contours = []
        self._undo = UndoStack()
        #
        self._dirty = False
        self.load()
    #
    @property
    def dirty(self):
        return self._dirty
    @dirty.setter
    def dirty(self, st):
        self._dirty = st
    #
    def loadSegmentedContours(self, csv_file):
        self._undo.clear()
        flat = rpeutil.import_contours(self.w, self.h, csv_file)
        self.contours = _unflatten_contours(flat)
        self._dirty = True
        self.save()
    #
    def loadReshapeContours(self, rs_tif, rs_csv):
        self._undo.clear()
        rs_mask = imread(rs_tif)
        self.contours = _unflatten_contours(rpeutil.import_reshape_contours(rs_mask, rs_csv))
        self._dirty = True
        self.save()
    #
    @property
    def empty(self):
        return len(self.contours) == 0
    #
    def deleteAll(self):
        self._undo.clear()
        self.contours = []
        self._dirty = True
        self.save()
    #
    def countContours(self):
        return len(self.contours)
    #
    def getContourList(self, z=None):
        return ContourList(self)
    #
    def getAnnotations(self):
        self.save()
        res = []
        for idx, cont in enumerate(self.contours):
            id = idx + 1
            all_x = [x for x, _ in cont]
            all_y = [y for _, y in cont]
            res.append({
                "shape_attributes" : {
                    "name" : "polygon",
                    "all_points_x" : all_x,
                    "all_points_y" : all_y},
                "region_attributes" : {
                    "cell" : id,
                    },}
                )
        return [res]
    def addContour(self, cont):
        if len(cont) >= 5:
            _idx = len(self.contours)
            cont = self.validateContour(cont, -1)
            if len(cont) >= 5:
                self.contours.append(cont)
        else:
            _idx = -1
        return _idx
    #
    def validateContour(self, cont, idx=-1):
        flat_slices = []
        for _idx, _cont in enumerate(self.contours):
            if _idx != idx:
                flat_slices.extend(flattenContour(_cont))
        flat_cont = flattenContour(cont)
        res = rpeutil.validate_contour(self.w, self.h, flat_slices, flat_cont)
        return _unflatten_contours(res)[0]
    #
    def get_flat_contours(self):
        flat_slices = []
        for cont in self.contours:
            flat_slices.extend(flattenContour(cont))
        return flat_slices
    def getEmptyMask(self):
        return np.empty(shape=(self.h, self.w), dtype=np.uint8)
    #
    def _filter_contours_at_border(self):
        flat_slices = self.get_flat_contours()
        at_bord_list = rpeutil.contours_at_border(self.w, self.h, flat_slices)
        if len(at_bord_list) == 0:
            return
        at_bord = set(at_bord_list)
        upd = [cont for idx, cont in enumerate(self.contours) if not idx in at_bord]
        # print ('At border:', len(at_bord_list), ' of ', len(self.contours))
        self.contours = upd
    #
    def saveSegmented(self, csv_file, tif_file):
        ch_mask = self.getEmptyMask()
        flat_slices = self.get_flat_contours()
        rpeutil.export_2d_contours(ch_mask, flat_slices, csv_file)
        print ('Write:', tif_file)
        imwrite(tif_file, ch_mask, compress=6, photometric='minisblack')
    #
    def savePickled(self, pkl_file):
        with open(pkl_file, 'wb') as fo:
            pickle.dump(self.contours, fo, pickle.HIGHEST_PROTOCOL)
    #
    def loadPickled(self, pkl_file):
        self._undo.clear()
        with open(pkl_file, 'rb') as fi:
            contours = pickle.load(fi)
        try:
            for cont in contours:
                for x, y in cont:
                    pass
        except Exception:
            raise RuntimeError('Invalid or incompatible file format')
        self.contours = contours
        self._filter_contours_at_border()
        self._dirty = True
        self.save()
    #
    def load(self):
        try:
            with open(self.cache_fpath, 'rb') as fi:
                self.contours = pickle.load(fi)
            self._filter_contours_at_border()
        except Exception:
            pass
    #
    def save(self):
        if self._dirty:
            try:
                self.savePickled(self.cache_fpath)
            except Exception:
                pass
            self._dirty = False
    #



# 3D segmentation

class aCell(object):
    def __init__(self, id, nfr):
        self.id = id
        self.nfr = nfr
        #
        self.zmap = tuple([[] for z in range(self.nfr)])
        self.bnd_valid = False
        self.xmin = self.xmax = self.ymin = self.ymax = self.zmin = self.zmax = None
    #
    def _calc_boundary(self):
        self.xmin = self.xmax = self.ymin = self.ymax = self.zmin = self.zmax = None
        for z in range(self.nfr):
            if len(self.zmap[z]) == 0: continue
            if self.zmin is None:
                self.zmin = z
            self.zmax = z
            for cont in self.zmap[z]:
                for x, y in cont:
                    if self.xmin is None:
                        self.xmin = self.xmax = x
                        self.ymin = self.ymax = y
                    else:
                        if x < self.xmin: self.xmin = x
                        if x > self.xmax: self.xmax = x
                        if y < self.ymin: self.ymin = y
                        if y > self.ymax: self.ymax = y
        self.bnd_valid = not self.zmin is None
    #
    def intersection(self, bnd):
        if not self.bnd_valid:
            self._calc_boundary()
            if not self.bnd_valid:
                return 0.
        _zmin = _zmax = None
        if len(bnd) == 4:
            _xmin, _xmax, _ymin, _ymax = bnd
        elif len(bnd) == 6:
            _xmin, _xmax, _ymin, _ymax, _zmin, _zmax = bnd
        else:
            raise ValueError('parameter must be of size 4 or 6')
        xl = max(_xmin, self.xmin)
        xr = min(_xmax, self.xmax)
        if xl > xr:
            return 0.
        yt = max(_ymin, self.ymin)
        yb = min(_ymax, self.ymax)
        if yt > yb:
            return 0.
        if _zmin is None:
            return (xr-xl+1.)*(yb-yt+1.)
        zs = max(_zmin, self.zmin)
        ze = min(_zmax, self.zmax)
        if zs > ze:
            return 0.
        return (xr-xl+1.)*(yb-yt+1.)*(ze-zs+1.)
    #
    def count(self):
        return sum(len(self.zmap[z]) for z in range(self.nfr))
    #
    def validate(self):
        nc = 0
        for z in range(self.nfr):
            if len(self.zmap[z]) == 0: continue
            self.zmap[z][:] = [cont for cont in self.zmap[z] if cont]
            nc += len(self.zmap[z])
        self.bnd_valid = False
        return nc > 0
    #
    def getContours(self, z):
        return [cont[:] for cont in self.zmap[z]]
    def setContours(self, z, contours):
        self.zmap[z][:] = [cont[:] for cont in contours]
        self.bnd_valid = False
    def addContour(self, z, cont):
        self.zmap[z].append(cont[:])
        self.bnd_valid = False
        return len(self.zmap[z]) - 1
    def delContourAt(self, z, idx):
        del self.zmap[z][idx]
        self.bnd_valid = False
    def replaceContourAt(self, z, idx, cont):
        self.zmap[z][idx] = cont[:]
        self.bnd_valid = False
    def deleteContours(self, z):
        self.zmap[z][:] = []
        self.bnd_valid = False
    def listContours(self, z):
        for idx, cont in enumerate(self.zmap[z]):
            if not cont is None:
                yield (self.id, idx, cont[:])
    #
    def moveContoursTo(self, other, idx_list):
        for z, idx in idx_list:
            other.addContour(z, self.zmap[z][idx])
            self.zmap[z][idx] = None
        for z in range(self.nfr):
            if len(self.zmap[z]) == 0: continue
            self.zmap[z][:] = [cont for cont in self.zmap[z] if not cont is None]
        self.bnd_valid = False
        other.bnd_valid = False
    #
    def join(self, other):
        for z in range(self.nfr):
            if len(other.zmap[z]) == 0: continue
            self.zmap[z].extend(other.zmap[z])
            other.zmap[z][:] = []
        self.bnd_valid = False
    #
    def flat(self):
        res = []
        for z in range(self.nfr):
            for cont in self.zmap[z]:
                if cont is None:
                    cont = []
                res.append(z)
                res.append(int(self.id))
                res.append(len(cont))
                for x, y in cont:
                    res.append(int(x+0.5))
                    res.append(int(y+0.5))
        return res
    #

class zContourList(object):
    def __init__(self, _mgr, _z):
        self._mgr = _mgr
        self._z = _z
        #
        self._ided = []
        self._undo = _mgr.undo_cache[_z]
        #
        self.reload()
    #
    def reload(self):
        self._ided[:] = []
        for id in sorted(self._mgr.cells.keys()):
            cell = self._mgr.cells[id]
            for id_cont in cell.listContours(self._z):
                self._ided.append(id_cont)
    #
    def update(self, upd_ided):
        for cell in self._mgr.cells.values():
            cell.deleteContours(self._z)
        for _id, _, _cont in upd_ided:
            self._mgr.cells[_id].addContour(self._z, _cont)
        self.reload()
        self._mgr.dirty = True
    #
    def __len__(self):
        return len(self._ided)
    #
    def __getitem__(self, idx):
        if isinstance(idx, slice):
            raise KeyError('Slicing not supported')
        _, _, _cont = self._ided[idx]
        return _cont
    #
    def __setitem__(self, idx, val):
        if isinstance(idx, slice):
            raise KeyError('Slicing not supported')
        _id, _idx, _cont = self._ided[idx]
        _cont[:] = self._mgr.validateContour(val, self._z, _id, _idx)
        if len(_cont) <= 3:
            self._mgr.cells[_id].delContourAt(self._z, _idx)
            self.reload()
        else:
            self._mgr.cells[_id].replaceContourAt(self._z, _idx, _cont)
        self._mgr.dirty = True
    #
    def __delitem__(self, idx):
        if isinstance(idx, slice):
            raise KeyError('Slicing not supported')
        _id, _idx, _ = self._ided[idx]
        self._mgr.cells[_id].delContourAt(self._z, _idx)
        self.reload()
        self._mgr.dirty = True
    #
    def __iter__(self):
        for (_, _, _cont) in self._ided:
            yield _cont
    #
    def append(self, val):
        _cont = val[:]
        _id, _idx, cont = self._mgr.addSlice(self._z, _cont)
        if not _id is None:
            self._mgr.selected_id = _id
            self._ided.append((_id, _idx, cont))
            self._mgr.dirty = True
    #
    def simplify(self, idx, simplifyContour=None):
        _id, _idx, _cont = self._ided[idx]
        _cont[:] = simplifyContour(_cont)
    #
    def delMultiple(self, idxlist):
        upd_ided = [item for idx, item in enumerate(self._ided) if not idx in idxlist]
        self.update(upd_ided)
    #
    def find(self, pt):
        for idx, (_, _, _cont) in enumerate(self._ided):
            if wn_PnPoly(pt, _cont) != 0:
                return idx
        return -1
    #
    def findIntersecting(self, cont):
        flat_cont = flattenContour(cont)
        flat_slices = []
        for _, _, _cont in self._ided:
            flat_slices.extend(flattenContour(_cont))
        return rpeutil.intersecting_slices(flat_slices, flat_cont)
        #return [idx for idx, (_, _, _cont) in enumerate(self._ided) if isIntersected(cont, _cont)]
    #
    def indexOf(self, cont):
        for idx, (_, _, _cont) in enumerate(self._ided):
            if sameContour(cont, _cont):
                return idx
        return -1
    #
    def undo_empty(self):
        return self._undo.is_empty()
    #
    def push_undo(self, to_remove, to_add):
        return self._undo.push_undo(to_remove, to_add)
    #
    def undo(self):
        to_remove, to_add = self._undo.pop_undo()
        rc = False
        if not to_remove is None:
            idxlist = []
            for cont in to_remove:
                idx = self.indexOf(cont)
                if idx >= 0:
                    idxlist.append(idx)
            if len(idxlist) > 0:
                self.delMultiple(idxlist)
                rc = True
        if not to_add is None:
            for cont in to_add:
                self.append(cont)
            rc = True
        return rc
    #
    @property
    def selectedCellId(self):
        if self._mgr.selected_id in self._mgr.cells:
            return self._mgr.selected_id
        return None
    @selectedCellId.setter
    def selectedCellId(self, id):
        if id in self._mgr.cells:
            self._mgr.selected_id = id
        else:
            self._mgr.selected_id = None
    #
    def deleteSelectedCell(self):
        if self._mgr.deleteSelected():
            self.reload()
            return True
        return False
    #
    def selectCell(self, idx):
        if idx >= 0:
            try:
                _id, _idx, _cont = self._ided[idx]
                if self._mgr.selected_id != _id:
                    self._mgr.selected_id = _id
                    return True
                return False
            except Exception:
                pass
        if not self._mgr.selected_id is None:
            self._mgr.selected_id = None
            #print ('No selected cell ID')
            return True
        return False
    #
    def selectedIndices(self):
        if self._mgr.selected_id is None:
            return set()
        return set([idx for idx, (_id, _idx, _cont) in enumerate(self._ided) \
                    if _id == self._mgr.selected_id])
    #


class ContourManager3d(object):
    def __init__(self, w, h, nfr, cache_fpath):
        self.w = w
        self.h = h
        self.nfr = nfr
        self.cache_fpath = cache_fpath
        #
        self.cells = {}
        self._selected_id = None
        self.undo_cache = defaultdict(UndoStack)
        #
        self._dirty = False
        self.load()
    #
    @property
    def dirty(self):
        return self._dirty
    @dirty.setter
    def dirty(self, st):
        self._dirty = st
    #
    @property
    def selected_id(self):
        return self._selected_id
    @selected_id.setter
    def selected_id(self, id):
        self._selected_id = id
    #
    def deleteSelected(self):
        if self.selected_id in self.cells:
            del self.cells[self.selected_id]
            self.selected_id = None
            return True
        return False
    #
    def loadSegmentedCells(self, csv_file):
        self.undo_cache.clear()
        self._selected_id = None
        self.cells = {}
        flat = rpeutil.import_3d_contours(self.w, self.h, self.nfr, csv_file)
        n_flat = len(flat)
        i = 0
        cnt = 0
        while i < n_flat:
            z = flat[i]
            id = flat[i+1]
            sz = flat[i+2]
            i += 3
            cont = []
            for n in range(sz):
                cont.append((flat[i], flat[i+1]))
                i += 2
            if id in self.cells:
                cell = self.cells[id]
            else:
                self.cells[id] = cell = aCell(id, self.nfr)
            cell.addContour(z, cont)
            cnt += 1
        # print ('loadSegmentedCells(): decoded', cnt, 'contours from', csv_file)
        self._dirty = True
        self.save()
    #
    @property
    def empty(self):
        return len(self.cells) == 0
    #
    def deleteAll(self):
        self.undo_cache.clear()
        self.cells = {}
        self._selected_id = None
        self._dirty = True
        self.save()
    #
    def countContours(self):
        return sum(cell.count() for cell in self.cells.values())
    #
    def getContourList(self, z):
        return zContourList(self, z)
    #
    def getAnnotations(self):
        self.save()
        res = [[] for z in range(self.nfr)]
        for id, cell in self.cells.items():
            for z in range(self.nfr):
                for cont in cell.zmap[z]:
                    all_x = [x for x, _ in cont]
                    all_y = [y for _, y in cont]
                    res[z].append({
                        "shape_attributes" : {
                            "name" : "polygon",
                            "all_points_x" : all_x,
                            "all_points_y" : all_y},
                        "region_attributes" : {
                            "cell" : id,
                            "frame" : z,
                            },}
                        )
        return res
    #
    def nextCellId(self):
        if len(self.cells) == 0:
            return 0
        return max(self.cells.keys()) + 1
    def addSlice(self, z, cont):
        xl, xr, yt, yb = contourBoundary(cont)
        if xl is None:
            return None, None, None
        zs = max(0, z - 2)
        ze = min(self.nfr - 1, z + 2)
        bnd = (xl, xr, yt, yb, zs, ze)
        #
        flat_cells = []
        for cell in self.cells.values():
            if cell.intersection(bnd) > 0:
                flat_cells.extend(cell.flat())
        flat_cont = flattenContour(cont)
        res = rpeutil.cell_for_contour(self.w, self.h, self.nfr, z, flat_cells, flat_cont)
        best_id = res[0]
        res = _unflatten_contours(res[1:])
        if len(res) > 0:
            cont = res[0]
            if len(cont) < 5:
                return None, None, None

        if len(res) > 1 and len(res[1]) > 0:
            #print ('Cell ID to split:', best_id)
            old_cell = self.cells[best_id]
            best_id = self.nextCellId()
            best_cell = aCell(best_id, self.nfr)
            self.cells[best_id] = best_cell
            old_cell.moveContoursTo(best_cell, res[1])
        elif len(res) > 2 and len(res[2]) > 0:
            #print ('Cell IDs to join:', res[2])
            for keep_id, join_id in res[2]:
                self.cells[keep_id].join(self.cells[join_id])
        #
        if best_id < 0:
            best_id = self.nextCellId()
            best_cell = aCell(best_id, self.nfr)
            self.cells[best_id] = best_cell
        else:
            best_cell = self.cells[best_id]
        #   
        _idx = best_cell.addContour(z, cont)
        return best_id, _idx, cont
    #
    def validateContour(self, cont, z, _id, _idx):
        flat_slices = []
        for id, cell in self.cells.items():
            for idx, _cont in enumerate(cell.zmap[z]):
                if id == _id and idx == _idx:
                    continue
                flat_slices.extend(flattenContour(_cont))
        flat_cont = flattenContour(cont)
        res = rpeutil.validate_contour(self.w, self.h, flat_slices, flat_cont)
        return _unflatten_contours(res)[0]
    #
    def _optimize_cells(self):
        self._selected_id = None
        upd = {}
        next_id = 0
        for cell in self.cells.values():
            if not cell.validate(): continue
            cell.id = next_id
            upd[next_id] = cell
            next_id += 1
        # print ('Optimize:', len(self.cells), '->', len(upd))
        self.cells = upd
    #
    def get_flat_cells(self):
        #self._optimize_cells()
        flat_cells = []
        for id in sorted(self.cells.keys()):
            cell = self.cells[id]
            flat_cells.extend(cell.flat())
        return flat_cells
    def getEmptyMask(self):
        return np.empty(shape=(self.nfr, self.h, self.w), dtype=np.uint8)
    #
    def _filter_cells_at_border(self):
        flat_cells = self.get_flat_cells()
        at_bord_ids = rpeutil.cells_at_border(self.w, self.h, flat_cells)
        if len(at_bord_ids) == 0:
            return
        upd = {}
        at_bord_ids = set(at_bord_ids)
        next_id = 0
        for cell in self.cells.values():
            if not cell.validate(): continue
            if cell.id in at_bord_ids: continue
            cell.id = next_id
            upd[next_id] = cell
            next_id += 1
        # print ('At border:', len(at_bord_ids), ' of ', len(self.cells))
        self.cells = upd
    #
    def saveSegmented(self, csv_file, tif_file, validate=True, separate=True):
        ch_mask = self.getEmptyMask()
        flat_cells = self.get_flat_cells()
        rpeutil.export_3d_contours(ch_mask, flat_cells, csv_file, validate, separate)
        print ('Write:', tif_file)
        imwrite(tif_file, ch_mask, compress=6, photometric='minisblack')
    #
    def savePickled(self, pkl_file):
        with open(pkl_file, 'wb') as fo:
            pickle.dump(self.cells, fo, pickle.HIGHEST_PROTOCOL)
    def loadPickled(self, pkl_file):
        self.undo_cache.clear()
        self._selected_id = None
        with open(pkl_file, 'rb') as fi:
            cells = pickle.load(fi)
            for idx, cell in cells.items():
                if cell.nfr != self.nfr:
                    raise RuntimeError('Invalid or incompatible file format')
                break
        self.cells = cells
        self._filter_cells_at_border()
        self._dirty = True
        self.save()
    #
    def load(self):
        try:
            with open(self.cache_fpath, 'rb') as fi:
                self.cells = pickle.load(fi)
            self._filter_cells_at_border()
        except Exception:
            pass
    #
    def save(self):
        if self._dirty:
            self._optimize_cells()
            try:
                self.savePickled(self.cache_fpath)
            except Exception:
                pass
            self._dirty = False
    #

if __name__ == '__main__':

    from src.RPE_Mask_RCNN_Jan22.rpefs import RpeTiffStack
    stack_path = r'C:\RPE_Data\Z01_suite\P1-W3-ZO1_D02_F006.rpe.json'
    tifstk = RpeTiffStack(stack_path)
    csv_file = tifstk.cache.getSegmName('DNA', ext='.csv')
    print (csv_file)
    pkl_file = tifstk.cache.getCacheName('DNA', ext='.pkl')
    print (pkl_file)
    mgr = ContourManager3d(tifstk.width, tifstk.height, tifstk.n_frames, pkl_file)
    print ('#cells:', len(mgr.cells))
    print ('#contours:', mgr.countContours())

    import json
    start_ts = datetime.datetime.now()
    annMap = mgr.getAnnotations()
    elapsed = datetime.datetime.now() - start_ts
    with open('annotations.json', 'w') as fo:
        json.dump(annMap, fo)
    print('getAnnotations(), elapsed:', str(elapsed))
    #mgr.loadSegmentedCells(csv_file)
    #lst = mgr.getContourList(12)
    #print ('at z=12: ', len(lst))
    #mgr.save()


