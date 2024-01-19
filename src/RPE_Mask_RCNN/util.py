__all__ = ('datadir', 'timing', 'UndoEntry', 'UndoStack', 'slice_area')

import sys, os, datetime
import csv
import enum
from collections import namedtuple

# Decorator for functions accessing data via relative paths
if hasattr(sys, '_MEIPASS'):
    # When distributed as single exe, temporarily switch to the _MEIPASS prefix,
    #    call decorated function, then switch back. This allows functions accessing data
    #    via paths relative to '.' work the same way
    def datadir(func):
        def wrapper(*args, **kwarg):
            cwd = os.getcwd()
            try:
                os.chdir(sys._MEIPASS)
                rc = func(*args, **kwarg)
            finally:
                os.chdir(cwd)
            return rc
        return wrapper
    #
    def timing(func):
        return func
else:
    # Decorator for functions accessing data via relative paths
    # If distributed as single dir (or while in dev sandbox), do nothing.
    def datadir(func):
        return func
    #
    # Profiling decorator for non-single-exe env only
    def timing(func):
        def wrapper(*args, **kwarg):
            start_ts = datetime.datetime.now()
            rc = func(*args, **kwarg)
            print(func.__name__+'() done in:', str(datetime.datetime.now()-start_ts))
            return rc
        return wrapper

UndoEntry = namedtuple('UndoEntry', ['to_remove', 'to_add'])

class UndoStack(object):
    def __init__(self, maxundo=100):
        self.maxundo = maxundo
        #
        self.buf = []
    #
    def clear(self):
        self.buf[:] = []
    #
    def is_empty(self):
        return len(self.buf) == 0
    #
    def push_undo(self, to_remove, to_add):
        if to_remove:
            to_remove = tuple(to_remove)
        else:
            to_remove = None
        if to_add:
            to_add = tuple(to_add)
        else:
            to_add = None
        if to_remove is None and to_add is None:
            return False
        self.buf.insert(0, UndoEntry(to_remove, to_add))
        if len(self.buf) > self.maxundo:
            self.buf[self.maxundo:] = []
        return True
    #
    def pop_undo(self):
        if len(self.buf) == 0:
            return None, None
        ent = self.buf.pop(0)
        return ent.to_remove, ent.to_add
    #

def slice_area(asize, tsize, minovl=200):
    aw = int(asize[1])
    ah = int(asize[0])
    w = int(tsize[1])
    h = int(tsize[0])
    #
    res = []
    if aw < w or ah < h:
        return res
    #
    nxt = aw // w + 1
    xovl = (nxt*w - aw) // (nxt - 1)
    while xovl < minovl:
        nxt += 1
        xovl = (nxt*w - aw) // (nxt - 1)
    #
    nyt = ah // h + 1
    yovl = (nyt*h - ah) // (nyt - 1)
    while yovl < minovl:
        nyt += 1
        yovl = (nyt*h - ah) // (nyt - 1)
    #
    hw = w//2
    hh = h//2
    xmax = aw - 1
    ymax = ah - 1
    #
    xstep = aw // nxt
    ystep = ah // nyt
    #
    for yc in range(nyt):
        y0 = yc * ystep + ystep//2 - hh
        if y0 < 0: y0 = 0
        y1 = y0 + h - 1
        if y1 > ymax:
            y1 = ymax
            y0 = y1 - h + 1
        for xc in range(nxt):
            x0 = xc * xstep + xstep//2 - hw
            if x0 < 0: x0 = 0
            x1 = x0 + w - 1
            if x1 > xmax:
                x1 = xmax
                x0 = x1 - w + 1
            res.append((x0, x1, y0, y1))
    return res
