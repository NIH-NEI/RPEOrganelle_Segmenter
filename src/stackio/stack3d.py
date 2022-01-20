import datetime
import os
import warnings

import numpy as np
from aicsimageio.writers.ome_tiff_writer import OmeTiffWriter
from tifffile import tifffile

from src.stackio import types, experimentalparams, Channel


class Stack():
    """
    Class for Stack object and its functions
    """
    alphabets = ["B", "C", "D", "E", "F", "G"]
    # alphabets = ["D"]

    def __init__(self, alphabet=None, channelname=None):
        if alphabet is None:
            warnings.warn(
                "if alphabet is None all acceptable alphabets will be tested. This may slow speed down significantly")
        self.reps = experimentalparams.generate_repinfo(_alphabets=alphabet)
        self.filepath = None
        self.savepath = None  # "../Results/2021/Sept24/actb/"
        self.dirs = None
        self.channels = ["C01", "C02", "C03", "C04"]
        self.Fs = ["F001", "F002", "F003", "F004", "F005", "F006"]
        self.channelname = channelname
        self.channel = Channel.channel(self.channelname)

    def findmesfile(self, flist):
        """
        locate and return the first encountered .mes file in a list of files. This file contains information about
        stacks that was missing from their filenames.

        :param flist: list of files
        :return: filename of .mes file; or None if no file found
        """
        for f in flist:
            if f.__contains__(".mes"):
                return f
        print("No .mes file found")
        return None

    def getdirectorylist(self, filepath: types.PathLike = ""):
        return os.listdir(filepath)

    def imagestostack(self, tiff_files, w, rep, foo):
        """

        :param tiff_files: list of
        :param w:
        :param rep:
        :param foo:
        :return:
        """
        flag = 0
        ims = []
        for coo in self.channels:
            # selectedfilelist = self.returnfilescontaining(listofallfiles=tiff_files, w, rep, foo, coo)
            cims = []
            for f in tiff_files:
                if self.filehasparameters(f, w, rep, foo, coo):
                    if w in f and rep in f and foo in f and coo in f:
                        flag = 1
                        # print(w, rep, foo, coo, f)
                        try:
                            cim = tifffile.imread(os.path.join(self.dpath, f))
                            cims.append(cim)
                        except Exception as e:
                            print("could not read", f, flush=True)
            cimages = np.asarray(cims)
            ims.append(cimages)
            del cims
        return flag, ims

    def filehasparameters(self, filename, w_in, r_in, f_in, c_in):
        """
        obtains the parameters from individual image files
        TODO: finalize and test

        :param filename: 
        :param w_in: 
        :param r_in: 
        :param f_in: 
        :param c_in: 
        :return: 
        """
        w, r, w_, r_, filename = experimentalparams.getwrprestack(filename)
        if (w == w_in) and (r == r_in):  # and (f == f_in) and (c == c_in):
            pass

    def returnstacksfromlist(self, filepath, savepath):
        """

        :param filepath: file path
        :param savepath: save path
        :return:
        """
        self.filepath = filepath
        self.savepath = savepath
        self.dirs = os.listdir(filepath)

        for dirname in self.dirs:
            dpath = os.path.join(self.filepath, dirname)
            if os.path.isdir(dpath):
                files = [f for f in os.listdir(dpath) if os.path.isfile(os.path.join(dpath, f))]
                mesfile = self.findmesfile(files)
                mesweekid = mesfile[:-4].split("-")[3:5]
                # print(mesweekid, mesfile)

                mesweekid = "-".join(["P1", mesweekid[1], mesweekid[0]])
                fsavepath = os.path.join(self.savepath, mesweekid)
                print(fsavepath, mesweekid, len(files))
                tiff_files = [f for f in files if f.__contains__(".tif")]
                tiff_files = sorted(tiff_files)
                if not os.path.exists(fsavepath):
                    os.mkdir(fsavepath)
                # for w in experimentalparams.WS:
                w = mesweekid[1]
                # done = self.getcompletedfilelist()
                # if w in done:
                #     continue
                for rep in self.reps:
                    for foo in self.Fs:
                        omefilename = "_".join([mesweekid, rep, foo + ".ome.tif"])
                        if os.path.exists(os.path.join(fsavepath, omefilename)):
                            print(os.path.join(fsavepath, omefilename) + " already done")
                        else:
                            try:
                                start_ts = datetime.datetime.now()
                                for coo in self.channels:
                                    flag, ims = self.imagestostack(tiff_files, w, rep, foo, coo)
                                    if flag:
                                        images = np.asarray(ims)
                                        expanded = np.expand_dims(images, 0)
                                        # print(expanded.shape) # TCZYX

                                        # expanded = expanded.transpose(0, 2, 1, 3, 4)
                                        # print(images.shape, expanded.shape, cimages.shape, expanded.transpose(0, 2, 1, 3, 4).shape)
                                        # print(expanded.shape) # TCZYX
                                        # raise Exception

                                        OmeTiffWriter.save(expanded, os.path.join(fsavepath, omefilename))
                                        # writer = omeTifWriter.OmeTifWriter(os.path.join(fsavepath, omefilename),
                                        #                                    overwrite_file=False)
                                        # writer.save(expanded)
                                        print("saved file: ", omefilename, expanded.shape)
                                        end_ts = datetime.datetime.now()
                                        print("Calculated in ", end_ts - start_ts, "@ :", end_ts)
                                        # raise Exception
                            except Exception as e:
                                print(e)
