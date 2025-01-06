import datetime
import os
import warnings

import numpy as np
from aicsimageio.writers.ome_tiff_writer import OmeTiffWriter
from tifffile import tifffile

from gfpseg.stackio import types, experimentalparams, Channel


class Stack():
    """
    Class for Stack object and its functions
    """
    alphabets = ["B", "C", "D", "E", "F", "G"]

    # alphabets = ["D"]

    def __init__(self, alphabet: str = None, channelname: str = None):
        """
        Args:
            alphabet:alphabet (e.g. A) used in naming the file e.g. P1-W2-*GFP*_*A04*_F003
            channelname: name of channel
        Returns:

        """
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
        # self.channel = Channel.channel(self.channelname)

    def findmesfile(self, flist: types.PathLike):
        """
        locate and return the first encountered .mes file in a list of files. This file contains information about
        stacks that was missing from their filenames.
        Args:
            flist: list of files
        Returns:
            filename of .mes file; or None if no file found
        """
        for f in flist:
            if f.__contains__(".mes"):
                return f
        print("No .mes file found")
        return None

    def getdirectorylist(self, filepath: types.PathLike = ""):
        """
        returns subfiles and directories.

        Args:
            filepath: filepath
        Returns:
            list of subfiles/subdirectories
        """
        return os.listdir(filepath)

    def imagestostack(self, tiff_files, dpath, rep, foo):
        """
        Convert a list of tiff files to a stack
        Args:
            tiff_files: list of tiff files
            rep: replicate
            foo: FOV
        Returns:

        """
        flag = 0
        ims = []
        for coo in self.channels:
            # selectedfilelist = self.returnfilescontaining(listofallfiles=tiff_files, w, rep, foo, coo)
            cims = []
            for f in tiff_files:
                if self.filehasparameters(f, rep, foo, coo):
                    # print(rep, foo, coo, f)

                    # if rep in f and foo in f and coo in f:
                    flag = 1
                    try:
                        cim = tifffile.imread(os.path.join(dpath, f))
                        cims.append(cim)
                    except Exception as e:
                        print("could not read", f, flush=True)
            cimages = np.asarray(cims)
            ims.append(cimages)
            del cims
        return flag, ims

    def filehasparameters(self, filename, r_in, f_in, c_in):
        """
        check for input parameters from individual image files
        Args:
            filename: filename
            w_in: Week no.
            r_in: replicate no.
            f_in: FOV no.
            c_in: Well plate no.
        Returns: 
            True if all parameters are part of filename. False otherwise
        """
        # only works for properly named files
        if not filename.endswith("tif"):
            return False
        if filename.__contains__("CH") or filename.__contains__("CMOS"):
            return False
        try:
            r, f, c = experimentalparams.get_rfcw(filename)
        except:
            print("file naming issue")

        if (r == r_in) and (f == f_in) and (c == c_in):
            return True
        else:
            return False

    def generatestacksfromdirs(self, filepath: types.PathLike, savepath: types.PathLike, writestack: bool=True):
        """
        generates stacks from list and saves then at savepath.
        Dir structure:
            savepath
                |
                |--W1 folder----|
                |--W2 folder    |--list of files for respective week.
                :
        Args:
            filepath: file path
            savepath: save path
        Returns:
            
        """
        self.filepath = filepath
        self.savepath = savepath
        self.dirs = os.listdir(filepath)

        for dirname in self.dirs:  # dirname + mesweekid for week
            dpath = os.path.join(self.filepath, dirname)
            if os.path.isdir(dpath):
                files = [f for f in os.listdir(dpath) if os.path.isfile(os.path.join(dpath, f))]
                mesfile = self.findmesfile(files)
                if not mesfile: continue
                mesweekid = mesfile[:-4].split("-")[3:5]
                print(mesweekid, mesfile)
                mesweekid = "-".join(["P1", mesweekid[1], mesweekid[0]])
                fsavepath = os.path.join(self.savepath, mesweekid)
                print(fsavepath, mesweekid, len(files))
                tiff_files = [f for f in files if f.__contains__(".tif")]
                tiff_files = sorted(tiff_files)
                # w = mesweekid[1]
                if not os.path.exists(fsavepath):
                    os.mkdir(fsavepath)
                # done = self.getcompletedfilelist()
                # if w in done:
                #     continue
                # for w in experimentalparams.WS:
                for rep in self.reps:
                    for foo in self.Fs:
                        omefilename = "_".join([mesweekid, rep, foo + ".ome.tif"])
                        if os.path.exists(os.path.join(fsavepath, omefilename)):
                            print(os.path.join(fsavepath, omefilename) + " already done")
                        else:
                            try:
                                start_ts = datetime.datetime.now()
                                flag, ims = self.imagestostack(tiff_files, dpath, rep, foo)
                                if flag:
                                    images = np.asarray(ims)
                                    expanded = np.expand_dims(images, 0)
                                    # note: if stack dimensions are in a different order uncomment and modify the following code to fix it
                                    # print(expanded.shape) # TCZYX
                                    # expanded = expanded.transpose(0, 2, 1, 3, 4)
                                    # print(images.shape, expanded.shape)
                                    # print(expanded.shape) # TCZYX

                                    if writestack:
                                        OmeTiffWriter.save(expanded, os.path.join(fsavepath, omefilename), compression='zlib')
                                        print("saved file: ", omefilename, expanded.shape)
                                        end_ts = datetime.datetime.now()
                                        print("Calculated in ", end_ts - start_ts, "@ :", end_ts)

                            except Exception as e:
                                print(e)


if __name__ == "__main__":
    channelname = "ZO1"
    s = Stack(alphabet=Channel.channel.getrepalphabet(channelname=channelname))
    filepath = "../data/prestack/ZO1/"
    savepath = "../Results/2022/May27/ZO1/"
    s.generatestacksfromlist(filepath=filepath, savepath=savepath)
