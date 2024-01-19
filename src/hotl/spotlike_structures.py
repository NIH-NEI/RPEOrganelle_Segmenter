"""
uncomment the methods to be tested
Algo:
1. get parameters for segmentations
2. generate set of minarea parameters
3. generate overlays with imagej
"""
import os
from src.stackio.Channel import channel

####################################dsp####################################
#
# from src.segmentation.segment_GFP import segmentdsp
# maindirpath = "../Results/../spotlike/Stacks/DSP/"
# savedirpath = "../Results/../spotlike/Segmented/DSP/"
#
# assert os.path.exists(maindirpath)
# assert os.path.exists(savedirpath)
# minareas = [0]
# # usechannels = ["dsp"]#, "slc25a17", "gja110"]
#
# methods = ["both","3dspot", "2dspot"]
# # methods = ["2dvessellness"]
# # usechanneldirnames  = [channel.dirnames[c] for c in usechannels]
# print(os.listdir(maindirpath))
# for f in os.listdir(maindirpath):
#     fpath = os.path.join(maindirpath, f)
#     for method in methods:
#         for minarea in minareas:
#             # savepath = os.path.join(savedirpath, f"{f}_{method}_minarea{minarea}")
#             print(fpath)
#             # print(channel.getdefaultparams("dsp")[0])
#             params, _ = channel.getdefaultparams("dsp")
#             segmentdsp(fpath, savedirpath, channel="dsp", params=params, minarea=minarea,
#                         method=method)
####################################dsp####################################

####################################gja1####################################

# from src.segmentation.segment_GFP import segmentgja
# maindirpath = "../Results/../spotlike/Stacks/GJA1/"
# savedirpath = "../Results/../spotlike/Segmented/GJA1/"
#
# assert os.path.exists(maindirpath)
# assert os.path.exists(savedirpath)
# minareas = [0]
#
# methods = ["both","3dspot", "2dspot"]
# # methods = ["2dvessellness"]
# # usechanneldirnames  = [channel.dirnames[c] for c in usechannels]
# print(os.listdir(maindirpath))
# for f in os.listdir(maindirpath):
#     fpath = os.path.join(maindirpath, f)
#     for method in methods:
#         for minarea in minareas:
#             # savepath = os.path.join(savedirpath, f"{f}_{method}_minarea{minarea}")
#             print(fpath)
#             # print(channel.getdefaultparams("dsp")[0])
#             params, _ = channel.getdefaultparams("gja1")
#             segmentgja(fpath, savedirpath, channel="gja1", params=params, minarea=minarea,
#                         method=method)
####################################gja1####################################

####################################SLC####################################

# from src.segmentation.segment_GFP import segmentslc
# maindirpath = "../Results/../spotlike/Stacks/SLC/"
# savedirpath = "../Results/../spotlike/Segmented/SLC/"
#
# assert os.path.exists(maindirpath)
# assert os.path.exists(savedirpath)
# minareas = [0]
# # usechannels = ["dsp"]#, "slc25a17", "gja110"]
#
# methods = ["both","3dspot", "2dspot"]
# # methods = ["2dvessellness"]
# # usechanneldirnames  = [channel.dirnames[c] for c in usechannels]
# print(os.listdir(maindirpath))
# for f in os.listdir(maindirpath):
#     fpath = os.path.join(maindirpath, f)
#     for method in methods:
#         for minarea in minareas:
#             # savepath = os.path.join(savedirpath, f"{f}_{method}_minarea{minarea}")
#             print(fpath)
#             # print(channel.getdefaultparams("dsp")[0])
#             params, _ = channel.getdefaultparams("slc25a17")
#             segmentslc(fpath, savedirpath, channel="slc25a17", params=params, minarea=minarea,
#                         method=method)
####################################SLC####################################
from src.segmentation.segment_GFP import segmentstgal
maindirpath = "../Results/../Stacks_hotl/ST6GAL1/"
savedirpath = "../Results/../Hotl_phases/test/ST6GAL1/"

assert os.path.exists(maindirpath)
assert os.path.exists(savedirpath)
minareas = [0]
# usechannels = ["dsp"]#, "slc25a17", "gja110"]

methods = ["nomo","nothin", ""]
# methods = ["2dvessellness"]
# usechanneldirnames  = [channel.dirnames[c] for c in usechannels]
print(os.listdir(maindirpath))
for f in os.listdir(maindirpath):
    fpath = os.path.join(maindirpath, f)
    for method in methods:
        for minarea in minareas:
            # savepath = os.path.join(savedirpath, f"{f}_{method}_minarea{minarea}")
            print(fpath)
            # params, _ = channel.getdefaultparams("ST6GAL1")
            params = {"s3params": [[3.75/3, 0.1]], "topothin": 1}
            segmentstgal(fpath, savedirpath, channel="ST6GAL1", params=params, minarea=minarea, method=method)