"""
uncomment the methods to be tested
Algo:
1. get parameters for segmentations
2. generate set of minarea parameters
3. generate overlays with imagej
"""
import os
from src.stackio.Channel import channel

####################################ACTB####################################

from src.segmentation.segment_GFP import segmentactb
maindirpath = "C:/Users/satheps/PycharmProjects/Results/2022/Apr1/filament/Stacks/ACTB/"
savedirpath = "C:/Users/satheps/PycharmProjects/Results/2022/Apr1/filament/Segmented/ACTB/"

assert os.path.exists(maindirpath)
assert os.path.exists(savedirpath)
minareas = [0]
# usechannels = ["actb"]#, "ctnnb1", "myh10"]

methods = ["both","3dvessellness", "2dvessellness"]
# methods = ["2dvessellness"]
# usechanneldirnames  = [channel.dirnames[c] for c in usechannels]
print(os.listdir(maindirpath))
for f in os.listdir(maindirpath):
    fpath = os.path.join(maindirpath, f)
    for method in methods:
        for minarea in minareas:
            # savepath = os.path.join(savedirpath, f"{f}_{method}_minarea{minarea}")
            print(fpath)
            # print(channel.getdefaultparams("actb")[0])
            params, _ = channel.getdefaultparams("actb")
            segmentactb(fpath, savedirpath, channel="actb", params=params, minarea=minarea,
                        method=method)
####################################ACTB####################################

####################################MYH####################################

# from src.segmentation.segment_GFP import segmentmyh
# maindirpath = "C:/Users/satheps/PycharmProjects/Results/2022/Apr1/filament/Stacks/MYH/"
# savedirpath = "C:/Users/satheps/PycharmProjects/Results/2022/Apr1/filament/Segmented/MYH/"
#
# assert os.path.exists(maindirpath)
# assert os.path.exists(savedirpath)
# minareas = [0]
#
# methods = ["both","3dvessellness", "2dvessellness"]
# # methods = ["2dvessellness"]
# # usechanneldirnames  = [channel.dirnames[c] for c in usechannels]
# print(os.listdir(maindirpath))
# for f in os.listdir(maindirpath):
#     fpath = os.path.join(maindirpath, f)
#     for method in methods:
#         for minarea in minareas:
#             # savepath = os.path.join(savedirpath, f"{f}_{method}_minarea{minarea}")
#             print(fpath)
#             # print(channel.getdefaultparams("actb")[0])
#             params, _ = channel.getdefaultparams("myh10")
#             segmentmyh(fpath, savedirpath, channel="myh10", params=params, minarea=minarea,
#                         method=method)
####################################MYH####################################

####################################CTNNB####################################

# from src.segmentation.segment_GFP import segmentctnnb
# maindirpath = "C:/Users/satheps/PycharmProjects/Results/2022/Apr1/filament/Stacks/CTNNB/"
# savedirpath = "C:/Users/satheps/PycharmProjects/Results/2022/Apr1/filament/Segmented/CTNNB/"
#
# assert os.path.exists(maindirpath)
# assert os.path.exists(savedirpath)
# minareas = [0]
# # usechannels = ["actb"]#, "ctnnb1", "myh10"]
#
# methods = ["both","2dspot", "2dvessellness"]
# # methods = ["2dvessellness"]
# # usechanneldirnames  = [channel.dirnames[c] for c in usechannels]
# print(os.listdir(maindirpath))
# for f in os.listdir(maindirpath):
#     fpath = os.path.join(maindirpath, f)
#     for method in methods:
#         for minarea in minareas:
#             # savepath = os.path.join(savedirpath, f"{f}_{method}_minarea{minarea}")
#             print(fpath)
#             # print(channel.getdefaultparams("actb")[0])
#             params, _ = channel.getdefaultparams("ctnnb1")
#             segmentctnnb(fpath, savedirpath, channel="ctnnb1", params=params, minarea=minarea,
#                         method=method)
####################################CTNNB####################################
