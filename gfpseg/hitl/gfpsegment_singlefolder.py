import os
from gfpseg.stackio.Channel import channel
from gfpseg.segmentation.segment_GFP import segmentslc


maindirpath = "../Results/2022/../bkgcorr_tests/bkgcorr_all/"
savedirpath = "../Results/2022/../bkgcorr_tests/bkgcorr_seg/"

assert os.path.exists(maindirpath)
assert os.path.exists(savedirpath)
minareas = [0]
# usechannels = ["dsp"]#, "slc25a17", "gja110"]

methods = ["both","3dspot", "2dspot"]
# methods = ["3dspot"]
# methods = ["2dvessellness"]
# usechanneldirnames  = [channel.dirnames[c] for c in usechannels]
print(os.listdir(maindirpath))
for f in os.listdir(maindirpath):
    fpath = os.path.join(maindirpath, f)
    for method in methods:
        for minarea in minareas:
            # savepath = os.path.join(savedirpath, f"{f}_{method}_minarea{minarea}")
            print(fpath)
            # print(channel.getdefaultparams("dsp")[0])
            params, _ = channel.getdefaultparams("slc25a17")
            segmentslc(fpath, savedirpath, channel="slc25a17", params=params, minarea=minarea, method=method)