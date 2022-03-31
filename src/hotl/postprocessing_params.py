"""
Algo:
1. get parameters for segmentations
2. generate set of minarea parameters 
3. ensure they work
4. generate overlays with imagej - attempt with 
"""
import os

from src.stackio.Channel import channel

maindirpath = "C:/Users/satheps/PycharmProjects/Results/2021/Oct8/Stacks/"
savedirpath = "C:/Users/satheps/PycharmProjects/Results/2022/Apr1/postproc_p3_minarea_iter2/"

assert os.path.exists(maindirpath)
assert os.path.exists(savedirpath)
# gfpchannel = channel("sec61b")
minareas = [0, 1, 2, 3, 4, 5]
# usechannels = ["dsp", "gja1", "slc25a17"]
usechannels = ["slc25a17"]
params = {
    # "dsp":[],
    # "gja1":[],
    "slc25a17":[[1.03125,0.07]]

}
usechanneldirnames  = [channel.dirnames[c] for c in usechannels]
print(os.listdir(maindirpath))
for subdir in os.listdir(maindirpath):
    # print(subdir, channel.dirnames.values())
    assert subdir in channel.dirnames.values()
    if subdir in usechanneldirnames:
        usechannel = usechannels[usechanneldirnames.index(subdir)]
        subdirpath = os.path.join(maindirpath, subdir).replace("\\", "/")
        savesubdirpath = os.path.join(savedirpath, subdir).replace("\\", "/")
        if not os.path.exists(savesubdirpath):
            os.mkdir(savesubdirpath)
        for filename in os.listdir(subdirpath):
            rpefilename = os.path.join(subdirpath, filename).replace("\\", "/")
            print(rpefilename)
            for minarea in minareas:
                try:
                    channel.segmentchannel(filename=rpefilename, savepath=savesubdirpath+"/",params=params[usechannel], channelname=usechannel, minarea=minarea)
                except KeyError as k:
                    print("keyerror, using current default parameters")
                    channel.segmentchannel(filename=rpefilename, savepath=savesubdirpath + "/", channelname=usechannel, minarea=minarea)