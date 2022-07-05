import os

from src.stackio.Channel import channel

channelname = "tjp1"

path_stackfolders = '../Results/2022/May27/ZO1/'
savedir = '../Results/2022/final_segmentations/ZO1/'
assert os.path.exists(path_stackfolders)
assert os.path.exists(savedir)

print(os.listdir(path_stackfolders))
# same as listdir ignoring directories
flist = [f for f in os.listdir(savedir) if os.path.isfile(os.path.join(savedir, f))]
# abspath for list
fpathlist = [os.path.abspath(os.path.join(savedir, fp)).replace("\\", "/") for fp in flist]
print(len(fpathlist))
for subdir in os.listdir(path_stackfolders):

    subdirpath = os.path.join(path_stackfolders, subdir).replace("\\", "/")
    gfpfiles = [f for f in os.listdir(subdirpath) if os.path.isfile(os.path.join(subdirpath, f))]

    for filename in gfpfiles:
        flag = 1
        rpefilename = os.path.join(subdirpath, filename).replace("\\", "/")
        print(filename)
        # print(rpefilename.split("."))
        # print(rpefilename.split(".")[-3])
        basepath = os.path.abspath(savedir + filename.split(".")[-3]).replace("\\", "/")
        for i in range(len(fpathlist)):
            if basepath in fpathlist[i]:
                print(f"basepath already exists. Ignore segmentation [{basepath} in {fpathlist[i]}]", flush=True)
                flag = 0
        if flag:
            print("basepath: ", basepath, flush=True)
            print("rpefilename: ", rpefilename, flush=True)
            channel.segmentchannel(filename=rpefilename, savepath=savedir, channelname=channelname)
