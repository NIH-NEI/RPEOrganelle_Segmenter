import os

from src.stackio.Channel import channel

path_stackfolders = '../Results/2021/july9/lampstacks/'
savedir = '../Results/2022/final_segmentations/LAMP1/'
# path_stackfolders = '../Results/2021/mar31/sec61stacks/'
channelname = "lamp1"
# path_stackfolders = '../Results/2021/july16/TOMstacks/'
# savedir = '../Results/2022/final_segmentations/TOM20/'
# channelname = "fbl"
# path_stackfolders = '../Results/2021/Sept3/FBLstacks/'
# savedir = '../Results/2022/final_segmentations/FBL/'
# channelname = "rab5"
# path_stackfolders = '../Results/2021/Sept10/rab5/'
# savedir = '../Results/2022/final_segmentations/RAB5/'
assert os.path.exists(path_stackfolders)
assert os.path.exists(savedir)

print(os.listdir(path_stackfolders))
print(os.listdir(savedir))
flist = [f for f in os.listdir(savedir) if os.path.isfile(os.path.join(savedir, f))]
fpathlist = [os.path.abspath(os.path.join(savedir, fp)).replace("\\", "/") for fp in flist]
for subdir in os.listdir(path_stackfolders):
    flag = 1
    subdirpath = os.path.join(path_stackfolders, subdir).replace("\\", "/")
    gfpfiles = [f for f in os.listdir(subdirpath) if os.path.isfile(os.path.join(subdirpath, f))]

    for filename in gfpfiles:
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
