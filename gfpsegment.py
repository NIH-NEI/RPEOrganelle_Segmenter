import os
import argparse
from src.stackio.Channel import channel


def segment_all_gfpfiles(path_stackfolders, channelname, savedir):
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


def main():
    parser = argparse.ArgumentParser(
        description='Segment GFP channels from given stacks')
    parser.add_argument('-c', '--channelname', required=True,
                        metavar='name of channel',
                        default='tjp1',
                        help='Name of gfp organelle. See readme file for supported channels')
    parser.add_argument('-p', '--path_stackfolders', required=True, action="store_true",
                        help="Path to where the folders containing the stacks are stored.",
                        default='../Results/../ZO1/')
    parser.add_argument('-s', '--savedir', required=True,
                        help="Directory where segmented files should be saved",
                        default='../Results/../final_segmentations/ZO1/')

    args = parser.parse_known_args()

    assert os.path.exists(args.path_stackfolders)
    assert os.path.exists(args.savedir)

    segment_all_gfpfiles(args.path_stackfolders, args.channelname, args.savedir)


if __name__ == "__main__":
    main()
