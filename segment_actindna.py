import os

from tensorflow.python.client import device_lib

from src.RPE_Mask_RCNN_Feb22 import generate_meta,predict
# from src.RPE_Mask_RCNN_Jan22 import predict


def get_available_gpus():
    local_device_protos = device_lib.list_local_devices()
    return [x.name for x in local_device_protos if x.device_type == 'GPU']


#
print("GPUS", get_available_gpus())
gpucount = '1'  # needs string input for argparse
imgspergpu = '1'  # ?

alldirs = 'E:/Pushkar_backup/Results/data/prestack/'
# path_stackfolders = '../Results/2021/july9/lampstacks/'

todo = ['SEC61', 'LAMP1']
# "E02"

# for dir in os.listdir(alldirs):  # for multichannel
#     print(dir)
#     if dir in todo:
#         datadir = os.path.join(alldirs, dir)
#
#         metadata_arguments = [None, datadir]
#         generate_meta.generatemeta(metadata_arguments)
#         modelweights = 'C:/Users/satheps/PycharmProjects/Results/2022/Jan21/model_weights'
#         predict_arguments = [None, '-d', datadir, '-w', modelweights, 'rpefile', datadir, '-g', gpucount, '-i',
#                              imgspergpu]
#         predict.predict(predict_arguments)


# imgdir = 'C:/Users/satheps/PycharmProjects/Results/2021/july9/lampstacks/'
imgdir = 'C:/Users/satheps/PycharmProjects/Results/2021/mar31/sec61stacks/'
# savedir = '../Results/2022/Mar4/channels/lamp1/segmented/'
# savedir = '../../Results/2022/Mar4/channels/sec61b/segmented/'
datadir = imgdir
# metadata_arguments = [None, datadir]
# generate_meta.generatemeta(metadata_arguments)
for r, readdir in enumerate(os.listdir(datadir)):
    # if r==0:
    #     continue
    modelweights = 'C:/Users/satheps/PycharmProjects/Results/2022/Jan21/model_weights'
    subdirname = os.path.join(datadir, readdir)
    # for ometiffile in os.listdir(subdirname):
    #     if ometiffile.__contains__(".ome.tif"):
    #         filepath = os.path.join(subdirname, ometiffile).replace("\\", "/")
    #         print(filepath)
    predict_arguments = [None, '-d', datadir, '-w', modelweights, 'rpefile', readdir, '-n']
    predict.predict(predict_arguments)
