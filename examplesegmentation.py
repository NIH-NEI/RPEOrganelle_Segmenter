import os

from tensorflow.python.client import device_lib

from src.RPE_Mask_RCNN_Jan22 import generate_meta
from src.RPE_Mask_RCNN_Jan22 import predict


def get_available_gpus():
    local_device_protos = device_lib.list_local_devices()
    return [x.name for x in local_device_protos if x.device_type == 'GPU']


print("GPUS", get_available_gpus())
gpucount = '1'  # needs string input for argparse
imgspergpu = '1'  # ?

# alldirs = 'E:/Pushkar_backup/Results/data/prestack/'
alldirs = '../Results/2022/Jan21/TOM_stack_18img'

completed = ['ACTB', 'CETN2', 'CTNNB1']
for dir in os.listdir(alldirs):
    print(dir)
    if dir not in completed:
        datadir = os.path.join(alldirs, dir)
        rpefiledir = datadir

        metadata_arguments = [None, datadir]
        generate_meta.generatemeta(metadata_arguments)
        modelweights = 'C:/Users/satheps/PycharmProjects/Results/2022/Jan21/model_weights'
        predict_arguments = [None, '-d', datadir, '-w', modelweights, 'rpefile', rpefiledir, '-g', gpucount, '-i',
                             imgspergpu]
        predict.predict(predict_arguments)
