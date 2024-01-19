import os

from tensorflow.python.client import device_lib

from src.RPE_Mask_RCNN import generate_meta, predict


def get_available_gpus():
    local_device_protos = device_lib.list_local_devices()
    return [x.name for x in local_device_protos if x.device_type == 'GPU']


meta = False
print("GPUS", get_available_gpus())
gpucount = '1'  # needs string input for argparse
imgspergpu = '1'  # ?

imgdir = '../Results/.././laminstacks/'

datadir = imgdir
for r, readdir in enumerate(os.listdir(datadir)):
    modelweights = '../Results/../../model_weights/'
    subdirname = os.path.join(datadir, readdir)
    if meta:
        metadata_arguments = [None, datadir]
        generate_meta.generatemeta(metadata_arguments)

    # for ometiffile in os.listdir(subdirname):
    #     if ometiffile.__contains__(".ome.tif"):
    #         filepath = os.path.join(subdirname, ometiffile).replace("\\", "/")
    #         print(filepath)
    # predict_arguments = [None, '-d', datadir, '-w', modelweights, 'rpefile', subdirname, '-n'] # for no assembly
    predict_arguments = [None, '-d', datadir, '-w', modelweights, 'rpefile', subdirname]
    print(predict_arguments)
    predict.predict(predict_arguments)
