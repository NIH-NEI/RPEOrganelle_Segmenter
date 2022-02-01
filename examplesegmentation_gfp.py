import os

# from tensorflow.python.client import device_lib
# from src.RPE_Mask_RCNN_Jan22 import generate_meta,predict
from src.stackio.Channel import channel


# def get_available_gpus():
#     local_device_protos = device_lib.list_local_devices()
#     return [x.napip install PyWaveletsme for x in local_device_protos if x.device_type == 'GPU']

#
# print("GPUS", get_available_gpus())
gpucount = '1'  # needs string input for argparse
imgspergpu = '1'  # ?

imgdir = '../Results/2021/July16/TOMstacks/'
savedir = '../Results/2022/Jan28/TOM/segmented/'
# mdir = '../Results/2022/Jan21/TOM_stack_18img/'
# tempdir = '../Results/2021/july16/TOMstacks' # generate meta doesnt work with stacks
# tempdir = '../Results/2022/Jan21/TOM/'  # generate meta doesnt work with stacks
print(os.listdir(imgdir))
for subdir in os.listdir(imgdir):
    subdirpath = os.path.join(imgdir,subdir).replace("\\","/")
    for filename in os.listdir(subdirpath):
        rpefilename = os.path.join(subdirpath, filename).replace("\\","/")
        print(rpefilename)
        gfpchannel = channel("tom20")
        params = {'f2params': [[2.75 / 3, 0.1]]}
        gfpchannel.segmentchannel(filename=rpefilename, savepath=savedir, params=params)




    # usefunction[dirname](mainfilepath, savepath, params)

# metadata_arguments = [None, tempdir]
# generate_meta.generatemeta(metadata_arguments)
#
# modelweights = 'C:/Users/satheps/PycharmProjects/Results/2022/Jan21/model_weights'
# predict_arguments = [None, '-d', tempdir, '-w', modelweights, 'rpefile', tempdir, '-g', gpucount, '-i',
#                  imgspergpu]
# predict.predict(predict_arguments)

exit()



alldirs = 'E:/Pushkar_backup/Results/data/prestack/'

todo = ['TOM']
"E02"

for dir in os.listdir(alldirs): # for multichannel
    print(dir)
    if dir in todo:
        datadir = os.path.join(alldirs, dir)
        rpefiledir = datadir

        metadata_arguments = [None, datadir]
        # generate_meta.generatemeta(metadata_arguments)
        modelweights = 'C:/Users/satheps/PycharmProjects/Results/2022/Jan21/model_weights'
        predict_arguments = [None, '-d', datadir, '-w', modelweights, 'rpefile', datadir, '-g', gpucount, '-i',
                             imgspergpu]
        predict.predict(predict_arguments)
