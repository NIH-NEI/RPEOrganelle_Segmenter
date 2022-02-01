import os
import sys
import json
import random
import datetime
import shutil
import numpy as np
import skimage.draw
from imgaug import augmenters as iaa

# Suppress warning messages from TensorFlow, since they tend to overwhelm the output
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

# Root directory of the project
ROOT_DIR = os.path.abspath(os.path.dirname(__file__))
DEFAULT_DATA_DIR = ROOT_DIR

# Import Mask RCNN
if not ROOT_DIR in sys.path:
    sys.path.append(ROOT_DIR)  # To find local version of the library
from mrcnn.config import Config
from mrcnn import model as modellib, utils

############################################################
#  Configurations
############################################################

class RpeTrainConfig(Config):
    # Give the configuration a recognizable name
    NAME = "DNA"
    
    TRAIN_VAL_SPLIT = 0.75
    LAYERS_TO_TRAIN = 'all'
    NO_OF_EPOCHS = 10
    INITIAL_EPOCH = 0

    # Limit # of multiprocess workers if > 0
    # Good for shared systems with limits on threads/processes
    MAX_WORKERS = 0

    # NUMBER OF GPUs to use. When using only a CPU, this needs to be set to 1.
    GPU_COUNT = 1
    # Approximately one 1024x1024 (upscaled) image per 6GB of GPU memory.
    # If you have <6GB GPU memory, it is best to disable GPU.
    IMAGES_PER_GPU = 1

    # Number of classes (including background)
    NUM_CLASSES = 1 + 1  # Background + DNA/Actin

    # Number of training steps per epoch
    STEPS_PER_EPOCH = 200
    VALIDATION_STEPS = 50

    # Skip detections with < 90% confidence
    DETECTION_MIN_CONFIDENCE = 0.9

    LEARNING_RATE = 0.001
    LEARNING_MOMENTUM = 0.9

    # Maximum number of ground truth instances to use in one image
    MAX_GT_INSTANCES = 1000

    # Input image resizing
    # Generally, use the "square" resizing mode for training and predicting
    # and it should work well in most cases. In this mode, images are scaled
    # up such that the small side is = IMAGE_MIN_DIM, but ensuring that the
    # scaling doesn't make the long side > IMAGE_MAX_DIM. Then the image is
    # padded with zeros to make it a square so multiple images can be put
    # in one batch.
    # Available resizing modes:
    # none:   No resizing or padding. Return the image unchanged.
    # square: Resize and pad with zeros to get a square image
    #         of size [max_dim, max_dim].
    # pad64:  Pads width and height with zeros to make them multiples of 64.
    #         If IMAGE_MIN_DIM or IMAGE_MIN_SCALE are not None, then it scales
    #         up before padding. IMAGE_MAX_DIM is ignored in this mode.
    #         The multiple of 64 is needed to ensure smooth scaling of feature
    #         maps up and down the 6 levels of the FPN pyramid (2**6=64).
    # crop:   Picks random crops from the image. First, scales the image based
    #         on IMAGE_MIN_DIM and IMAGE_MIN_SCALE, then picks a random crop of
    #         size IMAGE_MIN_DIM x IMAGE_MIN_DIM. Can be used in training only.
    #         IMAGE_MAX_DIM is not used in this mode.
    IMAGE_RESIZE_MODE = "square"
    IMAGE_MIN_DIM = 800
    IMAGE_MAX_DIM = 1024

#
def find_weights_path(wpath, chname):
    wpath = os.path.abspath(wpath)
    if os.path.isfile(wpath):
        return wpath
    if not os.path.isdir(wpath):
        return None
    prefix = 'mask_rcnn_' + chname.lower() + '_'
    epoch = -1
    res = None
    for fn in os.listdir(wpath):
        if not fn.startswith(prefix):
            continue
        fpath = os.path.join(wpath, fn)
        if not os.path.isfile(fpath):
            continue
        try:
            bn, _ = os.path.splitext(fn)
            bn = bn[len(prefix):].replace('-', '_').split('_')[0]
            _epoch = int(bn)
            if _epoch > epoch:
                epoch = _epoch
                res = fpath
        except Exception:
            pass
    return res
    #

def epoch_from_path(weights_path):
    try:
        fn = os.path.basename(weights_path)
        bn, _ = os.path.splitext(fn)
        epoch = int(bn.split('_')[-1])
        assert epoch >= 0
        return epoch
    except Exception:
        return 0
    #
    
def backup_path(fpath):
    if os.path.isfile(fpath):
        bpath = fpath+'.bak'
        if os.path.isfile(bpath):
            os.remove(bpath)
        os.rename(fpath, bpath)
        return bpath
    return None
    #

############################################################
#  Dataset
############################################################

def read_dataset(dataset_dir, val_split=0.75):
    all_data = {}
    for fn in os.listdir(dataset_dir):
        if not fn.endswith('_annotations_via.json'):
            continue
        jpath = os.path.join(dataset_dir, fn)
        try:
            with open(jpath, 'r') as fi:
                annotations = json.load(fi)
            for k, a in annotations.items():
                if a['regions']:
                    all_data[k] = a
            print ('Read:', jpath)
        except Exception:
            raise
            continue
    klist = sorted(all_data.keys())
    random.seed(1234)
    random.shuffle(klist)
    random.seed()
    val_idx = int(len(klist)*val_split)
    val_data = {}
    for k in klist[val_idx:]:
        val_data[k] = all_data.pop(k)
    print ('Train dataset size:', len(all_data), ' Validation dataset size:', len(val_data))
    return all_data, val_data

class RpeDataset(utils.Dataset):
    def __init__(self, dataset_dir, channel, subset, annotations):
        assert subset in ('train', 'val')
        super(RpeDataset, self).__init__()
        #
        self.dataset_dir = dataset_dir
        self.channel = channel
        self.subset = subset
        self.annotations = annotations
        #
        self.add_class(self.channel, 1, self.channel)
        #
        for a in annotations.values():
            if type(a['regions']) is dict:
                polygons = [r['shape_attributes'] for r in a['regions'].values()]
            else:
                polygons = [r['shape_attributes'] for r in a['regions']]
            if len(polygons) == 0: continue
            image_path = os.path.join(dataset_dir, a['filename'])
            try:
                fa = a['file_attributes']
                height = int(fa['height'])
                width = int(fa['width'])
            except Exception:
                image = skimage.io.imread(image_path)
                height, width = image.shape[:2]
            #
            self.add_image(
                self.channel,
                image_id=a['filename'],  # use file name as a unique image id
                path=image_path,
                width=width, height=height,
                polygons=polygons)
                
        #
        self.prepare()
    #
    def load_mask(self, image_id):
        """Generate instance masks for an image.
        Returns:
        masks: A byte array of shape [height, width, instance count] with
            one mask per instance.
        class_ids: a 1D array of class IDs of the instance masks.
        """
        # If not a balloon dataset image, delegate to parent class.
        image_info = self.image_info[image_id]
        if image_info["source"] != self.channel:
            return super(self.__class__, self).load_mask(image_id)

        # Convert polygons to a bitmap mask of shape
        # [height, width, instance_count]
        info = self.image_info[image_id]
        mask = np.zeros([info["height"], info["width"], len(info["polygons"])],
                        dtype=np.uint8)
        for i, p in enumerate(info["polygons"]):
            # Get indexes of pixels inside the polygon and set them to 1
            rr, cc = skimage.draw.polygon(p['all_points_y'], p['all_points_x'])
            mask[rr, cc, i] = 1

        # Return mask, and array of class IDs of each instance. Since we have
        # one class ID only, we return an array of 1s
        return mask, np.ones([mask.shape[-1]], dtype=np.int32)
    #
    def image_reference(self, image_id):
        """Return the path of the image."""
        info = self.image_info[image_id]
        if info["source"] == self.channel:
            return info["path"]
        else:
            super(self.__class__, self).image_reference(image_id)
    #
    def load_image(self, image_id):
        """Load the specified image and return a [H,W,3] Numpy array.
        """
        # Load image
        image = skimage.io.imread(self.image_info[image_id]['path'])
#         if image.ndim == 2:
#             image.reshape([image.shape[0], image.shape[1], 1])
        # If grayscale. Convert to RGB for consistency.
        if image.ndim != 3:
            image = skimage.color.gray2rgb(image)
        # If has an alpha channel, remove it for consistency
        if image.shape[-1] == 4:
            image = image[..., :3]
        return image

############################################################
#  Training
############################################################

def train(model, config, train_data, val_data):
    """Train the model."""
    
    channel = config.NAME
    # Training dataset.
    dataset_train = RpeDataset(dataset_dir, channel, 'train', train_data)
    # Validation dataset
    dataset_val = RpeDataset(dataset_dir, channel, 'val', val_data)

    # Image augmentation
    # http://imgaug.readthedocs.io/en/latest/source/augmenters.html
    augmentation = iaa.SomeOf((0, 2), [
        iaa.Fliplr(0.5),
        iaa.Flipud(0.5),
        iaa.OneOf([iaa.Affine(rotate=90),
                   iaa.Affine(rotate=180),
                   iaa.Affine(rotate=270)]),
    ])

    import keras
    class LfCallback(keras.callbacks.Callback):
        def on_batch_end(self, batch, logs=None):
            print('')

    # print("Train network heads")
    print ('Train %s at LR = %f, LM = %f' % (config.LAYERS_TO_TRAIN, config.LEARNING_RATE, config.LEARNING_MOMENTUM))
    model.train(dataset_train, dataset_val,
                learning_rate=config.LEARNING_RATE,
                epochs=model.epoch+config.NO_OF_EPOCHS,
                augmentation=augmentation,
                custom_callbacks=[LfCallback()],
                layers=config.LAYERS_TO_TRAIN)
    #

def parse_resp_file(resp_path):
    global DEFAULT_DATA_DIR
    fpath = os.path.abspath(resp_path)
    with open(fpath, 'rt') as fi:
        lines = fi.read().split('\n')
    arglist = []
    for _line in lines:
        line = _line.strip()
        if len(line) == 0 or line.startswith('#'):
            continue
        idx = line.find('=')
        if idx < 0:
            idx = line.find(' ')
        if idx < 0:
            arglist.append(line)
            continue
        arglist.append(line[:idx])
        arglist.append(line[idx+1:])
    DEFAULT_DATA_DIR = os.path.dirname(fpath)
    return arglist

if __name__ == '__main__':
    
    if len(sys.argv) > 1 and sys.argv[1].startswith('@'):
        try:
            arglist = parse_resp_file(sys.argv[1][1:])
        except Exception as ex:
            print ('Error parsing response file:', str(ex))
            sys.exit(1)
    else:
        arglist = sys.argv[1:]
    
    import argparse

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Train Mask R-CNN to segment an RPE channel.')
    parser.add_argument('channel', nargs=1,
            metavar='RPE_CHANNEL',
            help='RPE Microscope Channel (DNA, Actin, etc.)')
    parser.add_argument('dataset', nargs=1,
            metavar='/path/to/train/data',
            help='RPE Train Data Directory (same as in "Export Mask_RCNN"), can be relative to {data-dir}')
    parser.add_argument('-d', '--data-dir', required=False,
            metavar="/data/directory",
            default=DEFAULT_DATA_DIR,
            help='Directory where to look for data: model_weights, logs, etc. Default: "%s"' % \
                DEFAULT_DATA_DIR)
    parser.add_argument('-l', '--logs', required=False,
            metavar="/path/to/logs",
            default='logs',
            help='Logs and checkpoints directory, can be relative to {data-dir}. Default: "logs/".')
    parser.add_argument('-w', '--weights', required=False,
            metavar="/path/to/weights",
            default='model_weights',
            help='Path to weights .h5 file, directory containing .h5, or "coco".\n'+\
                'Can be relative to {data-dir}. Default: "model_weights/".')
    parser.add_argument('-W', '--workers', required=False, type=int,
            metavar="<max_workers>",
            default=0,
            help="If >0, limit the number of multiprocessing workers Mask_RCNN can create.\n"+\
                "May be needed on shared systems with limits on threads/processes.")
    parser.add_argument('-g', '--gpu-count', required=False, type=int,
            metavar="1..4",
            default=1,
            help="Number of GPUs to use (1..4); CPU only: 1, default: 1")
    parser.add_argument('-i', '--gpu-images', required=False, type=int,
            metavar="IMAGES_PER_GPU",
            default=RpeTrainConfig.IMAGES_PER_GPU,
            help="Number of images per GPU (default: %d, each image requires ~4GB GPU memory)" \
                % RpeTrainConfig.IMAGES_PER_GPU)
    parser.add_argument('-t', '--train-layers', required=False,
            metavar="LAYERS_TO_TRAIN",
            default=RpeTrainConfig.LAYERS_TO_TRAIN,
            help="Model layers to train: 'heads', 'all', '3+', '4+', '5+' (default: %s)" \
                % RpeTrainConfig.LAYERS_TO_TRAIN)
    parser.add_argument('-e', '--epochs', required=False, type=int,
            metavar="NO_OF_EPOCHS",
            default=RpeTrainConfig.NO_OF_EPOCHS,
            help="Number of training epochs (default: %d)" % RpeTrainConfig.NO_OF_EPOCHS)
    parser.add_argument('-r', '--rate', required=False, type=float,
            metavar="LEANING_RATE",
            default=RpeTrainConfig.LEARNING_RATE,
            help="Learning rate (default: %f)" % RpeTrainConfig.LEARNING_RATE)
    parser.add_argument('-m', '--momentum', required=False, type=float,
            metavar="LEARNING_MOMENTUM",
            default=RpeTrainConfig.LEARNING_MOMENTUM,
            help="Learning momentum (default: %f)" % RpeTrainConfig.LEARNING_MOMENTUM)
    parser.add_argument('-D', '--disable-gpu', required=False, action="store_true",
            help="Disable GPU(s) (this script will crash if GPUs have insufficient memory).")

    args = parser.parse_args(arglist)
    args.data_dir = os.path.abspath(os.path.join(DEFAULT_DATA_DIR, args.data_dir))
    
#     print(args)
#     sys.exit(0)
    
    if args.disable_gpu:
        print ('Disabling GPUs (if any).')
        os.environ['CUDA_VISIBLE_DEVICES']='-1'

    RpeTrainConfig.NAME = args.channel[0]
    # assert RpeTrainConfig.NAME in ('DNA', 'Actin')
    
    args.data_dir = os.path.abspath(args.data_dir)
    args.logs = os.path.abspath(os.path.join(args.data_dir, args.logs))

    dataset_dir = os.path.join(args.data_dir, args.dataset[0], 'Mask_RCNN', RpeTrainConfig.NAME)
    dataset_dir = os.path.abspath(dataset_dir)
    if not os.path.isdir(dataset_dir):
        msg = '''ERROR: "/path/to/train/data" parameter must point to a directory containing
        <..>/Mask_RCNN/%s/*.tif|*.png and *_annotations_via.json''' % (RpeTrainConfig.NAME,)
        print (msg)
        sys.exit(1)

    RpeTrainConfig.INITIAL_EPOCH = 0
    # Select weights file to load
    if args.weights.lower() == "coco":
        weights_path = os.path.join(args.data_dir, 'model_weights', 'mask_rcnn_coco.h5')
        # Download weights file
        if not os.path.exists(weights_path):
            print ('Downloading COCO weights ->', weights_path)
            utils.download_trained_weights(weights_path)
    else:
        args.weights = os.path.join(args.data_dir, args.weights)
        weights_path = find_weights_path(args.weights, RpeTrainConfig.NAME)
        RpeTrainConfig.INITIAL_EPOCH = epoch_from_path(weights_path)
    weights_dir = os.path.dirname(weights_path)
        
    logs_path = args.logs
    if not os.path.isdir(logs_path):
        try:
            os.mkdir(logs_path)
        except Exception:
            print ('ERROR: No access to the logs directory %s' % (logs_path,))
            sys.exit(1)
    
    RpeTrainConfig.MAX_WORKERS = int(args.workers)
    RpeTrainConfig.IMAGES_PER_GPU = int(args.gpu_images)
    assert RpeTrainConfig.IMAGES_PER_GPU > 0, 'Images per GPU must be >= 1'
    RpeTrainConfig.LAYERS_TO_TRAIN = args.train_layers
    RpeTrainConfig.NO_OF_EPOCHS = args.epochs
    RpeTrainConfig.LEARNING_RATE = args.rate
    RpeTrainConfig.LEARNING_MOMENTUM = args.momentum
    
    train_data, val_data = read_dataset(dataset_dir, val_split=RpeTrainConfig.TRAIN_VAL_SPLIT)
    batch_size = RpeTrainConfig.GPU_COUNT * RpeTrainConfig.IMAGES_PER_GPU
    RpeTrainConfig.STEPS_PER_EPOCH = len(train_data) // batch_size
    assert RpeTrainConfig.STEPS_PER_EPOCH > 4, 'Train dataset too small or empty'
    RpeTrainConfig.VALIDATION_STEPS = len(val_data) // batch_size
    print ('Training steps per epoch: %d; Validation steps: %d' % (RpeTrainConfig.STEPS_PER_EPOCH, RpeTrainConfig.VALIDATION_STEPS))
    
    config = RpeTrainConfig()
    # Create model
    model = modellib.MaskRCNN(mode="training", config=config,
                              model_dir=logs_path)
    if args.weights.lower() == "coco":
        print('Load COCO weights from:', weights_path, '...')
        # Exclude the last layers because they require a matching
        # number of classes
        model.load_weights(weights_path, by_name=True, exclude=[
            "mrcnn_class_logits", "mrcnn_bbox_fc",
            "mrcnn_bbox", "mrcnn_mask"])
    else:
        print('Load weights from:', weights_path, '...')
        model.load_weights(weights_path, by_name=True)
    model.epoch = RpeTrainConfig.INITIAL_EPOCH
    
    print ('Training epochs %d .. %d.' % (model.epoch+1, model.epoch+config.NO_OF_EPOCHS))
    
    train(model, config, train_data, val_data)
    
    checkpoint_path = model.checkpoint_path.format(epoch=model.epoch)
    print ('Checkpoint path:', checkpoint_path)
    if os.path.isfile(checkpoint_path):
        model_fn = os.path.basename(checkpoint_path)
        model_path = os.path.join(weights_dir, model_fn)
        print ('Copying model weights from:', checkpoint_path, 'to:', model_path)
        backup_path(model_path)
        shutil.copyfile(checkpoint_path, model_path)
        #print ('Saving model weights to:', model_path)
        #model.keras_model.save_weights(model_path, overwrite=True)

    print ('All good, exiting(0).')
    sys.exit(0)
