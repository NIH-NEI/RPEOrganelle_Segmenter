import os
import sys
import datetime

DEFAULT_DATA_DIR = os.path.abspath(os.path.dirname(__file__))
ALL_CHANNELS = ['DNA', 'Actin']

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

def parse_channel_list(channel):
    if channel == 'all':
        return ALL_CHANNELS
    return [ch.strip() for ch in channel.split(',')]

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
    
    parser = argparse.ArgumentParser(
            description='Make Mask R-CNN predictions on DNA or Actin channel of RPE stacks.')
    parser.add_argument('-c', '--channel', required=False,
            metavar='DNA|Actin|all',
            default='all',
            help='RPE Microscope Channel (DNA, Actin or all), default: "all".')
    parser.add_argument('-d', '--data-dir', required=False,
            metavar="/data/directory",
            default=DEFAULT_DATA_DIR,
            help="Directory where to look for data: model_weights, logs, etc.")
    parser.add_argument('-w', '--weights', required=False,
            metavar="/path/to/weights",
            default='model_weights',
            help='Path to weights .h5 file or directory containing .h5; default: "{data-dir}/model_weights"')
    parser.add_argument('-t', '--tilesize', required=False, type=int,
            metavar="<tilesize>",
            default=768,
            help="Image tile size (must be the same as the one used for training the model)")
    parser.add_argument('-g', '--gpu-count', required=False, type=int,
            metavar="1..4",
            default=1,
            help="Number of GPUs to use (1..4); CPU only: 1, default: 1")
    parser.add_argument('-i', '--gpu-images', required=False, type=int,
            metavar="<images_per_gpu>",
            default=1,
            help="Number of images per GPU (each image requires ~4GB GPU memory)")
    parser.add_argument('-r', '--raw', required=False, action="store_true",
            help="Write raw predictions (before 3d assembly) to .raw.tif")
    parser.add_argument('-n', '--no-adjust', required=False, action="store_true",
            help="Do not perform border adjustment at post-processing")
    parser.add_argument('-R', '--recurse', required=False, action="store_true",
            help="Recurse into subdirectories (1 level)")
    parser.add_argument('-D', '--disable-gpu', required=False, action="store_true",
            help="Disable GPU(s) if this script crashes due to insufficient GPU memory.")
    parser.add_argument('rpefile', nargs='+',
            help='RPE Meta File .rpe.json, OME TIFF .ome.tif (or directory containing .rpe.json/.ome.tif)')

    args = parser.parse_args(arglist)
    args.data_dir = os.path.abspath(os.path.join(DEFAULT_DATA_DIR, args.data_dir))
    
    #assert args.channel in ('all', 'DNA', 'Actin'), 'Channel must be "DNA", "Actin" or "all".'
    assert int(args.gpu_images), 'Number of images per GPU must be > 0.'
    assert int(args.tilesize), 'Image tile size must be >= 256.'
    assert int(args.gpu_count) >= 1 and int(args.gpu_count) <= 4, 'GPU count must be between 1 and 4.'
    
    if args.weights is None:
        args.weights = 'model_weights'
    args.weights = os.path.join(args.data_dir, args.weights)
    if args.channel == 'all':
        assert os.path.isdir(args.weights), '"--weights" must be a directory if "--channel=all"'
        
    #channels = ['DNA', 'Actin', 'DNA_f3', 'Actin_f3'] if args.channel == 'all' else [args.channel]
    channels = parse_channel_list(args.channel)
        
    #print(args)
    #sys.exit(0)
    from src.RPE_Mask_RCNN_Feb22 import predictors

    if args.disable_gpu:
        args.gpu_count = 1
        print ('Disabling GPUs (if any).')
        os.environ['CUDA_VISIBLE_DEVICES']='-1'

    start_ts = datetime.datetime.now()
    n_stacks = 0
    for channel in channels:
        weights_path = find_weights_path(args.weights, channel)
        if weights_path is None:
            print ('WARNING: Could not find model weights for %s in %s' % (channel, args.weights))
            print ('Skipping', channel)
            continue
        segmenter = predictors.get_segmenter(channel, weights_path, args)
        if segmenter is None:
            print ('WARNING: No segmenter found for channel %s' % (channel,))
            print ('Skipping', channel)
            continue
        for tifstk in predictors.iter_stacks(args.data_dir, args.rpefile, channel, args.recurse):
            print('Segmenting %s from %s ...' % (channel, tifstk.fpath))
            segmenter.segment(tifstk)
            n_stacks += 1
    elapsed = str(datetime.datetime.now() - start_ts)
    print('Processed %d stack channel(s) in %s. Exiting(0).' % (n_stacks, elapsed))

    sys.exit(0)


def predict(argv):
    if len(argv) > 1 and argv[1].startswith('@'):
        try:
            arglist = parse_resp_file(sys.argv[1][1:])
        except Exception as ex:
            print('Error parsing response file:', str(ex))
            sys.exit(1)
    else:
        arglist = argv[1:]

    import argparse

    parser = argparse.ArgumentParser(
        description='Make Mask R-CNN predictions on DNA or Actin channel of RPE stacks.')
    parser.add_argument('-c', '--channel', required=False,
                        metavar='DNA|Actin|all',
                        default='all',
                        help='RPE Microscope Channel (DNA, Actin or all), default: "all".')
    parser.add_argument('-d', '--data-dir', required=False,
                        metavar="/data/directory",
                        default=DEFAULT_DATA_DIR,
                        help="Directory where to look for data: model_weights, logs, etc.")
    parser.add_argument('-w', '--weights', required=False,
                        metavar="/path/to/weights",
                        default='model_weights',
                        help='Path to weights .h5 file or directory containing .h5; default: "{data-dir}/model_weights"')
    parser.add_argument('-t', '--tilesize', required=False, type=int,
                        metavar="<tilesize>",
                        default=768,
                        help="Image tile size (must be the same as the one used for training the model)")
    parser.add_argument('-g', '--gpu-count', required=False, type=int,
                        metavar="1..4",
                        default=1,
                        help="Number of GPUs to use (1..4); CPU only: 1, default: 1")
    parser.add_argument('-i', '--gpu-images', required=False, type=int,
                        metavar="<images_per_gpu>",
                        default=1,
                        help="Number of images per GPU (each image requires ~4GB GPU memory)")
    parser.add_argument('-r', '--raw', required=False, action="store_true",
                        help="Write raw predictions (before 3d assembly) to .raw.tif")
    parser.add_argument('-n', '--no-adjust', required=False, action="store_true",
                        help="Do not perform border adjustment at post-processing")
    parser.add_argument('-R', '--recurse', required=False, action="store_true",
                        help="Recurse into subdirectories (1 level)")
    parser.add_argument('-D', '--disable-gpu', required=False, action="store_true",
                        help="Disable GPU(s) if this script crashes due to insufficient GPU memory.")
    parser.add_argument('rpefile', nargs='+',
                        help='RPE Meta File .rpe.json, OME TIFF .ome.tif (or directory containing .rpe.json/.ome.tif)')

    args = parser.parse_args(arglist)
    args.data_dir = os.path.abspath(os.path.join(DEFAULT_DATA_DIR, args.data_dir))

    # assert args.channel in ('all', 'DNA', 'Actin'), 'Channel must be "DNA", "Actin" or "all".'
    assert int(args.gpu_images), 'Number of images per GPU must be > 0.'
    assert int(args.tilesize), 'Image tile size must be >= 256.'
    assert int(args.gpu_count) >= 1 and int(args.gpu_count) <= 4, 'GPU count must be between 1 and 4.'

    print(f"rpefile: {args.rpefile},    gpu count: {args.gpu_count},    gpu images: {args.gpu_images}")

    if args.weights is None:
        args.weights = 'model_weights'
    args.weights = os.path.join(args.data_dir, args.weights)
    if args.channel == 'all':
        assert os.path.isdir(args.weights), '"--weights" must be a directory if "--channel=all"'

    # channels = ['DNA', 'Actin', 'DNA_f3', 'Actin_f3'] if args.channel == 'all' else [args.channel]
    channels = parse_channel_list(args.channel)

    # print(args)
    # sys.exit(0)

    from src.RPE_Mask_RCNN_Feb22 import predictors

    if args.disable_gpu:
        args.gpu_count = 1
        print('Disabling GPUs (if any).')
        os.environ['CUDA_VISIBLE_DEVICES'] = '-1'

    start_ts = datetime.datetime.now()
    n_stacks = 0
    for channel in channels:
        weights_path = find_weights_path(args.weights, channel)
        if weights_path is None:
            print('WARNING: Could not find model weights for %s in %s' % (channel, args.weights))
            print('Skipping', channel)
            continue
        segmenter = predictors.get_segmenter(channel, weights_path, args)
        # print(channel, channels, segmenter, segmenter is None)
        if segmenter is None:
            print('WARNING: No segmenter found for channel %s' % (channel,))
            print('Skipping', channel)
            continue
        # print(args.data_dir, args.rpefile, channel, args.recurse)

        for tifstk in predictors.iter_stacks(args.data_dir, args.rpefile, channel, args.recurse):
            print('Segmenting %s from %s ...' % (channel, tifstk.fpath))
            segmenter.segment(tifstk)
            n_stacks += 1
    elapsed = str(datetime.datetime.now() - start_ts)
    print('Processed %d stack channel(s) in %s. Exiting(0).' % (n_stacks, elapsed))

    # sys.exit(0)
