# Software Suite for Detecting and Segmenting RPE Cell Walls and Nuclei.
*based on [Mask R-CNN implementation by Waleed Abdulla (2017)](https://github.com/matterport/Mask_RCNN).*

This software suite is intended to train Mask R-CNN model to detect and segment RPE channels "DNA" and "Actin".

The input data for training is exported from RPE Map Segmentation tool. The format of annotations is compatible with [VGG Image Annotator (VIA)](https://gitlab.com/vgg/via), even though DNA images are 16-bit grayscale TIFFs, which can not be loaded into browser.

The input data for predictions is "metadata" file(s) `*.rpe.json`, generated by RPE Map Segmentation tool. They must contain at least 3 reference channels, **DNA**, **Actin** and **Membrane**. The DNA training/prediction part only uses the DNA channel, but the Actin part uses all three.

Directory containing data (training data set, logs, model weights, input data for predictions) is separate from the directory containing the code. By default it is the same, i.e. if you place `model_weights` directory next to the main python scripts and use absolute paths to specify input `*.rpe.json files`, you don't have to specify `--data-dir` command line option. However this is not possible when running scripts in a Docker container. Keeping all data in a single base directory separate from the code simplifies Docker's command line (you only have to mount one host directory to a container's directory).

Both training and prediction scripts accept response file instead of the command line. A response file is a text file that lists all the desired command line arguments, one argument per line. A few examples of response files are supplied in the sample data distribution. When the first command line argument starts with **@**, it is treated as the path to a response file (without @), which is parsed instead of the command line, and the rest of the actual command line (if any) is ignored. When a response file is given, the default base directory is changed to the directory containing the response file. This can be changed by specifying `--base-dir /another/base/dir` within the response file.


## Requirements

This software suite requires Python 3.6, TensorFlow 1.15 and Keras 2.1. To take advantage of a GPU it also needs CUDA 10.0 and CuDNN 7.6.5. The prediction script uses a native code extension, which is pre-compiled for Windows. On Linux system it can be easily built, but make sure the **swig** package is installed in the system. To (re)compile it under Windows, **swig** must be installed and accessible via system's **PATH**, and also a C++ compiler (such as MS Visual C++) is available. The Python dependencies are listed in [requirements.txt](./requirements.txt).

Since this suite has a very specific set of requirements, it is generally a good idea to run it in an isolated environment, such as Miniconda VE, Docker container or locally installed/sourced CUDA/CuDNN and Python/dependencies. See below for step-by-step instructions in each case.


## Entry points

[train.py](./train.py) -- training script. The input is a training data directory exported by RPE Map Segmentation tool. The output is `mask_rcnn_dna_NNNN.h5` and `mask_rcnn_actin_NNNN.h5` files in the `model_weights` directory. Note that the output files are the last versions of the model weights, available only if the script finishes normally. If it is aborted, no new files are written into `model_weights`, but you can find the intermediate versions (at check points after each completed epoch) in the `logs` directory. If you want to continue training after an abnormal termination, you need to copy these files manually from `logs/XXXXXXXXX` to `model_weights`.

```
# Examples:

python train.py DNA RPE_Training --data-dir C:\rpemrcnn --epochs 15 --rate 0.0005
#     Train the DNA model using dataset in "C:\rpemrcnn\RPE_Training"
# for 15 epochs at learning rate 0.0005.

python train.py Actin C:\RPE_Training --weights coco --train-layers heads --epochs 5
#     Train the Actin model head layers starting from the MS COCO weights
# (download if necessary) for 5 epochs using dataset in C:\RPE_Training.
# The "logs" and "model_weights" directories are in the same directory as "train.py".
```

[predict.py](./predict.py) -- prediction script. The input is one or more `*.rpe.json` meta-files or directories containing these meta-files. The output is pairs of `*_RPE.csv`/`*_RPE.tif` files stored in the `Predicted` directory next to the corresponding `*.rpe.json`.

```
# Examples:

python predict.py C:\rpemrcnn\W2\P1-W2-ZO1_D02_F006.rpe.json
#     Segment both DNA and Actin of "C:\rpemrcnn\W2\P1-W2-ZO1_D02_F006.rpe.json"
# using latest model weights from ".\model_weights" directory.

python predict.py W1\P1-W1-ZO1_D04_F004.rpe.json --data-dir C:\rpemrcnn --channel DNA
#     Segment the DNA channel of stack "C:\rpemrcnn\W1\P1-W1-ZO1_D04_F004.rpe.json"
# using latest DNA model weights from "C:\rpemrcnn\model_weights\mask_rcnn_dna_NNNN.h5".
# The output is "C:\rpemrcnn\W1\Predicted\P1-W1-ZO1_D04_F004_DNA_RPE.csv" and
# "C:\rpemrcnn\W1\Predicted\P1-W1-ZO1_D04_F004_DNA_RPE.tif".

python predict.py @C:\rpemrcnn\predict.txt
#     Segment whatever is specified in the response file "C:\rpemrcnn\predict.txt"
# with whatever options are present in the same file.
```


## Miniconda Virtual Environment (Windows)

1. Download and install [Miniconda or Anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html).
2. Check out (or unzip) *RPE_Mask_RCNN* to a directory, e.g. `C:\ML\RPE_Mask_RCNN`. Run *Anaconda Prompt*, cd to this directory.
3. Create virtual environment (do this once, next time skip to the next step):

	`conda env create --file conda-environment.yml`

4. Activate virtual environment:

	`conda activate Mask_RCNN`

5. (optional) Run Python, and at the prompt type `import rpesegm`. If there are no error messages, the native extension is OK, you can skip the next step. If not, you need to build the native extension.
6. To build the native extension you will need [Microsoft Visual C++](https://www.microsoft.com/en-us/download/developer-tools.aspx) and [Swig](http://www.swig.org). Swig must be accessible via system's `PATH`. If both tools are installed in your system, simply type:

	`build_swig.bat`

7. Now you can run training and prediction scripts (see examples above). For a quick help on command line arguments, type:

	`python train.py --help`

	`python predict.py --help`

To delete the Virtual environment, deactivate it (if it is activated):

`conda deactivate`

then type:

`conda env remove --name Mask_RCNN`
   

## Docker container

1. Download and install [Docker or Docker Desktop](https://www.docker.com/).

2. Check out (or unzip) *RPE_Mask_RCNN* to a directory, e.g. `C:\ML\RPE_Mask_RCNN`, open Command Prompt, cd to this directory.

3. Create a Docker image (replace *rpe-mask-rcnn-cuda-10* with any other name if you want):

	`docker build . -t rpe-mask-rcnn-cuda-10`
	
4. Run scripts in the container, assuming your data files (model weights, training data, source data, response files, etc., are in the directory `C:\rpemrcnn`:

```
docker run -it -v "C:\\rpemrcnn:/rpemrcnn" --rm rpe-mask-rcnn-cuda-10 predict.py -d /rpemrcnn W1
# Note that "-v" option maps host directory "C:\rpemrcnn" to container directory "/rpemrcnn",
# this way making all the data inside this directory accessible to the ML scripts
# running within the container.
```

`docker run -it -v "C:\\rpemrcnn:/rpemrcnn" --rm rpe-mask-rcnn-cuda-10 train.py @/rpemrcnn/train.txt`

5. Export Docker image to a `*.tar.gz`. On Windows system, you will need to download [gzip](https://www.gnu.org/software/gzip/) and put it in a directory accessible via system PATH environment variable.

	`docker save -o rpe-mask-rcnn-cuda-10.tar rpe-mask-rcnn-cuda-10`
	
	`gzip rpe-mask-rcnn-cuda-10.tar`

6. Import Docker image (copied from another machine):

	`gzip -d rpe-mask-rcnn-cuda-10.tar.gz`
	
	`docker load -i rpe-mask-rcnn-cuda-10.tar`


## NIH HPC Biowulf

Biowulf is a 105,000+ processor Linux cluster managed by [NIH HPC Group](https://hpc.nih.gov/). It is accessible only within NIH intranet, including VPN. To run anything on Biowulf you need a [NIH HPC account](https://hpc.nih.gov/docs/accounts.html).

Software needed to connect to Biowulf:
* SFTP client, e.g. [WinSCP](https://winscp.net/eng/download.php); connect to `helix.nih.gov`.
* [NX No-Machine client](https://hpc.nih.gov/docs/nx.html); connect to `biowulf.nih.gov`.

The home directory `/home/$USER` in Biowulf has a very limited capacity, good for storing configurations and settings, but nothing else. Anything bigger must go to `/data/$USER`. You can create symbolic links in `/home/$USER/xxx` pointing to `/data/$USER/xxx` if you need something accessible via the **~** directory.

There is a useful [Biowulf User Guide](https://hpc.nih.gov/docs/userguide.html) about how to run stuff in Biowulf.

Biowulf has various versions of CUDA/CuDNN libraries, as well as various Python versions, that can be sourced via the `module` command. Unfortunately, **TensorFlow 1** is not supported, only **TensorFlow 2**, and creating an isolated environment (free from Python package conflicts) for Mask R-CNN can be a challenge. One way to handle this is to build Python 3.6 from sources locally.

### Build clean local Python 3.6 environment.

1. Start `NX`, open Terminal, `cd /data/$USER`

2. `mkdir python`

3. `cd python`

4. `wget https://www.python.org/ftp/python/3.6.14/Python-3.6.14.tgz`

5. `tar xzvf Python-3.6.14.tgz`

6. `find /data/$USER/python -type d|xargs chmod 0755`

7. `cd Python-3.6.14/`

8. `./configure --prefix=/data/$USER/python`

9. `make && make install`

10. `export PATH=/data/$USER/python/Python-3.6.14:$PATH`

This last command updates your **PATH** so that typing `python` will open the local (clean) version of Python. It needs to be repeated every time you log in or open a new Biowulf session. You can add it to your `~/bash_profile` as the last command, in which case it will be executed automatically every time you log in to Biowulf. Alternatively, you can create a separate file with a short name, e.g. `~\py36`, and source it every time you log in with the command `source ~\py36`.

11. `python -m pip install --no-cache-dir wheel`

12. `python -m pip install --upgrade pip`


### Install RPE Mask R-CNN dependencies.

1. Copy RPE_Mask_RCNN to `/data/$USER/RPE_Mask_RCNN` and rpercnn, to `/data/$USER/rpercnn`, using an SFTP client such as WinSCP.

2. `cd /data/$USER/RPE_Mask_RCNN`

3. `python -m pip install --no-cache-dir -r requirements.txt -t packages`

4. `export PYTHONPATH=/data/$USER/RPE_Mask_RCNN:/data/$USER/RPE_Mask_RCNN/packages`

This command makes Python packages installed in `RPE_Mask_RCNN/packages` available to Python scripts. It needs to be executed before running scripts using Mask_RCNN. The best place to do that is, perhaps, inside shell scripts submitted as Biowulf jobs.

5. `module load gcc/7.5.0`

6. `sh build_swig.sh`

This one builds the native extension. It is needed only for [predict.py](./predict.py), but not for [train.py](./train.py).

### Run scripts in an interactive session.

1. Start `NX`, open Terminal, `cd /data/$USER/RPE_Mask_RCNN`

2. List GPU partitions, available nodes, CPUs, GPUs: `freen | grep -E 'Partition|----|gpu'`

Note which *gpu* partition has lots of available GPUs. Let's assume it's "k80".

3. `sinteractive --gres=gpu:k80:1 --cpus-per-task=12 --mem 65536`

It may take a few minutes until the job is granted, and you are back at the command prompt. Notice the *Job ID* just displayed.

4. `newwall --jobid {Job ID} --time 1-00:00:00` -- Substitute {Job ID} with the actual Job ID assigned by the `sinteractive` command. This command extends the wall time from 8 hours (default) to 1 day.

5. `module load cuDNN/7.6.5/CUDA-10.0`

6. `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/CUDA/10.0.130/lib64` -- You may need to check the actual CUDA 10.0 version: `ls /usr/local/CUDA/`, use the one that starts with `10.0.`.

7. `export PYTHONPATH=/data/$USER/RPE_Mask_RCNN:/data/$USER/RPE_Mask_RCNN/packages`

8. Run the scripts (see examples above).

9. When done, just type `exit`. This will terminate the interactive session and cancel the job {Job ID}, but not close the terminal.


### Submit Biowulf jobs.

1. Create shell scripts, example:

*train_actin.sh*
```
#!/bin/bash
module load cuDNN/7.6.5/CUDA-10.0
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/CUDA/10.0.130/lib64
export PYTHONPATH=/data/$USER/RPE_Mask_RCNN:/data/$USER/RPE_Mask_RCNN/packages
python train.py @/data/$USER/rpemrcnn/train_actin.txt
```

Note this script uses a response file, which you can modify before submitting the job instead of modifying the script itself.

2. `chmod 0755 *.sh`

3. Modify the response file, if necessary.

4. Make sure directory `/data/$USER/rpemrcnn/logs` exists; if not, create it: `mkdir /data/$USER/rpemrcnn/logs`.

5. `sbatch --partition=gpu --cpus-per-task=12 --gres=gpu:k80:1 --mem 65536 --output=/data/$USER/rpemrcnn/logs/train_actin.log train_actin.sh`

This command displays {Job ID} right away (and nothing else), so modify the wall time in the next step.

6. `newwall --jobid 12345 --time 3-00:00:00`

The default wall time for batch jobs is 2 hours, so usually you have to extend it (up to 10-00:00:00).

7. Watch the output: `tail -f -n 20 /data/$USER/rpemrcnn/logs/train_actin.log`.

You can submit multiple jobs at the same time. There is a limit on jobs/GPUs per user, etc. -- if you exceed that limit, new jobs will be waiting in the queue until an old job finishes.

To list active jobs, `sacct`. To cancel a job, type `scancel {Job ID}`.
