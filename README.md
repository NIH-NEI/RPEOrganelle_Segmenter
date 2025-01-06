## RPEOrganelle_Segmenter

For Segmentation of GFP channels in confocal stacks from the RPE-map project. For segmentation of Actin and DNA
channels, code is reused
from [RPE_Segmentation](https://github.com/NIH-NEI/RPE_Segmentation) repository. It is strongly recommended to use it
for Actin and DNA segmentation for the latest version.

Also includes Human in the loop (HITL) process for parameter selection to improve segmentations. The first set of
parameters should be found
by approximation based on number of pixels.

## Instructions

## System requirements and Installation.

Requirements depend upon the size of data being used. It is recommended (but not necessary) that you have at least 16 Gb
ram and a few hundred Gb of free space to accommodate the data. Since segmented files are usually smaller in size when
compared to original intensity images, a rule of thumb is to keep the same amount of disk space free as that occupied by
the intensity images.

1. Download and install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
   or [Anaconda](https://www.anaconda.com/products/individual).


2. Clone **RPEOrganelle_Segmenter** to a local directory `<path>/RPEOrganelle_Segmenter`.
   (Where `<path>` is the local directory of your choosingRPEOrganelle_Segmenter). You can do this by navigating to the
   location and using the
   following command in your command prompt or terminal:

   ```
   git clone https://github.com/NIH-NEI/RPEOrganelle_Segmenter
   ```

<ul>
 Alternatively, you may simply download zip under the code button on the same webpage.
</ul>

3. Run Anaconda Prompt/Terminal or in your IDE, cd to `<path>/RPEOrganelle_Segmenter`.


4. Create Conda Virtual Environment (do this once on the first run):

   ```
   conda env create -f conda-environment.yml
   ```

5. Activate the Virtual Environment:

   ```
   conda activate RPEOrganelle_Segmenter
   ```

## Data

Data used for the experiment can be found at [Deepzoomweb: RPEmap](https://isg.nist.gov/deepzoomweb/data/RPEmap). For
this repository, download the 'intensity z-stack images' hosted in the above link. Create a main folder, and create a
folder for each cell line. Inside each folder should be subfolders for different weeks. The code will go through each
subfolder, segment each file and save them in the specified save directory.

# Segmentation of gfp channels

Given paths for stacks of a channel, and location for saving files, run gfpsegment.py

```
(<path-to-environment>/)python <path-to_project>/gfpsegment.py --channelname <channelname> --path_stackfolders <path-to-root-directory> --savedir <path-to-root-save-directory. 
```

For channelname, use the appropriate value from one of the following channels. Using a channelname will automatically
choose the optimized segmentation function with the final selected parameters. Note that these parameters and
segmentation functions are not universal, but may serve as a good starting point if your magnification and resolution is
similar. For 'path_stackfolder', use the folder for each cell line, which contains subfolders for each week (or you may
choose other variable).


> Supported Channels:
>> tom20, pxn, sec61b, tuba1b, lmnb1, fbl, actb, dsp, lamp1, tjp1, myh10, st6gal1, lc3b, cetn2, slc25a17, rab5, gja1,
> > ctnnb1, (dna, actin).

To access help information for the function, use:

```
python gfpsegment.py --help
```

## Contributors

* *Pushkar Sathe*
* *Andrei Volkov* (RPE_segmentation repository, version compatibility, API access from external Python projects)
