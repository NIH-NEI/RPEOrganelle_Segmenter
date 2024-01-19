## RPEOrganelle_Segmenter

For handling confocal stacks from the RPE-map project; 
Segmentation of GFP channels as well as Actin and DNA channels are implemented. 

Human in the loop process for parameter selection to improve segmentations. The first set of parameters should be found by approximation based on number of pixels.
 


## Instructions

# 1. Install:
```
git clone https://github.com/RPEGoogleMap/RPEOrganelle_Segmenter
```

# 2. Analyze set of segmented channels:
Given paths for stacks of a channel, and location for saving files, run gfpsegment.py
```python
(<path-to-environment>/)python <path-to_project>/gfpsegment.py
```
