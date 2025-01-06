import os
from os import listdir
from os.path import isfile, join

from gfpseg.segmentation.segment_GFP import segmentlaminstacks, segmentlampstacks, segmentsec61tacks, \
    segmenttom, segmentstgal, segmentfbl, segmentmyh, segmentrab5, segmenttub, segmentdsp, segmentpxn, segmentslc, \
    segmentactb, segmentcetn2, segmentctnnb, segmentgja, segmentlc3b, segmenttjp
from gfpseg.stackio import Channel

# maindirpath = "../Results/../Stacks/"
maindirpath = "../Results/../Stacks_hotl/"
savedirpath = "../Results/../stgaltest/"
maindirs = listdir(maindirpath)
savedirs = listdir(savedirpath)

temp = Channel

phase = 4


def getscaleparameters(mean_pixel_dimension=None, frac_deviation=0.5):
    """
    By default frac_deviation = 0.5 for phase 1, 0.25 for phase 2, 0.125 for phase 3

    :param mean_pixel_dimension: measurement of number of pixels (e.g. for width of an organelle)
    :param frac_deviation: +/- fraction to generate new set. e.g. 0.5  gives mean*(1 +/- frac_deviation)
    :return: list of 3 values with mean and mean*(1 +/- frac_deviation)
    """
    scaleparam = mean_pixel_dimension / 3
    scaleparams = [scaleparam * (1 - frac_deviation), scaleparam, scaleparam * (1 + frac_deviation)]
    return scaleparams


def returnscalecutofflists(scales, cutoffs, paramname="f2params", secondparamvals=None):
    """
    scales: list of scales
    cutoffs: list of cutoffs
    paramname = must be f2/s2/f3/s3+params.
    """
    paramlist = None
    # print("parametar name", paramname)
    if isinstance(paramname, list):
        if len(paramname) == 2:
            print(f"list of params found: {scales}, {cutoffs}, {secondparamvals}")
            paramslist = [{paramname[0]: [[scale, cutoff]], paramname[1]: spval} for scale in scales
                          for cutoff in cutoffs for spval in secondparamvals]
    else:
        paramslist = [{paramname: [[scale, cutoff]]} for scale in scales for cutoff in cutoffs]
    return paramslist


# replace directory names with names of channels in segmented files if they are different
channels = {
    "LaminB": "LMNB",
    "LAMP1": "LAMP1",
    "Sec61": "SEC",
    "ST6GAL1": "ST6GAL1",
    "TOM20": "TOM",
    "FBL": "FBL",
    "myosin": "MYH",
    "RAB5": "RAB5",
    "TUBA": "TUBA",
    "DSP": "DSP",
    "SLC": "SLC25A17",
    "PXN": "PXN",
    "GJA1": "GJA1",
    "CTNNB": "CTNNB",
    "ACTB": "ACTB",
    "CETN2": "CETN2",
    "LC3B": "LC3B",
    "ZO1": "TJP1"
}

usefunction = {
    "LaminB": segmentlaminstacks,
    "LAMP1": segmentlampstacks,
    "Sec61": segmentsec61tacks,
    "TOM20": segmenttom,
    "ST6GAL1": segmentstgal,
    "FBL": segmentfbl,
    "myosin": segmentmyh,
    "RAB5": segmentrab5,
    "TUBA": segmenttub,
    "DSP": segmentdsp,
    "SLC": segmentslc,
    "PXN": segmentpxn,
    "GJA1": segmentgja,
    "CTNNB": segmentctnnb,
    "ACTB": segmentactb,
    "CETN2": segmentcetn2,
    "LC3B": segmentlc3b,
    "ZO1": segmenttjp
}

fracdev = 0.5 / (2 ** (phase - 1))  # for generation of values for each phase

# PHASE 1
if phase == 1:
    cutoffparams = [0.01, 0.05, 0.1]  ## PHASE 1
    meanpixels = {
        "LaminB": 3,
        "LAMP1": 3.5,
        "Sec61": 3,
        "TOM20": 2.75,
        "ST6GAL1": 3.75,
        "FBL": 7.5,
        "myosin": 5,
        "RAB5": 3.25,
        "TUBA": 3,
        "DSP": 4.5,
        "SLC": 2.75,
        #     "PXN": , discarded - data quality
        "GJA1": 8.5,  # use membrane?
        "CTNNB": 3.25,  # use membrane
        "ACTB": 5,  # use observation
        "CETN2": 4.5,
        "LC3B": 4.5,
        "ZO1": 3
    }
    stgal_topothin_pvals = [1, 2, 3]
elif phase == 2:
    ## PHASE 2
    meanpixels = {
        #     "LaminB": , # ignored for phase 2
        "LAMP1": 3.5,
        "Sec61": 1.50,
        "TOM20": 2.75,
        "ST6GAL1": 5.625,
        "FBL": 7.5,
        "myosin": 2.5,
        "RAB5": 4.875,
        "TUBA": 4.5,
        "DSP": 2.25,
        "SLC": 2.75,
        # "PXN": , discarded - data quality
        "GJA1": 4.25,  # use membrane?
        "CTNNB": 1.625,  # use membrane
        "ACTB": 2.5,  # use observation
        "CETN2": 6.75,
        "LC3B": 4.5,
        "ZO1": 1.5

    }

    ## PHASE 2
    cutoffs = {
        #     "LaminB": , # ignored for phase 2
        "LAMP1": [0.03, 0.05, 0.07],
        "Sec61": [0.05, 0.075, 0.1],
        "TOM20": [0.075, 0.1, 0.125],
        "ST6GAL1": [0.075, 0.1, 0.125],
        "FBL": [0.03, 0.05, 0.07],
        "myosin": [0.01, 0.03, 0.05],
        "RAB5": [0.05, 0.075, 0.1],
        "TUBA": [0.01, 0.03, 0.05],
        "DSP": [0.01, 0.03, 0.05],
        "SLC": [0.05, 0.075, 0.1],
        #     "PXN": ,discarded - data quality
        "GJA1": [0.01, 0.03, 0.05],
        "CTNNB": [0.03, 0.05, 0.07],  # [0.01,0.03, 0.05],
        "ACTB": [0.03, 0.05, 0.07],  # [0.01, 0.03,0.05] ,
        "CETN2": [0.05, 0.075, 0.1],
        "LC3B": [0.075, 0.1, 0.125],
        "ZO1": [0.075, 0.1, 0.125]
    }
    stgal_topothin_pvals = [1]

elif phase == 3:
    stgal_topothin_pvals = [1]

    ## PHASE 3
    meanpixels = {
        #     "LaminB": , # ignored for phase 3
        "LAMP1": 3.5,
        "Sec61": 1.875,
        "TOM20": 3.4375,
        "ST6GAL1": 5.625,
        "FBL": 7.5,
        "myosin": 1.875,  # 3.125,
        "RAB5": 6.09375,
        "TUBA": 4.5,
        "DSP": 1.6875,
        "SLC": 2.75,
        # "PXN": , discarded - data quality
        "GJA1": 4.25,  # use membrane?
        "CTNNB": 2.03125,  # 1.625, #use membrane
        "ACTB": 2.5,  # 1.875, # use observation
        "CETN2": 5.0625,
        "LC3B": 3.375,
        "ZO1": 1.5

    }

    # phase3
    cutoffs = {
        #     "LaminB": , # ignored for phase 2
        "LAMP1": [0.05, 0.07, 0.09],
        "Sec61": [0.05, 0.075, 0.1],
        "TOM20": [0.075, 0.1, 0.125],
        "ST6GAL1": [0.05, 0.075, 0.1],
        "FBL": [0.01, 0.03, 0.05],
        "myosin": [0.01, 0.02, 0.03],  # [0.01, 0.03, 0.05],
        "RAB5": [0.03, 0.05, 0.07],
        "TUBA": [0.01, 0.03, 0.05],
        "DSP": [0.01, 0.03, 0.05],
        "SLC": [0.03, 0.05, 0.07],
        #     "PXN": ,discarded - data quality
        "GJA1": [0.01, 0.03, 0.05],
        "CTNNB": [0.04, 0.05, 0.06],  # [0.01, 0.03, 0.05],
        "ACTB": [0.04, 0.05, 0.06],  # [0.01, 0.03, 0.05],
        "CETN2": [0.075, 0.1, 0.125],
        "LC3B": [0.075, 0.1, 0.125],
        "ZO1": [0.02, 0.06, 0.1]
    }

elif phase == 4:
    stgal_topothin_pvals = [0]
    meanpixels = {
        "ST6GAL1": 5.625,
        "LAMP1": 3.5,
        "Sec61": 1.875,
        "TOM20": 3.4375,
        "FBL": 7.5,
        "myosin": 1.875,  # 3.125,
        "RAB5": 6.09375,
        "TUBA": 4.5,
        "DSP": 1.6875,
        "SLC": 2.75,
        # "PXN": , discarded - data quality
        "GJA1": 4.25,  # use membrane?
        "CTNNB": 2.03125,  # 1.625, #use membrane
        "ACTB": 2.5,  # 1.875, # use observation
        "CETN2": 5.0625,
        "LC3B": 3.375,
        "ZO1": 1.5

    }
    # temp
    cutoffs = {
        #     "LaminB": , # ignored for phase 2
        "LAMP1": [0.05, 0.07, 0.09],
        "Sec61": [0.05, 0.075, 0.1],
        "TOM20": [0.075, 0.1, 0.125],
        "FBL": [0.01, 0.03, 0.05],
        "myosin": [0.01, 0.02, 0.03],  # [0.01, 0.03, 0.05],
        "RAB5": [0.03, 0.05, 0.07],
        "TUBA": [0.01, 0.03, 0.05],
        "DSP": [0.01, 0.03, 0.05],
        "SLC": [0.03, 0.05, 0.07],
        #     "PXN": ,discarded - data quality
        "GJA1": [0.01, 0.03, 0.05],
        "CTNNB": [0.04, 0.05, 0.06],  # [0.01, 0.03, 0.05],
        "ACTB": [0.04, 0.05, 0.06],  # [0.01, 0.03, 0.05],
        "CETN2": [0.075, 0.1, 0.125],
        "LC3B": [0.075, 0.1, 0.125],
        "ZO1": [0.02, 0.06, 0.1],
        "ST6GAL1": [0.075, 0.0875, 0.1],
    }
paramtypes = {
    "LaminB": ["f2params", "useclosing"],
    "LAMP1": "s2params",  # dont use filament? change lamp code
    "Sec61": "f2params",
    "TOM20": "f2params",
    "ST6GAL1": ["s3params", "topothin"],
    "FBL": "s2params",
    "myosin": "f2params",
    "RAB5": "s2params",
    "TUBA": "f3params",
    "DSP": "s3params",
    "SLC": "s3params",
    "PXN": "f3params",
    "GJA1": "s3params",
    "CTNNB": "both",
    "ACTB": "f2params",
    "CETN2": "s3params",
    "LC3B": "s3params",
    "ZO1": "f2params"
}

dictofparams = {}

for i, dirname in enumerate(maindirs):
    print(i, dirname, end="\t")
    if dirname == "PXN" or dirname == "LaminB":
        continue
    print(i, dirname)
    #     print(i, channels[dirname])
    cstates = None
    # spvals = None
    if dirname == "LaminB":
        spvals = [True, False]  # useclosing
    elif dirname == "ST6GAL1":
        pass
        spvals = stgal_topothin_pvals  # no longer using topological thin
    elif dirname == "PXN":
        continue
    if phase > 1:
        cutoffparams = cutoffs[dirname]  # PHASE 2 onwards
    dictofparams[dirname] = returnscalecutofflists(getscaleparameters(meanpixels[dirname], frac_deviation=fracdev),
                                                   cutoffparams, paramname=paramtypes[dirname], secondparamvals=spvals)
# dochannels = ["ZO1","ST6GAL1"]
dochannels = ["ST6GAL1"]
# dochannels = ["ZO1"]

for i, dirname in enumerate(maindirs):
    if dirname == "PXN" or dirname == "LaminB":  # or dirname == "ST6GAL1":
        continue
    if dirname not in dochannels:
        continue
    mainpath = join(maindirpath, dirname)
    savepath = savedirpath + "/" + dirname + "/"
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    files = [f for f in listdir(mainpath) if isfile(join(mainpath, f))]
    print(i, channels[dirname], files)
    #     break
    if i >= 0:
        for f in files:
            mainfilepath = mainpath + "/" + f
            for params in dictofparams[dirname]:
                #                 print(channels[dirname] in f)
                # print(params)
                print(params, type(params))
                usefunction[dirname](mainfilepath, savepath, params=params, minarea=0)
