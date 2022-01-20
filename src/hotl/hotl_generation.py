import os
from os import listdir
from os.path import isfile, join
from src.stackio import Channel
from src.segmentation.segment_GFP import segmentlaminstacks, segmentlampstacks, segmentsec61tacks, \
    segmenttomstacks, segmentstgal, segmentfbl, segmentmyh, \
    segmentrab5, segmenttub, segmentdsp, segmentpxn, segmentslc, segmentactb, segmentcetn2, \
    segmentctnnb, segmentgja, segmentlc3b

maindirpath = "C:/Users/satheps/PycharmProjects/Results/2021/Oct8/Stacks/"
savedirpath = "C:/Users/satheps/PycharmProjects/Results/2021/Dec3/Stacks_seg/"

maindirs = listdir(maindirpath)
savedirs = listdir(savedirpath)

temp = Channel

# assert maindirs == savedirs

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
    #     print("PNAME", paramname)
    if isinstance(paramname, list):
        if len(paramname) == 2:
            print("list of params found")
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
    "LC3B": "LC3B"
}

usefunction = {
    "LaminB": segmentlaminstacks,
    "LAMP1": segmentlampstacks,
    "Sec61": segmentsec61tacks,
    "TOM20": segmenttomstacks,
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
    "LC3B": segmentlc3b
}

## PHASE 1
# cutoffparams = [0.01, 0.05, 0.1] ## PHASE 1
# meanpixels ={
#     "LaminB": 3,
#     "LAMP1": 3.5,
#     "Sec61": 3,
#     "TOM20": 2.75,
#     "ST6GAL1": 3.75,
#     "FBL": 7.5,
#     "myosin": 5,
#     "RAB5": 3.25,
#     "TUBA": 3,
#     "DSP": 4.5,
#     "SLC": 2.75,
# #     "PXN": , needs checking by Davide
#     "GJA1": 8.5, #use membrane?
#     "CTNNB": 3.25,#use membrane
#     "ACTB": 5, # use observation
#     "CETN2": 4.5,
#     "LC3B": 4.5
# }
## PHASE 2
meanpixels = {
    #     "LaminB": , # ignored for phase 2
    "LAMP1": 1.75,
    "Sec61": 1.50,
    "TOM20": 2.75,
    "ST6GAL1": 3.75,
    "FBL": 7.5,
    "myosin": 2.5,
    "RAB5": 4.875,
    "TUBA": 4.5,
    "DSP": 2.25,
    "SLC": 2.75,
    #     "PXN": , needs checking by Davide
    "GJA1": 4.25,  # use membrane?
    "CTNNB": 1.625,  # use membrane
    "ACTB": 2.5,  # use observation
    "CETN2": 2.25,
    "LC3B": 4.5
}

cutoffs = {
    #     "LaminB": , # ignored for phase 2
    "LAMP1": [0.03, 0.05, 0.07],
    "Sec61": [0.05, 0.075, 0.1],
    "TOM20": [0.075, 0.1, 0.125],
    "ST6GAL1": [0.1, 0.15, 0.2],
    "FBL": [0.03, 0.05, 0.07],
    "myosin": [0.01, 0.03, 0.05],
    "RAB5": [0.05, 0.075, 0.1],
    "TUBA": [0.01, 0.03, 0.05],
    "DSP": [0.01, 0.03, 0.05],
    "SLC": [0.05, 0.075, 0.1],
    #     "PXN": , needs checking by Davide
    "GJA1": [0.01, 0.03, 0.05],  # use membrane?
    "CTNNB": [0.01, 0.03, 0.05],  # use membrane
    "ACTB": [0.01, 0.03, 0.05],  # use observation
    "CETN2": [0.05, 0.075, 0.01],
    "LC3B": [0.075, 0.1, 0.125]
}

# f2_param1=[[0.1, 0.01], [0.2, 0.01], [0.4, 0.01], [0.8, 0.01], [1.6, 0.01]]
# f2_param2=[[0.1, 0.01], [0.25, 0.01], [0.5, 0.01], [1, 0.01], [1.5, 0.01], [2, 0.2], [3,0.5]]
# f2_param3=[[0.5, 0.01]]
# laminf2params = [f2_param1,f2_param2,f2_param3]
# closingstates = [ True, False]
# laminparams = [{"f2params": laminf2, "useclosing":closingstate} for laminf2 in laminf2params for closingstate in closingstates]

# s2_param1 = [[4, 0.12], [2, 0.09], [1, 0.02]]
# s2_param2 = [[5,0.09], [2.5,0.07], [1,0.01]]
# s2_param3 = [[2.5,0.09], [1.25,0.07], [0.5,0.01]]
# lamps2s = [s2_param1, s2_param2, s2_param3]
# lampf2scales = [0.5, 0.75, 1]
# lampf2cutoffs = [0.15]
# lampparams = [{"s2params":sparam,"f2params":[[scale,cutoff]]} for sparam in lamps2s for scale in lampf2scales for cutoff in lampf2cutoffs]


# stgaltopo =[1.6, 0.8]
# stgals3scale = [0.8, 1.2, 1.6]
# stgals3cutoff = [0.02, 0.2]
# stgalparams = [{"topothin":topo ,"s3params":[[scale,cutoff]]} for cutoff in stgals3cutoff for scale in stgals3scale for topo in stgaltopo]

paramtypes = {
    "LaminB": ["f2params", "useclosing"],
    "LAMP1": "s2params",  # dont use filament? change lamp code
    "Sec61": "f2params",
    "TOM20": "f2params",
    "ST6GAL1": ["s3params", "topothin"],
    "FBL": "s2params",
    "myosin": "f3params",
    "RAB5": "s2params",
    "TUBA": "f3params",
    "DSP": "s3params",
    "SLC": "s3params",
    "PXN": "f3params",
    "GJA1": "s3params",
    "CTNNB": "s2params",
    "ACTB": "f3params",
    "CETN2": "s3params",
    "LC3B": "s3params"
}

dictofparams = {}
phase = 2
for i, dirname in enumerate(maindirs):
    if dirname == "PXN" or dirname == "LaminB":
        continue
    #     print(i, channels[dirname])
    cstates = None
    spvals = None
    if dirname == "LaminB":
        spvals = [True, False]  # useclosing
    elif dirname == "ST6GAL1":
        #         spvals = [1.6, 0.8] # phase 1
        spvals = [0.5, 0.8, 1.1]  # phase 2
    elif dirname == "PXN":
        continue
    if phase > 1:
        cutoffparams = cutoffs[dirname]  # PHASE 2 onwards
    dictofparams[dirname] = returnscalecutofflists(
        getscaleparameters(meanpixels[dirname], frac_deviation=0.25), cutoffparams,
        paramname=paramtypes[dirname], secondparamvals=spvals)

for i, dirname in enumerate(maindirs):
    print(len(dictofparams[dirname]), i, dirname, dictofparams[dirname])

for i, dirname in enumerate(maindirs):
    if dirname == "PXN" or dirname == "LaminB":
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
                print(params)
                usefunction[dirname](mainfilepath, savepath, params)
