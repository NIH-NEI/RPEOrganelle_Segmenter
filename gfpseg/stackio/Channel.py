import warnings

from gfpseg.segmentation.segment_GFP import segmentlaminstacks, segmentlampstacks, segmentsec61tacks, segmenttom, \
    segmentstgal, segmentfbl, segmentmyh, segmentrab5, segmenttub, segmentdsp, segmentpxn, segmentslc, segmentactb, \
    segmentcetn2, segmentctnnb, segmentgja, segmentlc3b, segmenttjp


class channel():
    allchannelnames = ["dna", "actin", "membrane", "tom20", "pxn", "sec61b", "tuba1b",
                       "lmnb1", "fbl", "actb", "dsp", "lamp1", "tjp1", "myh10", "st6gal1",
                       "lc3b", "cetn2", "slc25a17", "rab5", "gja1", "ctnnb1"]
    organellestructure = {
        "dna": "Nucleus",  # check
        "actin": "Actin Filaments",
        "membrane": "Cell membrane",  # check
        "tom20": "Mitochondria",
        "pxn": "Matrix adhesions",
        "sec61b": "Endoplasmic reticulum",
        "tuba1b": "Microtubules",
        "lmnb1": "Nuclear Envelope",
        "fbl": "Nucleolus",
        "actb": "Actin Filaments",
        "dsp": "Desmosomes",
        "lamp1": "Lysosome",
        "tjp1": "Tight Junctions",
        "myh10": "Actomyosin bundles",
        "st6gal1": "Golgi Apparatus",
        "lc3b": "Autophagosomes",
        "cetn2": "Centrioles",
        "slc25a17": "Peroxisomes",
        "rab5": "Endosomes",
        "gja1": "Gap Junctions",
        "ctnnb1": "Adherens Junctions"
    }

    channelprotein = {
        "dna": "DNA",  # "Hoechst stain(nucleus)"
        "actin": "Beta-actin",
        "membrane": "Membrane",
        "tom20": "tom20",
        "pxn": "Paxillin",
        "sec61b": "Sec61 beta",
        "tuba1b": "Alpha Tubulin",
        "lmnb1": "Lamin B1",
        "fbl": "Fibrillarin",
        "actb": "Beta-actin",
        "dsp": "Desmoplakin",
        "lamp1": "LAMP-1",
        "tjp1": "Tight Junction Protein Z-01",
        "myh10": "Non-muscle myosin heavy chain IIB",
        "st6gal1": "Sialyltransferase 1",
        "lc3b": "Autophagy-related protein LC3 B",
        "cetn2": "Centrin-2",
        "slc25a17": "Peroxisomal membrane protein",
        "rab5": "Ras-related protein Rab-5A",
        "gja1": "Connexin-43",
        "ctnnb1": "Beta-catenin"
    }

    rep_alphabet = {
        "dna": "default",
        "actin": "default",
        "membrane": "default",
        "tom20": "E",
        "pxn": "F",
        "sec61b": "G",
        "tuba1b": "C",
        "lmnb1": "F",
        "fbl": "G",
        "actb": "D",
        "dsp": "B",
        "lamp1": "B",
        "tjp1": "D",
        "myh10": "F",
        "st6gal1": "C",
        "lc3b": "C",
        "cetn2": "G",
        "slc25a17": "E",
        "rab5": "E",
        "gja1": "G",
        "ctnnb1": "F"
    }

    defaultparameters = {
        "lmnb1": {"scale": [None], "cutoff": [None], "final": False, "minarea": 2},  # FINAL
        "lamp1": {"scale": [3.0625 / 3], "cutoff": [0.05], "final": True, "minarea": 3},  # FINAL
        "sec61b": {"scale": [0.625], "cutoff": [0.075], "final": True, "minarea": 4},  # FINAL
        "tom20": {"scale": [3.4375 / 3], "cutoff": [0.12], "final": True, "minarea": 2},  # FINAL
        "st6gal1": {"scale": [1.875], "cutoff": [0.0875], "final": True, "minarea": 0},  # FINAL
        "fbl": {"scale": [2.1875], "cutoff": [0.05], "final": True, "minarea": 5},  # FINAL
        "myh10": {"scale": [1.640625 / 3], "cutoff": [0.01], "final": True, "minarea": 0},
        "rab5": {"scale": [1.77734375], "cutoff": [0.07], "final": True, "minarea": 4},  # FINAL
        "tuba1b": {"scale": [None], "cutoff": [None], "final": False, "minarea": 4},
        "dsp": {"scale": [0.5625, 1.8984375 / 3], "cutoff": [0.03, 0.01], "final": True, "minarea": 0},
        "slc25a17": {"scale": [2.75 / 3, 2.40625/3], "cutoff": [0.05, 0.03], "final": True, "minarea": 0},  # FINAL
        "pxn": {"scale": [None], "cutoff": [None], "final": False, "minarea": 4},  # NOT USED
        "gja1": {"scale": [1.59375, 1.0625], "cutoff": [0.03, 0.01], "final": True, "minarea": 0},  # FINAL
        "ctnnb1": {"scale": [1.77345 / 3], "cutoff": [0.05], "final": True, "minarea": 0},  # FINAL
        "actb": {"scale": [2.1875 / 3], "cutoff": [0.05], "final": True, "minarea": 0},  # FINAL
        "cetn2": {"scale": [1.6875], "cutoff": [0.075], "final": True, "minarea": 2},  # FINAL
        "lc3b": {"scale": [1.125], "cutoff": [0.075], "final": True, "minarea": 4},  # FINAL
        "tjp1": {"scale": [0.5], "cutoff": [0.1], "final": True, "minarea": 0},
    }
    dirnames = {
        "lmnb1": "LaminB",
        "lamp1": "LAMP1",
        "sec61b": "Sec61",
        "st6gal1": "ST6GAL1",
        "tom20": "TOM20",
        "fbl": "FBL",
        "myh10": "myosin",
        "rab5": "RAB5",
        "tuba1b": "TUBA",
        "dsp": "DSP",
        "slc25a17": "SLC",
        "pxn": "PXN",
        "gja1": "GJA1",
        "ctnnb1": "CTNNB",
        "actb": "ACTB",
        "cetn2": "CETN2",
        "lc3b": "LC3B"}
    usefunction = {
        "lmnb1": segmentlaminstacks,
        "lamp1": segmentlampstacks,
        "sec61b": segmentsec61tacks,
        "tom20": segmenttom,
        "st6gal1": segmentstgal,
        "fbl": segmentfbl,
        "myh10": segmentmyh,
        "rab5": segmentrab5,
        "tuba1b": segmenttub,
        "dsp": segmentdsp,
        "slc25a17": segmentslc,
        "pxn": segmentpxn,
        "gja1": segmentgja,
        "ctnnb1": segmentctnnb,
        "actb": segmentactb,
        "cetn2": segmentcetn2,
        "lc3b": segmentlc3b,
        "tjp1": segmenttjp,
    }

    def __init__(self, inputchannelname=None):
        ### Scales and cutoffs from phase 3 human-in-the-loop parameter selection
        if channel.validchannelname(inputchannelname):
            self.channelname = inputchannelname
            self.channelprotein = channel.getproteinname(inputchannelname)
            self.organellestructurename = channel.getorganellestructurename(inputchannelname)
            self.repalphabet = channel.getrepalphabet(inputchannelname)

        else:
            raise Exception(
                f"Invalid Channel name:{inputchannelname}. Name must be one of {self.allchannelnames}")

        self.directory = channel.dirnames[inputchannelname]

    @staticmethod
    def getallallchannelnames():
        return channel.allchannelnames

    @staticmethod
    def getminarea(channelname=None):
        if channelname is None:
            print("please input a channel name")
        return channel.defaultparameters[channelname]["minarea"]

    @staticmethod
    def getproteinname(channelname=None):
        if channelname is None:
            print("please input a channel name")
        return channel.channelprotein[channelname]

    @staticmethod
    def getorganellestructurename(channelname=None):
        if channelname is None:
            print("please input a channel name")
        return channel.organellestructure[channelname]

    @staticmethod
    def validchannelname(channelname=None):
        if channelname is None:
            print("please input a channel name")
            # channelname = self.channelname
        return channelname.lower() in channel.allchannelnames

    @staticmethod
    def getrepalphabet(channelname=None):
        if channelname is None:
            print("please input a channel name")
        return channel.rep_alphabet[channelname]

    def setdirectoryname(self, dname):
        """
        Use to set directory name in case it is different from channelname
        :param dname:
        :return:
        """
        if dname is not None:
            self.directory = dname

    @staticmethod
    def getdefaultparams(channelname):
        if channel.defaultparameters[channelname]["final"] is not True:
            warnings.warn("Note: parameters are not yet final!")
        # parameterkeys = ["scale", "cutoff", "topothin"]
        parameterkeys = ["scale", "cutoff"]
        assert len(channel.defaultparameters[channelname]["scale"]) == len(
            channel.defaultparameters[channelname]["cutoff"]), "unequal number of parameters"
        nkeys =len(channel.defaultparameters[channelname]["scale"])
        threshparams = []
        for i in range(nkeys):
            subparams = []
            for subkey in parameterkeys:
                if subkey in channel.defaultparameters[channelname].keys():
                    subparams.append(channel.defaultparameters[channelname][subkey][i])
            threshparams.append(subparams)
        # threshparams = [[channel.defaultparameters[channelname][subkey][i] for subkey in parameterkeys if
        #                  subkey in channel.defaultparameters[channelname].keys() for i in
        #                  range(len(channel.defaultparameters[channelname][subkey]))]]
        minarea = channel.defaultparameters[channelname]["minarea"]
        return threshparams, minarea

    @staticmethod
    def segmentchannel(filename, savepath, params=None, channelname=None, minarea=None, **kwargs):
        if channelname is None:
            print("Input channel name!")
        if params is None:
            ch_params, ch_minarea = channel.getdefaultparams(channelname)
            params = ch_params
            print(f"Using default parameters:{params}. To use different parameters  input them as arguments")
        if minarea is not None:
            ch_minarea = minarea
        # print(self.usefunction[self.channelname])
        return channel.usefunction[channelname](fpath=filename, savepath=savepath, params=params,
                                                channel=channelname, minarea=ch_minarea, **kwargs)
