from src.segmentation.segment_GFP import segmentlaminstacks, segmentlampstacks, segmentsec61tacks, \
    segmenttomstacks, segmentstgal, segmentfbl, segmentmyh, \
    segmentrab5, segmenttub, segmentdsp, segmentpxn, segmentslc, segmentactb, segmentcetn2, \
    segmentctnnb, segmentgja, segmentlc3b


class channel():
    def __init__(self, inputchannelname=None):
        self.allchannelnames = ["dna", "actin", "membrane", "tom20", "pxn", "sec61b", "tuba1b",
                                "lmnb1", "fbl", "actb", "dsp", "lamp1", "tjp1", "myh10", "st6gal1",
                                "lc3b", "cetn2", "slc25a17", "rab5", "gja1", "ctnnb1"]
        self.organellestructure = {
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
        self.rep_alphabet = {
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
        self.channelprotein = {
            "dna": "",
            "actin": "Beta-actin",
            "membrane": "",
            "tom20": "tom20",
            "pxn": "Paxillin",
            "sec61b": "",
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
            "gja1": "Connxin-43",
            "ctnnb1": "Beta-catenin"
        }
        if self.validchannelname(inputchannelname):
            self.channelname = inputchannelname
            self.channelprotein = self.getproteinname(inputchannelname)
            self.organellestructurename = self.getorganellestructurename(inputchannelname)
            self.repalphabet = self.getrepalphabet(inputchannelname)

        else:
            raise Exception(
                f"Invalid Channel name:{inputchannelname}. Name must be one of {self.allchannelnames}")
        self.usefunction = {
            "lmnb1": segmentlaminstacks,
            "lamp1": segmentlampstacks,
            "sec61b": segmentsec61tacks,
            "tom20": segmenttomstacks,
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
            "tjp1": 'not assigned yet',
        }

        self.minarea = {
            "dna": 4,
            "actin": 4,
            "membrane": 4,
            "tom20": 4,
            "pxn": 4,
            "sec61b": 4,
            "tuba1b": 4,
            "lmnb1": 4,
            "fbl": 4,
            "actb": 4,
            "dsp": 4,
            "lamp1": 4,
            "tjp1": 4,
            "myh10": 4,
            "st6gal1": 4,
            "lc3b": 4,
            "cetn2": 4,
            "slc25a17": 4,
            "rab5": 4,
            "gja1": 4,
            "ctnnb1": 4
        }
        self.directory = None

    def getallallchannelnames(self):
        return self.allchannelnames

    def getminarea(self, key):
        return self.minarea[key]

    def getproteinname(self, key):
        return self.channelprotein[key]

    def getorganellestructurename(self, key):
        return self.organellestructure[key]

    def validchannelname(self, key):
        return key.lower() in self.allchannelnames

    def getrepalphabet(self, key):
        return self.rep_alphabet[key]

    def setdirectoryname(self, dname):
        """
        Use to set directory name in case it is different from channelname
        :param dname:
        :return:
        """
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
        self.directory = dname
