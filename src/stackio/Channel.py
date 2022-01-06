from src.segmentation.segment_GFP import segmentlaminstacks, segmentlampstacks, segmentsec61tacks, segmenttomstacks, segmentstgal, segmentfbl, segmentmyh, \
    segmentrab5, segmenttub, segmentdsp, segmentpxn, segmentslc, segmentactb, segmentcetn2, segmentctnnb, segmentgja, segmentlc3b


class channel():
    def __init__(self, key=None):
        self.allchannelnames = ["DNA", "Actin", "Membrane", "TOM20", "PXN", "SEC61B", "TOM20", "TUBA1B", "LMNB1", "FBL", "ACTB", "DSP", "LAMP1", "TJP1", "MYH10", "ST6GAL1", "LC3B", "CETN2",
                                "SLC25A17",
                                "RAB5", "GJA1", "CTNNB1"]

        self.organellestructure = {
            "DNA": "Nucleus",  # check
            "Actin": "Actin Filaments",
            "Membrane": "Cell membrane",  # check
            "TOM20": "Mitochondria",
            "PXN": "Matrix adhesions",
            "SEC61B": "Endoplasmic reticulum",
            "TUBA1B": "Microtubules",
            "LMNB1": "Nuclear Envelope",
            "FBL": "Nucleolus",
            "ACTB": "Actin Filaments",
            "DSP": "Desmosomes",
            "LAMP1": "Lysosome",
            "TJP1": "Tight Junctions",
            "MYH10": "Actomyosin bundles",
            "ST6GAL1": "Golgi Apparatus",
            "LC3B": "Autophagosomes",
            "CETN2": "Centrioles",
            "SLC25A17": "Peroxisomes",
            "RAB5": "Endosomes",
            "GJA1": "Gap Junctions",
            "CTNNB1": "Adherens Junctions"
        }

        self.channelprotein = {
            "DNA": "",
            "Actin": "Beta-actin",
            "Membrane": "",
            "TOM20": "Tom20",
            "PXN": "Paxillin",
            "SEC61B": "",
            "TUBA1B": "Alpha Tubulin",
            "LMNB1": "Lamin B1",
            "FBL": "Fibrillarin",
            "ACTB": "Beta-actin",
            "DSP": "Desmoplakin",
            "LAMP1": "LAMP-1",
            "TJP1": "Tight Junction Protein Z-01",
            "MYH10": "Non-muscle myosin heavy chain IIB",
            "ST6GAL1": "Sialyltransferase 1",
            "LC3B": "Autophagy-related protein LC3 B",
            "CETN2": "Centrin-2",
            "SLC25A17": "Peroxisomal membrane protein",
            "RAB5": "Ras-related protein Rab-5A",
            "GJA1": "Connxin-43",
            "CTNNB1": "Beta-catenin"
        }
        if self.validchannelname(key):
            self.channelname = key
            self.channelprotein = self.getproteinname(key)
            self.organellestructurename = self.getorganellestructurename(key)
        else:
            raise Exception(f"Invalid Channel name:{key}. Name must be one of {self.allchannelnames}")
        self.usefunction = {
            "LaminB": segmentlaminstacks,
            "LAMP1": segmentlampstacks,
            "Sec61": segmentsec61tacks,
            "TOM20": segmenttomstacks,
            "ST6GAL1": segmentstgal,
            "FBL": segmentfbl,
            "MYH10": segmentmyh,
            "RAB5": segmentrab5,
            "TUBA": segmenttub,
            "DSP": segmentdsp,
            "SLC25A17": segmentslc,
            "PXN": segmentpxn,
            "GJA1": segmentgja,
            "CTNNB1": segmentctnnb,
            "ACTB": segmentactb,
            "CETN2": segmentcetn2,
            "LC3B": segmentlc3b,
            "TJP1": 'not assigned yet',
        }

    def getallallchannelnames(self):
        return self.allchannelnames

    def getproteinname(self, key):
        return self.channelprotein[key]

    def getorganellestructurename(self, key):
        return self.organellestructure[key]

    def validchannelname(self, key):
        return key in self.allchannelnames
