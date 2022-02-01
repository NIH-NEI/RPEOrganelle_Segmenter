USEDTREATMENTS = 2
USEDWEEKS = 4
TREATMENT_TYPES = ["PGE2", "HPI4"]
WS = ["W1", "W2", "W3", "W4", "W5", "W6"]

XSCALE, YSCALE, ZSCALE = 0.21666666666666667, 0.21666666666666667, 0.500
VOLUMESCALE = XSCALE * YSCALE * ZSCALE
AREASCALE = XSCALE * YSCALE
T = 1
C = 4
Z = 27
Y = 1078
X = 1278
setstackshape = (T, C, Z, X, Y)


def generate_repinfo(_alphabets=None):
    """
    Input one of the alphabets based on specific file information
    :param _alphabets:
    :return: list of replicate ids
    """
    allalphabets = ["B", "C", "D", "E", "F", "G"]
    if _alphabets is None:
        _alphabets = allalphabets
    else:
        assert (_alphabets in allalphabets, f"Alphabets must be from {allalphabets}")
    _repnums = ["02", "03", "04", "05", "06", "07", "08", "09", "10", "11"]
    reps = []
    for a in _alphabets:
        for repnum in _repnums:
            reps.append(a + repnum)
    return reps


# def getwr(df, af, lf):
#     basestringdna = "_".join(df.split("_")[:-2])
#     basestringactin = "_".join(af.split("_")[:-3])
#     basesstringlmp = "_".join(lf.split("_")[:-1])
#
#     print(basestringdna, basestringactin, basesstringlmp)
#     assert basestringdna == basestringactin == basesstringlmp
#     s1, r, _ = basestringdna.split("_")
#     w = s1.split("-")[1]
#     w_ = WS.index(w)
#     r_ = int(r[1:]) - 2
#     return w, r, w_, r_, basestringdna

def getwrprestack(fname):  # TODO: test, implement remaining part
    s1, r, _ = fname.split("_")
    w = s1.split("-")[1]
    w_ = WS.index(w)
    r_ = int(r[1:]) - 2
    return w, r, w_, r_, fname


def findtreatment(r):  # TODO: check with getwr_3channel for inconsistencies
    """
    Returns the type of treatement based on replicate id
    :param r: replicate id ( converted to 0-9 range)
    :return: treatment id
    """
    assert (r < 10, "error")
    treatment = None
    if r < 5:
        treatment = 0
    else:
        treatment = 1
    return treatment

#
# def checkcellconditions(cellvals):
#     [centroid, vol, xspan, yspan, zspan, maxferet, minferet] = cellvals
#     satisfied_conditions = True
#     if (zspan < 1) or (xspan < 2) or (yspan < 2):
#         satisfied_conditions = False
#     if vol >= 100000 or vol <= 50:  # 50 = 2130,in cu micron
#         satisfied_conditions = False
#     return satisfied_conditions
#
minarea = {
            "dna": 4,
            "actin": 4,
            "membrane": 4,
            "tom20": 2,
            "pxn": 4,
            "sec61b": 4,
            "tuba1b": 4,
            "lmnb1": 4,
            "fbl": 5,
            "actb": 4,
            "dsp": 4,
            "lamp1": 3,
            "tjp1": 4,
            "myh10": 4,
            "st6gal1": 4,
            "lc3b": 4,
            "cetn2": 4,
            "slc25a17": 4,
            "rab5": 4,
            "gja1": 5,
            "ctnnb1": 5
        }
