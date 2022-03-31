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
        assert _alphabets in allalphabets, f"Alphabets must be from {allalphabets}"

        if isinstance(_alphabets,str):
            _alphabets=[_alphabets]
        # elif isinstance(_alphabets,list):
        #     assert len(_alphabets) == 1 ,"input must contain only 1 object"
        else:
            print("input alphabet must be a list or string")
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

def get_rfcw(fname,checkw=False):
    split = fname.split("_")
    if len(split)==3:
        s1, r, fc = split
    elif len(split)==4:
        s1, s2, r, fc = split
    c = fc.split(".tif")[0][-3:]
    f = fc.split(".tif")[0][-16:-12]
    # print("f,c",f,c)
    # raise Exception
    # print(fname.split("_"))
    # r_ = int(r[1:]) - 2
    if checkw:
        w = s1.split("-")[1]
        # w_ = WS.index(w)
        return r,f ,c, w
    return r,f,c


def findtreatment(r):  # TODO: check with getwr_3channel for inconsistencies
    """
    Returns the type of treatement based on replicate id
    :param r: replicate id ( converted to 0-9 range)
    :return: treatment id
    """
    assert r < 10, "error"
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

