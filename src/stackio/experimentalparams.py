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

    Args:
        _alphabets: alphabets that occur in the _ middle _ part of filename
    Returns:
        list of replicate ids
    """
    allalphabets = ["B", "C", "D", "E", "F", "G"]
    if _alphabets is None:
        _alphabets = allalphabets
    else:
        assert _alphabets in allalphabets, f"Alphabets must be from {allalphabets}"

        if isinstance(_alphabets, str):
            _alphabets = [_alphabets]
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


def get_rfcw(fname, checkw=False):
    """
    Parse filename to obtain values of replicate, FOV, Channel and week
    Args:
        fname: file name
        checkw: check for week number also in name
    Returns:
         replicate, fov no, wellplate no., (week no) : r,f,c,(w)
    """
    s1, s2, r, fc = None, None, None, None
    split = fname.split("_")
    if len(split) == 3:
        s1, r, fc = split
    elif len(split) == 4:
        s1, s2, r, fc = split
    c = fc.split(".tif")[0][-3:]
    f = fc.split(".tif")[0][-16:-12]
    # print("f,c",f,c)
    if checkw:
        w = s1.split("-")[1]
        return r, f, c, w
    return r, f, c


def findtreatment(r):
    """
    Returns the type of treatement based on replicate id

    :param r: replicate id (IMP: must be converted to 0-9 range)
    :return: treatment id
    """
    assert r < 10, "error"
    treatment = None
    if r < 5:
        treatment = 0
    else:
        treatment = 1
    return treatment
