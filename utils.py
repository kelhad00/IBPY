import os
import ffmpeg  # pip ffmpeg-python

##check duration
def check_duration(eaf):
    max_dur = 0
    names = eaf.get_tier_names()
    for n in names:
        tier = eaf.get_annotation_data_for_tier(n)
        if len(tier) < 1:
            continue
        if tier[-1][1] > max_dur:
            max_dur = tier[-1][1]
    return max_dur


def seconds_to_hhmmss(seconds):
    """convert seconds format to HH:MM:SS for vizualization"""
    hours = seconds // (60 * 60)
    seconds %= 60 * 60
    minutes = seconds // 60
    seconds %= 60
    return "%02i:%02i:%02i" % (hours, minutes, seconds)


def ind_to_val(lst_ind, lst_vals):
    """transform the list of indices to list of corresponding values"""
    return [lst_vals[i] for i in lst_ind]


def keep_pairs(paths_lst, patt1, patt2):
    """keep only in paths_lst the pairs based on matching patt1 and patt2.
    """
    ref = set(paths_lst)
    pairs = []
    for l in paths_lst:
        if os.path.basename(l).split("_")[1] == patt1:
            temp = l.replace(patt1, patt2)
            if temp in ref:
                pairs.extend([l, temp])
    return pairs


class AttributeGenerator:
    """create attributes from strings.
    """
    def _add(self, name, val):
        setattr(self, "_" + name, val)
        setattr(
            self,
            name,
            property(self._getter("_" + name), self._setter("_" + name, val)),
        )
        setattr(self, name, val)

    def _getter(self, attr):
        # add condition if needed
        return getattr(self, attr)

    def _setter(self, attr, val):
        # add condition if needed
        return setattr(self, attr, val)

    @property
    def _list_attr(self):
        """Get list of attributes added."""
        lst = [l[1:] for l in self.__dict__ if l.startswith("_")]
        return lst

def split_segments(file_path, dest_path, label, lst_vals, to=False):
    """split one file into several based on lst_vals.

    Args:
        file_path (str): path to file to split.
        dest_path (str): path to directory where to save.
        label (str): label corresponding to lst_vals.
        lst_vals (list): (start, stop, val)
    """

    for strt, stp, val in lst_vals:
        if not to:
            stp = stp - strt
        strt = str(strt / 1000)  # from ms to sec
        stp = str(stp / 1000)
        out_name, ext = os.path.splitext(os.path.basename(file_path))
        #TODO: add more flexibility to naming
        out_name = '_'.join([out_name, label, val, strt.replace('.', ''), stp.replace('.','')])
        out_name += ext
        out_name = os.path.join(dest_path, out_name)
        ffmpeg.input(file_path, **{"ss": strt, "t": stp}).output(out_name).run()


def ffmpeg_convert(infile, outfile):
    ffmpeg.input(infile).output(outfile).run()

def overlapping_dct_from_indices_to_vals(dct_inds, lstA, lstB):
    """Convert dictionary of indices to dictionary of values.
    
    Args:
        dct_inds (dict): {indA:[indB]} (output of get_overlapping_segments_ind output).
        lstA (list): labels of which the indices in dct_inds correspond to.
        lstB (list): labels of which the indices in dct_inds correspond to.
    
    Returns:
        dict: {val: [vals]}
    """

    dct_vals = {}
    for indA, B in dct_inds.items():
        dct_vals[lstA[indA]] = [lstB[indB] for indB in B]
    return dct_vals



class TierNavigator:
    def __init__(self, tier, ind = 0):
        self.tier = tier
        self.ind = ind
        self.len = len(tier)
    
    def inc(self, val = 1):
        if not (0 <= (self.ind + val) < len(self.tier)):
            raise ValueError("Index would become out of range.")
        self.ind+=val

    @property
    def next_ind(self):
        if (self.ind+1) >= len(self.tier):
            print("WARNING:No further element in tier.")
            return None
        return self.ind+1

    @property
    def stt(self):
        return self.tier[self.ind][0]

    @property
    def stp(self):
        return self.tier[self.ind][1]

    @property
    def next_stt(self):
        if (self.ind+1) > len(self.tier):
            print("WARNING:Index out of range.")
            return None
        return self.tier[self.ind+1][0]

    @property
    def next_stp(self):
        if (self.ind+1) > len(self.tier):
            print("WARNING:Index out of range.")
            return None
        return self.tier[self.ind+1][1]
    
    @property
    def prev_stt(self):
        if (self.ind-1) <= 0:
            print("WARNING:Index out of range.")
            return None
        return self.tier[self.ind-1][0]

    @property
    def prev_stp(self):
        if (self.ind-1) <= 0:
            print("WARNING:Index out of range.")
            return None
        return self.tier[self.ind-1][1]


def slice_and_split(seg, width, max_dur, overlap_perc = 0.33):
    """[summary]

    Args:
        seg ([type]): [description]
        width ([type]): [description]
        max_dur ([type]): [description]
        overlap_perc (float, optional): windows overlap ratio. Defaults to 0.3.

    Returns:
        [type]: [description]
    """
    stt, stp, val = seg
    seg_width = stp-stt
    final = []
    step = round(width*(1-overlap_perc))
    for i in range(stt, stp, step):
        if i+width>max_dur:
            return final
        print('start:{}, stop:{}, step:{}, val:{}'.format(stt, stp, step, val))
        final.append((i, i+width, val))
    return final

def pad_and_split(seg, width, max_dur, overlap_perc = 0.33):
    dur = seg[1]-seg[0]
    final = []
    step = round(width*(1-overlap_perc))
    stt = seg[1]-width
    stp = seg[0]+width
    for t in range(stt, stp, step):
        if (t < 0) or ((t+dur) > max_dur):# window before 0 or after total duration, skip it.
            continue
        if t >= seg[1]: # if window start>than segment length
            break
        final.append((t, t+width, seg[2]))
    return final

def adjust_segments_len(seg_lst, width, max_dur, overlap_perc=0.33):
    """Generate segments of fix lengths.

    Args:
        seg_lst (list of tuples): [(stt,stp,lab)]
        width (numeric): desired length in milliseconds.
        max_dur (numeric): total duration considered for entire tier in milliseconds.
        overlap_perc (float): percentage of window overlap.

    Returns:
        [type]: [description]
    """
    final_segs = []
    for seg in seg_lst:
        if seg[1]-seg[0]>width:
            final_segs.extend(slice_and_split(seg, width, max_dur, overlap_perc))
        elif seg[1]-seg[0]<width:
            final_segs.extend(pad_and_split(seg, width, max_dur, overlap_perc))
        else:
            final_segs.append(seg)
    return final_segs

def intersections():
    #same as union but modify which ones to keep.
    pass

def create_none_tier(tier, max_dur, stt_time = 0, lab = ''):
    """return segments not present in tier.

    Args:
        tier (list of tuples): [(stt,stp,lab)]
        max_dur (numeric): total duration considered for entire tier in milliseconds.
        stt_time (int, optional): tier starting time considered. Defaults to 0.
        lab (str, optional): label to add to the union segments tuples. Defaults to ''.

    Returns:
        list of tuples: [(stt,stp,lab)] tier of segments not present in input tier.
    """
    tier = tier[:] # avoid alteraction of input due to referencing
    tier.insert(0, (stt_time, stt_time, lab))
    tier.append((max_dur, max_dur, lab))
    tier = TierNavigator(tier)
    final = []
    while tier.next_ind:
        if tier.next_stt - tier.stp > 0:
            final.append((tier.stp, tier.next_stt, lab))
        tier.inc()
    return final

def union_tiers(tier1, tier2, lab = '', margin = 0):#margin of error 50 ms: in case of annotation error
    """return the union of the tiers' segments.

    To test :
    test1 = [(1, 10,'v'), (12, 20, 'v'), (22, 25, 'v')]
    test2 = [(3, 5, 'd'), (8, 13, 'c'), (17, 20, 'c'), (26, 27, 'f')]

    test1 = [(1,3,'a'), (7,8,'a')]
    test2 = [(4,5,'v')]

    Args:
        tier1 (list of tuples): [(stt,stp,lab)]
        tier2 (list of tuples): [(stt,stp,lab)]
        lab (str, optional): label to add to the union segments tuples. Defaults to ''.
        margin (int, optional): margin of annotation error in ms. Defaults to 50.

    Returns:
        list of tuples: [(stt,stp,lab)] union of input tiers segments.
    """
    tiers = tier1[:]
    tiers.extend(tier2)
    tiers = sorted(tiers, key = lambda x: (x[0], x[1]))
    tiers = TierNavigator(tiers)
    final = []
    while tiers.next_ind:
        stt = tiers.stt
        stp = tiers.stp
        while stp >= (tiers.stt+margin):
            if tiers.stp > (stp+margin):
                stp = tiers.stp
            else:
                if tiers.next_ind:
                    tiers.inc()
                else:
                    break
        final.append((stt, stp, lab))
    for i in range(tiers.ind, tiers.len-1):
        stt, stp, _ = tiers.tier[i]
        final.append((stt, stp, lab))
    return final
