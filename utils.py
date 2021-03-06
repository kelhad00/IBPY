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
    """Keep only in paths_lst the pairs based on matching patt1 and patt2.
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