import re
import pandas as pd

def generic_convert(contig):
    contig = re.sub(r'^chr(?=[0-9XY]+$)', '', contig)
    return re.sub(r'^chrM$', 'MT', contig)


def b37_convert(contig, b37_map, build='B37'):
    try:
        r = b37_map.loc[b37_map['HG19'] == contig, build].reset_index(drop=True).at[0]
        return r
    except:
        return contig