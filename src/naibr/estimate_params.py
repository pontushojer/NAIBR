from __future__ import print_function, division
import os, sys, subprocess, pysam, collections, time, gc
import multiprocessing as mp
from .utils import *
import numpy as np
from .global_vars import *


def estimate_delta(read_list):
    counter = 0
    gaps = []
    for b, l in read_list:
        counter += 1
        if counter > 10000:
            break
        l.sort(key=lambda x: x[0])
        qname = ""
        prev_end = 0
        for (chrom, start, end, hap, qname2, mapq) in l:
            if qname2 != qname and prev_end > 0 and start - prev_end < max_d:
                gaps.append(start - prev_end)
            qname = qname2
            prev_end = end
    mean = np.mean(gaps)
    gaps.sort()
    d = max(1000, int(gaps[len(gaps) - int((len(gaps) / 100) * 5)] / 1000) * 1000)
    return d
