from __future__ import print_function, division
import os
import collections
import itertools
from .utils import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pylab import rc
from scipy.stats import norm
from scipy.stats import gamma
import scipy.special as special
import scipy.optimize as optimize
import numpy as np
import mpmath

from .global_vars import *

mpl.use("Agg")


class NegBin(object):
    def __init__(self, p=0.1, r=10):
        self.nbin = np.frompyfunc(self._nbin, 3, 1)
        self.p = p
        self.r = r

    @staticmethod
    def _nbin(k, p, r):
        return mpmath.gamma(k + r) / (mpmath.gamma(k + 1) * mpmath.gamma(r)) * np.power(1 - p, r) * np.power(p, k)

    @staticmethod
    def mle(par, data, sm):
        p = par[0]
        r = par[1]
        n = len(data)
        f0 = sm / (r + sm) - p
        f1 = np.sum(special.psi(data + r)) - n * special.psi(r) + n * np.log(r / (r + sm))
        return np.array([f0, f1])

    def fit(self, data, p=None, r=None):
        if p is None or r is None:
            av = np.average(data)
            va = np.var(data)
            r = (av * av) / (va - av)
            p = (va - av) / (va)
        sm = np.sum(data) / len(data)
        x = optimize.fsolve(self.mle, np.array([p, r]), args=(data, sm))
        self.p = x[0]
        self.r = x[1]

    def pdf(self, k):
        return self.nbin(k, self.p, self.r).astype("float64")


def plot_distribution(p, distr, xlab, ylab, title):
    fname = "_".join(title.split(" "))
    nbins = 50
    fig, ax = plt.subplots()
    n, bins, patches = plt.hist(distr, nbins, density=True, facecolor="blue", alpha=0.70)
    rc("axes", linewidth=1)
    y = [p(b) for b in bins]
    plt.plot(bins, y, color="r", linewidth=5)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    x0 = max(n)
    xmax = max(bins)
    plt.axis([0, max(bins), 0, max(max(y), max(n))])
    fig.savefig(os.path.join(DIR, fname + ".pdf"), format="pdf")
    plt.close("all")
    return


def linked_reads(reads, chrom):
    reads.sort(key=lambda x: x[0])

    current_linkedread = [0, 0, 0, 0]
    linkedreads = []
    for start, end, hap, mapq in reads:
        if current_linkedread[0] == 0 or start - current_linkedread[2] > MAX_LINKED_DIST:
            if current_linkedread[0] != 0 and current_linkedread[3] >= MIN_READS and current_linkedread[2] - current_linkedread[1] >= MIN_LEN:
                linkedreads.append(current_linkedread)
            current_linkedread = [chrom, start, end, 1]
        else:
            current_linkedread[2] = max(current_linkedread[2], end)
            current_linkedread[3] += 1
    if current_linkedread[0] != 0 and current_linkedread[3] >= MIN_READS and current_linkedread[2] - current_linkedread[1] >= MIN_LEN:
        linkedreads.append(current_linkedread)
    return linkedreads


def get_overlap(barcode_linkedreads, barcode_overlap):
    for linkedread1, linkedread2 in itertools.combinations(barcode_linkedreads, 2):
        if linkedread1[0] > linkedread2[0] or linkedread1[1] > linkedread2[1]:
            linkedread1, linkedread2 = linkedread2, linkedread1
        chr1, start1, end1, count1 = linkedread1
        chr2, start2, end2, count2 = linkedread2
        index1 = {roundto(start1, MAX_LINKED_DIST), roundto(end1, MAX_LINKED_DIST)}
        index2 = {roundto(start2, MAX_LINKED_DIST), roundto(end2, MAX_LINKED_DIST)}
        for id1 in index1:
            for id2 in index2:
                barcode_overlap[(chr1, id1, chr2, id2)] += 1


def get_distributions(reads_by_barcode):
    linkedreads = []
    linkedreads_by_barcode = collections.defaultdict(list)
    barcode_overlap = collections.defaultdict(int)

    for key, reads in reads_by_barcode.items():
        chrom, barcode = key
        barcode_linkedreads = linked_reads(reads, chrom)
        linkedreads += barcode_linkedreads
        linkedreads_by_barcode[barcode] += barcode_linkedreads

    for barcode, barcode_linkedreads in linkedreads_by_barcode.items():
        if len(barcode_linkedreads) > 1:
            get_overlap(barcode_linkedreads, barcode_overlap)

    if len(linkedreads) < 100:
        return None, None, None

    p_rate = get_rate_distr(linkedreads)
    p_len = get_length_distr(linkedreads)
    return p_len, p_rate, barcode_overlap


def get_length_distr(linkedreads):
    lengths = [x[2] - x[1] for x in linkedreads]
    lengths.sort()
    assert len(lengths) >= 100
    assert np.var(lengths) != 0
    b = NegBin()
    b.fit(lengths)
    p = b.pdf
    pp = lambda x: max(1e-20, float(p([x])[0]))
    ## poisson distribution
    # p = scipy.stats.poisson(np.mean(lengths)).pmf
    # pp = lambda x: max(1e-20,float(p([x])[0]))
    # plot_distribution(pp,lengths,'Linked-read lengths (bp)','Frequency','Linked-read length distribution')
    return pp


def get_rate_distr(linkedreads):
    rate = [x[3] / float(x[2] - x[1]) for x in linkedreads]
    rate.sort()
    if len(rate) > 10:
        rate = rate[int(len(rate) / 10) : int(len(rate) / 10 * 9)]
    alpha, loc, beta = scipy.stats.gamma.fit(rate)
    p = scipy.stats.gamma(alpha, loc, beta).cdf
    pp = lambda x: max(1e-20, float(p([max(x, 1e-6)])[0] - p([max(x, 1e-6) - 1e-6])[0]))
    ## normal distribution
    # mu,sig = scipy.stats.norm.fit(rate)
    # p = scipy.stats.norm(mu,sig).cdf
    # pp = lambda x: max(1e-20,float(p([max(x,1e-6)])[0]-p([max(x,1e-6)-1e-6])[0]))
    # plot_distribution(pp,rate,'Per-molecule sequencing rate','Frequency','Per-molecule sequencing rate')
    return pp
