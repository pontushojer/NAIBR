import collections
import itertools
import logging
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import mpmath
import numpy as np
import scipy.optimize as optimize
import scipy.special as special
import scipy.stats as scipy_stats
from pylab import rc

from .utils import roundto

mpl.use("Agg")

logger = logging.getLogger(__name__)


class NegBin:
    def __init__(self, p=0.1, r=10):
        self.nbin = np.frompyfunc(self._nbin, 3, 1)
        self.p = p
        self.r = r

    @staticmethod
    def _nbin(k, p, r):
        return (
            mpmath.gamma(k + r)
            / (mpmath.gamma(k + 1) * mpmath.gamma(r))
            * np.power(1 - p, r)
            * np.power(p, k)
        )

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


def plot_distribution(p, distr, xlab, ylab, title, directory):
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
    plt.axis([0, max(bins), 0, max(max(y), max(n))])
    fig.savefig(os.path.join(directory, fname + ".pdf"), format="pdf")
    plt.close("all")
    return


def linked_reads(reads, chrom, configs):
    # Calculate the distance between neighbouring reads
    pair_dists = reads["start"][1:] - reads["end"][:-1]

    # Get indexes where the distance exceeds the MAX_LINKED_DIST
    breaks = np.where(pair_dists > configs.MAX_LINKED_DIST)[0] + 1

    linkedreads = []
    # Split reads using the indexes into groups of candidate linked reads
    for reads_group in np.split(reads, breaks):
        nr_reads = len(reads_group)
        start = reads_group["start"].min()
        end = reads_group["end"].max()
        mapqs = list(reads_group["mapq"])
        haps = list(reads_group["hap"])
        if nr_reads >= configs.MIN_READS and end - start >= configs.MIN_LEN:
            linkedreads.append((chrom, start, end, mapqs, haps))

    return linkedreads


def get_internal_overlap(barcode_linkedreads, barcode_overlap, configs):
    """Tally up the internal overlapp of normalized positions with linked reads
    for a barcode"""
    for linkedread in barcode_linkedreads:
        chrom, start, end, *_ = linkedread
        if end - start < configs.MAX_LINKED_DIST:
            continue

        norm_start = roundto(start, configs.MAX_LINKED_DIST)
        norm_end = roundto(end, configs.MAX_LINKED_DIST) + configs.MAX_LINKED_DIST
        positions = list(range(norm_start, norm_end, configs.MAX_LINKED_DIST))
        for s, e in itertools.combinations(positions, 2):
            barcode_overlap[(chrom, s, chrom, e)] += 1


def get_pairwise_overlap(barcode_linkedreads, barcode_overlap, configs):
    """Tally up the pair-wise overlapp of normalized positions between linked reads
    within a barcode"""
    for linkedread1, linkedread2 in itertools.combinations(barcode_linkedreads, 2):
        if linkedread1[0] > linkedread2[0] or linkedread1[1] > linkedread2[1]:
            linkedread1, linkedread2 = linkedread2, linkedread1
        chr1, start1, end1, *_ = linkedread1
        chr2, start2, end2, *_ = linkedread2
        norm_start1 = roundto(start1, configs.MAX_LINKED_DIST)
        norm_end1 = roundto(end1, configs.MAX_LINKED_DIST) + configs.MAX_LINKED_DIST
        norm_start2 = roundto(start2, configs.MAX_LINKED_DIST)
        norm_end2 = roundto(end2, configs.MAX_LINKED_DIST) + configs.MAX_LINKED_DIST

        # Connect each normalized position covered by each linked read between the
        # two.
        for id1 in range(norm_start1, norm_end1, configs.MAX_LINKED_DIST):
            for id2 in range(norm_start2, norm_end2, configs.MAX_LINKED_DIST):
                barcode_overlap[(chr1, id1, chr2, id2)] += 1


def get_linked_reads(reads_by_barcode, configs):
    linkedreads_by_barcode = collections.defaultdict(list)
    for (chrom, barcode), reads in reads_by_barcode.items():
        barcode_linkedreads = linked_reads(reads, chrom, configs)

        if barcode_linkedreads:
            linkedreads_by_barcode[barcode].extend(barcode_linkedreads)

    return linkedreads_by_barcode


def get_distributions(reads_by_barcode, configs):
    linkedreads_by_barcode = get_linked_reads(reads_by_barcode, configs)
    p_len, p_rate, barcode_overlap = get_linkedread_distributions(linkedreads_by_barcode, configs)
    return p_len, p_rate, barcode_overlap, linkedreads_by_barcode


def get_linkedread_distributions(linkedreads_by_barcode, configs):
    linkedreads = []
    barcode_overlap = collections.defaultdict(int)

    for barcode_linkedreads in linkedreads_by_barcode.values():
        linkedreads.extend(barcode_linkedreads)

        get_internal_overlap(barcode_linkedreads, barcode_overlap, configs)
        if len(barcode_linkedreads) > 1:
            get_pairwise_overlap(barcode_linkedreads, barcode_overlap, configs)

    if len(linkedreads) < 100:
        logger.warning("Too few linked reads to estimate distributions")
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
    # poisson distribution
    # p = scipy_stats.poisson(np.mean(lengths)).pmf
    # pp = lambda x: max(1e-20,float(p([x])[0]))
    # plot_distribution(pp,lengths,'Linked-read lengths (bp)','Frequency','Linked-read length distribution')
    return pp


def get_rate_distr(linkedreads):
    rate = [len(x[3]) / float(x[2] - x[1]) for x in linkedreads]
    rate.sort()
    if len(rate) > 10:
        rate = rate[int(len(rate) / 10) : int(len(rate) / 10 * 9)]
    alpha, loc, beta = scipy_stats.gamma.fit(rate)
    p = scipy_stats.gamma(alpha, loc, beta).cdf
    pp = lambda x: max(1e-20, float(p([max(x, 1e-6)])[0] - p([max(x, 1e-6) - 1e-6])[0]))
    # normal distribution
    # mu,sig = scipy_stats.norm.fit(rate)
    # p = scipy_stats.norm(mu,sig).cdf
    # pp = lambda x: max(1e-20,float(p([max(x,1e-6)])[0]-p([max(x,1e-6)-1e-6])[0]))
    # plot_distribution(pp,rate,'Per-molecule sequencing rate','Frequency','Per-molecule sequencing rate')
    return pp
