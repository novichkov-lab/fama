"""
Distribution of Read Insert Sizes Model

Reference:
Kristoffer Sahlin, Mattias Franberg, Lars Arvestad
Correcting bias from stochastic insert size in read pair data — applications to
structural variation detection and genome assembly
http://biorxiv.org/content/biorxiv/early/2015/08/04/023929.full.pdf

Original version:
https://github.com/ksahlin/GetDistr/blob/master/getdistr/assemblymodule/lib_est.py
"""
import os
import math
import bisect
import random
import json

# import numpy
from statsmodels.distributions.empirical_distribution import ECDF
# from rpy2 import robjects
# from rpy2.robjects import packages

from mathstats.normaldist.normal import MaxObsDistr
from mathstats.normaldist.truncatedskewed import param_est

try:
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_palette("husl", desat=.6)
    sns.set(font_scale=1.6)
    plt.rcParams.update({'font.size': 22})
except ImportError:
    pass

EMPIRICAL_BINS = 500
SAMPLE_SIZE = 2**32  # for estimating true full read pair distribution
NORMAL_QUANTILE_TWO_SIDED_95 = 1.96


class Parameters(object):
    """docstring for Parameters"""
    def __init__(self):
        super(Parameters, self).__init__()
        self.mu = None
        self.sigma = None
        self.adjusted_mu = None
        self.adjusted_sigma = None
        self.min_isize = None
        self.max_isize = None
        self.read_length = None
        self.pval = None
        self.total_basepairs = None
        self.ess_ratio = None

        self.scaf_lengths = {}
        self.outfolder = None
        self.plots = None
        self.stats_output = None

        self.nr_reads = 0
        self.nr_mapped = 0
        self.nr_proper_mapped = 0


def plot_isize(isizes, outfile):
    plt.hist(isizes, bins=100)
    plt.ylabel('frequency')
    plt.xlabel('fragment size')
    plt.title('Insert size distribution')
    plt.legend()
    plt.savefig(outfile)
    plt.close()
    plt.clf()


def adjust_insert_size_dist(mean_insert, std_dev_insert, insert_list):
    k = MaxObsDistr(len(insert_list), 0.95)
    # for ins in insert_list:
    #    if (ins > mean_insert + k * std_dev_insert or ins < mean_insert - k * std_dev_insert):
    #        print ins
    filtered_list = list(filter(
        (lambda x: (
            x < mean_insert + k * std_dev_insert and x > mean_insert - k * std_dev_insert
            )),
        insert_list
        ))
    result = False
    if len(insert_list) > len(filtered_list):
        result = True
    return(result, filtered_list)


def is_proper_aligned_unique_innie(read):
    mapped_ok = not (read.is_unmapped or read.mate_is_unmapped)
    orientation_ok = (
        read.is_reverse and not read.mate_is_reverse and read.tlen < 0
    ) or (
        not read.is_reverse and read.mate_is_reverse and read.tlen > 0
    )
    quality_ok = read.mapq > 10 and not read.is_secondary
    same_ref = read.rname == read.mrnm

    return mapped_ok and orientation_ok and quality_ok and same_ref

    # return not read.is_unmapped and not read.mate_is_unmapped and \
    #            ( (read.is_reverse and not read.mate_is_reverse  and read.tlen < 0 ) or \
    #              (not read.is_reverse and read.mate_is_reverse  and read.tlen > 0 ) ) \
    #              and  and read.mapq > 10 and not read.is_secondary and read.rname == read.mrnm


class LibrarySampler(object):
    """docstring for LibrarySampler"""
    def __init__(self, read_data, param):
        super(LibrarySampler, self).__init__()
        self.read_data = read_data
        self.param = param
        self.mu_sophisticated = 0.0
        self.lib_file = os.path.join(param.outfolder, 'library_info.txt')
        self.stats_file = os.path.join(param.outfolder, 'library_stats.txt')

        self.true_distr = None
        self.adjustedECDF_no_gap = None
        self.adjusted_mean = None
        self.adjusted_stddev = None

        self.sample_distribution()

    def get_correct_ecdf(self):
        read_len = int(self.read_length)
        softclipps = 0  # read_len #int(self.read_length*0.6)

        x_min = self.min_isize  # max(2*(read_len-softclipps) , int(self.mean - 5*self.stddev) )
        x_max = self.max_isize  # int(self.mean + 5*self.stddev)
        stepsize = int(max(1, (x_max - x_min) / EMPIRICAL_BINS))
        # print (x_max - x_min)/float(EMPIRICAL_BINS)
        cdf_list = [self.full_ECDF(x_min) * self.get_weight(
            int(round(x_min+stepsize/2.0, 0)),
            read_len,
            softclipps
            )]  # [ self.full_ECDF( 2*(read_len-softclipps)) * self.get_weight(2*\
        # (read_len-softclipps), gap_coordinates, read_len, softclipps) ]

        # create weigted (true) distribution
        # print('x_min:',x_min,',x_max:',x_max,',stepsize:',stepsize)
        for x_val in range(x_min + stepsize, x_max, stepsize):
            increment_area = self.get_weight(
                int((x_val + stepsize / 2)), read_len, softclipps
                ) * (self.full_ECDF(x_val) - self.full_ECDF(x_val - stepsize))
            cdf_list.append(cdf_list[-1] + increment_area)

        tot_cdf = cdf_list[-1]

        cdf_list_normalized = list(map(lambda x: x / float(tot_cdf), cdf_list))

        # Now create a weighted sample

        self.true_distr = [
            bisect.bisect(
                cdf_list_normalized, random.uniform(0, 1)
                )
            * stepsize + x_min for i in range(50000)
            ]

        # initialization of no gap true distribution
        self.adjustedECDF_no_gap = self.true_distr

        n = len(self.true_distr)
        self.adjusted_mean = sum(self.true_distr)/float(len(self.true_distr))
        self.adjusted_stddev = (sum(list(map((
            lambda x: x ** 2 - 2 * x * self.adjusted_mean + self.adjusted_mean ** 2
            ), self.true_distr))) / (n - 1)) ** 0.5

    def sample_distribution(self):
        # pilot sample
        read_lengths = []
        # max_tlen = 0
        # bam_filtered = ifilter(lambda r: is_proper_aligned_unique_innie(r), self.bamfile)
        isize_list = []
        mcmc_dict = {}
        reference_list = []
        reference_name_list = []
        # nr_reads = 0
        # nr_mapped = 0
        # nr_proper_mapped = 0
        # for sample_nr,read in enumerate(self.bamfile):
        reads_processed = 0
        for read_id in self.read_data:
            # add do insert size distribution calculation if proper pair
            # if is_proper_aligned_unique_innie(read) and not  read.is_reverse:
            self.param.nr_proper_mapped += 2  # add the read plus its mate since the mate
            # does not enter here
            # assert read.tlen > 0
            assert self.read_data[read_id]['tlen'] > 0
            read_lengths.append(self.read_data[read_id]['rlen'])
            isize_list.append(self.read_data[read_id]['tlen'])
#            if read.tid in mcmc_dict:
#                mcmc_dict[str(reads_processed)].append(self.read_data[read_id]['tlen'])
#            else:
            mcmc_dict[str(reads_processed)] = [self.read_data[read_id]['tlen']]
            reference_list.append(self.read_data[read_id]['ref_len'])
            reference_name_list.append(self.read_data[read_id]['ref_name'])
            # if abs(read.tlen) > max_tlen:
            #    max_tlen = abs(read.tlen)
            reads_processed += 1
            if reads_processed >= SAMPLE_SIZE:
                break
        # print(reads_processed, 'reads processed')
        # for sample_nr,read in enumerate(bam_filtered):
        #         add do insert size distribution calculation if proper pair
        #    if is_proper_aligned_unique_innie(read) and not  read.is_reverse:
        #        assert read.tlen > 0
        #        read_lengths.append(read.rlen)
        #        isize_list.append(read.tlen)
        #        # if abs(read.tlen) > max_tlen:
        #        #    max_tlen = abs(read.tlen)
        #    if sample_nr >= SAMPLE_SIZE:
        #        break

        # for read, mate_pos in fb.proper_read_isize(
        #      self.bamfile, self.param.lib_min, self.param.ligetdistr.assemblymodule.b_max
        #    ):
        #         sample_nr += 1
        #         ## add do insert size distribution calculation if proper pair
        #        if read.tlen >= 0:
        #    #if is_proper_aligned_unique_innie(read) and read.is_read1:
        #        read_lengths.append(read.rlen)
        #        isize_list.append(read.tlen)
        #        # if abs(read.tlen) > max_tlen:
        #        #    max_tlen = abs(read.tlen)
        #    if sample_nr >= SAMPLE_SIZE:
        #        break

#        self.bamfile.reset()
        # max_tlen = max_tlen+1000
        self.read_length = sum(read_lengths)/len(read_lengths)
        # print('Read length', str(self.read_length))

        # sample proper reads

        # isize_list = []
        # for sample_nr,read in enumerate(proper_read_isize_iter(
        #    self.bampath, self.read_length, max_tlen
        #    )):
        #        isize_list.append(read)
        #    if sample_nr > SAMPLE_SIZE:
        #        break
        params = dict()
        params["sample-nr"] = reads_processed  # sample_nr

        # print('isize_list before filter', str(len(isize_list)))
        # isize_list = filter(lambda x: 0 < x - 2*self.read_length,isize_list)
        # isize_list = list(filter(lambda x: 0 < x - 2*self.read_length,isize_list))
        # print('after:', str(len(isize_list)))

        n_isize = len(isize_list)
        # print ('n_isize', n_isize)
        mean_isize = sum(list(isize_list))/n_isize
        std_dev_isize = (sum(list(map((
            lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2
            ), isize_list))) / (n_isize - 1)) ** 0.5
        params["mu-raw"] = mean_isize
        params["sd-raw"] = std_dev_isize
        extreme_obs_occur = True
        while extreme_obs_occur:
            # print 'HERE!!'
            # print(mean_isize, std_dev_isize, str(len(isize_list)))
            extreme_obs_occur, filtered_list = adjust_insert_size_dist(
                mean_isize, std_dev_isize, isize_list
                )
            n_isize = float(len(filtered_list))
            mean_isize = sum(filtered_list) / n_isize
            std_dev_isize = (sum(list(map((
                lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2
                ), filtered_list))) / (n_isize - 1)) ** 0.5
            isize_list = filtered_list

        self.min_isize, self.max_isize = min(isize_list), max(isize_list)

        # filter outliers
        for ref in mcmc_dict:
            ref_isizes = mcmc_dict[ref]
            mcmc_dict[ref] = list(filter(
                lambda x: self.min_isize <= x <= self.max_isize, ref_isizes
                ))

        params["mu-filtered"] = mean_isize
        params["sd-filtered"] = std_dev_isize
        params["min-isize"] = self.min_isize
        params["max-isize"] = self.max_isize
        params["read-length"] = self.read_length

        self.nobs = n_isize
        self.mean = mean_isize
        self.stddev = std_dev_isize
        self.full_ECDF = ECDF(isize_list)
        self.adjustedECDF_no_gap = None
        self.adjustedECDF_no_gap = self.get_correct_ecdf()
        params["mu-adjusted"] = self.adjusted_mean
        params["sd-adjusted"] = self.adjusted_stddev

        samples = min(SAMPLE_SIZE, len(isize_list))
        # ess = self.effectiveSampleSize(mcmc_dict) #isize_list[:samples]) # mcmc_dict ) #
        # self.ess_ratio = ess / float(sum(map(lambda x: len(mcmc_dict[x]), mcmc_dict)))
        params["ess"] = 1  # self.ess_ratio
        reference_lengths = map(lambda x: int(x), reference_list)
        # reference_lengths = map(lambda x: int(x), self.bamfile.lengths)
        # ref_list = zip(self.bamfile.references, reference_lengths)
#        ref_list = zip(reference_name_list, reference_lengths)

        total_basepairs = sum(reference_lengths)
        self.param.total_basepairs = total_basepairs
        params["genome-length"] = total_basepairs
#        params["contigs"] = []
#        for ref, length in list(ref_list):
#            params["contigs"].append( { "name" : ref, "length" : length } )
        with open(self.lib_file, 'w') as out_f:
            json.dump(params, out_f, sort_keys=True, indent=4, separators=(',', ': '))
            out_f.close()

        params = dict()
        if self.param.nr_reads:
            reads = dict()
            reads["total"] = self.param.nr_reads
            reads["mapped"] = self.param.nr_mapped
            reads["properly-mapped"] = self.param.nr_proper_mapped
            reads["mapped-percentage"] = self.param.nr_mapped/float(self.param.nr_reads)
            reads["properly-mapped-percentage"] = self.param.nr_proper_mapped/float(
                self.param.nr_reads
                )
            reads["coverage"] = self.param.nr_reads/float(total_basepairs)
            reads["coverage-mapped"] = self.param.nr_mapped/float(total_basepairs)
            reads["coverage-properly-mapped"] = self.param.nr_proper_mapped/float(total_basepairs)
            params["reads"] = reads

        info = dict()
        info["proper-samples"] = samples
        info["ess-proper-samples"] = 1  # ess
        info["ess-ratio"] = 1  # self.ess_ratio
        coverage = self.read_length*samples*2/float(total_basepairs)
        info["mean-coverage-proper"] = coverage
        inner_span_coverage = coverage * (
            self.mean - 2 * self.read_length
            ) / (2 * self.read_length) if self.mean - 2 * self.read_length > 0 else coverage
        info["average-theoretical-inner-span-coverage"] = inner_span_coverage
        info["mu-full-lib"] = self.mean
        info["sd-full-lib"] = self.stddev
        info["mu-empirical"] = self.adjusted_mean
        info["sd-empirical"] = self.adjusted_stddev
        mu_naive = self.mean + self.stddev**2/float(self.mean - 2*self.read_length+1)
        sigma_naive = math.sqrt(
            abs(self.stddev ** 2 - self.stddev ** 4 / (self.mean - 2 * self.read_length + 1) ** 2)
            )
        info["mu-naive"] = mu_naive
        info["sd-naive"] = sigma_naive
        mu_sophisticated = param_est.mean_given_d(
            self.mean, self.stddev, self.read_length, total_basepairs, total_basepairs, 0
            )
        print('Average insert size:', str(mu_sophisticated))
        sigma_sophisticated = param_est.stddev_given_d(
            self.mean, self.stddev, self.read_length, total_basepairs, total_basepairs, 0
            )
        info["mu-sophisticated"] = mu_sophisticated
        self.mu_sophisticated = mu_sophisticated
        info["sd-sophisticated"] = sigma_sophisticated
        theoretical_margin_of_error = NORMAL_QUANTILE_TWO_SIDED_95*self.stddev / math.sqrt(
            inner_span_coverage
            )
        info["theoretical-error-margin-two-sided-95"] = theoretical_margin_of_error
        params["extra-info"] = info
        with open(self.stats_file, 'w') as out_f:
            json.dump(params, out_f, sort_keys=True, indent=4, separators=(',', ': '))
            out_f.close()

        # if self.param.plots:
            # outfile = os.path.join(self.param.plotfolder, 'isize.eps')
            # plot_isize(isize_list, outfile)
            # outfile = os.path.join(self.param.plotfolder, 'fitted_params_isize.eps')
            # fit.main(isize_list, outfile)

    def get_weight(self, val, read_len, softclipps):
        return val - (2 * (read_len - softclipps) - 1)

    def lowest_bound_biological(self):
        pass

    def plot(self):
        pass

    # def effectiveSampleSize_by_sampling(self, data):
        # """
        # Effective sample size by comparing emperical sample error to,
        # theoretical sample error
        # """

    # def effectiveSampleSize(self, data):
        # """
        # Effective sample size, as computed by EffectiveSize in coda, R.
        # returns a python float, the ess
        # """
        # robjects.packages.importr("coda")
        # r = robjects.r

        # ess_dict = {}
        # for ref in data.keys():
        #     if len(data[ref]) < 100:
        #         continue
        #     array_data = numpy.array(data[ref])
        #    normalizedData = array_data - array_data.mean()
        #     print array_data.mean()
        #     print array_data
        #    rData = robjects.IntVector(normalizedData)
        #    mcmc_r = r.mcmc(rData)
        #     print  mcmc_r
        #     print r.effectiveSize(mcmc_r)
        #     print list(r.effectiveSize(mcmc_r))
        #    ess_dict[ref] = list(r.effectiveSize(mcmc_r))[0]

        # print ess_dict
        # print sum(map(lambda x: ess_dict[x], ess_dict))
        # print sum(map(lambda x: len(data[x]), data))
        # # sum_ess = 0
        # # for ref in data:
        # #        nr_samples = len(data[ref])
        # #        sum_ess += nr_samples * ess_dict[ref] sum(ess_vector)
        # ess = sum(map(lambda x: ess_dict[x], ess_dict))

        # tot_samples = data #[item for sublist in data.values() for item in sublist]
        # array_data_tot = numpy.array( tot_samples )
        # normalizedData = array_data_tot - array_data_tot.mean()
        # print array_data_tot.mean()
        # rData = robjects.IntVector(normalizedData)
        # mcmc_r = r.mcmc(rData)
        # ess_tot = list(r.effectiveSize(mcmc_r))[0]
        # print 'Global:', ess_tot

        #   return ess
        # number_of_samples = len(data)
        # max_batches = number_of_samples / 10000
        # ess_list = []
        # for nr_batches in [1, max_batches]: #range(1, max_batches+1):
        #    chunk_size = number_of_samples / nr_batches
        #    print chunk_size
        #    ess_list.append([])
        #    for chunk in chunks(data,chunk_size):
        #        # dont look at the small rest chunk
        #        if len(chunk) == chunk_size:
        #            array_data = numpy.array(chunk)
        #            normalizedData = array_data - array_data.mean()
        #            rData = robjects.IntVector(normalizedData)
        #            mcmc_r = r.mcmc(rData)
        #            ess_list[-1].append( list(r.effectiveSize(mcmc_r))[0] )

        # print ess_list
        # sum_ess = 0

        # for i, ess_vector in enumerate(ess_list):
        #    sum_ess +=    sum(ess_vector)
        # avg_ess = sum_ess / len(ess_list)
        # return avg_ess

        # array_data = numpy.array(data)
        # normalizedData = array_data - array_data.mean()
        # rData = robjects.IntVector(normalizedData)
        # mcmc_r = r.mcmc(rData)
        # return list(r.effectiveSize(mcmc_r))[0]

# def chunks(l, n):
#    """ Yield successive n-sized chunks from l.
#    """
#    for i in range(0, len(l), n):
#        yield l[i:i+n]


def get_lib_est(read_data, outfolder):
    param = Parameters()
    param.plots = False
    param.outfolder = outfolder
    param.stats_output = os.path.join(param.outfolder, 'stats_output.txt')
    # param.lib_min = args.lib_min
    # param.lib_max = args.lib_max

    if not os.path.exists(param.outfolder):
        os.makedirs(param.outfolder)
    return LibrarySampler(read_data, param).mu_sophisticated
