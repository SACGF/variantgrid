#!/usr/bin/env python3

"""

17th Oct, 2018
David Lawrence (david.lawrence@sa.gov.au)

- Optimisations - use numpy maths rather than Python
- Python2/3 compatible (requires "six" compatibility library)

===============================================

8th May, 2018
Paul Wang (paul.wang@sa.gov.au)

Not much work was done the last time this was looked at.

After discussion with Emily, we decided to implement both error estimation and random starting seeds together by:

1) multiple iterations
2) each iteration will use a sub-sample of the input variants (default ~ 80%, say)
3) then each coefficient will have a range of values, described by mean and stddev, this then provides the error for each coefficient
4) need to check whether there are multiple local minima. Select two coefficients with largest sigma and plot their distribution (2D scatter)

Additional arguments:
- iterations
- sampling ratio (default = 0.8)
- processes (for multiprocessing)

--siginfo is not being used. Maybe consider stripping it off.


===============================================
28th November 2017
Paul Wang (paul.wang@sa.gov.au)

This is an upgrade to add some features:
- multiple iterations using different random starting seed
- provide some sort of estimate of how good the fit is

Maybe also add:
- different output format? (e.g. HTML?)
- some simple filters for the input VCF
- multiprocessing if running lots of iterations is going to be slow

Additional arguments:
- iterations

===============================================
12th May 2016
Paul Wang (ppswang@gmail.com)

The prototype (calculate_signatures.py) is done,
so this is a production version, which does things in Classes


Calculate cancer mutation signatures

"Signatures of mutational processes in human cancer"
Alexandrov et al.
Nature 500, 415-421 (22 August 2013) doi:10.1038/nature12477
http://www.nature.com/nature/journal/v500/n7463/full/nature12477.html

Signature data downloaded from:
http://cancer.sanger.ac.uk/cosmic/signatures
http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt


Usage:

calculate_cancer_mutation_signatures.py \
  --invcf somatic_SNVs.vcf \
  --outdir output_director \
  --reference genome_reference.fa \
  --sigdata signature_data \
  --siginfo signature_information

Cancer signature data and information

1) Signature data is directly downloaded from Sanger, however, sometimes may need
to be sorted by "Substitution Type" and "Trinucleotide"

2) Signature information is copied from the Sanger page, and formatted into ConfigParser-readable format
e.g.

[Signature 1]
Cancer types = Found in all cancer types and in most cancer samples.
Proposed aetiology = Result of an endogenous mutational process initiated by spontaneous deamination of 5-methylcytosine.
Additional mutational features = Associated with small numbers of small insertions and deletions in most tissue types.
Comments = Correlates with age of cancer diagnosis.


"""
# ===================================

import multiprocessing as mp
import random
from argparse import RawTextHelpFormatter
from collections import defaultdict, namedtuple
from functools import partial

import configargparse
import numpy as np
import vcf
from scipy.optimize import minimize

from library.genomics.fasta_wrapper import FastaFileWrapper

VERSION = "1.0.1"


def handle_args():
    parser = configargparse.ArgParser(description="""Calculate cancer mutation signatures

Usage:
calculate_signatures.py \
  --invcf somatic_SNVs.vcf \
  --out output_prefix \
  --reference genome_reference.fa \
  --sigdata signature_data

""", formatter_class=RawTextHelpFormatter)
    parser.add_argument("--invcf", "-i", required=True, help="Input VCF file")
    parser.add_argument("--outprefix", "-o", required=True, help="Output prefix")

    parser.add_argument("--config", "-c", required=False, is_config_file=True, help="Configuration file")
    parser.add_argument("--reference", "-R", required=True, help="Genome reference fasta file")
    parser.add_argument("--sigdata", "-s", required=True, help="Cancer mutation signature data file")
    parser.add_argument("--siginfo", "-g", required=False,
                        help="Cancer mutation signature information. If not provided, information about"
                        "signatures will not be included in the output.")

    parser.add_argument("--iterations", "-N", required=False, default=100, type=int,
                        help="Number of repeats with different random starting seed values and subsample of input data. Default=100")
    parser.add_argument("--processes", "-P", required=False, default=1, type=int,
                        help="Number of parallel processes to use. Default=1")
    parser.add_argument("--sampling", "-S", required=False, default=0.8, type=float,
                        help="Sampling fraction of total variants. Default=0.8")

    parser.add_argument("--precision", "-p", required=False, default="float",
                        choices=["float", "float64", "float128"], help="Precision for Numpy calculations")
    parser.add_argument("--minimization", "-m", default="LS", choices=["LS", "LA"], help="Minimise by least-squares (LS) or least-absolute-values (LA)")
    return parser.parse_args()

# ===================================


TOTAL_MUTATION_TYPES = 96

comp_nucs = {"A": "T",
             "C": "G",
             "G": "C",
             "T": "A"}

MinimizationResultsTuple = namedtuple("MinimizationResultsTuple", "full bootstrapped")


def invert_muttype(m):
    # muttype is context_var, e.g. ACA_G
    inv = comp_nucs[m[2]] + comp_nucs[m[1]] + comp_nucs[m[0]] + "_" + comp_nucs[m[4]]
    return inv


class MutationSignatures:

    def __init__(self,
                 invcf=None, reference=None,
                 sigdatafile=None,
                 iterations=100, processes=1, sampling=0.8,
                 precision="float", minimization="LS"):
        """ Set data by either passing invcf and reference - to read from VCF
            or calculate yourself and set via UpdateFrequencyData / UpdateMutIndexList
        """
        self.minimization = minimization
        self.PROCESSES = processes
        self.ITERATIONS = iterations
        self.SAMPLING = sampling
        self.version = VERSION
        self.Results = None

        NUMPY_DTYPE_TOL = {
            "float": (np.float, None),
            "float64": (np.float64, 1e-15),
            "float128": (np.float128, 1e-30)}
        (self.dtype, self.tolerance) = NUMPY_DTYPE_TOL[precision]

        self.sigdata = None
        self.NUMBER_OF_SIGNATURES = 0
        self.SIGNATURE_NAMES = {}
        self.MUTATION_TYPES = []
        self.mut_index = {}
        self.ParseSignatureData(sigdatafile)

        # --------------------------
        # get observed data in VCF file
        self.ObservedRawFrequency = None   # this is now
        self.RawMutIndexList = []    # list of observed mutations (by index)
        self.TOTAL_VARS = 0

        if invcf:
            self.ParseVCF(reference, invcf)

    # ---------------------------------------------------------------------------
    def calculate_signatures(self):
        if self.TOTAL_VARS == 0:
            msg = "calculate_signatures() called with no variants. You need to set it either by passing invcf/reference to the constructor or calculating yourself and calling UpdateMutIndexList()?"
            raise ValueError(msg)

        sampling_size = max(int(self.TOTAL_VARS * self.SAMPLING), 1)

        # first set is the full data
        normalised_freq = self.ObservedRawFrequency / sum(self.ObservedRawFrequency)
        input_data_set = [(self.dtype, self.tolerance, self.minimization, normalised_freq, self.sigdata)] * self.ITERATIONS

        # these are bootstrapped data
        for _ in range(self.ITERATIONS):
            subset = random.sample(self.RawMutIndexList, sampling_size)
            sub_freqarray = mutIdxList_to_freq_array(subset, dtype=self.dtype, normalised=True)
            input_data_set.append((self.dtype, self.tolerance, self.minimization, sub_freqarray, self.sigdata))

        if self.PROCESSES > 1:
            with mp.Pool(self.PROCESSES) as pool:
                results = pool.map(process_worker, input_data_set)
                pool.join()
        else:
            results = []
            for x in input_data_set:
                r = process_worker(x)
                results.append(r)

        self.Results = MinimizationResultsTuple(full=results[:self.ITERATIONS],
                                                bootstrapped=results[self.ITERATIONS:])

    # ---------------------------------------------------------------------------
    # parse the signature data
    def ParseSignatureData(self, sig_data_file):
        with open(sig_data_file, "r", encoding="utf-8") as fin:
            # how many columns are there? how many signatures?
            header = fin.readline()
            for idx, col_name in enumerate(header.strip().split("\t")[3:]):
                self.SIGNATURE_NAMES[idx] = col_name
            self.NUMBER_OF_SIGNATURES = len(self.SIGNATURE_NAMES)

            # read data and fill self.sigdata
            self.sigdata = np.zeros((self.NUMBER_OF_SIGNATURES, 96), dtype=self.dtype)

            row_ct = 0
            for line in fin:
                bits = line.strip().split("\t")
                var = bits[0][2]
                context = bits[1]
                key = f"{context}_{var}"
                self.MUTATION_TYPES.append(key)
                self.mut_index[key] = row_ct

                for idx, val_ in enumerate(bits[3:]):
                    freq = self.dtype(val_)
                    self.sigdata[idx][row_ct] = freq
                row_ct += 1

    def ParseVCF(self, reference_fasta, vcf_filename):
        """ previous this function would update self.observed_freq, which contains normalised frequencies
            but this is no longer useful, as it won't allow random subsampling, so need to keep it as counts """

        reference = FastaFileWrapper(reference_fasta)
        invcf = vcf.Reader(filename=vcf_filename)

        raw_mut_index_list = []
        for var in invcf:
            if not var.is_snp:
                continue
            start_ = var.POS - 1
            end_ = var.POS + 1
            context = reference[var.CHROM][start_ - 1:end_]
            mutkey = f"{context}_{var.ALT[0]}"
            if mutkey not in self.mut_index:
                invkey = invert_muttype(mutkey)
                mutidx = self.mut_index[invkey]
            else:
                mutidx = self.mut_index[mutkey]
            raw_mut_index_list.append(mutidx)

        self.UpdateMutIndexList(raw_mut_index_list)

    def UpdateFrequencyData(self, newFreqData):
        self.ObservedRawFrequency = newFreqData
        self.RawMutIndexList = freqArray_to_mut_index_list(newFreqData)
        self.TOTAL_VARS = len(self.RawMutIndexList)

    def UpdateMutIndexList(self, newMutList):
        self.RawMutIndexList = sorted(newMutList)
        self.ObservedRawFrequency = mutIdxList_to_freq_array(newMutList)
        self.TOTAL_VARS = len(self.RawMutIndexList)

    def get_mutation_type_labels(self):
        mutation_type_labels = []
        for l in self.MUTATION_TYPES:
            mutation_type_labels.append(f"{l[0]} ({l[1]}>{l[4]}) {l[2]}")
        return mutation_type_labels

# ==============================================================================


def mutIdxList_to_freq_array(mutIdxList, dtype=np.int, normalised=False):
    outputarray = np.asarray([0] * TOTAL_MUTATION_TYPES, dtype=dtype)
    idx_freq = defaultdict(int)

    for midx in mutIdxList:
        idx_freq[midx] += 1
    for midx in idx_freq:
        outputarray[midx] = idx_freq[midx]
    total_ = sum(outputarray)

    # if output_type is float type, then it should be normalised
    # but we'll still use a flag to determine this
    if normalised:
        outputarray = outputarray / total_
    return outputarray


def freqArray_to_mut_index_list(freq_array):
    mut_index_list = []
    for s_idx, s_freq in enumerate(freq_array):
        if s_freq == 0:
            continue
        mut_index_list += [s_idx] * s_freq
    return mut_index_list


# input data is a normalised frequency array
def process_worker(indata):
    dtype_, tol_, minimization_, freqdata_, sigdata_ = indata
    constraints = ({'type': 'eq', 'fun': lambda x: np.sum(x) - 1})
    boundaries = tuple([(0, 1)] * len(sigdata_))

    # create a partial function based on the input sub-sampled data to optimise
    partial_fn = partial(generalised_diff_function,
                         observed_freq=freqdata_,
                         minimization=minimization_,
                         sigdata=sigdata_)
    # start with a random seed
    x0 = np.random.rand(len(sigdata_))
    x0 = x0 / np.sum(x0)
    x0 = np.array(x0, dtype=dtype_)

    minimization_data = minimize(partial_fn, x0,
                                 bounds=boundaries,
                                 constraints=constraints,
                                 tol=tol_,
                                 method="SLSQP")  # , jac=jacobian)
    fitdata = recombine(minimization_data.x, sigdata_)
    diffdata = fitdata - freqdata_

    minimization_data["fitdata"] = fitdata
    minimization_data["diffdata"] = diffdata
    minimization_data["ls_sum_diff"] = (diffdata ** 2).sum()
    minimization_data["la_sum_diff"] = np.abs(diffdata).sum()

    return minimization_data


def recombine(coeffs, sigdata):
    """calculates and returns the linear combination of coeffs * sigdata
    """
    return (coeffs * sigdata.T).sum(axis=1)


def generalised_diff_function(coeffs, observed_freq, minimization="LS", sigdata=None):
    """This is the generalised difference function to minimize
    """
    rec = recombine(coeffs, sigdata)
    diff = rec - observed_freq
    if minimization == "LS":
        diff_sum = (diff ** 2).sum()
    elif minimization == "LA":
        diff_sum = np.abs(diff).sum()
    return diff_sum


def reorder_list(original_list, new_index_order):
    new_list = []
    for x in new_index_order:
        new_list.append(original_list[x])
    return new_list

    # =================================================================
    # the Jacobian, probably not needed for now
    # def jacobian(self, coeffs):
    #     grad = [0]*n_sigs
    #     for x in range(n_sigs):
    #         for k in range(96):
    #             grad[x] = 2 * sigdata[x][k] * (coeffs[x]*sigdata[x][k]-freq_norm[k])
    #     grad = np.asarray(grad)
    #     return grad

# --------------------------------------------------------------


if __name__ == '__main__':
    args = handle_args()
    ms = MutationSignatures(invcf=args.invcf,
                            reference=args.reference,
                            sigdatafile=args.sigdata,
                            iterations=args.iterations,
                            processes=args.processes,
                            sampling=args.sampling,
                            precision=args.precision,
                            minimization=args.minimization,)

    ms.calculate_signatures()
    outprefix = args.outprefix
    #ms.plot_data(outprefix + ".signatures.png")
    #ms.plot_residue_distribution(outprefix + ".residues.png")
    #ms.plot_datafit(outprefix+".datafit.png")
