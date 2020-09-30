import argparse
import glob
from itertools import chain
import drug_learning.two_dimensions.Input.fingerprints as fp
import drug_learning.two_dimensions.Input.applicability_domain as ad
import drug_learning.two_dimensions.Helpers.split as sp
import drug_learning.two_dimensions.Helpers.parallel as pl
import drug_learning.two_dimensions.Errors.errors as er


def get_parser():
    parser = argparse.ArgumentParser(description="Specify the fingerprints and the file format to which the user wants the sdf file to be converted.")
    parser.add_argument("infile",
                        type=str,
                        help="Input sdf file(s)")

    parser.add_argument('-n', "--nworkers",
                        dest = "n_workers",
                        action ='store',
                        default = 1,
                        type = int,
                        help = "Number of workers to parallelize the sdf transform into fingerprints")

    parser.add_argument('-mo', '--morgan',
                        dest ="morgan",
                        action = "store_const",
                        const = fp.MorganFP,
                        default = None,
                        help = "Convert molecules to Morgan fingerprint")

    parser.add_argument('-ma', "--maccs",
                        dest = "maccs",
                        action = "store_const",
                        const = fp.MACCS_FP,
                        default = None,
                        help = "Convert molecules to MACCS fingerprint")

    parser.add_argument('-rd', "--rdkit",
                        dest = "rdkit",
                        action = "store_const",
                        const = fp.RDkitFP,
                        default = None,
                        help = "Convert molecules to RDkit fingerprint")

    parser.add_argument('-md', "--mordred",
                        dest = "mordred",
                        action ='store_const',
                        const = fp.MordredFP,
                        default = None,
                        help = "Convert molecules to Mordred fingerprint")

    parser.add_argument('-urd', "--unfolded_rdkit",
                        dest = "unfolded_rdkit",
                        action ='store_const',
                        const = fp.UnfoldedRDkitFP,
                        default = None,
                        help = "Convert molecules to Unfolded RDkit fingerprint")

    parser.add_argument('-voc', "--vocabulary",
                        dest = "voc",
                        type = str,
                        default = None,
                        help = "Vocabulary for unfolded rdkit fingerprint")

    parser.add_argument('-c', "--csv",
                        dest = "to_csv",
                        action = "store_true",
                        help = "Save output to csv")

    parser.add_argument('-pq', "--parquet",
                        dest = "to_parquet",
                        action = "store_true",
                        help = "Save output to parquet")

    parser.add_argument('-f', "--feather",
                        dest = "to_feather",
                        action ="store_true",
                        help = "Save output to feather")

    parser.add_argument('-hdf', "--hdf",
                        dest = "to_hdf",
                        action = "store_true",
                        help = "Save output to hdf")

    parser.add_argument('-pk', "--pickle",
                        dest = "to_pickle",
                        action ="store_true",
                        help = "Save output to pickle")


    subparsers = parser.add_subparsers(dest = "split")

    parser_split = subparsers.add_parser("split", help = "Split the input sdf files into chunks" )

    parser_split.add_argument('-ch', "--chunk",
                        dest = "n_chunks",
                        action ='store',
                        default = 1000,
                        type = int,
                        help = "Number of molecules in each chunk.")

    parser_split.add_argument('-nsp', "--nworkers-sp",
                        dest = "nw_sp",
                        action ='store',
                        default = 1,
                        type = int,
                        help = "Number of workers to parallelize the split into chunks.")

    return parser

def sdf_to_fingerprint(input_file, fp_list, urdkit_voc, format_dict):
    for fp_class in fp_list:
       if fp_class:
            if fp_class == fp.UnfoldedRDkitFP:
                fingerprint = fp_class(urdkit_voc)
            else:
                fingerprint = fp_class()
            fingerprint.fit(input_file)
            fingerprint.transform()
            fingerprint.save(**format_dict)

def main():
    parser = get_parser()
    opt = parser.parse_args()

    input_sdfs = glob.glob(opt.infile)

    if opt.split:
        input_sdfs = pl.parallelize(sp.split_in_chunks, input_sdfs, opt.nw_sp, n_chunks = opt.n_chunks)
        input_sdfs = list(chain.from_iterable(input_sdfs))

    fp_options = [opt.morgan, opt.maccs, opt.rdkit, opt.mordred, opt.unfolded_rdkit]

    format_options = {'to_csv' : opt.to_csv, 'to_parquet' : opt.to_parquet, 'to_feather' : opt.to_feather,
                     'to_hdf' : opt.to_hdf, 'to_pickle' : opt.to_pickle}

    pl.parallelize(sdf_to_fingerprint, input_sdfs, opt.n_workers, fp_list = fp_options, urdkit_voc = opt.voc, format_dict = format_options )
                              

if __name__ == "__main__":
    main()
