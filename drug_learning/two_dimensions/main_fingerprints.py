import argparse
import glob
from itertools import chain
import drug_learning.two_dimensions.Input.applicability_domain as ad
import drug_learning.two_dimensions.Helpers.split as sp
import drug_learning.two_dimensions.Helpers.parallel as pl
import drug_learning.two_dimensions.Helpers.sdf_to_fingerprint as sf
import drug_learning.two_dimensions.Errors.errors as er
import drug_learning.two_dimensions.fingerprint_argparse as fp_arg
import drug_learning.two_dimensions.split_argparse as sp_arg

def get_parser():
    parser = argparse.ArgumentParser(description="Specify the fingerprints and the file format to which the user wants the sdf file to be converted.")
    fp_arg.parse_args(parser)
    sp_arg.parse_args(parser)
    return parser

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

    pl.parallelize(sf.sdf_to_fingerprint, input_sdfs, opt.n_workers, fp_list = fp_options, urdkit_voc = opt.voc, format_dict = format_options )

if __name__ == "__main__":
    main()
