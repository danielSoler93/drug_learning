import argparse
import glob
from itertools import chain
import ray
import os
import drug_learning.two_dimensions.Input.applicability_domain as ad
import drug_learning.two_dimensions.Helpers.split as sp
import drug_learning.two_dimensions.Helpers.sdf_to_fingerprint as sf
import drug_learning.two_dimensions.Helpers.parallel as pl
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

    if opt.ray:
        ray.init(address=os.environ.get("address"), _redis_password=os.environ.get("password"))

    if opt.split:
        if opt.ray:
            results = [pl.split_in_chunks_ray.remote(f, opt.n_chunks) for f in input_sdfs]
            input_sdfs = ray.get(results)
        else:
            input_sdfs = [sp.split_in_chunks(f, opt.n_chunks) for f in input_sdfs]

        input_sdfs = list(chain.from_iterable(input_sdfs))

    fp_options = [opt.morgan, opt.maccs, opt.rdkit, opt.mordred, opt.unfolded_rdkit]

    format_options = {'to_csv' : opt.to_csv, 'to_parquet' : opt.to_parquet, 'to_feather' : opt.to_feather,
                     'to_hdf' : opt.to_hdf, 'to_pickle' : opt.to_pickle}

    if any(fp_options):
        if opt.ray:
            results = [pl.sdf_to_fingerprint_ray.remote(f, fp_options, format_options, opt.voc) for f in input_sdfs]
            ray.get(results)
        else:
            for f in input_sdfs: sf.sdf_to_fingerprint(f, fp_options, format_options, opt.voc)

if __name__ == "__main__":
    main()
