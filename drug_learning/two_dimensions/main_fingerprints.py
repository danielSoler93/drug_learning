import argparse
import drug_learning.two_dimensions.Input.fingerprints as fp

def parse_arguments():
    parser = argparse.ArgumentParser(description="Specifying the fingerprints to which the user wants the sdf file to be converted.")

    parser.add_argument('-i', '--input',
                        dest="infile",
                        action="store",
                        default = ".",
                        help="Input sdf file")

    parser.add_argument('-mo', '--morgan',
                        dest ="morgan",
                        default = False,
                        help = "Convert molecules to Morgan fingerprint",
                        action = "store_true")

    parser.add_argument('-ma', "--maccs",
                        dest = "maccs",
                        default = False,
                        help = "Convert molecules to MACCS fingerprint",
                        action = "store_true")

    parser.add_argument('-rd', "--rdkit",
                        dest = "rdkit",
                        help = "Convert molecules to RDkit fingerprint",
                        default = False,
                        action = "store_true")

    parser.add_argument('-md', "--mordred",
                        dest = "mordred",
                        help = "Convert molecules to Mordred fingerprint",
                        default = False,
                        action ='store_true')

    options = parser.parse_args()

    return options

def main():
    options = parse_arguments()

    if options.morgan:
        morgan_fps = fp.MorganFP()
        morgan_fps.fit(options.infile)
        morgan_fps.transform()
        morgan_fps.save(to_csv=False, to_parquet=False, to_feather=False, to_hdf=False, to_pickle=True)

    if options.maccs:
        maccs_fps = fp.MACCS_FP()
        maccs_fps.fit(options.infile)
        maccs_fps.transform()
        maccs_fps.save(to_csv=False, to_parquet=False, to_feather=False, to_hdf=False, to_pickle=True)

    if options.rdkit:
        rdkit_fps = fp.RDkit_FP()
        rdkit_fps.fit(options.infile)
        rdkit_fps.transform()
        rdkit_fps.save(to_csv=False, to_parquet=False, to_feather=False, to_hdf=False, to_pickle=True)

    if options.mordred:
        mordred_fps = fp.MordredFP()
        mordred_fps.fit(options.infile)
        mordred_fps.transform()
        mordred_fps.save(to_csv=False, to_parquet=False, to_feather=False, to_hdf=False, to_pickle=True)

if __name__ == "__main__":
    main()
