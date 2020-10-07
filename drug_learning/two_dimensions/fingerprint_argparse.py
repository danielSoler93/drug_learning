import drug_learning.two_dimensions.Input.fingerprints as fp

def parse_args(parser):
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
