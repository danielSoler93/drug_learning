import argparse
import drug_learning.two_dimensions.Input.fingerprints as fp
import drug_learning.two_dimensions.Errors.errors as er

def parse_arguments():
    parser = argparse.ArgumentParser(description="Specifying the fingerprints to which the user wants the sdf file to be converted.")

    parser.add_argument('-i', '--input',
                        dest="infile",
                        action="store",
                        help="Input sdf file")

    parser.add_argument('-mo', '--morgan',
                        dest ="morgan",
                        action = "store_true",
                        default = False,
                        help = "Convert molecules to Morgan fingerprint")

    parser.add_argument('-ma', "--maccs",
                        dest = "maccs",
                        action = "store_true",
                        default = False,
                        help = "Convert molecules to MACCS fingerprint")

    parser.add_argument('-rd', "--rdkit",
                        dest = "rdkit",
                        action = "store_true",
                        default = False,
                        help = "Convert molecules to RDkit fingerprint")

    parser.add_argument('-md', "--mordred",
                        dest = "mordred",
                        action ='store_true',
                        default = False,
                        help = "Convert molecules to Mordred fingerprint")

    parser.add_argument('-urd', "--unfolded_rdkit",
                        dest = "unfolded_rdkit",
                        action ='store_true',
                        default = False,
                        help = "Convert molecules to Unfolded RDkit fingerprint")

    parser.add_argument('-voc', "--vocabulary",
                        dest = "voc",
                        type = str,
                        default = None,
                        help = "Vocabulary for unfolded rdkit fingerprint")

    parser.add_argument('-c', "--csv",
                        dest = "to_csv",
                        action = "store_true",
                        default = False,
                        help = "Save output to csv")

    parser.add_argument('-pq', "--parquet",
                        dest = "to_parquet",
                        action = "store_true",
                        default = True,
                        help = "Save output to parquet")

    parser.add_argument('-f', "--feather",
                        dest = "to_feather",
                        action ='store_true',
                        default = False,
                        help = "Save output to feather")

    parser.add_argument('-hdf', "--hdf",
                        dest = "to_hdf",
                        action = "store_true",
                        default = False,
                        help = "Save output to hdf")

    parser.add_argument('-pk', "--pickle",
                        dest = "to_pickle",
                        action ='store_true',
                        default = False,
                        help = "Save output to pickle")


    options = parser.parse_args()

    return options

def main():
    opt = parse_arguments()
    if opt.morgan:
        morgan_fps = fp.MorganFP()
        morgan_fps.fit(opt.infile)
        morgan_fps.transform()
        morgan_fps.save(to_csv=opt.to_csv, to_parquet=opt.to_parquet, to_feather=opt.to_feather, to_hdf=opt.to_hdf, to_pickle=opt.to_pickle)

    if opt.maccs:
        maccs_fps = fp.MACCS_FP()
        maccs_fps.fit(opt.infile)
        maccs_fps.transform()
        maccs_fps.save(to_csv=opt.to_csv, to_parquet=opt.to_parquet, to_feather=opt.to_feather, to_hdf=opt.to_hdf, to_pickle=opt.to_pickle)

    if opt.rdkit:
        rdkit_fps = fp.RDkitFP()
        rdkit_fps.fit(opt.infile)
        rdkit_fps.transform()
        rdkit_fps.save(to_csv=opt.to_csv, to_parquet=opt.to_parquet, to_feather=opt.to_feather, to_hdf=opt.to_hdf, to_pickle=opt.to_pickle)

    if opt.mordred:
        mordred_fps = fp.MordredFP()
        mordred_fps.fit(opt.infile)
        mordred_fps.transform()
        mordred_fps.save(to_csv=opt.to_csv, to_parquet=opt.to_parquet, to_feather=opt.to_feather, to_hdf=opt.to_hdf, to_pickle=opt.to_pickle)

    if opt.unfolded_rdkit:
        if not opt.voc:
            raise er.NotVocabularyUnfolded("Vocabulary (--voc) must be pass to use unfolded rdkit fingerprints")
        mordred_fps = fp.UnfoldedRDkitFP(opt.voc)
        mordred_fps.fit(opt.infile)
        mordred_fps.transform()
        mordred_fps.save(to_csv=False, to_parquet=True, to_feather=False, to_hdf=False, to_pickle=False)

if __name__ == "__main__":
    main()
