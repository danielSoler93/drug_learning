import os
import drug_learning.two_dimensions.Input.fingerprints as fp


def main(ligand_sdf):
    dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    input = os.path.join(dir, ligand_sdf)
    morgan_fps = fp.MorganFP()
    morgan_fps.fit(input)
    morgan_fps.structures
    morgan_fps.transform()
    morgan_fps.save(to_csv=False, to_parquet=False, to_feather=False, to_hdf=False, to_pickle=True)


if __name__ == "__main__":
    dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    input = os.path.join(dir, "datasets/ligands.sdf")
    main(input)
