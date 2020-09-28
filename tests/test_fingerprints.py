import os
import pytest
import drug_learning.two_dimensions.Input.fingerprints as fp

DIR = os.path.dirname(__file__)
MOLECULE = os.path.join(DIR, "data/test_molecule.sdf")
NAME = "core,42852325"

def test_fp_fit(molecule = MOLECULE, length = 1, name = NAME):
    fingerprint = fp.bc.Fingerprint()
    fingerprint.fit(molecule)
    assert len(fingerprint.structures) == length
    assert fingerprint.structures[0].GetProp("_Name") == name

@pytest.mark.parametrize( "fp_class, shape",
                            [
                            (fp.MorganFP, (1, 2048)),
                            (fp.MACCS_FP, (1, 167)),
                            (fp.RDkitFP, (1, 2048)),
                            (fp.MordredFP, (1, 1613)) # mirar error
                            ]
                        )
def test_fingerprint_transform(get_fingerprint, fp_class, shape, molecule = MOLECULE):
    fingerprint = get_fingerprint(fp_class, molecule)
    assert fp.np.shape(fingerprint.features) == shape
    assert len(fingerprint.mol_names) == len(fingerprint.features)
    assert len(fingerprint.columns) == fp.np.shape(fingerprint.features)[1]

@pytest.mark.parametrize( "csv, parquet, feather, hdf, pickle, extension",
                            [
                            (True, False, False, False, False, ".csv"),
                            (False, True, False, False, False, ".gzip"),
                            (False, False, True, False, False,".ftr"),
                            (False, False, False, True, False, ".hdf5"),
                            (False, False, False, False, True, ".pkl"),
                            ]
                        )
def test_fingerprint_save(get_fingerprint, csv, parquet, feather, hdf, pickle, extension, molecule = MOLECULE):
    fingerprint = get_fingerprint(fp.MorganFP, molecule)
    fingerprint.save(to_csv= csv, to_parquet= parquet, to_feather= feather, to_hdf= hdf, to_pickle= pickle)
    file = fingerprint.filename + fingerprint.fp_name + extension
    assert os.path.exists(file)
    os.remove(file)
