import os
import pytest
import numpy as np
import drug_learning.two_dimensions.main_fingerprints as mn

DIR = os.path.dirname(__file__)
MOLECULE1 = os.path.join(DIR, "data/test_ad_molecule1.sdf")
MOLECULE2 = os.path.join(DIR, "data/test_ad_molecule2.sdf")

def test_AD_fit(get_fingerprint, result = [0., 0.]):
    AD = mn.ad.ApplicabilityDomain()
    fingerprint = get_fingerprint(mn.fp.MorganFP, MOLECULE1)
    AD.fit(fingerprint.features)
    assert np.all(AD.thresholds == np.array(result))

@pytest.mark.parametrize( "molecule, other_molecule, result",
                            [
                            (MOLECULE1, MOLECULE1, [2,2] ),
                            (MOLECULE1, MOLECULE2, [0,0] ),
                            ]
                        )
def test_AD_predict(get_fingerprint, molecule, other_molecule, result):
    AD = mn.ad.ApplicabilityDomain()
    fingerprint = get_fingerprint(mn.fp.MorganFP, molecule)
    other_fingerprint = get_fingerprint(mn.fp.MorganFP, other_molecule)
    AD.fit(fingerprint.features)
    AD.predict(other_fingerprint.features)
    assert np.all(AD.n_insiders == np.array(result))
