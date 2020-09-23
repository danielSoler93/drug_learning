import pytest

@pytest.fixture
def get_fingerprint():
    def _get_fingerprint(fp_class, input_sdf):
        fingerprint = fp_class()
        fingerprint.fit(input_sdf)
        features = fingerprint.transform()
        return fingerprint
    return _get_fingerprint
