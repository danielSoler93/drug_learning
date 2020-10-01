import drug_learning.two_dimensions.Errors.errors as ce
import drug_learning.two_dimensions.Input.base_class as bc

def test_output_format_error():
    try:
        fp = bc.Fingerprint()
        fp.save()
    except ce.NotOutputFormat:
        return True
    return False
        
