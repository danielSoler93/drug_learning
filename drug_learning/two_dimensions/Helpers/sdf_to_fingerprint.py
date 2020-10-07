import drug_learning.two_dimensions.Input.fingerprints as fp

def sdf_to_fingerprint(input_file, fp_list, format_dict, urdkit_voc=None):
    for fp_class in fp_list:
       if fp_class:
            if fp_class == fp.UnfoldedRDkitFP:
                fingerprint = fp_class(urdkit_voc)
            else:
                fingerprint = fp_class()
            fingerprint.fit(input_file)
            fingerprint.transform()
            fingerprint.save(**format_dict)
