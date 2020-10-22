from multiprocessing import Pool
from functools import partial
import ray
import drug_learning.two_dimensions.Helpers.split as sp
import drug_learning.two_dimensions.Helpers.sdf_to_fingerprint as sf

def parallelize(func, iterable, n_workers, **kwargs):
    f = partial(func, **kwargs)
    if n_workers > 1:
        pool = Pool(n_workers)
        return pool.map(f, iterable)
        pool.close()
        pool.join()
    else:
        return list(map(f, iterable))

@ray.remote
def split_in_chunks_ray(sdf_file, n_chunks):
    split_files = sp.split_in_chunks(sdf_file, n_chunks)
    return split_files

@ray.remote
def sdf_to_fingerprint_ray(input_file, fp_list, format_dict, urdkit_voc):
    sf.sdf_to_fingerprint(input_file, fp_list, format_dict, urdkit_voc)
