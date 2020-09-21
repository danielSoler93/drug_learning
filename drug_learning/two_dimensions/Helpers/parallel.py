from multiprocessing import Pool
from functools import partial

def parallelize(func, iterable, n_workers, **other_args):
    pool = Pool(n_workers)
    f = partial(func, **other_args)
    return pool.map(f, iterable)
    pool.close()
    pool.join()
