from multiprocessing import Pool
from functools import partial

def parallelize(func, iterable, n_workers, **kwargs):
    pool = Pool(n_workers)
    f = partial(func, **kwargs)
    return pool.map(f, iterable)
    pool.close()
    pool.join()
