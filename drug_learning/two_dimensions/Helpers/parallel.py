from multiprocessing import Pool
from functools import partial

def parallelize(func, iterable, n_workers, **kwargs):
    f = partial(func, **kwargs)
    if n_workers > 1:
        pool = Pool(n_workers)
        return pool.map(f, iterable)
        pool.close()
        pool.join()
    else:
        return list(map(f, iterable))
