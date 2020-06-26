import multiprocessing  
import json
from frap import FRAP
import numpy as np

def loguniform(low=0, high=1, size=None):
    return np.exp(np.random.uniform(low, high, size))

def uniform(low=0, high=1, size=None):
    return np.random.uniform(low, high, size)

def run(frap: FRAP, **params):

    # Parameter ranges
    D       = loguniform(low=params["lb_D"], high=params["ub_D"])
    c0      = uniform(low=params["lb_c"], high=params["ub_c"])
    alpha   = uniform(low=params["lb_alpha"], high=params["ub_alpha"])

    # Variance ranges
    a = loguniform(low=params["lb_a"], high=params["ub_a"])
    b = loguniform(low=params["lb_b"], high=params["ub_b"])

    # Other 
    mobile_fraction = 1
    beta            = 1
    gamma           = 0

    sys_params = {"D": D,
                  "c0": c0,
                  "alpha": alpha,
                  "a": a,
                  "b": b,
                  "mobile_fraction": mobile_fraction,
                  "beta": beta,
                  "gamma": gamma}

    return frap.generate(**sys_params)

if __name__ == "__main__":

    with open('configs/test.json') as f:
        config = json.load(f)

    frap = FRAP(**config["experiment"])

    pool = multiprocessing.Pool(processes=4)
    pool.map(run(frap, **config["ranges"]), range(10))
    pool.close()
    pool.join()   
    print('done')
