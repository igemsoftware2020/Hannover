from math import sqrt
from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt

def brownian(x0, n, dt, delta):
    x0 = np.asarray(x0)

    # For each element of x0, generate a sample of n numbers from a
    # normal distribution.
    r = norm.rvs(size=x0.shape + (n,), scale=delta*sqrt(dt))
    out = np.empty(r.shape)

    np.cumsum(r, axis=-1, out=out)

    # Add the initial condition.
    out += np.expand_dims(x0, axis=-1)
    return out


# The Wiener process parameter.
delta = 3
# Total time.
T = 10.0
# Number of steps.
N = 500
# Time step size
dt = T/N
# Number of realizations to generate.
m = 1
# Create an empty array to store the realizations.
x = np.empty((m,N+1))
# Initial values of x.
x[:, 0] = 0

out = brownian(x[:,0], N, dt, delta)
print(out)