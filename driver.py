#!/usr/bin/env python

import numpy as np
import scipy
from scipy.optimize import linprog
from parseprovider import parseprovider
from master import masterproblem
from computeQUBO import computeQUBO
from gen_Q import gen_Q

# Step 1 => Parse file and get the input parameters

input_file = 'data/edisp_10.txt';
N,D,cost,xmin,xmax = parseprovider(input_file)

# Step 2 call the QUBO problem
status = [1] * N
# Assume it gives commitment u
result = masterproblem(status,N,D,cost,xmin,xmax)

x = result.x

Q = gen_Q(xmax)

computeQUBO(Q)
    




