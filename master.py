#!/usr/bin/env python

import numpy
import scipy
from scipy.optimize import linprog

def masterproblem(u,N,D,cost,xmin,xmax):

    Aeq = [1] * N
    beq = D
    Aineq = []
    bineq = []

    for i in range(0,N):
        bineq.append(-u[i]*xmin[i])
        bineq.append(u[i]*xmax[i])
        A_1 = [0] * N
        A_2 = [0] * N
        A_1[i] = -1.0
        A_2[i] =  1.0
        Aineq.append(A_1)
        Aineq.append(A_2)

    res = linprog(cost,A_ub=Aineq,b_ub=bineq,A_eq=[Aeq],b_eq=[beq])

    return res
# Run master problem
#res = master(u,N,D,c,xmin,xmax)

              


