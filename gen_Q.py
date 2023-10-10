import numpy as np

# def gen_Q(xs, cost_coeffs, delta_f):
#     n = len(xs)
#     c_tildes = cost_coeffs * xs
    
#     Q = np.zeros((n, n))

#     for i in range(n):
#         for j in range(n):
#             if i == j:
#                 Q[i][j] = -delta_f
#             else:
#                 Q[i][j] = c_tildes[i] * c_tildes[j]
#
#     return Q

def gen_Q(x_maxes):
    n = len(x_maxes)
    
    Q = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            if i == j:
                Q[i][j] = 0
            else:
                Q[i][j] = x_maxes[i] * x_maxes[j]
                
    # Rescale the matrix such that max value is 20
    max_value = np.max(Q)
    if max_value != 0:
        scale_factor = 20.0 / max_value
        Q *= scale_factor

    return Q
