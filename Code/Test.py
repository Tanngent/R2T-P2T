# -*- coding: utf-8 -*-
import getopt
import math
import sys
import random
import cplex
import numpy as np
import time
import gc
import multiprocessing
from multiprocessing.sharedctypes import Value
from ctypes import c_double
manager = multiprocessing.Manager()
gc.enable()

def LapNoise(sensitivity, epsilon):
    scale = sensitivity / epsilon
    u = random.uniform(0, 1)
    if u < 0.5:
        return scale * math.log(2 * u)
    else:
        return -scale * math.log(2 * (1 - u))

def runP4T(inp = "../Information/TPCH/Q18_0.txt", global_sens = 5000, utility = 1000, bins = 10, beta = 0.1, eps_h = 100):
    eps_q = 10000
    beta_h = beta / 2
    beta_q = beta / 2
    tau_q = 0

    real_query_result = 0
    size_dic = {}
    input_file = open(inp,'r')
    for line in input_file.readlines():
        elements = line.split()
        value = float(elements[0])
        entity = int(elements[1])
        real_query_result += value
        if entity in size_dic.keys():
            size_dic[entity] += value
        else:
            size_dic[entity] = value
    #print(size_dic)
    unoised_hist = {}
    for i in range(bins):
        unoised_hist[int(global_sens / bins * (i+1))] = 0
    for k, v in size_dic.items():
        v_bin = int(math.ceil(v / (global_sens / bins)) * (global_sens / bins))
        unoised_hist[v_bin] += 1
    #print(unoised_hist)
    noised_hist = {k: max(0, v + LapNoise(1, eps_h) + 1/eps_h*math.log(bins/2/beta_h)) for k, v in unoised_hist.items()}
    #print(noised_hist)
    for i in range(bins):
        tau = int(global_sens / bins * (i+1))
        bias = sum([max(k - tau, 0) * v for k, v in noised_hist.items()])
        t_tau = utility - bias
        if t_tau > 0:
            #print(tau / t_tau * math.log(1 / beta_q))
            if eps_q > tau / t_tau * math.log(1 / beta_q):
                eps_q = tau / t_tau * math.log(1 / beta_q)
                tau_q = tau
                #print(tau, eps_q)
    clipped_query = sum([min(v, tau_q) for k, v in size_dic.items()])
    release = clipped_query + LapNoise(tau_q, eps_q)
    print(release, real_query_result - release, eps_h + eps_q)
    return release, real_query_result - release, eps_h + eps_q

if __name__ == "__main__":
	runP4T()