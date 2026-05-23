# -*- coding: utf-8 -*-
import math
import random
import numpy as np
from R2T_PARAM import runR2T
from R2TSJF_PARAM import runR2TSJF
from P4T import runP4TPrivRelax

def LapNoise(sensitivity, epsilon):
    scale = sensitivity / epsilon
    u = random.uniform(0, 1)
    if u < 0.5:
        return scale * math.log(2 * u)
    else:
        return -scale * math.log(2 * (1 - u))
    
def R2T_Translate(global_sens, alpha, beta):
    eps = 4 * math.log(global_sens) * math.log(global_sens / beta) * global_sens / alpha
    # print(f"R2T_Translate: eps={eps}")
    return eps

def runP4TSJF(inp = "../Information/TPCH/Q18_0.txt", global_sens = 50000, utility = 100000, bins = 10, beta = 0.1, eps_h = 1):
    eps_q = 10000
    beta_h = beta / (2*bins)
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
    # print(unoised_hist)
    noise_scale = 2 * sum(unoised_hist) + 1
    print(noise_scale)
    noised_hist = {k: max(0, v + LapNoise(noise_scale, eps_h) + noise_scale/eps_h*math.log(bins/2/beta_h)) for k, v in unoised_hist.items()}
    # print(noised_hist)
    # errors = []
    for i in range(bins):
        tau = int(global_sens / bins * (i+1))
        bias = sum([max(k - tau, 0) * v for k, v in noised_hist.items()])
        # print(bias, 2 * tau ** 2 / eps_h ** 2, bias + 2 * tau ** 2 / eps_h ** 2)
        # errors[i] = bias + 2 * tau ** 2 / eps_h ** 2
        t_tau = utility - bias
        if t_tau > 0:
            #print(tau / t_tau * math.log(1 / beta_q))
            if eps_q > tau / t_tau * math.log(1 / beta_q):
                eps_q = tau / t_tau * math.log(1 / beta_q)
                tau_q = tau
                #print(tau, eps_q)
    clipped_query = sum([min(v, tau_q) for k, v in size_dic.items()])
    release = clipped_query + LapNoise(tau_q, eps_q)
    # print(f"#bins={bins} utility={utility} release={release} diff={real_query_result - release} eps={eps_h + eps_q} truncation= {tau_q}")
    return release, abs(real_query_result - release), eps_h + eps_q, tau_q

def rcond_y_given_x_lap_paper(x, eps1, eps2, size=1):
    """
    Sample the second noise V2 = y given the first noise x and privacy levels eps1 < eps2.
    """
    assert eps2 > eps1 > 0, "Require eps2 > eps1 > 0"
    
    p1 = (eps1 / eps2) * np.exp(-(eps2 - eps1) * abs(x))
    p2 = (eps2 - eps1) / (2 * eps2)
    p3 = (eps1 + eps2) / (2 * eps2) * (1 - np.exp(-(eps2 - eps1) * abs(x)))
    p4 = (eps2 - eps1) / (2 * eps2) * np.exp(-(eps2 - eps1) * abs(x))

    # print(p1,p2,p3,p4)
    samples = []

    for _ in range(size):
        r = np.random.rand()
        # Case 1
        if r < p1:
            # print('case1', x)
            y = x
        
        # Case 2
        elif r < p1 + p2:
            z = np.random.exponential(scale=1/(eps1+eps2)) #* (eps1+eps2)
            z = -z  # z <= 0
            # print('case2', np.sign(x) * z)
            y = np.sign(x) * z
        
        # Case 3
        elif r < p1 + p2 + p3:
            z = np.random.exponential(scale=1/(eps2 - eps1)) #* ((eps2 - eps1))
            # truncate to [0, |x|]
            z = min(z, abs(x))
            # print('case3', np.sign(x) * z)
            y = np.sign(x) * z
        
        # Case 4
        else:
            z = np.random.exponential(scale=1/(eps1+eps2)) #* (eps1+eps2)
            # shift to >= |x|
            z = z + abs(x)
            # print('case4', np.sign(x) * z)
            y = np.sign(x) * z
        samples.append(y)
    return samples

def runP4TSJFPrivRelax(inp = "../Information/TPCH/Q18_0.txt", global_sens = 5000, bins = 10, alpha = 10, beta = 0.1, eps_max = 20, kappa = 2):
    # starting historgram parameters
    eps_h = eps_max / math.pow(kappa, math.ceil(math.log(eps_max, kappa)))
    beta_h = beta / (2*bins)

    # starting query parameters
    eps_q = eps_max + 1
    beta_q = beta / 2

    # calculate real results and distrubution
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
    
    # calcualte unoised histogram
    unoised_hist = {}
    for i in range(bins):
        unoised_hist[int(global_sens / bins * (i+1))] = 0
    for k, v in size_dic.items():
        v_bin = int(math.ceil(v / (global_sens / bins)) * (global_sens / bins))
        unoised_hist[v_bin] += 1
    #print(unoised_hist)

    # calcualte starting noised histogram
    noise_scale = 2 * sum(unoised_hist) + 1
    print(noise_scale)
    noised_hist = {k: max(0, v + LapNoise(noise_scale, eps_h) + noise_scale/eps_h*math.log(bins/2/beta_h)) for k, v in unoised_hist.items()}

    # finding suitable eps_h
    while True:
        #print(eps_h, eps_q)
        for i in range(bins):
            tau = int(global_sens / bins * (i+1)) # for each tau
            bias = sum([max(k - tau, 0) * v for k, v in noised_hist.items()])
            t_tau = alpha - bias # leftover budget
            if t_tau > 0:
                #print(tau / t_tau * math.log(1 / beta_q))
                if eps_q > tau / t_tau * math.log(1 / beta_q): # find min eps_q
                    eps_q = tau / t_tau * math.log(1 / beta_q)
                    tau_q = tau
        #print(noised_hist)
        if eps_h + eps_q < eps_max:
            break # found suitable combination
        if eps_h > eps_max:
            raise  RuntimeError("no answer found")
        noised_hist = {k: unoised_hist[k] + rcond_y_given_x_lap_paper(v - unoised_hist[k], eps_h, eps_h * kappa)[0] for k, v in noised_hist.items()}
        eps_h = eps_h * kappa
        
    clipped_query = sum([min(v, tau_q) for k, v in size_dic.items()])
    release = clipped_query + LapNoise(tau_q, eps_q)
    #print(f"#bins={bins} alpha={alpha} release={release} diff={real_query_result - release} eps={eps_h + eps_q} truncation={tau_q}")
    print(real_query_result - release, eps_h + eps_q)
    return release, abs(real_query_result - release), eps_h + eps_q, tau_q

if __name__ == "__main__":
    """
    bins = [10, 20, 100, 500, 5000, 50000, 500000, 5000000, 50000000]
    # bins=[10]
    
    inp = "../Information/TPCH/Q18_0.txt"
    
    global_sens = 5000
    utility = 150000
    beta = 0.01
    times = 10
    results = []
    
    for bin in bins:
        if bin > global_sens:
            continue
        temp = []
        for t in range(times):
            release, diff, eps, tau_q = runP4T(inp=inp, bins=bin, global_sens = global_sens, utility = utility, beta = beta, eps_h = 0.5)
            temp.append((release, diff, eps, tau_q))
        avg_release = sum(x[0] for x in temp) / len(temp)
        avg_diff = sum(x[1] for x in temp) / len(temp)
        avg_eps = sum(x[2] for x in temp) / len(temp)
        avg_tau_q = sum(x[3] for x in temp) / len(temp)
        results.append((bin, avg_release, avg_diff, avg_eps, avg_tau_q))
    
    print("bin, release, diff, eps, tau_q", results)
    
    R2T_results = []
    for t in range(times):
        eps = R2T_Translate(global_sens = global_sens, alpha = utility, beta = beta) 
        R2T_results.append(eps)
    avg_R2T_eps = sum(R2T_results) / len(R2T_results)
    print("R2T eps:", avg_R2T_eps)
    """

    global_sens = 5000
    alpha = 10000
    beta = 0.01
    kappa = 2
    
    bins = 20
    r2t_eps = R2T_Translate(global_sens, alpha, beta)

    print("R2T_SJF")
    runR2TSJF(global_sensitivity=global_sens, epsilon=r2t_eps, beta=beta)
    print("R2T")
    runR2T(global_sensitivity=global_sens, epsilon=r2t_eps, beta=beta)
    print("P4T_SJF")
    runP4TSJFPrivRelax(global_sens=global_sens, bins=bins, alpha=alpha, beta=beta, kappa=kappa)
    print("P4T")
    runP4TPrivRelax(global_sens=global_sens, bins=bins, alpha=alpha, beta=beta, kappa=kappa)