import math
import random
import numpy as np
import matplotlib.pyplot as plt

def LapNoise(sensitivity, epsilon):
    scale = sensitivity / epsilon
    u = random.uniform(0, 1)
    if u < 0.5:
        return scale * math.log(2 * u)
    else:
        return -scale * math.log(2 * (1 - u))

def testHistDiff(inp = "../Information/TPCH/Q18_0.txt", global_sens = 5000, bins = 10, alpha = 10, beta = 0.1, eps_max = 20, kappa = 2, override_size_dic = None):
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
    print(max(size_dic, key=size_dic.get))

    if override_size_dic is not None:
        size_dic = override_size_dic
    
    # calcualte unoised histogram
    unoised_hist = {}
    for i in range(bins):
        unoised_hist[int(global_sens / bins * (i+1))] = 0
    for k, v in size_dic.items():
        v_bin = int(math.ceil(v / (global_sens / bins)) * (global_sens / bins))
        unoised_hist[v_bin] += 1
    print(unoised_hist)
    # calcualte starting noised histogram
    noised_hist = {k: max(0, v + LapNoise(1, eps_h) + 1/eps_h*math.log(bins/2/beta_h)) for k, v in unoised_hist.items()}
    print(noised_hist)
    t, r, u, n = [], [], [], []
    for i in range(bins):
        tau = int(global_sens / bins * (i+1)) # for each tau
        real_truncation_error = 0
        for k, v in size_dic.items():
            real_truncation_error += max(0, v - tau)
        unoised_truncation_error = sum([max(k - tau, 0) * v for k, v in unoised_hist.items()])
        noised_truncation_error = sum([max(k - tau, 0) * v for k, v in noised_hist.items()])
        t.append(tau)
        r.append(real_truncation_error)
        u.append(unoised_truncation_error)
        n.append(noised_truncation_error)
    return t, r, u, n

if __name__ == "__main__":

    import sys
    print(sys.prefix)

    global_sens = 50000
    s = np.random.beta(4, 4, 10000)
    #s = np.random.beta(8, 2, 10000)
    #s = np.random.beta(2, 8, 10000)
    #s = np.random.rand(10000)
    override_size_dic = {i: max(1, min(v*global_sens, global_sens)) for i, v in enumerate(s)}
    plt.hist(override_size_dic.values(), bins=100)
    plt.savefig(f'dist.png')
    plt.clf()

    bins = [10, 20, 30, 40, 50]
    for bin in bins:
        t, r, u, n = testHistDiff(bins=bin, global_sens=global_sens, override_size_dic=override_size_dic)
        #print(r)
        #print(u)
        #print(n)
        #if bin == 100:
            #plt.plot(t, r, label ='real')
        #plt.plot(t, u, '-.', label ='unoised')
        plt.plot(t, [i - j for i, j in zip(n, r)], '-.', label =f'noised_{bin}')
    plt.legend(ncol=3)
    plt.savefig(f'my_plot_normal.png')