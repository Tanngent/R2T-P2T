# -*- coding: utf-8 -*-
import getopt
import math
import sys
import random
import time



def ReadInput(input_file_path):
    real_query_result = 0
    size_dic = {}
    input_file = open(input_file_path,'r')
    for line in input_file.readlines():
        elements = line.split()
        value = float(elements[0])
        entity = int(elements[1])
        real_query_result += value
        if entity in size_dic.keys():
            size_dic[entity] += value
        else:
            size_dic[entity] = value
    return size_dic, real_query_result
    
    
    
def LapNoise():
    a = random.uniform(0,1)
    b = math.log(1/(1-a))
    c = random.uniform(0,1)
    if c>0.5:
        return b
    else:
        return -b



def RunAlgorithm(size_dic, global_sensitivity, epsilon, beta):
    base = 5.5
    max_i = int(math.log(global_sensitivity,base))+1
    if max_i<=1:
        max_i+=1
    max_res1 = -10000000
    for i in range(1,max_i+1):
        tau = math.pow(base,i)
        t_res2 = LP(size_dic, tau)+LapNoise()*math.pow(base,i)/epsilon*max_i
        t_res1 = t_res2 - math.pow(base,i)/epsilon*max_i*math.log(max_i/beta,2.9718)
        if t_res1>max_res1:
            max_res1 = t_res1
    return max_res1
    
    

def LP(size_dic, tau):
    res = 0
    for element in size_dic.keys():
        res += min(tau,size_dic[element])
    return res



def runR2T(input_file_path = "../Information/TPCH/Q18_0.txt", global_sensitivity = 5000, utility = 1000, beta = 0.1):
    base = 5.5
    epsilon = 4 * math.log(global_sensitivity, base) * math.log(math.log(global_sensitivity, base) / beta) * global_sensitivity / utility
    size_dic, real_query_result = ReadInput(input_file_path)
    res = RunAlgorithm(size_dic, global_sensitivity, epsilon, beta)
    print(res, real_query_result - res, epsilon)
    return res, real_query_result - res, epsilon

if __name__ == "__main__":
	runR2T()