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

gc.enable()
class Optimizer(cplex.callbacks.SimplexCallback):
    def __call__(self):
        value = self.get_objective_value()
        if self.LP_type == 0:
            self.primals[self.tau_index] = value
            # Update the threshold for early stop
            if self.global_max.value <= self.factor + value:
                self.global_max.value = self.factor + value
            # If arrive the stop condition for the approximate algorithm
            # If the dual/primal has already stopped
            # If arrive the early stop
            if (self.stop_duals[self.tau_index] == 1 or 
                abs(self.duals[self.tau_index] - value) <= self.tau or 
                self.early_stops[self.tau_index] == 1):
                self.abort()
        else:
            self.duals[self.tau_index] = value
            if value + self.factor < self.global_max.value:
                self.early_stops[self.tau_index] = 1
                self.abort()
            if (self.stop_primals[self.tau_index] == 1 or 
                abs(self.primals[self.tau_index] - value) <= self.tau):
                self.abort()

def LapNoise(sensitivity, epsilon):
    scale = sensitivity / epsilon
    u = random.uniform(0, 1)
    if u < 0.5:
        return scale * math.log(2 * u)
    else:
        return -scale * math.log(2 * (1 - u))

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

def ThresholdRunAlgorithm(base, downward_sensitivity, assigned_of_tau_ids, LP_type, 
                          tilde_Q_tau, entities, connections, aggregation_values, 
                          approximate_factor, primals, duals, stop_primals, 
                          stop_duals, early_stops, global_max):
    """Process assigned tau values."""
    for i in assigned_of_tau_ids:
        tau = math.pow(base, i)
        if tau < downward_sensitivity:
            LPSolver(tau, LP_type, i, tilde_Q_tau[i], entities, connections, 
                    aggregation_values, approximate_factor, primals, duals, 
                    stop_primals, stop_duals, early_stops, global_max)


def LPSolver(tau, LP_type, tau_index, factor, entities, connections, 
             aggregation_values, approximate_factor, primals, duals, 
             stop_primals, stop_duals, early_stops, global_max):
    """Solve the LP for a given tau value."""
    num_constraints = len(entities)
    num_variables = len(connections)
    
    # Set the obj
    cpx = cplex.Cplex()
    cpx.objective.set_sense(cpx.objective.sense.maximize)
    
    # Set variables
    obj = np.ones(num_variables)
    ub = np.zeros(num_variables)
    for i in range(num_variables):
        ub[i] = aggregation_values[i]
    cpx.variables.add(obj=obj, ub=ub)
    
    # Set the right hand side and the sign
    rhs = np.ones(num_constraints) * tau
    senses = "L" * num_constraints
    cpx.linear_constraints.add(rhs=rhs, senses=senses)
    
    # Set the coefficients
    cols = []
    rows = []
    vals = []
    for i in range(num_variables):
        for j in connections[i]:
            cols.append(i)
            rows.append(j)
            vals.append(1)
    cpx.linear_constraints.set_coefficients(zip(rows, cols, vals))
    
    cpx.set_log_stream(None)
    cpx.set_error_stream(None)
    cpx.set_warning_stream(None)
    cpx.set_results_stream(None)
    
    # Set the optimizer
    if LP_type == 0:
        cpx.parameters.lpmethod.set(1)
    else:
        cpx.parameters.lpmethod.set(2)
    
    optimizer = cpx.register_callback(Optimizer)
    optimizer.threshold = tau * approximate_factor
    optimizer.tau = tau
    optimizer.LP_type = LP_type
    optimizer.tau_index = tau_index
    optimizer.factor = factor
    optimizer.primals = primals
    optimizer.duals = duals
    optimizer.stop_primals = stop_primals
    optimizer.stop_duals = stop_duals
    optimizer.early_stops = early_stops
    optimizer.global_max = global_max
    
    cpx.solve()
    
    # If the get the feasible solution, update the information
    if LP_type == 0 and cpx.solution.get_status() == cpx.solution.status.optimal:
        primals[tau_index] = cpx.solution.get_objective_value()
        if global_max.value <= factor + primals[tau_index]:
            global_max.value = factor + primals[tau_index]
        stop_primals[tau_index] = 1
        if global_max.value <= factor + primals[tau_index]:
            global_max.value = factor + primals[tau_index]
    elif LP_type == 1 and cpx.solution.get_status() == cpx.solution.status.optimal:
        duals[tau_index] = cpx.solution.get_objective_value()
        if global_max.value <= factor + duals[tau_index]:
            global_max.value = factor + duals[tau_index]
        stop_duals[tau_index] = 1
        if global_max.value <= factor + duals[tau_index]:
            global_max.value = factor + duals[tau_index]

def ReadInput(input_file_path):
    """Read input file and return all data structures."""
    # Store the ids of entities
    entities = []
    # The connections between entities and join results
    connections = []
    # The aggregation values of join results
    aggregation_values = []
    # The dictionary to store the tuples' sensitivities
    entities_sensitivity_dic = {}
    # The dictionary to re-id entities
    id_dic = {}
    # The number of base table tuples
    id_num = 0
    # Collect the DS
    downward_sensitivity = 0
    
    input_file = open(input_file_path, 'r')
    for line in input_file.readlines():
        elements = line.split()
        connection = []
        # The first value is the aggregation value
        aggregation_value = float(elements[0])
        # For each entity contribution to that join result
        for element in elements[1:]:
            element = int(element)
            # Re-order the IDs
            if element in id_dic.keys():
                element = id_dic[element]
            else:
                entities.append(id_num)
                id_dic[element] = id_num
                element = id_num
                id_num += 1
            # Update the entity's sensitivity
            if element in entities_sensitivity_dic.keys():
                entities_sensitivity_dic[element] += aggregation_value
            else:
                entities_sensitivity_dic[element] = aggregation_value
            # Update the DS
            if downward_sensitivity <= entities_sensitivity_dic[element]:
                downward_sensitivity = entities_sensitivity_dic[element]
            connection.append(element)
        connections.append(connection)
        aggregation_values.append(aggregation_value)
    
    input_file.close()
    real_query_result = sum(aggregation_values)
    
    return entities, connections, aggregation_values, downward_sensitivity, real_query_result

def runP4TPrivRelax(inp = "../Information/TPCH/Q18_0.txt", global_sens = 5000, bins = 20, alpha = 150000, beta = 0.1, eps_max = 20, kappa = 2):
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
    print(noised_hist)
    # finding suitable eps_h
    while True:
        #print(eps_h, eps_q)
        for i in range(bins):
            tau = int(global_sens / bins * (i+1)) # for each tau
            bias = sum([max(k - tau, 0) * v for k, v in noised_hist.items()])
            t_tau = alpha - bias # leftover budget
            #print(t_tau)
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
    
    primals = [0]
    duals = [10 * real_query_result]
    stop_duals = [0]
    stop_primals = [0]
    early_stops = [0]
    global_max = Value(c_double, -global_sens / eps_q, lock=True)
    entities, connections, aggregation_values, downward_sensitivity, real_query_result = ReadInput(inp)
    LPSolver(tau_q, 1, 0, 0, entities, connections, aggregation_values, 0, primals, duals, stop_primals, stop_duals, early_stops, global_max)
    if tau_q > downward_sensitivity:
        clipped_query = real_query_result
    elif stop_duals[0] == 1:
        clipped_query = duals[0]
    else:
        clipped_query = primals[0]
    release = clipped_query + LapNoise(tau_q, eps_q)
    #print(f"#bins={bins} alpha={alpha} release={release} diff={real_query_result - release} eps={eps_h + eps_q} truncation={tau_q}")
    print(real_query_result - release, eps_h + eps_q)
    return release, abs(real_query_result - release), eps_h + eps_q, tau_q

if __name__ == "__main__":
    runP4TPrivRelax()