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


def LapNoise():
    """Generate Laplace noise."""
    a = random.uniform(0, 1)
    b = math.log(1 / (1 - a))
    c = random.uniform(0, 1)
    if c > 0.5:
        return b
    else:
        return -b


def RunAlgorithm(entities, connections, aggregation_values, downward_sensitivity, 
                 real_query_result, global_sensitivity, epsilon, beta, 
                 approximate_factor, processor_num):
    """Run the R2T algorithm."""
    manager = multiprocessing.Manager()
    
    # Base of Log function
    # Besides, one notice is that we also optimize the base used in Theorem 5.1.
    # In the paper, we use 2 while it can be proven the optimal one is e (by optimizing the error bound in Theorem 5.1).
    # Here, we use 2e for its better practical performance.
    # Please note any theoretical result in the paper will not be affected as long as the base is still a constant number.
    base = 5.5
    # The number of all tau's
    max_i = int(math.log(global_sensitivity, base)) + 1
    if max_i <= 1:
        max_i += 1
    
    # Used to store the results
    # The Q(I,tau)
    Q_tau = {}
    # The Q(I,tau)-factor+Noise
    tilde_Q_tau = {}
    # The Q(I,tau)+Noise
    hat_Q_tau = {}
    
    # Initialize the variables for optimizations
    primals = manager.dict()
    duals = manager.dict()
    stop_primals = manager.dict()
    stop_duals = manager.dict()
    early_stops = manager.dict()
    global_max = Value(c_double, -global_sensitivity / epsilon, lock=True)
    
    # Assign the tau's
    arrangement_of_tau_ids = []
    for i in range(processor_num):
        arrangement_of_tau_ids.append([])
    
    # Indicate which processor assigned to one given tau/i
    j = 0
    # Enumerate i inversely
    for ii in range(max_i):
        i = max_i - ii
        tau = math.pow(base, i)
        arrangement_of_tau_ids[j].append(i)
        j = (j + 1) % processor_num
        Q_tau[i] = 0
        # Add the noise
        hat_Q_tau[i] = LapNoise() * tau / epsilon * max_i * (approximate_factor + 1)
        # Add the noise and the factor
        tilde_Q_tau[i] = hat_Q_tau[i] - tau / epsilon * max_i * math.log(max_i / beta) * (approximate_factor + 1)
        primals[i] = 0
        duals[i] = 10 * real_query_result
        stop_duals[i] = 0
        stop_primals[i] = 0
        early_stops[i] = 0
    
    # Create the processors
    threads = []
    for i in range(processor_num):
        # Try to make parameters locally
        threads.append(multiprocessing.Process(
            target=ThresholdRunAlgorithm,
            args=(base, downward_sensitivity, arrangement_of_tau_ids[i], 1, tilde_Q_tau,
                  entities, connections, aggregation_values, approximate_factor,
                  primals, duals, stop_primals, stop_duals, early_stops, global_max)
        ))
        threads[i].start()
    
    for i in range(processor_num):
        threads[i].join()
    
    # Obtain the tau with maximum tilde_Q_tau
    max_ind = 1
    max_val = 0
    for i in range(1, max_i + 1):
        tau = math.pow(base, i)
        if tau >= downward_sensitivity:
            Q_tau[i] = real_query_result
        else:
            if stop_duals[i] == 1:
                Q_tau[i] = duals[i]
            else:
                Q_tau[i] = primals[i]
        hat_Q_tau[i] += Q_tau[i]
        tilde_Q_tau[i] += Q_tau[i]
        # The the the LP is early stoped, the result should not be counted
        if early_stops[i] == 1:
            continue
        if tilde_Q_tau[i] > max_val:
            max_val = tilde_Q_tau[i]
            max_ind = i
    
    final_res = tilde_Q_tau[max_ind]
    return final_res


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


def runR2T(input_file_path = "../Information/TPCH/Q18_0.txt", global_sensitivity = 5000, utility = 150000, epsilon = 1, beta = 0.01):
    # The number of processor
    processor_num = 10
    # The approximate factor
    approximate_factor = 0
    
    # Two processors for one task: primal and dual
    processor_num = int(processor_num)
    if processor_num < 1:
        processor_num = 1
    
    start = time.time()
    entities, connections, aggregation_values, downward_sensitivity, real_query_result = ReadInput(input_file_path)
    res = RunAlgorithm(entities, connections, aggregation_values, downward_sensitivity, 
                      real_query_result, global_sensitivity, epsilon, beta, 
                      approximate_factor, processor_num)
    print(real_query_result - res, epsilon)
    return real_query_result - res


if __name__ == "__main__":
    times = 10
    err = 0
    for t in range(times):
        err += abs(runR2T())
    
    print("average error:", err/times)