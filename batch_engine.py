# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 12:03:10 2020

@author: Lio
"""

# Import

import chirality_sim as cs
import numpy as np
run_num = 0

# Define

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 50, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    Source: https://stackoverflow.com/a/34325723
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()

def batch_recur(folder, batch, variables, ranges, reps=1, depth=0):
    """
    Recursively run a batch of simulations.
    
    parameters:
    
    variables (list):
    A list of string literals referring to the parameters that should vary in the 
      batch. The strings should be keys in the params dictionary.
    
    ranges (list):
    A list of lists containing the values each corresponding parameter should take
      during the batch.
    
    reps (int, optional):
    The number of times to run each combination of parameter values. Defaults to 1.
    
    depth (int, optional):
    The current layer of recursion, used within the function to keep track of its
      current status. Defaults to 0, don't change it.
    """
    # for each value that this particular parameter should take,
    if ranges == []:
        raise ValueError("Need to choose at least one parameter to be dynamic!!")
    for value in ranges[depth]:
        # Set the parameter to that value.
        params[variables[depth]] = value
        # Run the sim with these settings reps times.
        for n in range(reps):
            # But hang on, is it time to run the sim yet?
            if depth < len(variables)-1:
                # If we haven't recursed enough to set every variable, we need to go deeper.
                batch_recur(folder, batch, variables, ranges, reps, depth+1)
            else:
                # But if we have set everything, we can run the sim and step back up a level.
                simulate(folder,batch)
                # Update progress bar
                global run_num
                run_num += 1
                printProgressBar(run_num,n_runs)


def simulate(folder,batch):
    # Run sim

    sim = cs.Simulation(params)
    sim.run()
    
    # Export data and params
    
    sim.export(f"{folder}/{batch}")


##################################################################

# Set base params

params = {"BASE_BOND_BREAK_PROBABILITY":0.5, 
          "HOMOCHIRAL_BREAK_FACTOR":0.5,
          "HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR":0.8,
          "LENGTH_FACTOR":0.6,
          "N":40,
          "LAMBDA":6,
          "HOMOCHIRAL_BREAK_FACTOR_LEFT":0.5,
          "HOMOCHIRAL_BREAK_FACTOR_RIGHT":0.5,
          "HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR_LEFT":0.5,
          "HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR_RIGHT":0.5,
          "POOF_CHANCE":0.33,
          "BOND_PROB":0.33,
          "POOL_SIZE":500,
          "ITERATIONS":500,
          "METHOD":'standard',
          "POISSON_FACTOR":1.3,
          "FUSION":False,
          "REFILL":["refill normal"]}


# Set which parameters to change

variable_params = ["HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR", "BASE_BOND_BREAK_PROBABILITY"]

# Choose what values to give them

var_ranges = [[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]]

# Name folder and batch

folder = "data"
batch = "hnif-bbbp"

# How many times to run each combination

reps = 1


#################################################################

# Figure out how many runs it's gonna take
run_lens = [len(r) for r in var_ranges]
n_combs = np.prod(run_lens)
n_runs = n_combs * (reps**len(run_lens))
print(f"{n_runs} total runs")
printProgressBar(run_num,n_runs)

# Here it goes 
batch_recur(folder,batch,variable_params,var_ranges,reps)























