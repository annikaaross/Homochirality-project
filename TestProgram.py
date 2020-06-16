#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 11:53:02 2020

@author: anni
"""
import Reactable
import Polymer
import Monomer
import globalvariables
import random
import numpy as np
import copy
from google.colab import widgets
import matplotlib.pyplot as plt
import jdc
import more_itertools

POOL_SIZE = 100
ITERATIONS = 100


a = Reactables(make_pool(POOL_SIZE))

print("Running.............|\n", end="")

monomercounts = []

for n in range(ITERATIONS):
    a.iterate(POOL_SIZE)
    monomercounts.append(len(a.get_free_monomers())+len(a.get_bound_monomers()))
    if n % (ITERATIONS / 10) == 0:
        print("##",end="")

print("\nDone.")

plots = ["parameters","polymer ee spread","ee heatmap","leftright chirality",
         "homochiral chain lengths", "homochirality vs length"]
tb = widgets.TabBar(plots)

with tb.output_to("parameters"):
    a.visualize_parameters()

with tb.output_to("polymer ee spread"):
    a.plot_signed_ee_spread()

with tb.output_to("ee heatmap"):
    a.plot_ee_heatmap()

with tb.output_to("leftright chirality"):
    a.plot_leftrighthomochirality()

with tb.output_to("homochiral chain lengths"):
    a.plot_homochiral_chain_lengths()

with tb.output_to("homochirality vs length"):
    a.plot_homochirality_vs_length()
