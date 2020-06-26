# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 11:42:54 2020

@author: Lio
"""
import numpy as np
import pandas as pd
import more_itertools

# Shorthands #
L = True
R = False

# Column labels #
Type = "Type"                   # Monomer or polymer
Length = "Length"               # How many monomers in the thing
nLefts = "#Lefts"               # How many left monomers in the thing
nRights = "#Rights"             # How many right monomers in the thing
Sequence = "Sequence"           # The string sequence of the thing
nLhomo = "#LeftHomochiral"      # The number of left homochiral bonds in the polymer
nRhomo = "#RightHomochiral"     # The number of right homochiral bonds in the polymer
nHomo = "#Homochiral"           # The overall number of homochiral bonds in the polymer
sEE = "Signed ee"               # The signed enantiomeric excess of the polymer (+ if more True)
pcHomo = "%Homochirality"       # The proportion of bonds in the polymer that are homochiral
pcLhomo = "%LeftHomochirality"  # The proportion of bonds in the polymer that are left homochiral
pcRhomo = "%RightHomochirality" # The proportion of bonds in the polymer that are right homochiral




lookup = {}


hist = [[(False,), (True,), (False,), (False,), (True,), (True,), (False, True), (True, False)], 
        [(False, True, True, False), (True,), (False,), (False,), (False,), (False,), (True,), (True, True)], 
        [(True, True), (False,), (True,), (False,), (False, True, True), (False,), (False, False), (True, False), (False,)], 
        [(True, True, False, True, True), (False,), (False,), (False, False), (False,), (True,), (True, False), (False, False)], 
        [(False, False, True, False), (True,), (True,), (False, False), (False,), (False,), (False,), (False, True, True, False, True, True)], 
        [(False, False), (False, False, True, False, True), (False,), (False,), (False,), (False,), (True,), (True,), (False, True, True, False, True, True)], 
        [(False, True, True, False, True, True, True), (True,), (False, False, True, False, True), (False,), (False,), (False,), (False, False), (True, False)], 
        [(False, False, True, False, True), (True,), (False, True, True, False, True, True, True), (False,), (False, False), (False,), (False,), (False,), (True, False), (False,)], 
        [(False,), (False, False, False), (False,), (True,), (False, True, True, False, True, True, True), (False, False, True, False, True), (True, False), (False,), (False,)], 
        [(False,), (True,), (False, False, False), (True, False, False, True, True, False, True, True, True), (False,), (False, False), (False, False, False, True, False, True)]]

stats = []


def parse_history(history):
    """ Create an array of plottable information from the history log. """

    for iteration in history:

        iter_stats = pd.DataFrame()
        for item in iteration:
            if len(item) == 1: # It's a monomer
                
                #Translate
                sequence = standard_form(item)
                
                #Is it in the lookup table?
                if sequence not in lookup:
                    
                    # Get the info
                    lr = count_LR(item)
                    
                    # And put it in the lookup table
                    lookup[sequence]={Type: 'Monomer', 
                                     Length:1, 
                                     nLefts:lr[0], 
                                     nRights:lr[1],
                                     Sequence:sequence}
                    
                # Now that the data is searchable...
                # Log it
                new_log = pd.Series(lookup.get(sequence))
                iter_stats = iter_stats.append(new_log, ignore_index=True)

            elif len(item) > 1: # It's a polymer
                
                #Translate
                sequence = standard_form(item)
                
                #Is it in the lookup table?
                if sequence not in lookup:

                    # Get basic info
                    length = len(item)
                    lr = count_LR(item)
                    bonds = homochiral_bond_counts(item)
                    total_homos = bonds[0]+bonds[1]
                    signed_ee = (lr[0]-lr[1])/(lr[0]+lr[1])
                    homochirality = total_homos/(length-1)
                    lhomochirality = bonds[0]/(length-1)
                    rhomochirality = bonds[1]/(length-1)
                    
                    # And put it in the lookup table
                    lookup[sequence] = {Type: 'Polymer', 
                                     Length: length, 
                                     nLefts:lr[0], 
                                     nRights:lr[1],
                                     nLhomo:bonds[0],
                                     nRhomo:bonds[1],
                                     nHomo:total_homos,
                                     sEE:signed_ee,
                                     pcHomo:homochirality,
                                     pcLhomo:lhomochirality,
                                     pcRhomo:rhomochirality,
                                     Sequence:sequence}
                
            
                # Now that the data is searchable...
                # Log it
                new_log = pd.Series(lookup.get(sequence))
                iter_stats = iter_stats.append(new_log, ignore_index=True)
                
            else:
                raise ValueError("There's something with length 0 in your history.")
        stats.append(iter_stats)

    
# Functions for getting the info we want to log
    
# Number of right and left monomers in a polymer (works on monomers too)
def count_LR(log):
    """ Return (n_Lefts, n_Rights) """
    if L and not R:
        return (sum(log), len(log)-sum(log))
    elif R and not L:
        return (len(log)-sum(log), sum(log))
    else:
        raise ValueError("Your L's and R's are screwed up somehow.")
    

def eAnd(*args): # From https://stackoverflow.com/q/2770434
    return [all(tuple) for tuple in zip(*args)]
    


# HOMOCHIRALITY SEQUENCES #
def homochirality_sequence(log):
  """ Return a boolean list of bonds within a given logged polymer, True if the bond is
      homochiral. """
  return [i[0]==i[1] for i in list(more_itertools.pairwise(log))]

def _is_L_True(log):
    return log[:-1]

def T_homochirality_sequence(h_seq,L_seq):
    """ Return a boolean list of bonds within a polymer, true if the bond is 
    a homochiral bond between two True monomers. Parameters are the outputs of 
    homochirality_sequence() and _is_L_True(). """
    return eAnd(h_seq,L_seq)

# HOMOCHIRAL BOND COUNTS #
    
def homochiral_bond_counts(log):
    """ Return (number of left homochiral bonds, number of right homochiral bonds) """
    homo = homochirality_sequence(log)
    west_true = log[:-1]
    true_homochiral = eAnd(homo,west_true)
    west_false = [not m for m in west_true]
    false_homochiral = eAnd(homo,west_false)
    if L and not R:
        return (sum(true_homochiral), sum(false_homochiral))
    elif R and not L: 
        return (sum(false_homochiral), sum(true_homochiral))
    else:
        raise ValueError("Your L's and R's are screwed up somehow.")


def get_polymer_chain_lengths(polylog):
  """ Return the lengths of the homochiral chains in a given polymer log. """
  count = 1
  lengths = [] 
  for n in range(1,len(polylog)):
    if polylog[n] == polylog[n-1]:
      count += 1
    else:
      lengths.append(count)
      count = 1
  lengths.append(count)
  return lengths
    
def standard_form(poly):
  if len(poly) == 0:
    return ''
  return ''.join(['L' if m else 'R' for m in poly])
    
    

##############################################################################
#%% MAIN %%#
##############################################################################

parse_history(hist)
