#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 16:55:42 2020

@author: anni
"""

#the bond break probbility applied to all bonds
BASE_BOND_BREAK_PROBABILITY = 1
#if a bond itself is homochiral the base bond fator is multiplied by
#this factor
HOMOCHIRAL_BREAK_FACTOR = 0.9
#this factor is incorporated for each additional neighboring homochiral bond
# becomes less and less influential as the bond gets further away
HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR = 0.3
#decreases the break probability based on the length of the polymer
LENGTH_FACTOR = 0.9

#dictionary of poisson values
LOOKUP_TABLE={}
#used in Poisson distribution
LAMBDA=6
#################################################################
#these are the same factors as before but with on for left and  #
#one for right homochiral bonds.                                #
#used in biaseastbondbreakability method which can be choosen to#
#be used over the regular eastbondbreakability method in the    #
#self.reset_break_probability() method                          #
#################################################################
HOMOCHIRAL_BREAK_FACTOR_LEFT = 0.3
HOMOCHIRAL_BREAK_FACTOR_RIGHT = 0.9
HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR_LEFT = 0.3
HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR_RIGHT = 0.3


#These are the liklihoods that a monomer will choose to poof or bond, respectively.
#Used in the Reactables class.
POOF_CHANCE = 0.3333
BOND_PROB = 0.3333