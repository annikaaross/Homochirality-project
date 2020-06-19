#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 16:58:31 2020

@author: anni
"""

import globalvariables
import random
import numpy as np
import copy

#from google.colab import widgets
import matplotlib.pyplot as plt
import jdc
import more_itertools

# These are our Monomers. The individuals in our dating show. They are the best thing a girl
# can be in this worl, beautiful little fools. They know that they exist and they know
# what kind of monomer they are, left or right(thats handedness not political party 
#alignment) What more can you ask for?


class Monomer:
    
    
# Constructor
    def __init__(self, **kwargs):
        """
        The constructor method for Monomer objects. Assigns default handedness.
        Accepts argument hand = bool with keyword.
        Considering other syntax for passing in handedness.
        """
        self._handedness = kwargs['hand'] if 'hand' in kwargs else random.choice([True,False])
        self._eastbrkprob = kwargs['brkprob'] if 'brkprob' in kwargs else -1
  
    def __repr__(self):
        return str(f"{self.get_handedness()}-monomer")

    # Getters and Setters
    def get_handedness(self):
        """Getter method for Monomer handedness."""
        return self._handedness

    def get_eastbrkprob(self):
        """Getter method for east bond brk probability"""
        return self._eastbrkprob

    def set_eastbrkprob(self, newbrk):
        """Setter method for east bond brk probability"""
        self._eastbrkprob = newbrk

    def reset_eastbrkprob(self) :
        """reset east brk probability to -1 
        to be called when single monomer is broken from polymer"""
        self.set_eastbrkprob(-1)
   
    def get_chirality(self):
        """ Return an expression of handedness compatible with the Polymer notation
        of chirality (n_lefts, n_rights).
        """
        if self.get_handedness():
            return (1,0)
        else:
            return (0,1)

    def get_length(self):
        """ Return 1. """
        return 1

    def get_ee(self):
        """ Return the enantiomeric excess of the monomer. This value is always 1. """
        return 1

    def get_signed_ee(self):
        """ Return the signed enantiomeric excess of the monomer.
        Returns 1 if left handed, -1 if right handed.
        """
        if self.get_handedness():
            return 1
        else:
            return -1

    def generate_old_log(self):
        """ Return a log tuple in the format
        (length, signed ee, total # homochiral bonds)
         """
        return (self.get_length(), 
        self.get_signed_ee())
        #sum(self.get_leftrighthomochiralbonds()))

    def generate_log(self):
        """ Return a log compatible with the reactables history. """
        return (self.get_handedness(),)
