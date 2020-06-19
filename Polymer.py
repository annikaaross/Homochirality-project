#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 17:33:09 2020

@author: anni
"""
import Monomer
import globalvariables
import random
import numpy as np
import copy

#from google.colab import widgets
import matplotlib.pyplot as plt
import jdc
import more_itertools

#This is our Polymer class, aka the relationships in this dating show. They contain a list
# of all the monomers that are inside of them. It is kept in a specific order. They are 
#pretty open books. They can tell us how many monomers they contain, they add new members
#(either monomers or whole other polymers), they can remove entire groups of members 
#after a certain index, they can check which of their bonds are homochiral and use this
# to calculate probability of a certain bond breaking(the weak areas in their relationship).
# They can use this to then tell the producers if and where they are planning on breaking 
#off members of they ranks. Like all good relationships, they use purely mathematics and 
#chance to determine these heart breaking decisions. They can also tell us how many left 
#handed monomers and how many right handed monomers they contain.
class Polymer:
    
    #Constructor
    def __init__(self, monomers = []):
        """ Constructor method for Polymer.
        In addition to setting the monomers list, also calculates all break probabilities
        of the bonds it contains using reset_break_probability.
        """
        self._monomers = monomers
        self.reset_break_probability()

    def __str__(self):
        """ str function for Polymer class. """
        return f"\nA polymer of length {self.get_length()}: {self.get_sequence()}"
  
    def __repr__(self):
        """ repr function for Polymer class. """
        return f"{self.get_sequence()}-polymer"
    #Getters and Setters
    def get_monomers(self):
        """ Return list of monomers in the polymer. """
        return self._monomers


    def get_length(self):
        """ Return the number of monomers in the polymer. """
        return len(self._monomers)


    def get_sequence(self):
        """ Return the sequence of monomers in the polymer in a human-readable (and loggable) format. """
        sequence = []
        for monomer in self.get_monomers():
            sequence.append(monomer.get_handedness())
        return sequence


    def set_monomers(self, new_monomers):
        """ Set the monomers in the polymer to new_monomers. """
        self._monomers = new_monomers


    def get_chirality(self):
        """ Return the chiralities of the monomers within the polymer.

        RETURN: tuple
        (number_lefts, number_rights)
        """
        #Get the monomer handednesses
        sequence = self.get_sequence()
        #Count the number of trues
        n_true = sum(sequence)
        #The number of falses is the total minus the number of trues
        #Return that.
        return (n_true, len(sequence) - n_true)

    def get_signed_ee(self):
        lefts, rights = self.get_chirality()
        return (lefts - rights) / self.get_length()

    def get_ee(self):
        lefts, rights = self.get_chirality()
        return abs((lefts - rights) / self.get_length())

    def get_leftrighthomochiralbonds(self):
        """
        returns  # of left and right homochiral bonds in tuple (l,r)
        """
        l=0
        r=0
        for n in range(self.get_length()):
            if (self.easthomochiralbiascheck(n) == 3):
                l+=1
            elif (self.easthomochiralbiascheck(n) == 5):
                r+=1
        tuple1=(l,r)
        return tuple1
  
    def generate_old_log(self):
        """ Return a log tuple in the format
        (length, signed ee, total # homochiral bonds)
         """
        return (self.get_length(), 
        self.get_signed_ee())
          #sum(self.get_leftrighthomochiralbonds()))
  
    def generate_log(self):
        """ Return a log tuple formatted simply as the list of boolean handednesses of the monomers contained
        in the polymer. Will be heavier than holding data about the polymer like in the other log format,
        but should allow access to more information and more kinds of information.
        """
        return tuple(self.get_sequence())
    
    #FUnctionality Methods
    def append(self, other):
        """
        Adds a reactable (monomer or polymer) to the end of the Polymer.
        If the added reactable is a polymer, append() deletes it after taking its monomers.
        """
        #Check the type of thing to add
        if isinstance(other, Monomer):
            #Monomers can just be appended to the list
            self._monomers.append(other)
        elif isinstance(other, Polymer):
            #For polymers we need to extend rather than append, to avoid nested lists
            self._monomers.extend(other.get_monomers())
            #Then once the monomers are safely in their new polymer, the old one is deleted
            del other

    def brkloc(self):

        #makes a list of all the indices of monomers in the polymer except 
        #for the last one
        indices = []
        for n in range(self.get_length() - 1):
            indices.append(n)
            #puts indices in random order
        random.shuffle(indices)
        #tests if each indices will break
        for index in indices:
        #gets the breakprobability for monomer at n index
            brkprob = self._monomers[n].get_eastbrkprob()
            rand = random.random()
            # checks if the random number is less than the brkprob then the index of
            #the monomer is returned from the method
            #indicating this polymer breaks at the east bond of this monomer
            if(brkprob > rand):
                return n
        return None

    def removeright(self,location):
        """
        Removes all monomers to the right (east) of a given index
        """
        #creates a list of monomers in the polymer to the left, and including,
        #the monomer at the location(index) provided
        newList = self._monomers[0:location+1]

        #makes a list of the monoers to the rightof the monomer at the 
        #index provided
        removed = self._monomers[location+1:]

        #sets the polymers monomers to newList(the monomers to the left
        #of  and the index provided)
        self.set_monomers(newList)
        #print("after breaking ractable is:",self)
        #print("after breaking removed is:",removed)

        #returns the monomers that were removed, everything the right
        #of the index provided
        return removed
    
    #Homochirality Handling
    def easthomochiralcheck(self,numbermonomer):
        """
        Takes the index of a monomer within the Polymer and returns whether its east bond is homochiral.
        """
        #First check that the index isn't out of bounds 
        if (0 > numbermonomer or numbermonomer >= self.get_length()-1):
            #Should this raise an IndexError instead of returning False?
            return False
        #get the handedness of this monomer and its east neighbor, and return whether or not they're equal (bool)
        return self._monomers[numbermonomer].get_handedness() == self._monomers[numbermonomer+1].get_handedness()

    def easthomochiralbiascheck (self, numbermonomer):
        """ this method can be used in place of the east homochiral check method
        returns false if not homochiral, 3 if homochiral left, 5 if homochiral right"""
        #First check that the index isn't out of bounds 
        if (0 > numbermonomer or numbermonomer >= self.get_length()-1):
            #Should this raise an IndexError instead of returning False?
            return False
        #get the handedness of this monomer and its east neighbor, and return whether or not they're equal (bool)
        if (self._monomers[numbermonomer].get_handedness() == self._monomers[numbermonomer+1].get_handedness()):
            if self._monomers[numbermonomer].get_handedness():
                return 3
            elif not self._monomers[numbermonomer].get_handedness():
                return 5
        return False

    def eastbondbreakprobability(self,numbermonomer):
        """
        Takes the index of a monomer within the Polymer
        Returns the probability that the monomer's east bond will break
        returns -3 if the monomer has no east bond
        """
        ##############################
        #Override to return an equal break chance for every bond in the polymer
        #Simply comment out this line to get the regular function
        #return (0.3)
        ##############################
        #Initialize vars
        brk = -3
        length = self.get_length()
        #First check if the monomer is the last one in the polymer
        if (numbermonomer >= length-1):
        #If so, it doesn't have an east bond to break, so the probability is 0
            return brk
        #With that out of the way, we set brk to the base probability of breaking a bond
        #set to base bondbreak  multiplied by length factor to the power of 
        #the length of the polymer(breaking becomes less likely as the
        #polymer becomes longer)
        brk = globalvariables.BASE_BOND_BREAK_PROBABILITY * (globalvariables.LENGTH_FACTOR**self.get_length())

        #check if the east bond is homochiral
        if (self.easthomochiralcheck(numbermonomer)):
            brk *= globalvariables.HOMOCHIRAL_BREAK_FACTOR
            j = numbermonomer+1

        #going right to check for neighboring homochirality
        while(self.easthomochiralcheck(j)):
            brk *= 1 - (globalvariables.HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR**abs(j - numbermonomer))
            j += 1

        #going left to check for neighboring homochirality
        j=numbermonomer-1
        while(self.easthomochiralcheck(j)):
            brk *= 1-(globalvariables.HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR**abs(numbermonomer - j))
            j -= 1
        return brk

    def biaseastbondbreakprobability(self,numbermonomer):
        """
        ***left-right sensitivity****
        Takes the index of a monomer within the Polymer
        Returns the probability that the monomer's east bond will break
        returns -3 if the monomer has no east bond
        """
        ##############################
        #Override to return an equal break chance for every bond in the polymer
        #Simply comment out this line to get the regular function
        #return (0.3)
        ##############################
        #Initialize vars
        brk = -3
        length = self.get_length()
        #First check if the monomer is the last one in the polymer
        if (numbermonomer >= length-1):
            #If so, it doesn't have an east bond to break, so the probability is 0
            return brk
        #With that out of the way, we set brk to the base probability of breaking a bond
        #set to base bondbreak  multiplied by length factor to the power of 
        #the length of the polymer(breaking becomes less likely as the
        #polymer becomes longer)
        brk = globalvariables.BASE_BOND_BREAK_PROBABILITY * (globalvariables.LENGTH_FACTOR**self.get_length())

        #Check if the east bond is not homorchiral
        if (not self.easthomochiralbiascheck(numbermonomer)):
            return brk

        #check if the east bond is homochiral left
        if (self.easthomochiralbiascheck(numbermonomer)== 3):
            brk *= globalvariables.HOMOCHIRAL_BREAK_FACTOR_LEFT
            j = numbermonomer+1

            #going right to check for neighboring lefts
            while(self.easthomochiralcheck(j)):
                brk *= 1 - (globalvariables.HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR_LEFT**abs(j - numbermonomer))
                j += 1

            #going left to check for neighboring lefts
            j=numbermonomer-1
            while(self.easthomochiralcheck(j)):
                brk *= 1-(globalvariables.HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR_LEFT**abs(numbermonomer - j))
                j -= 1

        #checks if bond is homochiral right
        elif (self.easthomochiralbiascheck(numbermonomer) == 5):
            brk *= globalvariables.HOMOCHIRAL_BREAK_FACTOR_RIGHT
            j = numbermonomer+1

            #going right checking for neighboring rights
            while(self.easthomochiralcheck(j)):
                brk *= 1 - (globalvariables.HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR_RIGHT**abs(j - numbermonomer))
                j += 1

            #going left checking for neighboring lefts
            j=numbermonomer-1
            while(self.easthomochiralcheck(j)):
                brk *= 1-(globalvariables.HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR_RIGHT**abs(numbermonomer - j))
                j -= 1
        return brk

    def poissonbreakprobability(self,numbermonomer):
        """
        Takes the index of a monomer within the Polymer
        Returns the probability that the monomer's east bond will break
        based on Poisson distribution
        returns -3 if the monomer has no east bond
        """
        #Initialize brk
        brk = -3
        length = self.get_length()
        #First check if the monomer is the last one in the polymer
        if (numbermonomer >= length-1):
            #If so, it doesn't have an east bond to break, so the probability is 0
            return brk
        #With that out of the way, we set brk to the base probability of breaking a bond
        #set to base bondbreak  multiplied by length factor to the power of 
        #the length of the polymer(breaking becomes less likely as the
        #polymer becomes longer)
        brk = globalvariables.BASE_BOND_BREAK_PROBABILITY * (globalvariables.LENGTH_FACTOR**self.get_length())
        #keeps track of the number of homochiral bonds
        poisson=0 
        #check if the east bond is homochiral
        if (self.easthomochiralcheck(numbermonomer)):
            poisson+=1
            brk *= 1-self.poissonequation(poisson)
            j = numbermonomer+1
            #going right to check for neighboring homochirality
            while(self.easthomochiralcheck(j)):
                poisson+=1
                brk *= 1 - (self.poissonequation(poisson))
                j += 1
                #going left to check for neighboring homochirality
            j=numbermonomer-1
            while(self.easthomochiralcheck(j)):
                poisson+=1
                brk *= 1-(self.poissonequation(poisson))
                j -= 1
        return brk

    def poissonequation(self,k):
        #checks if value has already been calulated
        if (k not in globalvariables.LOOKUP_TABLE):
            #calculates
            p= np.exp(-globalvariables.LAMBDA)*((globalvariables.LAMBDA**k)/(np.math.factorial(k)))
            #adds to dictionary
            globalvariables.LOOKUP_TABLE[k]=p
        return globalvariables.LOOKUP_TABLE[k]



    def reset_break_probability(self):
        """ recalculates break probabilities for all monomers in polymer 
        can be altered to incorporate left right homochiral sensitivity"""
        #if the polymer is of length one(meaning it is about to be deleted
        #and simply recorded as a monomer) it sets its monomers 
        #brk probability to =-1
        if self.get_length() == 1 :
            self._monomers[0].set_eastbrkprob(-1)
            #otherwise goes through indices and recalcultes break prob
        else:
            for n in range(self.get_length()):
                ####################################################################
                #LEFT RIGHT SENSITIVITY CHOICE
                #comment out biaseastbondbreak to get no left/right sensitivity.   #
                #comment out eastbonbreakprobability to get left/right sensirtivity#
                ####################################################################
                #self._monomers[n].set_eastbrkprob(self.biaseastbondbreakprobability(n))
                self._monomers[n].set_eastbrkprob(self.eastbondbreakprobability(n))
                #calculates poisson brk probabilities with poisson distribution
                #self._monomers[n].set_eastbrkprob(self.poissonbreakprobability(n))
    
  
    
    
    