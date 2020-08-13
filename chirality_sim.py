# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 13:38:00 2020

@author: Lio

Local importable homochirality simulator

Not really a thing that works atm because of our reliance on global vars.
"""


#%%### IMPORTS ### 

import random
import numpy as np
import copy
import more_itertools
import pandas as pd
import uuid
import datetime
from collections import namedtuple

#%%### DEFINE RULES CLASS ###

Rules = namedtuple("Rules",['BASE_BOND_BREAK_PROBABILITY',
                            'HOMOCHIRAL_BREAK_FACTOR',
                            'HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR',
                            'LENGTH_FACTOR',
                            'N',
                            'HOMOCHIRAL_BREAK_FACTOR_LEFT',
                            'HOMOCHIRAL_BREAK_FACTOR_RIGHT',
                            'HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR_LEFT',
                            'HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR_RIGHT',
                            'POOF_CHANCE',
                            'BOND_PROB',
                            'POOL_SIZE',
                            'ITERATIONS',                        
                            'METHOD',
                            'LAMBDA',
                            'POISSON_FACTOR',
                            'FUSION',
                            'REFILLRANDOM',
                            'REFILLPERCENT',
                            'REFILLNUMBER',
                            'REFILLNUMBERDECREASE',
                            'REFILLNORMAL',
                            'REFILLPERCENTDECREASE'])




#%% GLOBALS

#initialize all global parameters - values will be replaced so they are unimportant
BASE_BOND_BREAK_PROBABILITY = 0.9
HOMOCHIRAL_BREAK_FACTOR = 0.9
HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR = 0.3
LENGTH_FACTOR = 0.6
N=40
LAMBDA=6
HOMOCHIRAL_BREAK_FACTOR_LEFT = 0.0
HOMOCHIRAL_BREAK_FACTOR_RIGHT = 0.9
HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR_LEFT = 0.3
HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR_RIGHT = 0.3
POOF_CHANCE = 0.3333
BOND_PROB = 0.3333
POOL_SIZE = 10
ITERATIONS = 10
METHOD = 'standard'
POISSON_FACTOR = 1.3
FUSION = False
REFILL = ["no refill"]
rands = []
poisson_dict = {}
break_prob_lookup_table = {}

L = True
R = False


#%%### HELPERS ###


def make_pool(n):
    """
    Return a list of n new monomers.
    """
    monomers = []
    for i in range(n):
        monomers.append(_Monomer())
    return monomers

def eAnd(*args): # From https://stackoverflow.com/q/2770434
    return [all(tuple) for tuple in zip(*args)]

def get_rand():

    """A function to provide uniform random numbers"""

    if len(rands)<= 10:
        #refills the array if contains less than 10
        randsarray = np.random.rand(1000)
        rands.extend(randsarray.tolist())
    #return last number and deletes it from the array
    return rands.pop()

def poissonequation(k):

    """Given a value k, returns the result of the poisson equation for k"""

    #checks if value has already been calculated
    if k not in poisson_dict:
        #otherwise calculates the values
        p = np.exp(-LAMBDA)*((LAMBDA**k)/(np.math.factorial(k)))
        #adds to dictionary
        poisson_dict[k] = 1- 1.3*(p)
    return poisson_dict[k]
  
def dict_to_csv(params,filename):
    with open(filename, 'w') as f:
        for key in params.keys():
              f.write("%s,%s\n"%(key,params[key]))
    f.close()
    
    
#%%### MONOMER ###
    
class _Monomer:

    def __init__(self, **kwargs):
        """
        The constructor method for Monomer objects. Assigns default handedness.
        Accepts argument hand = bool with keyword.
        Unless provided, assigns easbrkprob to -1 (indicating it is not in a polymer).
        Initializes age, how many iterations the monomer has survived, as 0.
        """

        self._handedness = kwargs['hand'] if 'hand' in kwargs else random.choice([True,False])
        self._eastbrkprob = kwargs['brkprob'] if 'brkprob' in kwargs else -1
        self._age = 0
        self._uuid = uuid.uuid4()
  
    def __repr__(self):
        """repr function for Monomer Class"""
        return str(f"{self.get_handedness()}-monomer")
    
    def get_handedness(self):
        """Getter method for Monomer handedness."""
        return self._handedness
    
    def get_age(self):
        """Getter method for Monomer age"""
        return self._age
    
    
    def get_eastbrkprob(self):
        """Getter method for east bond brk probability"""
        return self._eastbrkprob
    
    def set_eastbrkprob(self, newbrk):
        """Setter method for east bond brk probability"""
        self._eastbrkprob = newbrk
    
    def set_age(self, newage):
        """setter method for age"""
        self._age = newage
    
    def reset_eastbrkprob(self) :
        """ reset east brk probability to -1 
            to be called when single monomer is broken from polymer"""
        self.set_eastbrkprob(-1)
    
    def get_id(self):
        return self._uuid
    
    def get_sequence(self):
        return (self.get_handedness(),)
    
    def birthday(self):
        """ages the monomer up"""
        self._age += 1
    

#%%### POLYMER ###
        
class _Polymer:
    
    def __init__(self, monomers = []):
        """ 
        Constructor method for Polymer. 
        Sets list of monomers it contains.
        Keeps track of the age of itself (number of iterations it has survived)
        """
        self._monomers = monomers
        self._age = 0
        self._uuid = uuid.uuid4()

    def __str__(self):
        """ str function for Polymer class. """
        return f"\nA polymer of length {self.get_length()}: {self.get_sequence()}"

    def __repr__(self):
        """ repr function for Polymer class. """
        return f"{self.get_sequence()}-polymer"
    
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
    
    def get_age(self):
        """
        Getter method for the age of the polymer
        """
        return self._age
    
    
    def set_age(self, newage):
        """
        Setter method for the age of the polymer
        """
        self._age = newage
    
    def set_monomers(self, new_monomers):
        """ Set the monomers in the polymer to new_monomers. """
        self._monomers = new_monomers
    
    def birthday(self):
        """
        Ages up the polymer and all of the monomers in the polymer
        """
        #ages up itself
        self._age += 1
        #ages up all its monomers
        for n in range(self.get_length()):
            self.get_monomers()[n].birthday()
    
    def get_id(self):
        return self._uuid
    
    def append(self, other):
        """
        Adds a reactable (monomer or polymer) to the end of the Polymer.
        If the added reactable is a polymer, this method deletes it after taking its monomers.
        """
        #Check the type of the thing to be added
    
        if isinstance(other, _Monomer):
            #Monomers can just be appended to the list
            self._monomers.append(other)
    
        elif isinstance(other, _Polymer):
            #For polymers we need to extend rather than append, to avoid nested lists
            self._monomers.extend(other.get_monomers())
            #Then once the monomers are safely in their new polymer, the old one is deleted
            del other
    
    def brkloc(self):
    
        """This method randomly goes through the bonds of a polymer
        and asks each if it will brk, returns first that will, or None if none break"""
    
        #list of all the indices of monomers in the polymer except for the last one
        indices = []
        for n in range(self.get_length() - 1):
            indices.append(n)
    
        #puts indices in random order
        random.shuffle(indices)
    
        #tests if each index will break
        for index in indices:
            #gets the breakprobability for monomer at n index
            brkprob = self._monomers[index].get_eastbrkprob()
            #pulls a random number
            rand = get_rand()
            # checks if the random< brkprob of index
            if(brkprob > rand):
                #if so returns index
                return index
        #if none break, return none
        return None
    
    def removeright(self, location):
        """
        Removes all monomers to the right (east) of a given index
        and returns a list of the monomers that were removed
        """
        #list of the monomers at the, and to the left of
        #the location provided
        newList = self._monomers[0:location+1]
    
        #list of the monomers to the right of
        #the location provided
        removed = self._monomers[location+1:]
    
        #resets monomers to only include those to the left
        self.set_monomers(newList)
    
        #returns the monomers that were removed
        return removed
    
# HOMOCHIRALITY HANDLING #
            
    def easthomochiralcheck(self,numbermonomer):
        """
        Takes the index of a monomer within the Polymer and returns whether its east bond is homochiral.
        """
        #First check that the index isn't out of bounds 
        if (0 > numbermonomer or numbermonomer >= self.get_length()-1):
            return False
        #get the handedness of this monomer and its east neighbor, and return whether or not they're equal (bool)
        return self._monomers[numbermonomer].get_handedness() == self._monomers[numbermonomer+1].get_handedness()
    
    def easthomochiralbiascheck (self, numbermonomer):
        """ this method can be used in place of the east homochiral check method
        returns false if not homochiral, 3 if homochiral left, 5 if homochiral right"""
    
        #First check that the index isn't out of bounds or if not homochiral
        if (0 > numbermonomer or numbermonomer >= self.get_length()-1 or not self.easthomochiralcheck(numbermonomer)):
            return False
    
        #otherwise check if they are homochiral left or right
        elif self._monomers[numbermonomer].get_handedness():
            #means left homochiral
            return 3
        elif not self._monomers[numbermonomer].get_handedness():
            #means right homochiral
            return 5
        
        
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
    
        #First check if the monomer is the last one in the polymer
        if (numbermonomer >= self.get_length()-1):
            #If so, it doesn't have an east bond to break, so the probability is 0(return -3)
            return -3
    
        #now we initialize brk probability (brk)
        brk = (BASE_BOND_BREAK_PROBABILITY) * (LENGTH_FACTOR**(self.get_length()/N))
        #check if the east bond is homochiral
        if (self.easthomochiralcheck(numbermonomer)):
            #if so multiply it by homochiral break factor (shrinks probability)
            brk *= HOMOCHIRAL_BREAK_FACTOR
            #goes through method which checks and calculates benefit of all homochiral neighbors
            brk = self.checkforhomochiralneighbors(numbermonomer,brk,HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR)
        #end def: returns break probability of monomers east bond
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
    
        #First check if the monomer is the last one in the polymer
        if (numbermonomer >= self.get_length()-1):
            #If so, it doesn't have an east bond to break, so the probability is 0 (returns -3)
            return -3
    
        #initialize the brk probability
        brk = (BASE_BOND_BREAK_PROBABILITY) * (LENGTH_FACTOR**(self.get_length()/N))
    
        #check if the east bond is homochiral left
        if (self.easthomochiralbiascheck(numbermonomer) == 3):
            brk *= HOMOCHIRAL_BREAK_FACTOR_LEFT
    
            #run through function that recalculates brk based on benfits of homochiral neighbors
            brk = self.checkforhomochiralneighbors(numbermonomer, brk, HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR_LEFT)
    
        #otherwise checks if bond is homochiral right
        elif (self.easthomochiralbiascheck(numbermonomer) == 5):
            brk *= HOMOCHIRAL_BREAK_FACTOR_RIGHT
    
            #run through function that recalculates brk based on benfits of homochiral neighbors
            brk = self.checkforhomochiralneighbors(numbermonomer,brk, HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR_RIGHT)
    
        return brk
    
    def checkforhomochiralneighbors(self, numbermonomer, brk, neighborfactor):
    
        """
        helps all the versions of the eastbond break calculating functions
        takes in a brk probability a number monomer and the desired neighbor improvement factor
        and calculates how the neighbors a bond will increase the strength of the bond
        it then returns the new brk prob to the eastbond brk function
        """
        brk = brk
        j = numbermonomer + 1
        #going right to check for homochiral neighbors
        while(self.easthomochiralcheck(j)):
            #calculates decrease to brk prob for a neighbor of that distance for the bond
            brk *= 1 - (neighborfactor**abs(j - numbermonomer))
            j += 1
    
        #going left to check for homochiral neighbors
        j=numbermonomer-1
        while(self.easthomochiralcheck(j)):
            #calculates decrease to brk prob for a neighbor of that distance for the bond
            brk *= 1-(neighborfactor**abs(numbermonomer - j))
            j -= 1
    
        return brk
    
    
    def poissonbreakprobability(self,numbermonomer):
        """
        Takes the index of a monomer within the Polymer
        Returns the probability that the monomer's east bond will break
        based on Poisson distribution
        returns -3 if the monomer has no east bond
        """
    
        #First check if the monomer is the last one in the polymer
        if (numbermonomer >= self.get_length()-1):
            #If so, it doesn't have an east bond to break, so the probability is 0(return -3)
            return -3
    
        #With that out of the way, we initialize brk to the base probability of breaking a bond
        brk = (BASE_BOND_BREAK_PROBABILITY) * (LENGTH_FACTOR**(self.get_length()/N))
    
        #keeps track of the number of the total consecutive homochiral bonds
        homochiralcount=0 
    
        #check if the east bond is homochiral
        if (self.easthomochiralcheck(numbermonomer)):
            #if so homochiral count increases by 1
            homochiralcount +=1
            #brk gets benefited by the poisson distribution of homochiral count 
            brk *= poissonequation(homochiralcount)
    
            #set j to monomer to right
            j = numbermonomer+1
    
            #going right to check for neighboring homochirality
            while(self.easthomochiralcheck(j)):
                #if homochiral increases count
                homochiralcount += 1
                #recalculates brk prob
                brk *= poissonequation(homochiralcount)
    
                #check next monomer
                j += 1
    
        #going left to check for neighboring homochirality
        j=numbermonomer-1
        homochiralcount = 1
        while(self.easthomochiralcheck(j)):
            homochiralcount += 1
            brk *= poissonequation(homochiralcount)
            j -= 1
    
        return brk
    
    
    def reset_break_probability(self, method):
        
        """ 
        resets the break probabilities of the monomers in a 
        sequence. 1)checks lookup table to see if similar polymer values
        have been caluclated, otherwise calculates values and adds
        polymer to lookup table 2)sets values to monomers.
    
        """
        #store self as list of booleans to be added to lookup table
        polyseq = self.get_monomers().copy()
        for n in range(self.get_length()):
            polyseq[n] = polyseq[n].get_handedness()
    
        ##########################
        # LEFT/RIGHT SENSITIVITY #
        ##########################
        # All stored sequences start with True. If the one being requested doesn't, invert so that it does
        if not method == "bias":
            if polyseq[0] == True:
                  polyseq = [not mono for mono in polyseq]
    
        #make polyseq a tuple able to be added as an index of a dictionary  
        polyseq = tuple(polyseq)
    
        # If the sequence is already in the table (dict), retrieves its previously calculated break probabilities
        if polyseq in break_prob_lookup_table:
            the_probs = break_prob_lookup_table[polyseq]
    
        else: # Otherwise calculate and add to the lookup table
            # Calculate the break probabilities
            the_probs = self.calculatebrkprob(method)
            #and add it to the lookup table
            break_prob_lookup_table[polyseq] = the_probs
    
        #set monomers eastbond brk probs to correlated values stored in the_probs
        for n in range (self.get_length()):
            self._monomers[n].set_eastbrkprob(the_probs[n])
    
        #end def: polymer break probabilities recalculated and added
        #to lookup table (if not already there)
    
    def calculatebrkprob(self, method):
    
        """
        generates and returns a list of break probabilities for the monomers at each index of
        a polymer. This method can be altered to calculate normally
        with left/right sensitivity or with a poisson distribution
    
        Calculation methods: [bias, standard, poisson]
    
        """
    
        #initialize brk_probs as an empty list
        brk_probs = []
    
        #checks if polymer is of length one -meaning it is lone monomer
        if self.get_length() == 1 :
            #if so adds -1 to list(brk prob value assigned to individual monomers)
            brk_probs.append(-1)
        
        else:
            #otherwise goes through all the indices of the polymer,calculates its brk probability, and appends it to brk_probs 
            for n in range(self.get_length()):
        
                #############################################################
                #LEFT RIGHT SENSITIVITY/POISSON EQUATIONS CHOICE.           #
                #############################################################
        
                if method == "bias":
                    brk_probs.append(self.biaseastbondbreakprobability(n))
        
                elif method == "standard":
                    brk_probs.append(self.eastbondbreakprobability(n))
        
                elif method == "poisson":
                    brk_probs.append(self.poissonbreakprobability(n))
        
                else:
                    raise ValueError(f"'{method}' is not a recognized calculation method. Please use one of ['bias', 'standard', 'poisson'].")
    
    
        #returns in order list of brk probs for the east bond of all the monomers in a polymer
        return brk_probs


#%%### REACTABLES ###
        
class _Reactables:

    def __init__(self, reactables = []):
        """
        Constructor class for Reactables object
        Accepts list of reactables or defaults to empty list
        Reactables functions as a bin for reactable objects monomer and polymer
        Methods consist mostly of list handling
        """
        self._reactables = reactables
        #The reactables bag also handles remembering its history
        self._history = []
        self._lookup = {}
        self._hist_stats = pd.DataFrame()


    def __str__(self):
        """
        str method for Reactables class
        """
        return str(self.get_reactables())
    
    def get_reactables(self):
        """
        Getter method for reactables in Reactables
        Returns a list of the objects currently in the Reactables bag
        """
        return self._reactables
    
    def get_history(self):
        """ Return the history of the reactables bag.
        """
        return self._history
    
    def set_history(self, new_hist):
        """ Set the history to new_hist. """
        self._history = new_hist
        
    def get_count(self):
        """
        Getter method for Reactables class
        Returns the number of objects in the Reactables bag
        """
        return len(self.get_reactables()) 
    
    
    def set_reactables(self, new_list):
        """
        Setter method for the Reactables class
        Sets the reactables bag to a new list
        Used mostly as a helper function in bonding
        """
        self._reactables = new_list
    
    
    def get_stats(self):
        return self._hist_stats
    
    def set_hist_stats(self,df):
        self._hist_stats = df
    
    def get_lookup(self):
        return self._lookup
    
    def add_stat_lookup(self, key, value):
        self._lookup[key] = value
    
    def add_log(self, log):
        (self._history).append(log.copy())
    


    def refill(self, pool_size):
    
        """Method that refills the pool size after each iteration"""
        
        #initialize a number to record number of monomers to be added to the bag
        number = 0
        #calculate the difference between the pool size and the actual size
        dearth = pool_size - self.get_count()
        #if reactant size is less than pool size
        if dearth  > 0:
        
            #if refill random method is selected
            if REFILL[0] == "refill random":
                #number = random percentage of the dearth
                number = int(get_rand() * dearth) 
        
            #if refill percent is selected
            elif REFILL[0] == "refill percent":
                #refill by selected percent of dearth 
                number = int(REFILL[1] * dearth)
        
            #if refill number is selected
            elif REFILL[0] == "refill number":
                #refill selected number of monomers
                number = REFILL[1]
        
            #if refill number decrease
            elif REFILL[0] == "refill number decrease":
                #refill recorded number
                number = REFILL[1]
                #if number is >0
                if REFILL[1] > 0:
                    #decrease by 1
                    REFILL[1] -= 1
            #if refill normal
            elif REFILL[0] == "refill normal":
                #refill dearth
                number = dearth
            #if refill percent decrease
            elif REFILL[0] == "refill percent decrease":
                #refill recorded percent of difference
                number = int(REFILL[1] * dearth)
                #if percent is greater than 0
                if REFILL[1]> 0:
                    #decrease by .01
                    REFILL[1] -= .01
                    
            else:
                raise ValueError("Bad refill method.")
            
            #loop through number of monomers to be added
            for n in range(number):
                #add a monomer to the reacables bag
                self.add([_Monomer()])
    
    def randomize_reactables(self):
        """
        Randomizes the order of the reactables list
        """
        #Get the current reactables
        reactables = self.get_reactables()
        #Shuffle them
        random.shuffle(reactables)
        #And reset the reactables list to the new sorted list
        self.set_reactables(reactables)
    
    def do_the_thing(self, method):
        """
        Handle a single iteration of the reactables.
    
        """
        #We keep track of bonding by storing the most recent reactable to choose to bond as the 'bachelor.'
        bachelor = None
        #We need a copy of the reactables because indices will change as soon as we start doing stuff
        reactables = copy.copy(self.get_reactables())
        #Iterate through each reactable
        for item in reactables:
            #Roll a random number in (0,1) to compare against our probabilities
            roll = get_rand()
            #If our reactable is a monomer...
            if isinstance(item, _Monomer):
                #We get our breaking and bonding chances (parameters now, may become functions later?)
                brk = POOF_CHANCE
                bond = BOND_PROB
                #If we roll to break the monomer
                if 0 < roll <= brk:
                    #Just delete it
                    del item
                #If we roll to bond the monomer
                elif brk < roll <= brk+bond:
                    #Check if there's a reactable waiting to bond
                    if bachelor == None:
                        #If not, make this monomer the bachelor. Someone will come along to bond later.
                        bachelor = item
                    else:
                        #If there is someone waiting, bond them together and reset the bachelor chair to empty.
                        self.bond_pair(bachelor, item, method)
                        bachelor= None
                else:
                    continue
            #If the reactable is a polymer instead...
            elif isinstance(item, _Polymer):
                #Choose whether the polymer will break or bond. It's 50-50 right now.
                if roll >= 0.33:
                    #This is the same bonding logic as for the monomer. It's not in a helper function because it needs to be able to access bachelor.
                    if bachelor == None:
                        bachelor = item
                    elif (type(bachelor) != _Polymer or type(item) != _Polymer):
                        self.bond_pair(bachelor, item, method)
                        bachelor = None
                    elif FUSION:
                        self.bond_pair(bachelor,item, method)
                        bachelor = None
                #Here's what happens if the polymer is chosen to check breaking
                else:
                    #It finds its break location
                    break_spot = item.brkloc()
                    #Which might turn out to be nowhere.
                    if break_spot != None:
                        #If somewhere does break, call break_polymer to handle the breaking
                        self.break_polymer(item,break_spot,method)
            else:
              #You never know what might end up in your reactables bag
              raise ValueError("This thing is neither a monomer or a polymer. What?!")
    
    def iterate(self,size,iteration_number,method):
        """
        Handles full sequence of iteration
        """
        self.randomize_reactables()
        self.do_the_thing(method)
        self.log(iteration_number)
        self.ageup()
        self.refill(size)
    
    def simulate(self, poolsize, iterations, method="standard"):
        for n in range(iterations):
            self.iterate(poolsize, n, method)
        self.parse_history()

# Functionality helpers #
        
    def add(self, new_reactables):
        """
        Adds a list of reactables to the reactables bag
        Built on extend(); argument must be iterable
        """
    
        self.get_reactables().extend(new_reactables)
        #print("reactables list is now",self.get_reactables())
    
    
    def subtract(self, removables):
        """
        Removes a list of reactables from the reactables bag
        Built on list comprehension; argument must be iterable
        """
        self.set_reactables([reactable for reactable in self.get_reactables() if not (reactable in removables)])
    
    
    def bond_pair(self, west, east, method):
        """
        Accepts two reactables from the bag
        Bonds them together, handling deletion of any emptied polymer
        """
    
        #this is my attempt to get the ages working properly the idea is:
    
        #if both are polymers or both are monomers:
        if (type(west) == type(east) and type(west) == _Polymer):
            #set the new age to the maximum age of the polymers
            age = max(west.get_age(), east.get_age())
      
        if (type(west) == type(east) and type(west) == _Monomer):
            #set age to 0 it is a brand new polymer
            age = 0
    
        #if only one is a polymer:
        elif isinstance(west, _Polymer) and isinstance(east, _Monomer):
            #set the age to the age of the polymer
            age = west.get_age()
        elif isinstance(west, _Monomer) and isinstance(east, _Polymer):
            # set the age to the age of the polymer
            age = east.get_age()
    
        #Take your two reactables and check if the west one is a monomer or a polymer 
        if isinstance(west, _Polymer):
            #If it's a polymer, all you've got to do is append the east item and then take it out of the bag
            west.append(east)
            self.subtract([east])
            #recalculate the breakprobabilities of the bonds in polymer
            west.reset_break_probability(method)
            #set new polymer to age that was calculated before
            west.set_age(age)
    
        if isinstance(west, _Monomer):
            #If west friend is a Monomer, then make a new polymer containing west
            newpoly = _Polymer([west])
            #Add the east things to it
            newpoly.append(east)
            #Remove both the west monomer and the east thing
            self.subtract([east])
            self.subtract([west])
            #And put the new polymer into the reactables bag
            self.add([newpoly])
            #recalculates break probability
            newpoly.reset_break_probability(method)
            #set new polymer to age that was calculated before
            newpoly.set_age(age)
    
    
    
    def break_polymer(self, polymer, brk_location, method):
        """
        breaks polymer at given location and creates a new polymer of the 
        monomers removed
        """
        age = polymer.get_age()
        #when polymer is made-breakprobabilities are calculated
        newPolymer = _Polymer(polymer.removeright(brk_location))
        polys = [newPolymer, polymer]
        #goes throught the two new polymers
        for poly in polys:
            #resets break probabilities in polymers
            poly.reset_break_probability(method)
            #sets age to that of the oldest polymer
            poly.set_age(age)
    
            #checks if they are of length 1(monomer)
            if poly.get_length() <= 1:
                #adds them to reactable as a single monomer
                self.add(poly.get_monomers())
                #subtracts polymer from reactable list
                self.subtract([poly])
                del poly
            #checks if the polymer is not in the reactables list
            elif (poly not in self.get_reactables()):
                #add polymer to reactable list
                self.add([poly])
    
    def ageup(self):
        """
        Method that ages up every reactable in the reactables bag
        """
        for reactable in self.get_reactables():
            reactable.birthday()
    
        
# History handling #
            
    def log(self, iteration):
        """ 
        Translate current reactables into a loggable list of tuples.
        This log keeps the full sequences of the polymers without digesting the data, 
        so it should be more flexible as we more forward with pulling new information
        out of this simulation. 
        There is one major piece of information lost, however, which is the identities
        of the individual polymers and monomers. 
        """
    
        # Column labels for History Handling #
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
        Iter = "Iteration"              # The iteration number at which the item is found
        Age = "Age"                     # The age of the item
        Id = "ID"                       # The ID of the item
        C_lens = "Chain lengths"        # A list of the chain lengths in the item
        Max_C = "Longest chain length"  # The length of the longest chain in the item
        Contains = "Contains"           # The monomers that the polymer contains
        
    
        for r in self.get_reactables():
            
            item = r.get_sequence()
            age = r.get_age()
            # Translate
            sequence = self.standard_form(item)
    
            #Is it in the lookup table?
            if sequence not in self.get_lookup():
    
                # If not, check what kind of data to gether
                if len(item) == 1: # It's a monomer
                    # Get the info
                    lr = self.count_LR(item)
                    # And put it in the lookup table
                    self.add_stat_lookup(sequence, {Type: 'Monomer', 
                                                    Length:1, 
                                                    nLefts:lr[0], 
                                                    nRights:lr[1],
                                                    Sequence:sequence})
                elif len(item) > 1: # It's a polymer
                    # Get basic info
                    length = len(item)
                    lr = self.count_LR(item)
                    bonds = self.homochiral_bond_counts(item)
                    total_homos = bonds[0]+bonds[1]
                    signed_ee = (lr[0]-lr[1])/(lr[0]+lr[1])
                    homochirality = total_homos/(length-1)
                    lhomochirality = bonds[0]/(length-1)
                    rhomochirality = bonds[1]/(length-1)
                    chains = self.chain_lengths(item)
                    max_c = max(chains)
                    # And put it in the lookup table
                    self.add_stat_lookup(sequence, {Type: 'Polymer', 
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
                                                    Sequence:sequence,
                                                    C_lens:chains,
                                                    Max_C:max_c})
                else:
                    raise ValueError("There's something with length 0 in your history.")
        
            # Now that the data is searchable...
            # Log it
            new_log = self.get_lookup().get(sequence)
            new_log[Iter] = iteration
            new_log[Age] = age
            new_log[Id] = r.get_id()
            if len(item) > 1:
                new_log[Contains] = self.get_contents(r)
            self.add_log(new_log)
    
    
    def parse_history(self):
        """ Create an array of plottable information from the history log. """
        parsed_hist = pd.DataFrame(self.get_history())
        self.set_hist_stats(parsed_hist)
    


# Low level info functions #
        
    # Functions for getting the info we want to log
        
    
    def get_contents(self,polymer):
        monomers = polymer.get_monomers()
        return [m.get_id() for m in monomers]
    
    
    # Number of right and left monomers in a polymer (works on monomers too)
    
    
    def count_LR(self, log):
        """ Return (n_Lefts, n_Rights) """
        if L and not R:
            return (sum(log), len(log)-sum(log))
        elif R and not L:
            return (len(log)-sum(log), sum(log))
        else:
            raise ValueError("Your L's and R's are screwed up somehow.")
        
    
    
    # HOMOCHIRALITY SEQUENCES #
    
    def homochirality_sequence(self, log):
        """ Return a boolean list of bonds within a given logged polymer, True if the bond is
          homochiral. """
        return [i[0]==i[1] for i in list(more_itertools.pairwise(log))]
    
    
    # HOMOCHIRAL BOND COUNTS #
        
    
    def homochiral_bond_counts(self,log):
        """ Return (number of left homochiral bonds, number of right homochiral bonds) """
        homo = self.homochirality_sequence(log)
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
    
    
    def chain_lengths(self,polylog):
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
        
    
    def standard_form(self,poly):
        if len(poly) == 0:
            return ''
        return ''.join(['L' if m else 'R' for m in poly])
        
    
#%%### SIMULATION ###
        
    
class Simulation:
    def __init__(self,params):
        self._parameters = params
        self.set_globals(self.get_params())
        self._results = None
        # Make reactables
        self._react_bag = _Reactables(make_pool(POOL_SIZE))
        self._timestamp = str(datetime.datetime.now()).replace(" ","_").replace(":",".")
        
        
    def set_globals(self,params):
        # Set the vars
        global BASE_BOND_BREAK_PROBABILITY
        global HOMOCHIRAL_BREAK_FACTOR
        global HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR
        global LENGTH_FACTOR
        global N
        global LAMBDA
        global HOMOCHIRAL_BREAK_FACTOR_LEFT
        global HOMOCHIRAL_BREAK_FACTOR_RIGHT
        global HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR_LEFT
        global HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR_RIGHT
        global POOF_CHANCE
        global BOND_PROB
        global POOL_SIZE
        global ITERATIONS
        global METHOD
        global POISSON_FACTOR
        global FUSION
        global REFILL
        
        BASE_BOND_BREAK_PROBABILITY = params["BASE_BOND_BREAK_PROBABILITY"]
        HOMOCHIRAL_BREAK_FACTOR = params["HOMOCHIRAL_BREAK_FACTOR"]
        HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR = params["HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR"]
        LENGTH_FACTOR = params["LENGTH_FACTOR"]
        N = params["N"]
        LAMBDA = params["LAMBDA"]
        HOMOCHIRAL_BREAK_FACTOR_LEFT = params["HOMOCHIRAL_BREAK_FACTOR_LEFT"]
        HOMOCHIRAL_BREAK_FACTOR_RIGHT = params["HOMOCHIRAL_BREAK_FACTOR_RIGHT"]
        HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR_LEFT = params["HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR_LEFT"]
        HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR_RIGHT = params["HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR_RIGHT"]
        POOF_CHANCE = params["POOF_CHANCE"]
        BOND_PROB = params["BOND_PROB"]
        POOL_SIZE = params["POOL_SIZE"]
        ITERATIONS = params["ITERATIONS"]
        METHOD = params["METHOD"]
        POISSON_FACTOR = params["POISSON_FACTOR"]
        FUSION = params["FUSION"]
        REFILL = params["REFILL"]
        

    def get_params(self):
        return self._parameters
    
    def _set_params(self,new_params):
        self._parameters = new_params
        
    def _set_results(self,results):
        self._results = results
        
    def get_results(self):
        return self._results
    
    def get_timestamp(self):
        return self._timestamp
        
        
    
    def run(self):
        # Set parameters
        self.set_globals(self.get_params())
        
        # Init lookup tables
        global L
        global R
        global N_RANDS
        global rands
        global break_prob_lookup_table
        global poisson_dict
        
        rands = []
        L = True
        R = False
        N_RANDS = 1000
        break_prob_lookup_table = {}
        poisson_dict = {}
       
        # Run simulation
        self._react_bag.simulate(POOL_SIZE, ITERATIONS, METHOD)
        # Store results
        self._set_results(self._react_bag.get_stats())

    def export(self, path):
        if not isinstance(self.get_results(), pd.DataFrame):
            raise RuntimeError("Can't export sim before running")
       
        # Export parameter file
        dict_to_csv(self.get_params(),f"{path}_{self.get_timestamp()}_params.csv")
        # Export data file
        self.get_results().to_csv(f"{path}_{self.get_timestamp()}_data.csv")
        
    
            
