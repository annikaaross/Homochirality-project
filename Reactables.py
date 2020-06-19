#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 09:52:40 2020

@author: anni
"""
import Polymer
import Monomer
import globalvariables
import random
import numpy as np
import copy
#from google.colab import widgets
import matplotlib.pyplot as plt
import jdc
import more_itertools

#This is the Reactables class aka the stylish little bag we keep our contestants 
#in. As all good dating shows, we define our members by the relationship they are in. 
#So a single monomer will be addd to the bag alone, but a monomer in a polymer will not be.
# Rather each polymer will be recorded once. Its like if you have charlie, ben, and jack. 
#Jack and Ben are in a relationship so they are jointly called Jen. If these three 
#were in the reactable bag they would be recognized as two things, charlie and 
#Jen(Ben+JAck). The reactable class is the heart and center of our game- 
#I mean dating- show. Its here we ask all the members just the right 
#questions to get them falling in love and breaking hearts, and then we are 
#curteous enough to do all the splitting and merging for them.

#Constructor
class Reactables:

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
        self._leftrighthistory = []


    def __str__(self):
        """
        str method for Reactables class
        """
        return str(self.get_reactables())
    
#Getters and Setters
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

    def get_leftrighthistory(self):
        return self._leftrighthistory

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


    def get_bond_chance(self, reactable1, reactable2):
        """
        Returns the probability that two given reactables in the reactables bag will bond
        """
        return globalvariables.BOND_PROB

    def get_free_monomers(self):
        return [r for r in self.get_reactables() if isinstance(r, Monomer)]

    def get_polymers(self):
        return [r for r in self.get_reactables() if isinstance(r, Polymer)]

    def get_bound_monomers(self):
        bound = []
        for p in self.get_polymers():
            bound.extend(p.get_monomers())
        return bound

    def record_leftrighthomochiral(self):
        """returns the number of left and right homochiral bonds, (left,right)"""
        l=0
        r=0
        for polymer in self.get_polymers() :
            tuple1=polymer.get_leftrighthomochiralbonds()
            l += tuple1[0]
            r += tuple1[1]
        tuple2=(l,r)
        self._leftrighthistory.append(tuple2)
        return tuple2
    
    #Functionality Methods####
    #Top Level Functionality
    def refill(self, pool_size):
        #Do we have less than the required number of things in the pool?
        dearth = pool_size - self.get_count()
        #If so,
        if dearth > 0:
            new = []
            #Make as many monomers as we need to get back up to the riquired amount
            for n in range(dearth):
                new.append(Monomer())
            #And add them to the reactables bag
            self.add(new)

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

    def do_the_thing(self):
        """Handle a single iteration of the reactables.

        """
        #We keep track of bonding by storing the most recent reactable to choose to bond as the 'bachelor.'
        bachelor = None
        #We need a copy of the reactables because indices will change as soon as we start doing stuff
        reactables = copy.copy(self.get_reactables())
        #Iterate through each reactable
        for item in reactables:
            #Roll a random number in (0,1) to compare against our probabilities
            roll = random.random()
            #If our reactable is a monomer...
            if isinstance(item, Monomer):
            #We get our breaking and bonding chances (parameters now, may become functions later?)
                brk = globalvariables.POOF_CHANCE
                bond = globalvariables.BOND_PROB
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
                        self.bond_pair(bachelor, item)
                        bachelor = None
                else:
                    continue
                #If the reactable is a polymer instead...
            elif isinstance(item, Polymer):
                #Choose whether the polymer will break or bond. It's 50-50 right now.
                if roll >= 0.5:
                    #This is the same bonding logic as for the monomer. It's not in a helper function because it needs to be able to access bachelor.
                    if bachelor == None:
                        bachelor = item
                    else:
                        self.bond_pair(bachelor, item)
                        bachelor = None
                #Here's what happens if the polymer is chosen to check breaking
                else:
                    #It finds its break location
                    break_spot = item.brkloc()
                    #Which might turn out to be nowhere.
                    if break_spot != None:
                    #If somewhere does break, call break_polymer to handle the breaking
                        self.break_polymer(item,break_spot)
            else:
            #You never know what might end up in your reactables bag
                raise ValueError("This thing is neither a monomer or a polymer. What?!")

    def iterate(self,size):
        """
        Handles full sequence of iteration
        """
        self.randomize_reactables()
        self.do_the_thing()
        self.log()
        self.record_leftrighthomochiral()
        self.refill(size)
        
    #functionality helpers
    def add(self, new_reactables):
        """
        Adds a list of reactables to the reactables bag
        Built on extend(); argument must be iterable
        """
        #print("list of new Reactbles", new_reactables)
        #print("list of get Reactables:",self.get_reactables())
        self.get_reactables().extend(new_reactables)
        #print("reactables list is now",self.get_reactables())


    def subtract(self, removables):
        """
        Removes a list of reactables from the reactables bag
        Built on list comprehension; argument must be iterable
        """
        self.set_reactables([reactable for reactable in self.get_reactables() if not (reactable in removables)])


    def bond_pair(self, west, east):
        """
        Accepts two reactables from the bag
        Bonds them together, handling deletion of any emptied polymer
        """
        #Take your two reactables and check if the west one is a monomer or a polymer 
        if isinstance(west, Polymer):
            #If it's a polymer, all you've got to do is append the east item and then take it out of the bag
            west.append(east)
            self.subtract([east])
            #recalculate the breakprobabilities of the bonds in polymer
            west.reset_break_probability()
        if isinstance(west, Monomer):
            #If west friend is a Monomer, then make a new polymer containing west
            newpoly = Polymer([west])
            #Add the east things to it
            newpoly.append(east)
            #Remove both the west monomer and the east thing
            self.subtract([east])
            self.subtract([west])
            #And put the new polymer into the reactables bag
            self.add([newpoly])
            #recalculates break probability
            newpoly.reset_break_probability()


    def break_polymer(self, polymer, brk_location):
        #breaks polymer at given location and creates a new polymer of the 
        #monomers removed
        #when polymer is made-breakprobabilities are calculated
        newPolymer = Polymer(polymer.removeright(brk_location))
        polys = [newPolymer, polymer]
        #resets break probabilities in polymer
        polymer.reset_break_probability()
        #goes throught the two new polymers
        for poly in polys:
            #checks if they are of length 1
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
    #Current State Data Methods
    def get_overall_chirality(self):
        """ Return the total quantity of left- and right-handed monomers in the reactable
        pool, including monomers within polymers.

        RETURN tuple
        (Left-handed_quantity, Right-handed_quantity)
        """
        left_count = 0
        right_count = 0
        for reactable in self.get_reactables():
            if isinstance(reactable, Monomer):
                if reactable.get_handedness():
                    left_count += 1
                else:
                    right_count += 1
            elif isinstance(reactable, Polymer):
                l, r = reactable.get_chirality()
                left_count += l
                right_count += r
            else:
                print(f"{reactable}, which is a {type(reactable)}")
                raise ValueError
        return (left_count, right_count)

    def get_polymer_chirality(self):
        """ Return the total quantity of left- and right-handed monomers in the polymers
        within the reactable pool.

        RETURN tuple
        (Left-handed_quantity, Right-handed_quantity)
        """
        left_count = 0
        right_count = 0
        for reactable in self.get_reactables():
            if isinstance(reactable, Polymer):
                l, r = reactable.get_chirality()
                left_count += l
                right_count += r
        return (left_count, right_count)

    def get_free_proportion(self):
        """Return the proportion of free monomers to total population"""
        return len(self.get_free_monomers())/self.unit_count()

    def unit_count(self):
        """return the total number of monomers in the bag, both free and bound"""
        bound = 0
        for p in self.get_polymers():
            bound += p.get_length()
        return len(self.get_free_monomers()) + bound

    def get_polymer_ee(self):
        ee = []
        for polymer in self.get_polymers():
            ee.append(polymer.get_signed_ee())
        return ee

    def avg_ee_by_length(self):
        """ Return a list containing the enantiomeric excess of all reactables sorted by size.
        That is:
        [(ee of monomers),(ee of 2-mers),(ee of 3-mers) etc]
        """
        output = []
        for n in range(1,self.max_length()+1):
            group = self.get_reactables_by_length(n)
            ees = []
            # print(f"For length {n} I got {len(group)} reactables.")
            if len(group) > 0:
                # eestr = ""
                for r in group:
                    ees.append(r.get_ee())
                    # eestr += f", {r.get_ee()}"
                output.append(np.mean(ees))
                # print(f"The ee's of this length were {eestr}.\nThe average ee was {np.mean(ees)}.")
            else:
                output.append(0)
        return output

    def get_reactables_by_length(self,length):
        if length <= 1:
            #return the monomers
            return self.get_free_monomers()
        else:
            return [p for p in self.get_polymers() if p.get_length() == length]

    def max_length(self):
        """ Return the length of the longest polymer in the reactables bag.
        """
        length = 0
        for p in self.get_polymers():
            if p.get_length() > length:
                length = p.get_length()
        return length
  #History Handling Methods
  # def old_log(self):
#   """ Translate current reactables list into a lighter-weight list of tuples 
#     describing the bag state, then logs that list into self._history.
#     Log format: (length,signed_ee)
#     e.g.  True-monomer would be (1,1)
#           [False, True, False, False]-polymer would be (4,-0.5)
#     """
#   log = []
#   for r in self.get_reactables():
#     log.append(r.generate_old_log())
#   (self._history).append(log)

    def log(self):
        """ Translate current reactables into a loggable list of tuples.
        This log keeps the full sequences of the polymers without digesting the data, 
        so it should be more flexible as we more forward with pulling new information
        out of this simulation. 
        There is one major piece of information lost, however, which is the identities
        of the individual polymers and monomers. """
        log = []
        for r in self.get_reactables():
            log.append(r.generate_log())
        (self._history).append(log)




# def get_iter_polymers_by_length(self, single_iteration_log):
#   # So we have something like this: [(t,f,f),(t),(f,f),(f),(t,t,t,f)]
#   polymers_by_length = []
#   lengths = [len(log) for log in single_iteration_log]
#   biggest = max(lengths)
#   for n in range(1, biggest+1):
#     polys = [log for log in single_iteration_log if len(log) == n]
#     if len(polys) > 0:
#       polymers_by_length.append(polys)
#     else:
#       polymers_by_length.append([])
#   return polymers_by_length

# def get_iter_ees_by_length(self,single_iter_log):
#   polymers_by_length = self.get_iter_polymers_by_length(single_iter_log)
#   for length_category in polymers_by_length:
#     for polymer in length_category:
#       ee = hist_get_polymer_ee(polymer)

#Low Level Info Functions
    def hist_get_polymer_homochirality_of_bonds(self,polylog):
        """ Return a boolean list of bonds within a given logged polymer, True if the bond is
        homochiral. """
        return [i[0]==i[1] for i in list(more_itertools.pairwise(polylog))]

    def hist_get_polymer_ee(self,polylog):
        """ Return the (signed) enantiomeric excess of a logged polymer. """
        length = len(polylog)
        n_True = sum(polylog)
        n_False = length - n_True
        return (n_True - n_False) / length

    def hist_get_polymers(self,iteration):
        """ Return a list of the polymer logs in an iteration. This is just the iteration
        but without the monomers. """
        return [r for r in iteration if len(r) > 1]

    def hist_count_longest_homochiral_chain(self,polylog):
        """ Return the length of the longest homochiral chain given the log of a polymer. """
        previous = None
        count = 1
        longest = 1
        for monomer in polylog:
            if monomer == previous:
                count += 1
            else:
                longest = max(count, longest)
                count = 1
            previous = monomer
        longest = max(count, longest)
        return longest

    def hist_get_polymer_chain_lengths(self,polylog):
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
    
    def hist_get_iteration_chain_lengths(self,iteration):
        """ Return the lengths of all the homochiral chains in polymers in a given iteration. """
        polymers = self.hist_get_polymers(iteration)
        chain_lengths = []
        for polymer in polymers:
            chain_lengths.extend(self.hist_get_polymer_chain_lengths(polymer))
        return chain_lengths

#Plotting Methods
    def visualize_parameters(self,stamp=""):
        fig,ax = plt.subplots()

        parameters = ['POOF_CHANCE','BOND_PROB','BASE_BOND_BREAK_PROBABILITY',
                'HOMOCHIRAL_BREAK_FACTOR','HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR',
                'LENGTH_FACTOR','HOMOCHIRAL_BREAK_FACTOR_LEFT',
                'HOMOCHIRAL_BREAK_FACTOR_RIGHT',
                'HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR_LEFT',
                'HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR_RIGHT']
        values = [globalvariables.POOF_CHANCE,globalvariables.BOND_PROB,globalvariables.BASE_BOND_BREAK_PROBABILITY,
            globalvariables.HOMOCHIRAL_BREAK_FACTOR,globalvariables.HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR,
            globalvariables.LENGTH_FACTOR,globalvariables.HOMOCHIRAL_BREAK_FACTOR_LEFT,
            globalvariables.HOMOCHIRAL_BREAK_FACTOR_RIGHT,
            globalvariables.HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR_LEFT,
            globalvariables.HOMOCHIRAL_NEIGHBOR_IMPROV_FACTOR_RIGHT]

        ypos = np.arange(len(parameters))

        ax.barh(ypos, values, align='center')
        ax.set_yticks(ypos)
        ax.set_yticklabels(parameters)
        ax.invert_yaxis()  # labels read top-to-bottom
        ax.set_xlabel('Value')
        ax.set_title(f"{stamp}Parameter values")



    def plot_signed_ee_spread(self,cmap='bone',stamp=""):
        #First make the polymerspread thingy
        poly_ee_hist = []
        for iteration in self.get_history():
            polys = [log for log in iteration if not isinstance(log,bool)]
            polydata = [log for log in polys if len(log) > 1] # All the polymers (not monomers) in the iteration
            poly_ees = [self.hist_get_polymer_ee(poly) for poly in polydata] # the EE of each polymer in the iteration
            poly_ee_hist.append(poly_ees) # Put the EE list in our plottable structure
        #Now we can make the plot
        fig,ax = plt.subplots()
        ys = []
        xs = []
        for n in range(len(poly_ee_hist)):
            y = poly_ee_hist[n]
            x = [n] * len(y)
            ys.extend(y)
            xs.extend(x)
        ax.hexbin(xs,ys,cmap=cmap)
        ax.set_title(f"{stamp}polymer spread")
        ax.set_xlabel("iteraions")
        ax.set_ylabel("ee")

    def plot_ee_heatmap(self,cmap='bone',stamp=""):
        #Need to figure out how to get ees_by_length from history
        ees_by_length = []
        for iteration in self.get_history():
            ee_data = []
            #get the polymers
            polys = [log for log in iteration if not isinstance(log,bool)]
            #get the length of the longest polymer
            lengths = [len(log) for log in polys]
            try:
                biggest = max(lengths)
            except:
                biggest = 0
            for n in range(1,biggest+1):
                ees = [abs(self.hist_get_polymer_ee(log)) for log in polys if len(log) == n]
                if len(ees) > 0:
                    ee_data.append(np.mean(ees))
                else:
                    ee_data.append(0)
        ees_by_length.append(ee_data)
        #Need to make a new ee array with rectangular dimensions
        ee = numpy_fillna(ees_by_length)
        fig,ax = plt.subplots(dpi=150)
        im = ax.imshow(ee,cmap=cmap)
        ax.set_title(f"{stamp}ee of polymers by length through iterations")
        ax.set_xlabel("polymer lengths")
        ax.set_ylabel("iterations")
        fig.tight_layout()
        plt.show()

    def plot_leftrighthomochirality(self):
        plt.figure()
        l=[]
        r=[]
        x=[]
        i=1
        for iteration in self.get_leftrighthistory():
            l.append(iteration[0])
            r.append(iteration[1])
            x.append(i)
            i+= 1
        plt.plot(x,l,label="left homochiral")
        plt.plot(x,r,label="right homochiral")
        plt.ylabel("Number of bonds")
        plt.xlabel("Iterations")
        plt.title("Left right homochiral bonds when LEFT="+str(globalvariables.HOMOCHIRAL_BREAK_FACTOR_LEFT) +  "and right=" + str(globalvariables.HOMOCHIRAL_BREAK_FACTOR_RIGHT))
        plt.legend()
        plt.show()

    def plot_homochiral_chain_lengths(self):
        """ Plot a chart of the mean length of homochiral chains present at each iteration.
        The plot generated will attempt to summarize the mean lengths of homochiral
        chains within polymers across iterations. It will also attempt to visualize the 
        varience with continuous lines plotting the standard deviation around that mean.
        
        The x axis of the resulting plot is iterations, and the y axis is length of 
        the homochiral chains. Both of these are discrete categories, but the plot will
        treat iterations as continuous to improve readability for large n.
        """
       
        history = self.get_history()
        means = []
        stdevs= []
        maxes = []
        for iteration in history:
            # Get the chain lengths in that iteration
            iter_data = self.hist_get_iteration_chain_lengths(iteration)
            means.append(np.mean(iter_data))
            stdevs.append(np.std(iter_data))
            maxes.append(max(iter_data))
        # Plot those data
        fig,ax = plt.subplots()
        ax.fill_between(np.arange(0,len(means)), [m + s for m, s in zip(means, stdevs)],[m - s for m, s in zip(means, stdevs)], alpha=0.2, label = "Means +- one standard deviation")
        ax.plot(means, 'k-', label = "Mean homochiral chain length")
        ax.plot(maxes, 'b.', label = "Max homochiral chain length")
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Homochiral chain length")
        ax.legend()

    def plot_homochirality_vs_length(self):
        """ Take all the polymers that ever existed in history and plot them on a
        scatter plot of length vs homochirality. Homochirality is the percent of their
        bonds that are homochiral. """
        #Get the history
        history = self.get_history()
        #Get the polymers
        x = []
        y = []
        for i in history:
            for p in self.hist_get_polymers(i):
                homochirality = sum(self.hist_get_polymer_homochirality_of_bonds(p))/len(p)
                x.append(len(p))
                y.append(homochirality)
        fig,ax = plt.subplots()
        ax.scatter(x,y)
        ax.set_xlabel("Polymer lengths")
        ax.set_ylabel("Homochirality")
  




