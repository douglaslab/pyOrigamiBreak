#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2017-08-10 01:32:19
# @Author  : Tural Aksel (turalaksel@gmail.com)
# @Link    : http://example.org
# @Version : $Id$

import numpy as np
import cadnano
import random
import os
import sys
import argparse
import utilities
import math
import copy

from cadnano.document import Document
from functools import reduce


class Sequence:
    def __init__(self):
        '''
        DNA sequence class
        '''
        self.dna           = None
        self.next_sequence = None

class OligoBreakSolution:
    def __init__(self):
        '''
        Break solution for an oligo
        '''
        
        self.start_break   = None
        self.final_break   = None
        self.oligo         = None
        self.break_paths   = None
        self.score         = None
        self.bad_list      = None
        self.breaks        = None
        self.scores        = None
        self.edges         = None

    def initialize(self):
        '''
        Prepare the lists
        '''
        self.breaks = [break_path.break_node for break_path in self.break_paths[::-1]]
        self.edges  = [break_path.break_edge for break_path in self.break_paths[::-1]]
        self.scores = [break_path.score      for break_path in self.break_paths[::-1]]

    def is_identical(self,other_solution,max_index=None):
        '''
        Compare the baths between two solutions
        '''
        #Determine maximum index
        max_i = len(self.breaks)-1
        if max_index:
            max_i = max_index

        #Identical variable
        identical = True

        #Breaks for the current path
        current_breaks = self.breaks[:max_i]

        #Other breaks
        other_breaks   = other_solution.breaks[:max_i]

        if not len(current_breaks) == len(other_breaks):
            return False

        #Pairwise comparison of elements
        for i in range(max_i):
            identical *= (current_breaks[i] == other_breaks[i])

        return identical
            

class GroupBreaksolution:
    def __init__(self):
        '''
        Break solutions for oligo group
        '''

        self.break_solutions = None
        self.total_score     = 0
        self.total_penalty   = 0

    def calculate_penalty(self):
        '''
        Calculate penalty for the oligo group solution
        '''
        self.total_score   = 0
        self.total_penalty = 0
        
        #Iterate over each solution
        for key in self.break_solutions: 
            #Get break solution
            break_solution    = self.break_solutions[key]
            
            #If the solution doesnt exist move to next break solution
            if not break_solution:
                continue

            #Update total score
            self.total_score += break_solution.score

            #Initialize the bad break list
            break_solution.bad_list = []

            #Iterate over breaks in breaks list
            for new_break in break_solution.breaks[:-1]:

                #Get neighbor key
                if new_break.neighbor_break:
                    neighbor_break     = new_break.neighbor_break
                    neighbor_oligo_key = neighbor_break.oligo.key
                    
                    if self.break_solutions[neighbor_oligo_key] and neighbor_break in self.break_solutions[neighbor_oligo_key].breaks:
                        break_solution.bad_list.append(new_break)

                        #Update total penalty score
                        self.total_penalty += 1

        #Divide penalty score by 2
        self.total_penalty = int(self.total_penalty/2)

class Scaffold:
    def __init__(self):
        '''
        Scaffold class
        '''

class OligoGroup:
    def __init__(self):
        '''
        Oligo group connected to each other through break point neighborship
        '''
        self.oligos                = None
        self.breaks                = None
        self.group_solutions       = None
        self.best_score_solution   = None
        self.best_penalty_solution = None

    def create_group_solutions(self,num_solutions=200):
        '''
        Create random group solutions
        '''
        #Initialize group solutions
        self.group_solutions = []

        for i in range(num_solutions):
            #Initialize group break solution
            new_group_solution = GroupBreaksolution()
            new_group_solution.break_solutions = {}
            
            for oligo in self.oligos:
                #Randomly pick one solution for an oligo
                if len(oligo.break_solutions) > 0:
                    new_group_solution.break_solutions[oligo.key]  = random.choice(oligo.break_solutions)
                else:
                     new_group_solution.break_solutions[oligo.key] = None

            #Calculate the penalties for each group solution
            new_group_solution.calculate_penalty()

            #Add solution to list
            self.group_solutions.append(new_group_solution)

    def sort_solutions(self):
        '''
        Sort solutions based on the penalty score
        '''
        #1. Sort based on total score
        self.group_solutions.sort(key=lambda x:x.total_score,reverse=True)
        
        #2. Assign best total score solution
        self.best_score_solution = self.group_solutions[0]

        #3. Sort based on total penalty
        self.group_solutions.sort(key=lambda x:x.total_penalty)

        #4. Assign best penalty solution
        self.best_penalty_solution = self.group_solutions[0]

    def print_solutions(self):
        '''
        Print solutions
        '''
        for solution in self.group_solutions:
            print(solution.total_score,solution.total_penalty)


class Strand:
    def __init__(self):
        '''
        Strand class
        '''
        self.origami      = None
        self.sequences    = []
        self.next_strand  = None
        self.oligo        = None
        self.length       = None
        self.final_strand = False

        #Possible break locations
        self.fwd_breaks  = None    #Fwd break locations along a strand
        self.rev_breaks  = None    #Rev break locations along a strand

        #Break rule
        self.break_rule = ['cross','long']
        self.LONG_STRAND_LENGTH = 21
        self.LONG_STRAND_STEP   = 7

    def apply_break_rule(self):
        #Get the break rule
        self.break_rule = self.origami.break_rule
        
        #initialize break location
        self.fwd_breaks = []
        self.rev_breaks = []
        self.all_breaks = []

        #1. Rule number 1 - Crossover only
        if 'cross' in self.break_rule:
            self.rev_breaks.append(-1)
        
        #2. Rule number 2 - Only 3 base away from 5p crossover
        if '3f' in self.break_rule:
            self.fwd_breaks.append(+2)

        #3. Rule number 3 - Only 3 base away from 3p crossover
        if '3r' in self.break_rule:
            self.rev_breaks.append(-4)

        #4. Rule number 4 - Only 3 base away from crossovers
        if '3' in self.break_rule:
            self.fwd_breaks.append(+2)
            self.rev_breaks.append(-4)
        
        #5. Rule number 5 - break long strands
        if 'long' in self.break_rule:
            if self.length >= self.LONG_STRAND_LENGTH:
                #Prepare the fwd and rev break locations
                fwd_breaks = list(range( self.LONG_STRAND_STEP-1,  self.length-self.LONG_STRAND_STEP, self.LONG_STRAND_STEP))
                rev_breaks = list(range(-self.LONG_STRAND_STEP-1, -self.length+self.LONG_STRAND_STEP,-self.LONG_STRAND_STEP))
                
                #Add break locations to list
                self.fwd_breaks += fwd_breaks
                self.rev_breaks += rev_breaks

        #6. For final strand, add the last position(-1)
        if self.final_strand:
            self.rev_breaks.append(-1)

        #Make the breaks array and sort them
        self.fwd_breaks = np.array(sorted(self.fwd_breaks),dtype=int)
        self.rev_breaks = np.array(sorted(self.rev_breaks),dtype=int) + int(self.length)

        #Combine the two arrays
        self.all_breaks = np.sort(np.unique(np.hstack((self.fwd_breaks,self.rev_breaks))))

class Oligo:
    def __init__(self):
        '''
        Oligo class
        '''
        self.origami             = None
        self.strands             = None
        self.type                = 'broken'
        self.crossovers          = []
        self.break_solutions     = []
        self.dont_break          = False
        self.start_break         = None
        self.final_break         = None 
        self.max_score           = 0
        self.initial_score       = 0
        self.end_to_end_edge     = None

    def get_initial_score(self):
        '''
        Get the initial score before the oligo is broken
        '''
        #Create break
        new_edge = BreakEdge()
        
        #Assign
        new_edge.origami = self.origami
        
        #Make the connection
        if self.circular:
            start_break = self.breaks[0]
            final_break = self.breaks[0]
        else:
            start_break = self.start_break
            final_break = self.final_break

        #Get break to break distance
        break_distance  = start_break.get_break_distance(final_break)

        #Assign edge length
        new_edge.edge_length = break_distance

        #Make the connection
        new_edge.make_connection(start_break,final_break)
        self.end_to_end_edge = new_edge

        #Assign the score
        self.initial_score = new_edge.edge_weight

    def keep_best_break_solutions(self):
        '''
        Keep best break solutions
        '''
        #If the solution list is empty return
        if len(self.break_solutions) == 0:
            return

        #Sort the list based on scores
        self.break_solutions.sort(key=lambda x: x.score, reverse=True)

        #Get the maximum score
        self.max_score = self.break_solutions[0].score

        #Find the first index that has a lower score
        max_valid_index = len(self.break_solutions)

        for i in range(len(self.break_solutions)):
            if self.break_solutions[i].score < self.max_score:
                max_valid_index = i
                break

        #Keep only the solutions with maximum score
        self.break_solutions = self.break_solutions[:max_valid_index] 

    def reset_break_order_ids(self,start_break,final_break):
        '''
        Renumber break ids for shortest path algorithm
        '''
        current_break = start_break

        #Break id counter
        id_counter = 0

        #Make the id of current break 0
        current_break.order_id = id_counter

        #Get the next break
        current_break = current_break.next_break
        
        #Iterate over each break
        while current_break:
            #Update id counter
            id_counter += 1
            
            #Update break id
            current_break.order_id = id_counter

            #If we reach the start break quit
            if current_break == final_break:

                #If start and final breaks are the same, make first break number 0
                if start_break == final_break:
                    current_break.order_id = 0
                break

            #Update current breaks
            current_break = current_break.next_break

    def generate_shortest_paths(self,num_solutions=100):
        '''
        Get the shortest paths for the oligo if only it allowed to break it
        '''
        
        #Get k-select parameter
        k_select = self.origami.autobreak.k_select

        #Show oligo being processed
        print('Oligo: ',self.key,'Number of breaks: %d'%(len(self.breaks)))

        if self.dont_break:
            self.break_solutions = []
            return

        #Initialize break solutions
        self.break_solutions = []
        
        if self.circular:
            if len(self.breaks) > 0:
                #Determine number of solutions per oligo
                self.num_solutions_per_oligo = math.ceil(1.0*num_solutions/len(self.breaks))
        else:
            self.num_solutions_per_oligo = num_solutions

        if self.circular:
            
            for current_break in self.breaks:
            
                #Reset break path variables
                self.reset_break_paths()
                #Reset break id numbers
                self.reset_break_order_ids(current_break,current_break)

                #Generate the shortest k-paths
                shortest_k_paths = current_break.get_k_shortest_paths(current_break,self.num_solutions_per_oligo,k_select)

                #Add solutions to solutions list
                self.break_solutions += shortest_k_paths

        else:
    
            #Reset break path variables
            self.reset_break_paths()
            #Reset break id numbers
            self.reset_break_order_ids(self.start_break,self.final_break)

            #Generate shortest paths
            shortest_k_paths = self.start_break.get_k_shortest_paths(self.final_break,self.num_solutions_per_oligo,k_select)

            #Add solutions to solutions list
            self.break_solutions += shortest_k_paths

        #If number of solutions is 1 then remove the worse results
        self.keep_best_break_solutions()

    def reset_break_paths(self):
        for current_break in self.breaks:
            current_break.reset_break_path()


class Crossover:
    def __init__(self):
        '''
        Crossover class
        '''
        self.neighbor   = None
        self.break_node = None

class CrossoverSet:
    def __init__(self):
        '''
        Crossover set class
        '''

class Origami:
    def __init__(self):
        '''
        Origami class
        '''

        self.json_input    = None
        self.oligos        = None
        self.oligo_map     = {}
        self.break_edge_map= {}
        self.oligo_groups  = None

        self.staples       = None
        self.scaffolds     = None
        self.idnums        = None 

    def set_sequence_file(self,sequence_file=None):
        '''
        Set sequence filename
        '''
        self.sequence_file = sequence_file

    def prepare_origami(self):
        '''
        List of commands to prepare origami for break
        '''
        #Get oligos
        self.get_oligos()

        #Reset oligos list
        self.reset_oligos()

        #Read sequence file
        self.read_sequence()

        #Read scaffolds and staples
        self.read_staples()

        #Generate crossovers
        self.generate_crossovers()
        
        #Generate sequences
        self.generate_sequences()

        #Generate break points
        self.generate_break_points()

        #Cluster break points
        self.cluster_oligo_groups()

    def get_cadnano_strand(vh,idx,direction):
        '''
        Get cadnano strand from vh, idx, direction information
        '''
        return self.part.getStrand(direction>0, vh, idx)

    def split_cadnano_strand(vh,idx,direction):
        '''
        Split cadnano strand at vh, idx, direction
        '''
        new_strand = self.part.getStrand(direction>0, vh, idx)
        new_strand.split(idx)

    def initialize(self,input_filename):
        #Initialize cadnano
        app = cadnano.app();
        self.doc = app.document = Document();

        #Assign cadnano input file
        self.json_input = input_filename

        #Read cadnano input file
        self.doc.readFile(self.json_input);

        #Assign part
        self.part = self.doc.activePart()

    def read_oligo(self, oligo, oligo_type='staple'):
        '''
        Read oligo from 5' to 3'
        '''
        #Oligo strand generator
        generator    = oligo.strand5p().generator3pStrand()

        #Create a null strand object for oligo
        previous_strand          = Strand()
        previous_strand.length   = 0
        previous_strand.distance = 0 

        #Create Oligo object
        new_oligo             = Oligo()
        new_oligo.type        = oligo_type
        new_oligo.circular    = oligo.isCircular() 
        new_oligo.null_strand = previous_strand
        new_oligo.length      = oligo.length()
        new_oligo.origami     = self

        #Get 5p strand
        strand5p              = oligo.strand5p()
        idx5p                 = strand5p.idx5Prime()
        vh                    = strand5p.idNum()
        direction             = -1 + 2*strand5p.isForward()

        #Assign key for oligo
        new_oligo.key         = (vh,idx5p,direction)

        #Add oligo to list
        self.oligos[oligo_type].append(new_oligo)
        
        #Add oligo to map
        self.oligo_map[new_oligo.key] = new_oligo
        
        #Get Strand parameters
        for strand in generator:
            #Create new strand
            new_strand                    = Strand()

            new_strand.vh                 = strand.idNum()
            new_strand.idx5p              = strand.idx5Prime()
            new_strand.idx3p              = strand.idx3Prime()
            new_strand.forward            = strand.isForward()
            new_strand.idxLow,new_strand.idxHigh = (new_strand.idx5p,new_strand.idx3p) if new_strand.forward else (new_strand.idx3p,new_strand.idx5p)
            new_strand.dna                = strand.sequence()
            new_strand.direction          = -1 + 2*strand.isForward()
            new_strand.complement_strands = strand.getComplementStrands()[::new_strand.direction]
            new_strand.length             = new_strand.direction*(new_strand.idx3p-new_strand.idx5p)+1
            new_strand.distance           = previous_strand.distance + previous_strand.length
            new_strand.origami            = self

            #Prepare the break points
            new_strand.apply_break_rule()

            #Make the strand connection
            previous_strand.next_strand = new_strand

            #Update previous trand 
            previous_strand = new_strand

        #If oligo is not circular make the last strand final strand
        if not new_oligo.circular:
            previous_strand.final_strand = True

    def generate_crossovers(self):
        '''
        Generate crossover objects
        '''

        #Initialize the crossovers for the origami
        self.crossovers = {}

        for oligo in self.oligos['staple']:
            #Initialize crossovers for the oligo
            oligo.crossovers = {}

            #Get current strand
            current_strand = oligo.null_strand.next_strand

            #Iterate over the strands
            while current_strand:
                if not current_strand.final_strand:

                    #Get the connected strand
                    next_strand   = current_strand.next_strand or oligo.null_strand.next_strand

                    new_crossover = Crossover()
                    new_crossover.vh5p           = current_strand.vh
                    new_crossover.vh3p           = next_strand.vh
                    new_crossover.idx            = current_strand.idx3p
                    new_crossover.direction      = current_strand.direction
                    new_crossover.key            = (new_crossover.vh5p, new_crossover.idx, new_crossover.direction)                    
                    new_crossover.oligo          = oligo
                    new_crossover.current_strand = current_strand
                    new_crossover.next_strand    = next_strand
                    
                    new_crossover.neighbor       = None
                    new_crossover.neighbor_key   = (new_crossover.vh3p,new_crossover.idx+new_crossover.direction,-new_crossover.direction)

                    #Add crossovers to the list
                    oligo.crossovers[new_crossover.key] = new_crossover
                    self.crossovers[new_crossover.key]  = new_crossover

                #Update current strand
                current_strand          = current_strand.next_strand

        #Check crossover neighbors
        for key in self.crossovers:

            #Get neighbor key
            neighbor_key = self.crossovers[key].neighbor_key

            #Make the neighbor connection
            if neighbor_key in self.crossovers:
                self.crossovers[key].neighbor = self.crossovers[neighbor_key]

    def generate_sequences(self):
        '''
        Generate sequences from strands
        '''

        for oligo in self.oligos['staple']:
            #Initialize sequence list for the oligo
            oligo.sequences = []

            current_strand = oligo.null_strand.next_strand

            #Make null sequence object
            previous_sequence = Sequence()
            current_strand.null_sequence = previous_sequence

            while current_strand:
                complement_strands = current_strand.complement_strands

                #Strand sequence positions
                current_strand.sequence_idxLows = [] 

                #Sequence array
                current_strand.sequences = []

                #Go through each complement strand to determine the boundaries
                for strand in complement_strands:
                    comp_idx5p   = strand.idx5Prime()
                    comp_idx3p   = strand.idx3Prime()
                    comp_forward = strand.isForward()

                    #Get the low and high indexes for complementary strand
                    comp_idxLow,comp_idxHigh = (comp_idx5p, comp_idx3p) if comp_forward else (comp_idx3p,comp_idx5p) 

                    #Make sequence object
                    new_sequence         = Sequence()
                    new_sequence.idNum   = current_strand.vh
                    new_sequence.idxLow  = max(current_strand.idxLow ,comp_idxLow)
                    new_sequence.idxHigh = min(current_strand.idxHigh,comp_idxHigh)
                    new_sequence.length  = new_sequence.idxHigh-new_sequence.idxLow+1
                    new_sequence.idx5p, new_sequence.idx3p = (new_sequence.idxLow,new_sequence.idxHigh) if current_strand.forward else (new_sequence.idxHigh,new_sequence.idxLow)
                    new_sequence.forward = current_strand.forward

                    #Assign string indexes
                    new_sequence.strLow  = current_strand.direction*(new_sequence.idx5p - current_strand.idx5p)
                    new_sequence.strHigh = current_strand.direction*(new_sequence.idx3p - current_strand.idx5p)

                    #Assign sequence distance from 5' end of oligo
                    new_sequence.distance= new_sequence.strLow + current_strand.distance

                    #Keep the low position
                    current_strand.sequence_idxLows.append(new_sequence.strLow)

                    #Get the sequence
                    new_sequence.dna = current_strand.dna[new_sequence.strLow:new_sequence.strHigh+1] 

                    #Add new sequence to strand sequence list
                    current_strand.sequences.append(new_sequence)
                    
                    #Add new sequence to oligo sequence list
                    oligo.sequences.append(new_sequence)

                    #Make the next sequence link
                    previous_sequence.next_sequence = new_sequence

                    #Update previuous sequence
                    previous_sequence = new_sequence

                #Make the sequence starting position numpy array
                current_strand.sequence_idxLows = np.array(current_strand.sequence_idxLows)

                #Update current strand
                current_strand = current_strand.next_strand

            #For circular oligos connect last sequence to first sequence
            if oligo.circular:
                previous_sequence.next_sequence = oligo.null_strand.next_strand.sequences[0]

    def generate_break_points(self):
        '''
        Generate break points
        '''
        
        #Initialize breaks
        self.breaks = []

        for oligo in self.oligos['staple']:
            #Initialize oligo breaks and break map
            oligo.breaks = []

            #Assign current strand
            current_strand = oligo.null_strand.next_strand

            #Set break id counter for graph algorithms
            order_id_counter = 0

            #Create a null break point representing the 5'-end of the nucleotide
            previous_break   = BreakNode()
            oligo.null_break = previous_break

            #Initiliaze null break point
            oligo.null_break.break_point = -1
            oligo.null_break.idx         = current_strand.idx5p + oligo.null_break.break_point
            oligo.null_break.vh          = current_strand.vh
            oligo.null_break.direction   = current_strand.direction
            oligo.null_break.distance    = 0
            oligo.null_break.strand      = current_strand
            oligo.null_break.key         = (oligo.null_break.vh,oligo.null_break.idx,oligo.null_break.direction)
            oligo.null_break.sequence    = current_strand.sequences[0]
            oligo.null_break.oligo       = oligo
            oligo.null_break.order_id    = order_id_counter

            #Add null oligo to break list
            oligo.breaks.append(oligo.null_break)

            #If the oligo is not circular make the first break node start node
            if not oligo.circular:
                oligo.start_break = oligo.null_break

            #Iterate over each strand
            while current_strand:
                
                #Iterate through all positions
                for break_position in current_strand.all_breaks:
                    #Update break id counter
                    order_id_counter += 1

                    #Find the sequence for the break point
                    sequence_id = np.searchsorted(current_strand.sequence_idxLows, break_position,side='right') - 1
                    
                    #Make a break
                    new_break   = BreakNode()

                    #Break position
                    new_break.break_point = break_position
                    new_break.idx         = current_strand.idx5p + current_strand.direction*break_position
                    new_break.vh          = current_strand.vh
                    new_break.direction   = current_strand.direction
                    new_break.distance    = current_strand.distance + break_position + 1
                    new_break.key         = (new_break.vh,new_break.idx,new_break.direction)
                    new_break.order_id    = order_id_counter

                    #Assign sequence to break object
                    new_break.sequence    = current_strand.sequences[sequence_id]

                    #Assign strand to new break
                    new_break.strand      = current_strand

                    #Assign oligo 
                    new_break.oligo       = oligo

                    #Check if the break is at a cross-over location
                    if new_break.key in self.crossovers:
                        new_break.crossover = self.crossovers[new_break.key]
                        self.crossovers[new_break.key].break_node = new_break

                    #Assign to previous break
                    previous_break.next_break = new_break

                    #Make the previous connection
                    new_break.previous_break  = previous_break

                    #Update previous break
                    previous_break = new_break

                    #Add break to break list
                    oligo.breaks.append(new_break)

                #Update current strand
                current_strand   = current_strand.next_strand

            #If oligo is circular, connect final node to null break's next break
            if oligo.circular:
                #Update the final node to first node forward connection
                previous_break.next_break = oligo.null_break.next_break
                
                #Update the final node to first node reverse connection
                if oligo.null_break.next_break:
                    oligo.null_break.next_break.previous_break = previous_break

                #Remove the first break node from list
                oligo.breaks.pop(0)

                #Assign final break to null break
                oligo.null_break  = previous_break
            else:
                #If the oligo is not circular make the final break node final node
                oligo.final_break = previous_break

            #Add breaks to origami list
            self.breaks += oligo.breaks

        #Connect the break boints and set break constraints based on connectivity
        self.connect_break_points()

    def connect_break_points(self):
        '''
        Connect break points
        '''
        for oligo in self.oligos['staple']:
            
            #Visit each break object
            for current_break in oligo.breaks:

                #Get the crossover for the current break
                current_crossover = current_break.crossover

                #Get neighbor break
                if current_crossover:
                    current_break.type = 'crossover'
                    #If a neighbor exists make the connection
                    if current_crossover.neighbor:
                        current_break.neighbor_break = current_crossover.neighbor.break_node
                    else:
                        current_break.dont_break = True
                else:
                    current_break.type = 'break'

    def cluster_oligo_groups(self):
        '''
        Cluster oligos based on break connectivity
        '''
        #Initialize oligo groups
        self.oligo_groups = []

        for current_break in self.breaks:
            if not current_break.visited:
                #Create new oligo group
                new_oligo_group = OligoGroup() 
                
                #Perform depth first search starting from current break
                visited_breaks  = current_break.depth_first_search()
                
                #Get the oligos for the visited breaks
                visited_oligos  = set([new_break.oligo for new_break in visited_breaks])
                
                #For all visited breaks make the oligo group assignment
                for new_break in visited_breaks:
                    new_break.oligo_group = new_oligo_group

                #For all visited oligos make the oligo group assignment
                for new_oligo in visited_oligos:
                    new_oligo.oligo_group = new_oligo_group

                #Assign breaks and oligos to oligo group
                new_oligo_group.breaks = visited_breaks
                new_oligo_group.oligos = visited_oligos

                #Add oligo group to the list
                self.oligo_groups.append(new_oligo_group)

    def read_staple(self, staple):
        '''
        Read staple from 5' to 3'
        '''
        self.read_oligo(staple,oligo_type='staple')

    def read_staples(self):
        '''
        Read oligos
        '''

        for staple in self.staples:
            self.read_staple(staple)

    def read_scaffolds(self):
        '''
        Read the scaffolds
        '''
        #Initialize idnums list
        self.idnums = []
        for scaffold in self.scaffolds:
            self.read_scaffold(scaffold)

    def reset_oligos(self):
        '''
        Reset oligos list
        '''
        self.oligos = {'scaffold':[],'staple':[]}

    def read_scaffold(self,scaffold):
        '''
        Read scaffold
        '''
        self.read_oligo(scaffold,oligo_type='scaffold')

    def get_oligos(self):
        '''
        Get oligos
        '''
        self.oligos = self.part.oligos()
        self.oligos = sorted(self.oligos, key=lambda x: x.length(), reverse=True)

        #Initialize the scaffolds and staples
        self.scaffolds = []
        self.staples   = []

        #Iterate over oligos
        for oligo in self.oligos:
            strand5p  = oligo.strand5p()
            vh        = strand5p.idNum()
            isForward = strand5p.isForward()

            #Scaffold criteria: (forward and even-number helix) or (reverse and odd-number helix)
            if int(isForward) + vh%2 == 1:
                self.scaffolds.append(oligo)
            else:
                self.staples.append(oligo)

    def read_sequence(self):
        '''
        Read sequence file
        '''

        if self.sequence_file and os.path.isfile(self.sequence_file):
            
            #Read sequence from file
            f = open(self.sequence_file) 
            self.scaffold_sequence = ''.join([line.strip() for line in f.readlines()])
            f.close()
            
            #Convert to upper case
            self.scaffold_sequence = self.scaffold_sequence.upper()

            #Apply sequence to scaffolds
            for scaffold in self.scaffolds:
                scaffold.applySequence(self.scaffold_sequence)
        else:
            #Assign random sequence
            for scaffold in self.scaffolds:

                scaffold.length
                scaffold.applySequence(self.scaffold_sequence)


    def get_coordinates(self, vh, index):
        '''
        Given a vh and a index, returns (x,y,z)
        for the sidechain pts and backbones fwd and rev
        '''

        #Need to reverse the sign of y-axis(it could be any axis) to make the coordinate system right-handed
        #Cadnano coordinate system is left-handed, not compatible with right-handed A-DNA
         
        axis_pts = self.part.getCoordinates(vh)[0][index]*(1,-1,1)
        fwd_pts  = self.part.getCoordinates(vh)[1][index]*(1,-1,1)
        rev_pts  = self.part.getCoordinates(vh)[2][index]*(1,-1,1)

        return {-1:rev_pts, 0:axis_pts, 1:fwd_pts}

class AutoStaple:
    def __init__(self):
        '''
        Auto staple class
        '''
        self.origami = None

    def decorate(self):
        '''
        Decorate the design
        '''
        self.part.potentialCrossoverMap(0)

class AutoBreak:
    def __init__(self):
        
        #Cadnano parameters
        self.origami                         = None
        self.json_input                      = None
        self.json_output                     = None
        
        #Constraints
        self.MAX_OLIGO_LENGTH                = 60
        self.UPPER_BOUND                     = 60
        self.LOWER_BOUND                     = 21

        #Local and global solutions
        self.NUM_OLIGO_SOLUTIONS             = 1000
        self.NUM_GLOBAL_SOLUTIONS            = 1000

        #CONSTANTS
        self.INFINITY                        = 1E6

        #Score parameters
        self.total_score                     = 0
        self.total_penalty                   = 0 

        #k-shortest path parameter
        self.k_select                        = 'best'

        #Break rule
        self.break_rule                      = ['cross','long']
        
        #Optimization function
        self.optim_args                      = [['14'],['glength','45','5']]

        #Length gaussian parameters
        self.optim_length_mean               = 45
        self.optim_length_tolerance          = 5

        #Tm gaussian parameters
        self.optim_Tm_mean                   = 60
        self.optim_Tm_tolerance              = 5

        #Max sequence gaussian parameters
        self.optim_maxseq_mean               = 14
        self.optim_maxseq_tolerance          = 2

        #Set function dictionary
        self.optim_funcs_dict                = {'14'     :self._optimize_14,
                                                'Tm'     :self._optimize_Tm,
                                                'glength':self._gauss_length,
                                                'gmaxseq':self._gauss_maxseq,
                                                'gTm'    :self._gauss_Tm}
        #Set optimization params 
        self.optim_params_dict              = {'14'     :[],
                                                'Tm'     :[],
                                                'glength':[self.optim_length_mean,self.optim_length_tolerance],
                                                'gmaxseq':[self.optim_maxseq_mean,self.optim_maxseq_tolerance],
                                                'gTm'    :[self.optim_Tm_mean    ,self.optim_Tm_tolerance]}



    def set_solution_nums(self,solutions_per_oligo=1000,global_solutions=1000):
        '''
        Set solution numbers
        '''
        self.NUM_OLIGO_SOLUTIONS  = solutions_per_oligo
        self.NUM_GLOBAL_SOLUTIONS = global_solutions 

    def set_k_select(self,k_parameter='best'):
        '''
        Set k-select value
        '''
        self.k_select = k_parameter

    def set_break_rule(self,new_break_rule=['cross','long']):
        '''
        Set break rule
        '''
        self.break_rule         = [rule for rule in new_break_rule]
        self.origami.break_rule = self.break_rule  

    def run_autobreak(self):
        '''
        Run basic autobreak protocol
        '''
        #Make break-break graph
        self.initialize()

        #Determine initial scores
        self.determine_initial_scores()

        #Break oligos
        self.create_oligo_solutions()

        #Create group solutions
        self.create_group_solutions()

        #Sort and print group solutions
        self.sort_group_solutions()

    def initialize(self):
        '''
        Initialize the connectivity maps
        '''
        
        for oligo in self.origami.oligos['staple']:

            #Check oligo length, if the length is within length limits dont break it
            if oligo.length <= self.MAX_OLIGO_LENGTH and not oligo.circular:
                oligo.dont_break = True
            
            #Visit each break object
            for current_break in oligo.breaks:
                #Initialize the break edges
                current_break.break_edges = []

                #Get next break
                next_break = current_break.next_break

                #Iterate over the breaks
                while next_break:
                    
                    #Determine break to break distance
                    break_distance = current_break.get_break_distance(next_break)
                    
                    if break_distance >= self.LOWER_BOUND and break_distance <= self.UPPER_BOUND:
                        #Create break
                        new_edge = BreakEdge()
                        
                        #Assign origami
                        new_edge.origami = self.origami

                        #Assign edge length
                        new_edge.edge_length = break_distance

                        #Make the connection
                        new_edge.make_connection(current_break,next_break)

                        #Set edge weight
                        new_edge.edge_weight = self.optimize(new_edge)

                        #Add directed edge to current break's edges
                        current_break.break_edges.append(new_edge)

                        #Add break edge to edge map
                        self.origami.break_edge_map[current_break.key+next_break.key] = new_edge

                    #Stop criteria
                    if break_distance > self.UPPER_BOUND or next_break==current_break:
                        break

                    next_break = next_break.next_break


    def determine_initial_scores(self):
        '''
        Determine starting scores
        '''
        for oligo in self.origami.oligos['staple']:
            oligo.get_initial_score()

    def create_oligo_solutions(self):
        '''
        Break oligos
        '''
        for oligo in self.origami.oligos['staple']:
            if not oligo.dont_break:
                oligo.generate_shortest_paths(self.NUM_OLIGO_SOLUTIONS)

    def create_group_solutions(self):
        '''
        Create oligo group solutions
        '''
        for oligo_group in self.origami.oligo_groups:
            oligo_group.create_group_solutions(self.NUM_GLOBAL_SOLUTIONS)

    def sort_group_solutions(self):
        '''
        Sort group solutions based on penalty score
        '''
        for oligo_group in self.origami.oligo_groups:
            oligo_group.sort_solutions()
            
    def combine_group_solutions(self):
        '''
        Combine group solutions
        '''
        self.best_score_solutions   = []
        self.best_penalty_solutions = []

        for oligo_group in self.origami.oligo_groups:
            
            #Add the best and best penalty solutions
            self.best_score_solutions.append(oligo_group.best_score_solution)
            self.best_penalty_solutions.append(oligo_group.best_penalty_solution)

    def break_oligos(self):
        '''
        Oligo breaking routine
        '''
        

        

    def set_optimization_func(self,func_args):
        '''
        Set optimization function
        '''
        self.optim_args         = func_args
        self.optim_args_funcs   = [function[0]  for function in func_args]
        self.optim_args_params  = [[int(x) for x in function[1:]] for function in func_args]
        self.optimize_func_list = []

        #Set optimization function parameters
        for i in range(len(self.optim_args_funcs)):
            func   = self.optim_args_funcs[i]
            params = self.optim_args_params[i]
            
            #Make the optimize function
            if func in self.optim_funcs_dict:
                self.optimize_func_list.append(self.optim_funcs_dict[func])

            #Assign function parameters
            if func in self.optim_params_dict and len(params) <= len(self.optim_params_dict[func]): 
                for j in range(len(params)):
                    self.optim_params_dict[func][j] = params[j]

    def optimize(self,edge):
        '''
        final optimization function
        '''
        return np.product(np.array([func(edge) for func in self.optimize_func_list])) 
    
    def _optimize_14(self,edge):
        '''
        Optimization function 14
        '''
        return edge.edge_has14
    
    
    def _optimize_Tm(self,edge):
        '''
        Optimization function Tm
        '''
        return edge.edge_Tm
    
    
    def _gauss_length(self,edge):
        '''
        Optimization function gauss length
        '''

        return np.exp(-(edge.edge_length-self.optim_params_dict['glength'][0])**2/self.optim_params_dict['glength'][1]**2)

    
    def _gauss_Tm(self,edge):
        '''
        Optimization function gauss Tm
        '''

        return np.exp(-(edge.edge_Tm-self.optim_params_dict['gTm'][0])**2/self.optim_params_dict['gTm'][1]**2)

    
    def _gauss_maxseq(self,edge):
        '''
        Optimization function gauss Tm
        '''

        return np.exp(-(edge.edge_maxseq-self.optim_params_dict['gmaxseq'][0])**2/self.optim_params_dict['gmaxseq'][1]**2)


class BreakEdge:
    def __init__(self):
        '''
        Break edge class
        '''
        self.origami       = None 
        self.current_break = None 
        self.next_break    = None
        self.edge_weight   = None
        self.edge_length   = None
        self.edge_maxseq   = None
        self.edge_Tm       = None
        self.LOW_TM        = 60

        #state parameters
        self.active        = True

        #Length parameters
        self.edge_has14    = None
        self.edge_num14    = None
        
        #Tm parameters
        self.edge_hasTm    = None
        self.edge_numTm    = None

        self.sequence_list = None

        #Loop parameter
        self.isloop        = False

    def set_edge_weight(self):
        '''
        Set edge weight
        '''
        self.edge_weight = self.origami.autobreak.optimize(self)

    def is_valid(self):
        '''
        Determine if edge is valid
        '''
        return self.active and not self.current_break.dont_break and not self.next_break.dont_break

    def make_connection(self,from_break,to_break):
        '''
        Make the connection between self and another edge
        '''
        #Initialize sequence and dna list
        self.sequence_list = []
        self.dna_list      = []

        #Set the break nodes
        self.current_break = from_break
        self.next_break    = to_break

        #Make the sequence list
        self.sequence_list.append(self.current_break.sequence)

        #Check if the breaks are in consecutive positions on the same sequence
        if not (self.current_break.sequence == self.next_break.sequence and \
            self.current_break.direction*(self.next_break.idx - self.current_break.idx) > 0):
            
            #Get forward link
            next_sequence = self.current_break.sequence.next_sequence

            #Iterate over all the sequences
            while next_sequence:
                self.sequence_list.append(next_sequence)

                #Stop if the sequence is same as the final sequence
                if next_sequence == self.next_break.sequence:
                    break
                
                #Update next sequence
                next_sequence = next_sequence.next_sequence

        #Iterate over sequence list
        if len(self.sequence_list) == 1:
            start_point  = self.current_break.break_point+1
            final_point  = self.next_break.break_point+1
            dna_sequence = self.current_break.sequence.dna[start_point:final_point]

            self.dna_list.append(dna_sequence)
        else:
            #1. Get the 5' sequence 
            start_point  = self.current_break.break_point+1
            final_point  = self.current_break.sequence.strHigh+1
            dna_sequence = self.current_break.sequence.dna[start_point:final_point]

            self.dna_list.append(dna_sequence)
            #2. Get the sequences in between
            for sequence in self.sequence_list[1:-1]:
                self.dna_list.append(sequence.dna)

            #3. Get the 3' sequence
            start_point  = self.next_break.sequence.strLow
            final_point  = self.next_break.break_point+1
            dna_sequence = self.next_break.sequence.dna[start_point:final_point]

            self.dna_list.append(dna_sequence)
        
        #Remove empty sequences
        self.dna_list = [dna for dna in self.dna_list if len(dna)]

        #Determine Tm
        self.Tm_list     = np.array([utilities.sequence_to_Tm(dna) for dna in self.dna_list])

        #Determine lengths
        self.length_list = np.array([len(dna) for dna in self.dna_list])

        #Determine the edge weights
        self.edge_Tm     = max(self.Tm_list)
        
        #Length parameters
        self.edge_maxseq = max(self.length_list)
        self.edge_num14  = np.sum(self.length_list >= 14)  
        self.edge_has14  = int(self.edge_num14 > 0)

        #Tm parameters
        self.edge_numTm  = np.sum(self.Tm_list >= self.LOW_TM)
        self.edge_hasTm  = np.sum(self.edge_numTm > 0)

        #Set edge weight
        self.set_edge_weight()
        
        #Assign loop parameters and loop edge
        if from_break == to_break:
            self.isloop                  = True
            self.current_break.loop_edge = self

class BreakPath:
    def __init__(self,break_node,break_edge=None,score=0):
        '''
        Break path object
        '''
        self.break_node = break_node
        self.break_edge = break_edge
        self.score      = score

class BreakNode:
    def __init__(self):
        '''
        Break node class
        '''
        self.crossover        = None
        self.next_break       = None
        self.previous_break   = None 
        self.neighbor_break   = None
        self.connected_breaks = None
        self.break_edges      = None
        self.edge_nodes       = None
        self.type             = None
        self.sequence         = None
        self.crossover        = None
        self.loop_edge        = None

        #State parameters
        self.visited          = False
        self.active           = True
        self.break_state      = None  #broken for break,  not-broken for no break
        self.dont_break       = False #If set to True, keep the break not-broken

        #Cluster info
        self.oligo_group      = None

        #Graph parameter
        self.score               = 0
        self.best_path_nodes     = None
        self.shortest_paths      = None
        self.k_potential_paths   = None
        self.k_shortest_paths    = None
        self.traverse_path       = None
        self.shortest_score      = 0
        self.order_id            = None

    def reset_break_path(self):
        self.order_id         = -1
        self.traverse_path    = []
        self.best_path_nodes  = []
        self.score            = 0
        self.visited          = False
        self.shortest_paths   = []

    def reset_k_paths(self):
        self.k_shortest_paths = []
        self.k_potential_paths= []

    def is_break_edge_possible(self,other_break):
        '''
        Check if break edge is possible from self to to another break node for a circular oligo
        '''
        order_difference = self.get_break_order_difference(other_break)

        #Criteria for a proper connection (no loop) in the right direction 
        return order_difference > 0 or (not order_difference ==0 and order_difference == -self.order_id)

    def get_loop_edge(self):
        '''
        Return loop edge
        '''

        return self.loop_edge

    def get_break_distance(self,other_break):
        '''
        Break to break distance
        '''
        return (other_break.distance - self.distance)%self.oligo.length + self.oligo.length*(self==other_break)*self.oligo.circular

    def get_break_order_difference(self,other_break):
        '''
        Return break to order id difference
        '''
        return other_break.order_id - self.order_id


    def get_valid_edge_nodes(self):
        '''
        Get break nodes connected by edges
        '''
        #Get the connected break nodes
        return [break_edge.next_break for break_edge in self.break_edges if not break_edge.next_break.dont_break] 


    def get_valid_edges(self):
        '''
        Get edges that lead to break nodes that can be broken
        '''

        #Get the connected break nodes
        return [break_edge for break_edge in self.break_edges if break_edge.is_valid()]

    def get_k_shortest_paths(self,final_break,k_num=10,k_select='best'):
        '''
        Get k-shortest path results
        '''
        #Initialize k-shortest paths
        self.k_shortest_paths = []

        #1. Get the shortest paths
        shortest_paths = self.get_shortest_paths(final_break,num_solutions=1)

        #If the here is no path found return empty list
        if len(shortest_paths) == 0:
            return self.k_shortest_paths

        #2.Add best path to k-path list
        self.k_shortest_paths  = [shortest_paths[0]]
        self.k_potential_paths = []
        num_k_solutions        = 1 

        #Check the final score of the path
        if shortest_paths[0].score == 0:
            return self.k_shortest_paths

        #3. Make the paths
        while num_k_solutions < k_num:
            
            #Get the last best path
            last_solution = self.k_shortest_paths[-1]


            #Iterate through the edges
            for i in range(len(last_solution.edges[:-1])):
                #Inactive edge list
                inactive_edges = []

                for j in range(len(self.k_shortest_paths)):
                    
                    if last_solution.is_identical(self.k_shortest_paths[j],i):
                        inactive_edges.append(last_solution.edges[i])
                        
                        #Make the edge inactive
                        last_solution.edges[i].active = False
                
                #Get last break
                last_break = last_solution.breaks[i]

                #Reset break paths and find solution
                self.oligo.reset_break_paths()
                self.oligo.reset_break_order_ids(last_break,final_break)

                #Get sub solution
                sub_solution = last_solution.breaks[i].get_shortest_paths(final_break,num_solutions=1)
                
                #If there is no solution continue
                if len(sub_solution) == 0:
                    continue

                sub_solution[0].initialize()

                #Make the edges active again
                for edge in inactive_edges:
                    edge.active = True

                #Get final score from last solution
                pre_score = last_solution.scores[:i][-1] if len(last_solution.scores[:i]) > 0 else 0

                #Update subsolution scores
                sub_solution[0].scores = [score+pre_score for score in sub_solution[0].scores]

                #Create new solution
                new_solution = OligoBreakSolution()
                new_solution.start_break  = self
                new_solution.final_break  = final_break
                new_solution.edges        = last_solution.edges[:i]  + sub_solution[0].edges
                new_solution.breaks       = last_solution.breaks[:i] + sub_solution[0].breaks 
                new_solution.scores       = last_solution.scores[:i] + sub_solution[0].scores
                new_solution.score        = new_solution.scores[-1]

                #Add potential solution if it doesnt exist in k-shortest paths
                potential_solution_exists = False
                for solution in self.k_shortest_paths:
                    if solution.is_identical(new_solution):
                        potential_solution_exists = True
                        break

                #Add new solution to potential paths
                if not potential_solution_exists:
                    self.k_potential_paths.append(new_solution)

            #Check if the potential path list is empty, if empty quit
            if len(self.k_potential_paths) == 0:
                break

            #Sort the potential paths and add the best one to k-shortest paths list
            self.k_potential_paths.sort(key=lambda x:x.score, reverse=True)

            if k_select == 'best':
                solution_index = 0
            else:
                #Add random one to the k-shortest paths
                solution_index = random.randint(0, len(self.k_potential_paths)-1)
            
            #Add item to k-shortest path list
            self.k_shortest_paths.append(self.k_potential_paths[solution_index])
            
            #Update num solutions
            num_k_solutions += 1

            #Remove the result from potential paths list
            self.k_potential_paths.pop(solution_index)

        return self.k_shortest_paths

    def get_shortest_paths(self,final_break,num_solutions=10):
        '''
        Find the shortest path between current and final break points
        '''

        #Initialize the set and stack
        stack = [self]
        
        while stack:
            #Pop the break node
            new_break = stack.pop(0)

            #If current node is final break, quit
            if final_break.visited and new_break == final_break:
                break
            
            #Add break to visited list
            new_break.visited = True
            
            #Get valid break edges
            valid_break_edges = new_break.get_valid_edges()

            #Update the scores for connected breaks
            for break_edge in valid_break_edges:
                
                #If id difference is in wrong direction and if it is a loop, discard the edge
                if not new_break.is_break_edge_possible(break_edge.next_break):
                    continue

                #Determine the new score
                new_score = new_break.score + break_edge.edge_weight

                #Make new break Path object
                new_break_path = BreakPath(new_break, break_edge, new_score)

                #If new score is higher than the previous one, make a new list
                if new_score > break_edge.next_break.score:
                    break_edge.next_break.best_path_nodes = [new_break_path]
                    break_edge.next_break.score           = new_score
                #Otherwise add to existng path list
                elif new_score == break_edge.next_break.score:
                    break_edge.next_break.best_path_nodes.append(new_break_path)

                #Add next break to connected breaks list
                if not break_edge.next_break.visited: stack.append(break_edge.next_break)

        #Finally compare with the loop connection
        if self == final_break and self.loop_edge and self.loop_edge.is_valid():
            #Make a loop break path object
            loop_break_path = BreakPath(self, self.loop_edge, self.loop_edge.edge_weight)
            
            if self.loop_edge.edge_weight > final_break.score:
                final_break.best_path_nodes = [loop_break_path]
                final_break.score           =  self.loop_edge.edge_weight
            elif self.loop_edge.edge_weight == final_break.score:
                final_break.best_path_nodes.append(loop_break_path)
        
        #Return best paths
        final_break.shortest_paths = final_break.traverse_best_paths(self,num_solutions)
        
        return final_break.shortest_paths

    def traverse_best_paths(self, start_break, num_solutions=1):
        
        #Initialize shortest path
        self.shortest_paths = []

        #Previous break on path
        stack = self.best_path_nodes

        #Make the traverse path for starting node empty
        self.traverse_path = []

        #Make final node break path object
        final_break_path = BreakPath(self,None,self.score)

        #Initalize the traverse path for the initial nodes
        for break_node_path in stack:
            break_node_path.break_node.traverse_path = self.traverse_path + [final_break_path]

        while stack:

            #Check number of shortest paths
            if len(self.shortest_paths) == num_solutions:
                break

            #Pop the break path
            new_break_path = stack.pop()
            
            #Get break node
            new_break = new_break_path.break_node

            if new_break == start_break:
                new_break_solution              = OligoBreakSolution()
                new_break_solution.start_break  = start_break
                new_break_solution.final_break  = self
                new_break_solution.break_paths  = new_break.traverse_path + [new_break_path]  
                new_break_solution.score        = self.score

                #Initialize the solution
                new_break_solution.initialize()

                self.shortest_paths.append(new_break_solution)
                continue

            #Get the new path nodes
            new_path_nodes = new_break.best_path_nodes

            #Iterate through each node
            for new_path_node in new_path_nodes:
                new_path_node.break_node.traverse_path = new_break.traverse_path + [new_break_path]

            #Extend the group with new breaks
            stack.extend(new_path_nodes)
        
        return self.shortest_paths

    def get_connected_breaks(self):
        '''
        Return connected breaks
        '''
        #Initialize connected breaks
        self.connected_breaks = []
        
        #1. Add the directly connected next break
        if self.next_break:
            self.connected_breaks.append(self.next_break)

        #2. Add the directly connected previous break
        if self.previous_break:
            self.connected_breaks.append(self.previous_break)

        #3. Add the neighbor break
        if self.neighbor_break:
            self.connected_breaks.append(self.neighbor_break)

        return self.connected_breaks

    def depth_first_search(self):
        '''
        Depth-first-search graph traverse algorithm to find connected components
        '''
        #Initialize the set and stack
        visited, stack = set(), [self]
        while stack:
            #Pop the break node
            new_break = stack.pop()
            if new_break not in visited:
                #Add break to visited list
                visited.add(new_break)

                #Make the break visited 
                new_break.visited = True
                
                #Get the connected breaks
                connected_breaks = new_break.get_connected_breaks()

                #Extend the group with new breaks
                stack.extend(set(connected_breaks) - visited)
        
        return list(visited)

#Parse functions
def parse_break_rule(break_rule):
    '''
    Parse break rule
    '''
    return break_rule.split('.')

def parse_optim_function(function_input):
    '''
    Parse optimizatiotion function input
    '''
    #Get the function groups
    groups = function_input.split('.')

    #Get the functions and its parameters
    functions = [group.split(':') for group in groups]

    return functions

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i",   "--input",    type=str,     help="Cadnano json file" )
    parser.add_argument("-s",   "--sequence", type=str,     help="Sequence file in txt", default=None )
    parser.add_argument("-r",   "--rule",     type=str,     help="Break rule",default='cross.long')
    parser.add_argument("-f",   "--func",     type=str,     help="Optimization function", default='14.glength:45:10')
    parser.add_argument("-o",   "--osol",     type=int,     help="Solutions per oligo", default=1000)
    parser.add_argument("-g",   "--gsol",     type=int,     help="Global solutions", default=1000)
    parser.add_argument("-k",   "--kselect",  type=str,     help="selection method for k-shortest path", choices=['best','random'],default='best')

    args = parser.parse_args()

    #Assign the parameters
    input_filename          = args.input
    sequence_filename       = args.sequence 
    break_rule              = parse_break_rule(args.rule)
    optimization_func       = parse_optim_function(args.func)
    sols_per_oligo          = args.osol
    global_solutions        = args.gsol
    k_selection_method      = args.kselect

    #Check if input file exists
    if not os.path.isfile(input_filename):
        sys.exit('Input file does not exist!')

    #Create Origami object
    new_origami = Origami()

    #Create Autobreak object and make cross assignments
    new_autobreak = AutoBreak()

    #Cross assignments
    new_origami.autobreak = new_autobreak
    new_autobreak.origami = new_origami

    #Assign break rule
    new_autobreak.set_break_rule(break_rule)

    #Set optimization function
    new_autobreak.set_optimization_func(optimization_func)

    #Set solution numbers
    new_autobreak.set_solution_nums(sols_per_oligo,global_solutions)

    #Set k-shortest path selection method
    new_autobreak.set_k_select(k_selection_method)

    #Initialize origami object
    new_origami.initialize(input_filename)

    #Set sequence filename
    new_origami.set_sequence_file(sequence_filename)

    #Prepare origami for autobreak
    new_origami.prepare_origami()

    #Run autobreak
    new_autobreak.run_autobreak()


if __name__ == "__main__":
  main()








