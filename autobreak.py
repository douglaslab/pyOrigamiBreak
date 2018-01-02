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
import matplotlib
import argparse

from cadnano.document import Document
from functools import reduce
from math import ceil
from matplotlib import cm


#Sequence enthalpy/entropy values at 37C-1M NaCl from SANTALUCIA-HICKS paper
#Units are in kcal/mol for Enthalpy and cal/mol.Kelvin
TUPLE_ENTHALPY = {  'AA'  : -7.6,
                    'TT'  : -7.6,
                    
                    'AT'  : -7.2,
                    
                    'TA'  : -7.2,
                    
                    'CA'  : -8.5,
                    'TG'  : -8.5,
                    
                    'GT'  : -8.4,
                    'AC'  : -8.4,
                    
                    'CT'  : -7.8,
                    'AG'  : -7.8,

                    'GA'  : -8.2,
                    'TC'  : -8.2,

                    'CG'  :-10.6,
                    
                    'GC'  : -9.8,
                    
                    'GG'  : -8.0,
                    'CC'  : -8.0,
                    
                    'INIT':  0.2,
                    'PENL':  2.2
}

TUPLE_ENTROPY  = {  'AA'  : -21.3,
                    'TT'  : -21.3,
                    
                    'AT'  : -20.4,
                    
                    'TA'  : -21.3,
                    
                    'CA'  : -22.7,
                    'TG'  : -22.7,

                    'GT'  : -22.4,
                    'AC'  : -22.4,

                    'CT'  : -21.0,
                    'AG'  : -21.0,

                    'GA'  : -22.2,
                    'TC'  : -22.2,

                    'CG'  : -27.2,
                    
                    'GC'  : -24.4,
                    
                    'GG'  : -19.9,
                    'CC'  : -19.9,

                    'INIT': -5.7 ,
                    'PENL':  6.9
}

class AutoBreakSolution:
    def __init__(self):
        self.oligo_type            = None
        self.oligo_num             = None
        self.constraints_satisfied = None
        self.valid                 = None
        self.solution_path         = []
        self.crossover_path        = []
        self.score                 = []          

class Structure:
    def __init__(self):
        self.segments              = []
        self.oligos                = []
        self.xyz                   = []

class AutoBreak:
    def __init__(self):
        
        #Cadnano file parameters
        self.doc                             = None
        self.part                            = None
        self.json_input                      = None
        self.json_output_autobreak           = None
        self.json_output_longstrands_broken  = None
        
        #Cadnano oligo/strand information
        self.scaffold                = None
        self.sequence_file           = ''
        self.scaffold_sequence       = ''
        self.oligos                  = []
        self.strands                 = []
        self.strand_lengths          = []
        self.strand_sequences        = []
        self.strand_Tm               = []
        self.crossover_nums_to_names = []
        self.crossover_state_matrix  = []

        #Crossover - interaction info
        self.oligo_oligo_interaction    = {}
        self.crossover_names_to_states  = {}
        self.crossover_constraints      = []

        self.crossover_names_to_nums     = {}
        self.crossover_names_to_neighbor = {}

        #Constraints
        self.UPPER_BOUND               = 60
        self.LOWER_BOUND               = 23
        self.BREAK_LONG_STRANDS        = False

        #Lookup tablesspot
        self.DISTANCE_MATRIX           = []
        self.UPPER_BOUND_MATRIX        = []
        self.LOWER_BOUND_MATRIX        = []
        self.SLENGTH_CONSTRAINT_MATRIX = []
        self.STRAND_MATRIX             = []

        #Thermodynamics constants
        self.SCAF_CONC                 = 10E-9
        self.STAP_CONC                 = 100E-9
        self.MAGNESIUM                 = 10E-3 
        self.R                         = 0.593/298.15

        #Local and global solutions
        self.local_solutions           = [[],[]]
        self.local_solutions_complete  = True
        self.pruned_local_solutions    = {}
        self.best_global_solutions     = {}
        self.best_global_solution      = {}
        self.NUM_SOLUTIONS_PER_OLIGO   = 100
        self.NUM_GLOBAL_ITERATIONS     = 50

        #CONSTANTS
        self.INFINITY                   = 1E6

    def crossover_position(self,helix1):
        '''
        Returns the position of crossover
        '''
        
        if helix1%2 == 0:
            position = 'right'
        else:
            position = 'left'

        return position

    def break_broken_long_strands(self):
        '''
        Break 21 and 28 length strands for broken oligos
        '''
        break_points = []
        for i in range(len(self.strand_lengths[0])):
            if sum(self.strand_lengths[0][i]) > self.LOWER_BOUND:
                strands_21 = np.nonzero(np.array(self.strand_lengths[0][i]) == 21)[0]
                strands_28 = np.nonzero(np.array(self.strand_lengths[0][i]) == 28)[0]

                #For every 28 long strands, mid-point is the break point
                for strand_28 in strands_28:
                    #Get the crossover locations
                    crossover_1 = [int(x) for x in self.crossover_nums_to_names[0][i][strand_28-1].split('-')]
                    crossover_2 = [int(x) for x in self.crossover_nums_to_names[0][i][strand_28  ].split('-')]

                    helix_num   = crossover_1[2]
                    mid_break   = int((crossover_2[1] - crossover_1[3])/27*13+crossover_1[3])

                    break_points.append([helix_num,mid_break])

                #For every 21 strand, find the best break point
                for strand_21 in strands_21:
                    left_14s  = np.nonzero(np.array(self.strand_lengths[0][i][:strand_21]) >= 14)[0]
                    right_14s = np.nonzero(np.array(self.strand_lengths[0][i][strand_21:]) >= 14)[0]
                    
                    #If there is no 14 in the oligo, skip
                    if len(left_14s) == 0 or len(right_14s) == 0:
                        continue

                    #Get the crossover locations
                    crossover_1 = [int(x) for x in self.crossover_nums_to_names[0][i][strand_21-1].split('-')]
                    crossover_2 = [int(x) for x in self.crossover_nums_to_names[0][i][strand_21  ].split('-')]

                    helix_num   = crossover_1[2]

                    right_break = int((crossover_2[1] - crossover_1[3])/20*13+crossover_1[3])
                    left_break  = int((crossover_2[1] - crossover_1[3])/20*6 +crossover_1[3])
                    
                    if len(left_14s) == 0:
                        break_points.append([helix_num,left_break])
                    elif len(right_14s) == 0:
                        break_points.append([helix_num,right_break])
                    elif strand_21 - max(left_14s) > min(right_14s) - strand_21:
                        break_points.append([helix_num,left_break])
                    else:
                        break_points.append([helix_num,right_break])

        #Split the strands
        for break_point in break_points:
            strand = self.get_strand(break_point[0],break_point[1])
            strand.split(break_point[1])

    def break_circular_long_strands(self):
        '''
        Break 21 and 28 length strands for circular oligos
        '''
        break_points = []
        for i in range(len(self.strand_lengths[1])):
            if sum(self.strand_lengths[1][i]) > self.LOWER_BOUND:
                
                strands_21 = np.nonzero(np.array(self.strand_lengths[1][i]) == 21)[0]
                strands_28 = np.nonzero(np.array(self.strand_lengths[1][i]) == 28)[0]

                #Number of strands per oligo
                num_strands = len(self.crossover_nums_to_names[1][i])

                #For every 28 long strands, mid-point is the break point
                for strand_28 in strands_28:
                    #Get the crossover locations
                    crossover_1 = [int(x) for x in self.crossover_nums_to_names[1][i][(strand_28-1)%num_strands].split('-')]
                    crossover_2 = [int(x) for x in self.crossover_nums_to_names[1][i][strand_28%num_strands].split('-')]

                    helix_num   = crossover_1[2]
                    mid_break   = int((crossover_2[1] - crossover_1[3])/27*13+crossover_1[3])

                    break_points.append([helix_num,mid_break])

                #For every 21 strand, find the best break point
                for strand_21 in strands_21:
                    left_14s  = np.nonzero(np.array(self.strand_lengths[1][i][:strand_21]+self.strand_lengths[1][i][strand_21:]) >= 14)[0]
                    right_14s = np.nonzero(np.array(self.strand_lengths[1][i][strand_21:]+self.strand_lengths[1][i][:strand_21]) >= 14)[0]
                    
                    #If there is no 14 in the oligo, skip
                    if len(left_14s) == 0 or len(right_14s) == 0:
                        continue

                    #Get the crossover locations
                    crossover_1 = [int(x) for x in self.crossover_nums_to_names[1][i][(strand_21-1)%num_strands].split('-')]
                    crossover_2 = [int(x) for x in self.crossover_nums_to_names[1][i][strand_21%num_strands].split('-')]

                    helix_num   = crossover_1[2]

                    right_break = int((crossover_2[1] - crossover_1[3])/20*13+crossover_1[3])
                    left_break  = int((crossover_2[1] - crossover_1[3])/20*6+crossover_1[3])
                    
                    if len(left_14s) == 0:
                        break_points.append([helix_num,left_break])
                    elif len(right_14s) == 0:
                        break_points.append([helix_num,right_break])
                    elif strand_21 - max(left_14s) > min(right_14s) - strand_21:
                        break_points.append([helix_num,left_break])
                    else:
                        break_points.append([helix_num,right_break])
        
        #Split the strands
        for break_point in break_points:
            strand = self.get_strand(break_point[0],break_point[1])
            strand.split(break_point[1])

    def states_to_matrix(self):
        '''
        Convert crossover states variables to state matrix
        '''
        for crossover_name in self.crossover_names_to_ids.keys():
            oligo_type,oligo_num,crossover_num = self.crossover_names_to_ids[crossover_name]
            self.state_matrix[oligo_type][oligo_num][crossover_num] = self.crossover_states[crossover_name]

    def parse_oligos(self):
        '''
        Parse oligos
        '''
        oligos_sorted_by_length = sorted(self.part.oligos(), key=lambda x: x.length(), reverse=True)
        
        self.scaffold                       = oligos_sorted_by_length[0]
        self.oligos                         = oligos_sorted_by_length[1:]
        self.strand_lengths                 = [[],[]]
        self.strand_sequences               = [[],[]]
        self.crossover_nums_to_names        = [[],[]]
        self.crossover_state_matrix         = [[],[]]
        self.crossover_names_to_neighbors   = {}
        self.crossover_names_to_nums        = {}
        self.crossover_names_to_states      = {}

        #If there is a scaffold sequence entered apply the sequence
        if len(self.scaffold_sequence) > 0:
            self.apply_sequence()

        #Oligo number
        oligo_num                = [0,0] 

        for oligo in self.oligos:

            #Determine oligo type: 0 for broken, 1 for circular
            oligo_type = int(oligo.isCircular())

            unbroken_length  = oligo.length()
            strand_lengths   = oligo.getStrandLengths()
            current_strand5p = oligo.strand5p()
            start_strand     = current_strand5p

            #Get the sequence
            strand_sequences = []
            
            #Keep track of strand and crossover number
            strand_num       = 0
            crossover_num    = 0

            #Crossovers
            crossovers       = []

            while current_strand5p is not None and (current_strand5p != start_strand or strand_num == 0):

                if current_strand5p.connection3p() is not None:
                    helix1    = current_strand5p.idNum()
                    position1 = current_strand5p.idx3Prime()
                    helix2    = current_strand5p.connection3p().idNum()
                    position2 = current_strand5p.connection3p().idx5Prime()
                    
                    #Cross over info
                    crossover    = [helix1,position1,helix2,position2]
                    if helix1%2 == 0:
                        neighbor  = [helix2,position2-1,helix1,position1-1]
                    else:
                        neighbor  = [helix2,position2+1,helix1,position1+1]

                    #Determine the ids
                    crossover_name = '-'.join([str(x) for x in crossover]) 
                    neighbor_name  = '-'.join([str(x) for x in neighbor])

                    #Assign crossover name to id
                    self.crossover_names_to_nums[crossover_name]      = [oligo_type,oligo_num[oligo_type],crossover_num]

                    #Assign neigbor
                    self.crossover_names_to_neighbors[crossover_name] = neighbor_name

                    #Assign to crossover states
                    self.crossover_names_to_states[crossover_name] = 1 

                    #Add to crossover list
                    crossovers.append(crossover_name)
                    crossover_num += 1

                    #Add the strand sequence to strand sequence list
                    strand_sequences.append(current_strand5p.sequence())

                current_strand5p = current_strand5p.connection3p()
                strand_num += 1

            self.strand_sequences[oligo_type].append(strand_sequences)
            self.strand_lengths[oligo_type].append(strand_lengths)
            self.crossover_nums_to_names[oligo_type].append(crossovers)
            self.crossover_state_matrix[oligo_type].append(np.ones(len(crossovers),dtype=int))

            #Update oligo numbers
            oligo_num[oligo_type] += 1


    def sequence_to_Tm(self,strand):
        #Count the number of A:T on the ends
        end_AT = (strand[0] == 'T' or strand[0] == 'A') + (strand[-1] == 'T' or strand[-1] == 'A')
        
        dH     = TUPLE_ENTHALPY['INIT'] + sum([TUPLE_ENTHALPY[strand[i:i+2]] for i in range(len(strand)-2)]) + end_AT*TUPLE_ENTHALPY['PENL']
        dS     = TUPLE_ENTROPY['INIT'] + sum([TUPLE_ENTROPY[strand[i:i+2]] for i in range(len(strand)-2)]) + end_AT*TUPLE_ENTROPY['PENL']
        
        #Equation from Santalucia paper
        Tm_1MNaCl = 1.0*dH/(self.R*np.log(self.STAP_CONC-0.5*self.SCAF_CONC)+dS/1000) - 273.15

        #Fraction of G/C
        fGC    = 1.0*sum([(c == 'G' or c == 'C') for c in strand])/len(strand)
        
        #Equation 16 from salt correction paper: Richard Owczarzy et al, Biochemistry, 2008 
        ln_Mg     = np.log(self.MAGNESIUM)
        Nbp       = len(strand) 
        
        a         = 3.92E-5
        b         = -9.11E-6
        c         = 6.26E-5
        d         = 1.42E-5
        e         =-4.82E-4
        f         = 5.25E-4
        g         = 8.31E-5

        Tm_Mg_inv = 1.0/Tm_1MNaCl + a + b*ln_Mg + fGC*(c+d*ln_Mg) + 1.0/(2*(Nbp-1))*(e+f*ln_Mg+g*ln_Mg**2)
        Tm_Mg     = 1.0/Tm_Mg_inv
        
        return Tm_Mg

    def calculate_strands_Tm(self):
        '''    
        Calculate strand Tm
        '''
        if self.scaffold_sequence == '':
            return

        self.strand_Tm = [[],[]]
        for oligo_type in range(len(self.strand_sequences)):
            for oligo_num in range(len(self.strand_sequences[oligo_type])):
                Tm_values = []
                for strand_num in range(len(self.strand_sequences[oligo_type][oligo_num])):
                    sequence = self.strand_sequences[oligo_type][oligo_num][strand_num]
                    Tm_values.append(self.sequence_to_Tm(sequence))
                self.strand_Tm.append(Tm_values)

    def assign_all_break(self):
        '''
        Breaks all crossovers  - Assign 0 for the break
        '''
        for crossover_name in self.crossover_names_to_nums.keys():
            self.crossover_names_to_states[crossover_name] = 0

    def update_strands(self):
        '''
        Update the strands based on the crossover configuration
        '''
        for crossover_name in self.crossover_names_to_states.keys():
            if self.crossover_names_to_states[crossover_name] == 0:
                crossover = [int(x) for x in crossover_name.split('-')]
                strand5p = self.get_strand(crossover[0],crossover[1])
                strand3p = self.get_strand(crossover[2],crossover[3])

                self.part.removeXover(strand5p, strand3p)

    def calculate_lookup_tables(self):
        '''
        Calculates the lookup tables for auto staple breaking
        '''

        #Lookup tables
        self.DISTANCE_MATRIX           = [[],[]]
        self.UPPER_BOUND_MATRIX        = [[],[]]
        self.LOWER_BOUND_MATRIX        = [[],[]]
        self.FINAL_CONSTRAINTS         = [[],[]]
        self.EDGE_WEIGHT_MATRIX        = [[],[]]

        #Broken oligos
        for i in range(len(self.strand_lengths[0])):
            strand_lengths             = self.strand_lengths[0][i]

            #Max length constraints
            max_length_7                = np.array(np.array(strand_lengths) >= 7 ,dtype=int) 
            max_length_14               = np.array(np.array(strand_lengths) >= 14,dtype=int)
            max_length_21               = np.array(np.array(strand_lengths) >= 21,dtype=int)
            max_length_28               = np.array(np.array(strand_lengths) >= 28,dtype=int)
            
            #Initialize the arrays
            self.DISTANCE_MATRIX[0].append([])
            self.LOWER_BOUND_MATRIX[0].append([])
            self.UPPER_BOUND_MATRIX[0].append([])
            self.FINAL_CONSTRAINTS[0].append([])
            self.EDGE_WEIGHT_MATRIX[0].append([])
            
            self.DISTANCE_MATRIX[0][i]           = np.zeros((len(strand_lengths),len(strand_lengths)))
            self.LOWER_BOUND_MATRIX[0][i]        = np.zeros((len(strand_lengths),len(strand_lengths)))
            self.UPPER_BOUND_MATRIX[0][i]        = np.zeros((len(strand_lengths),len(strand_lengths)))
            self.EDGE_WEIGHT_MATRIX[0][i]        = np.zeros((len(strand_lengths),len(strand_lengths)))

            for j in range(len(strand_lengths)):
                self.DISTANCE_MATRIX[0][i][j,j:]           = np.cumsum(strand_lengths[j:])
                self.EDGE_WEIGHT_MATRIX[0][i][j,j:]        = np.cumsum(max_length_7[j:]) + np.cumsum(max_length_14[j:]) + np.cumsum(max_length_21[j:]) + np.cumsum(max_length_28[j:])       

            self.DISTANCE_MATRIX[0][i]           = np.array(self.DISTANCE_MATRIX[0][i],dtype=int)
            self.LOWER_BOUND_MATRIX[0][i]        = np.array(self.DISTANCE_MATRIX[0][i] >= self.LOWER_BOUND,dtype=int)
            self.UPPER_BOUND_MATRIX[0][i]        = np.array(self.DISTANCE_MATRIX[0][i] <= self.UPPER_BOUND,dtype=int)
            self.FINAL_CONSTRAINTS[0][i]         = np.array(self.LOWER_BOUND_MATRIX[0][i]*self.UPPER_BOUND_MATRIX[0][i] > 0,dtype=int)
            self.EDGE_WEIGHT_MATRIX[0][i]        = np.array(self.EDGE_WEIGHT_MATRIX[0][i] > 0, dtype = int)

            #Apply the constraints on the edge weight matrix
            self.EDGE_WEIGHT_MATRIX[0][i]        = self.EDGE_WEIGHT_MATRIX[0][i]*self.FINAL_CONSTRAINTS[0][i]

            #Make all the 0 edges -INFINITY
            valid_row = np.nonzero(self.EDGE_WEIGHT_MATRIX[0][i]==0)[0]
            valid_col = np.nonzero(self.EDGE_WEIGHT_MATRIX[0][i]==0)[1]
            self.EDGE_WEIGHT_MATRIX[0][i][valid_row,valid_col] = -self.INFINITY
        
        #Circular oligos
        for i in range(len(self.strand_lengths[1])):
            strand_lengths        = self.strand_lengths[1][i] 
            double_strand_lengths = strand_lengths + strand_lengths

            #Max length constraints
            max_length_7                     = np.array(np.array(strand_lengths) >= 7,dtype=int) 
            max_length_14                    = np.array(np.array(strand_lengths) >= 14,dtype=int)
            max_length_21                    = np.array(np.array(strand_lengths) >= 21,dtype=int)
            max_length_28                    = np.array(np.array(strand_lengths) >= 28,dtype=int)  

            double_max_length_7              = np.hstack((max_length_7,max_length_7))
            double_max_length_14             = np.hstack((max_length_14,max_length_14))
            double_max_length_21             = np.hstack((max_length_21,max_length_21))
            double_max_length_28             = np.hstack((max_length_28,max_length_28))

            #Initialize the arrays
            self.DISTANCE_MATRIX[1].append([])
            self.LOWER_BOUND_MATRIX[1].append([])
            self.UPPER_BOUND_MATRIX[1].append([])
            self.FINAL_CONSTRAINTS[1].append([])
            self.EDGE_WEIGHT_MATRIX[1].append([])

            for j in range(len(strand_lengths)):
                #Distance matrix
                edge_distances = np.cumsum(double_strand_lengths[j:j+len(strand_lengths)])
                self.DISTANCE_MATRIX[1][i].append(np.hstack((edge_distances[-j:],edge_distances[:-j])))
                
                #Edge weight matrix edges
                edge_weights  = np.cumsum(double_max_length_7[j:j+len(max_length_7)]) + np.cumsum(double_max_length_14[j:j+len(max_length_14)]) + np.cumsum(double_max_length_21[j:j+len(max_length_21)]) + np.cumsum(double_max_length_28[j:j+len(max_length_28)])
                self.EDGE_WEIGHT_MATRIX[1][i].append(np.hstack((edge_weights[-j:],edge_weights[:-j])))

            self.DISTANCE_MATRIX[1][i]           = np.array(self.DISTANCE_MATRIX[1][i],dtype=int)
            self.LOWER_BOUND_MATRIX[1][i]        = np.array(self.DISTANCE_MATRIX[1][i] >= self.LOWER_BOUND,dtype=int)
            self.UPPER_BOUND_MATRIX[1][i]        = np.array(self.DISTANCE_MATRIX[1][i] <= self.UPPER_BOUND,dtype=int)
            self.FINAL_CONSTRAINTS[1][i]         = np.array(self.LOWER_BOUND_MATRIX[1][i]*self.UPPER_BOUND_MATRIX[1][i] > 0, dtype=int)
            self.EDGE_WEIGHT_MATRIX[1][i]        = np.array(self.EDGE_WEIGHT_MATRIX[1][i],dtype=int)

            #Apply the constraints on the edge weight matrix
            self.EDGE_WEIGHT_MATRIX[1][i]        = self.EDGE_WEIGHT_MATRIX[1][i]*self.FINAL_CONSTRAINTS[1][i]
            
            #Make all the 0 edges -INFINITY
            valid_row = np.nonzero(self.EDGE_WEIGHT_MATRIX[1][i]==0)[0]
            valid_col = np.nonzero(self.EDGE_WEIGHT_MATRIX[1][i]==0)[1]
            
            self.EDGE_WEIGHT_MATRIX[1][i][valid_row,valid_col] = -self.INFINITY

    def parse_crossover_constraints(self):
        '''
        Read crossover constraints and build a oligo connectivity map
        '''
        
        #Determine which oligo interacts with which oligo
        self.oligo_oligo_interaction  = {}
        self.crossover_constraints    = {}

        for crossover_name_1 in self.crossover_names_to_neighbors.keys():
            crossover_name_2 = self.crossover_names_to_neighbors[crossover_name_1]

            crossover_id_1   = self.crossover_names_to_nums[crossover_name_1]
            oligo_type_1,oligo_num_1,crossover_num_1 = crossover_id_1
            
            if not crossover_name_2 in self.crossover_names_to_nums.keys():
                self.crossover_constraints[(oligo_type_1,oligo_num_1,crossover_num_1)] = 1
                continue
            crossover_id_2 = self.crossover_names_to_nums[crossover_name_2]
            oligo_type_2,oligo_num_2,crossover_num_2 = crossover_id_2
            
            if not (oligo_type_1,oligo_num_1) in self.oligo_oligo_interaction.keys():
                self.oligo_oligo_interaction[(oligo_type_1,oligo_num_1)]  = {}
            if not (oligo_type_2,oligo_num_2) in self.oligo_oligo_interaction[(oligo_type_1,oligo_num_1)].keys():
                self.oligo_oligo_interaction[(oligo_type_1,oligo_num_1)][(oligo_type_2,oligo_num_2)] = []

            self.oligo_oligo_interaction[(oligo_type_1,oligo_num_1)][(oligo_type_2,oligo_num_2)] += [(crossover_num_1,crossover_num_2)]

    def shortest_path_to_result(self,path):
        '''
        Convert shortest path to result
        '''

        crossover_list = []

        #Retrieve the final score and crossover numbers
        crossover_i, crossover_j, score = path[-1]
        crossover_list          = [(crossover_j,score)] + crossover_list
        for new_crossover_i,new_crossover_j,new_score in path[-2::-1]:
            if new_crossover_j == crossover_i:
                crossover_list          = [(new_crossover_j,new_score)] + crossover_list
                crossover_i,crossover_j = new_crossover_i,new_crossover_j
        
        return crossover_list


    def find_local_solutions(self):
        '''
        Find local solutions
        '''
        self.local_solutions = [[],[]]

        for i in range(len(self.FINAL_CONSTRAINTS[0])):
            best_paths = self.k_shortest_path_broken(self.EDGE_WEIGHT_MATRIX[0][i],-1,self.EDGE_WEIGHT_MATRIX[0][i].shape[0],num_solutions=self.NUM_SOLUTIONS_PER_OLIGO)
            self.local_solutions[0].append(best_paths)

        for i in range(len(self.FINAL_CONSTRAINTS[1])):
            best_paths                  = []
            num_crossovers              = len(self.strand_lengths[1][i])
            #Determine the number of solutions for each subsolution
            num_solutions_per_crossover = ceil(1.0*self.NUM_SOLUTIONS_PER_OLIGO/num_crossovers)  
            for start in range(len(self.FINAL_CONSTRAINTS[1][i])):
                best_paths += self.k_shortest_path_circular(self.EDGE_WEIGHT_MATRIX[1][i],start,start,num_solutions=num_solutions_per_crossover)
            self.local_solutions[1].append(best_paths)
        
    def find_global_solutions(self):
        '''
        Find global solutions
        '''
        #Initialize the best solution lists
        
        self.best_global_solutions  = {}   #Dictionary of best solutions
        self.best_global_solution   = {}   #Dictionary of best solution

        #1. pick the best local solutions
        for main_key in self.pruned_local_solutions.keys():
            oligo_type,oligo_num = main_key
            for solution_num in self.pruned_local_solutions[main_key].keys():
                solution_list      = [(oligo_type,oligo_num,solution_num)]
                solution_keys      = [main_key]
                score_list         = [self.local_solutions[oligo_type][oligo_num][solution_num]['score']]
                current_solution   = self.pruned_local_solutions[main_key][solution_num]
                self.traverse_pruned_local_solutions(current_solution,solution_list,solution_keys,score_list,random_solution=False)
                global_key         = sorted(solution_keys,key=lambda x:1000*x[0]+x[1])
                global_key         = '_'.join([str(x[0])+'-'+str(x[1]) for x in global_key])
                if not global_key in self.best_global_solutions.keys():
                    self.best_global_solutions[global_key] = []
                penalty = self.check_solution_constraints(solution_list)
                self.best_global_solutions[global_key].append([solution_list,penalty,sum(score_list)])
        
        #2. pick random solutions
        for itr in range(self.NUM_GLOBAL_ITERATIONS):
            for main_key in self.pruned_local_solutions.keys():
                oligo_type,oligo_num = main_key
                for solution_num in self.pruned_local_solutions[main_key].keys():
                    solution_list      = [(oligo_type,oligo_num,solution_num)]
                    solution_keys      = [main_key]
                    score_list         = [self.local_solutions[oligo_type][oligo_num][solution_num]['score']]
                    current_solution   = self.pruned_local_solutions[main_key][solution_num]
                    self.traverse_pruned_local_solutions(current_solution,solution_list,solution_keys,score_list,random_solution=True)
                    global_key         = sorted(solution_keys,key=lambda x:1000*x[0]+x[1])
                    global_key         = '_'.join([str(x[0])+'-'+str(x[1]) for x in global_key])
                    penalty = self.check_solution_constraints(solution_list)
                    self.best_global_solutions[global_key].append([solution_list, penalty,sum(score_list)])
        
        #3. Pick the best solution for each cluster
        for global_key in self.best_global_solutions.keys():
            self.best_global_solutions[global_key] = sorted(self.best_global_solutions[global_key],key=lambda x:x[1]*1E10-x[2])
            self.best_global_solution[global_key]  = self.best_global_solutions[global_key][0]

    def print_global_solutions(self,num_solutions=10):
        '''
        Print the scores for the global solutions
        '''
        for global_key in self.best_global_solutions.keys():
            solution = 0
            oligo_names = global_key.replace('-','.')
            oligo_names = oligo_names.replace('_',':')
            print('Best solutions for domain with oligos {0}'.format(oligo_names))
            while solution < len(self.best_global_solutions[global_key]) and solution < num_solutions:
                penalty_score      = self.best_global_solutions[global_key][solution][1]
                optimization_score = self.best_global_solutions[global_key][solution][2]
                print('Solution {0}: penalty score= {1}, optimization score= {2}'.format(solution,penalty_score,optimization_score))
                solution+=1

    def check_solution_constraints(self,solution_list):
        '''
        Check the constraints for a solution
        '''
        crossover_name_list = []
        penalty           = 0
        for solution in solution_list:
            oligo_type,oligo_num,solution_num = solution
            solution_path  = self.local_solutions[oligo_type][oligo_num][solution_num]['path']
            crossover_path = [cross for cross,score in solution_path]
            
            #Broken oligo
            crossover_list = [crossover for crossover in crossover_path[1:-1]]
            
            #Circular oligo
            if oligo_type == 1:
                crossover_list = [crossover for crossover in crossover_path[1:]]

            #Crossover id list
            crossover_name_list += [self.crossover_nums_to_names[oligo_type][oligo_num][crossover] for crossover in crossover_list]
            
        for crossover_name in crossover_name_list:
            if self.crossover_names_to_neighbors[crossover_name] in crossover_name_list:
                penalty += 1

        return penalty


    def apply_best_solution(self):
        '''
        Apply best solution to crossover list
        '''
        for key in self.best_global_solution.keys():
            #Get the best solution results
            solution_list,penalty, score = self.best_global_solution[key]
            for solution in solution_list:
                oligo_type,oligo_num,solution_num = solution
                solution_path  = self.local_solutions[oligo_type][oligo_num][solution_num]['path']
                crossover_path = [cross for cross,score in solution_path]
                
                #Broken oligo
                crossover_list = [crossover for crossover in crossover_path[1:-1]]
                
                #Circular oligo
                if oligo_type == 1:
                    crossover_list = [crossover for crossover in crossover_path[1:]]
                
                #Update crossover list
                for crossover in crossover_list:
                    crossover_name = self.crossover_nums_to_names[oligo_type][oligo_num][crossover]
                    self.crossover_names_to_states[crossover_name] = 0

    def traverse_pruned_local_solutions(self,current_solution,solution_list,solution_keys,score_list,random_solution=False):
        '''
        Traverse final local solutions to combine local solutions
        '''
        for main_key in current_solution.keys():
            if not main_key in solution_keys and len(current_solution[main_key]):
                solution_key = 0
                
                #Find random solution
                if random_solution:
                    solution_key = random.randint(0,len(current_solution[main_key])-1)
                solution       = current_solution[main_key][solution_key]
                oligo_type,oligo_num,solution_num = solution
                next_solution  = self.pruned_local_solutions[main_key][solution_num]
                solution_list += [solution]
                solution_keys += [main_key]
                score_list    += [self.local_solutions[oligo_type][oligo_num][solution_num]['score']]
                self.traverse_pruned_local_solutions(next_solution,solution_list,solution_keys,score_list,random_solution)
    
    def check_local_solutions(self):
        '''
        Check all location solutions and 
        throw a warning if a solution doesn't exist for at least one oligo
        ''' 
        self.local_solutions_complete = True

        for oligo_type in range(len(self.local_solutions)):
            for oligo_num in range(len(self.local_solutions[oligo_type])):
                solutions = []
                for solution in range(len(self.local_solutions[oligo_type][oligo_num])):
                    solutions.append(self.local_solutions[oligo_type][oligo_num][solution]['valid'])
                if sum(solutions) == 0:
                    self.local_solutions_complete = False
                    return 0

    def apply_constraints_to_local_solutions(self):
        '''
        Find global solution from local solutions based on the neighborhood constraints
        '''

        #1. iterate over the hard constraints
        for key in self.crossover_constraints.keys():
            oligo_type,oligo_num,crossover = key
            oligo_solutions = self.local_solutions[oligo_type][oligo_num]
            
            #Iterate over all the solutions
            for solution in oligo_solutions:
                crossover_path = [cross for cross,score in solution['path']]
                if not crossover in crossover_path:
                    solution['satisfied'][0][key] = 1
                else:
                    solution['satisfied'][0][key] = 0
        
        #2. Check internal contraints for each oligo solution
        for oligo_type in range(len(self.local_solutions)):
            for oligo_num in range(len(self.local_solutions[oligo_type])):
                oligo_solutions  = self.local_solutions[oligo_type][oligo_num]
                pruned_solutions = [] 
                for solution in oligo_solutions:

                    crossover_path        = [cross for cross,score in solution['path']]
                    constraints_satisfied = True                    
                    for crossover_key_1 in self.oligo_oligo_interaction.keys():
                        oligo_type_1,oligo_num_1 = crossover_key_1
                        for crossover_key_2 in self.oligo_oligo_interaction[crossover_key_1].keys():
                            oligo_type_2,oligo_num_2 = crossover_key_2
                            for crossover_1,crossover_2 in self.oligo_oligo_interaction[crossover_key_1][crossover_key_2]:
                                if oligo_type == oligo_type_1 and oligo_type_1 == oligo_type_2 and oligo_num == oligo_num_1 and oligo_num_1 == oligo_num_2:
                                    if crossover_1 in crossover_path and crossover_2 in crossover_path:
                                        constraints_satisfied = False
                            
                    #Internally consistent solution
                    if constraints_satisfied:
                        pruned_solutions.append(solution)
                #Pruned solution list
                self.local_solutions[oligo_type][oligo_num] = pruned_solutions
        
        #3. Iterate over the neighborhood constraints
        for crossover_key_1 in self.oligo_oligo_interaction.keys():
            
            #Retrieve the oligo-1 info 
            oligo_type_1, oligo_num_1  = crossover_key_1
            oligo_solutions_1          = self.local_solutions[oligo_type_1][oligo_num_1]

            for crossover_key_2 in self.oligo_oligo_interaction[crossover_key_1].keys():
                #If the constraint is within the oligo, skip the constraint
                if crossover_key_1 == crossover_key_2:
                    continue

                #Retrieve the oligo-2 info
                oligo_type_2, oligo_num_2 = crossover_key_2
                oligo_solutions_2         = self.local_solutions[oligo_type_2][oligo_num_2]
                
                for crossover_1,crossover_2 in self.oligo_oligo_interaction[crossover_key_1][crossover_key_2]:
                    for i in range(len(oligo_solutions_1)):
                        #Retrieve the solution
                        solution_1 = oligo_solutions_1[i]
                        
                        #Initialize the constraint list
                        if not (oligo_type_2,oligo_num_2) in solution_1['satisfied'][1].keys():
                            solution_1['satisfied'][1][(oligo_type_2,oligo_num_2)] = {}
                        solution_1['satisfied'][1][(oligo_type_2,oligo_num_2)][crossover_2] = []
                        
                        #Build the crossover path
                        crossover_path_1 = [crossover for crossover,score in solution_1['path']]
                        for j in range(len(oligo_solutions_2)):
                            #Retrieve the solution
                            solution_2 = oligo_solutions_2[j]

                            #Build the crossover path
                            crossover_path_2 = [crossover for crossover,score in solution_2['path']]
                            
                            if not(crossover_1 in crossover_path_1 and crossover_2 in crossover_path_2):
                                solution_1['satisfied'][1][(oligo_type_2,oligo_num_2)][crossover_2].append((oligo_type_2,oligo_num_2,j))
        
        #4. Evaluate results
        for oligo_type in range(len(self.local_solutions)):
            for oligo_num in range(len(self.local_solutions[oligo_type])):
                for solution in range(len(self.local_solutions[oligo_type][oligo_num])):

                    #1. Initialize the valid property 
                    self.local_solutions[oligo_type][oligo_num][solution]['valid'] = 1
                    
                    #2. Check the hard constraints first
                    for key in self.local_solutions[oligo_type][oligo_num][solution]['satisfied'][0]:
                        self.local_solutions[oligo_type][oligo_num][solution]['valid'] *= (self.local_solutions[oligo_type][oligo_num][solution]['satisfied'][0][key] == 1)

                    #3. Check if there is at least one solution for each constraint
                    for key_1 in self.local_solutions[oligo_type][oligo_num][solution]['satisfied'][1].keys():
                        
                        combined_solutions = []
                        for key_2 in self.local_solutions[oligo_type][oligo_num][solution]['satisfied'][1][key_1]:
                            combined_solutions.append(self.local_solutions[oligo_type][oligo_num][solution]['satisfied'][1][key_1][key_2])
                        combined_solutions = reduce(lambda x,y:set(x).intersection(y),combined_solutions)
                        
                        #Keep the intersection list
                        self.local_solutions[oligo_type][oligo_num][solution]['satisfied'][1][key_1] = list(combined_solutions)
                        
                        #If the solution is empty, assign the solution invalid
                        if len(self.local_solutions[oligo_type][oligo_num][solution]['satisfied'][1][key_1]) == 0:
                            self.local_solutions[oligo_type][oligo_num][solution]['valid'] = 0
                        
        #5. Initialize the final solution list 
        self.pruned_local_solutions = {}
        for oligo_type in range(len(self.local_solutions)):
            for oligo_num in range(len(self.local_solutions[oligo_type])):
                for solution in range(len(self.local_solutions[oligo_type][oligo_num])):
                    if self.local_solutions[oligo_type][oligo_num][solution]['valid'] == 1:
                        #Initialize the dictionary
                        if not (oligo_type,oligo_num) in self.pruned_local_solutions.keys():
                            self.pruned_local_solutions[(oligo_type,oligo_num)] = {}

                        self.pruned_local_solutions[(oligo_type,oligo_num)][solution] = {}
                        
                        #Keep the best results
                        for key in self.local_solutions[oligo_type][oligo_num][solution]['satisfied'][1].keys():
                            self.pruned_local_solutions[(oligo_type,oligo_num)][solution][key] = []

        #6. Fill in the final solution list
        for main_key_1 in self.pruned_local_solutions.keys():
            oligo_type,oligo_num = main_key_1
            for solution in self.pruned_local_solutions[main_key_1].keys():
                for main_key_2 in self.pruned_local_solutions[main_key_1][solution].keys():
                    for possible_solution in self.local_solutions[oligo_type][oligo_num][solution]['satisfied'][1][main_key_2]:
                        possible_oligo_type,possible_oligo_num,possible_solution_num = possible_solution
                        if (possible_oligo_type,possible_oligo_num) in self.pruned_local_solutions.keys() and possible_solution_num in self.pruned_local_solutions[(possible_oligo_type,possible_oligo_num)].keys():
                            self.pruned_local_solutions[(oligo_type,oligo_num)][solution][main_key_2].append(possible_solution)

    def shortest_path_broken(self,matrix,start,final,initial_score=0):
        '''
        Apply shortest path between start and final crossovers
        '''
        best_scores = [initial_score]
        best_paths  = [(None,start,initial_score)]
        for column in range(start+1,final):
            
            scores = np.array(matrix[start+1:column+1,column])+np.array(best_scores)

            max_score = np.amax(scores)
            best_path = np.argmax(scores)+start+1
            best_scores += [max_score]
            best_paths  += [(best_path-1,column,max_score)]

        return best_paths

    def k_shortest_path_broken(self,matrix,start,final,num_solutions=10):
        '''
        Find the k-shortest paths for the broken oligos
        '''
        best_k          = []
        shortest_path   = self.shortest_path_broken(matrix,start,final)
        shortest_result = self.shortest_path_to_result(shortest_path) 
        
        #Add the shortest to best-k
        best_k.append({'path':shortest_result,'score':shortest_result[-1][1],'satisfied':[{},{}]})
        for k in range(num_solutions):
            potential_k = []
            for i in range(len(best_k[-1]['path'])-1):
                new_matrix           = np.copy(matrix)
                root_path            = best_k[-1]['path'][:i+1]
                from_cross,to_cross  = best_k[-1]['path'][i][0],best_k[-1]['path'][i+1][0]
                new_matrix[from_cross+1,to_cross] = -self.INFINITY
                #If root path exists in the other best paths remove the forward edge
                for r in range(len(best_k)-1):
                    if root_path == best_k[r]['path'][:i+1]:
                        from_cross_other,to_cross_other  = best_k[r]['path'][i][0],best_k[r]['path'][i+1][0]
                        new_matrix[from_cross_other+1,to_cross_other] = -self.INFINITY
                new_shortest_result  = root_path[:-1]+self.shortest_path_to_result(self.shortest_path_broken(new_matrix,from_cross,final,best_k[-1]['path'][i][1]))
                potential_k.append({'path':new_shortest_result,'score':new_shortest_result[-1][1],'satisfied':[{},{}]})

            #Sort potential candidates
            potential_k = (sorted(potential_k,key=lambda x:x['score'])[::-1])
            
            #Add the best solution 
            if potential_k[0]['path'] in [solution['path'] for solution in best_k] or potential_k[0]['score'] < 0:
                return best_k
            else:
                best_k.append(potential_k[0])

        return best_k

    def k_shortest_path_circular(self,matrix,start,final,num_solutions=10):
        '''
        Find the k-shortest paths for the circular oligos
        '''
        best_k          = []
        shortest_path   = self.shortest_path_circular(matrix,start,final)
        shortest_result = self.shortest_path_to_result(shortest_path)

        #Add the shortest to best-k
        if len(shortest_result) > 1 and shortest_result[0][0] == shortest_result[-1][0]:
            best_k.append({'path':shortest_result,'score':shortest_result[-1][1],'satisfied':[{},{}]})
        else:
            return []

        for k in range(num_solutions):
            potential_k = []
            for i in range(len(best_k[-1]['path'])-1):
                new_matrix           = np.copy(matrix)
                root_path            = best_k[-1]['path'][:i+1]
                from_cross,to_cross  = best_k[-1]['path'][i][0],best_k[-1]['path'][i+1][0]
                new_matrix[(from_cross+1)%matrix.shape[0],to_cross] = -self.INFINITY
                
                #If root path exists in the other best paths remove the forward edge
                for r in range(len(best_k)-1):
                    if root_path == best_k[r]['path'][:i+1]:
                        from_cross_other,to_cross_other  = best_k[r]['path'][i][0],best_k[r]['path'][i+1][0]
                        new_matrix[(from_cross_other+1)%matrix.shape[0],to_cross_other] = -self.INFINITY
                if from_cross < final:
                    new_shortest_result  = root_path[:-1]+self.shortest_path_to_result(self.shortest_path_broken(new_matrix,from_cross,final,best_k[-1]['path'][i][1]))
                else:
                    new_shortest_result  = root_path[:-1]+self.shortest_path_to_result(self.shortest_path_circular(new_matrix,from_cross,final,best_k[-1]['path'][i][1]))
                
                #If the solution is circular add to potential solution list
                if len(new_shortest_result) > 1 and new_shortest_result[0][0] == new_shortest_result[-1][0]:
                    potential_k.append({'path':new_shortest_result,'score':new_shortest_result[-1][1],'satisfied':[{},{}]})
            #Sort potential candidates
            potential_k = (sorted(potential_k,key=lambda x:x['score'])[::-1])
            
            #Add the best solution 
            if len(potential_k) == 0 or potential_k[0]['path'] in [solution['path'] for solution in best_k] or potential_k[0]['score'] < 0:
                return best_k
            else:
                best_k.append(potential_k[0])

        return best_k

    def shortest_path_circular(self,matrix,start,final,initial_score=0):
        '''
        Apply shortest path between start and final crossovers
        '''
        best_scores   = [initial_score]
        best_paths    = [(None,start,initial_score)]
        matrix_length = matrix.shape[0]
        for column in range(start+1,final+matrix_length+1):
            
            if column < matrix_length:
                previous_scores = matrix[start+1:column+1,column%matrix_length]
            else:
                previous_scores = np.hstack((matrix[start+1:matrix_length,column%matrix_length],matrix[:column-matrix_length+1,column%matrix_length]))
            scores = np.array(previous_scores)+np.array(best_scores)

            max_score    = np.amax(scores)
            best_path    = (np.argmax(scores)+start+1)%matrix_length
            best_scores += [max_score]
            best_paths  += [((best_path-1)%matrix_length,column%matrix_length,max_score)]

        return best_paths

    def break_21_strands(self):
        '''
        Break strands that are exactly 21 base long flanked by crossovers
        '''

    def color_by_penalty(self,color_map='viridis'):
        '''
        Color oligos by the penalty score
        '''
        
        flag_color  = matplotlib.colors.to_hex('red')
        valid_color = matplotlib.colors.to_hex('gray')

    def color_by_segments(self,color_map='jet'):
        '''
        Color by the segments
        '''

        #Number of segments
        num_segments = len(self.best_global_solution.keys())
        
        if num_segments == 0:
            return None

        #Prepare the color map
        cmap      = cm.get_cmap(color_map, num_segments)
        hexcolors = [] 
        for i in range(cmap.N):
            rgb = cmap(i)[:3]       # will return rgba, we take only first 3 so we get rgb
            hexcolors.append(matplotlib.colors.rgb2hex(rgb))
        
        #Retrieve the segment info
        segments = list(self.best_global_solution.keys())

        #Crossover colors
        crossover_colors = []

        for segment_id in range(len(segments)):
            #Oligo color
            oligo_color = hexcolors[segment_id]

            #Get the segment info
            segment = segments[segment_id]

            #Retrieve oligo information from the solutions 
            oligo_tuplets_str = [ oligo.split('-') for oligo in segment.split('_')]
            oligo_tuplets_int = [ (int(oligo_type),int(oligo_num)) for oligo_type,oligo_num in oligo_tuplets_str]
            
            #Retrieve all the crossovers
            crossover_names = []
            for oligo_type,oligo_num in oligo_tuplets_int:
                crossover_names += [self.crossover_nums_to_names[oligo_type][oligo_num][crossover_num] for crossover_num in range(len(self.crossover_nums_to_names[oligo_type][oligo_num]))]

            #Get the crossover neighbors
            neighbor_names = []
            for crossover_name in crossover_names:
                neighbor_name = self.crossover_names_to_neighbors[crossover_name]
                if not neighbor_name in crossover_names:
                    neighbor_names += [self.crossover_names_to_neighbors[crossover_name]]
            
            #Combine the crossover lists
            crossover_names += neighbor_names

            #Crossover colors
            crossover_colors += [(crossover_name,oligo_color) for crossover_name in crossover_names]
            
            #Assign the colors
            for crossover_name,oligo_color in crossover_colors:
                crossover = [int(x) for x in crossover_name.split('-')]

                #Get the strands
                strand_1 = self.get_strand(crossover[0],crossover[1])
                strand_2 = self.get_strand(crossover[2],crossover[3])

                #Apply colors
                strand_1.oligo().applyColor(oligo_color)
                strand_2.oligo().applyColor(oligo_color)

    def color_by_oligo_length(self):
        '''
        Color by the length
        '''
        cmap      = cm.get_cmap(color_map, 60)
        hexcolors = [] 
        for i in range(cmap.N):
            rgb = cmap(i)[:3]       # will return rgba, we take only first 3 so we get rgb
            hexcolors.append(matplotlib.colors.rgb2hex(rgb))

    def color_by_max_strand_length(self):
        '''
        Color by the length
        '''
        cmap      = cm.get_cmap(color_map, 60)
        hexcolors = [] 
        for i in range(cmap.N):
            rgb = cmap(i)[:3]       # will return rgba, we take only first 3 so we get rgb
            hexcolors.append(matplotlib.colors.rgb2hex(rgb))

    def get_strand(self,idnum,idx):
        '''
        Get strand from idnum and idx
        '''

        direction = (idnum+1)%2
        ss        = self.part.getStrandSets(idnum)[direction]
        strand    = ss.getStrand(idx)

        return strand

    def check_connection(self,crossover):
        '''
        Checks if the crossover is broken or not
        Returns:
            1 if connected 
            0 if broken
        '''

        strand1 = self.get_strand(crossover[0],crossover[1])
        strand2 = self.get_strand(crossover[2],crossover[3])

        if strand1 != None:
            return int(strand1.connection3p() == strand2)
        else:
            return 0

    def read_json(self,json_file):
        '''
        Read Cadnano json file
        '''
        if os.path.isfile(json_file):
            app = cadnano.app()

            self.json_input = json_file 
            self.doc        = app.document = Document()

            self.doc.readFile(json_file)
            
            self.part      = self.doc.activePart()
        else:
            sys.exit(-1)

    def read_sequence(self,txt_file):
        '''
        Read sequence file
        '''
        if os.path.isfile(txt_file):
            self.sequence_file     = txt_file
            f = open(txt_file) 
            self.scaffold_sequence = ''.join([line.strip() for line in f.readlines()])
            f.close()
        else:
            sys.exit(-1)

    def apply_sequence(self):
        '''
        Apply the sequence
        '''
        if len(self.scaffold_sequence) > 0:
            self.scaffold.applySequence(self.scaffold_sequence)

    def auto_break_rule_3(self):
        '''
        Auto staple breaking with Rule 3 using google-ortools
        '''

        path, ext = os.path.splitext(self.json_input)
        self.json_output_autobreak                = path+'_autobreak_rule_3'+ext

        if self.BREAK_LONG_STRANDS:
            self.json_output_longstrands_broken   = path+'_longstrands_broken'+ext
            self.parse_oligos()
            self.break_broken_long_strands()
            self.break_circular_long_strands()
            self.doc.writeToFile(self.json_output_longstrands_broken,legacy=True)
            print('Long strands are broken')

        self.parse_oligos()
        print('Oligos are imported')

        self.calculate_strands_Tm()

        self.calculate_lookup_tables()
        print('Look up tables are prepared')
        self.parse_crossover_constraints()

        self.find_local_solutions()
        print('Local solutions found')

        self.apply_constraints_to_local_solutions()
        print('Constraints applied to local solutions')

        self.check_local_solutions()
        print('Local solutions are checked')

        if not self.local_solutions_complete:
            print('Warning: Solutions are not found for some oligos. Increase the number of solutions per oligo')
        
        self.find_global_solutions()
        print('Global solutions found')
        
        self.print_global_solutions()

        self.apply_best_solution()
        self.update_strands()
        print('Best solution is applied to DNA structure')

        self.color_by_segments()
        self.doc.writeToFile(self.json_output_autobreak,legacy=True)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i",   "--input",         type=str,            help="Cadnano json file" )
    parser.add_argument("-s",   "--sequence",      type=str,            help="Sequence file in txt", default=None )
    
    parser.add_argument("-bl",   "--breaklong",     action='store_true', help="Break long (14,21,28 base long strands)")
    parser.add_argument("-ubl",  "--upperboundlen", type=int,            help="Oligo length upper bound",default=60 )
    parser.add_argument("-lbl",  "--lowerboundlen", type=int,            help="Oligo length lower bound",default=23 )

    parser.add_argument("-ns",   "--nsolsperoligo" , type=int,            help="Number of solutions per oligo",default=100 )
    parser.add_argument("-ni",   "--nitr"          , type=int,            help="Number of global iterations",default=50 )

    args = parser.parse_args()
    
    #Assign the parameters
    input_filename          = args.input
    sequence_filename       = args.sequence 
    break_long_strands      = args.breaklong
    num_solutions_per_oligo = args.nsolsperoligo 
    num_global_iterations   = args.nitr
    upper_bound_length      = args.upperboundlen
    lower_bound_length      = args.lowerboundlen

    #Check if input file exists
    if not os.path.isfile(input_filename):
        sys.exit('Input file does not exist!')

    new_autobreak = AutoBreak()
    new_autobreak.read_json(input_filename)
    
    if not sequence_filename == None and os.path.isfile(sequence_filename):
        new_autobreak.read_sequence(sequence_filename)

    #Assign the optimization parameters
    new_autobreak.BREAK_LONG_STRANDS      = break_long_strands
    new_autobreak.NUM_SOLUTIONS_PER_OLIGO = num_solutions_per_oligo
    new_autobreak.NUM_GLOBAL_ITERATIONS   = num_global_iterations
    new_autobreak.LOWER_BOUND             = lower_bound_length
    new_autobreak.UPPER_BOUND             = upper_bound_length

    #Apply Rule 3
    new_autobreak.auto_break_rule_3()

if __name__ == "__main__":
  main()








