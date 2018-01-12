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

from cadnano.document import Document
from functools import reduce

class Sequence:
    def __init__(self):
        '''
        DNA sequence class
        '''
        self.dna           = None
        self.next_sequence = None

class BreakSolution:
    def __init__(self):
        '''
        Break solution for an oligo
        '''
        
        self.oligo = None

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
        self.fwd_breaks  = np.array([+2])       #Fwd break locations along a strand
        self.rev_breaks  = np.array([-4,-1])    #Rev break locations along a strand

        #Break rule
        self.brake_rules = ['crossover','3']
        self.LONG_STRAND_LENGTH = 21
        self.LONG_STRAND_STEP   = 7

    def break_rule(self,rule_names=['crossover']):
        self.break_rule = rule_names
        
        #initialize break location
        self.fwd_breaks = []
        self.rev_breaks = []

        #1. Rule number 1 - Crossover only
        if 'crossover' in self.break_rule:
            self.rev_breaks.append(-1)
        
        #2. Rule number 2 - Only 3 base away from crossovers
        if '3' in self.break_rule:
            self.fwd_breaks.append(+2)
            self.rev_breaks.append(-4)
        
        #3. Rule number 3 - break long strands
        if 'long' in self.break_rule:
            if self.length >= self.LONG_STRAND_LENGTH:
                
                #Prepare the fwd and rev break locations
                fwd_breaks = list(range( self.LONG_STRAND_STEP,  self.length-self.LONG_STRAND_STEP, self.LONG_STRAND_STEP))
                rev_breaks = list(range(-self.LONG_STRAND_STEP, -self.length+self.LONG_STRAND_STEP,-self.LONG_STRAND_STEP))
                
                #Add break locations to list
                self.fwd_breaks += fwd_breaks
                self.rev_breaks += rev_breaks

        #4. For final strand, add the last position(-1)
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
        self.strands    = None
        self.type       = 'broken'
        self.crossovers = []

class Crossover:
    def __init__(self):
        '''
        Crossover class
        '''

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
        self.staples       = None
        self.scaffolds     = None
        self.idnums        = None 

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

        self.oligos[oligo_type].append(new_oligo)
        
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

            #Prepare the break points
            new_strand.break_rule()

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
        
        for oligo in self.oligos['staple']:
            #Initialize oligo breaks
            oligo.breaks = []

            #Assign current strand
            current_strand = oligo.null_strand.next_strand

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

            #Add null oligo to break list
            oligo.breaks.append(oligo.null_break)

            #If the oligo is not circular make the first break node start node
            if not oligo.circular:
                oligo.start_break = oligo.null_break

            #Iterate over each strand
            while current_strand:
                
                #Iterate through all positions
                for break_position in current_strand.all_breaks:

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

                    #Assign sequence to break object
                    new_break.sequence    = current_strand.sequences[sequence_id]

                    #Assign strand to new break
                    new_break.strand      = current_strand

                    #Check if the break is at a cross-over location
                    if new_break.key in self.crossovers:
                        new_break.crossover = self.crossovers[new_break.key]

                    #Assign to previous break
                    previous_break.next_break = new_break

                    #Update previous break
                    previous_break = new_break

                    #Add break to break list
                    oligo.breaks.append(new_break)

                #Update current strand
                current_strand   = current_strand.next_strand

            #If oligo is circular, connect final node to null break's next break
            if oligo.circular:
                previous_break.next_break = oligo.null_break.next_break

                #Remove the first break node from list
                oligo.breaks.pop(0)

                #Assign final break to null break
                oligo.null_break  = previous_break
            else:
                #If the oligo is not circular make the final break node final node
                oligo.final_break = previous_break

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

    def read_sequence(self,sequence_file):
        '''
        Read sequence file
        '''
        if os.path.isfile(sequence_file):
            self.sequence_file     = sequence_file
            
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
            sys.exit(-1)

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
        self.UPPER_BOUND                     = 60
        self.LOWER_BOUND                     = 21
        self.BREAK_LONG_STRANDS              = False

        #Local and global solutions
        self.NUM_OLIGO_SOLUTIONS             = 100
        self.NUM_GLOBAL_SOLUTIONS            = 50

        #CONSTANTS
        self.INFINITY                        = 1E6

    def initialize(self):
        '''
        Initialize the connectivity maps
        '''
        
        for oligo in self.origami.oligos['staple']:
            
            #Visit each break object
            for current_break in oligo.breaks:
                #Initialize the break edges
                current_break.break_edges = []

                #Get next break
                next_break = current_break.next_break

                #Iterate over the breaks
                while next_break:
                    
                    #Determine break to break distance
                    break_break_distance = (next_break.distance - current_break.distance)%oligo.length + oligo.length*(next_break==current_break)
                    
                    if break_break_distance >= self.LOWER_BOUND and break_break_distance <= self.UPPER_BOUND:
                        #Create break
                        new_edge = BreakEdge()
                        
                        #Assign edge length
                        new_edge.edge_length = break_break_distance
                        
                        #Make the connection
                        new_edge.make_connection(current_break,next_break)

                        
                        #Add directed edge to current break's edges
                        current_break.break_edges.append(new_edge)

                    if break_break_distance > self.UPPER_BOUND or next_break==current_break:
                        break

                    next_break = next_break.next_break

class BreakEdge:
    def __init__(self):
        '''
        Break edge class
        '''
        self.current_break = None 
        self.next_break    = None
        self.edge_weight   = None
        self.edge_length   = None
        self.edge_Tm       = None
        
        #state parameters
        self.active        = True

        #Length parameters
        self.edge_has14    = None
        self.edge_num14    = None
        
        #Tm parameters
        self.edge_has60    = None
        self.edge_num60    = None

        self.edge_has55    = None
        self.edge_num55    = None

        self.edge_has50    = None
        self.edge_num50    = None

        self.sequence_list = None

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
        self.edge_num14  = np.sum(self.length_list >= 14)  
        self.edge_has14  = int(self.edge_num14 > 0)

        self.edge_num10  = np.sum(self.length_list >= 10)  
        self.edge_has10  = int(self.edge_num10 > 0)

        #Tm parameters
        self.edge_num60  = np.sum(self.Tm_list >= 60)
        self.edge_has60  = np.sum(self.edge_num60 > 0)

        self.edge_num55  = np.sum(self.Tm_list >= 55)
        self.edge_has55  = np.sum(self.edge_num55 > 0)

        self.edge_num50  = np.sum(self.Tm_list >= 50)
        self.edge_has50  = np.sum(self.edge_num50 > 0)

class BreakNode:
    def __init__(self):
        '''
        Break node class
        '''
        self.crossover  = None
        self.next_break = None
        self.type       = None
        self.sequence   = None
        self.crossover  = None

        #State parameters
        self.active     = True
        self.state      = False  #True for break, false for no-break

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i",   "--input",         type=str,            help="Cadnano json file" )
    parser.add_argument("-s",   "--sequence",      type=str,            help="Sequence file in txt", default=None )

    args = parser.parse_args()
    
    #Assign the parameters
    input_filename          = args.input
    sequence_filename       = args.sequence 

    #Check if input file exists
    if not os.path.isfile(input_filename):
        sys.exit('Input file does not exist!')

    #Create Origami object
    new_origami = Origami()

    #Initialize origami object
    new_origami.initialize(input_filename)

    #Get oligos
    new_origami.get_oligos()

    #Reset oligos list
    new_origami.reset_oligos()

    #Read sequence file
    if not sequence_filename == None and os.path.isfile(sequence_filename):
        new_origami.read_sequence(sequence_filename)

    #Read scaffolds and staples
    new_origami.read_staples()

    #Generate crossovers
    new_origami.generate_crossovers()
    
    #Generate sequences
    new_origami.generate_sequences()

    #Generate break points
    new_origami.generate_break_points()

    #Create Autobreak object and make assignments
    new_autobreak = AutoBreak()
    new_origami.autobreak = new_autobreak
    new_autobreak.origami = new_origami

    #Make break-break graph
    new_autobreak.initialize()


if __name__ == "__main__":
  main()








