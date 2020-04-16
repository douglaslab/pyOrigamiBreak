#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2018-05-21 13:40:55
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import random
import os
import matplotlib
import math
import cadnano
import numpy as np
import sys

from tqdm import tqdm
from matplotlib import cm
from cadnano.document import Document

import ctk_autobreak as autobreak
import utilities

matplotlib.use('TkAgg')


class Sequence:
    def __init__(self):
        '''
        DNA sequence class
        '''
        self.dna           = None
        self.next_sequence = None

    def assign_scaffold_positions(self):
        '''
        Assign scaffold positions
        '''
        self.scaffoldPos = []
        for idx in range(self.idx5p, self.idx3p+self.direction, self.direction):
            # Get scaffold position
            scaffold_positions = self.origami.get_scaffold_positions(self.idNum, idx)

            # Add new positions to list
            self.scaffoldPos.extend(scaffold_positions)


class OligoGroup:
    def __init__(self):
        '''
        Oligo group connected to each other through break point neighborship
        '''
        self.key                   = None
        self.oligos                = None
        self.breaks                = None
        self.group_solutions       = None
        self.best_score_solution   = None
        self.best_penalty_solution = None

    def remove_incomplete_solutions(self):
        '''
        Remove incomplete group solutions
        '''
        self.group_solutions = list(filter(lambda x: x.complete, self.group_solutions))

    def shuffle_oligos(self):
        '''
        Shuffle oligos
        '''
        random.shuffle(self.oligos)

    def sort_oligos_by_length(self, reverse=True):
        '''
        Sort oligos by length
        '''
        self.oligos.sort(key=lambda x: x.length, reverse=reverse)

    def reset_temp_neighbor_constraints(self):
        '''
        Reset temporary neighbor constraints
        '''
        for oligo in self.oligos:
            oligo.reset_temp_neighbor_constraints()

    def combine_oligo_solutions(self, num_global_solutions=200):
        '''
        Create random group solutions
        '''
        # Initialize group solutions
        self.group_solutions = []

        for i in range(num_global_solutions):
            # Initialize group break solution
            new_group_solution = autobreak.GroupBreaksolution()
            new_group_solution.break_solutions = {}
            new_group_solution.origami = self.origami

            for oligo in self.oligos:
                # Randomly pick one solution for an oligo
                if len(oligo.break_solutions) > 0:
                    new_group_solution.break_solutions[oligo.key]  = random.choice(oligo.break_solutions)
                else:
                    new_group_solution.break_solutions[oligo.key] = None

            # Calculate the penalties for each group solution
            new_group_solution.calculate_penalty()

            # Add solution to list
            self.group_solutions.append(new_group_solution)

    def create_stepwise_oligo_solutions(self, num_oligo_solutions=100, num_global_solutions=500,
                                        pick_method='random', shuffle_oligos=True, verbose=False):
        '''
        Create stepwise oligo solutions
        '''

        # Initialize group solutions
        self.group_solutions = []

        # Number of solutions
        for i in tqdm(range(num_global_solutions), desc='Global loop', leave=False,
                      dynamic_ncols=True, bar_format='{l_bar}{bar}', 
                      file=self.origami.tqdm_output_file):

            # Reset temporary neighbor constraints
            self.reset_temp_neighbor_constraints()

            # Initialize group break solution
            new_group_solution = autobreak.GroupBreaksolution()
            new_group_solution.break_solutions = {}
            new_group_solution.origami = self.origami

            # Shuffle oligos
            if shuffle_oligos:
                self.shuffle_oligos()

            # Iterate over every oligo
            for oligo in tqdm(self.oligos, desc='Oligo loop', leave=False,
                              dynamic_ncols=True, bar_format='{l_bar}{bar}',
                              file=self.origami.tqdm_output_file):

                # If oligo has dont break flag, skip
                if oligo.dont_break:
                    continue

                # 1. Create shortest paths
                oligo.generate_shortest_paths(num_oligo_solutions, verbose=verbose)

                # 2. Remove penalized solutions
                oligo.remove_penalized_solutions()

                # 3. Pick a solution
                chosen_solution = oligo.pick_break_solution(pick_method)

                # 4. Assign solution
                new_group_solution.break_solutions[oligo.key] = chosen_solution

                # 5. Apply temporary neighbor constraints
                if chosen_solution:
                    chosen_solution.apply_temp_neighbor_constraints()

            # Calculate the penalties for each group solution
            new_group_solution.calculate_penalty()

            # Print new group solution
            if verbose:
                new_group_solution.print_solution()

            # Add solution to list
            self.group_solutions.append(new_group_solution)

    def sort_solutions(self, filter_incomplete=True):
        '''
        Sort solutions based on the penalty score
        '''
        # 1. Filter incomplete solutions
        if filter_incomplete:
            self.group_solutions = list(filter(lambda x: x.complete, self.group_solutions))

        # 2. Sort based on total score
        self.group_solutions.sort(key=lambda x: x.total_score, reverse=True)

        # 3. Assign best total score solution
        if len(self.group_solutions) > 0:
            self.best_score_solution = self.group_solutions[0]
        else:
            self.best_score_solution = None

        # 4. Sort based on total penalty
        self.group_solutions.sort(key=lambda x: x.total_penalty)

        # 5. Assign best penalty solution
        if len(self.group_solutions) > 0:
            self.best_penalty_solution = self.group_solutions[0]
        else:
            self.best_penalty_solution = None

    def print_solutions(self):
        '''
        Print solutions
        '''
        for group_solution in self.group_solutions:
            group_solution.print_solution()


class Nucleotide:
    def __init__(self):
        '''
        Nucleotide class
        '''
        self.vh        = None
        self.idx       = None
        self.direction = None
        self.key       = None
        self.dsDNA     = False

        self.next_nucleotide     = None
        self.previous_nucleotide = None


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

        # Possible break locations
        self.fwd_breaks  = None    # Fwd break locations along a strand
        self.rev_breaks  = None    # Rev break locations along a strand

        # Break rule
        self.break_rule = ['xstap', 'all3']
        self.LONG_STRAND_LENGTH = 21
        self.LONG_STRAND_STEP   = 7

    def get_inserts(self, idx_a, idx_b):
        '''
        Get inserts between two idx values on a strand
        '''
        idx_low, idx_high = (idx_a, idx_b) if idx_b > idx_a else (idx_b, idx_a)

        # Get insertion length between two idx
        return self.cadnano_strand.insertionLengthBetweenIdxs(idx_low, idx_high)

    def apply_break_rule(self):
        # Get the break rule
        self.break_rule = self.origami.break_rule

        # initialize break location
        self.fwd_breaks = []
        self.rev_breaks = []
        self.all_breaks = []

        # 1. Rule number 1 - Crossover only
        if 'xstap' in self.break_rule:
            self.rev_breaks.append(-1)

        # 2. Rule number 2 - Only 3 base away from 5p crossover
        if '3f' in self.break_rule:
            self.fwd_breaks.append(+2)

        # 3. Rule number 3 - Only 3 base away from 3p crossover
        if '3r' in self.break_rule:
            self.rev_breaks.append(-4)

        # 4. Rule number 4 - Only 3 base away from crossovers
        if '3' in self.break_rule:
            self.fwd_breaks.append(+2)
            self.rev_breaks.append(-4)

        # 5. Rule number 5 - break long strands
        if 'long' in self.break_rule:
            if self.length >= self.LONG_STRAND_LENGTH:
                # Prepare the fwd and rev break locations
                fwd_breaks = list(range(self.LONG_STRAND_STEP-1, self.length-self.LONG_STRAND_STEP,
                                  self.LONG_STRAND_STEP))
                rev_breaks = list(range(-self.LONG_STRAND_STEP-1, -self.length+self.LONG_STRAND_STEP,
                                  -self.LONG_STRAND_STEP))

                # Add break locations to list
                self.fwd_breaks += fwd_breaks
                self.rev_breaks += rev_breaks

        # 6. For final strand, add the last position(-1)
        if self.final_strand:
            self.rev_breaks.append(-1)

        # 7. Check if the rule is all2
        if 'all2' in self.break_rule:
            self.fwd_breaks.extend(np.arange(1, int(self.length)-2))

        # 8. Check if the rule is all3
        if 'all3' in self.break_rule:
            self.fwd_breaks.extend(np.arange(2, int(self.length)-3))

        # Make the breaks array and sort them
        self.fwd_breaks = np.array(sorted(self.fwd_breaks), dtype=int)
        self.rev_breaks = np.array(sorted(self.rev_breaks), dtype=int) + int(self.length)

        # All values need to be between 0 and length
        self.fwd_breaks = list(filter(lambda x: x >= 0 and x < self.length, self.fwd_breaks))
        self.rev_breaks = list(filter(lambda x: x >= 0 and x < self.length, self.rev_breaks))

        # Convert the lists to int array
        self.fwd_breaks = np.array(sorted(self.fwd_breaks), dtype=int)
        self.rev_breaks = np.array(sorted(self.rev_breaks), dtype=int)

        # Combine the two arrays
        self.all_breaks = np.sort(np.unique(np.hstack((self.fwd_breaks, self.rev_breaks))))

        # Adjust break positions for the inserts and skips
        self.adjust_break_positions()

    def adjust_break_positions(self):
        '''
        Adjust relative break positions taking insers/skips into account
        '''
        self.all_breaks_adjusted = []

        for current_break in self.all_breaks:
            # Get idx for break
            break_idx3p = self.idx5p + self.direction*current_break

            # Get number of inserts/skips between 5p and the break position
            num_inserts = self.get_inserts(self.idx5p, break_idx3p)

            # Add adjusted break to list
            self.all_breaks_adjusted.append(current_break+num_inserts)


class Oligo:
    def __init__(self):
        '''
        Oligo class
        '''
        self.cadnano_oligo       = None
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

        # Break solutions
        self.break_solutions       = []
        self.chosen_solution       = None
        self.best_score_solution   = None
        self.best_penalty_solution = None

    def break_in_half(self):
        '''
        Break the oligo in half
        '''

        # Get the median break
        median_break_num  = round(0.5*len(self.breaks))
        self.median_break = self.breaks[median_break_num]

        # Get the first break in the list
        start_break       = self.breaks[0]

        # If the oligo is circular break the first and middle breaks
        if self.circular:
            start_break.break_cadnano()
            self.median_break.break_cadnano()
        else:
            self.median_break.break_cadnano()

    def pick_break_solution(self, method='best'):
        '''
        Pick a break solution
        '''

        # Sort solutions by score
        self.sort_solutions_by_score()

        if len(self.break_solutions) == 0:
            return None

        if method == 'best':
            return self.break_solutions[0]
        else:
            return random.choice(self.break_solutions)

    def sort_solutions_by_score(self):
        '''
        Sort solutions by score
        '''
        self.break_solutions.sort(key=lambda x: x.score, reverse=True)

    def get_initial_score(self):
        '''
        Get the initial score before the oligo is broken
        '''
        # Create break
        new_edge = autobreak.BreakEdge()

        # Assign origami
        new_edge.origami = self.origami

        # Assign autobreak
        new_edge.autobreak = self.origami.autobreak

        # Make the connection
        if self.circular:
            start_break = self.breaks[0]
            final_break = self.breaks[0]
        else:
            start_break = self.start_break
            final_break = self.final_break

        # Get break to break distance
        break_distance  = start_break.get_break_distance(final_break)

        # Assign edge length
        new_edge.edge_length = break_distance

        # Make the connection
        new_edge.make_connection(start_break, final_break)
        self.end_to_end_edge = new_edge

        # Set edge weight
        new_edge.edge_weight = self.origami.autobreak.optimize(new_edge)

        # Assign the score
        self.initial_score = new_edge.edge_weight
        self.folding_prob  = new_edge.edge_prob
        self.Tf            = new_edge.edge_Tf
        self.dsDNA_length  = sum(new_edge.dsDNA_length_list)

    def color_by_folding_prob(self, color_map='bwr'):
        '''
        Color by the segments
        '''

        # Prepare the color map
        cmap      = cm.get_cmap(color_map, 1000)

        # Get color index
        cmap_index = int(1000*self.folding_prob[0])

        # Get RGB value
        self.rgb = cmap(cmap_index)[:3]

        # Get hex-color
        self.hexcolor = matplotlib.colors.rgb2hex(self.rgb)

        # Apply the color
        self.cadnano_oligo.applyColor(self.hexcolor)

    def color_by_Tf(self, color_map='bwr', min_Tf=30, max_Tf=50):
        '''
        Color by the segments
        '''

        # Prepare the color map
        cmap      = cm.get_cmap(color_map, max_Tf-min_Tf+1)

        # Get color index
        if self.Tf <= min_Tf:
            cmap_index = 0
        elif self.Tf >= max_Tf:
            cmap_index = max_Tf - min_Tf
        else:
            cmap_index = int(self.Tf-min_Tf)

        # Get RGB value
        self.rgb = cmap(cmap_index)[:3]

        # Get hex-color
        self.hexcolor = matplotlib.colors.rgb2hex(self.rgb)

        # Apply the color
        self.cadnano_oligo.applyColor(self.hexcolor)

    def remove_penalized_solutions(self):
        '''
        Remove solutions with non-zero self penalty score
        '''
        self.break_solutions = list(filter(lambda x: x.self_penalty == 0, self.break_solutions))

    def keep_best_break_solutions(self):
        '''
        Keep best break solutions
        '''
        # If the solution list is empty return
        if len(self.break_solutions) == 0:
            return

        # Sort the list based on scores
        self.break_solutions.sort(key=lambda x: x.score, reverse=True)

        # Get the maximum score
        self.max_score = self.break_solutions[0].score

        # Remove the solutions that have less score than the max score
        self.break_solutions = list(filter(lambda x: x.score == self.max_score, self.break_solutions))

    def reset_temp_neighbor_constraints(self):
        '''
        Reset temporary neighbor constraints
        '''
        for current_break in self.breaks:
            current_break.dont_break_temp = False

    def reset_break_order_ids(self, start_break, final_break):
        '''
        Renumber break ids for shortest path algorithm
        '''
        current_break = start_break

        # Break id counter
        id_counter = 0

        # Make the id of current break 0
        current_break.order_id = id_counter

        # Get the next break
        current_break = current_break.next_break

        # Iterate over each break
        while current_break:
            # Update id counter
            id_counter += 1

            # Update break id
            current_break.order_id = id_counter

            # If we reach the start break quit
            if current_break == final_break:

                # If start and final breaks are the same, make first break number 0
                if start_break == final_break:
                    current_break.order_id = 0
                break

            # Update current breaks
            current_break = current_break.next_break

    def generate_shortest_paths(self, num_solutions=1, verbose=False):
        '''
        Get the shortest paths for the oligo if only it allowed to break it
        '''

        # Get k-select parameter
        k_select = self.origami.autobreak.k_select

        # Show oligo being processed - use tdqm
        if verbose:
            tqdm.write('Processing oligo:%-15s Number of breaks:%-3d' % (self.key, len(self.breaks)),
                        file=self.origami.tqdm_output_file)

        if self.dont_break:
            self.break_solutions = []
            return

        # Initialize break solutions
        self.break_solutions = []

        if self.circular:
            if len(self.breaks) > 0:
                # Determine number of solutions per oligo
                self.num_solutions_per_oligo = math.ceil(1.0*num_solutions/len(self.breaks))
        else:
            self.num_solutions_per_oligo = num_solutions

        if self.circular:
            for current_break in tqdm(self.breaks, desc='Shortest path loop ',
                                      dynamic_ncols=True, bar_format='{l_bar}{bar}',
                                      file=self.origami.tqdm_output_file):

                # Check the constraints. If not allowed to break, skip
                if current_break.dont_break or current_break.dont_break_temp:
                    continue

                # Reset break path variables
                self.reset_break_paths()
                # Reset break id numbers
                self.reset_break_order_ids(current_break, current_break)

                # Generate the shortest k-paths
                shortest_k_paths = current_break.get_k_shortest_paths(current_break, self.num_solutions_per_oligo,
                                                                      k_select)
                # Add solutions to solutions list
                self.break_solutions += shortest_k_paths

        else:
            # Reset break path variables
            self.reset_break_paths()
            # Reset break id numbers
            self.reset_break_order_ids(self.start_break, self.final_break)

            # Generate shortest paths
            shortest_k_paths = self.start_break.get_k_shortest_paths(self.final_break, self.num_solutions_per_oligo,
                                                                     k_select)

            # Add solutions to solutions list
            self.break_solutions += shortest_k_paths

        # Calculate self penalty scores for the solutions
        for break_solution in self.break_solutions:
            break_solution.calculate_self_penalty()

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

        self.json_input     = None
        self.oligos         = {'scaffold': [], 'staple': []}
        self.oligo_map      = {}
        self.break_edge_map = {}
        self.oligo_groups   = None

        self.cadnano_oligos = None
        self.staples        = None
        self.scaffolds      = None
        self.idnums         = None

        # Structure parameter
        self.num_crossovers = 0

        # Keeps whether very long staples exist
        self.very_long_staples_exist      = True
        self.dont_break_very_long_staples = False

        # DNA Sequence parameters
        self.sequence_offset    = None
        self.corrected_offset   = 0
        self.sequence_start_pos = None
        self.current_start_pos  = None

        # tqdm output file
        self.tqdm_output_file   = None
        self.std_output_file    = None

    def set_std_output_file(self, filename=None):
        '''
        Set tqdm output file
        '''
        self.std_output_file = filename

    def set_tqdm_output_file(self, filename=None):
        '''
        Set tqdm output file
        '''
        self.tqdm_output_file = filename

    def warn_circular_scaffold(self):
        '''
        Warn user about adding a break in a scaffold
        '''
        cadnano_oligos = self.part.oligos()
        cadnano_oligos = sorted(cadnano_oligos, key=lambda x: x.length(), reverse=True)
        if cadnano_oligos[0].isCircular():
            print('Warning: Scaffold is circular.' +
                  ' It is recommended to add a break in scaffold before applying autobreak.')

    def circularize_scaffold(self):
        '''
        Circularize scaffold
        '''

        strand5p = self.scaffolds[0].strand5p()
        strand3p = self.scaffolds[0].strand3p()

        # Set sequence start position before circularizing
        if self.sequence_start_pos is None:
            self.sequence_start_pos = (strand5p.idNum(), strand5p.idx5Prime())

        # Connect strand5p and strand3p
        if self.circularize and not self.scaffolds[0].isCircular():
            if(strand5p.idNum() == strand3p.idNum() and
               abs(strand5p.idx5Prime() - strand3p.idx3Prime()) == 1):
                strand3p.merge(strand3p.idx3Prime())
            else:
                strand3p.setConnection3p(strand5p)

            # Read cadnano oligos again
            self.get_oligos()

        # Get the new start position for the scaffold
        strand5p = self.scaffolds[0].strand5p()
        strand3p = self.scaffolds[0].strand3p()

        # Get current start position
        self.current_start_pos = (strand5p.idNum(), strand5p.idx5Prime())

        # Set circularize oligo False so that this function is excuted only once
        self.circularize = False

    def determine_num_crossovers(self):
        '''
        Determine total number of crossovers
        '''
        self.num_crossovers = 0
        for oligo in self.oligos['staple']:
            self.num_crossovers += oligo.num_crossovers

    def assign_scaffold_positions(self):
        '''
        Assign scaffold positions for the strands
        '''
        for oligo in self.oligos['staple']:

            # Assign current strand
            current_strand = oligo.null_strand.next_strand

            # Iterate over each strand
            while current_strand:
                direction = current_strand.direction
                vh        = current_strand.vh
                idx5p     = current_strand.idx5p
                idx3p     = current_strand.idx3p

                # Keep the positions
                current_strand.scaffoldPos  = []

                for idx in range(idx5p, idx3p+direction, direction):
                    # Get scaffold position
                    scaffold_positions = self.get_scaffold_positions(vh, idx)

                    # Add new positions to list
                    current_strand.scaffoldPos.extend(scaffold_positions)

                # Update current strand
                current_strand = current_strand.next_strand

    def set_dont_break_oligos(self, maximum_length=0):
        '''
        Set dont break oligos
        '''
        for oligo in self.oligos['staple']:
            if oligo.length < maximum_length:
                oligo.dont_break = True

    def set_sequence_offset(self):
        '''
        Set sequence offset based on the difference between
        scaffold_start_pos and current_start_pos

        '''
        # Determine the sequence distance
        # between scaffold_start_pos (original) and current_start_pos(after circularize)

        if self.sequence_offset is None:
            sequence_distance = self.get_scaffold_distance(self.sequence_start_pos,
                                                           self.current_start_pos)

            # Assign the sequence offset
            self.sequence_offset = sequence_distance % self.scaffolds[0].length()

    def set_circularize(self, circularize=False):
        '''
        Set circularize (scaffold) parameter
        '''
        self.circularize = circularize

    def set_sequence_start_pos(self, key):
        '''
        Set sequence start position in (vh, idx)
        '''
        if key != (-1, -1):
            self.sequence_start_pos = key

    def set_sequence_start_offset(self, offset=0):
        '''
        Set sequence offset
        '''
        if offset >= 0:
            self.sequence_start_offset = offset
        else:
            self.sequence_start_offset = 0

    def get_scaffold_distance(self, current_key, next_key):
        '''
        Get scaffold distance
        Distance is between the minimum positions in the positions list
        '''
        current_vh, current_idx = current_key
        next_vh, next_idx = next_key
        current_positions = self.get_scaffold_positions(current_vh, current_idx)
        next_positions = self.get_scaffold_positions(next_vh, next_idx)

        if len(current_positions) == 0 or len(next_positions) == 0:
            return None
        else:
            return min(next_positions) - min(current_positions)

    def get_scaffold_positions(self, vh, idx):
        '''
        Get scaffold position
        '''
        # Make the key
        key = (vh, idx)

        if key in self.key2scaffold:
            if self.key2scaffold[key] == -1:
                return []
            else:
                return self.key2scaffold[key]
        else:
            return [None]

    def get_current_nucleotide(self, key):
        '''
        Get nucleotide from nucleotide map
        '''
        if key in self.nucleotide_map:
            return self.nucleotide_map[key]
        else:
            return None

    def get_next_nucleotide(self, key):
        '''
        Get next nucleotide
        '''
        if key in self.nucleotide_map:
            return self.nucleotide_map[key].next_nucleotide
        else:
            return None

    def is_dsDNA(self, vh, idx):
        '''
        Check if it is dsDNA
        '''
        forward_strand = self.get_cadnano_strand(vh, idx, True)
        reverse_strand = self.get_cadnano_strand(vh, idx, False)

        if forward_strand and reverse_strand:
            return True
        else:
            return False

    def color_by_Tm(self):
        '''
        Color oligos by Tm
        '''

    def set_dont_break_very_long_staples(self, value=False):
        '''
        Set break very long staples
        '''
        self.dont_break_very_long_staples = value

    def sort_staples_by_length(self):
        '''
        Sort oligos by length
        '''
        self.oligos['staple'].sort(key=lambda x: x.length)

    def sort_staples_by_key(self):
        '''
        Sort staples by key
        '''
        self.oligos['staple'].sort(key=lambda x: x.key)

    def sort_breaks_by_key(self):
        '''
        Sort breaks by key
        '''
        self.breaks.sort(key=lambda x: x.key)

    def set_sequence_file(self, sequence_file=None):
        '''
        Set sequence filename
        '''
        self.sequence_file = sequence_file

    def break_very_long_staples(self):
        '''
        Break very long oligos
        '''
        if self.dont_break_very_long_staples:
            return

        for oligo in self.oligos['staple']:
            if len(oligo.breaks) > self.autobreak.max_num_breaks:
                oligo.break_in_half()

    def is_there_very_long_staples(self):
        '''
        Returns whether there are any very long staples
        '''
        return self.very_long_staples_exist

    def check_very_long_staples(self):
        '''
        Check if there is any long staples
        '''
        self.very_long_staples_exist = False
        for oligo in self.oligos['staple']:
            if len(oligo.breaks) > self.autobreak.max_num_breaks:
                self.very_long_staples_exist = True
                break

    def prepare_origami(self):
        '''
        List of commands to prepare origami for break
        '''
        # Get oligos
        self.get_oligos()

        # Circularize scaffold
        self.circularize_scaffold()

        # Read scaffolds
        self.read_scaffolds()

        # Build scaffold map
        self.build_scaffold_map()

        # Determine the proper sequence offset
        self.set_sequence_offset()

        # Reset oligos list
        self.reset_oligos()

        # Read sequence file
        self.read_sequence()

        # Read scaffolds
        self.read_scaffolds()

        # Read scaffolds and staples
        self.read_staples()

        # Build scaffold map
        self.build_scaffold_map()

        # Apply sequence
        self.apply_sequence(self.sequence_offset)

        # Assign strand sequences
        self.assign_strands_dna()

        # Assign scaffold positions
        self.assign_scaffold_positions()

        # Sort staple by length
        self.sort_staples_by_length()

        # Sort staples by key
        self.sort_staples_by_key()

        # Generate staple crossovers
        self.generate_staple_crossovers()

        # Generate scaffold crossovers
        self.generate_scaffold_crossovers()

        # Link crossovers
        self.link_crossovers()

        # Build nucleotide map
        self.build_nucleotide_map()

        # Generate dsDNA sequences
        self.generate_dsDNA_sequences()

        # Generate ssDNA sequences
        self.generate_ssDNA_sequences()

        # Connect the sequences
        self.connect_sequences()

        # Determine total crossovers
        self.determine_num_crossovers()

        # Apply break rules
        self.apply_break_rules()

        # Generate break points
        self.generate_break_points()

        # Connect the break boints and set break constraints based on connectivity
        self.connect_break_points()

        # Apply crossover rule
        self.apply_cross_rule()

    def set_cadnano_sequence_offset(self):
        '''
        Set cadnano sequenceOffset
        '''
        self.part.setSequenceOffset(self.corrected_offset)

    def get_cadnano_strand(self, vh, idx, direction):
        '''
        Get cadnano strand from vh, idx, direction information
        '''
        return self.part.getStrand(direction > 0, vh, idx)

    def split_cadnano_strand(self, vh, idx, direction):
        '''
        Split cadnano strand at vh, idx, direction
        '''
        new_strand = self.part.getStrand(direction > 0, vh, idx)

        # Break only if it is a valid strand
        if new_strand:
            new_strand.split(idx)

    def remove_cadnano_crossover(self, vh, idx, direction):
        '''
        Remove cadnano crossover
        '''
        strand5p = self.part.getStrand(direction > 0, vh, idx)
        strand3p = strand5p.connection3p()

        self.part.removeXover(strand5p, strand3p)

    def initialize(self, input_filename):
        # Initialize cadnano
        app = cadnano.app()
        self.doc = app.document = Document()

        # Assign cadnano input file
        self.json_input = input_filename

        # Read cadnano input file
        self.doc.readFile(self.json_input)

        # Assign part
        self.part = self.doc.activePart()

    def split_scaffold(self):
        '''
        Split scaffold at the starting position - (direction)
        '''

        (vh, idx) = self.sequence_start_pos

        # Determine direction
        if vh % 2 == 0:
            direction =  1
        else:
            direction = -1

        # Break the scaffold at start position
        self.split_cadnano_strand(vh, idx-direction, direction)

    def assign_strands_dna(self):
        '''
        Assign strand sequences from cadnano part
        '''
        for oligo in self.oligos['staple']:

            # Assign current strand
            current_strand = oligo.null_strand.next_strand

            while current_strand:
                current_strand.dna = current_strand.cadnano_strand.sequence()

                # Current cadnano2.5 version doesnt assign sequences for strands
                # with no scaffold in the same virtual helix
                # Assign the proper empty length sequence in case it is not assigned

                if len(current_strand.dna) != current_strand.totalLength:
                    current_strand.dna = ' '*current_strand.totalLength

                # Update current strand
                current_strand = current_strand.next_strand

    def read_oligo(self, oligo, oligo_type='staple'):
        '''
        Read oligo from 5' to 3'
        '''
        # Oligo strand generator
        generator    = oligo.strand5p().generator3pStrand()

        # Create a null strand object for oligo
        previous_strand             = Strand()
        previous_strand.length      = 0
        previous_strand.totalLength = 0
        previous_strand.distance    = 0

        # Create Oligo object
        new_oligo               = Oligo()
        new_oligo.type          = oligo_type
        new_oligo.circular      = oligo.isCircular()
        new_oligo.null_strand   = previous_strand
        new_oligo.length        = oligo.length()
        new_oligo.origami       = self
        new_oligo.cadnano_oligo = oligo

        # Get 5p strand
        strand5p              = oligo.strand5p()
        idx5p                 = strand5p.idx5Prime()
        idx3p                 = strand5p.idx3Prime()
        vh                    = strand5p.idNum()
        direction             = -1 + 2*strand5p.isForward()

        # Assign key for oligo
        new_oligo.key         = (vh, idx5p, direction)

        # Add oligo to list
        self.oligos[oligo_type].append(new_oligo)

        # Add oligo to map
        self.oligo_map[new_oligo.key] = new_oligo

        # Get Strand parameters
        for strand in generator:
            # Create new strand
            new_strand                    = Strand()

            new_strand.cadnano_strand     = strand
            new_strand.vh                 = strand.idNum()
            new_strand.idx5p              = strand.idx5Prime()
            new_strand.idx3p              = strand.idx3Prime()
            new_strand.forward            = strand.isForward()
            new_strand.idxLow, new_strand.idxHigh = ((new_strand.idx5p, new_strand.idx3p)
                                                     if new_strand.forward else (new_strand.idx3p, new_strand.idx5p))
            new_strand.direction          = -1 + 2*strand.isForward()
            new_strand.complement_strands = strand.getComplementStrands()[::new_strand.direction]
            new_strand.length             = new_strand.direction*(new_strand.idx3p-new_strand.idx5p)+1
            new_strand.totalLength        = strand.totalLength()
            new_strand.distance           = previous_strand.distance + previous_strand.totalLength
            new_strand.origami            = self

            # Prepare the insert/skip list
            new_strand.inserts            = []
            for idx in range(idx5p, idx3p + direction, direction):
                # Initialize insert size
                insert_size = 0
                if strand.hasInsertionAt(idx):
                    insert_size = strand.insertionLengthBetweenIdxs(idx, idx)
                new_strand.inserts.append(insert_size)

            # Make a numpy array
            new_strand.inserts = np.array(new_strand.inserts)

            # Make the strand connection
            previous_strand.next_strand = new_strand

            # Update previous trand
            previous_strand = new_strand

        # If oligo is not circular make the last strand final strand
        if not new_oligo.circular:
            previous_strand.final_strand = True

    def generate_staple_crossovers(self):
        '''
        Generate staple crossover objects
        '''

        # Initialize the crossovers for the origami
        self.crossovers = {}

        # Read staple crossovers
        for oligo in self.oligos['staple']:
            # Initialize crossovers for the oligo
            oligo.crossovers = {}

            # Get current strand
            current_strand = oligo.null_strand.next_strand

            # Iterate over the strands
            while current_strand:
                if not current_strand.final_strand:

                    # Get the connected strand
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
                    new_crossover.type           = 'staple'
                    new_crossover.neighbor       = None
                    new_crossover.neighbor_key   = (new_crossover.vh3p, new_crossover.idx+new_crossover.direction,
                                                    -new_crossover.direction)

                    # Add crossovers to the list
                    oligo.crossovers[new_crossover.key] = new_crossover
                    self.crossovers[new_crossover.key]  = new_crossover

                # Update current strand
                current_strand          = current_strand.next_strand

    def generate_scaffold_crossovers(self):
        '''
        Generate scaffold crossovers
        '''
        # Read scaffold crossovers
        for oligo in self.oligos['scaffold']:
            # Initialize crossovers for the oligo
            oligo.crossovers = {}

            # Get current strand
            current_strand = oligo.null_strand.next_strand

            # Iterate over the strands
            while current_strand:
                if not current_strand.final_strand:

                    # Get the connected strand
                    next_strand   = current_strand.next_strand or oligo.null_strand.next_strand

                    new_crossover = Crossover()
                    new_crossover.vh5p           = current_strand.vh
                    new_crossover.vh3p           = next_strand.vh
                    new_crossover.idx            = current_strand.idx3p+current_strand.direction
                    new_crossover.direction      = -current_strand.direction
                    new_crossover.key            = (new_crossover.vh5p, new_crossover.idx, new_crossover.direction)
                    new_crossover.oligo          = oligo
                    new_crossover.current_strand = current_strand
                    new_crossover.next_strand    = next_strand
                    new_crossover.type           = 'scaffold'
                    new_crossover.neighbor       = None
                    new_crossover.neighbor_key   = (new_crossover.vh3p, new_crossover.idx+new_crossover.direction,
                                                    -new_crossover.direction)

                    # Add crossovers to the list
                    oligo.crossovers[new_crossover.key] = new_crossover
                    self.crossovers[new_crossover.key]  = new_crossover

                # Update current strand
                current_strand          = current_strand.next_strand

    def link_crossovers(self):
        '''
        Link crossovers
        '''
        # Check crossover neighbors
        for key in self.crossovers:

            # Get neighbor key
            neighbor_key = self.crossovers[key].neighbor_key

            # Make the neighbor connection
            if neighbor_key in self.crossovers:
                self.crossovers[key].neighbor = self.crossovers[neighbor_key]

    def generate_dsDNA_sequences(self):
        '''
        Generate sequences from strands
        '''

        for oligo in self.oligos['staple']:
            # Initialize sequence list for the oligo
            oligo.sequences = []

            current_strand = oligo.null_strand.next_strand

            while current_strand:
                # Make null sequence object
                previous_sequence = Sequence()
                current_strand.null_sequence               = previous_sequence
                current_strand.null_sequence.next_sequence = None

                # Get complementary strands
                complement_strands = current_strand.complement_strands

                # Strand sequence positions
                current_strand.sequence_idxLows = []

                # Sequence array
                current_strand.sequences = []

                # Go through each complement strand to determine the boundaries
                for strand in complement_strands:
                    comp_idx5p   = strand.idx5Prime()
                    comp_idx3p   = strand.idx3Prime()
                    comp_forward = strand.isForward()

                    # Get the low and high indexes for complementary strand
                    comp_idxLow, comp_idxHigh = (comp_idx5p, comp_idx3p) if comp_forward else (comp_idx3p, comp_idx5p)

                    # Make sequence object
                    new_sequence         = Sequence()
                    new_sequence.idNum   = current_strand.vh
                    new_sequence.idxLow  = max(current_strand.idxLow, comp_idxLow)
                    new_sequence.idxHigh = min(current_strand.idxHigh, comp_idxHigh)
                    new_sequence.length  = new_sequence.idxHigh-new_sequence.idxLow+1
                    new_sequence.idx5p, new_sequence.idx3p = ((new_sequence.idxLow, new_sequence.idxHigh)
                                                              if current_strand.forward
                                                              else (new_sequence.idxHigh, new_sequence.idxLow))
                    new_sequence.forward   = current_strand.forward
                    new_sequence.type      = 'dsDNA'
                    new_sequence.direction = current_strand.direction
                    new_sequence.strand    = current_strand
                    new_sequence.origami   = self

                    # Assign string indexes
                    new_sequence.strLow  = (current_strand.direction*(new_sequence.idx5p - current_strand.idx5p) +
                                            current_strand.get_inserts(new_sequence.idx5p, current_strand.idx5p))
                    new_sequence.strHigh = (current_strand.direction*(new_sequence.idx3p - current_strand.idx5p) +
                                            current_strand.get_inserts(new_sequence.idx3p, current_strand.idx5p))

                    # Get scaffold positions
                    new_sequence.assign_scaffold_positions()

                    # Get the total length for the sequence
                    new_sequence.totalLength = new_sequence.strHigh - new_sequence.strLow + 1

                    # Assign sequence distance from 5' end of oligo
                    new_sequence.distance = new_sequence.strLow + current_strand.distance

                    # Get the sequence
                    new_sequence.dna = current_strand.dna[new_sequence.strLow:new_sequence.strHigh+1]

                    # Make the next sequence link
                    previous_sequence.next_sequence = new_sequence

                    # Make the previous sequence link
                    new_sequence.previous_sequence = previous_sequence

                    # Update previous sequence
                    previous_sequence = new_sequence

                # Make the sequence starting position numpy array
                current_strand.sequence_idxLows = np.array(current_strand.sequence_idxLows)

                # Update current strand
                current_strand = current_strand.next_strand

    def generate_ssDNA_sequences(self):
        '''
        Genereate ssDNA sequences for oligo
        '''
        for oligo in self.oligos['staple']:

            # Initialize sequence list for the oligo
            oligo.sequences = []

            # Get the first strand
            current_strand = oligo.null_strand.next_strand

            while current_strand:

                # Initialize the start and final idx for the sequence
                start_idx    = current_strand.idx5p
                final_idx    = -1

                # Get current sequence
                current_sequence = current_strand.null_sequence.next_sequence

                # Get previous sequence
                previous_sequence = current_strand.null_sequence

                # Strand sequence positions
                current_strand.sequence_idxLows = []

                # Sequence array
                current_strand.sequences = []

                while current_sequence:

                    # Get final idx
                    final_idx = current_sequence.idx5p - current_sequence.direction

                    # Criteria to make a new ssDNA sequence
                    if current_sequence.direction*(final_idx - start_idx) > 0:
                        # Make sequence object
                        new_sequence         = Sequence()
                        new_sequence.idNum   = current_strand.vh
                        new_sequence.idxLow  = min(start_idx, final_idx)
                        new_sequence.idxHigh = max(start_idx, final_idx)
                        new_sequence.length  = new_sequence.idxHigh-new_sequence.idxLow+1
                        new_sequence.idx5p, new_sequence.idx3p = (start_idx, final_idx)
                        new_sequence.forward   = current_strand.forward
                        new_sequence.type      = 'ssDNA'
                        new_sequence.direction = current_strand.direction
                        new_sequence.strand    = current_strand
                        new_sequence.origami   = self

                        # Assign string indexes
                        new_sequence.strLow  = (current_strand.direction*(new_sequence.idx5p - current_strand.idx5p) +
                                                current_strand.get_inserts(new_sequence.idx5p, current_strand.idx5p))
                        new_sequence.strHigh = (current_strand.direction*(new_sequence.idx3p - current_strand.idx5p) +
                                                current_strand.get_inserts(new_sequence.idx3p, current_strand.idx5p))

                        # Assign scaffold positions
                        new_sequence.scaffoldPos = []

                        # Get the total length for the sequence
                        new_sequence.totalLength = new_sequence.strHigh - new_sequence.strLow + 1

                        # Assign sequence distance from 5' end of oligo
                        new_sequence.distance = new_sequence.strLow + current_strand.distance

                        # Get the sequence
                        new_sequence.dna = current_strand.dna[new_sequence.strLow:new_sequence.strHigh+1]

                        # Make the link between the sequences
                        previous_sequence.next_sequence = new_sequence
                        new_sequence.previous_sequence  = previous_sequence

                        new_sequence.next_sequence = current_sequence
                        current_sequence.previous_sequence = new_sequence

                        # Keep the low position
                        current_strand.sequence_idxLows.append(new_sequence.strLow)

                        # Add new sequence to strand sequence list
                        current_strand.sequences.append(new_sequence)

                        # Add new sequence to oligo sequence list
                        oligo.sequences.append(new_sequence)

                    # Add the current sequence values to the lists
                    current_strand.sequence_idxLows.append(current_sequence.strLow)
                    current_strand.sequences.append(current_sequence)
                    oligo.sequences.append(current_sequence)

                    # Update start idx
                    start_idx = current_sequence.idx3p + current_sequence.direction

                    # Update sequences
                    previous_sequence = current_sequence
                    current_sequence  = current_sequence.next_sequence

                # Make the 3p terminal ssDNA sequence
                final_idx   = current_strand.idx3p
                if current_strand.direction*(final_idx - start_idx) > 0:
                    # Make sequence object
                    new_sequence         = Sequence()
                    new_sequence.idNum   = current_strand.vh
                    new_sequence.idxLow  = min(start_idx, final_idx)
                    new_sequence.idxHigh = max(start_idx, final_idx)
                    new_sequence.length  = new_sequence.idxHigh-new_sequence.idxLow+1
                    new_sequence.idx5p, new_sequence.idx3p = (start_idx, final_idx)
                    new_sequence.forward   = current_strand.forward
                    new_sequence.type      = 'ssDNA'
                    new_sequence.direction = current_strand.direction
                    new_sequence.strand    = current_strand

                    # Assign string indexes
                    new_sequence.strLow  = (current_strand.direction*(new_sequence.idx5p - current_strand.idx5p) +
                                            current_strand.get_inserts(new_sequence.idx5p, current_strand.idx5p))
                    new_sequence.strHigh = (current_strand.direction*(new_sequence.idx3p - current_strand.idx5p) +
                                            current_strand.get_inserts(new_sequence.idx3p, current_strand.idx5p))

                    # Assign scaffold positions
                    new_sequence.scaffoldPos = []

                    # Get the total length for the sequence
                    new_sequence.totalLength = new_sequence.strHigh - new_sequence.strLow + 1

                    # Assign sequence distance from 5' end of oligo
                    new_sequence.distance = new_sequence.strLow + current_strand.distance

                    # Get the sequence
                    new_sequence.dna = current_strand.dna[new_sequence.strLow:new_sequence.strHigh+1]

                    # Make the link between the sequences
                    previous_sequence.next_sequence = new_sequence
                    new_sequence.previous_sequence  = previous_sequence

                    new_sequence.next_sequence = None

                    # Keep the low position
                    current_strand.sequence_idxLows.append(new_sequence.strLow)

                    # Add new sequence to strand sequence list
                    current_strand.sequences.append(new_sequence)

                    # Add new sequence to oligo sequence list
                    oligo.sequences.append(new_sequence)

                # Make the sequence starting position numpy array
                current_strand.sequence_idxLows = np.array(current_strand.sequence_idxLows)

                # Update current strand
                current_strand = current_strand.next_strand

    def connect_sequences(self):
        '''
        Make the connection between strand sequences
        '''
        for oligo in self.oligos['staple']:

            oligo.num_crossovers = 0

            # Get the first strand
            current_strand = oligo.null_strand.next_strand

            # Assign previous  and current sequences
            previous_sequence = current_strand.null_sequence

            while current_strand:
                # Get current sequence
                current_sequence = current_strand.null_sequence.next_sequence

                # Iterate over the sequence
                while current_sequence:
                    # Update sequence numbers
                    oligo.num_crossovers += int(current_sequence.type == 'dsDNA')

                    # Make the links
                    previous_sequence.next_sequence = current_sequence
                    current_sequence.previous_sequence = previous_sequence

                    # Update previous and current sequence
                    previous_sequence = current_sequence
                    current_sequence = current_sequence.next_sequence

                # Update current strand
                current_strand = current_strand.next_strand

            # For circular oligos connect last sequence to first sequence
            if oligo.circular:
                previous_sequence.next_sequence = oligo.null_strand.next_strand.sequences[0]

                # For circular oligos add one more crossover
                oligo.num_crossovers += 1

    def update_sequences_dna(self):
        '''
        Update the sequences dna
        '''
        for oligo in self.oligos['staple']:

            # Iterate over the sequences
            for current_sequence in oligo.sequences:
                current_sequence.dna = current_sequence.strand.dna[current_sequence.strLow:current_sequence.strHigh+1]

    def apply_break_rules(self):
        '''
        Apply break rules
        '''
        for oligo in self.oligos['staple']:
            # Assign current strand
            current_strand = oligo.null_strand.next_strand

            while current_strand:
                # Apply the break rules
                current_strand.apply_break_rule()

                # Update current strand
                current_strand   = current_strand.next_strand

    def apply_cross_rule(self):
        '''
        Apply cross rule
        '''
        # If xscaf is not in the rule list disable scaffold crossover break points
        if 'xscaf' not in self.break_rule:
            self.disable_scaffold_crossovers()

        # If xstap is not in the rule list disable staple crossover break points
        if 'xstap' not in self.break_rule:
            self.disable_staple_crossovers()

    def determine_longrange_breaks(self):
        '''
        Determine long range contacts
        '''
        self.long_range_breaks = []

    def build_nucleotide_map(self):
        '''
        Make nucleotide map
        '''
        # Initialize nucleotide map
        self.nucleotide_map = {}

        for oligo in self.oligos['staple']:

            # Make a dummy nucleotide
            previous_nucleotide = Nucleotide()
            oligo.null_nucleotide = previous_nucleotide

            # Get current strand
            current_strand = oligo.null_strand.next_strand

            while current_strand:
                idx5p = current_strand.idx5p
                idx3p = current_strand.idx3p
                direction = current_strand.direction
                for idx in range(idx5p, idx3p+direction, direction):
                    # Make a new nucleotide
                    new_nucleotide = Nucleotide()
                    new_nucleotide.vh = current_strand.vh
                    new_nucleotide.idx = idx
                    new_nucleotide.direction = direction
                    new_nucleotide.key = (new_nucleotide.vh, new_nucleotide.idx,
                                          new_nucleotide.direction)
                    new_nucleotide.dsDNA = self.is_dsDNA(new_nucleotide.vh,
                                                         new_nucleotide.idx)
                    # Assign new nucleotide
                    self.nucleotide_map[new_nucleotide.key] = new_nucleotide

                    # Make the links
                    previous_nucleotide.next_nucleotide = new_nucleotide
                    new_nucleotide.previous_nucleotide = previous_nucleotide

                    # Update previous nucleotide
                    previous_nucleotide = new_nucleotide

                # Update current strand
                current_strand = current_strand.next_strand

    def generate_break_points(self):
        '''
        Generate break points
        '''

        # Initialize breaks
        self.breaks = []

        for oligo in self.oligos['staple']:
            # Initialize oligo breaks and break map
            oligo.breaks = []

            # Assign current strand
            current_strand = oligo.null_strand.next_strand

            # Set break id counter for graph algorithms
            order_id_counter = 0

            # Create a null break point representing the 5'-end of the nucleotide
            previous_break   = autobreak.BreakNode()
            oligo.null_break = previous_break

            # Initiliaze null break point
            oligo.null_break.break_point          = -1
            oligo.null_break.break_point_adjusted = -1
            oligo.null_break.idx                  = (current_strand.idx5p +
                                                     current_strand.direction*oligo.null_break.break_point)
            oligo.null_break.vh                   = current_strand.vh
            oligo.null_break.direction            = current_strand.direction
            oligo.null_break.distance             = 0
            oligo.null_break.strand               = current_strand
            oligo.null_break.key                  = (oligo.null_break.vh, oligo.null_break.idx,
                                                     oligo.null_break.direction)
            oligo.null_break.current_nucleotide   = self.get_current_nucleotide(oligo.null_break.key)
            oligo.null_break.next_nucleotide      = self.get_next_nucleotide(oligo.null_break.key)
            oligo.null_break.sequence             = (current_strand.sequences[0] if len(current_strand.sequences) > 0
                                                     else None)
            oligo.null_break.oligo                = oligo
            oligo.null_break.order_id             = order_id_counter
            oligo.null_break.origami              = self
            oligo.null_break.dsDNA                = oligo.null_break.is_dsDNA()
            oligo.null_break.insert               = current_strand.get_inserts(oligo.null_break.idx,
                                                                               oligo.null_break.idx)

            # Decide if the break cant be broken
            if not oligo.null_break.dsDNA or oligo.null_break.insert == -1:
                oligo.null_break.dont_break = True

            # Add null oligo to break list
            oligo.breaks.append(oligo.null_break)

            # If the oligo is not circular make the first break node start node
            if not oligo.circular:
                oligo.start_break = oligo.null_break

            # Iterate over each strand
            while current_strand:
                # Iterate through all positions
                for i in range(len(current_strand.all_breaks)):
                    # Get break position
                    break_position = current_strand.all_breaks[i]

                    # Get adjusted break position
                    break_position_adjusted = current_strand.all_breaks_adjusted[i]

                    # Update break id counter
                    order_id_counter += 1

                    # Find the sequence for the break point
                    sequence_id = np.searchsorted(current_strand.sequence_idxLows, break_position, side='right') - 1

                    # Make a break
                    new_break   = autobreak.BreakNode()

                    # Break position
                    new_break.break_point          = break_position
                    new_break.break_point_adjusted = break_position_adjusted
                    new_break.idx                  = current_strand.idx5p + current_strand.direction*break_position
                    new_break.vh                   = current_strand.vh
                    new_break.direction            = current_strand.direction
                    new_break.distance             = current_strand.distance + break_position_adjusted + 1
                    new_break.key                  = (new_break.vh, new_break.idx, new_break.direction)
                    new_break.current_nucleotide   = self.get_current_nucleotide(new_break.key)
                    new_break.next_nucleotide      = self.get_next_nucleotide(new_break.key)
                    new_break.order_id             = order_id_counter
                    new_break.origami              = self
                    new_break.dsDNA                = new_break.is_dsDNA()
                    new_break.insert               = current_strand.get_inserts(new_break.idx, new_break.idx)
                    new_break.crossover            = None
                    new_break.location             = 'internal'

                    # If the break point is not located between dsDNA segments, dont break
                    if not new_break.dsDNA or new_break.insert != 0:
                        new_break.dont_break = True

                    # Assign sequence to break object
                    new_break.sequence    = current_strand.sequences[sequence_id] if sequence_id >= 0 else None

                    if new_break.sequence is None:
                        print(new_break.key, current_strand.length, break_position, current_strand.sequence_idxLows)

                    # Assign strand to new break
                    new_break.strand      = current_strand

                    # Assign oligo
                    new_break.oligo       = oligo

                    # Check if the break is at a cross-over location
                    if new_break.key in self.crossovers:
                        new_break.crossover = self.crossovers[new_break.key]
                        self.crossovers[new_break.key].break_node = new_break

                    # Assign to previous break
                    previous_break.next_break = new_break

                    # Make the previous connection
                    new_break.previous_break  = previous_break

                    # Update previous break
                    previous_break = new_break

                    # Add break to break list
                    oligo.breaks.append(new_break)

                # Update current strand
                current_strand   = current_strand.next_strand

            # If oligo is circular, connect final node to null break's next break
            if oligo.circular:
                # Update the final node to first node forward connection
                previous_break.next_break = oligo.null_break.next_break

                # Update the final node to first node reverse connection
                if oligo.null_break.next_break:
                    oligo.null_break.next_break.previous_break = previous_break

                # Remove the first break node from list
                oligo.breaks.pop(0)

                # Assign final break to null break
                oligo.null_break  = previous_break
            else:
                # If the oligo is not circular make the final break node final node
                oligo.final_break = previous_break

                # Make the start and final break valid for breaking
                oligo.start_break.dont_break = False
                oligo.final_break.dont_break = False

                # Set the terminus break types
                oligo.start_break.location = 'terminus'
                oligo.final_break.location = 'terminus'

            # Add breaks to origami list
            self.breaks += oligo.breaks

    def disable_staple_crossovers(self):
        '''
        Disable all crossover break nodes
        '''
        for key in self.crossovers:
            if self.crossovers[key].break_node is None or self.crossovers[key].type == 'scaffold':
                continue
            if self.crossovers[key].break_node.location == 'internal':
                self.crossovers[key].break_node.dont_break = True
            else:
                self.crossovers[key].break_node.dont_break = False

    def disable_scaffold_crossovers(self):
        '''
        Disable all crossover break nodes
        '''
        for key in self.crossovers:
            if self.crossovers[key].break_node is None or self.crossovers[key].type == 'staple':
                continue
            if self.crossovers[key].break_node.location == 'internal':
                self.crossovers[key].break_node.dont_break = True
            else:
                self.crossovers[key].break_node.dont_break = False

    def connect_break_points(self):
        '''
        Connect break points
        '''
        for oligo in self.oligos['staple']:

            # Visit each break object
            for current_break in oligo.breaks:

                # Get the crossover for the current break
                current_crossover = current_break.crossover

                # Get neighbor break
                if current_crossover:
                    # Determine break type for the strategy to break
                    if current_crossover.type == 'staple':
                        current_break.type = 'crossover'
                    else:
                        current_break.type = 'break'

                    # If a neighbor exists make the connection
                    if current_crossover.neighbor:
                        current_break.neighbor_break = current_crossover.neighbor.break_node
                    # If there is a neighbor and if it is an internal scaffold crossover
                    elif current_break.location == 'internal':
                        current_break.dont_break = True
                else:
                    current_break.type = 'break'

    def cluster_oligo_groups(self):
        '''
        Cluster oligos based on break connectivity
        '''
        # Initialize oligo groups
        self.oligo_groups = []

        # Group key
        group_key = 0

        # Sort breaks by key
        self.sort_breaks_by_key()

        for current_break in self.breaks:
            if not current_break.visited:
                # Create new oligo group
                new_oligo_group = OligoGroup()

                # Perform depth first search starting from current break
                visited_breaks  = current_break.depth_first_search()

                # Get the oligos for the visited breaks
                visited_oligos  = set([new_break.oligo for new_break in visited_breaks])

                # For all visited breaks make the oligo group assignment
                for new_break in visited_breaks:
                    new_break.oligo_group = new_oligo_group

                # For all visited oligos make the oligo group assignment
                for new_oligo in visited_oligos:
                    new_oligo.oligo_group = new_oligo_group

                # Assign breaks and oligos to oligo group
                new_oligo_group.breaks = list(visited_breaks)
                new_oligo_group.oligos = list(visited_oligos)

                # Sort group breaks and oligos by key
                new_oligo_group.breaks.sort(key=lambda x: x.key)
                new_oligo_group.oligos.sort(key=lambda x: x.key)

                # Assign group key
                new_oligo_group.key    = group_key

                # Assign origami to oligo group
                new_oligo_group.origami = self

                # Update oligo group key
                group_key += 1

                # Add oligo group to the list
                self.oligo_groups.append(new_oligo_group)

    def read_staple(self, staple):
        '''
        Read staple from 5' to 3'
        '''
        self.read_oligo(staple, oligo_type='staple')

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
        # Initialize idnums list
        self.idnums = []
        for scaffold in self.scaffolds:
            self.read_scaffold(scaffold)

    def reset_oligos(self):
        '''
        Reset oligos list
        '''
        self.oligos    = {'scaffold': [], 'staple': []}
        self.oligo_map = {}

    def read_scaffold(self, scaffold):
        '''
        Read scaffold
        '''
        self.read_oligo(scaffold, oligo_type='scaffold')

    def build_scaffold_map(self):
        '''
        Build scaffold map
        '''
        # Initialize scaffold maps
        self.scaffold2key = {}
        self.key2scaffold = {}

        # Get the first scaffold
        scaffold = self.oligos['scaffold'][0]

        # Position counter
        position = 0

        # Assign current strand
        current_strand = scaffold.null_strand.next_strand

        while current_strand:
            # Get the position parameters
            idx5p = current_strand.idx5p
            idx3p = current_strand.idx3p
            direction = current_strand.direction
            vh = current_strand.vh

            for idx in range(idx5p, idx3p+direction, direction):
                key = (vh, idx)
                num_inserts = current_strand.get_inserts(idx, idx)
                position += (num_inserts+1)

                # Assign the values
                if position not in self.scaffold2key:
                    self.scaffold2key[position] = key

                if num_inserts == -1:
                    self.key2scaffold[key] = -1
                else:
                    self.key2scaffold[key] = list(range(position-num_inserts, position+1))[::-1]

            # Update current strand
            current_strand = current_strand.next_strand

    def get_oligos(self):
        '''
        Get oligos
        '''
        self.cadnano_oligos = self.part.oligos()

        # In place sort based on vh and idx5p and oligo direction
        self.cadnano_oligos = sorted(self.cadnano_oligos,
                                     key=lambda x: (x.length(),
                                                    x.strand5p().idNum(),
                                                    x.strand5p().idx5Prime(),
                                                    int(x.strand5p().isForward())), reverse=True)

        # Initialize the scaffolds and staples
        self.scaffolds = []
        self.staples   = []

        # Iterate over oligos
        for oligo in self.cadnano_oligos:
            strand5p  = oligo.strand5p()
            vh        = strand5p.idNum()
            isForward = strand5p.isForward()

            # Scaffold criteria: (forward and even-number helix) or (reverse and odd-number helix)
            if int(isForward) + vh % 2 == 1:
                self.scaffolds.append(oligo)
            else:
                self.staples.append(oligo)

    def read_sequence(self):
        '''
        Read sequence file
        '''

        # 1. Check the popular scaffolds list
        # 2. Check if the file exists
        # 3. Generate random sequence

        if self.sequence_file in utilities.SCAFFOLD_SEQUENCES:
            print("Sequence code %s exists." % (self.sequence_file))
            self.scaffold_sequence = utilities.SCAFFOLD_SEQUENCES[self.sequence_file]

        elif self.sequence_file and os.path.isfile(self.sequence_file):
            # Read sequence from file
            f = open(self.sequence_file)
            self.scaffold_sequence = ''.join([line.strip() for line in f.readlines()])
            f.close()

            # Convert to upper case
            self.scaffold_sequence = self.scaffold_sequence.upper()

            # If sequence length is less than scaffold length quit
            if len(self.scaffold_sequence) < self.scaffolds[0].length():
                sys.exit('SCAFFOLD SEQUENCE IS SHORTER THAN SCAFFOLD LENGTH!')
        elif self.sequence_file == 'allT':
            # Make sequence allT
            self.scaffold_sequence = utilities.generate_nT(self.scaffolds[0].length())
        elif self.sequence_file == 'allA':
            # Make sequence allA
            self.scaffold_sequence = utilities.generate_nA(self.scaffolds[0].length())
        elif self.sequence_file == 'allC':
            # Make sequence allT
            self.scaffold_sequence = utilities.generate_nC(self.scaffolds[0].length())
        elif self.sequence_file == 'allG':
            # Make sequence allA
            self.scaffold_sequence = utilities.generate_nG(self.scaffolds[0].length())
        elif len(self.scaffolds) > 0:
            self.scaffold_sequence = utilities.generate_random_sequence(self.scaffolds[0].length())

    def apply_sequence(self, offset=0):
        '''
        Apply sequence to scaffold
        '''
        self.sequence_offset = offset
        if len(self.scaffolds) > 0:
            self.scaffolds[0].applySequence(self.scaffold_sequence[offset:] +
                                            self.scaffold_sequence[:offset])

    def get_coordinates(self, vh, index):
        '''
        Given a vh and a index, returns (x,y,z)
        for the sidechain pts and backbones fwd and rev
        '''

        # Need to reverse the sign of y-axis(it could be any axis) to make the coordinate system right-handed
        # Cadnano coordinate system is left-handed, not compatible with right-handed A-DNA

        axis_pts = self.part.getCoordinates(vh)[0][index]*(1, -1, 1)
        fwd_pts  = self.part.getCoordinates(vh)[1][index]*(1, -1, 1)
        rev_pts  = self.part.getCoordinates(vh)[2][index]*(1, -1, 1)

        return {-1: rev_pts, 0: axis_pts, 1: fwd_pts}
