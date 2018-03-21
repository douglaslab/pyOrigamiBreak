#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2017-08-10 01:32:19
# @Author  : Tural Aksel (turalaksel@gmail.com)
# @Link    : http://example.org
# @Version : $Id$

from tqdm import tqdm
from cadnano.document import Document
from matplotlib import cm
from operator import itemgetter
from shutil import copyfile

import numpy as np
import cadnano
import random
import os
import sys
import argparse
import utilities
import math
import glob
import pandas
import openpyxl

import matplotlib
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


class OligoBreakSolution:
    def __init__(self):
        self.breaks = None
        self.edges  = None

    def get_cvs_rows(self):
        '''
        Prepare cvs writer object

        COLUMNS
        1.  Oligo key
        2.  Oligo Group key
        3.  Break 1 key
        4.  Break 1 type
        5.  Break 2 key
        6.  Break 2 type
        7.  Length
        8.  Sequence
        9.  Edge weight
        10.  ProbFolding
        11.  Log-Probfolding
        12.  Tf
        13. maxTm
        14. maxSeq
        15. has14
        16. dGtotal
        17. dGintrin
        18. dGloop
        19. dGconc
        '''

        # Initialize rows
        cvs_writer_rows = []
        for edge in self.edges[:-1]:
            cvs_writer_rows.append(edge.get_cvs_row_object())

        return cvs_writer_rows

    def break_oligo_solution(self):
        '''
        Break oligo solution
        '''
        for current_break in self.breaks[:-1]:
            current_break.break_cadnano()

    def apply_temp_neighbor_constraints(self):
        '''
        Apply neighbor constraints
        '''
        for new_break in self.breaks[:-1]:

            # Get neighbor break
            if new_break.neighbor_break:
                neighbor_break                 = new_break.neighbor_break
                neighbor_break.dont_break_temp = True

    def print_solution(self):
        '''
        Print break solution
        '''
        if self.breaks:
            tqdm.write('Break path:\n'+'->\n'.join(["(%3d.%3d.%3d)" %
                       (current_break.key) for current_break in self.breaks]))
            tqdm.write('Edge length/weight: '+'->'.join(['(%d / %.1e)' %
                       (edge.edge_length, edge.edge_weight) for edge in self.edges[:-1]]))

    def reset_temp_neighbor_constraints(self):
        '''
        Reset neighbor constraints
        '''
        for new_break in self.breaks[:-1]:

            # Get neighbor break
            if new_break.neighbor_break:
                neighbor_break                 = new_break.neighbor_break
                neighbor_break.dont_break_temp = False

    def calculate_self_penalty(self):
        '''
        Calculate self penalty score
        '''
        self.self_penalty = 0
        for new_break in self.breaks[:-1]:

            # Get neighbor break
            if new_break.neighbor_break:

                neighbor_break     = new_break.neighbor_break

                if neighbor_break and neighbor_break in self.breaks:

                    # Update total penalty score
                    self.self_penalty += 1

        # Divide penalty by 2
        self.self_penalty = int(1.0*self.self_penalty/2)

    def initialize(self):
        '''
        Prepare the lists
        '''
        self.breaks = [break_path.break_node for break_path in self.break_paths[::-1]]
        self.edges  = [break_path.break_edge for break_path in self.break_paths[::-1]]
        self.scores = [break_path.score for break_path in self.break_paths[::-1]]

    def is_identical(self, other_solution, max_index=None):
        '''
        Compare the baths between two solutions
        '''
        # Determine maximum index
        max_i = len(self.breaks)-1
        if max_index:
            max_i = max_index

        # Identical variable
        identical = True

        # Breaks for the current path
        current_breaks = self.breaks[:max_i]

        # Other breaks
        other_breaks   = other_solution.breaks[:max_i]

        if not len(current_breaks) == len(other_breaks):
            return False

        # Pairwise comparison of elements
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
        self.complete        = True

    def get_cvs_rows(self):
        cvs_writer_rows = []
        for key in self.break_solutions:
            # Get break solution
            break_solution    = self.break_solutions[key]

            # If the solution doesnt exist move to next break solution
            if not break_solution:
                continue

            # Extend the list with the row objects
            cvs_writer_rows.extend(break_solution.get_cvs_rows())

        return cvs_writer_rows

    def break_group_solution(self):
        '''
        Break group solution
        '''
        for key in self.break_solutions:
            # Get break solution
            break_solution    = self.break_solutions[key]

            # If the solution doesnt exist move to next break solution
            if not break_solution:
                continue

            # Perform break
            break_solution.break_oligo_solution()

    def print_solution(self):
        '''
        Print group solution
        '''
        if self.break_solutions:
            tqdm.write('Complete:%-5s TotalScore:%-5.2f - TotalCrossoverPenalty:%-3d' %
                       (self.complete, self.total_score, self.total_penalty))
            # Print the solutions
            for oligo_key in self.break_solutions:
                # Print solution for oligo
                tqdm.write('Solution for oligo: (%d,%d,%d)' % oligo_key)
                if self.break_solutions[oligo_key]:
                    self.break_solutions[oligo_key].print_solution()

    def calculate_penalty(self, verbose=False):
        '''
        Calculate penalty for the oligo group solution
        '''
        self.total_score   = 0
        self.total_penalty = 0
        self.complete      = True

        # Iterate over each solution
        for key in self.break_solutions:
            # Get break solution
            break_solution    = self.break_solutions[key]

            # If the solution doesnt exist move to next break solution
            if not break_solution:
                if verbose:
                    tqdm.write('SOLUTION DOESNT EXIST for oligo (%d,%d,%d)' % (key[0], key[1], key[2]))
                self.complete = False
                continue

            # Update total score
            self.total_score += break_solution.score

            # Initialize the bad break list
            break_solution.bad_list = []

            # Iterate over breaks in breaks list
            for new_break in break_solution.breaks[:-1]:

                # Get neighbor key
                if new_break.neighbor_break:
                    neighbor_break     = new_break.neighbor_break
                    neighbor_oligo_key = neighbor_break.oligo.key

                    # Get the neighbor information
                    if(neighbor_oligo_key in self.break_solutions and self.break_solutions[neighbor_oligo_key] and
                       neighbor_break in self.break_solutions[neighbor_oligo_key].breaks):
                        break_solution.bad_list.append(new_break)

                        # Update total penalty score
                        self.total_penalty += 1

        # Divide penalty score by 2
        self.total_penalty = int(self.total_penalty/2)


class CompleteBreakSolution:
    '''
    Complete break solution class
    '''

    def __init__(self):

        self.temperature_celcius = None
        self.group_solutions     = {}
        self.total_prob          = 1
        self.total_score         = 0
        self.total_penalty       = 0
        self.sequence_offset     = 0
        self.complete            = True
        self.cvs_writer_rows     = []

        '''
        COLUMNS
        1.  Oligo key
        2.  Oligo Group key
        3.  Break 1 key
        4.  Break 1 type
        5.  Break 1 location
        5.  Break 2 key
        6.  Break 2 type
        7.  Break 2 location
        8.  Length
        9.  Sequence
        10. Edge weight
        11. ProbFolding
        12. Log-Probfolding
        13. Tf
        14. maxTm
        15. maxSeq
        16. has14
        17. dGtotal
        18. dGintrin
        19. dGloop
        20. dGconc
        '''
        self.cvs_header = ['OligoKey',
                           'OligoGroupkey',
                           'StartBreakKey',
                           'StartBreakType',
                           'StartBreakLocation',
                           'EndBreakKey',
                           'EndBreakType',
                           'EndBreakLocation',
                           'Length',
                           'Sequence',
                           'EdgeWeight',
                           'ProbFold',
                           'LogProbFold',
                           'Tf',
                           'maxTm',
                           'maxSeqLength',
                           'Has14',
                           'dGtotal',
                           'dGintrin',
                           'dGLoop',
                           'dGconc']

    def calculate_total_score(self):
        '''
        Calculate total score for the best solutions
        '''
        self.total_prob    = 1.0
        self.total_score   = 0
        self.total_penalty = 0
        self.complete      = True

        for key in self.group_solutions:
            # Break group solution
            if self.group_solutions[key]:
                self.total_prob    *= np.exp(self.group_solutions[key].total_score)
                self.total_score   += self.group_solutions[key].total_score
                self.total_penalty += self.group_solutions[key].total_penalty
                self.complete      *= self.group_solutions[key].complete
            else:
                self.complete       = False

    def get_cvs_rows(self):
        '''
        Get cvs writer rows
        '''
        cvs_writer_rows = []
        for key in self.group_solutions:
            # Check if the solution exists
            if self.group_solutions[key]:
                cvs_writer_rows.extend(self.group_solutions[key].get_cvs_rows())

        return cvs_writer_rows

    def get_summary_rows(self):
        '''
        Get summary rows
        '''
        summary_rows = [['Temperature', self.temperature_celcius],
                        ['TotalProb', self.total_prob],
                        ['TotalScore', self.total_score],
                        ['TotalPenalty', self.total_penalty],
                        ['Complete', int(self.complete)],
                        ['SequenceOffset', self.sequence_offset]]

        return summary_rows

    def export_staples(self, filename):
        '''
        Export the staples and its scores into an excel file
        '''
        # Get summary rows
        summary_rows = np.array(self.get_summary_rows())

        # Get the cvs rows
        cvs_rows = np.array(self.get_cvs_rows())

        # If the data arrays are empty, quit
        if len(summary_rows) == 0 or len(cvs_rows) == 0:
            return

        # Load workbook
        book = openpyxl.load_workbook(filename)
        writer      = pandas.ExcelWriter(filename, engine='openpyxl')
        writer.book = book

        # Assign sheet number
        sheet_number = self.sequence_offset

        # Create data frames
        summary_frame = pandas.DataFrame(summary_rows)

        # Create staples frames
        staples_frame  = pandas.DataFrame(cvs_rows)

        # Write summary data
        summary_frame.to_excel(writer, sheet_name=str(sheet_number), header=None, index=False)

        # Write staples data
        staples_frame.to_excel(writer, sheet_name=str(sheet_number), header=self.cvs_header, index=False, startrow=10)

        # Save writer and close
        writer.save()
        writer.close()


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

    def sort_oligos_by_length(self):
        '''
        Sort oligos by length
        '''
        self.oligos.sort(key=lambda x: x.length, reverse=True)

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
            new_group_solution = GroupBreaksolution()
            new_group_solution.break_solutions = {}

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
                      dynamic_ncols=True, bar_format='{l_bar}{bar}'):

            # Reset temporary neighbor constraints
            self.reset_temp_neighbor_constraints()

            # Initialize group break solution
            new_group_solution = GroupBreaksolution()
            new_group_solution.break_solutions = {}

            # Shuffle oligos
            if shuffle_oligos:
                self.shuffle_oligos()

            # Iterate over every oligo
            for oligo in tqdm(self.oligos, desc='Oligo loop', leave=False,
                              dynamic_ncols=True, bar_format='{l_bar}{bar}'):

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
        new_edge = BreakEdge()

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
            tqdm.write('Processing oligo:%-15s Number of breaks:%-3d' % (self.key, len(self.breaks)))

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
                                      dynamic_ncols=True, bar_format='{l_bar}{bar}'):

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
        self.oligos         = None
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
        self.sequence_offset = 0

    def circularize_scaffold(self):
        '''
        Circularize scaffold
        '''

        strand5p = self.scaffolds[0].strand5p()
        strand3p = self.scaffolds[0].strand3p()

        # Set sequence start position before circularizing
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

    def set_sequence_offset(self, key):
        '''
        Set sequence offset for scaffold from position key (vh,idx)
        '''
        # Determine the sequence offset
        self.sequence_offset = self.sequence_start_offset

        if key in self.key2scaffold:
            self.sequence_offset += (-self.key2scaffold[key][0]+1)

    def set_circularize(self, circularize=False):
        '''
        Set circularize (scaffold) parameter
        '''
        self.circularize = circularize

    def set_sequence_start_pos(self, key):
        '''
        Set sequence start position in (vh, idx)
        '''
        self.sequence_start_pos = key

    def set_sequence_start_offset(self, offset=0):
        '''
        Set sequence offset
        '''
        if offset >= 0:
            self.sequence_start_offset = offset
        else:
            self.sequence_start_offset = 0

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

        # Set sequence offset
        self.set_sequence_offset(self.sequence_start_pos)

        # Apply sequence
        self.apply_sequence(self.sequence_offset)

        # Assign strand sequences
        self.assign_strands_dna()

        # Assign scaffold positions
        self.assign_scaffold_positions()

        # Sort staple by length
        self.sort_staples_by_length()

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

    def assign_strands_dna(self):
        '''
        Assign strand sequences from cadnano part
        '''
        for oligo in self.oligos['staple']:

            # Assign current strand
            current_strand = oligo.null_strand.next_strand

            while current_strand:
                current_strand.dna = current_strand.cadnano_strand.sequence()

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
            previous_break   = BreakNode()
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
            oligo.null_break.sequence             = current_strand.sequences[0] if len(current_strand.sequences) > 0 else None
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
                    new_break   = BreakNode()

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

                # Assign group key
                new_oligo_group.key    = group_key

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
        self.oligos = {'scaffold': [], 'staple': []}

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
        self.cadnano_oligos = sorted(self.cadnano_oligos, key=lambda x: x.length(), reverse=True)

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

        # Cadnano parameters
        self.origami                         = None
        self.json_input                      = None

        # Constraints
        self.MAX_OLIGO_LENGTH                = 60
        self.UPPER_BOUND                     = 60
        self.LOWER_BOUND                     = 21

        # Local and global solutions
        self.NUM_OLIGO_SOLUTIONS             = 100
        self.NUM_GLOBAL_SOLUTIONS            = 1

        # Score parameters
        self.total_score                     = 0
        self.total_penalty                   = 0

        # k-shortest path parameter
        self.k_select                        = 'best'

        # Break rule
        self.break_rule                      = ['xstap', 'all3']

        # Maximum number of breaks allowed per oligo for autobreak
        self.max_num_breaks                  = 10000000

        # Optimization parameters
        self.optim_shuffle_oligos            = True
        self.optim_pick_method               = 'best'

        # Output file
        self.json_legacy_output             = 'out_legacy.json'
        self.json_cn25_output               = 'out_cn25.json'
        self.output_directory                = '.'

        # Optimization function
        self.optim_args                      = [['dG:50']]

        # Maxseq parameter
        self.optim_maxseq_length             = 14

        # Length gaussian parameters
        self.optim_length_mean               = 45
        self.optim_length_tolerance          = 5

        # Tm gaussian parameters
        self.optim_Tm_mean                   = 60
        self.optim_Tm_tolerance              = 5

        # Temperature parameters for dG optimization
        self.optim_temperature_celcius       = 50
        self.optim_temperature_kelvin        = self.optim_temperature_celcius+273.15

        # Max sequence gaussian parameters
        self.optim_maxseq_mean               = 14
        self.optim_maxseq_tolerance          = 2

        # Set optimization score functions
        self.optim_score_functions           = ['sum']

        # Structure factor
        self.optim_structure_factor         = 1.0

        # Set function dictionary
        self.optim_funcs_dict                = {'14': self._optimize_14,
                                                '16': self._optimize_16,
                                                'Tm': self._optimize_Tm,
                                                'dG': self._optimize_dG,
                                                'structure': self._optimize_structure,
                                                'length': self._optimize_length,
                                                'maxseq': self._optimize_maxseq,
                                                'glength': self._gauss_length,
                                                'gmaxseq': self._gauss_maxseq,
                                                'gTm': self._gauss_Tm}
        # Set optimization params
        self.optim_params_dict              = {'14': [],
                                               '16': [],
                                               'Tm': [],
                                               'dG': [self.optim_temperature_celcius],
                                               'structure': [self.optim_structure_factor],
                                               'length': [],
                                               'maxseq': [self.optim_maxseq_length],
                                               'glength': [self.optim_length_mean, self.optim_length_tolerance],
                                               'gmaxseq': [self.optim_maxseq_mean, self.optim_maxseq_tolerance],
                                               'gTm': [self.optim_Tm_mean, self.optim_Tm_tolerance]}

        # Verbose output
        self.verbose_output                 = False

        # Solutions containers
        self.complete_solutions             = {}
        self.best_complete_solution         = None

        # Sequence parameter
        self.best_sequence_offset           = 0

        # Permutation parameter
        self.permute_sequence               = False

        # Excel file that stores the results
        self.results_excel_file             = None
        self.autobreak_excel_file           = None
        self.summary_excel_file             = None

    def get_results_summary(self):
        '''
        Get results summary
        '''
        self.results_summary = []

        # Print the scores
        for complete_solution in self.sorted_complete_solutions:
            self.results_summary.append([complete_solution.sequence_offset,
                                         complete_solution.total_prob,
                                         complete_solution.total_score,
                                         complete_solution.total_penalty])

        return self.results_summary

    def write_results_summary(self):
        '''
        Write results summary
        '''

        # Get results summary
        results_summary = np.array(self.get_results_summary())

        # If the array is empty, quit
        if len(results_summary) == 0:
            sys.exit('SOLUTION DOESNT EXIST')

        # Create writer
        writer      = pandas.ExcelWriter(self.summary_excel_file, engine='openpyxl')

        # Create data frames
        summary_frame  = pandas.DataFrame(results_summary)

        # Create summary header
        summary_header = ['SequenceOffset', 'TotalProb', 'TotalScore', 'TotalPenalty']

        # Write summary data
        summary_frame.to_excel(writer, sheet_name='summary', header=summary_header, index=False)

        # Save writer and close
        writer.save()
        writer.close()

    def write_results(self, sequence_offset=0):
        '''
        Write Solution results to a sheet in results excel file
        '''
        if sequence_offset in self.complete_solutions:
            self.complete_solutions[sequence_offset].export_staples(self.results_excel_file)

    def write_best_result(self):
        '''
        Write best result
        '''
        if not self.write_all_results:
            self.best_complete_solution.export_staples(self.results_excel_file)

    def create_results_excel_file(self):
        '''
        Create resuls excel file
        '''
        wb = openpyxl.Workbook()
        wb.save(self.results_excel_file)

    def set_write_all_results(self, write_all_results=False):
        '''
        Set write all results
        '''
        self.write_all_results = write_all_results

    def set_temperature_parameter(self):
        '''
        Set temperature parameters
        '''
        self.optim_temperature_celcius = self.optim_params_dict['dG'][0]
        self.optim_temperature_kelvin  = self.optim_temperature_celcius + 273.15

    def set_permute_sequence(self, permute=False):
        '''
        Set permute sequence
        '''
        self.permute_sequence = permute

    def color_oligos_by_folding_prob(self):
        '''
        Color oligos by folding prob
        '''
        for oligo in self.origami.oligos['staple']:
            oligo.color_by_folding_prob()

    def color_oligos_by_Tf(self):
        '''
        Color oligos by folding prob
        '''
        for oligo in self.origami.oligos['staple']:
            oligo.color_by_Tf()

    def preprocess_optim_params(self):
        '''
        Preprocess optimization parameters


        RULE1. If there is no cross in break rule,
               make oligo solution and global solution number 1
               since Optimization yields the best result

        RUL2. If oligos are not shuffled, make num oligo solutions
              and global solutions 1

        '''

        # RULE 1
        if 'xscaf' not in self.break_rule and 'xstap' not in self.break_rule:
            self.NUM_OLIGO_SOLUTIONS  = 1
            self.NUM_GLOBAL_SOLUTIONS = 1

        # RULE 2
        if not self.optim_shuffle_oligos:
            self.NUM_OLIGO_SOLUTIONS  = 1
            self.NUM_GLOBAL_SOLUTIONS = 1

    def set_output_directory(self, input_filename):
        '''
        Set output directory
        '''
        # Split input file
        head, tail       = os.path.split(os.path.abspath(input_filename))
        root, ext        = os.path.splitext(tail)

        # List existing output directories
        potential_directories = list(filter(lambda x: os.path.isdir(x),
                                            glob.glob(head+'/'+root+'_autobreak_[0-9][0-9][0-9]')))

        # Get the number extensions
        number_extensions = [int(x[-3:]) for x in potential_directories]

        # Get the counter
        output_counter = 1
        if len(number_extensions) > 0:
            output_counter = max(number_extensions)+1

        self.output_directory = head+'/'+root+"_autobreak_%03d" % (output_counter)

        # Make directory
        os.mkdir(self.output_directory)

        # Copy input file to output directory
        copyfile(input_filename, self.output_directory+'/'+tail)

    def copy_sequence_file(self):
        '''
        Copy sequence file to output directory
        '''
        if self.origami.sequence_file and os.path.isfile(self.origami.sequence_file):
            head, tail = os.path.split(self.origami.sequence_file)
            copyfile(self.origami.sequence_file, self.output_directory+'/'+tail)
        elif self.origami.sequence_file not in utilities.SCAFFOLD_SEQUENCES:
            f = open(self.output_directory+'/random_sequence.txt', 'w')
            f.write(self.origami.scaffold_sequence)
            f.close()

    def write_input_args(self):
        '''
        Write input arguments in a txt file
        '''
        self.args_file = self.output_directory+'/args.txt'
        f = open(self.args_file, 'w')
        for key in self.args_dict:
            f.write('%-8s: %-8s\n' % (key, self.args_dict[key]))

        f.close()

    def define_output_files(self):
        '''
        Define json output
        '''

        # Split input file
        head, tail       = os.path.split(self.origami.json_input)
        root, ext        = os.path.splitext(tail)

        # Output files
        self.json_start_output  = self.output_directory+'/'+root+'_start.json'

        self.json_legacy_output = self.output_directory+'/'+root+'_autobreak_legacy.json'
        self.json_cn25_output   = self.output_directory+'/'+root+'_autobreak_cn25.json'

        # Results excel file
        self.results_excel_file   = self.output_directory+'/'+root+'.xlsx'
        self.autobreak_excel_file = self.output_directory+'/'+root+'_autobreak.xlsx'
        self.summary_excel_file   = self.output_directory+'/'+root+'_summary.xlsx'

    def write_part_to_json(self, filename, legacy_option=True):
        '''
        Write part to json
        '''
        self.origami.doc.writeToFile(filename, legacy=legacy_option)

    def write_start_part_to_json(self):
        '''
        Write the starting point for autobreak
        '''
        self.origami.doc.writeToFile(self.json_start_output, legacy=True)

    def write_final_part_to_json(self):
        '''
        Write cadnano part to json
        '''
        self.origami.doc.writeToFile(self.json_legacy_output, legacy=True)
        self.origami.doc.writeToFile(self.json_cn25_output, legacy=False)

    def set_pick_method(self, pick_method='random'):
        '''
        Set solution picking method
        '''
        self.optim_pick_method = pick_method

    def set_maximum_breaks(self, maximum_num_breaks=100000):
        '''
        Set maximum number of breaks parameter
        '''
        self.max_num_breaks = maximum_num_breaks

    def set_oligo_shuffle_parameter(self, oligo_shuffle_parameter=False):
        '''
        Set oligo shuffle parameter
        '''
        self.optim_shuffle_oligos = oligo_shuffle_parameter

    def set_solution_nums(self, solutions_per_oligo=1000, global_solutions=1000):
        '''
        Set solution numbers
        '''
        self.NUM_OLIGO_SOLUTIONS  = solutions_per_oligo
        self.NUM_GLOBAL_SOLUTIONS = global_solutions

    def set_k_select(self, k_parameter='best'):
        '''
        Set k-select value
        '''
        self.k_select = k_parameter

    def set_break_rule(self, new_break_rule=['xstap', 'all2']):
        '''
        Set break rule
        '''
        self.break_rule         = [rule for rule in new_break_rule]
        self.origami.break_rule = self.break_rule

    def set_verbose_output(self, verbose=False):
        '''
        Set verbose output
        '''
        self.verbose_output = verbose

    def shift_scaffold_sequence(self, offset=0):
        '''
        Shift scaffold sequence and update autobreak graphs
        '''
        # Apply the offset and shift sequence
        self.origami.apply_sequence(offset)

        # Update strand and sequence dna
        self.origami.assign_strands_dna()
        self.origami.update_sequences_dna()

        # Update graph edge weights
        self.update_edge_weights()

    def permute_scaffold_sequence(self):
        '''
        Permute scaffold sequence
        '''
        # Set final offset
        final_offset = 1
        if self.permute_sequence:
            final_offset = len(self.origami.scaffold_sequence)

        for offset in tqdm(range(0, final_offset), desc='Permutation loop', leave=False,
                           dynamic_ncols=True, bar_format='{desc}: {percentage:3.2f}%|'+'{bar}'):
            # Shift the sequence to the offset
            self.shift_scaffold_sequence(offset)

            # Run autobreak
            self.run_autobreak()

            # Write results
            if self.write_all_results:
                self.write_results(offset)

    def run_autobreak(self):
        '''
        Run basic autobreak protocol
        '''

        # Create stepwise group solutions
        self.create_stepwise_group_solutions()

        # Sort and print group solutions
        self.sort_group_solutions()

        # Combine group solutions
        self.combine_group_solutions()

    def determine_oligo_scores(self):
        '''
        Determine oligo scores from design
        '''
        # Make break-break graph
        self.origami.prepare_origami()

        # Determine initial scores
        self.determine_initial_scores()

    def create_independent_group_solutions(self):
        '''
        Create independent group solution
        No temporary neighbor constraints are imposed during run
        '''

        # Break oligos
        self.create_oligo_solutions()

        # Create group solutions
        self.combine_oligo_solutions()

    def create_stepwise_group_solutions(self):
        '''
        Main function for solution determination
        '''
        for oligo_group in tqdm(self.origami.oligo_groups, desc='Main loop ',
                                dynamic_ncols=True, bar_format='{l_bar}{bar}'):

            # Sort oligos by length
            oligo_group.sort_oligos_by_length()

            # Create solutions via stepwise approach
            oligo_group.create_stepwise_oligo_solutions(self.NUM_OLIGO_SOLUTIONS, self.NUM_GLOBAL_SOLUTIONS,
                                                        self.optim_pick_method, self.optim_shuffle_oligos,
                                                        verbose=self.verbose_output)

            # Remove incomplete solutions
            oligo_group.remove_incomplete_solutions()

    def update_edge_weights(self):
        '''
        Update edge weights
        '''
        for oligo in self.origami.oligos['staple']:
            # Visit each break object
            for current_break in oligo.breaks:
                # Iterate over the break edges
                for break_edge in current_break.break_edges:
                    # Update edge weights
                    break_edge.update_connection()

    def initialize(self):
        '''
        Initialize the connectivity maps
        '''

        for oligo in self.origami.oligos['staple']:
            # Check oligo length, if the length is within length limits dont break it
            if oligo.length < self.LOWER_BOUND:
                oligo.dont_break = True

            # Visit each break object
            for current_break in oligo.breaks:
                # Initialize the break edges
                current_break.break_edges = []

                # Get next break
                next_break = current_break.next_break
                # Iterate over the breaks
                while next_break:

                    # Determine break to break distance
                    break_distance = current_break.get_break_distance(next_break)

                    if break_distance >= self.LOWER_BOUND and break_distance <= self.UPPER_BOUND:
                        # Create break
                        new_edge = BreakEdge()

                        # Assign origami
                        new_edge.origami = self.origami

                        # Assign autobreak
                        new_edge.autobreak = self

                        # Assign edge length
                        new_edge.edge_length = break_distance

                        # Make the connection
                        new_edge.make_connection(current_break, next_break)

                        # Make loop edge
                        new_edge.make_loop_edge()

                        # Set edge weight
                        new_edge.edge_weight = self.optimize(new_edge)

                        # Add directed edge to current break's edges
                        current_break.break_edges.append(new_edge)

                        # Add break edge to edge map
                        self.origami.break_edge_map[current_break.key+next_break.key] = new_edge

                    # Stop criteria
                    if break_distance > self.UPPER_BOUND or next_break == current_break:
                        break

                    next_break = next_break.next_break

    def reset_temp_neighbor_constraints(self):
        '''
        Reset temporary neighbor constraints
        '''
        for oligo in self.origami.oligos['staple']:
            oligo.reset_temp_neighbor_constraints()

    def determine_initial_scores(self):
        '''
        Determine starting scores
        '''
        for oligo in self.origami.oligos['staple']:
            oligo.get_initial_score()

    def export_initial_scores(self):
        '''
        Export initial scores to final excel file
        '''

        # Create a dummy Complete Break Solution object
        self.final_complete_break_solution = CompleteBreakSolution()

        # Store final cvs rows
        self.final_cvs_rows = []

        # Save total score and prob
        total_score = 0
        total_prob  = 1.0
        for oligo in self.origami.oligos['staple']:
            self.final_cvs_rows.append(oligo.end_to_end_edge.get_cvs_row_object())

            # Add scores
            total_score += oligo.initial_score
            total_prob  *= oligo.folding_prob

        # Prepare summary data
        self.final_summary_data = [['TotalProb', total_prob],
                                   ['TotalScore', total_score]]

        # Check if the data arrays are empty
        if len(self.final_cvs_rows) == 0:
            return

        # Results header
        cvs_header = self.final_complete_break_solution.cvs_header

        # Write the results
        writer      = pandas.ExcelWriter(self.autobreak_excel_file,
                                         engine='openpyxl')

        # Create data frames
        summary_frame = pandas.DataFrame(np.array(self.final_summary_data))

        # Create staples frames
        staples_frame  = pandas.DataFrame(np.array(self.final_cvs_rows))

        # Write summary data
        summary_frame.to_excel(writer, sheet_name=str('final'), header=None, index=False)

        # Write staples data
        staples_frame.to_excel(writer, sheet_name=str('final'), header=cvs_header, index=False, startrow=10)

        # Save writer and close
        writer.save()
        writer.close()

    def create_oligo_solutions(self):
        '''
        Break oligos
        '''
        for oligo in self.origami.oligos['staple']:
            if not oligo.dont_break:
                oligo.generate_shortest_paths(self.NUM_OLIGO_SOLUTIONS, verbose=self.verbose_output)
                oligo.remove_penalized_solutions()

    def combine_oligo_solutions(self):
        '''
        Create oligo group solutions
        '''
        for oligo_group in self.origami.oligo_groups:
            oligo_group.combine_oligo_solutions(self.NUM_GLOBAL_SOLUTIONS)

    def sort_group_solutions(self):
        '''
        Sort group solutions based on penalty score
        '''
        for oligo_group in self.origami.oligo_groups:
            oligo_group.sort_solutions()

            # If verbose output print solutions
            if self.verbose_output:
                oligo_group.print_solutions()

    def combine_group_solutions(self):
        '''
        Combine group solutions
        '''

        # Make new complete solution
        new_complete_solution = CompleteBreakSolution()
        new_complete_solution.group_solutions = {}
        new_complete_solution.sequence_offset = self.origami.sequence_offset

        # Iterate over all the oligo groups
        for oligo_group in self.origami.oligo_groups:

            # Add the best and best penalty solutions
            new_complete_solution.group_solutions[oligo_group.key]  = oligo_group.best_score_solution

        # Calculate the total score and penalty for the complete solution
        new_complete_solution.calculate_total_score()

        # Assign the best complete solution if it is complete
        if new_complete_solution.complete:
            self.complete_solutions[self.origami.sequence_offset] = new_complete_solution

    def compare_complete_solutions(self):
        '''
        Compare complete solutions
        '''
        # Sort complete solutions
        self.sorted_complete_solutions = sorted(list(self.complete_solutions.values()),
                                                key=lambda solution: solution.total_score, reverse=True)

        # Print the scores
        for complete_solution in self.sorted_complete_solutions:
            # Print total score and crossover penalty for the best solution
            tqdm.write('Complete solutions: Offset: %-5d' % (complete_solution.sequence_offset) +
                       ' - TotalProb:%-5.2f' % (complete_solution.total_prob) +
                       ' - TotalScore:%-5.2f' % (complete_solution.total_score) +
                       ' - TotalCrossoverPenalty:%-3d' % (complete_solution.total_penalty))

        # Assign best solution
        if len(self.sorted_complete_solutions) > 0:
            self.best_complete_solution = self.sorted_complete_solutions[0]

    def set_best_sequence_offset(self):
        '''
        Set best sequence offset
        '''
        # Best sequence offset
        if self.best_complete_solution:
            self.best_sequence_offset = self.best_complete_solution.sequence_offset

        # Apply the offset and shift sequence
        self.origami.apply_sequence(self.best_sequence_offset)

    def break_best_complete_solution(self):
        '''
        Oligo breaking routine
        '''
        if self.best_complete_solution is None:
            sys.exit('SOLUTION DOESNT EXIST!')

        for key in self.best_complete_solution.group_solutions:
            # Break group solution
            if self.best_complete_solution.group_solutions[key]:
                self.best_complete_solution.group_solutions[key].break_group_solution()

    def set_score_func(self, func_args):
        '''
        Set optimization score functions
        '''
        self.optim_score_functions = func_args

    def set_optimization_func(self, func_args):
        '''
        Set optimization function
        '''
        self.optim_args         = func_args
        self.optim_args_funcs   = [function[0] for function in func_args]
        self.optim_args_params  = [[float(x) for x in function[1:]] for function in func_args]
        self.optimize_func_list = []

        # Set optimization function parameters
        for i in range(len(self.optim_args_funcs)):
            func   = self.optim_args_funcs[i]
            params = self.optim_args_params[i]

            # Make the optimize function
            if func in self.optim_funcs_dict:
                self.optimize_func_list.append(self.optim_funcs_dict[func])

            # Assign function parameters
            if func in self.optim_params_dict and len(params) <= len(self.optim_params_dict[func]):
                for j in range(len(params)):
                    self.optim_params_dict[func][j] = params[j]

    def optimize(self, edge):
        '''
        final optimization function
        '''
        score = 0
        # Get the score for each function type
        score_list = np.array([func(edge) for func in self.optimize_func_list])

        if 'sum' in self.optim_score_functions:
            score += np.sum(score_list)
        if len(score_list) > 1 and 'product' in self.optim_score_functions:
            score += np.product(score_list)
        return score

    def _optimize_structure(self, edge):
        '''
        Optimization function for structure
        '''
        return edge.edge_structure*self.optim_params_dict['structure'][0]

    def _optimize_dG(self, edge):
        '''
        Optimization function dG
        '''
        return edge.edge_logprob

    def _optimize_14(self, edge):
        '''
        Optimization function 14
        '''
        return edge.edge_has14

    def _optimize_16(self, edge):
        '''
        Optimization function 16
        '''
        return edge.edge_has16

    def _optimize_length(self, edge):
        '''
        Optimization function 16
        '''
        return edge.edge_length

    def _optimize_maxseq(self, edge):
        '''
        Optimization function for N
        '''
        return int(edge.edge_maxseq >= self.optim_params_dict['maxseq'][0])

    def _optimize_Tm(self, edge):
        '''
        Optimization function Tm
        '''
        return edge.edge_maxTm

    def _gauss_length(self, edge):
        '''
        Optimization function gauss length
        '''

        return np.exp(-(edge.edge_length -
                      self.optim_params_dict['glength'][0])**2/self.optim_params_dict['glength'][1]**2)

    def _gauss_Tm(self, edge):
        '''
        Optimization function gauss Tm
        '''

        return np.exp(-(edge.edge_maxTm-self.optim_params_dict['gTm'][0])**2/self.optim_params_dict['gTm'][1]**2)

    def _gauss_maxseq(self, edge):
        '''
        Optimization function gauss Tm
        '''

        return np.exp(-(edge.edge_maxseq -
                      self.optim_params_dict['gmaxseq'][0])**2/self.optim_params_dict['gmaxseq'][1]**2)


class BreakEdge:
    def __init__(self):
        '''
        Break edge class
        '''
        self.autobreak     = None
        self.origami       = None
        self.current_break = None
        self.next_break    = None
        self.edge_weight   = None
        self.edge_length   = None
        self.edge_maxseq   = None
        self.edge_Tm       = None
        self.LOW_TM        = 60

        # state parameters
        self.active        = True

        # Length parameters
        self.edge_has14    = None
        self.edge_num14    = None

        self.edge_has16    = None
        self.edge_num16    = None

        # Tm parameters
        self.edge_hasTm    = None
        self.edge_numTm    = None

        self.sequence_list = None

        # Loop parameter
        self.isloop        = False

    def get_cvs_row_object(self):
        '''
        Return cvs row object

        COLUMNS
        1.  Oligo key
        2.  Oligo Group key
        3.  Break 1 key
        4.  Break 1 type
        5.  Break 1 location
        6.  Break 2 key
        7.  Break 2 type
        8.  Break 2 location
        9.  Length
        10.  Sequence
        11.  Edge weight
        12.  ProbFolding
        13.  Log-Probfolding
        14.  Tf
        15. maxTm
        16. maxSeq
        17. has14
        18. dGtotal
        19. dGintrin
        20. dGloop
        21. dGconc

        '''

        # Check if the oligo group exists
        # Doesn't exist if the clustering is not performed

        oligo_group_key = -1
        if self.current_break.oligo_group:
            oligo_group_key = self.current_break.oligo_group.key

        cvs_row =      ['.'.join([str(x) for x in self.current_break.oligo.key]),
                        oligo_group_key,
                        '.'.join([str(x) for x in self.current_break.key]),
                        self.current_break.type,
                        self.current_break.location,
                        '.'.join([str(x) for x in self.next_break.key]),
                        self.next_break.type,
                        self.next_break.location,
                        self.edge_length,
                        ''.join([str(x) for x in self.ssDNA_seq_list]),
                        self.edge_weight,
                        self.edge_prob,
                        self.edge_logprob,
                        self.edge_Tf,
                        self.edge_maxTm,
                        self.edge_maxseq,
                        int(self.edge_has14),
                        self.dG_total,
                        np.sum(self.dG_intrin_list),
                        np.sum(self.dG_inter_list),
                        self.dG_conc]

        return cvs_row

    def set_edge_weight(self):
        '''
        Set edge weight
        '''
        self.edge_weight = self.origami.autobreak.optimize(self)

    def is_valid(self):
        '''
        Determine if edge is valid
        '''
        return self.active and (not self.current_break.dont_break and not self.next_break.dont_break and not
                                self.current_break.dont_break_temp and not self.next_break.dont_break_temp)

    def make_connection(self, from_break, to_break):
        '''
        Make the connection between two break nodes
        '''

        # Set the break nodes
        self.current_break = from_break
        self.next_break    = to_break

        # Assign the edge weights
        self.update_connection()

    def update_connection(self):
        '''
        Update edge weights
        '''

        # Get the temperature parameter for dG optimization
        temperature_kelvin = self.autobreak.optim_temperature_kelvin

        # Initialize sequence and dna list
        self.sequence_list  = []
        self.ssDNA_seq_list = []
        self.dsDNA_seq_list = []

        self.ssDNA_pos_list = []
        self.dsDNA_pos_list = []

        # Make the sequence list
        self.sequence_list.append(self.current_break.sequence)

        # Check if the breaks are in consecutive positions on the same sequence
        if not (self.current_break.sequence == self.next_break.sequence and
                self.current_break.direction*(self.next_break.idx - self.current_break.idx) > 0):

            # Get forward link
            next_sequence = self.current_break.sequence.next_sequence

            # Iterate over all the sequences
            while next_sequence:
                self.sequence_list.append(next_sequence)

                # Stop if the sequence is same as the final sequence
                if next_sequence == self.next_break.sequence:
                    break

                # Update next sequence
                next_sequence = next_sequence.next_sequence

        # Iterate over sequence list
        if len(self.sequence_list) == 1:
            start_point   = self.current_break.break_point_adjusted+1
            final_point   = self.next_break.break_point_adjusted+1
            dna_sequence  = self.current_break.strand.dna[start_point:final_point]
            position_list = self.current_break.strand.scaffoldPos[start_point:final_point]

            if len(position_list) > 0:
                self.ssDNA_pos_list.append(position_list)
            if len(dna_sequence) > 0:
                self.ssDNA_seq_list.append(dna_sequence)
        else:
            # 1. Get the 5' sequence
            start_point   = self.current_break.break_point_adjusted+1
            final_point   = self.current_break.sequence.strHigh+1
            dna_sequence  = self.current_break.strand.dna[start_point:final_point]
            position_list = self.current_break.strand.scaffoldPos[start_point:final_point]

            if len(position_list) > 0:
                self.ssDNA_pos_list.append(position_list)
            if len(dna_sequence) > 0:
                self.ssDNA_seq_list.append(dna_sequence)

            # 2. Get the sequences in between
            for sequence in self.sequence_list[1:-1]:
                self.ssDNA_pos_list.append(sequence.scaffoldPos)
                self.ssDNA_seq_list.append(sequence.dna)

            # 3. Get the 3' sequence
            start_point   = self.next_break.sequence.strLow
            final_point   = self.next_break.break_point_adjusted+1
            dna_sequence  = self.next_break.strand.dna[start_point:final_point]
            position_list = self.next_break.strand.scaffoldPos[start_point:final_point]

            if len(position_list) > 0:
                self.ssDNA_pos_list.append(position_list)
            if len(dna_sequence) > 0:
                self.ssDNA_seq_list.append(dna_sequence)

        # Remove empty sequences
        self.dsDNA_seq_list = [dna.strip() for dna in self.ssDNA_seq_list if len(dna.strip()) > 0]
        self.dsDNA_pos_list = [list(filter(lambda x: x, pos_list)) for pos_list in self.ssDNA_pos_list]

        # Replace all empty characters in ssDNA seq list with ?
        self.ssDNA_seq_list = [dna.replace(" ", "?") for dna in self.ssDNA_seq_list]

        # Remove empty lists from positions list
        self.dsDNA_pos_list = list(filter(lambda x: len(x), self.dsDNA_pos_list))

        # Determine mean values for the positions
        self.dsDNA_mean_pos_list = [np.round(np.mean(x)) for x in self.dsDNA_pos_list]

        # Determine Tm
        self.Tm_list       = np.array([utilities.sequence_to_Tm(dna) for dna in self.dsDNA_seq_list])

        # 1. Determine intrinsic free energies
        self.dG_intrin_list = []
        self.dH_intrin_list = []
        self.dS_intrin_list = []

        for dna in self.dsDNA_seq_list:
            (dGintrin,
             dHintrin,
             dSintrin) = utilities.sequence_to_dG(dna, temperature_kelvin)

            # Add free energies to the lists
            self.dG_intrin_list.append(dGintrin)
            self.dH_intrin_list.append(dHintrin)
            self.dS_intrin_list.append(dSintrin)

        # Make the lists arrays
        self.dG_intrin_list = np.array(self.dG_intrin_list)
        self.dH_intrin_list = np.array(self.dH_intrin_list)
        self.dS_intrin_list = np.array(self.dS_intrin_list)

        # Get scaffold length
        scaffold_length   = self.origami.oligos['scaffold'][0].length
        scaffold_circular = self.origami.oligos['scaffold'][0].circular

        # 2. Determine interfacial coupling energies
        self.dG_inter_list = []
        self.dS_inter_list = []

        for i in range(len(self.dsDNA_mean_pos_list)-1):
            dGinter, dSinter = utilities.position_to_loop_dG(self.dsDNA_mean_pos_list[i],
                                                             self.dsDNA_mean_pos_list[i+1],
                                                             scaffold_length, scaffold_circular,
                                                             temperature_kelvin)
            self.dG_inter_list.append(dGinter)
            self.dS_inter_list.append(dSinter)

        # Make the lists arrays
        self.dG_inter_list = np.array(self.dG_inter_list)
        self.dS_inter_list = np.array(self.dS_inter_list)

        # 3. Get free energy due to concentration
        dGconc, dSconc = utilities.conc_to_dG(temperature_kelvin)

        self.dG_conc = dGconc
        self.dS_conc = dSconc

        # Get RT values
        self.RT = 0.593/298.15*temperature_kelvin

        # 4. Get total energies
        self.dG_total = np.sum(self.dG_intrin_list) + np.sum(self.dG_inter_list) + self.dG_conc
        self.dS_total = np.sum(self.dS_intrin_list) + np.sum(self.dS_inter_list) + self.dS_conc
        self.dH_total = np.sum(self.dH_intrin_list)

        # 5. Determine probabilities
        self.edge_prob = np.exp(-self.dG_total/self.RT)/(1.0+np.exp(-self.dG_total/self.RT))

        # 6. Determine log-probabilities
        self.edge_logprob = np.log(self.edge_prob)

        # 7. Determine Tf in Celcius
        self.edge_Tf = self.dH_total/self.dS_total - 273.15

        # Determine lengths
        self.ssDNA_length_list = np.array([len(dna) for dna in self.ssDNA_seq_list])
        self.dsDNA_length_list = np.array([len(dna) for dna in self.dsDNA_seq_list])

        # Determine the edge weights
        self.edge_maxTm  = max(self.Tm_list)

        # Length parameters
        self.edge_maxseq = max(self.dsDNA_length_list)
        self.edge_num14  = np.sum(self.dsDNA_length_list >= 14)
        self.edge_has14  = int(self.edge_num14 > 0)

        self.edge_num16  = np.sum(self.dsDNA_length_list >= 16)
        self.edge_has16  = int(self.edge_num16 > 0)

        # Tm parameters
        self.edge_numTm  = np.sum(self.Tm_list >= self.LOW_TM)
        self.edge_hasTm  = np.sum(self.edge_numTm > 0)

        # Determine structure score
        self.edge_num_cross = len(self.dsDNA_seq_list)
        self.edge_structure = self.edge_num_cross**2.0

        # Set edge weight
        self.set_edge_weight()

    def make_loop_edge(self):
        '''
        Make loop edge
        '''
        if self.current_break == self.next_break:
            self.isloop                  = True
            self.current_break.loop_edge = self


class BreakPath:
    def __init__(self, break_node, break_edge=None, score=0):
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

        # Nucleotide parameters
        self.current_nucleotide = None
        self.next_nucleotide    = None

        # State parameters
        self.visited          = False
        self.active           = True
        self.break_state      = None   # broken for break,  not-broken for no break
        self.dont_break       = False  # If set to True, keep the break not-broken
        self.dont_break_temp  = False  # Transient version of dont_break. If set to True, keep the break not-broken

        # Cluster info
        self.oligo_group      = None

        # Graph parameter
        self.score               = -utilities.INFINITY
        self.best_path_node      = None
        self.best_path_nodes     = None
        self.shortest_paths      = None
        self.k_potential_paths   = None
        self.k_shortest_paths    = None
        self.traverse_path       = None
        self.shortest_score      = 0
        self.order_id            = None

    def is_dsDNA(self):
        '''
        Determine if the break node is located on dsDNA
        '''
        # Initialize dsDNA parameters
        current_dsDNA = True
        next_dsDNA = True

        if self.current_nucleotide:
            current_dsDNA = self.current_nucleotide.dsDNA
        else:
            current_dsDNA = False

        if self.next_nucleotide:
            next_dsDNA = self.next_nucleotide.dsDNA
        else:
            next_dsDNA = False

        return current_dsDNA and next_dsDNA

    def reset_break_path(self):
        self.order_id         = -1
        self.traverse_path    = []
        self.best_path_nodes  = []
        self.best_path_node   = None
        self.score            = -utilities.INFINITY
        self.visited          = False
        self.shortest_paths   = []

    def reset_k_paths(self):
        self.k_shortest_paths  = []
        self.k_potential_paths = []

    def is_break_edge_possible(self, other_break):
        '''
        Check if break edge is possible from self to to another break node for a circular oligo
        '''
        order_difference = self.get_break_order_difference(other_break)

        # Criteria for a proper connection (no loop) in the right direction
        return order_difference > 0 or (not order_difference == 0 and order_difference == -self.order_id)

    def get_loop_edge(self):
        '''
        Return loop edge
        '''

        return self.loop_edge

    def get_break_distance(self, other_break):
        '''
        Break to break distance
        '''
        distance = other_break.distance - self.distance
        if distance <= 0 and self.oligo.circular:
            return distance+self.oligo.length
        else:
            return distance

    def get_break_order_difference(self, other_break):
        '''
        Return break to order id difference
        '''
        return other_break.order_id - self.order_id

    def get_valid_edge_nodes(self):
        '''
        Get break nodes connected by edges
        '''
        # Get the connected break nodes
        return [break_edge.next_break for break_edge in self.break_edges if not break_edge.next_break.dont_break]

    def get_valid_edges(self):
        '''
        Get edges that lead to break nodes that can be broken
        '''

        # Get the connected break nodes
        return [break_edge for break_edge in self.break_edges if break_edge.is_valid()]

    def get_k_shortest_paths(self, final_break, k_num=10, k_select='best'):
        '''
        Get k-shortest path results
        '''
        # Initialize k-shortest paths
        self.k_shortest_paths = []

        # 1. Get the shortest paths
        shortest_path = self.get_shortest_path(final_break)

        # If the here is no path found return empty list
        if shortest_path is None:
            return self.k_shortest_paths

        # 2.Add best path to k-path list
        self.k_shortest_paths  = [shortest_path]
        self.k_potential_paths = []
        num_k_solutions        = 1

        # Check the final score of the path
        if shortest_path.score == 0:
            return self.k_shortest_paths

        # 3. Make the paths
        while num_k_solutions < k_num:

            # Get the last best path
            last_solution = self.k_shortest_paths[-1]

            # Iterate through the edges
            for i in range(len(last_solution.edges[:-1])):
                # Inactive edge list
                inactive_edges = []

                for j in range(len(self.k_shortest_paths)):

                    if last_solution.is_identical(self.k_shortest_paths[j], i):
                        inactive_edges.append(last_solution.edges[i])

                        # Make the edge inactive
                        last_solution.edges[i].active = False

                # Reset break paths and find solution
                self.oligo.reset_break_paths()
                self.oligo.reset_break_order_ids(self, final_break)

                # Get sub solution
                potential_solution = self.get_shortest_path(final_break)

                # Make the edges active again
                for edge in inactive_edges:
                    edge.active = True

                # If there is no solution continue
                if potential_solution is None:
                    continue

                # Get the first solution in the list
                new_solution = potential_solution

                # Add potential solution if it doesnt exist in k-shortest paths
                potential_solution_exists = False
                for solution in self.k_shortest_paths:
                    if solution.is_identical(new_solution):
                        potential_solution_exists = True
                        break

                # Add new solution to potential paths
                if not potential_solution_exists:
                    self.k_potential_paths.append(new_solution)

            # Check if the potential path list is empty, if empty quit
            if len(self.k_potential_paths) == 0:
                break

            # Sort the potential paths and add the best one to k-shortest paths list
            self.k_potential_paths.sort(key=lambda x: x.score, reverse=True)

            if k_select == 'best':
                solution_index = 0
            else:
                # Add random one to the k-shortest paths
                solution_index = random.randint(0, len(self.k_potential_paths)-1)

            # Add item to k-shortest path list
            self.k_shortest_paths.append(self.k_potential_paths[solution_index])

            # Update num solutions
            num_k_solutions += 1

            # Remove the result from potential paths list
            self.k_potential_paths.pop(solution_index)

        return self.k_shortest_paths

    def get_shortest_path(self, final_break):
        '''
        Find the shortest path between current and final break points
        '''
        # Initialize the set and stack
        stack = [self]

        while stack:
            # Pop the break node
            new_break = stack.pop(0)

            # If current node is final break, quit
            if final_break.visited and new_break == final_break:
                break

            # Get valid break edges
            valid_break_edges = new_break.get_valid_edges()

            # Update the scores for connected breaks
            for break_edge in valid_break_edges:
                # If id difference is in wrong direction and if it is a loop, discard the edge
                if not new_break.is_break_edge_possible(break_edge.next_break):
                    continue

                # Determine the new score
                if new_break.order_id == 0:
                    new_score = break_edge.edge_weight
                else:
                    new_score = new_break.score + break_edge.edge_weight

                # Update the score based on the existence of a neighbor crossover
                if new_break.neighbor_break in new_break.best_path_nodes:
                    new_score += -utilities.INFINITY*utilities.INFINITY

                # Make new break Path object
                new_break_path = BreakPath(new_break, break_edge, new_score)

                # If new score is higher than the previous one, make a new list
                if new_score > break_edge.next_break.score:
                    break_edge.next_break.best_path_node  = new_break_path
                    break_edge.next_break.best_path_nodes = new_break.best_path_nodes + [new_break]
                    break_edge.next_break.score           = new_score

                # Add next break to connected breaks list
                if not break_edge.next_break.visited:
                    stack.append(break_edge.next_break)

                # Make the next break visited
                break_edge.next_break.visited = True

        # Finally compare with the loop connection
        if self == final_break and self.loop_edge and self.loop_edge.is_valid():
            # Make a loop break path object
            loop_break_path = BreakPath(self, self.loop_edge, self.loop_edge.edge_weight)

            if self.loop_edge.edge_weight > final_break.score:
                final_break.best_path_node = loop_break_path
                final_break.score          =  self.loop_edge.edge_weight

        # Return best path
        final_break.shortest_path = final_break.traverse_best_path(self)

        return final_break.shortest_path

    def traverse_best_path(self, start_break):

        # Initialize shortest path
        self.shortest_path = None

        # Check if the final break has a path node
        if self.best_path_node is None:
            return self.shortest_path

        # Make the traverse path for starting node empty
        self.traverse_path = []

        # Make final node break path object
        final_break_path = BreakPath(self, None, self.score)

        # Initalize the traverse path for the initial node
        self.best_path_node.break_node.traverse_path = self.traverse_path + [final_break_path]

        # Assign new break path
        new_break_path = self.best_path_node

        # Assign new break
        new_break = new_break_path.break_node

        while new_break:

            if new_break == start_break:
                new_break_solution              = OligoBreakSolution()
                new_break_solution.start_break  = start_break
                new_break_solution.final_break  = self
                new_break_solution.break_paths  = new_break.traverse_path + [new_break_path]
                new_break_solution.score        = self.score

                # Initialize the solution
                new_break_solution.initialize()

                self.shortest_path = new_break_solution
                break

            # Get the new path nodes
            new_path_node = new_break.best_path_node

            # Check if there is a path from new break in reverse direction
            if new_path_node is None:
                break

            # Iterate through each node
            new_path_node.break_node.traverse_path = new_break.traverse_path + [new_break_path]

            # Update new break path
            new_break_path = new_path_node

            # Update new break
            new_break = new_path_node.break_node

        return self.shortest_path

    def get_connected_breaks(self):
        '''
        Return connected breaks
        '''
        # Initialize connected breaks
        self.connected_breaks = []

        # 1. Add the directly connected next break
        if self.next_break:
            self.connected_breaks.append(self.next_break)

        # 2. Add the directly connected previous break
        if self.previous_break:
            self.connected_breaks.append(self.previous_break)

        # 3. Add the neighbor break
        if self.neighbor_break:
            self.connected_breaks.append(self.neighbor_break)

        return self.connected_breaks

    def break_cadnano(self):
        '''
        Break the breaknode
        '''
        if self.type == 'crossover':
            self.origami.remove_cadnano_crossover(self.vh, self.idx, self.direction)
        else:
            self.origami.split_cadnano_strand(self.vh, self.idx, self.direction)

    def depth_first_search(self):
        '''
        Depth-first-search graph traverse algorithm to find connected components
        '''
        # Initialize the set and stack
        visited, stack = set(), [self]
        while stack:
            # Pop the break node
            new_break = stack.pop()
            if new_break not in visited:
                # Add break to visited list
                visited.add(new_break)

                # Make the break visited
                new_break.visited = True

                # Get the connected breaks
                connected_breaks = new_break.get_connected_breaks()

                # Extend the group with new breaks
                stack.extend(set(connected_breaks) - visited)

        return list(visited)


# Parse functions

def parse_break_rule(break_rule):
    '''m
    Parse break rule
    '''
    return break_rule.split('.')


def parse_score_function(score_function):
    '''m
    Parse break rule
    '''
    return score_function.split('.')


def parse_optim_function(function_input):
    '''
    Parse optimizatiotion function input
    '''
    # Get the function groups
    groups = function_input.split('.')

    # Get the functions and its parameters
    functions = [group.split(':') for group in groups]

    return functions


def parse_sequence_position(key_input):
    '''
    Parse sequence position
    '''
    return tuple([int(x) for x in key_input.split('.')])


def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i",   "--input",    type=str,
                        help="Cadnano json file")

    parser.add_argument("-read",   "--read",  action='store_true',
                        help="Read-only to determine oligo scores")

    parser.add_argument("-rule",   "--rule",     type=str, default='xstap',
                        help="Break rule")

    parser.add_argument("-score",   "--score",     type=str, default='sum',
                        help="Optimization score function")

    parser.add_argument("-func",   "--func",     type=str, default='dG:50',
                        help="Optimization function")

    parser.add_argument("-seq",   "--sequence", type=str, default=None,
                        help="Sequence file in txt")

    parser.add_argument("-pos",   "--position", type=str, default=None,
                        help="Sequence start position")

    parser.add_argument("-offset", "--offset", type=int, default=None,
                        help="Sequence offset")

    parser.add_argument("-sort", "--sort", action='store_true',
                        help="Sort oligos for stepwise optimization.")

    parser.add_argument("-nsol",   "--nsol",     type=int,
                        help="Number of solutions", default=10)

    parser.add_argument("-dontb",   "--dontbreak",  type=int,
                        help="Dont break oligos less than the length specified", default=0)

    parser.add_argument("-v",   "--verbose",  action='store_true',
                        help="Verbose output")

    parser.add_argument("-circularize", "--circularize",  action='store_true',
                        help="Circularize scaffold")

    parser.add_argument("-permute",   "--permute",  action='store_true',
                        help="Permute sequence")

    parser.add_argument("-writeall",   "--writeall",  action='store_true',
                        help="Write all results")

    parser.add_argument("-seed",   "--seed",  type=int, default=0,
                        help="Random seed")

    args = parser.parse_args()

    # Check if the required arguments are passed to the code
    if args.input is None:
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Assign the parameters
    input_filename          = args.input
    sequence_filename       = args.sequence
    read_only               = args.read
    break_rule              = parse_break_rule(args.rule)
    score_func              = parse_score_function(args.score)
    optimization_func       = parse_optim_function(args.func)
    global_solutions        = args.nsol
    dontbreak_less_than     = args.dontbreak
    verbose_output          = args.verbose
    random_seed             = args.seed
    permute_sequence        = args.permute
    write_all_results       = args.writeall
    start_offset            = args.offset
    shuffle_oligos          = not args.sort
    circularize             = args.circularize

    # Create args dictionary
    args_dict = {'input': args.input,
                 'sequence': args.sequence,
                 'read': args.read,
                 'rule': args.rule,
                 'score': args.score,
                 'func': args.func,
                 'nsol': args.nsol,
                 'dontb': args.dontbreak,
                 'verbose': args.verbose,
                 'seed': args.seed,
                 'permute': args.permute,
                 'writeall': args.writeall,
                 'sort': args.sort,
                 'circularize': args.circularize}

    # Check if offset and/or offset is provided
    if args.position or args.offset:
        permute_sequence = False

    # Initialize start position and sequence offsets
    start_offset = -1
    start_pos    = (-1, -1)
    if args.offset:
        start_offset = args.offset
    if args.position:
        start_pos    = parse_sequence_position(args.position)

    # Check if input file exists
    if not os.path.isfile(input_filename):
        sys.exit('Input file does not exist!')

    # Set random seed in order to get the same results with the same set of parameters
    random.seed(random_seed)

    # Create Origami object
    new_origami = Origami()

    # Create Autobreak object and make cross assignments
    new_autobreak = AutoBreak()

    # Store args
    new_autobreak.args_dict = args_dict

    # Cross assignments
    new_origami.autobreak = new_autobreak
    new_autobreak.origami = new_origami

    # Assign break rule
    new_autobreak.set_break_rule(break_rule)

    # Set solution numbers
    new_autobreak.set_solution_nums(1, global_solutions)

    # Set optimization  functions
    new_autobreak.set_optimization_func(optimization_func)

    # Set score functions
    new_autobreak.set_score_func(score_func)

    # Set permute sequence
    new_autobreak.set_permute_sequence(permute_sequence)

    # Set oligo shuffle parameter
    new_autobreak.set_oligo_shuffle_parameter(shuffle_oligos)

    # Preprocess optimization parameters
    new_autobreak.preprocess_optim_params()

    # Set verbose parameter
    new_autobreak.set_verbose_output(verbose_output)

    # Set output directory
    new_autobreak.set_output_directory(input_filename)

    # Set write all flag
    new_autobreak.set_write_all_results(write_all_results)

    # Set optimization temperature
    new_autobreak.set_temperature_parameter()

    # Initialize origami object
    new_origami.initialize(input_filename)

    # Set sequence filename
    new_origami.set_sequence_file(sequence_filename)

    # Set sequence start position
    new_origami.set_sequence_start_pos(start_pos)

    # Set sequence start offset
    new_origami.set_sequence_start_offset(start_offset)

    # Set circularize
    new_origami.set_circularize(circularize)

    # Define json output
    new_autobreak.define_output_files()

    # Write input arguments
    new_autobreak.write_input_args()

    # Check if it is a read-only or autobreak run
    if not read_only:
        # Prepare origami for autobreak
        new_origami.prepare_origami()

        # Cluster staples
        new_origami.cluster_oligo_groups()

        # Set dont break flag for oligos
        new_origami.set_dont_break_oligos(dontbreak_less_than)

        # Make break-break graph
        new_autobreak.initialize()

        # Create results excel file
        new_autobreak.create_results_excel_file()

        # Run the permutation protocol
        new_autobreak.permute_scaffold_sequence()

        # Compare complete solutions
        new_autobreak.compare_complete_solutions()

        # Write results summary
        new_autobreak.write_results_summary()

        # Write best result
        new_autobreak.write_best_result()

        # Break best solutions
        new_autobreak.break_best_complete_solution()

    # Read only the scores
    new_autobreak.determine_oligo_scores()

    # Color oligos by folding prob
    new_autobreak.color_oligos_by_Tf()

    # Export initial scores
    new_autobreak.export_initial_scores()

    # Write result to json
    new_autobreak.write_final_part_to_json()

    # Write sequence file to output directory
    new_autobreak.copy_sequence_file()


if __name__ == "__main__":
    main()
