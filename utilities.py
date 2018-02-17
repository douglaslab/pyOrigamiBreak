#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2018-01-08 20:33:39
# @Author  : Tural Aksel (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import numpy as np

# CONSTANTS
INFINITY = 1E20

# Sequence enthalpy/entropy values at 37C-1M NaCl from SANTALUCIA-HICKS paper
# https://doi.org/10.1146/annurev.biophys.32.110601.141800 - DOI:10.1146/annurev.biophys.32.110601.141800

# Units are in kcal/mol for Enthalpy and cal/mol.Kelvin

# Thermodynamics constants
SCAF_CONC  = 10E-9
STAP_CONC  = 100E-9

MAGNESIUM  = 12.5E-3
TRIS       = 40E-3
R          = 0.593/298.15

TUPLE_ENTHALPY = {'AA': -7.6,
                  'TT': -7.6,

                  'AT': -7.2,
                  'TA': -7.2,

                  'CA': -8.5,
                  'TG': -8.5,

                  'GT': -8.4,
                  'AC': -8.4,

                  'CT': -7.8,
                  'AG': -7.8,

                  'GA': -8.2,
                  'TC': -8.2,

                  'CG': -10.6,
                  'GC': -9.8,

                  'GG': -8.0,
                  'CC': -8.0,

                  'INIT': 0.2,
                  'PENL': 2.2}

TUPLE_ENTROPY  = {'AA': -21.3,
                  'TT': -21.3,

                  'AT': -20.4,
                  'TA': -21.3,

                  'CA': -22.7,
                  'TG': -22.7,

                  'GT': -22.4,
                  'AC': -22.4,

                  'CT': -21.0,
                  'AG': -21.0,

                  'GA': -22.2,
                  'TC': -22.2,

                  'CG': -27.2,
                  'GC': -24.4,

                  'GG': -19.9,
                  'CC': -19.9,

                  'INIT': -5.7,
                  'PENL':  6.9}


def generate_random_sequence(oligo_length):
    '''
    Generate random DNA sequence
    '''
    return ''.join(np.random.choice(['A', 'C', 'G', 'T'], int(oligo_length)))


def generate_nA(oligo_length):
    '''
    Generate AAAAA...AAA
    '''
    return 'A'*int(oligo_length)


def generate_nC(oligo_length):
    '''
    Generate CCCCC...CCC
    '''
    return 'C'*int(oligo_length)

def get_min_distance(from_position, to_position, oligo_length):
    '''
    Get min distance between two positions on a circular oligo
    '''
    first_distance  = (to_position - from_position) % oligo_length 
    second_distance = oligo_length - first_distance

    return min(first_distance, second_distance) 

def end_to_end_distance(num_bases):
    '''
    Estimate end to end distance - Units are in nm
    Model from https://www.sciencedirect.com/science/article/pii/S0022283698918307?via%3Dihub
    https://doi.org/10.1006/jmbi.1998.1830

    L_base and L_p from https://www.nature.com/articles/nature14860 (doi:10.1038/nature14860)
    '''

    L_base    = 0.6
    L_contour = num_bases*L_base

    L_p       = 0.9

    return 2*L_p*L_contour*(1-1.0*L_p/L_contour*(1-np.exp(-1.0*L_contour/L_p)))


def distance_to_loop_dG(distance_square):
    '''
    Calculate loop entropy
    distance square is in nm2
    '''
    # Equation from https://www.nature.com/articles/nature14860 (doi:10.1038/nature14860)
    # Divide by 1000 to convert cal/(mol.K) to kcal/(mol.K)
    effective_concentration = 1.0/6.02E23*(3.0/(2*np.pi*distance_square*1E-18))**(3/2.0)/1000
    dSloop = R*np.log(effective_concentration)
    dG37   = -310.15*dSloop
    dG50   = -323.15*dSloop

    return (dG37, dG50, dSloop)

def position_to_loop_dG(from_position, to_position, oligo_length):
    # Get number of minimum number of bases between to location
    base_distance   = get_min_distance(from_position, to_position, oligo_length)

    #Determine end to end distance^2
    distance_square = end_to_end_distance(base_distance)

    # Get the free energies
    energies = distance_to_loop_dG(distance_square)

    return energies

def sequence_to_Tm(sequence):
    '''
    Convert sequence to Tm
    '''

    if len(sequence) == 0:
        return None

    # Count the number of A:T on the ends
    end_AT = (sequence[0] == 'T' or sequence[0] == 'A') + (sequence[-1] == 'T' or sequence[-1] == 'A')

    # If sequence length is 1, check only the first base
    if len(sequence) == 1:
        end_AT = (sequence[0] == 'T' or sequence[0] == 'A')

    dH     = (TUPLE_ENTHALPY['INIT'] +
              sum([TUPLE_ENTHALPY[sequence[i:i+2]] for i in range(len(sequence)-1)]) +
              end_AT*TUPLE_ENTHALPY['PENL'])
    dS    = (TUPLE_ENTROPY['INIT'] +
             sum([TUPLE_ENTROPY[sequence[i:i+2]] for i in range(len(sequence)-1)]) +
             end_AT*TUPLE_ENTROPY['PENL'])

    # Equation from Santalucia paper (https://doi.org/10.1146/annurev.biophys.32.110601.141800)
    # DOI:10.1146/annurev.biophys.32.110601.141800
    Tm_1MNaCl = 1.0*dH/(R*np.log(STAP_CONC-0.5*SCAF_CONC)+dS/1000) - 273.15

    # Fraction of G/C
    fGC    = 1.0*sum([(c == 'G' or c == 'C') for c in sequence])/len(sequence)

    # Equation 16 from salt correction paper:
    # Richard Owczarzy et al, Biochemistry, 2008 (http://pubs.acs.org/doi/abs/10.1021/bi702363u - DOI:10.1021/bi702363u)
    ln_Mg     = np.log(MAGNESIUM)
    Nbp       = len(sequence)

    a = 3.92E-5
    b = -9.11E-6
    c = 6.26E-5
    d = 1.42E-5
    e = -4.82E-4
    f = 5.25E-4
    g = 8.31E-5

    # If the sequence is less than 2 bp long, dont incorporate Mg
    if Nbp > 1:
        Tm_Mg_inv = 1.0/Tm_1MNaCl + a + b*ln_Mg + fGC*(c+d*ln_Mg) + 1.0/(2*(Nbp-1))*(e+f*ln_Mg+g*ln_Mg**2)
        Tm_Mg     = 1.0/Tm_Mg_inv
    else:
        Tm_Mg = Tm_1MNaCl

    return Tm_Mg

def conc_to_dG():
    # Entropy factor due to concenteration
    dSconc = R*np.log(STAP_CONC-0.5*SCAF_CONC)
    dG37 = -310.15*dSconc
    dG50 = -323.15*dSconc

    return (dG37, dG50, dSconc)

def sequence_to_dG(sequence):
    '''
    Convert sequence to dG
    '''
    if len(sequence) == 0:
        return None

    # Count the number of A:T on the ends
    end_AT = (sequence[0] == 'T' or sequence[0] == 'A') + (sequence[-1] == 'T' or sequence[-1] == 'A')

    # If sequence length is 1, check only the first base
    if len(sequence) == 1:
        end_AT = (sequence[0] == 'T' or sequence[0] == 'A')

    dHtotal = (TUPLE_ENTHALPY['INIT'] +
               sum([TUPLE_ENTHALPY[sequence[i:i+2]] for i in range(len(sequence)-1)]) +
               end_AT*TUPLE_ENTHALPY['PENL'])
    dSbase  = (TUPLE_ENTROPY['INIT'] +
               sum([TUPLE_ENTROPY[sequence[i:i+2]] for i in range(len(sequence)-1)]) +
               end_AT*TUPLE_ENTROPY['PENL'])

    # Determine salt correction
    # Salt correction from https://www.nature.com/articles/nature14860(doi:10.1038/nature14860)
    dSsalt  = 0.368*np.log(0.5*TRIS + 3.3*np.sqrt(MAGNESIUM))*(len(sequence)-1)

    # Get total entropy
    dStotal = (dSbase + dSsalt)/1000.0

    # Determine dG at two temperatures
    dG37 = dHtotal - 310.15*dStotal
    dG50 = dHtotal - 323.15*dStotal

    return (dG37, dG50, dHtotal, dStotal)
