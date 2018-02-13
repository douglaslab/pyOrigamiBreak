#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2018-01-08 20:33:39
# @Author  : Tural Aksel (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import numpy as np
import contextlib

# CONSTANTS
INFINITY = 1E6

# Sequence enthalpy/entropy values at 37C-1M NaCl from SANTALUCIA-HICKS paper
# https://doi.org/10.1146/annurev.biophys.32.110601.141800 - DOI:10.1146/annurev.biophys.32.110601.141800

# Units are in kcal/mol for Enthalpy and cal/mol.Kelvin

# Thermodynamics constants
SCAF_CONC  = 10E-9
STAP_CONC  = 100E-9

MAGNESIUM  = 10E-3
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


def sequence_to_Tm(sequence):
    '''
    Convert sequence to Tm
    '''

    if len(sequence) == 0:
        return None
    elif len(sequence) == 1:
        return -273.15

    # Count the number of A:T on the ends
    end_AT = (sequence[0] == 'T' or sequence[0] == 'A') + (sequence[-1] == 'T' or sequence[-1] == 'A')

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

    Tm_Mg_inv = 1.0/Tm_1MNaCl + a + b*ln_Mg + fGC*(c+d*ln_Mg) + 1.0/(2*(Nbp-1))*(e+f*ln_Mg+g*ln_Mg**2)
    Tm_Mg     = 1.0/Tm_Mg_inv

    return Tm_Mg