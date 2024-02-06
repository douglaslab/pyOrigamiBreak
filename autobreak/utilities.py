import numpy as np

# Units are in kcal/mol for Enthalpy and [cal/mol•°Kelvin]

# Physical constants
INFINITY     = 1E20     # Close enough
R_SI_units   = 8.31447  # [J/mol•°K]
cal_per_J    = 1/4.184  # [cal/J]
kcal_per_cal = 0.001    # [kcal/cal]
R            = R_SI_units*cal_per_J*kcal_per_cal  # 0.0019872 [kcal/mol•°K]

# Default concentrations
SCAF_CONC    =  10.0E-9  #  10.0 nM (Typical folding reaction)
STAP_CONC    = 100.0E-9  # 100.0 nM (Typical folding reaction)
Mg_CONC      =  12.5E-3  #  12.5 mM (Assumed in Dunn et al 2015)
TRIS_CONC    =  40.0E-3  #  40.0 mM (Assumed in Dunn et al 2015)

# Sequence enthalpy/entropy values at 37°C and 1 M NaCl from SantaLucia & Hicks 2004
# https://doi.org/10.1146/annurev.biophys.32.110601.141800 
SantaLucia2004Table1 = {
    'dH': {  # Nearest-neighbor enthalpy values ∆H° [kcal/mol]
        'AA':  -7.6, 'TT': -7.6, 'AT': -7.2, 'TA': -7.2,
        'CA':  -8.5, 'TG': -8.5, 'GT': -8.4, 'AC': -8.4,
        'CT':  -7.8, 'AG': -7.8, 'GA': -8.2, 'TC': -8.2,
        'CG': -10.6, 'GC': -9.8, 'GG': -8.0, 'CC': -8.0,
        'Initiation': 0.2, 'Terminal-AT-Penalty': 2.2
    },
    'dS': {  # Nearest-neighbor entropy values ∆S° [cal/mol•°K]
        'AA': -21.3, 'TT': -21.3, 'AT': -20.4, 'TA': -21.3,
        'CA': -22.7, 'TG': -22.7, 'GT': -22.4, 'AC': -22.4,
        'CT': -21.0, 'AG': -21.0, 'GA': -22.2, 'TC': -22.2,
        'CG': -27.2, 'GC': -24.4, 'GG': -19.9, 'CC': -19.9,
        'Initiation': -5.7, 'Terminal-AT-Penalty': 6.9
    }
}

# Count number of 'A' or 'T' chars at the start and end of sequence
def count_Terminal_ATs(sequence):
    count = 0
    if sequence[0] in ['A', 'T']:  # First char
        count += 1
    if len(sequence) > 1 and sequence[-1] in ['A', 'T']:  # Last char
        count += 1
    return count

# Calculate ∆H, per SantaLucia 2004, Equation 1 and Table 1
def get_dH_SantaLucia2004(sequence):
    NN_dH_table = SantaLucia2004Table1['dH']
    end_AT = count_Terminal_ATs(sequence)
    dH = (
        NN_dH_table['Initiation'] +
        sum([NN_dH_table[sequence[i:i+2]] for i in range(len(sequence)-1)]) +
        end_AT*NN_dH_table['Terminal-AT-Penalty']
    )
    return dH

# Calculate ∆S, per SantaLucia 2004 Equation 1 and Table 1
def get_dS_SantaLucia2004(sequence):
    NN_dS_table = SantaLucia2004Table1['dS']
    end_AT = count_Terminal_ATs(sequence)
    dS = (
        NN_dS_table['Initiation'] +
        sum([NN_dS_table[sequence[i:i+2]] for i in range(len(sequence)-1)]) +
        end_AT*NN_dS_table['Terminal-AT-Penalty']
    )
    return dS

# Calculate Tm, per SantaLucia 2004 Equation 3
def get_Tm_SantaLucia2004(dH, dS, scaf_conc=SCAF_CONC, stap_conc=STAP_CONC):
    R_cal = R*1000  # Rescale R to [cal/mol•°K]
    x = 1           # Assume duplexes are self-complementary
    CT = stap_conc-0.5*scaf_conc  # Total molar strand concentration
    Tm_1M_NaCl =  dH * 1000/(dS + (R_cal * np.log(CT/x))) - 273.15
    return Tm_1M_NaCl

# Owczarzy 2008: Magnesium-corrected Tm

# Calculate magnesium-corrected Tm using Equation 16 from Owczarzy et al., Biochemistry, 2008
# Predicting Stability of DNA Duplexes in Solutions Containing Magnesium and Monovalent Cations
# http://pubs.acs.org/doi/abs/10.1021/bi702363u
def get_Tm_Mg_Owczarzy2008(sequence, seq_length, Tm_1M_NaCl, Mg_conc=Mg_CONC):
    # Empirical values from Table 2, units [1/°K]
    a =  3.92E-5
    b = -9.11E-6
    c =  6.26E-5
    d =  1.42E-5
    e = -4.82E-4
    f =  5.25E-4
    g =  8.31E-5

    # Fraction of G or C in sequence
    fGC = 1.0 * sum([nucleotide in ['G', 'C'] for nucleotide in sequence]) / seq_length
    ln_Mg = np.log(Mg_conc)

    Tm_Mg_inv = 1.0/Tm_1M_NaCl + a + b*ln_Mg + fGC*(c+d*ln_Mg) + 1.0/(2*(seq_length-1))*(e+f*ln_Mg+g*ln_Mg**2)
    Tm_Mg     = 1.0/Tm_Mg_inv

    return Tm_Mg


# Dunn 2015 salt-corrected dS entropy value 
# Combine constants into a single factor
EQN12_FACTOR = 0.368*np.log(0.5*TRIS_CONC + 3.3*np.sqrt(Mg_CONC))
def get_salt_corrected_dS_Dunn2015(seq_length):
    # Dunn et al. Equation 12
    # dSsalt = 0.368*(seq_length-1)*np.log(0.5*TRIS_CONC + 3.3*np.sqrt(Mg_CONC))
    dSsalt  = EQN12_FACTOR*(seq_length-1)
    return dSsalt

# Calculate the shortest distance in nucleotides between start and end index on a (typically circular) scaffold
def get_min_scaffold_distance(scaffold_length, start_index, end_index, is_scaffold_circular=True):
    '''Get min distance between two positions on a circular oligo'''
    forward_distance  = (end_index - start_index) % scaffold_length 
    reverse_distance = scaffold_length - forward_distance
    if is_scaffold_circular:
        return min(forward_distance, reverse_distance)
    else:
        return abs(end_index - start_index)

# Given a ssDNA length, calculate the end-to-end distance in [nm^2]
def end_to_end_distance(num_bases):
    '''
    Estimate end-to-end distance for DNA of length num_bases

    Model is from Rivetti et al. 1998, Eqn. 5
    https://www.sciencedirect.com/science/article/pii/S0022283698918307

    Constants are from Dunn et al. 2015, Eqn 12 
    https://www.nature.com/articles/nature14860 
    - Contour length L = 0.6 [nm]
    - Kuhn length λ_ss = 1.8 [nm]
    - Persistence length = λ_ss/2 = 0.9 [nm]
    '''
    L_base    = 0.6  # [nm]
    L_contour = num_bases*L_base
    L_p       = 0.9
    return 2*L_p*L_contour*(1-1.0*L_p/L_contour*(1-np.exp(-1.0*L_contour/L_p)))


# Calculate loop-closure entropy for given end-to-end distance
def distance_to_loop_dG(distance_square, temperature_kelvin=323.15):
    '''Calculate loop entropy for distance_square, which has units of nm^2'''
    # Equation from https://www.nature.com/articles/nature14860 (doi:10.1038/nature14860)
    # Divide by 1000 to convert cal/(mol.K) to kcal/(mol.K)
    effective_concentration = 1.0/6.02E23*(3.0/(2*np.pi*distance_square*1E-18))**(3/2.0)/1000
    dSloop = R*np.log(effective_concentration)
    dGloop = -temperature_kelvin*dSloop

    return (dGloop, dSloop)


# Given a start and end index
def position_to_loop_dG(scaffold_length, start_index, end_index, scaffold_circular=True, 
                        temperature_kelvin=323.15):
    # Get number of minimum number of bases between to location
    base_distance  = get_min_scaffold_distance(scaffold_length, start_index, end_index, scaffold_circular)

    # Determine end to end distance^2
    distance_square = end_to_end_distance(base_distance)

    # Get the free energies
    energies = distance_to_loop_dG(distance_square, temperature_kelvin)

    return energies


# Calculate melting temperature [°C] for given sequence
def sequence_to_Tm(sequence):
    seq_length = len(sequence)

    if seq_length == 0:
        return

    # Determine SantaLucia dH and dS
    dH = get_dH_SantaLucia2004(sequence)
    dS = get_dS_SantaLucia2004(sequence)

    # Equation from Santalucia paper (https://doi.org/10.1146/annurev.biophys.32.110601.141800)
    Tm_1M_NaCl = get_Tm_SantaLucia2004(dH, dS)

    # Perform Owczarzy Mg++ correction for oligomers
    if seq_length > 1:
        Tm_Mg = get_Tm_Mg_Owczarzy2008(sequence, seq_length, Tm_1M_NaCl)
        return Tm_Mg
    else:
        return Tm_1M_NaCl 


# Entropy factor due to concentration
def conc_to_dG(temperature_kelvin, stap_conc=STAP_CONC, scaf_conc=SCAF_CONC):
    dSconc = R*np.log(stap_conc-0.5*scaf_conc)
    dGconc = -temperature_kelvin*dSconc
    return (dGconc, dSconc)


# Calculate changes in free energy, enthalpy, and entropy for hybridization of given sequence
def sequence_to_dG_dH_dS(sequence, temperature_kelvin=323.15):
    seq_length = len(sequence)

    if seq_length == 0:
        return

    # Determine SantaLucia dH and dS
    dHtotal = get_dH_SantaLucia2004(sequence)
    dSbase  = get_dS_SantaLucia2004(sequence)

    # Determine Dunn2015 salt correction
    dSsalt  = get_salt_corrected_dS_Dunn2015(seq_length)

    # Get total entropy
    dStotal = (dSbase + dSsalt)/1000.0

    # Determine dG at the specified temperature 
    dGtotal = dHtotal - temperature_kelvin*dStotal

    return (dGtotal, dHtotal, dStotal)


# SEQUENCE GENERATORS

def generate_random_sequence(oligo_length):
    '''Generate random DNA sequence'''
    return ''.join(np.random.choice(['A', 'C', 'G', 'T'], int(oligo_length)))


def generate_nA(oligo_length):
    '''Generate AAAAA...AAA'''
    return 'A'*int(oligo_length)


def generate_nT(oligo_length):
    '''Generate TTTTT...TTT'''
    return 'T'*int(oligo_length)


def generate_nC(oligo_length):
    '''Generate CCCCC...CCC'''
    return 'C'*int(oligo_length)


def generate_nG(oligo_length):
    '''Generate GGGGG...GGG'''
    return 'G'*int(oligo_length)


# PARSING FUNCTIONS

def parse_break_rule(break_rule):
    '''Parse break rule'''
    return break_rule.split('.')


def parse_score_function(score_function):
    '''Parse break rule'''
    return score_function.split('.')


def parse_optim_function(function_input):
    '''Parse optimizatiotion function input'''
    # Get the function groups
    groups = function_input.split('.')
    # Get the functions and its parameters
    functions = [group.split(':') for group in groups]
    return functions


def parse_sequence_position(key_input):
    '''Parse sequence position'''
    return tuple([int(x) for x in key_input.split('.')])


# FILE COMPRESSION
def zip_directory(input_dir_path, output_zip_path):
    '''
    Compresses the contents of the specified directory into a zip file.

    :param input_dir_path: Path to the directory to be zipped.
    :param output_zip_path: Path to the output zip file.
    '''
    import os, zipfile
    with zipfile.ZipFile(output_zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(input_dir_path):
            for file in files:
                file_path = os.path.join(root, file)
                arcname = os.path.relpath(file_path, start=input_dir_path)
                zipf.write(file_path, arcname)

