#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2018-03-03 15:44:55
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

from cadnano.document import Document

import os
import sys
import cadnano
import argparse
import glob
import utilities
import openpyxl
import numpy as np


# GLOBALS
WHITE_FONT = openpyxl.styles.Font(color=openpyxl.styles.colors.WHITE)


class Project:
    '''
    Origami Project class that keeps multiple structures
    '''
    def __init__(self):
        self.json_files     = []
        self.structures     = []
        self.plates         = []
        self.date           = ''
        self.author         = ''
        self.stocks         = {}
        self.oligos_list    = []
        self.oligos_dict    = {}

        self.stock_counts   = {}    # Number of stocks
        self.stock_keys     = []
        self.total_bases    = 0     # Total number of bases
        self.num_structures = 0
        self.wb             = None  # Workbook object

    def assign_color_counters(self):
        '''
        Assign color counters
        '''

        for key in self.oligos_dict:
            self.oligos_dict[key].colorctr = self.color_counts[self.oligos_dict[key].color]

    def assign_sortkeys(self):
        '''
        Assign sorkeys
        '''
        for key in self.oligos_dict:
            self.oligos_dict[key].sortkey  = (self.oligos_dict[key].colorctr,
                                              self.oligos_dict[key].color,
                                              self.oligos_dict[key].bitseq)

    def assign_bitseqs(self):
        '''
        Assign bitseqs
        '''
        for key in self.oligos_dict:
            self.oligos_dict[key].bitseq   = ''.join([str(x) for x in self.oligos_dict[key].bitlist])

    def prepare_oligos_list(self):
        '''
        Prepare oligos list
        '''
        self.oligos_list = []
        for key in self.oligos_dict:

            # Add staple to list
            self.oligos_list.append(self.oligos_dict[key])

        # 4. Sort staples based on length
        self.oligos_list.sort(key=lambda x: x.length)

        # 5. Sort staples based on the sort key
        self.oligos_list.sort(key=lambda x: x.sortkey, reverse=True)

    def count_colors(self):
        '''
        Count colors
        '''
        self.color_counts = {}
        for key in self.oligos_dict:
            color = self.oligos_dict[key].color
            if color not in self.color_counts:
                self.color_counts[color] = 0
            self.color_counts[color] += 1

    def add_structure(self, new_structure):
        '''
        Add new structure
        '''
        self.structures.append(new_structure)

        # Reset bits for the new structure
        new_structure.reset_bits(self.num_structures)

        # Add staples to main list
        for key in new_structure.oligos_dict:

            if key not in self.oligos_dict:
                self.oligos_dict[key] = new_structure.oligos_dict[key]

        # Assign project
        new_structure.project = self

    def assign_bit_id(self, new_structure):
        '''
        Add new staples to staples list
        '''

        # Assign the bit for the staples
        for key in self.oligos_dict:
            if key in new_structure.oligos_dict:
                self.oligos_dict[key].bitlist[new_structure.structure_id] = 1

    def read_sequence(self, sequence_file):
        '''
        Read sequence file
        '''

        # 1. Check the popular scaffolds list
        # 2. Check if the file exists

        # Initialize scaffold sequence
        self.scaffold_sequence = ''

        if sequence_file in utilities.SCAFFOLD_SEQUENCES:
            self.scaffold_sequence = utilities.SCAFFOLD_SEQUENCES[sequence_file]

        elif os.path.isfile(sequence_file):
            # Read sequence from file
            f = open(sequence_file)
            self.scaffold_sequence = ''.join([line.strip() for line in f.readlines()])
            f.close()

            # Convert to upper case
            self.scaffold_sequence = self.scaffold_sequence.upper()

        return self.scaffold_sequence

    def parse_input_files(self, input_fnames):
        '''
        Parse input files
        '''

        self.json_files = []
        for input_fname in input_fnames:
            self.json_files += glob.glob(input_fname)

        # Filter out any non-json file
        self.json_files = list(filter(lambda x: os.path.splitext(x)[1] == '.json', self.json_files))

        # Sort json files
        self.json_files.sort()

        return self.json_files

    def _create_workbook_96well(self):
        '''
        Create workbook for 96well plate output format
        '''
        self.wb = openpyxl.Workbook()

    def calculate_bases(self):
        '''
        Calculate total number of bases
        '''
        self.total_bases = 0

        for current_oligo in self.oligos_list:
            # Update total bases
            self.total_bases += current_oligo.length

    def prepare_stocks(self):
        '''
        Prepare stocks
        '''
        self.stocks       = {}
        self.stock_counts = {}
        self.stock_keys   = []

        # Initialize stock id
        stock_id = 0

        for current_oligo in self.oligos_list:

            if current_oligo.sortkey not in self.stock_keys:
                # Create new stock
                new_stock = Stock()
                new_stock.key         = current_oligo.sortkey
                new_stock.count       = 0
                new_stock.stock_id    = stock_id
                new_stock.oligos_list = [current_oligo]
                new_stock.color       = current_oligo.color

                self.stocks[current_oligo.sortkey]       = new_stock
                self.stock_counts[current_oligo.sortkey] = 0
                self.stock_keys.append(current_oligo.sortkey)

                # Update stock id
                stock_id += 1

            # Update oligo counts for each stock
            self.stock_counts[current_oligo.sortkey] += 1
            self.stocks[current_oligo.sortkey].count += 1
            self.stocks[current_oligo.sortkey].oligos_list.append(current_oligo)

    def add_json_files(self, json_files):
        '''
        Add json files
        '''
        self.json_files = json_files

        # Get number of structures
        self.num_structures = len(self.json_files)

    def _prepare_structure_sheet_96well(self):
        '''
        Prepare structures
        '''
        for i in range(self.num_structures):
            # Get structure
            current_structure = self.structures[i]

            stock_row = [self.json_files[i]]

            # Starting H2O volume
            water = 300
            for j in range(len(self.stock_keys)):
                current_stock_key = self.stock_keys[j]

                # Get membership key
                membership_key = current_stock_key[2]

                # Check stock membership
                if membership_key[i] == '1':
                    stock_row.append('%d (%d)' % (j+1, self.stocks[current_stock_key].count))
                    water -= self.stocks[current_stock_key].count

                    # Add stock to structure
                    if not current_stock_key in current_structure.stocks:
                        current_structure.stocks[current_stock_key] = self.stocks[current_stock_key]
            # Add water
            stock_row.append('H2O (%d)' % (water))
            self.ws_str.append(stock_row)

    def _color_structure_sheet_96well(self):
        '''
        Color structrue sheet for 96 well plate
        '''
        # Get maximum row and column for structrue sheet
        maxRow = self.ws_str.max_row
        maxCol = self.ws_str.max_column

        for colNum in range(2, maxCol + 1, 1):
            for rowNum in range(2, maxRow + 1):
                # Get cell
                current_cell = self.ws_str.cell(row=rowNum, column=colNum)
                if current_cell.value and current_cell.value[0] != 'H':

                    # Get stock id and stock color
                    stock_id    = int(current_cell.value.split('(')[0]) - 1
                    stock_color = openpyxl.styles.colors.Color(rgb=self.stock_keys[stock_id][1][1:])

                    # Set color
                    current_cell.fill = openpyxl.styles.fills.PatternFill(patternType='solid', fgColor=stock_color)

                    # Set font
                    current_cell.font = WHITE_FONT

    def _write_structure_sheet_96well(self):
        '''
        Write structure sheet for 96 well plate format output
        '''
        self.ws_str = self.wb.active
        self.ws_str.title = "Structures"

        # Update column width
        self.ws_str.column_dimensions['A'].width = 20

        # Add header row
        self.ws_str.append(('Structure', 'Stocks (ul)'))

        # Calculate bases and stocks
        self.calculate_bases()
        self.prepare_stocks()

        # Prepare structures sheet
        self._prepare_structure_sheet_96well()

        # Color cells for structures sheet
        self._color_structure_sheet_96well()

        # Append total base count
        self.ws_str.append(['Total bases', self.total_bases])

    def _write_plate_sheets_96well(self, plate_header=''):

        # Header for the plate sheets
        self.plate_column_header = ['Plate Name', 'Well Position', 'Sequence Name', 'Sequence',
                                    'Stock Id', 'Oligo Start', 'Oligo End',
                                    'Oligo Length', 'Oligo Color', 'Oligo Membership']

        # Initialize Plate/Well numbers
        plate_id  = 1
        row_id    = 1
        col_id    = 0
        seq_id    = 0
        stock_id  = 1

        # Stock id/key
        prev_stock_key    = self.oligos_list[0].sortkey
        current_stock_key = None

        # Store sheet labels
        self.sheet_labels = [self.ws_str.title]

        # Save plate header
        self.plate_header = plate_header

        # Get plate label and assign the first startwell
        plate_label = 'Plate%s-%d' % (self.plate_header, plate_id)
        self.stocks[prev_stock_key].startwell = (plate_label, 'A1')

        # Aggregate sheet
        self.ws_ALL = self.wb.create_sheet(title='PlatesAll')

        # Add header row for the new sheet
        self.ws_ALL.append(self.plate_column_header)

        # Set column widths
        self.ws_ALL.column_dimensions['A'].width = 20
        self.ws_ALL.column_dimensions['D'].width = 60
        self.ws_ALL.column_dimensions['J'].width = 20

        # Initialize the current plate
        current_plate = None

        for current_oligo in self.oligos_list:
            # Get current stock key
            current_stock_key = current_oligo.sortkey

            # Update sequence id
            seq_id += 1

            if current_stock_key == prev_stock_key:
                # Update column id
                col_id += 1
            else:
                stock_id += 1
                col_id    = 1
                row_id   += 1

            # Update columns if column id is over 12
            if col_id > 12:
                col_id = 1
                row_id += 1

            # Update rows if row id is over 8
            if row_id > 8:
                row_id = 1
                plate_id += 1

            # Get plate id
            plate_label = 'Plate%s-%d' % (self.plate_header, plate_id)

            # Get well id
            well_id  = chr(64+row_id)+str(col_id)

            # Update start and endwells
            if current_stock_key == prev_stock_key:
                self.stocks[current_stock_key].endwell   = (plate_label, well_id)
            else:
                self.stocks[current_stock_key].startwell = (plate_label, well_id)
                self.stocks[current_stock_key].endwell   = (plate_label, well_id)

            # Check if plate is in the list
            if plate_label not in self.sheet_labels:

                self.sheet_labels.append(plate_label)
                ws_current = self.wb.create_sheet(title=plate_label)

                # Add header row for the new sheet
                ws_current.append(self.plate_column_header)

                # Set column widths
                ws_current.column_dimensions['A'].width = 20
                ws_current.column_dimensions['D'].width = 60
                ws_current.column_dimensions['J'].width = 20

                # Create new 96well plate
                current_plate = Plate()
                current_plate.oligos = {}
                current_plate.worksheet   = ws_current
                current_plate.plate_label = plate_label

                # Add current plate to plates list
                self.plates.append(current_plate)

            # Row
            current_plate.oligos[well_id] = current_oligo

            # Update current oligo entries
            current_oligo.plate96_well_id = well_id
            current_oligo.plate96_seq_id  = seq_id
            current_oligo.plate96         = current_plate
            current_oligo.stock           = self.stocks[current_stock_key]

            row = [plate_label, well_id, seq_id, current_oligo.sequence,
                   stock_id, current_oligo.startkey, current_oligo.finishkey,
                   current_oligo.length, current_oligo.color, current_oligo.bitseq]

            # Write the row
            ws_current.append(row)
            self.ws_ALL.append(row)

            # Update stock keys
            prev_stock_key = current_stock_key

    def _color_structures_sheet_96well(self):
        # Color cells for structures sheet
        for current_sheet in self.wb:
            if not current_sheet.title[:5] == 'Plate':
                continue
            maxRow = current_sheet.max_row
            for rowNum in range(2, maxRow + 1):
                current_cell = current_sheet.cell(row=rowNum, column=5)
                if current_cell.value:
                    stock_id    = int(current_cell.value) - 1
                    stock_color = openpyxl.styles.colors.Color(rgb=self.stock_keys[stock_id][1][1:])
                    current_cell.fill = openpyxl.styles.fills.PatternFill(patternType='solid', fgColor=stock_color)
                    current_cell.font = WHITE_FONT

    def _write_stocks_sheet_96well(self):
        # Write stocks sheet
        self.ws_stocks = self.wb.create_sheet(title='Stocks')
        # Append header
        self.ws_stocks.append(['Stock Id', 'Start Plate', 'Start Well', 'End Plate', 'End Well', '#Oligos'])

        for current_stock_key in self.stock_keys:
            self.ws_stocks.append([self.stocks[current_stock_key].stock_id+1,
                                   self.stocks[current_stock_key].startwell[0],
                                   self.stocks[current_stock_key].startwell[1],
                                   self.stocks[current_stock_key].endwell[0],
                                   self.stocks[current_stock_key].endwell[1],
                                   self.stocks[current_stock_key].count])

        # Change the colors
        maxRow = self.ws_stocks.max_row
        for rowNum in range(2, maxRow + 1):
            current_cell = self.ws_stocks.cell(row=rowNum, column=1)
            if current_cell.value:
                stock_id    = int(current_cell.value) - 1
                stock_color = openpyxl.styles.colors.Color(rgb=self.stock_keys[stock_id][1][1:])
                current_cell.fill = openpyxl.styles.fills.PatternFill(patternType='solid', fgColor=stock_color)
                current_cell.font = WHITE_FONT

    def _order_sheets_96well(self):
        # Order the sheets

        sheet_order = np.arange(len(self.wb._sheets), dtype=int)
        sheet_order[2:] -= 1
        sheet_order[1] = len(self.wb._sheets) - 1
        self.wb._sheets = [self.wb._sheets[i] for i in sheet_order]

    def save_sheets_96well(self, out_fname):
        # Save excel file (workbook)
        self.wb.save(out_fname)

    def write_oligos_96well(self, out_fname, plate_header=''):
        '''
        Export oligos
        '''
        self._create_workbook_96well()
        self._write_structure_sheet_96well()
        self._write_plate_sheets_96well(plate_header)
        self._color_structures_sheet_96well()
        self._write_stocks_sheet_96well()
        self._order_sheets_96well()
        self.save_sheets_96well(out_fname)


class Structure:
    '''
    Structure class that keeps the oligos to order and the Cadnano input file
    '''
    def __init__(self):
        self.structure_name = ''
        self.structure_id   = ''
        self.json_file      = ''
        self.oligos_list    = ''
        self.oligos_dict    = ''
        self.project        = None
        self.stocks         = {}

    def reset_bits(self, num_structures):
        '''
        Reset bits
        '''
        for key in self.oligos_dict:
            self.oligos_dict[key].bitlist = [0]*num_structures

    def read_oligos(self, cadnano_part, scaffold_sequence, add_T=False):
        self.cadnano_oligos = cadnano_part.oligos()
        self.cadnano_oligos = sorted(self.cadnano_oligos, key=lambda x: x.length(), reverse=True)

        # Apply sequence to scaffold
        self.cadnano_oligos[0].applySequence(scaffold_sequence)

        # Initialize the scaffolds and staples
        self.oligos_dict = {}

        # Empty character
        empty_ch = '?'
        if add_T:
            empty_ch = 'T'

        # Iterate over oligos
        for oligo in self.cadnano_oligos[1:]:
            strand5p  = oligo.strand5p()
            strand3p  = oligo.strand3p()

            # Make new oligo
            new_oligo           = Oligo()
            new_oligo.length    = oligo.length()
            new_oligo.vh5p      = strand5p.idNum()
            new_oligo.vh3p      = strand3p.idNum()
            new_oligo.idx5p     = strand5p.idx5Prime()
            new_oligo.idx3p     = strand3p.idx3Prime()
            new_oligo.color     = oligo.getColor()
            new_oligo.sequence  = ''

            # Get the strand generator
            generator           = oligo.strand5p().generator3pStrand()
            for strand in generator:
                if strand.totalLength() != len(strand.sequence()):
                    new_oligo.sequence += strand.length()*empty_ch
                else:
                    new_oligo.sequence += strand.sequence().replace(' ', empty_ch)

            new_oligo.key       = '-'.join([str(new_oligo.vh5p),
                                            str(new_oligo.idx5p),
                                            str(new_oligo.vh3p),
                                            str(new_oligo.idx3p),
                                            new_oligo.sequence])

            new_oligo.startkey  = '%d[%d]' % (new_oligo.vh5p, new_oligo.idx5p)
            new_oligo.finishkey = '%d[%d]' % (new_oligo.vh3p, new_oligo.idx3p)

            self.oligos_dict[new_oligo.key] = new_oligo

        return self.oligos_dict


class Plate:
    '''
    Plate class
    '''
    def __init__(self):
        self.plate_label = ''
        self.plate_id    = ''
        self.plate_size  = 96
        self.oligos      = {}
        self.worksheet   = None

        self.row_length  = 8
        self.col_length  = 12

    def set_dimensions(self, plate_size=384):
        # Assign plate dimensions
        self.plate_size = plate_size

        if self.plate_size == 384:
            self.row_length = 16
            self.col_length = 24
        elif self.plate_size == 96:
            self.row_length = 8
            self.col_length = 12

    def advance_row_order(self, row_id, column_id):
        '''
        Advance row order
        '''
        # Get new row id
        new_column_id = (column_id+1) % self.col_length
        new_row_id    = row_id

        # Get new column id
        if column_id+1 > self.col_length:
            new_row_id = row_id + 1

        # Check new row id
        if new_row_id > self.row_length:
            return False
        else:
            return (new_row_id, new_column_id)

    def advance_col_order(self, row_id, column_id):
        '''
        Advance col order
        '''

    def get_well_id(self, row_id, column_id):
        '''
        Get well id
        '''


class Stock:
    '''
    Origami stock class
    '''
    def __init__(self):
        self.key         = ''
        self.color       = ''
        self.count       = 0
        self.stock_id    = ''
        self.stock_name  = ''
        self.oligos_list = []
        self.oligos_dict = {}

        self.startwell   = ()
        self.endwell     = ()


class Oligo:
    def __init__(self):
        self.key      = None
        self.color    = None
        self.colorctr = None
        self.sequence = None
        self.vh5p     = None
        self.idx5p    = None
        self.vh3p     = None
        self.idx3p    = None
        self.bitlist  = None
        self.bitseq   = None
        self.sortkey  = None

        # Plate and stock parameters
        self.plate96_well_id  = None
        self.plate96_seq_id   = None
        self.plate96          = None

        self.plate384_well_id = None
        self.plate384_seq_id  = None
        self.plate384         = None

        self.stock            = None


def parse_input_files(input_fnames):
    '''
    Parse input files
    '''

    json_files = []
    for input_fname in input_fnames:
        json_files += glob.glob(input_fname)

    # Filter out any non-json file
    json_files = list(filter(lambda x: os.path.splitext(x)[1] == '.json', json_files))

    # Sort json files
    json_files.sort()

    return json_files


def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i",   "--input",    type=str, nargs='+', required=True,
                        help="Cadnano json file(s)")

    parser.add_argument("-seq",   "--seq", type=str,
                        help="Scaffold sequence file")

    parser.add_argument("-header", "--header", type=str, default='',
                        help="Plate header")

    parser.add_argument("-offset", "--offset", type=int, default=0,
                        help="Sequence offset")

    parser.add_argument("-addT", "--addT", action='store_true',
                        help="Replace ? with T")

    args = parser.parse_args()

    # Check if the required arguments are passed to the code
    if args.input is None:
        parser.print_help()
        sys.exit('Input file does not exist!')

    # Get json files
    json_files = parse_input_files(args.input)

    # Check if the json file list is empty
    if len(json_files) == 0:
        sys.exit('Input file does not exist!')

    # Create a project
    new_project = Project()

    # Get scaffold sequence
    scaffold_sequence      = new_project.read_sequence(args.seq)

    # Add json files
    new_project.add_json_files(json_files)

    # Check if sequence file exists
    if len(scaffold_sequence) == 0:
        sys.exit('Scaffold sequence is not available!')

    # Apply the offset
    new_project.scaffold_sequence = (new_project.scaffold_sequence[args.offset:] +
                                     new_project.scaffold_sequence[:args.offset])

    # Define output file
    head, tail         = os.path.split(os.path.abspath(json_files[0]))
    xlsx_output_96well  = head+'/oligos_96well.xlsx'
    xlsx_output_384well = head+'/oligos_384well.xlsx'
    echo_output_384well = head+'/echo_384well.cvs'

    # Initialize cadnano
    app = cadnano.app()
    doc = app.document = Document()

    # 1. Iterate over all the json files
    for i in range(len(json_files)):

        # Get json file
        json_input = json_files[i]

        # Define json output
        head, ext  = os.path.splitext(os.path.abspath(json_input))

        # Json output
        json_output = head+'_export.json'

        # Read cadnano input file
        doc.readFile(json_input)

        # Get part
        part = doc.activePart()

        # Get staples
        # Create structure object
        new_structure = Structure()
        new_structure.oligos_dict  = new_structure.read_oligos(part, new_project.scaffold_sequence, args.addT)
        new_structure.structure_id = i

        # Write json output
        doc.writeToFile(json_output, legacy=False)

        # Add structure to list
        new_project.add_structure(new_structure)

        # Assign the bit for the staples
        new_project.assign_bit_id(new_structure)

    # 2. Count colors
    new_project.count_colors()

    # 3. Assign the color counters for the staples and make bitseq
    new_project.assign_color_counters()
    new_project.assign_bitseqs()
    new_project.assign_sortkeys()
    new_project.prepare_oligos_list()

    # 7. Export oligos for 96 well plate format
    new_project.write_oligos_96well(xlsx_output_96well, args.header)


if __name__ == "__main__":
    main()
