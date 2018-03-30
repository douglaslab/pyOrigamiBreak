#!/usr/bin/env python
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


def write_oligos(oligos, stocks, json_files, fname, plate_header=''):
    '''
    Export oligos
    '''

    # Create white font for colored backgrounds
    white_font = openpyxl.styles.Font(color=openpyxl.styles.colors.WHITE)

    # Header for the plate sheets
    plate_column_header = ['Plate Name', 'Well Position', 'Sequence Name', 'Sequence',
                           'Stock Id', 'Oligo Start', 'Oligo End',
                           'Oligo Length', 'Oligo Color', 'Oligo Membership']

    wb = openpyxl.Workbook()

    # Current sheet
    ws1 = wb.active
    ws1.title = "Structures"

    # Get number of structures
    num_structures = len(json_files)

    # Update column width
    ws1.column_dimensions['A'].width = 20

    # Add header row
    ws1.append(('Structure', 'Stocks'))

    # 1. Write structures sheet
    for i in range(num_structures):
        stock_list = [json_files[i]]
        for j in range(len(stocks)):
            stock_id = stocks[j]
            if stock_id[2][i] == '1':
                stock_list.append(j+1)
        ws1.append(stock_list)

    # 2. Color cells for structures sheet
    maxRow = ws1.max_row
    maxCol = ws1.max_column
    for colNum in range(2, maxCol + 1, 1):
        for rowNum in range(2, maxRow + 1):
            current_cell = ws1.cell(row=rowNum, column=colNum)
            if current_cell.value:
                stock_id    = int(current_cell.value) - 1
                stock_color = openpyxl.styles.colors.Color(rgb=stocks[stock_id][1][1:])
                current_cell.fill = openpyxl.styles.fills.PatternFill(patternType='solid', fgColor=stock_color)
                current_cell.font = white_font

    # 3. Initialize Plate/Well numbers
    plate_id  = 1
    row_id    = 1
    col_id    = 0
    seq_id    = 0

    # Stock id/key
    prev_stock_key    = oligos[0].sortkey
    current_stock_key = None

    # Store sheet labels
    sheet_labels = [ws1.title]
    stock_id     = 1

    # Stock well positions
    stock_positions = {1: {'startwell': ('Plate%s-%d' % (plate_header, 1), 'A1'),
                           'endwell': ('Plate%s-%d' % (plate_header, 1), 'A1')}}

    # Aggregate sheet
    ws_ALL = wb.create_sheet(title='PlatesAll')

    # Add header row for the new sheet
    ws_ALL.append(plate_column_header)

    # Set column widths
    ws_ALL.column_dimensions['A'].width = 20
    ws_ALL.column_dimensions['D'].width = 60
    ws_ALL.column_dimensions['J'].width = 20

    for current_oligo in oligos:
        # Get current stock key
        current_stock_key = current_oligo.sortkey
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
        plate_label = 'Plate%s-%d' % (plate_header, plate_id)

        # Get well id
        well_id  = chr(64+row_id)+str(col_id)

        # Update start and endwells
        if current_stock_key == prev_stock_key:
            stock_positions[stock_id]['endwell'] = (plate_label, well_id)
        else:
            stock_positions[stock_id] = {}
            stock_positions[stock_id]['startwell'] = (plate_label, well_id)
            stock_positions[stock_id]['endwell']   = (plate_label, well_id)

        # Check if plate is in the list
        if plate_label not in sheet_labels:
            sheet_labels.append(plate_label)
            ws_current = wb.create_sheet(title=plate_label)

            # Add header row for the new sheet
            ws_current.append(plate_column_header)

            # Set column widths
            ws_current.column_dimensions['A'].width = 20
            ws_current.column_dimensions['D'].width = 60
            ws_current.column_dimensions['J'].width = 20

        # Row
        row = [plate_label, well_id, seq_id, current_oligo.sequence,
               stock_id, current_oligo.startkey, current_oligo.finishkey,
               current_oligo.length, current_oligo.color, current_oligo.bitseq]

        # Write the row
        ws_current.append(row)
        ws_ALL.append(row)

        # Update stock keys
        prev_stock_key = current_stock_key

    # 4. Color cells for structures sheet
    for current_sheet in wb:
        if not current_sheet.title[:5] == 'Plate':
            continue
        maxRow = current_sheet.max_row
        for rowNum in range(2, maxRow + 1):
            current_cell = current_sheet.cell(row=rowNum, column=5)
            if current_cell.value:
                stock_id    = int(current_cell.value) - 1
                stock_color = openpyxl.styles.colors.Color(rgb=stocks[stock_id][1][1:])
                current_cell.fill = openpyxl.styles.fills.PatternFill(patternType='solid', fgColor=stock_color)
                current_cell.font = white_font

    # 5. Write the Stocks sheet
    ws_stocks = wb.create_sheet(title='Stocks')
    # Append header
    ws_stocks.append(['Stock Id', 'Start Plate', 'Start Well', 'End Plate', 'End Well'])
    for stock_id in stock_positions:
        ws_stocks.append([stock_id, stock_positions[stock_id]['startwell'][0],
                          stock_positions[stock_id]['startwell'][1],
                          stock_positions[stock_id]['endwell'][0],
                          stock_positions[stock_id]['endwell'][1]])

    # Change the colors
    maxRow = ws_stocks.max_row
    for rowNum in range(2, maxRow + 1):
        current_cell = ws_stocks.cell(row=rowNum, column=1)
        if current_cell.value:
            stock_id    = int(current_cell.value) - 1
            stock_color = openpyxl.styles.colors.Color(rgb=stocks[stock_id][1][1:])
            current_cell.fill = openpyxl.styles.fills.PatternFill(patternType='solid', fgColor=stock_color)
            current_cell.font = white_font

    # Order the sheets
    sheet_order = np.arange(len(wb._sheets), dtype=int)
    sheet_order[2:] -= 1
    sheet_order[1] = len(wb._sheets) - 1
    wb._sheets = [wb._sheets[i] for i in sheet_order]

    wb.save(fname)


def read_sequence(sequence_file):
    '''
    Read sequence file
    '''

    # 1. Check the popular scaffolds list
    # 2. Check if the file exists

    # Initialize scaffold sequence
    scaffold_sequence = ''

    if sequence_file in utilities.SCAFFOLD_SEQUENCES:
        scaffold_sequence = utilities.SCAFFOLD_SEQUENCES[sequence_file]

    elif os.path.isfile(sequence_file):
        # Read sequence from file
        f = open(sequence_file)
        scaffold_sequence = ''.join([line.strip() for line in f.readlines()])
        f.close()

        # Convert to upper case
        scaffold_sequence = scaffold_sequence.upper()

    return scaffold_sequence


def read_oligos(cadnano_part, scaffold_sequence):
    cadnano_oligos = cadnano_part.oligos()
    cadnano_oligos = sorted(cadnano_oligos, key=lambda x: x.length(), reverse=True)

    # Apply sequence to scaffold
    cadnano_oligos[0].applySequence(scaffold_sequence)

    # Initialize the scaffolds and staples
    staples = {}

    # Iterate over oligos
    for oligo in cadnano_oligos[1:]:
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
            if strand.length() !=  len(strand.sequence()):
                new_oligo.sequence += strand.length()*'?'
            else:
                new_oligo.sequence += strand.sequence().replace(' ', '?')

        new_oligo.key       = '-'.join([str(new_oligo.vh5p),
                                        str(new_oligo.idx5p),
                                        str(new_oligo.vh3p),
                                        str(new_oligo.idx3p),
                                        new_oligo.sequence])

        new_oligo.startkey  = '%d[%d]' % (new_oligo.vh5p, new_oligo.idx5p)
        new_oligo.finishkey = '%d[%d]' % (new_oligo.vh3p, new_oligo.idx3p)

        staples[new_oligo.key] = new_oligo

    return staples


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

    # Get scaffold sequence
    scaffold_sequence = read_sequence(args.seq)

    # Check if sequence file exists
    if len(scaffold_sequence) == 0:
        sys.exit('Scaffold sequence is not available!')

    # Apply the offset
    scaffold_sequence = scaffold_sequence[args.offset:]+scaffold_sequence[:args.offset]

    # Define output file
    head, tail       = os.path.split(os.path.abspath(json_files[0]))
    xlsx_output = head+'/oligos_96well.xlsx'

    # Initialize cadnano
    app = cadnano.app()
    doc = app.document = Document()

    # All staples dictionary
    all_staples_dict = {}

    # 1. Iterate over all the json files
    for i in range(len(json_files)):
        # Get json file
        json_input = json_files[i]

        # Read cadnano input file
        doc.readFile(json_input)

        # Get part
        part = doc.activePart()

        # Get staples
        staples = read_oligos(part, scaffold_sequence)

        # Add bit key to staples and add to all staples list
        for key in staples:
            staples[key].bitlist = [0]*len(json_files)

            if key not in all_staples_dict:
                all_staples_dict[key] = staples[key]

        # Assign the bit for the staples
        for key in all_staples_dict:
            if key in staples:
                all_staples_dict[key].bitlist[i] = 1

    # 2. Count colors
    color_counts = {}
    for key in all_staples_dict:
        color = all_staples_dict[key].color
        if color not in color_counts:
            color_counts[color] = 0
        color_counts[color] += 1

    # 3. Assign the color counters for the staples and make bitseq
    all_staples_list = []
    for key in all_staples_dict:
        all_staples_dict[key].colorctr = color_counts[all_staples_dict[key].color]
        all_staples_dict[key].bitseq   = ''.join([str(x) for x in all_staples_dict[key].bitlist])
        all_staples_dict[key].sortkey  = (all_staples_dict[key].colorctr,
                                          all_staples_dict[key].color,
                                          all_staples_dict[key].bitseq)
        # Add staple to list
        all_staples_list.append(all_staples_dict[key])

    # 4. Sort staples based on length
    all_staples_list.sort(key=lambda x: x.length)

    # 5. Sort staples based on the sort key
    all_staples_list.sort(key=lambda x: x.sortkey, reverse=True)

    # 6. Get the stock names
    stocks  = []
    for oligo in all_staples_list:
        if oligo.sortkey not in stocks:
            stocks.append(oligo.sortkey)

    # 7. Export oligos
    write_oligos(all_staples_list, stocks, json_files, xlsx_output, args.header)


if __name__ == "__main__":
    main()
