#!/usr/bin/env python3

import cadnano
from cadnano.document import Document
from cadnano.proxies.cnenum import GridType

class Validator(object):
    """
    Analyzes likely scaffold locations and staple lengths for a Cadnano part.
    """
    def __init__(self, part: object, scaffold_sequence:str, unbroken_max_len:int) -> None:
        super(Validator, self).__init__()
        self.part = part
        self.part_props = part_props = part.getModelProperties().copy()
        self.vh_order = part_props['virtual_helix_order']
        self.vh_props, self.vh_origins = part.helixPropertiesAndOrigins()
        self.vh_radius = part.radius()
        self.max_vhelix_length = max(self.vh_props['length'])
        self.square_lattice = True if part_props['grid_type'] == GridType.SQUARE else False
        self.oligos = part.oligos()
        self.scaffold_sequence = scaffold_sequence
        self.scaffold_length = len(scaffold_sequence)
        self.unbroken_max_len = unbroken_max_len
        self.scaf_lengths = []
        self.scaf_starts = []

        self.attributes = {
            'input': {
                'unbroken_max_len': self.unbroken_max_len,
                'lattice_type': 'square' if self.square_lattice else 'honeycomb'
            },
            'result': {
                'circular_oligo_count': None,
                'likely_scaffolds': [],
                'likely_scaffold_count': None,
                'long_staple_count': None
            }
        }
    # end def

    def analyze(self):
        # scaffold counts
        scaf_count = 0
        long_staple_count = 0
        circular_oligo_count = 0
        likely_scaffolds = []
        warnings = []

        for o in self.oligos:
            # find scaffolds
            o_len = o.length()

            # tally long oligos
            if o_len > self.unbroken_max_len:
                # is it a scaffold?
                if o.strand5p().strandSet().isScaffold():
                    scaf_count += 1
                    s = o.strand5p()
                    # length:vh.idx.polarity
                    direction = 'fwd' if s.isForward() else 'rev'
                    scaf = f"{o_len}:{s.idNum()}[{s.idx5Prime()}][{direction}]"
                    likely_scaffolds.append(scaf)
                    if o_len != self.scaffold_length:
                        warnings.append(f" WARNING: Input sequence length {self.scaffold_length} does not match design scaffold length {o_len}.")
                else:
                    long_staple_count += 1

            # tally circular
            if o.isCircular():
                circular_oligo_count += 1

        self.attributes['result']['circular_oligo_count'] = circular_oligo_count
        self.attributes['result']['likely_scaffold_count'] = scaf_count
        self.attributes['result']['likely_scaffolds'] = likely_scaffolds
        self.attributes['result']['long_staple_count'] = long_staple_count

        print("Analyzing input")
        print(f" Lattice type:          {self.attributes['input']['lattice_type']}")
        print(f" Circular oligo count:  {circular_oligo_count}")
        print(f" Scaffold count:        {scaf_count}")
        if scaf_count != 1:
            print(" WARNING: designs with more than one scaffold are unsupported.")
        likely_scaffolds_str = ', '.join(likely_scaffolds)
        plural = 's' if scaf_count > 1 else ''
        print(f" Likely scaffold{plural}:       {likely_scaffolds_str}")
        if len(warnings) > 0:
            for warning in warnings:
                print(warning)

    # import yaml
    # def writeYamlOutput(self, output_file):
    #     print('writing', output_file)
    #     with open(output_file, 'w') as outfile:
    #         yaml.dump(self.attributes, outfile)
# end class
