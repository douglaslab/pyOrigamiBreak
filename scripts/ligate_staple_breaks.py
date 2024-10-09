"""
# Instructions:

1. Make a virtual environment and install cadnano2.5:

	python -m venv ~/virtualenvs/ligate
	source ~/virtualenvs/ligate/bin/activate
	pip install git+https://github.com/douglaslab/cadnano2.5.git

2. Copy this script and your Cadnano design (mydesign.json) into the same folder, e.g. your Desktop

	cd ~/Desktop
	python ligate_staple_breaks.py <filename>

3. Several "Joining ..." print statements should appear. 
   The design will be saved with a "_ligated" suffix, e.g. mydesign_ligated.json

4. Examine the ligated version of the file with Cadnano, make any necessary changes, and then run autobreak.

"""
import sys
import cadnano
from cadnano.document import Document

def ligate_adjacent_staples(input_file, output_file):
	"""
	Opens a cadnano design file, identifies and ligates adjacent staple ends,
	and saves the modified design to a new file.

	Args:
		input_file (str): Path to the input cadnano JSON file.
		output_file (str): Path to the output cadnano JSON file.
	"""

	# Load the cadnano design
	app = cadnano.app()
	doc = app.document = Document()
	doc.readFile(input_file)
	part = doc.activePart()

	linear_breaks = []
	crossover_breaks = []
	ends5p, ends3p = [], []
	dict5p = {}

	# Iterate through all staple oligos, collecting 5' and 3' ends
	for oligo in part.oligos():
		s5p = oligo.strand5p()
		if s5p.strandSet().isScaffold(): continue
		id5, fwd5, idx5 = s5p.idNum(), s5p.isForward(), s5p.idx5Prime()
		ends5p.append((id5, fwd5, idx5))
		s3p = oligo.strand3p()
		id3, fwd3, idx3 = s3p.idNum(), s3p.isForward(), s3p.idx3Prime()
		ends3p.append((id3, fwd3, idx3))


	# Find linear breaks
	for id5, fwd5, idx5 in ends5p:
		dict5p[f'{id5}.{idx5}'] = fwd5  # Store table for later crossover lookups
		for id3, fwd3, idx3 in ends3p:
			if id5 == id3 and fwd5 == fwd3:
				if fwd5 == False and idx5 == idx3 - 1:
					linear_breaks.append((id5, fwd5, idx5, idx3))
				elif fwd5 == True and idx5 == idx3 + 1:
					linear_breaks.append((id5, fwd5, idx5, idx3))

	# Find crossover nicks
	for id3, fwd3, idx3 in ends3p:
		per_neighbor_hits, pairs = part.potentialCrossoverMap(id3, idx3)
		for neighbor_id, hits in per_neighbor_hits.items():
			fwd_axis_hits, rev_axis_hits = hits
			neighbor_hits = fwd_axis_hits if fwd3 else rev_axis_hits
			potential_crossover_idxs = [i[0] for i in neighbor_hits]
			if idx3 in potential_crossover_idxs:
				if f'{neighbor_id}.{idx3}' in dict5p:
					crossover_breaks.append((id3, fwd3, idx3, neighbor_id))

	# Merge linear breaks
	for vh_num, is_fwd, idx5, idx3 in linear_breaks:
		fwd_str = 'fwd' if is_fwd else 'rev'
		print(f"Joining linear break at {vh_num}[{idx5}] and {vh_num}[{idx3}]")
		fwd_strand_set, rev_strand_set = part.getStrandSets(vh_num)  # Get reverse strands (staples)
		stap_ss = fwd_strand_set if is_fwd else rev_strand_set
		strand = stap_ss.getStrand(idx5)
		strand.merge(idx5)

	# Merge broken crossovers
	for id3, fwd3, idx3, neighbor_id in crossover_breaks:
		fwd_strand_set3, rev_strand_set3 = part.getStrandSets(id3)
		ss3 = fwd_strand_set3 if fwd3 else rev_strand_set3
		strand3 = ss3.getStrand(idx3)
		fwd_strand_set5, rev_strand_set5 = part.getStrandSets(neighbor_id)
		ss5 = rev_strand_set5 if fwd3 else fwd_strand_set5
		strand5 = ss5.getStrand(idx3)
		print(f"Joining crossover break at {id3}[{idx3}] and {neighbor_id}[{idx3}]")
		part.createXover(strand3, idx3, strand5, idx3, use_undostack=False)

	# Save the modified design
	doc.writeToFile(output_file, legacy=True)

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print("Usage: python ligate_staple_nicks.py <filename>")
		sys.exit(1)
	input_file = sys.argv[1]
	output_file = input_file.replace(".json", "_ligated.json")  # Add suffix to output file
	ligate_adjacent_staples(input_file, output_file)
	print(f"Ligated design saved to {output_file}")
