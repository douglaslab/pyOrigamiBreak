import os
import sys

def read_sequence(file_path):
    with open(file_path, 'r') as file:
        return file.read().strip()

def write_permutation(sequence, index, folder):
    file_name = f'{os.path.basename(folder)}-{index:04d}.txt'
    file_path = os.path.join(folder, file_name)
    with open(file_path, 'w') as file:
        file.write(sequence)

def generate_circular_permutations(sequence, folder):
    if not os.path.exists(folder):
        os.makedirs(folder)

    n = len(sequence)
    for i in range(n):
        new_sequence = sequence[i:] + sequence[:i]
        write_permutation(new_sequence, i+1, folder)

def main():
    if len(sys.argv) != 2:
        print("Usage: python make_scaffold_permutations.py <filename>")
        sys.exit(1)

    input_file = sys.argv[1]
    sequence = read_sequence(input_file)
    folder_name = os.path.splitext(os.path.basename(input_file))[0]
    generate_circular_permutations(sequence, folder_name)

if __name__ == "__main__":
    main()
