import os
import subprocess
import re

def delta_to_fasta(delta_filename, output_filename=None):
    if not os.path.exists(delta_filename):
        raise FileNotFoundError(f"Error: {delta_filename} does not exist.")

    if output_filename is None:
        assert len(delta_filename.rsplit('.')) == 2, 'Give the delta file a file extension so that a basename can be identified and used for the output filename, otherwise supply a filename of your choosing with the `output_filename` parameter.'
        delta_file_basename = delta_filename.rsplit('.')[0]
        fasta_filename = '{}.fasta'.format(delta_file_basename)
    else:
        fasta_filename = output_filename

    coords_file = subprocess.run(['show-coords', '-c', '-l', '-r', '-T', delta_filename],
                                  stdout=subprocess.PIPE).stdout.decode('utf-8')

    coords_lines = coords_file.split('\n')
    if len(coords_lines) < 5:
        raise ValueError("Error: Unexpected output from show-coords. Check if the delta file is valid.")

    coords_file_contents = coords_lines[4:-1]

    seq_names = None
    alignments = []
    fasta_strings = []

    for line in coords_file_contents:
        pct_identity = line.split('\t')[6]

        if seq_names != tuple(line.split('\t')[-2:]): 
            seq_names = tuple(line.split('\t')[-2:])
            first_seq_name, second_seq_name = seq_names

            aligns_file = subprocess.run(['show-aligns', delta_filename, first_seq_name, second_seq_name],
                                        stdout=subprocess.PIPE).stdout.decode('utf-8')
            alignments = alignments_from_aligns_file(aligns_file)

        if not alignments:
            raise ValueError(f"Error: No alignments found for {first_seq_name} vs {second_seq_name}.")

        alignment = alignments.pop(0)
        fasta_strings.append(sequences_lines_from_alignment(alignment, first_seq_name, second_seq_name, pct_identity))

    with open(fasta_filename, 'w') as fasta_file:
        fasta_file.write('\n\n'.join(fasta_strings))

    print(f'Finished: FASTA file written to {fasta_filename}')

def alignments_from_aligns_file(aligns_file):
    file_contents = '\n'.join(aligns_file.split('\n')[3:-3])
    alignments = file_contents.split('\n-- BEGIN alignment ')[1:]
    return alignments

def sequences_lines_from_alignment(alignment_string, first_seq_name='first_sequence', second_seq_name='second_sequence', pct_identity=None):
    split_alignment = alignment_string.split('\n\n', 2)  # Ensure at most 3 splits
    
    if len(split_alignment) != 3:
        raise ValueError(f"Error: Unexpected alignment format. Got {len(split_alignment)} sections instead of 3.")

    header, lines, footer = split_alignment

    coordinates_regex = r'\[ ([\+\-])1 ([0-9]+) \- ([0-9]+) \| ([\+\-])1 ([0-9]+) \- ([0-9]+) \]'
    match = re.findall(coordinates_regex, header)

    if not match:
        raise ValueError("Error: Could not extract coordinates from alignment header.")

    first_strand_dir, first_seq_start_coord, first_seq_end_coord, second_strand_dir, second_seq_start_coord, second_seq_end_coord = match[0]

    lines = lines.split('\n')
    lines = ["\n".join(lines[i:i+4]) for i in range(0, len(lines), 4)]

    first_sequence = ""
    second_sequence = ""

    for line in lines:
        try:
            _, first_sequence_part, second_sequence_part, _ = line.split('\n')
        except ValueError:
            raise ValueError(f"Error: Unexpected sequence line format in:\n{line}")

        remove_coord_numbers_regex = r'[0-9]+\s+(.*)'
        first_sequence_part = re.search(remove_coord_numbers_regex, first_sequence_part).group(1)
        second_sequence_part = re.search(remove_coord_numbers_regex, second_sequence_part).group(1)

        first_sequence += first_sequence_part
        second_sequence += second_sequence_part

    pct_identity_string = f' {pct_identity}% identity' if pct_identity else ''

    first_fasta_string = f'>{first_seq_name}:{first_seq_start_coord}-{first_seq_end_coord}({first_strand_dir}){pct_identity_string}\n{first_sequence}'
    second_fasta_string = f'>{second_seq_name}:{second_seq_start_coord}-{second_seq_end_coord}({second_strand_dir}){pct_identity_string}\n{second_sequence}'

    return f"{first_fasta_string}\n{second_fasta_string}"
  
if __name__ == "__main__":
    delta_filename = "l_dearborni_h1tg000019l_h2tg000034l.delta"
    output_filename = "l_dearborni_h1tg000019l_h2tg000034l.fasta"
    delta_to_fasta(delta_filename, output_filename)
