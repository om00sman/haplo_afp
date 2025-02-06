#!/bin/bash
#SBATCH --job-name=hap_align
#SBATCH --mail-user=omoosman@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/hb/home/omoosman/owen/zoarcoidei/analysis/realignment/l_dearborni/mummer/align_%A_%a.out
#SBATCH --error=/hb/home/omoosman/owen/zoarcoidei/analysis/realignment/l_dearborni/mummer/align_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30GB
#SBATCH --time=5-00:00:00

#load required modules
module load miniconda3
conda activate mummer

nucmer --coords -p l_dearborni_h1tg000019l_h2tg000034l /hb/home/omoosman/owen/zoarcoidei/data/assemblies/l_dearborni/l_dearborni_h1tg000019l.fasta /hb/home/omoosman/owen/zoarcoidei/data/assemblies/l_dearborni/l_dearborni_h2tg000034l.fasta

import os
import subprocess
import re

def delta_to_fasta(delta_filename, output_filename=None):
    if output_filename is None:
        assert len(delta_filename.rsplit('.')) == 2, 'Give the delta file a file extension so that a basename can be identified and used for the output filename, otherwise supply a filename of your choosing with the `output_filename` parameter.'
        delta_file_basename = delta_filename.rsplit('.')[0]

        fasta_filename = '{}.fasta'.format(delta_file_basename)
    else:
        fasta_filename = output_filename

    coords_file = subprocess.run(['conda', 'run', '-n', 'mummer', 'show-coords', '-c', '-l', '-r', '-T', delta_filename],
                              stdout=subprocess.PIPE).stdout.decode('utf-8')

    coords_file_contents = coords_file.split('\n')[4:-1]

    seq_names = None
    alignments = []

    fasta_strings = []

    for line in coords_file_contents:
        pct_identity = line.split('\t')[6]

        if seq_names != tuple(line.split('\t')[-2:]): # this if clause should also work correctly in base case
#            assert len(alignments) == 0
            seq_names = tuple(line.split('\t')[-2:])
            first_seq_name, second_seq_name = seq_names
            
            aligns_file = subprocess.run(['conda', 'run', '-n', 'mummer','show-aligns', delta_filename, first_seq_name, second_seq_name],
                                        stdout=subprocess.PIPE).stdout.decode('utf-8')
            alignments = alignments_from_aligns_file(aligns_file)

        alignment = alignments.pop(0)
        fasta_strings.append(sequences_lines_from_alignment(alignment, first_seq_name, second_seq_name, pct_identity))

    with open(fasta_filename, 'w') as fasta_file:
        fasta_file.write('\n\n'.join(fasta_strings))
    print('Finished, aligned sequences with coordinates and percent identities of matches should have been written in FASTA format to {}'.format(fasta_filename))

def alignments_from_aligns_file(aligns_file):
    file_contents = '\n'.join(aligns_file.split('\n')[3:-3]) # get rid of cruft in first 3 and last 3 lines of output
    alignments = file_contents.split('\n-- BEGIN alignment ')[1:]
    return alignments

def sequences_lines_from_alignment(alignment_string, first_seq_name='first_sequence', second_seq_name='second_sequence', pct_identity=None):
    header, lines, footer = alignment_string.split('\n\n')

    coordinates_regex = r'\[ ([\+\-])1 ([0-9]+) \- ([0-9]+) \| ([\+\-])1 ([0-9]+) \- ([0-9]+) \]'
    first_strand_dir, first_seq_start_coord, first_seq_end_coord, second_strand_dir, second_seq_start_coord, second_seq_end_coord = re.findall(coordinates_regex, header)[0]

    lines = lines.split('\n')
    lines = ["\n".join(lines[i:i+4]) for i in range(0, len(lines), 4)]

    first_sequence = ""
    second_sequence = ""

    for line in lines:
        _, first_sequence_part, second_sequence_part, _ = line.split('\n')
        remove_coord_numbers_regex = r'[0-9]+\s+(.*)'
        first_sequence_part = re.search(remove_coord_numbers_regex, first_sequence_part).group(1)
        second_sequence_part = re.search(remove_coord_numbers_regex, second_sequence_part).group(1)

        first_sequence += first_sequence_part
        second_sequence += second_sequence_part

    if pct_identity is not None:
        pct_identity_string = ' {}% identity'.format(pct_identity)
    else:
        pct_identity_string = ''

    first_fasta_string = '>{}:{}-{}({}){}\n{}'.format(first_seq_name, first_seq_start_coord, first_seq_end_coord,
                                                      first_strand_dir, pct_identity_string, first_sequence)
    second_fasta_string = '>{}:{}-{}({}){}\n{}'.format(second_seq_name, second_seq_start_coord, second_seq_end_coord,
                                                       second_strand_dir, pct_identity_string, second_sequence)

    fasta_string = '\n'.join([first_fasta_string, second_fasta_string])
    return fasta_string

delta_filename = "l_dearborni_h1tg000019l_h2tg000034l.delta"
output_filename = "l_dearborni_h1tg000019l_h2tg000034l.fasta"
