import scripts.bio_files_processor_scripts as bfp
import os


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta=None):
    """
    Creates fasta file with oneline sequences based on given fasta file with multiline sequences in the same directory.
    :param input_fasta: path to the multiline fasta file (str)
    :param output_fasta: name of oneline fasta file (str)
    """
    if output_fasta is None:
        output_fasta = f'oneline_{os.path.basename(os.path.realpath(input_fasta))}'
    if not output_fasta.endswith('.fasta'):
        output_fasta += '.fasta'

    location = os.path.dirname(os.path.realpath(input_fasta))
    seq = ''
    with open(input_fasta, mode='r') as fr, open(os.path.join(location, output_fasta), mode='w') as fw:
        for line in fr:
            if line.startswith('>'):
                if seq:
                    fw.write(seq + '\n')
                    seq = ''
                fw.write(line)
            else:
                seq += line.strip('\n')
        fw.write(seq)


def select_genes_from_gbk_to_fasta(*genes: str, input_gbk: str, n_before: int, n_after: int, output_fasta=None):
    """
    Creates fasta file with neighbor CDSs to given genes from GBK file and stores it
    in fasta_selected_from_gbk directory.
    :param genes: genes of interest that are used for neighbor CDSs search (str)
    :param input_gbk: path to GBK file (str)
    :param n_before: number of neighbor CDSs before gene of interest (int)
    :param n_after: number of neighbor CDSs after gene of interest (int)
    :param output_fasta: name of the output fasta file (str)
    """
    gene_list = genes
    cds_list = []
    gene_cds_dict = {}
    translation_dict = bfp.get_cds_translation_dict(input_gbk)

    with open(input_gbk, mode='r') as gbk:
        for line in gbk:
            if line.startswith('     CDS'):
                current_cds = line.replace(' ', '').strip('\n')
                cds_list.append(current_cds)
            if '/gene=' in line:
                current_gene = line.replace(' ', '').strip('\n').strip('/gene=').strip('"')
                gene_cds_dict[current_gene] = current_cds

        cds_of_interest = bfp.get_cds_of_interest(gene_list, gene_cds_dict, cds_list, n_before, n_after)

    os.makedirs('fasta_selected_from_gbk', exist_ok=True)

    if output_fasta is None:
        output_fasta = 'CDS_selected_from_gbk.fasta'
    if not output_fasta.endswith('.fasta'):
        output_fasta = output_fasta + '.fasta'

    bfp.get_fasta(output_fasta, cds_of_interest, translation_dict)


def change_fasta_start_pos(input_fasta: str, shift: int, output_fasta=None):
    """
    Change the starting position of a DNA sequence in a FASTA file.
    :param input_fasta: Path to the input FASTA file (str).
    :param shift: The number of positions to shift the sequence (int).
    :param output_fasta: Name of the output FASTA file (str).
    """
    with open(input_fasta, mode='r') as fa:
        lines = fa.readlines()
        seq_name = lines[0]
        seq = lines[1].strip()
        shifted_seq = f'{seq[shift:]}{seq[:shift]}\n'

    location = os.path.dirname(input_fasta)
    if output_fasta is None:
        output_fasta = f'Shifted{shift}.fasta'
    if not output_fasta.endswith('.fasta'):
        output_fasta = output_fasta + '.fasta'

    with open(os.path.join(location, output_fasta), mode='w') as sfa:
        sfa.write(seq_name)
        sfa.write(shifted_seq)





