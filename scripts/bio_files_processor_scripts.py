import os


def get_cds_translation_dict(input_gbk: str) -> dict:
    """
    Creates dict where CDSs are keys and gene names (if they are presented) and translation sequences
    are values from GBK file.
    :param input_gbk: path to GBK file (str)
    :return: dict where CDSs are keys and gene names (if they are presented) and translation sequences
    are values from GBK file (dict).
    """
    cds_dict = {}
    current_gene = ''
    current_seq = ''
    translation_start = False
    with open(input_gbk, mode='r') as gbk:
        for line in gbk:
            if line.startswith('     CDS'):
                translation_start = False
                if current_seq:
                    cds_dict[current_cds] = (current_gene, current_seq)
                    current_gene = ''
                    current_seq = ''
                current_cds = line.replace(' ', '').strip('\n')
            elif '/gene=' in line:
                current_gene = line.replace(' ', '').strip('\n').strip('/gene=').strip('"')
            elif '/translation' in line:
                current_seq = line.replace(' ', '').strip('\n').strip('/translation=')
                translation_start = True
            elif 'ORIGIN' in line:
                translation_start = False
            elif translation_start:
                current_seq += line.replace(' ', '').strip('\n')
    cds_dict[current_cds] = (current_gene, current_seq)
    return cds_dict


def get_fasta(output_fasta: str, cds_of_interest: list, translation_dict: dict):
    """
    Get FASTA file from given list of CDSs and their translations
    :param output_fasta: name for output FASTA file (str)
    :param cds_of_interest: list of cdc that is used for names
    :param translation_dict: dict where CDSs are keys and gene names (if they are presented) and translation sequences
    are values from GBK file (dict).
    """
    with open(os.path.join('fasta_selected_from_gbk', output_fasta), mode='w') as fasta:
        for cds in cds_of_interest:
            fasta.write(f'>{cds} gene: {translation_dict[cds][0]}\n')
            fasta.write(translation_dict[cds][1].replace('"', '') + '\n')

    print('FASTA file for neighbour CDSs of given genes is created ')