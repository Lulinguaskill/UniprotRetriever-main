from UniprotFrontEnd import *
from DBTools import *
from utilities import *
from PDBTools import *
from ThreeDimTools import *
from BasicClasses import *
from blast_cmd import run_blast
from blastUtils import *
from Threshold_optimize import *
from inflexion_points import *
from config_manager import get_config

import re
import builtins


def split_family_fasta(id_p, input_file):
    """ divide the family fasto into several fasta.txt in the protein/family/fasta directory """

    fasta_directory = f"data/proteins/{id_p}/family/fasta"
    os.makedirs(fasta_directory, exist_ok=True)
   
    for record in SeqIO.parse(input_file, "fasta"):
        
        id = re.search(r"sp\|([A-Za-z0-9]+)\.", record.id).group(1)        

        output_file = f"data/proteins/{id_p}/family/fasta/{id}.fasta.txt"
        
        with builtins.open(output_file, 'w') as f:
            f.write(f'>{record.description}\n')
            f.write(f'{record.seq}\n')


def family_fasta_blast(id_p, database_path, threshold):  
    """ run BLAST on all fasta files from the id_p"""

    xml_directory = f"data/proteins/{id_p}/family/xml"
    os.makedirs(xml_directory, exist_ok=True)

    input_path = f"data/proteins/{id_p}/family/fasta"
    files = os.listdir(input_path)

    fasta_files = [file for file in files if file.endswith('.fasta.txt')]

    for file in fasta_files:
        file_path = f"data/proteins/{id_p}/family/fasta/{file}"

        parts = os.path.basename(file_path)
        id = parts.replace(".fasta.txt", "")
        output_file = f"data/proteins/{id_p}/family/xml/{id}.xml"

        threading.Thread(target=run_blast, args=(file_path, database_path, output_file, threshold)).start()


def family_blast(id, family):
    """ run all blast from the family of the id"""

    family_directory = f"data/proteins/{id}/family"
    os.makedirs(family_directory, exist_ok=True)

    output_fasta = f"data/proteins/{id}/{id}_family.fasta.txt"
    sequences_to_fasta(family.sequences, output_fasta)

    split_family_fasta(id, output_fasta)

    family_fasta_blast(id, get_config('config.yml')['db_path'], 10)


def get_sequences_from_xml_family(id_p):

    """ retourne une liste, chaque élement est un dictionnaire des résultats de blast d'un élément de la famille"""

    xml_directory = f"data/proteins/{id_p}/family/xml"

    files = os.listdir(xml_directory)

    fasta_files = [file for file in files if file.endswith('.xml')]

    total_sequences = []

    for file in fasta_files:

        with builtins.open(f'data/proteins/{id_p}/family/xml/{file}', 'r') as f:

            blast_record = NCBIXML.read(f)

            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    sequence_data = {
                        "title": alignment.title,
                        "length": alignment.length,
                        "e_value": hsp.expect,
                        "sequence": alignment.hsps[0].sbjct
                    }
                    total_sequences.append(sequence_data)

    return(total_sequences)
    



def get_percent_new(sequences, total_sequences):

    """renvoie le pourcentage de nouvelles séquences dans total_sequences et la liste des nouvlles séquences"""

    sequences_name = []
    total_sequences_names = []

    for seq in sequences:
        title = seq['title']
        id = re.search(r"sp\|([A-Za-z0-9]+)\.", title).group(1)
        sequences_name.append(id)

    for seq in total_sequences:
        title = seq['title']
        id = re.search(r"sp\|([A-Za-z0-9]+)\.", title).group(1)
        total_sequences_names.append(id)

    N_tot = len(total_sequences_names)

    counter = 0

    new_sequences = []

    for i, name in enumerate(total_sequences_names):
        if name not in sequences_name:
            counter += 1
            new_sequences.append(total_sequences[i])
            
    
    return counter/N_tot, new_sequences
 
        







