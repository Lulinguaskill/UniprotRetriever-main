""" File for dealing with pdb files and their 3D representation in pymol"""

import os
import requests
from Bio.PDB import *
from RequestTools import *


def exec_pymol_instruction(pdb_file, instructions):
    import pymolPy3
    """ execute the list of instructions on the pdb file """
    pm = pymolPy3.pymolPy3(0)  # Opening pymol instance
    pm(f'load {pdb_file}')  # Charging the pdb file into pymol
    for instruction in instructions:
        pm(instruction)
    pm.__del__()


def get_pymol_instructions(script_name):
    """ return a list of pymol instructions from the associated script"""
    commands = []
    with open(f'data/pdb/pymol scripts/{script_name}.txt', 'r') as f:
        lines = f.readlines()
    for line in lines:
        if not line.startswith(('#', ' ')):
            commands.append(line)
    return commands


def _download_pdb_file(file_url, output_path):
    """ download the pdb file and stores it in output_path"""
    try:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        response = requests.get(file_url)
        if response.status_code == 200:
            with open(output_path, 'wb') as f:
                f.write(response.content)
            print(f"PDB file downloaded successfully to: {output_path}")
        else:
            print(f"Error downloading PDB file. Status code: {response.status_code}")
    except Exception as e:
        print(f"An error occurred: {str(e)}")


def get_atoms(pdb_file):
    """ returns the atoms line of a pdb file"""
    with open(pdb_file, 'r') as file:
        lines = file.readlines()
    ATOM_lines = [line.strip() for line in lines if line.startswith('ATOM')]
    return ATOM_lines


def get_pdb_file(prot):
    """ Download and return the pdb file associated to the uniprot it"""
    if f'{prot.prot_id}.pdb' in os.listdir(prot.path):
        print(f'Already have pdb file for {prot.prot_id}')
        return os.path.join(prot.path, f'{prot.prot_id}.pdb')
    else:
        print(f'Downloading pdb file for {prot.prot_id}')
        prot.pdb_path = os.path.join(prot.path, f'{prot.prot_id}.pdb')
        url = get_pdb_url_from_uniprot_id(prot.prot_id)
        output_path = prot.pdb_path
        _download_pdb_file(url, output_path)
        return output_path
