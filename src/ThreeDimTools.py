from utilities import *
from Bio.PDB import PDBParser
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import gaussian_kde
from math import cbrt, sqrt
from scipy.ndimage import gaussian_filter
from itertools import combinations
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl
import matplotlib.patches as mpatches


molar_masses = {
    'H': 1.00794,    # Hydrogen
    'C': 12.0107,    # Carbon
    'N': 14.0067,    # Nitrogen
    'O': 15.9994,    # Oxygen
    'S': 32.065,     # Sulfur
    'P': 30.973762,  # Phosphorus
}


def get_struct(pdb_file: str):
    """ return a Bio.PDB.Structure object of the read pdb file"""
    pdb_parser = PDBParser()  # Instance for parsing the pdb file
    struct = pdb_parser.get_structure("structure_id", pdb_file)
    return struct


def filter_atoms(atoms, type=None):
    """ Returns a list of atoms of a specified type /!\ Could be optimized using the element (last column PDB)"""
    selection = []
    for atom in atoms:
        if type is None:
            selection.append(list(atom.get_coord()))
        else:
            if atom.get_name() == type:
                selection.append(list(atom.get_coord()))
    return selection


def scatter_atoms(atoms_coords):
    """ Display a scatter plot of all the atoms /!\ Laggy """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for atom in atoms_coords:
        ax.scatter(atom[0], atom[1], atom[2], c='b', marker='o')
    plt.show()


def heatmap_atoms(atoms_coords):
    """ Creates a heatmap of atom density : projection of the 3D structure on a 2D plan"""
    x, y, z = zip(*atoms_coords)
    kde = gaussian_kde([x, y])  # Gaussian kernel estimation (wtf tho?)
    xx, yy = np.meshgrid(np.linspace(min(x), max(x), 100), np.linspace(min(y), max(y), 100))  # Meshgrid for 2D Surface
    zz = np.reshape(kde([xx.ravel(), yy.ravel()]), xx.shape)
    plt.imshow(zz, cmap='hot', extent=[min(x), max(x), min(y), max(y)], origin='lower')  # Heatmap creation
    plt.colorbar(label='Density')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('2D heatmap projection of alpha carbons')
    plt.show()


def split_hydro(struct):
    """ Returns two sets of alpha carbons coordinates depending on the hydrophobicity of the amino-acid"""
    hydrophobic_aa = {'ALA', 'VAL', 'LEU', 'ILE', 'PHE', 'TRP', 'MET', 'PRO', 'GLY'}  # Exact ?
    hydrophobic_residues, polar_residues = [], []
    for model in struct:  # Parsing the structure
        for chain in model:
            for residue in chain:
                residue_name = residue.get_resname()
                if residue_name in hydrophobic_aa:  # Check for hydrophobicity
                    hydrophobic_residues.append(residue)
                else:
                    polar_residues.append(residue)
    return extract_CA_from_residue(hydrophobic_residues), extract_CA_from_residue(polar_residues)


def extract_CA_from_residue(residues):
    """ Extract the alpha carbon coordinates from a set of residues (aa) """
    CA = []
    for residue in residues:
        atoms = residue.get_atoms()
        for atom in atoms:
            if atom.get_name() == 'CA':
                CA.append(list(atom.get_coord()))
    return CA


def extract_barycenters_from_residue(residues):
    """ extract the barycenter coordinate of a residue (amino-acid)"""
    barys = []
    for residue in residues:
        atoms = residue.get_atoms()
        barys.append(get_barycenter(atoms))
    return barys


def get_barycenter(atoms):
    """ return a point in space which is the weighted barycenter of a set of atom"""
    total_weight = 0
    weighted_coord = np.zeros(3)
    for atom in atoms:
        weight = molar_masses[atom.get_name()[0]]
        weighted_coord += atom.get_coord() * weight
        total_weight += weight
    return weighted_coord/total_weight


def scatter_hydro(hydrophobic, polar):
    """ Double scatter plot for displaying dichotomed data, here hydrophobicity"""
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x1, y1, z1 = zip(*hydrophobic)  # Extract coordinate set
    x2, y2, z2 = zip(*polar)
    ax.scatter(x1, y1, z1, c='red', label='Hydrophobic')
    ax.scatter(x2, y2, z2, c='blue', label='Hydrophilic')
    plt.legend()
    plt.show()


def count_atoms_in_sphere(coordinates, radius):
    """ returns the number of atoms in the radius vicinity of each atoms of the coordinates set /!\ use np.array"""
    num_atoms = len(coordinates)
    atom_counts = []
    for i in range(num_atoms):  # First loop on atoms
        atom_count = 0
        for j in range(num_atoms):  # Iterating on each atom to see if close of the first one
            if i != j:
                distance = np.linalg.norm(coordinates[i] - coordinates[j])
                if distance <= radius:
                    atom_count += 1
        atom_counts.append(atom_count)
    return atom_counts


def create_atom_matrix_CA(struct):
    """ Return a 2D matrix of the alpha carbon distances"""
    CA = extract_CA_from_residue(struct.get_residues())  # List of coordinates
    size = len(CA)
    matrix = np.zeros((size, size))
    for i, icarbon in enumerate(CA):
        for j, jcarbon in enumerate(CA):
            matrix[i, j] = np.linalg.norm(np.array(icarbon)-np.array(jcarbon))
    return matrix


def create_atom_matrix_bary(struct):
    """ Return a 2D matrix of the resiudues barycenter distances"""
    bary = extract_barycenters_from_residue(struct.get_residues())  # List of coordinates
    size = len(bary)
    matrix = np.zeros((size, size))
    for i, ibary in enumerate(bary):
        for j, jbary in enumerate(bary):
            matrix[i, j] = np.linalg.norm(np.array(ibary)-np.array(jbary))
    return matrix


def display_atom_matrix(matrix1, matrix2, matrix3, threshold, mask=100):
    """ 2D overview of a protein structure. Mask aims to hide to far away residues """
    matrix1 = np.ma.masked_where(matrix1 > mask, matrix1)
    matrix2 = np.ma.masked_where(matrix2 > mask, matrix2)
    fig, axs = plt.subplots(2, 2, figsize=(12, 6))
    im1 = axs[0, 0].imshow(matrix1, cmap='inferno')
    axs[0, 0].set_title('Proximity matrix of alpha carbons in MRCP-20')
    axs[0, 0].set_xlabel('Residue N°')
    axs[0, 0].set_ylabel('Residue N°')
    cbar1 = fig.colorbar(im1, ax=axs[0, 0], shrink=0.55)
    cbar1.set_label('Distance in Angstroms')

    im2 = axs[0, 1].imshow(matrix2, cmap='inferno')
    axs[0, 1].set_title('Proximity matrix of residues barycenter in MRCP-20')
    axs[0, 1].set_xlabel('Residue N°')
    axs[0, 1].set_ylabel('Residue N°')
    cbar2 = fig.colorbar(im2, ax=axs[0, 1], shrink=0.55)
    cbar2.set_label('Distance in Angstroms')

    difference_matrix = np.abs(matrix1 - matrix2)
    im3 = axs[1, 0].imshow(difference_matrix, cmap='inferno')
    axs[1, 0].set_title('Absolute difference matrix')
    axs[1, 0].set_xlabel('Residue N°')
    axs[1, 0].set_ylabel('Residue N°')
    cbar3 = fig.colorbar(im3, ax=axs[1, 0], shrink=0.55)
    cbar3.set_label('Difference in Angstroms')

    axs[1, 1].imshow(matrix3, cmap='binary', interpolation='nearest')
    axs[1, 1].set_title(f'Interactions with dmax = {threshold}')
    axs[1, 1].set_xlabel('Residue N°')
    axs[1, 1].set_ylabel('Residue N°')

    plt.tight_layout()
    plt.show()


def can_interact(struct, chain, index1, index2, threshold):
    """ checks if the amino acids 1 and 2 are closer than the threshold in the provided struct"""
    res1 = struct[0][chain][index1]
    res2 = struct[0][chain][index2]
    coords = extract_barycenters_from_residue([res1, res2])
    distance = np.linalg.norm(np.array(coords[0]) - np.array(coords[1]))
    return distance < threshold


def create_atom_interaction_matrix(struct, threshold):
    """ Return a 2D matrix of the interactions"""
    bary = extract_barycenters_from_residue(struct.get_residues())  # List of coordinates
    size = len(bary)
    matrix = np.zeros((size, size))
    for i, ibary in enumerate(bary):
        for j, jbary in enumerate(bary):
            if j > i + int(threshold/2):
                matrix[i, j] = np.linalg.norm(np.array(ibary)-np.array(jbary)) < threshold
            else:
                matrix[i, j] = False
    return matrix


def get_thresholded_matrix_render(matrix, threshold, prot_id):
    """ let the values beyond the threshold blank on the plot"""
    masked_data = np.ma.masked_where(matrix > threshold, matrix)  # Mask values beyond threshold
    plt.imshow(masked_data, cmap='inferno', interpolation='nearest', vmin=0, vmax=threshold)
    colorbar = plt.colorbar()
    colorbar.set_label('Distance in Å')
    plt.tick_params(axis='y', pad=20)
    rect1 = mpatches.Patch(color='#000000', label='β-sheet')
    rect2 = mpatches.Patch(color='#EC4776', label='α-helix')
    plt.xlabel('Residue N°')
    plt.ylabel('Residue N°')
    plt.title(f'Proximity matrix of alpha carbons in {prot_id}', y=1.06)
    plt.legend(handles=[rect1, rect2], loc='upper left', bbox_to_anchor=(1.05, 1.17))
    fig = plt.gcf()
    fig.canvas.draw()
    image_np = np.array(fig.canvas.renderer.buffer_rgba())
    image_pil = Image.fromarray(image_np)
    return image_pil


def build_2D_render(prot):
    struct = get_struct(prot.pdb_path)  # retrieving the protein structure
    distances_render = get_thresholded_matrix_render(create_atom_matrix_CA(struct), 30, prot.prot_id)  # 2D dij mpl render

    data = open_saved_json(prot.features_path)  # opening features file
    feat = extract_features(data)
    struct = create_features_strand(feat, factor=369/(len(list(struct.get_residues()))))  # Creating the features img

    merged = merge_feat_and_graph(struct, distances_render)  # Render with pasted features
    merged.save(prot.render_2D_path)
    print('Render done')
    return prot.render_2D_path
