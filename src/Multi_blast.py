import os

import numpy as np
from Bio.Blast import NCBIXML
import glob

from matplotlib import pyplot as plt

from Threshold_optimize import *
from inflexion_points import find_inflexion_points


def get_multi_blast_result(id):
    """Return a list of lists of sequences from XML files as a list of lists of dictionaries"""
    all_sequences = []
    directory = f'data/proteins/{id}/family/xml'
    for filename in os.listdir(directory):
        # print(filename)
        sequences = []
        with open(os.path.join(directory, filename), 'r') as f:
            blast_record = NCBIXML.read(f)
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    sequence_data = {
                        "title": alignment.title,
                        "length": alignment.length,
                        "e_value": hsp.expect,
                        "sequence": alignment.hsps[0].sbjct
                    }
                    sequences.append(sequence_data)
        all_sequences.append(sequences)
    return all_sequences


def multi_graph_Threshold(res, window_length=20):
    """ display the multi blast graph (mean, mean-sigma, mean+sigma) """
    x, y = zip(*[get_blast_xy(item) for item in res])
    y_array = np.array(y)
    x = np.array(x)
    y_mean = np.mean(y_array, axis=0)
    y_std = np.std(y_array, axis=0)

    dy, dy2, inflexion_up, inflexion_down, thresholds = find_inflexion_points(x[0], y_mean, window_length)
    # print(inflexion_up)

    # Tracer le graphique
    for i in range(len(y)):
        plt.semilogx(x[0, :], y[i], marker='.', alpha=0.03)
    plt.semilogx(x[0,:], y_mean, marker='.', color='blue', label='Mean')
    plt.semilogx(x[0, :], y_mean-y_std, marker='.', color='cyan', label='Mean - Std_deviation')
    plt.semilogx(x[0, :], y_mean+y_std, marker='.', color='magenta', label='Mean + Std_deviation')
    plt.plot(*zip(*inflexion_up), marker='o', markersize=8, color='red', linestyle='', label="Points d'inflexion vers le haut")
    # plt.twinx()
    # plt.semilogx(x[0, :], dy2, marker='.', color='green', label='dy2')
    plt.xlabel('Seuil e-value')
    plt.ylabel('Nombre de voisins moyen')
    plt.title('Nombre de voisins moyen trouvés en fonction du seuil e-value')
    plt.xscale('log')  # Échelle logarithmique pour l'axe x (optionnel)
    plt.grid(True)  # Activer la grille (optionnel)
    plt.legend()
    plt.show()  # Afficher le graphique
