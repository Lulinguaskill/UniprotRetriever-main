import numpy as np
from scipy.signal import find_peaks
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
from BasicClasses import *
from Threshold_optimize import *
from blastUtils import filter_blast_result



def find_inflexion_points(x, y, n=20):
    # calcul des points d'inflection (vers le haut et vers le bas) d'une courbe y=f(x)
    # en utilisant les pics (positifs ou négatifs) de la dérivé seconde
    y_smoothed = savgol_filter(y, window_length=n, polyorder=2)
    dy = np.gradient(y_smoothed)
    dy2 = np.gradient(dy)

    peaks, _ = find_peaks(dy2, height=max(dy2)/4, prominence=max(dy2)/3)  # ,  distance=2, width=2,
    # threshold=1, rel_height=0.5, wlen=1)

    inflexion_up = [(x[p], y[p]) for p in peaks if dy2[p] > 0]
    inflexion_down = [(x[p], y[p]) for p in peaks if dy2[p] < 0]

    thresholds = [t[0] for t in inflexion_up]

    return dy, dy2, inflexion_up, inflexion_down, thresholds


def find_thresholds(sequences):
    """" return a list of thresholds from a list of sequences after a blast"""
    x, y = get_blast_xy(sequences)
    _, _, _, _, thresholds = find_inflexion_points(x, y, 20)
    return thresholds

  
def graph(sequences):
    x, y = get_blast_xy(sequences)
    dy, dy2, inflexion_up, inflexion_down, thresholds = find_inflexion_points(x, y, 20)
    graph_threshold(x, y, dy, dy2, inflexion_up, inflexion_down)

    
def define_families(sequences, thresholds):
    """Retourne une liste où chaque élément est une classe PrteinFamily qui contient la séquence de la famille"""

    filtered_sequences = filter_blast_result(sequences, 10)
    families = []
    prev_treshold = 0

    for threshold in thresholds:
        family_sequence = [seq for seq in filtered_sequences if prev_treshold <= seq['e_value'] <= threshold]
        protein_family = ProteinFamily(family_sequence)
        prev_treshold = threshold
        families.append(protein_family)

    return families
