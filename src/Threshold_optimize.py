import matplotlib.pyplot as plt
import blastUtils
import blast_cmd
import numpy as np
from blastUtils import filter_blast_result
from PIL import Image
from inflexion_points import find_thresholds, find_inflexion_points
from utilities import segment_background

colors = ["#ac3931", "#e5d352", "#d9e76c", "#537d8d", "#482c3d"]


def get_blast_xy(sequences):
    voisins = []
    seuils = np.logspace(-150, 1, 300)

    for threshold in seuils:
        voisins.append(len(filter_blast_result(sequences, threshold)))
    return seuils, voisins


def dy_graph_threshold(seuils, voisins, dy, dy2, inflexion_up, inflexion_down):
    plt.semilogx(seuils, voisins, marker='.')  # Tracer le graphique avec des marqueurs circulaires
    for point in inflexion_up:
        plt.scatter(point[0], point[1], marker='o', color='red', label="Point d'inflexion")
    # for point in inflexion_down:
    #     plt.scatter(point[0], point[1], color='red', label="Point d'inflexion")
    # plt.twinx()
    # plt.plot(seuils, dy, 'y', label="Dérivée")
    plt.twinx()
    plt.plot(seuils, dy2, 'r', label="Dérivée seconde")
    plt.xlabel('Blast E-value')  # Libellé de l'axe x
    plt.ylabel('Number of hits')  # Libellé de l'axe y
    plt.title('Nombre de voisins trouvés en fonction du seuil e-value')  # Titre du graphique
    plt.xscale('log')  # Échelle logarithmique pour l'axe x (optionnel)
    plt.grid(True)  # Activer la grille (optionnel)
    plt.show()  # Afficher le graphique


def graph_threshold(seuils, voisins, output_path, prot_id, prot_name):
    plt.figure()
    #plt.plot(seuils, voisins, marker='.', color='black')  # Tracer le graphique avec des marqueurs circulaires
    plt.xlabel('Blast E-value')  # Libellé de l'axe x
    plt.ylabel('Number of hits')  # Libellé de l'axe y
    plt.title(f'Blast graph of {prot_id} : {prot_name}')  # Titre du graphique
    plt.xscale('log')  # Échelle logarithmique pour l'axe x (optionnel)
    plt.grid(True)  # Activer la grille (optionnel)
    t = find_inflexion_points(seuils, voisins, 20)[4]
    t.insert(0, min(seuils))  # ajout de la valeur à gauche du graphe
    segment_background(seuils, voisins, plt, t, colors)
    fig = plt.gcf()
    fig.canvas.draw()
    image_np = np.array(fig.canvas.renderer.buffer_rgba())
    image_pil = Image.fromarray(image_np)
    image_pil.save(output_path)
    return output_path


def colored_graph_threshold(seuils, voisins, thresholds, colors):
    plt.figure()
    plt.plot(seuils, voisins, marker='.', color='black')  # Tracer le graphique avec des marqueurs circulaires
    plt.xlabel('Seuil e-value')  # Libellé de l'axe x
    plt.ylabel('Nombre de voisins')  # Libellé de l'axe y
    plt.title('Nombre de voisins trouvés en fonction du seuil e-value')  # Titre du graphique
    plt.xscale('log')  # Échelle logarithmique pour l'axe x (optionnel)
    plt.grid(True)  # Activer la grille (optionnel)

    # Diviser les données en segments en fonction des seuils et colorer chaque segment
    for i in range(len(thresholds) - 1):
        lower_threshold = thresholds[i]
        upper_threshold = thresholds[i + 1]
        segment_indices = [index for index, seuil in enumerate(seuils) if lower_threshold <= seuil <= upper_threshold]
        segment_seuils = [seuils[index] for index in segment_indices]
        segment_voisins = [voisins[index] for index in segment_indices]
        segment_color = colors[i % len(colors)]  # Use modulo to cycle through colors
        plt.plot(segment_seuils, segment_voisins, marker='.', color=segment_color,
                 label=f'Threshold: {lower_threshold}-{upper_threshold}')
    plt.legend()
    plt.show()  # Afficher le graphique


def find_breakpoint_x(x_values, y_values):
    """
    Find the x value for the breakpoint in a curve.

    Args:
    x_values (list or array): List of x values.
    y_values (list or array): List of y values.

    Returns:
    float: The x value for the breakpoint.
    """
    # Calculate the differences in y values
    diffs = [y_values[i + 1] - y_values[i] for i in range(len(y_values) - 1)]

    # Find the index where the difference is maximum
    max_diff_index = diffs.index(max(diffs))

    # Return the corresponding x value
    return x_values[max_diff_index]
