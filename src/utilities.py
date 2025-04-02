""" File for utility functions of all kinds"""

import os
import shutil
import ctypes
import requests
from RequestTools import *
import re
import numpy as np
import matplotlib.colors as mcolors
from PIL import Image, ImageDraw
from PDBTools import *
import json
import win32clipboard
import io
from collections import Counter
import matplotlib.pyplot as plt
from pywaffle import Waffle
import matplotlib.patches as mpatches


def extract_protein_name(fasta: str):
    """
    :param fasta: content of fasta file
    :return: the common name of the protein
    """
    pattern = r'\|([^|]+)\|([^OS]+)OS=([^ ]+)'  # Regular expression searching for the name
    match = re.search(pattern, fasta)
    if match:
        identifier = match.group(2).strip()
        return ' '.join(identifier.split()[1:])
    else:
        return None


def extract_fasta_sequence(fasta):
    """ Return the amino acid sequence as a string from the fasta file """
    inverted_sequence = fasta[::-1]
    i = 0
    aa_sequence, char = inverted_sequence[i], inverted_sequence[i]
    while not char.isnumeric():  # while we don't reach the SV value
        aa_sequence += char
        i += 1
        char = inverted_sequence[i]
    aa_sequence = aa_sequence[::-1].replace('\n', '')  # Remove residuals line returns
    return aa_sequence


def extract_fasta_data(fasta):
    """ returns the name and specie from a fasta file """
    name_pattern = r'\s(.*?)(?=\sOS=)'  # Everything between the first space and os
    name = re.findall(name_pattern, fasta)[0][:-1]  # /!\ Handle the case where no find
    specie_pattern = r'OS=(.*?)\sOX='  # Everything between os and ox
    specie = re.findall(specie_pattern, fasta)[0]
    return name, specie


def clear_bin():
    """ Remove all temp files located in bin"""
    for file in os.listdir('bin'):
        if file != ".gitkeep":
            os.remove(os.path.join('bin', file))


def download_img(url: str, output_path: str):
    """ download an image from its url"""
    print(f'url : {url}')
    response = requests.get(url, headers={'User-Agent': 'Uniprot retriever guillaume.boulet@grenoble-inp.org'})
    with open(output_path, 'wb') as file:
        file.write(response.content)


def open_img(path: str):
    """open locally an image"""
    os.startfile(path)


def open_specie_img(specie: str):
    url = retrieve_image_url(specie)  # Get the URL of the image from the wikipedia page of the specie
    path = os.path.join('bin', f'{specie}.png')  # Output path from the image
    download_img(url, path)
    open_img(path)


def download_and_get_image_path(prot):
    """ find an image on wikipedia of the associated organism and returns its path."""
    path = os.path.join(prot.path, f'{prot.organism}.png')  # Output path from the image
    if 'no_photo_placeholder.txt' in os.listdir(prot.path):
        return None
    if not f'{prot.organism}.png' in os.listdir(prot.path):
        print('Downloading photo...')
        url = retrieve_image_url(prot.organism)  # Get the URL of the image from the wikipedia page of the specie
        if url is None:
            return url
        else:
            download_img(url, path)
            return path
    else:
        print('Already have photo...')
        return path


def center_window(x, y, sx, sy):
    """ return the geometry for centering the window on the screen"""
    pos_x = (sx - x) // 2
    pos_y = (sy - y) // 2
    scaleFactor = int(ctypes.windll.shcore.GetScaleFactorForDevice(0) / 100)
    return '{}x{}+{}+{}'.format(x, y, pos_x*scaleFactor, pos_y*scaleFactor)


def hex_to_rgb(hex_color):  # conversion hexadecimal color to rgb
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i + 2], 16) for i in (0, 2, 4))


def rgb_to_hex(rgb_color):  # same shit but opposite
    return '#%02x%02x%02x' % rgb_color


def get_color_gradient(color1, color2, n):
    """ Generate n colors between color 1 and color 2"""
    rgb_color1 = hex_to_rgb(color1)
    rgb_color2 = hex_to_rgb(color2)
    rgb_color1 = np.array(rgb_color1) / 255.0  # Normalization
    rgb_color2 = np.array(rgb_color2) / 255.0
    r_gradient = np.linspace(rgb_color1[0], rgb_color2[0], n)  # Gradient creation
    g_gradient = np.linspace(rgb_color1[1], rgb_color2[1], n)
    b_gradient = np.linspace(rgb_color1[2], rgb_color2[2], n)
    gradient = np.clip(np.vstack((r_gradient, g_gradient, b_gradient)), 0.0, 1.0).T  # combine rgb values
    gradient *= 255.0  # Conversion to RGB
    gradient = gradient.round().astype(int)  # Conversion in INT
    hex_colors = [rgb_to_hex(tuple(rgb)) for rgb in gradient]  # back in HEX

    return hex_colors


def get_hover_color(color):
    """ returns a slightly brighter color than color"""
    rgb_color = hex_to_rgb(color)
    brighter_rgb = tuple(
        min(255, int(component + 25)) for component in rgb_color)  # Brighten by addind a bit to all rgb
    return rgb_to_hex(brighter_rgb)


def remove_content_between_parentheses(input_string):
    pattern = r'\([^)]*\)'
    result = re.sub(pattern, '', input_string)
    return result


def concatenate_infos(names, infos):
    """ format data for display in ProteinInfoFrame"""
    return [f'{name} : {info}' for name, info in zip(names, infos)]


def get_constrained_image_dimensions(constraints, image_path):
    """ returns new image size to fit in the constraints"""
    with Image.open(image_path) as img:
        image_size = img.size
    if constraints[0] >= image_size[0] and constraints[1] >= image_size[1]:  # If image fits in constraints
        return image_size
    else:
        ratio = max(image_size[0] / constraints[0], image_size[1] / constraints[1])
        return image_size[0] / ratio, image_size[1] / ratio


def create_and_save_render(prot):
    """ retrieve the pdb file for the prod_id from alpha fold website, build the render and return the save path"""
    if not f'{prot.prot_id}_render.png' in os.listdir(prot.path):
        print('Building render...')
        pdb_file = get_pdb_file(prot)  # Path to the downloaded pdb file
        exec_pymol_instruction(pdb_file, get_pymol_instructions('basic_render_script'))
        output_path = move_temp_render(prot.prot_id)
        return output_path
    else:
        print('Already have render')
        return os.path.join(prot.path, f'{prot.prot_id}_render.png')


def move_temp_render(prot_id):
    """ move the temporary pymol render to the associated protein directory"""
    source = os.path.join('bin', f'{prot_id}_render.png')
    target_dir = os.path.join('data/proteins', str(prot_id))
    target_file = target_dir + f'/{prot_id}_render.png'
    os.rename(os.path.join('bin', 'temp.png'), source)
    try:
        shutil.move(source, target_dir)
    except shutil.Error:
        os.remove(target_file)
        shutil.move(source, target_dir)
    return target_file


def read_fasta_file(path):
    with open(path, 'r') as file:
        return file.read()


def write_fasta(prot, content):
    with open(os.path.join(prot.path, f'{prot.prot_id}.fasta.txt'), 'w') as file:
        file.write(content)


def auto_crop(image_path, threshold=0):
    """ Automatically crops the image to remove excess blank space around the object """
    img = Image.open(image_path)
    img_gray = img.convert('L')
    bbox = img_gray.getbbox()
    bbox_expanded = (bbox[0] - threshold, bbox[1] - threshold, bbox[2] + threshold, bbox[3] + threshold)
    cropped_img = img.crop(bbox_expanded)
    cropped_img.save(image_path)
    return image_path


def find_key_by_value(dictionary, value):
    for key, val in dictionary.items():
        if val == value:
            return key
    return None


def extract_features(features: dict):
    """ returns a list of [alpha-helix/ beta sheet, (start aa, end aa)] from a loaded uniprot json file"""
    structures = []
    for feature in features:
        feature_type = feature['type']
        start = feature['location']['start']['value']
        end = feature['location']['end']['value']
        structures.append((feature_type, (start, end)))
    return structures


def write_json(data, output):
    with open(output, 'w') as f:
        json.dump(data, f)


def create_features_strand(features, factor=1):
    """ Creates a PIL image of a protein features to paste along a graph"""
    image_height = 20  # Height of the bands
    image = Image.new("RGB", (369, image_height), (255, 255, 255))  # 369 determined by mpl constant
    draw = ImageDraw.Draw(image)

    for struct_type, (start, end) in features:
        rect_color = (25, 25, 25) if struct_type == 'Beta strand' else (236, 71, 118)
        draw.rectangle([(int(start * factor), 0), (int(end * factor), image_height)], fill=rect_color)

    return image


def open_saved_json(file_path):
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        return data
    except FileNotFoundError:
        return None
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON: {e}")
        return None


def merge_feat_and_graph(feat, graph):
    """ Place both vertical and horizontal features images on a mpl graph /!\ Square format for the mpl"""
    y_offset = 38
    x_offset = 108
    merged_img = graph.copy()
    merged_img.paste(feat, (x_offset, y_offset))  # Horizontal bar

    rotated_feat = feat.rotate(270, expand=True)
    rotated_img = Image.new("RGBA", merged_img.size)
    rotated_img.paste(rotated_feat, (x_offset - 21, y_offset + 20))
    merged_img = Image.alpha_composite(merged_img.convert("RGBA"), rotated_img)
    # merged_img.paste(Image.open('ressources/legend.png').resize((82, 43)), (550, 6))  # Manual legend

    return merged_img


def get_fasta_path(prot):
    if not f'{prot.prot_id}.fasta.txt' in os.listdir(prot.path):
        print('Retrieving fasta...')
        answer = retrieve_protein_entry(prot.prot_id)
        if isinstance(answer, int):
            print(f'Fasta Error : {answer}')
        else:
            write_fasta(prot, answer)
            fasta = answer
    else:
        print('Already have fasta')
        fasta = read_fasta_file(os.path.join(prot.path, f'{prot.prot_id}.fasta.txt'))
    return fasta


def download_and_save_features(prot):
    if not f'{prot.prot_id}_features.json' in os.listdir(prot.path):
        print('Retrieving features...')
        answer = get_structural_features(prot.prot_id)
        if isinstance(answer, int):
            print(f'Features Error : {answer}')
        else:
            prot.features_path = os.path.join(prot.path, f'{prot.prot_id}_features.json')
            write_json(answer, prot.features_path)
            return prot.features_path
    else:
        print('Already have features')
    return prot.features_path


def open_taxonomy(prot):
    try:
        return open_saved_json(prot.taxo_path)
    except FileNotFoundError:
        return None


def write_placeholder(prot):
    with open(os.path.join(prot.path, 'no_photo_placeholder.txt'), 'w') as file:
        file.write('Placeholder file to not trying to download a non-existing photo.')


def send_to_clipboard(image):
    """ copy the PIL image to the clipboard"""
    output = io.BytesIO()
    image.convert('RGB').save(output, 'BMP')
    data = output.getvalue()[14:]
    output.close()

    win32clipboard.OpenClipboard()
    win32clipboard.EmptyClipboard()
    win32clipboard.SetClipboardData(win32clipboard.CF_DIB, data)
    win32clipboard.CloseClipboard()
    print('Copied')


def get_local_proteins():
    return os.listdir('data/proteins')

  
def display_histogram(data, bins=10, xlabel="Value", ylabel="Frequency", title="Histogram"):
    """ Display a histogram of a series of values. """
    plt.hist(data, bins=bins)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.show()


def extract_seq_data(title):
    """ return the name and the organism of a protein in a sequence """
    pattern1 = r'RecName: Full=([^;]+)'
    pattern2 = r'\[([^\]]+)\]'
    name, specie = '', ''
    match1 = re.search(pattern1, title)
    match2 = re.search(pattern2, title)
    if match1:
        name = match1.group(1)
    if match2:
        specie = match2.group(1)
    return name, specie


def most_common_names_percentages(names, n=5):
    """ Returns the n most common names from a list with their occurrences as percentages. """
    total_names = len(names)
    name_counter = Counter(names)
    most_common = name_counter.most_common(n)
    most_common_percentages = {name: (count / total_names) * 100 for name, count in most_common}
    return most_common_percentages


def plot_waffle_chart(data, prot_name, output_path):
    plt.figure(
        FigureClass=Waffle,
        values=data,
        rows=10, columns=10,
        legend={'loc': 'upper left', 'bbox_to_anchor': (1, 1)}
    )

    plt.title(f'Family distribution of {prot_name} blast')
    fig = plt.gcf()
    fig.canvas.draw()
    image_np = np.array(fig.canvas.renderer.buffer_rgba())
    image_pil = Image.fromarray(image_np)
    image_pil.save(output_path)
    return output_path


def fill_majoritary_dictionnary(freq):
    present_freq = sum(val for val in freq.values())
    filling = 100 - present_freq
    freq['others'] = filling
    return freq


def segment_background(seuils, voisins, plt, thresholds, colors):
    """ Segments the background of the plot based on the provided thresholds and colors. """
    c = 0
    for i in range(len(thresholds) - 1):
        start = thresholds[i]
        end = thresholds[i + 1]
        color = colors[i % len(colors)]  # Use modulus to repeat colors if necessary
        # Change the background color of the section between start and end
        plt.axvspan(start, end, color=color, alpha=0.5)
        lower_threshold = thresholds[i]
        upper_threshold = thresholds[i + 1]
        segment_indices = [index for index, seuil in enumerate(seuils) if lower_threshold <= seuil <= upper_threshold]
        segment_seuils = [seuils[index] for index in segment_indices]
        segment_voisins = [voisins[index] for index in segment_indices]
        plt.plot(segment_seuils, segment_voisins, marker='.', color=color,
                 label=f'Threshold: {lower_threshold}-{upper_threshold}')
        c += 1
    handles = []
    used_colors = colors if c > len(colors) else colors[:c]
    for i, color in enumerate(used_colors, start=1):
        handle = mpatches.Rectangle((0, 0), 1, 1, color=color, alpha=0.5, label=f'Family {i}')
        handles.append(handle)
    plt.legend(handles=handles, loc='upper left')


def truncate_data(x, y, tolerance):
    """ Truncates data x and y starting from the index where y values are no longer constant. """
    start = y[0]
    max_value = tolerance*max(y)
    for index, value in enumerate(y):
        print(value, value-start, max_value)
        if abs(value-start) > max_value:
            return x[index:], y[index:]
    return x, y

