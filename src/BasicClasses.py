""" File defining classes to store the structure information of proteins, sequences, etc """
import threading

from DBTools import *
from RequestTools import *
from Bio.Seq import Seq
from utilities import *
from webbrowser import open
from ThreeDimTools import *
import random
from blastUtils import get_blast_result, blast
from Threshold_optimize import graph_threshold, get_blast_xy
from config_manager import *

db_path = get_config('config.yml')['db_path']


""">sp|P32387|RM41_YEAST Large ribosomal subunit protein uL23m 
OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=MRP20 PE=1 SV=1
MPRLTVGTKNMLYPLQKTLAVGSCKPEQVPIRSLASVVESSSKILDKSGSDREVDINVSE
KIYKWTKAGIEQGKEHFKVGGNKVYFPKARIILLRPNAKHTPYQAKFIVPKSFNKLDLRD
YLYHIYGLRAMNITTQLLHGKFNRMNLQTTRFREPQIKKMTIEMEEPFIWPEEPRPDENS
FWDSTTPDNMEKYREERLNCLGSDANKPGTAFDGVVGPYERVAQPFIPRFLKREIDNKRE
RHAAELQRADKLIALNRYIEDLH"""


class Protein:
    """ Class to represent a protein """

    def __init__(self, prot_id, ui=None):

        self.prot_id = prot_id  # Uniprot ID
        self.fasta = None  # Fasta file as str from uniprot
        self.name = None  # Bio name
        self.organism = None  # Specie with this protein
        self.sequence = None  # ProteinSequence object
        self.taxonomy = None  # Taxonomy object
        self.photo_path = None
        self.render_path = None
        self.blast_graph_path = None
        self.path = None
        self.taxo_path = None
        self.pdb_path = None
        self.features_path = None
        self.render_2D_path = None
        self.ui = ui  # Bind to UI for progressbar
        self.build_directory()
        self.get_fasta()
        self.build_necessary()

    def build_directory(self):
        if self.ui is not None:
            self.ui.progress()
        self.path = os.path.join("data/proteins", f"{self.prot_id}")
        try:
            os.mkdir(self.path)
        except FileExistsError:
            pass

    def download_photo(self):
        self.photo_path = download_and_get_image_path(self)
        if self.photo_path is None:
            print('No photo can be found')
            write_placeholder(self)
            return os.path.join('ressources', 'no_photo.png')
        else:
            return self.photo_path

    def build_render(self):
        self.render_path = create_and_save_render(self)
        return auto_crop(self.render_path)

    def get_features(self):
        self.features_path = os.path.join(self.path, f'{self.prot_id}_features.json')
        download_and_save_features(self)
        return self.features_path

    def get_fasta(self):
        if self.ui is not None:
            self.ui.progress()
        self.fasta = get_fasta_path(self)
        return self.fasta

    def build_necessary(self):
        """ extract the only information necessary to build the dashboard """
        if self.ui is not None:
            self.ui.progress()
        self.name, self.organism = extract_fasta_data(self.fasta)
        self.sequence = Seq(extract_fasta_sequence(self.fasta))

    def get_info(self):
        return self.name, self.organism, len(self.sequence)  # Rapid summary of the prot for display

    def build_taxonomy(self):
        """ build the Taxonomy class for treeview """
        self.taxo_path = os.path.join(self.path, f'{self.prot_id}_taxo.json')
        tax = open_taxonomy(self)
        if tax is None:  # If taxo file not already present
            self.taxonomy = Taxonomy(self)
            self.taxonomy.build_taxonomy()
        else:
            print('Already have taxonomy')
            self.taxonomy = Taxonomy(self, tax)

    def blast(self):
        """ run a blast without threshold """
        blast(self.prot_id, db_path, 10)
        x, y = get_blast_xy(get_blast_result(self.prot_id))
        x, y = truncate_data(x, y, 0)
        self.blast_graph_path = graph_threshold(x, y, os.path.join(self.path, f'{self.prot_id}_blastgraph.png'),
                                                self.prot_id, self.name)
        return self.blast_graph_path

    def get_2D_render(self):
        self.render_2D_path = os.path.join(self.path, f'{self.prot_id}_2D_render.png')
        if not f'{self.prot_id}_2D_render.png' in os.listdir(self.path):  # If render not present, build it
            self.pdb_path = get_pdb_file(self)
            self.get_features()
            build_2D_render(self)
        return self.render_2D_path

    def save(self):
        """ Add a line in the database :  TO DO : Check if protein ID isn't already present"""
        if self.sequence and self.taxonomy:
            return insert_protein([self.prot_id, self.sequence.get_aa_sequence(), self.taxonomy.get_string_taxonomy()])

    def open_uniprot_webpage(self, event=None):
        open(f'https://www.uniprot.org/uniprotkb/{self.prot_id}/entry')  # Allows to go on uniprot entry from UI

    def copy_2D_render_to_clipboard(self):
        send_to_clipboard(Image.open(self.render_2D_path))

    def save_as(self, output_path):
        Image.open(self.render_2D_path).save(output_path)

    def copy_blastrender_to_clipboard(self):
        send_to_clipboard(Image.open(self.blast_graph_path))


# noinspection PyMissingConstructor
class Taxonomy(Protein):
    """ Basic class to represent the taxonomy of a protein"""

    def __init__(self, protein: Protein, tax=None):
        self.master = protein
        self.taxonomy = None
        self.specie = None
        if tax is not None:
            self.taxonomy = tax

    def build_taxonomy(self):
        self.taxonomy = retrieve_taxonomy(self.master.prot_id)  # Hides a call to Uniprot API

    def get_taxonomy(self):
        """ Return the taxonomy as a list of species"""
        return self.taxonomy

    def get_specie(self):
        """Return the specie at the end of protein taxonomy"""
        if self.specie is None:
            self.specie = self.taxonomy.keys()[-1]
        return self.specie

    def save_descr(self, descriptions):
        write_json(descriptions, self.master.taxo_path)


class Task:
    """ Class for running task in background. Allows to display a loading animation and stop it when ended"""

    def __init__(self, target_func):
        self.target_func = target_func
        self.ui = None
        self.output = None

    def start(self, ui=None):
        if ui is not None:
            self.ui = ui
        threading.Thread(target=self._start).start()

    def _start(self):
        self.output = self.target_func()
        if self.ui is not None:
            self.ui.end(self.output)


class ProteinFamily:
    """ class for storing multiple proteins which can be a part of a section of a blast result"""

    def __init__(self, seq_list: list):
        self.sequences = seq_list
        self.mean = None  # e-value mean of the family
        self.aa_lengths = []  # List of aa lenghts
        self.range = None  # Range of the e-valyes
        self.species = []  # List of species associated with the family proteins
        self.names = []  # List of names associated with the family proteins
        self.names_percentages = None  # Dictionnary of protein names with percentages

    def get_mean(self):
        mean = 0
        for seq in self.sequences:
            mean += seq['e_value']
        self.mean = mean/len(self.sequences)
        return self.mean

    def get_range(self):
        self.range = self.sequences[-1]['e_value'] - self.sequences[0]['e_value']
        return self.range

    def get_aa_dispersion(self):
        for seq in self.sequences:
            self.aa_lengths.append(seq['length'])
        return self.aa_lengths

    def display_aa_dispersion_graph(self):
        display_histogram(self.aa_lengths, 10, 'Length', 'Number of proteins',
                          'Dispersion of protein lenght in the family')

    def display_family_data(self):
        print(f'Range : {self.range}')
        print(f'Mean : {self.mean}')
        self.display_aa_dispersion_graph()

    def build_name_and_speices(self):
        if not self.names:
            for seq in self.sequences:
                result = extract_seq_data(seq['title'])
                self.names.append(result[0])
                self.species.append(result[1])

    def get_random_seq(self, number):
        """ returns random sequences from the family"""
        return random.sample(self.sequences, number)

    def get_majoritary_names(self, number):
        self.names_percentages = most_common_names_percentages(self.names, number)
        return self.names_percentages

    def build_waffle(self):
        path = plot_waffle_chart(fill_majoritary_dictionnary(self.get_majoritary_names(10)),
                                 'Family 1', 'bin/temp_render.png')
        return path

