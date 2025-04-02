import matplotlib.pyplot as plt

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

""" Main for Uniprot Retriever. This file is made to be called from the root directory of the project, not from src.
 Consider change the interpreter working directory into ./UniprotRetriever """


if __name__ == '__main__':
    create_DB()  # Create DB if non-existent
    clear_bin()  # Emptying temp files
    app = MainApp()  # UI starts
