import yaml
from yaml import load, dump, Loader


def get_config(config_path):
    fichier = open(config_path, 'r')
    return load(fichier, Loader)


def save_config(config, config_path):
    fichier = open(config_path, 'w')
    dump(config, fichier)