""" File containing tools for request to Uniprot"""

import requests
import json  # Serializing
import wikipediaapi  # Accessing wikipedia website
from bs4 import BeautifulSoup  # Parsing HTML documents
import webbrowser


def retrieve_protein_entry(uniprot_id: str):
    """ Return the fasta data associated to the ID of the protein uniprot_id"""
    api_url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"  # UniProt API URL
    response = requests.get(api_url)  # Get the response from Uniprot
    if response.status_code == 200:  # If the request is successful
        return response.text  # Return the protein sequence
    else:
        #  print(f"Error retrieving protein entry. Status code: {response.status_code}")
        return response.status_code


def retrieve_taxonomy(uniprot_id: str):
    """ Return a list of species linked to the protein of uniprot_id"""
    print('Retreiving taxonomy...')
    api_url1 = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}/"   # UniProt API URL
    response = requests.get(api_url1)
    if response.status_code == 200:  # If the request is successful
        return response.json()['organism']['lineage']  # Conversion to dictionnary and access by organism and lineage keys
    else:
        #  print(f"Error retrieving protein taxonomy. Status code: {response.status_code}")
        return response.status_code


def retrieve_image_url(specie: str):
    headers = {'User-Agent': 'Uniprot retriever guillaume.boulet@grenoble-inp.org'}
    wikipedia_url = f'https://en.wikipedia.org/wiki/{specie}'
    html_content = requests.get(wikipedia_url, headers=headers).text  # HTML code of the wiki page
    parser = BeautifulSoup(html_content, 'html.parser')  # Parser object to inspect the HTML page
    infobox = parser.find('table', {'class': 'infobox'})  # Finding the table of the wikipedia page
    if infobox:  # If found
        image_div = infobox.find('img')  # Image div of the main page image
        if image_div and 'src' in image_div.attrs:  # If the image is found and it has source
            return 'https:' + image_div['src']  # Get the source
    return None


def check_for_prot(protein_id: str):
    """ Check if the protein exists on uniprot"""
    api_url = f"https://www.uniprot.org/uniprot/{protein_id}.fasta"  # UniProt API URL
    try:
        response = requests.get(api_url)
        return response.status_code == 200
    except requests.RequestException:  # Erreur
        return False


def get_pdb_url_from_uniprot_id(uniprot_id: str):
    """ explicit enough"""
    if '_' in uniprot_id:
        uniprot_id = uniprot_id.split('_')[0]
    print(uniprot_id)
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"  # Url of the protein on alphafold website
    response = requests.get(api_url)
    if response.status_code == 200:  # If the request is successful
        return response.json()[0]['pdbUrl']  # Conversion to dictionnary and access by the right key
    else:
        return response.status_code


def get_rank_data(rank: str):
    """ retrieve a short description of one and oly rank in a taxonomy"""
    headers = {'User-Agent': 'Uniprot retriever guillaume.boulet@grenoble-inp.org'}
    wiki = wikipediaapi.Wikipedia(headers['User-Agent'], 'en')
    page = wiki.page(rank)
    return page.summary


def get_taxonomy_data(taxo: list[str]):
    """ retrieve the wikipedia intro from all the rank in taxonomy. Made for faster execution better than N executions of get_rank_data"""
    descriptions = {}
    headers = {'User-Agent': 'Uniprot retriever guillaume.boulet@grenoble-inp.org'}
    wiki = wikipediaapi.Wikipedia(headers['User-Agent'], 'en')
    for rank in taxo:
        descriptions[rank] = wiki.page(rank).summary
    return descriptions


def get_structural_features(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json?fields=ft_strand%2Cft_helix"
    headers = {"Accept": "application/json"}

    try:
        response = requests.get(url, headers=headers)
        response.raise_for_status()  # Raise an exception for 4xx or 5xx status codes
        data = response.json()

        structural_features = data.get('features', [])
        return structural_features

    except requests.exceptions.RequestException as e:
        print(f"Error fetching data: {e}")
        return None


