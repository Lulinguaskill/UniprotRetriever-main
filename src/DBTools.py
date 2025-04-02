"""" File for managing local protein database"""
import os
import sqlite3
from sqlite3 import *


def create_DB():
    """Creates the protein database if not existing """
    try:
        if not 'proteins.db' in os.listdir('data'):  # Check if db is existing
            conn = connect('data/proteins.db')
            cursor = conn.cursor()
            with open('data/queries/create_query.sql', 'r') as query_file:  # Loading the query from SQL file
                query = query_file.read()
            cursor.execute(query)
            conn.commit()
            conn.close()
            print('Proteins database created.')
        else:
            print('Proteins database already existing.')
    except sqlite3.Error as e:  # Handles SQL errors
        print(f'Database error : {e}')


def insert_protein(line):
    """ Insert a protein line in the associated database """
    try:
        conn = connect('data/proteins.db')
        cursor = conn.cursor()
        with open('data/queries/insert_query.sql', 'r') as query_file:
            query = query_file.read()
        cursor.execute(query, line)
        conn.commit()
        conn.close()
        print(f'Line {line[0]} inserted.')
        return True
    except sqlite3.Error as e:
        print(f'Database error : {e}')
        return False
