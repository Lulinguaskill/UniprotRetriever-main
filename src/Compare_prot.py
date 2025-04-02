from PIL import Image, ImageDraw, ImageFont
import numpy as np
from BasicClasses import *
from Bio.Align import substitution_matrices

# test:

# Première famille
seq_list1 = [
    "KKTAYIAKQRTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL",
    "KMTAYGIGGKLSWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL"
]
family1 = ProteinFamily(seq_list1)
family1.mean = 1.5
family1.aa_lengths = [367, 367]
family1.range = (0.9, 2.1)
family1.species = ["Homo sapiens", "Mus musculus"]
family1.names = ["Protein A", "Protein B"]
family1.names_percentages = {"Protein A": 60, "Protein B": 40}

# Deuxième famille
seq_list2 = [
    "AAGSSVSADPALPEDLRFAFAWAGWLLLLFSGAVAAEAGSGEVEVFAQVAQFQVKGFGGRGAAVLGPDAAEAGVAVGAVAGVLGGTVAEALAIAG",
    "AGSSVSADPALPEDLRFAFAWAGWLLLLFSGAVAAEAGSGEVEVFAQVAQFQVKGFGGRGLAPPEQARERQVAQLFGLFGLGDDEEEEEEEEEEEE"
]
family2 = ProteinFamily(seq_list2)
family2.mean = 0.8
family2.aa_lengths = [310, 310]
family2.range = (0.5, 1.2)
family2.species = ["Escherichia coli", "Bacillus subtilis"]
family2.names = ["Protein X", "Protein Y"]
family2.names_percentages = {"Protein X": 70, "Protein Y": 30}

# Troisième famille
seq_list3 = [
    "MMSLHAGRGVVTFASMRDNFAGLVQLAWFDSRK"
    "MSLHAGRGVVTFASMRDNFAGLVQLAWFDSRKK"]
family3 = ProteinFamily(seq_list3)
family3.mean = 2.2
family3.aa_lengths = [467, 467]
family3.range = (1.8, 2.6)
family3.species = ["Arabidopsis thaliana", "Oryza sativa"]
family3.names = ["Protein M", "Protein N"]
family3.names_percentages = {"Protein M": 80, "Protein N": 20}

# Create a list of families
family_list = [family1, family2, family3]

prot_ref="KKTAYIAKQRTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL"
### FIN TEST

# demander à l'utilisateur de donner i et j pour choisir quelle famille comparer avec query prot.
def extraire_2prot_diff_families(family_list, i, j): #Prend 2 famille et return le nom et la sequence de 2 prot pris au hasard au sein de ces 2familles
    family1 = family_list[i - 1]
    family2 = family_list[j - 1]
    print(f"Comparing Family {i} with Family {j}:")
    sequence1 = random.choice(family1.sequences)
    sequence2 = random.choice(family2.sequences)
    # Obtenez les index des séquences aléatoires
    index1 = family1.sequences.index(sequence1)
    index2 = family2.sequences.index(sequence2)
    # Obtenez les noms correspondant aux séquences choisies aléatoirement
    name1 = family1.names[index1]
    name2 = family2.names[index2]
    print("Random Sequence from Family", i, ":", name1)
    print(sequence1)
    print("Random Sequence from Family", j, ":", name2)
    print(sequence2)
    print()  # Print a newline for clarity
    return name1, name2, sequence1, sequence2

#creation des sauvegarde des scores et des etats pour chaque déplacement
all_states = []
all_score = []

def indice_first_letter(prot): #return l'indice de la premiere lettre d'une chaine de carcatere
    i = 0
    while prot[i] == '-':
        i = i + 1
    return i
def indice_last_letter(prot): #return l'indice de la dernière lettre d'une chaine de carcatere
    indice = len(prot) - 1  # Indice du dernier caractère dans la chaîne
    while indice >= 0 and prot[indice] == '-':
        indice -= 1
    return indice

#creer un compteur de score qui s'indente dès que 2 caractere sont les mêmes et qu'ils sont adjacents
def compte_score(prot_i, prot_ref):
    # Initialiser le nombre de lettres côte à côte
    count = 0
    # Parcourir les deux chaînes simultanément jusqu'à ce que l'une d'elles soit entièrement parcourue
    ind = 0
    while ind < len(prot_i) and ind < len(prot_ref):
        # Vérifier si les caractères actuels sont les mêmes et différents de '-'
        if prot_i[ind] == prot_ref[ind] and prot_i[ind] != '-':
            # Augmenter le nombre de lettres côte à côte
            count += 1
            # Chercher combien de caractères identiques suivent à la suite
            while (ind+1 < len(prot_i) and ind+1 < len(prot_ref)) and prot_i[ind+1] == prot_ref[ind+1] and prot_i[ind+1] != '-':
                count += 1
                ind += 1
            return count
        # Passer au caractère suivant
        ind += 1
    return count

def save_state(prot, prot_fix):
    state = [prot, prot_fix]
    score = [compte_score(prot, prot_fix)]
    return state, score
#Trouve l'état avec le score le plus grand:
def best_align(all_states, all_score):
    ind_max = all_score.index(max(all_score))
    return all_states[ind_max][1], all_states[ind_max][0]
#Ajoute le score et l'état dans les listes pour les sauvegarder
def indentation_state_score(prot, prot_ref):
    state, score = save_state(prot, prot_ref)
    all_score.append(score)
    all_states.append(state)
    return all_states, all_score


def total_moves(prot, prot_ref): #Aligne REF et Prot_i
    while indice_last_letter(prot) != indice_first_letter(prot_ref):
        prot_ref = '-'+prot_ref
    while len(prot)!= len(prot_ref):
        prot = prot+'-'
    all_states=[]
    all_score = []
    while indice_first_letter(prot)!= len(prot_ref):
        prot='-'+prot
        all_states, all_score = indentation_state_score(prot,prot_ref)
    while len(prot)!= len(prot_ref):
        prot_ref = prot_ref+'-'
    align_ref_for_j, best_align_i = best_align(all_states, all_score)
    all_score.clear()
    all_states.clear()
    return (prot_ref, best_align_i)

#Aligne prot_i et prto_ref avec prot_j
def second_prot(prot_i,prot_j, prot_ref):
    prot_ref, prot_i = total_moves(prot_i,prot_ref)
    while indice_first_letter(prot_ref) < indice_last_letter(prot_j):
        prot_i = '-'+prot_i
        prot_ref = '-'+prot_ref
    while len(prot_j)!= len(prot_i):
        prot_j = prot_j+'-'
    while indice_first_letter(prot_j) != indice_last_letter(prot_ref):
        all_states,all_score = indentation_state_score(prot_j,prot_ref)
        prot_j='-'+prot_j
    align_ref_for_j,best_align_j = best_align(all_states, all_score)
    return (prot_ref, best_align_j, prot_i)

namei, namej, prot_i, prot_j = extraire_2prot_diff_families(family_list, 1, 2)
prot_ref, best_align_j, prot_i = second_prot(prot_i,prot_j, prot_ref)

def add_tirets(*proteins): #Ajoute les tirets après les prot pour avoir la même taille
    max_length = max(map(len, proteins))
    return tuple(p.ljust(max_length, '-') for p in proteins)

def cancel_tirets(prot1, prot2, prot3): #retire les tireits inuitles à gauche et à droite
    prot1, prot2, prot3 = add_tirets(prot1, prot2, prot3)
    first_indices = [indice_first_letter(prot) for prot in (prot1, prot2, prot3)]
    shift_amount = min(first_indices)
    prot1, prot2, prot3 = prot1[shift_amount:], prot2[shift_amount:], prot3[shift_amount:]

    prot_list = [prot1, prot2, prot3]
    last_indices = [indice_last_letter(prot) for prot in prot_list]
    longest_index = last_indices.index(max(last_indices))
    longest_prot = prot_list[longest_index]

    while longest_prot.endswith('-'):
        prot1, prot2, prot3 = prot1[:-1], prot2[:-1], prot3[:-1]
        longest_prot = longest_prot[:-1]

    return prot1, prot2, prot3

prot_i, best_align_j, prot_ref=cancel_tirets(prot_i,best_align_j,prot_ref)

def write(x,y,prot,draw,font,color,ind): #ecrit chaque caractère
    char_bbox = draw.textbbox((x, y), prot[ind], font=font)  # Obtenir la boîte englobante du caractère
    char_width = char_bbox[2] - char_bbox[0]  # Largeur du caractère
    draw.text((x, y), prot[ind], fill=color, font=font)
                      # Incrémenter x en fonction de la largeur du caractère
    return char_width

def dessin(prot_i,prot_ref,draw,x,y,l,font): #Dessine prot_ref et prot_i
    nbr= compte_score(prot_i, prot_ref)
    # Dessiner une ligne de text
    ind=0
    y2 = y+l
    x2 = x+3.2
    for char in prot_i:
        if char == prot_ref[ind] and char !='-': #Si element en commun on les écrit en rouge
            for p in range (ind,nbr+ind):
                char_width = write(x, y, prot_i, draw, font, "red", p)
                x += char_width
                char_width = write(x2, y2, prot_ref, draw, font, "red", p)
                x2 += char_width
                p= 1+p
            for k in range (ind+nbr,len(prot_i)): #quand les sequences communes sont finit on ecrit le reste en noir
                char_width = write(x, y, prot_i, draw, font, "black", k)
                x += char_width  # Incrémenter x en fonction de la largeur du caractère

                char_width = write(x2, y2, prot_ref, draw, font, "black",k)
                x2 += char_width  # Incrémenter x en fonction de la largeur du caractère
                k=k+1
            return

        else: #si pas communs : on écrit en noir
            char_width = write(x, y, prot_i, draw, font, "black", ind)
            x += char_width  # Incrémenter x en fonction de la largeur du caractère

            char_width = write(x2, y2, prot_ref, draw, font, "black", ind)
            x2 += char_width  # Incrémenter x en fonction de la largeur du caractère
        ind+=1

def dessin2(prot,prot_ref,draw,x,y,l,font):  #Dessine prot_j
    nbr = compte_score(prot, prot_ref)
    # Dessiner une ligne de text
    ind = 0
    y2 = y-l
    for char in prot:
        if char == prot_ref[ind] and char !='-': #Si element en commun on les écrit en rouge
            for p in range(ind , nbr+ind):
                char_width = write(x, y, prot, draw, font, "red", p)
                x += char_width  # Incrémenter x en fonction de la largeur du caractère

                x2 = x - char_width
                char_width2 = write(x2, y2, prot_ref, draw, font, "red", p)
                draw.rectangle([x2, y2, x2 + char_width2, y2 + 18], fill="white") #RECTANGLE BLANC sur prot_ref pour pas superposer les lettres
                char_width2 = write(x2, y2, prot_ref, draw, font, "red", p) #ecrit les lettres communes en rouge
                x2 += char_width2  # Incrémenter x en fonction de la largeur du caractère
                p= 1+p
            for k in range (ind+nbr,len(prot)): #quand les sequences communes sont finit on ecrit le reste en noir
                char_width = write(x, y, prot, draw, font, "black",k)
                x += char_width  # Incrémenter x en fonction de la largeur du caractère
                k=k+1
            return
        else:#si pas communs : on écrit en noir
            char_width = write(x, y, prot, draw, font, "black", ind)
            x += char_width  # Incrémenter x en fonction de la largeur du caractère
        ind+=1


def create_text_image(prot_i, prot_ref, best_align_j, font_size, output_file, font_path, nom_i, nom_ref,
                                 nom_j):
    # Calculer la longueur maximale des noms de protéines
    longueur_max_nom = max(len(nom_i), len(nom_ref), len(nom_j))

    # Calculer la largeur nécessaire pour les noms de protéines
    largeur_nom = longueur_max_nom * font_size

    # Calculer les coordonnées de départ pour les noms et les séquences
    x_nom = 10
    x_seq = largeur_nom   # Ajouter un peu de marge entre les noms et les séquences
    y = 10
    l = 18

    # Créer une nouvelle image avec un fond blanc
    largeur = largeur_nom + 10 * x_seq
    hauteur = 70
    image = Image.new("RGB", (largeur, hauteur), "white")
    draw = ImageDraw.Draw(image)

    # Définir la police
    font = ImageFont.truetype(font_path, font_size)

    # Écrire les noms de protéines
    draw.text((x_nom, y), nom_i, fill="black", font=font)
    draw.text((x_nom, y + 1 * l), nom_ref, fill="black", font=font)
    draw.text((x_nom, y + 2 * l), nom_j, fill="black", font=font)

    # Dessiner chaque ligne de texte pour les séquences
    dessin(prot_i, prot_ref, draw, x_seq, y, l, font)
    dessin2(best_align_j, prot_ref, draw, x_seq, y + 2 * l, l, font)

    # Enregistrer l'image au format PNG
    image.save(output_file)
    image.show()

# Spécifiez le chemin d'accès à votre police TrueType (.ttf)
font_path = "C:\\Users\\antoi\\PyCharm\\Uniprot\\UniprotRetriever\\Comfortaa,Courier_Prime,Zilla_Slab\\Courier_Prime\\CourierPrime-Regular.ttf"
print('res')
print(prot_i)
print(prot_ref)
print(best_align_j)
# Exemple d'utilisation
font_size = 20
output_file = "text_image.png"
create_text_image(prot_i,prot_ref,best_align_j, font_size, output_file,font_path, namei,"query protein",namej)








