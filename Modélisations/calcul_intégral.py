# Bibliothèques locales

from constantes import *
from base_données import *

# Bibliothèques externes

from scipy.integrate import quad
from scipy.interpolate import interp1d
from random import random
from os import system, name

import numpy as np
import matplotlib.pyplot as plt
import warnings 

# Déclarations de fonctions
    
def fonction_mathématique_corps_noir():
    
    def valeur_luminance_corps_noir (lambda_, T = T_0):
        """ Renvoie la fonction de luminance d'un corps noir de température T 
        et en fonction de la longueur d'onde lambda_ """
        term_1 = np.pi * (C_1 * (10 **6) ** 4) / lambda_ ** 5
        term_2 = np.exp(C_2 * (10 ** 6) / (T * lambda_)) - 1
        return term_1 / term_2
    
    return valeur_luminance_corps_noir

def produit_de_fonctions(fonction1, fonction2):
    
    def fonction_produit(x):
        return fonction1(x) * fonction2(x)
    
    return fonction_produit

def tableau_valeurs_fonction (fonction, x_min, x_max, nb_points):
    """ Renvoie un tableau de valeurs (discrétisation) d'une fonction sur un 
    intervalle [x_min ; x_max] avec un pas de delta """
    tab_values = []
    delta = (x_max - x_min) / nb_points
    x = x_min
    for _ in range(nb_points):
        x += delta
        tab_values.append(fonction(x))
    return tab_values

def intégrale(fonction, borne_inf, borne_sup): 
    return quad(fonction, borne_inf, borne_sup, limit = 10 ** 7, full_output = 1)

def spectre_transmission_CO2 (données_abscisses, données_ordonnées):
    """ Représentation graphique spectre transmission CO2 en fonction de la longueur d'onde """
    plt.plot(données_abscisses, données_ordonnées)
    plt.title("Taux de transmission du CO2 en fonction de la longueur d'onde")
    plt.xlabel("Longueur d'onde (en m)")
    plt.ylabel("Taux de transmission du CO2 (en %)")
    plt.grid(True)
    plt.show()
    return None

def spectre_luminance_corps_noir (données_abscisses, données_ordonnées):
    """ Représentation graphique luminance corps noir en fonction de la longuer d'onde """
    plt.scatter(données_abscisses, données_ordonnées, marker = '.')
    plt.title("")
    plt.xlabel("Longueur d'onde (en m)")
    plt.ylabel('Luminance spectrale (en kg.m^-1.s^-3)')
    # 'Exitance totale corps noir (kg.s^-3.K^-4)'
    plt.grid(True)
    plt.show()
    return None

def test_fonction_mathématique(abcisses, fonction):
    x = []
    y = []
    for valeur in abcisses:
        décalage = random()
        image_valeur = évaluer_fonction_interpolation(fonction, valeur)
        x.append(valeur + décalage)
        y.append(image_valeur)
    plt.scatter(x, y)
    plt.show()
    return None

def fonction_mathématique_interpolation(longueur_onde, taux_CO2):
    
    def valeur_interpolation(longueur_onde, taux_CO2):
        return interp1d(longueur_onde, taux_CO2, kind = 'linear', fill_value = 'extrapolate')
    
    return valeur_interpolation(longueur_onde, taux_CO2)

def évaluer_fonction_interpolation(fonction, valeur):
    return (fonction.__call__(valeur)).tolist()

def affichage_physique(paramètre, M_0, T):
    """ Affichage de résultat avec valeur, unité et incertitude associées """
    print(f"Exitance {paramètre}: M_0({T} K) = {M_0[0]} ± {M_0[1]} W.m^-2 = kg.s^-3.K^-4")
        
# Programme principal

T = T_0

system('cls' if name == 'nt' else 'clear')

""" Chargement données selon le modèle choisi """
longueur_onde, taux_CO2 = chargement_données()

""" Transformation vers fonctions mathématiques continues """
fonction_mathématique_taux_CO2_longueur_onde = fonction_mathématique_interpolation(longueur_onde, taux_CO2)

""" Calculs différents flux selon les cas """

warnings.filterwarnings('ignore') # Supprimer les warnings

M_0_théorique_sans_CO2_corps_noir = C_S * T ** 4
print(f'Exitance théorique sans CO2 à {T} K = {M_0_théorique_sans_CO2_corps_noir} W.m^-2 = kg.s^-3.K^-4')

intégrande_sans_taux_CO2 = fonction_mathématique_corps_noir()
M_0_sans_CO2 = intégrale(intégrande_sans_taux_CO2, 0, np.inf)
affichage_physique('sans CO2', M_0_sans_CO2, T)

intégrande_avec_taux_C02 = produit_de_fonctions(intégrande_sans_taux_CO2, fonction_mathématique_taux_CO2_longueur_onde)
M_0_avec_CO2 = intégrale(intégrande_avec_taux_C02, 0, np.inf)
affichage_physique('avec CO2', M_0_avec_CO2, T)

""" Affichage spectre tranmission CO2 """
# spectre_transmission_CO2(longueur_onde, taux_CO2)

""" Affichage spectre luminance spectrale corps noir """
lumiance_corps_noir_tableau = tableau_valeurs_fonction(intégrande_sans_taux_CO2, 
                                                       min(longueur_onde), max(longueur_onde),
                                                       len(longueur_onde))
# spectre_luminance_corps_noir(longueur_onde, lumiance_corps_noir_tableau)