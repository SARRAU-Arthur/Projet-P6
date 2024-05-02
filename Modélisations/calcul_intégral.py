# Bibliothèques locales

from constantes import *
from base_données import *

# Bibliothèques externes

from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.optimize import fmin
from random import random
from os import system, name

import numpy as np
import matplotlib.pyplot as plt

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

def tableau_valeurs_fonction (fonction, x_min, x_max, delta):
    """ Renvoie un tableau de valeurs (discrétisation) d'une fonction sur un 
    intervalle [x_min ; x_max] avec un pas de delta """
    tab_values = [None] * int(np.abs((x_max - x_min)) / delta)
    x = x_min
    for i in range (0,(np.size(tab_values))):
        x += delta
        tab_values[i] = fonction(x)
    return tab_values

def intégrale(fonction, borne_inf, borne_sup): 
    return quad(fonction, borne_inf, borne_sup, limit = 10 ** 7, full_output = 0)
    
def intégrale_coefficient_ensemble_discret(borne_inf, borne_sup, liste_intégrande, liste_variable_intégration):
    """ Calcul numérique de l'intégrale de borne_inf à borne_sup de la fonction en entrée 
    grâce à la méthode des trapèzes, utile pour intégrer sur une liste discrète de valeurs"""
    aire_trapèzes = liste_intégrande[0] * (liste_variable_intégration[0] - liste_variable_intégration[1]) 
    for i in range(1, len(liste_intégrande) - 1):
        aire_trapèzes += np.dot(liste_intégrande[i] - liste_intégrande[i - 1], 
                           liste_variable_intégration[i + 1] - liste_variable_intégration[i])
    luminance_corps_noir_fonction = lambda x: intégrande_luminance_corps_noir(x, T)
    aire_scipy = quad(luminance_corps_noir_fonction, borne_inf, borne_sup, full_output = 0)[0]
    return aire_trapèzes / 2 + aire_scipy

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

def fonction_mathématique_interpolation():
    
    def valeur_interpolation():
        
        longueur_onde, taux_CO2 = chargement_données()
        return interp1d(longueur_onde, taux_CO2, kind = 'linear', fill_value = 'extrapolate')
    
    return valeur_interpolation

def évaluer_fonction_interpolation(fonction, valeur):
    return (fonction.__call__(valeur)).tolist()
        
# Programme principal

T = T_0 #

system('cls' if name == 'nt' else 'clear')

""" Chargement données selon le modèle choisi """
# longueur_onde, taux_CO2 = chargement_données()

""" Transformation vers fonctions mathématiques continues """
fonction_mathématique_taux_CO2_longueur_onde = fonction_mathématique_interpolation()

""" Calcul intégral de la valeur avec incertitude de l'exitance totale du corps noir """
M_0 = intégrale(fonction_mathématique_taux_CO2_longueur_onde, 0, np.inf)
print(f"Exitance M_0 = {M_0[0]} ± {M_0[1]} kg.s^-3.K^-4")

""" Calculs flux """
intégrande_sans_taux_CO2 = fonction_mathématique_corps_noir()
intégrande_avec_taux_C02 = produit_de_fonctions(intégrande_sans_taux_CO2, 
                                  évaluer_fonction_interpolation(fonction_mathématique_taux_CO2_longueur_onde, 
                                                                 longueur_onde))
exitance_classique = intégrale(intégrande_sans_taux_CO2, 0, np.inf)
exitance_taux_CO2 = intégrale(intégrande_avec_taux_C02, 0, np.inf)
print('Flux théorique = ', C_S * T ** 4)
print('Flux sans CO2 = ', exitance_classique)
print('Flux avec CO2 = ', exitance_taux_CO2)

""" Affichage spectre CO2 """
# spectre_transmission_CO2(longueur_onde, taux_CO2)

""" Affichage de la luminance spectrale en fonction de la longueur d'onde """
# spectre_luminance_corps_noir()