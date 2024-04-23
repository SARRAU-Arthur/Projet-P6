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

# Déclarations de fonctions

def intégrande_luminance_corps_noir (lambda_, T):
    """ Renvoie la fonction de luminance d'un corps noir de température T 
    et en fonction de la longueur d'onde lambda_ """
    term_1 = np.pi * (C_1 * (10 **6) ** 4) / lambda_ ** 5
    term_2 = np.exp(C_2 * (10 ** 6) / (T * lambda_)) - 1
    return term_1 / term_2

def tableau_valeurs_fonction (fonction, x_min, x_max, delta):
    """ Renvoie un tableau de valeurs (discrétisation) d'une fonction sur un 
    intervalle [x_min ; x_max] avec un pas de delta """
    tab_values = [None] * int(np.abs((x_max - x_min)) / delta)
    x = x_min
    for i in range (0,(np.size(tab_values))):
        x += delta
        tab_values[i] = fonction(x)
    return tab_values

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
    plt.plot(données_abscisses,données_ordonnées)
    plt.title("Taux de transmission du CO2 en fonction de la longueur d'onde")
    plt.xlabel("Longueur d'onde (en m)")
    plt.ylabel("Taux de transmission du CO2 (en %)")
    plt.grid(True)
    plt.show()
    return ()

def spectre_luminance_corps_noir (données_abscisses, données_ordonnées):
    """ Représentation graphique luminance corps noir en fonction de la longuer d'onde """
    plt.plot(données_abscisses,données_ordonnées)
    plt.title("")
    plt.xlabel("Longueur d'onde (en m)")
    plt.ylabel('Luminance spectrale (en kg.m^-1.s^-3)')
    # 'Exitance totale corps noir (kg.s^-3.K^-4)'
    plt.grid(True)
    plt.show()
    return ()

def test_fonction_mathématique(abcisses, fonction):
    points = []
    for valeur in abcisses:
        décalage = random()
        points.append([valeur + décalage, fonction(valeur)])
    print(type(fonction(valeur)))
    plt.scatter(points[:][0], points[:][1])
    plt.show()
    return None

def fonction_interpolation(tableau_abcisses, tableau_ordonnées):
    return interp1d(tableau_abcisses, tableau_ordonnées)
        
# Programme principal

# system('cls' if name == 'nt' else 'clear')

""" Chargement données et grandeurs """
longueur_onde, taux_CO2 = chargement_données_HITRAN()


""" Transformation vers fonctions mathématiques continues """
fonction_taux_CO2_longueur_onde = fonction_interpolation(longueur_onde, taux_CO2)
print(type(fonction_taux_CO2_longueur_onde))
test_fonction_mathématique(longueur_onde, fonction_taux_CO2_longueur_onde)
aire_scipy = quad(fonction_taux_CO2_longueur_onde, min(longueur_onde), max(longueur_onde), 
                  limit = 1000, full_output = 0)[0]
print(aire_scipy)

""" Calcul flux """
luminance_corps_noir_tableau = intégrande_luminance_corps_noir(longueur_onde, T)
intégrande = taux_CO2 * luminance_corps_noir_tableau
exitance_classique = intégrale_coefficient_ensemble_discret(0, np.inf, longueur_onde, luminance_corps_noir_tableau)
exitance_taux_CO2 = intégrale_coefficient_ensemble_discret(0, np.inf, longueur_onde, intégrande)
print('Flux classique = ', exitance_classique)
print('Flux théorique = ', C_S * T ** 4)
print('Flux avec CO2 = ', exitance_taux_CO2)

""" Calcul intégral de la valeur avec incertitude de l'exitance totale du corps noir """
# M_0 =
# print(f"M_0 = {M_0[0]} ± {M_0[1]} kg.s^-3.K^-4")

""" Affichage spectre CO2 """
# spectre_transmission_CO2(longueur_onde, taux_CO2)

""" Affichage de la luminance spectrale en fonction de la longueur d'onde """
# spectre_luminance_corps_noir(longueur_onde, luminance_corps_noir_tableau)