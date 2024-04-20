# Import bibliothèques

from scipy.integrate import quad
from scipy.interpolate import interp1d
from random import random
from os import system, name
import numpy as np
import matplotlib.pyplot as plt

# Définition constantes d'après les constantes fondamentales en unité SI

h = 6.62607015e-34 # Constante de Plank (kg.m^2.s^-1)
c_0 = 2.998E8 # Célérité de la lumière (m.s^-1)
k_B = 1.380649E-23 # Constante de Boltzmann (kg.m^2.s^-2.K^-1)
R = 8.31446262 # Constante universelle des gaz parfaits (kg.m^2.s^-2.mol^−1.K^−1)
g = 9.80665 # Constante gravitationelle terrestre (m.s^-2)
n = 1 # Indice de réfraction
C_1 = 2 * h * (c_0 / n) ** 2 # Constante de Plank 1 (kg.s^−3.m^-4)
C_2 = (h * c_0) / (n * k_B) # Constante de Plank 2 (m.K)
C_S = 5.67E-8 # Constante de Stephan (kg.s^-3.K^-4)
T_S = 5.772E3 # Température soleil (K)
T_0 = 2.8815E2 # Température 15 degrés Celcius (K)

T = T_0

# Déclarations de fonctions

def chemin_acces (langue, lettre, extension):
    """ Raccourci pour aller chercher le fichier au format choisi 
    dans un dossier potentiellement différent de l'actuel """
    return f'../{langue}/{lettre}.{extension}'

def intégrande_luminance_corps_noir (lambda_, T):
    """ Renvoie la fonction de luminance d'un corps noir de température T 
    et en fonction de la longueur d'onde lambda_ """
    term_1 = np.pi * (C_1 * (10 **6) ** 4) / lambda_ ** 5
    term_2 = np.exp(C_2 * (10 ** 6) / (T * lambda_)) - 1
    return term_1 / term_2

def intégrande_luminance_corps_noir_discrétisation (longueur_onde, T):
    """ Renvoie la fonction de luminance d'un corps noir de température T 
    et en fonction de la longueur d'onde x """
    term_1 = np.pi * (C_1 * (10 **6) ** 4) / longueur_onde ** 5
    term_2 = np.exp(C_2 * (10 ** 6) / (T * longueur_onde)) - 1
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
    
def chargement_données(nom_fichier_avec_extension, delimiteur_csv):
    """ Chargement données de la base donnée en deux tableaux: 
    taux transmission CO2 (en %) et longueur d'onde (en m) """
    data = np.loadtxt(nom_fichier_avec_extension, delimiter = delimiteur_csv)
    taux_CO2 = []
    nombre_onde = 1E-2 * data[:-1,0] # Nan dernière ligne, on exclu pour éviter erreurs à l'exécution
    for i in range(0,np.size(data[:,0]) - 1):
        taux_CO2.append(np.mean(data[i,1:5]) * 100)
    return 1 / nombre_onde, taux_CO2

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
    # 'Exitance totale corps noir ((kg.s^-3.K^-4))'
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

system('cls' if name == 'nt' else 'clear')

""" Chargement données et grandeurs """
longueur_onde, taux_CO2 = chargement_données(chemin_acces('Bases de données','CO2 Absorption NIST','csv'), ';')
lambda_min = min(longueur_onde)
lambda_max = max(longueur_onde)

""" Transformation vers fonctions mathématiques continues """
fonction_taux_CO2_longueur_onde = fonction_interpolation(longueur_onde, taux_CO2)
print(type(fonction_taux_CO2_longueur_onde))
test_fonction_mathématique(longueur_onde, fonction_taux_CO2_longueur_onde)
aire_scipy = quad(fonction_taux_CO2_longueur_onde, min(longueur_onde), max(longueur_onde), 
                  limit = 1000, full_output = 0)[0]
print(aire_scipy)

""" Calcul flux """
luminance_corps_noir_tableau = intégrande_luminance_corps_noir_discrétisation(longueur_onde, T)
intégrande = taux_CO2 * luminance_corps_noir_tableau
exitance_classique = intégrale_coefficient_ensemble_discret(0, np.inf, longueur_onde, luminance_corps_noir_tableau)
exitance_taux_CO2 = intégrale_coefficient_ensemble_discret(0, np.inf, longueur_onde, intégrande)
print('Flux classique = ', exitance_classique)
print('Flux théorique = ', C_S * T ** 4)
print('Flux avec CO2 = ', exitance_taux_CO2)

""" Affichage spectre CO2 """
# spectre_transmission_CO2(longueur_onde, taux_CO2)

""" Calcul intégral de la valeur avec incertitude de l'exitance totale du corps noir """
# delta_adapté = (lambda_max - lambda_min) / len(longueur_onde)
# corps_noir = tableau_valeurs_fonction(intégrande_luminance_corps_noir(T), 
#                                                   lambda_min, lambda_max, 
#                                                   delta_adapté)
# M_0 = intégrale_trapèzes(intégrande_luminance_corps_noir(T), lambda_min, lambda_max)
# print(f"M_0 = {M_0[0]} ± {M_0[1]} kg.s^-3.K^-4")

""" Affichage de la luminance spectrale en fonction de la longueur d'onde """
# spectre_luminance_corps_noir(longueur_onde, luminance_corps_noir_tableau)

""" Test fonctionnement fonction taux absorption CO2 """
# fonction_taux_CO2 = fonction_taux_transmission_CO2(longueur_onde, taux_CO2)
# print(fonction_taux_CO2(lambda_min/2))

