from hapi import *

reportMissingImports = False

from constantes import *
from scipy.integrate import quad
import matplotlib.pyplot as plt

import numpy as np

def chemin_acces(langue, lettre, extension):
    """ Raccourci pour aller chercher le fichier au format choisi 
    dans un dossier potentiellement différent de l'actuel """
    return f'../{langue}/{lettre}.{extension}'

def chargement_données_NIST():
    """ Chargement données de la base donnée NIST en deux tableaux: 
    taux transmission CO2 (en %) et longueur d'onde (en m) """
    data = np.loadtxt(chemin_acces('Bases de données','CO2 Absorption NIST','csv'), delimiter = ';')
    taux_CO2 = []
    nombre_onde = 1E-2 * data[:-1,0] # Nan dernière ligne, on exclu pour éviter erreurs à l'exécution
    for i in range(0,np.size(data[:,0]) - 1):
        taux_CO2.append(np.mean(data[i,1:5]))
    return 1 / nombre_onde, taux_CO2

def chargement_données_HITRAN_complet_z_constant():
    """ Chargement données de la base donnée depuis le site HITRAN.org
    et redistribution dans un fichier csv en deux colonnes:
    taux transmission CO2 (pas en %) et nombre d'onde (en m) """
    db_begin('data') # Chargement données depuis site
    fetch('CO2', 2, 1, 500, 2100) # Accès aux données: (numéro molécule CO2 = 2) entre 500 et 2100 cm^-1
    nombre_onde, k_abs = absorptionCoefficient_Lorentz(SourceTables = 'CO2', 
                                         Diluent = {'air': 1.0},
                                         HITRAN_units = False)
    nombre_onde, transmittance = transmittanceSpectrum(nombre_onde, k_abs)
    nom_fichier = chemin_acces('Bases de données','CO2 Absorption z constant HITRAN','csv')
    with open(nom_fichier, mode = 'w', newline = '') as fichier_csv:
        for lignes in range(0,len(nombre_onde)):
            fichier_csv.write(str(nombre_onde[lignes]) + ' ; ' + str(transmittance[lignes]) + '\n')
    return None

def fonction_mathématiques_coefficients():
    
    def coefficient_a(z):
        a = [-6.5E-3, 0, 1E-3, 2.8E-3]
        return np.piecewise(z,
                [z < z_trop,
                    (z_trop <= z) & (z < z_strat1),
                    (z_strat1 <= z) & (z < z_strat2),
                    (z_strat2 <= z) & (z < z_meso),],
                [a[0],
                    a[1],
                    a[2],
                    a[3]])
    
    def coefficient_b(z):
        b = [288.15, 216.5, 196.5, 138.9]
        return np.piecewise(z,
                [z < z_trop,
                    (z_trop <= z) & (z < z_strat1),
                    (z_strat1 <= z) & (z < z_strat2),
                    (z_strat2 <= z) & (z < z_meso),],
                [b[0],
                    b[1],
                    b[2],
                    b[3]])
            
    return coefficient_a, coefficient_b

def fonction_mathématique_température_altitude():
    
    def température_altitude(z):
        a, b = fonction_mathématiques_coefficients()
        return a(z) * z + b(z)
    
    return température_altitude
    
def fonction_mathématique_quantité_matière_altitude():
    """ Renvoie la fonction de luminance d'un corps noir de température T 
        et en fonction de la longueur d'onde lambda_ """
    
    def quantité_matière_altitude(z):
        """ Renvoie l'image d'une valeur lamba_ donnée à travers la fonction de luminance 
        d'un corps noir de température T """
        T = fonction_mathématique_température_altitude()
        a, b = fonction_mathématiques_coefficients()
        term_1 = lambda z: P_0 * np.exp((-g * M) / (R * T(z))) if (z_trop <= z) & (z < z_strat1) \
                              else P_0 * (1 - a(z) * z / b(z)) ** (M * g) / (R * a(z))
        term_2 = k_B * T(z)
        return term_1(z) / term_2
    
    return quantité_matière_altitude

def graphique_quantité_matière_fonction_altitude():
    x = np.linspace(1, 1E4)
    y = fonction_mathématique_quantité_matière_altitude()
    plt.plot(x, y(x))
    plt.show()
    return None

def chargement_données_HITRAN_complet_fonction_z():
    """ Chargement données de la base donnée depuis le site HITRAN.org
    et redistribution dans un fichier csv en deux colonnes:
    taux transmission CO2 (pas en %) et nombre d'onde (en m) """
    db_begin('data') # Chargement données depuis site
    fetch('CO2', 2, 1, 500, 2100) # Accès aux données: (numéro molécule CO2 = 2) entre 500 et 2100 cm^-1
    nombre_onde, k_abs = absorptionCoefficient_Lorentz(SourceTables = 'CO2', 
                                         Diluent = {'air': 1.0},
                                         HITRAN_units = False)
    transmittance = k_abs * quad(fonction_mathématique_quantité_matière_altitude(), 0, h_max, \
                                 limit = 10 ** 7, full_output = 1)
    nom_fichier = chemin_acces('Bases de données', 'CO2 Absorption fonction z HITRAN', 'csv')
    with open(nom_fichier, mode = 'w', newline = '') as fichier_csv:
        for lignes in range(0,len(nombre_onde)):
            fichier_csv.write(str(nombre_onde[lignes]) + ' ; ' + str(transmittance[lignes]) + '\n')
    return None

def chargement_données_HITRAN(nom_fichier):
    """ Chargement données de la base donnée HITRAN en un tableau à 2 colonnes: 
    taux transmission CO2 (pas en %) et longueur d'onde (en m). Préciser dans <nom_fichier>
    si l'on souhaite utiliser les données avec z constant ou en fonction de z """
    data = np.loadtxt(chemin_acces('Bases de données', nom_fichier, 'csv'), delimiter = ';')
    taux_CO2 = []
    nombre_onde = 1E-2 * data[:,0]
    for i in range(0,np.size(data[:,0])):
        taux_CO2.append(data[i,1]) 
    return 1 / nombre_onde, taux_CO2

def chargement_données():
    return chargement_données_HITRAN('CO2 Absorption fonction z HITRAN')

x = np.linspace(0, h_max -1 )
T = fonction_mathématique_température_altitude()
plt.plot(T(x), x)
plt.show()

print(quad(fonction_mathématique_quantité_matière_altitude(), 0, z_trop + 1, limit = 10 ** 7, full_output = 1))
# chargement_données_HITRAN_complet_fonction_z()
