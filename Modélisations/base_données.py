from hapi import *

reportMissingImports = False

from constantes import *

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
    nombre_onde, coef = absorptionCoefficient_Lorentz(SourceTables = 'CO2', 
                                         Diluent = {'air': 1.0},
                                         HITRAN_units = False)
    nombre_onde, transmittance = transmittanceSpectrum(nombre_onde, coef)
    nom_fichier = chemin_acces('Bases de données','CO2 Absorption z constant HITRAN','csv')
    with open(nom_fichier, mode = 'w', newline = '') as fichier_csv:
        for lignes in range(0,len(nombre_onde)):
            fichier_csv.write(str(nombre_onde[lignes]) + ' ; ' + str(transmittance[lignes]) + '\n')
    return None

def fonction_mathématique_température_altitude():
    
    def température_altitude(z):
        z_trop = 1.0E4
        z_strat1 = 1.9E4
        z_strat2 = 3.2E4
        z_meso = 4.7E4
        
        def coefficients(z):
            a = (-1.4452E2, 2.19E2, 1.300E3, 3.41E2)
            b = (4.173913E4, 0, -2.65E4, -4.6068E4)
            return (z,
                    [z < z_trop,
                    (z >= z_strat1) & (z < z_strat2),
                    (z >= z_strat2) & (z < z_meso),
                    z >= z_meso],
                    [(a[0], b[0]),
                     (a[1], b[1]),
                     (a[2], b[2]),
                     (a[3], b[3])])
         
        return coefficients, \
                (z,
                [z < z_trop,
                    (z >= z_strat1) & (z < z_strat2),
                    (z >= z_strat2) & (z < z_meso),
                    z >= z_meso],
                [lambda z: coefficients[0][0] * z + coefficients[0][1],
                    lambda z: coefficients[1][0] * z + coefficients[1][1],
                    lambda z: coefficients[2][0] * z + coefficients[2][1],
                    lambda z: coefficients[3][0] * z + coefficients[3][1]])
    
    return température_altitude
    
def fonction_mathématique_quantité_matière_altitude():
    """ Renvoie la fonction de luminance d'un corps noir de température T 
        et en fonction de la longueur d'onde lambda_ """
    
    def quantité_matière_altitude (z):
        """ Renvoie l'image d'une valeur lamba_ donnée à travers la fonction de luminance 
        d'un corps noir de température T """
        T =  fonction_mathématique_température_altitude[2]
        term_1 = P_0 * (1 - coefficients(z) * z / T_T)
        term_2 = np.exp(C_2 * (10 ** 6) / (T * lambda_)) - 1
        return term_1 / term_2
    
    return quantité_matière_altitude

def chargement_données_HITRAN_complet_fonction_z():
    """ Chargement données de la base donnée depuis le site HITRAN.org
    et redistribution dans un fichier csv en deux colonnes:
    taux transmission CO2 (pas en %) et nombre d'onde (en m) """
    db_begin('data') # Chargement données depuis site
    fetch('CO2', 2, 1, 500, 2100) # Accès aux données: (numéro molécule CO2 = 2) entre 500 et 2100 cm^-1
    nombre_onde, coef = absorptionCoefficient_Lorentz(SourceTables = 'CO2', 
                                         Diluent = {'air': 1.0},
                                         HITRAN_units = False)
    nombre_onde, transmittance = transmittanceSpectrum(nombre_onde, coef)
    nom_fichier = chemin_acces('Bases de données','CO2 Absorption fonction z HITRAN','csv')
    with open(nom_fichier, mode = 'w', newline = '') as fichier_csv:
        for lignes in range(0,len(nombre_onde)):
            fichier_csv.write(str(nombre_onde[lignes]) + ' ; ' + str(transmittance[lignes]) + '\n')
    return None

def chargement_données_HITRAN():
    """ Chargement données de la base donnée HITRAN en deux tableaux: 
    taux transmission CO2 (pas en %) et longueur d'onde (en m) """
    data = np.loadtxt(chemin_acces('Bases de données','CO2 Absorption HITRAN','csv'), delimiter = ';')
    taux_CO2 = []
    nombre_onde = 1E-2 * data[:,0]
    for i in range(0,np.size(data[:,0])):
        taux_CO2.append(data[i,1]) 
    return 1 / nombre_onde, taux_CO2

def chargement_données():
    return chargement_données_HITRAN()
