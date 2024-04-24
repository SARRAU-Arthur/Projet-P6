from hapi import *

reportMissingImports = False

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
        taux_CO2.append(np.mean(data[i,1:5]) * 100)
    return 1 / nombre_onde, taux_CO2

def chargement_données_HITRAN_complet():
    """ Chargement données de la base donnée depuis le site HITRAN.org
    et redistribution dans un fichier csv en deux colonnes:
    taux transmission CO2 (en %) et nombre d'onde (en m) """
    db_begin('data') # Chargement données depuis site
    fetch('CO2', 2, 1, 500, 2100) # Accès aux données: (numéro molécule CO2 = 2) entre 500 et 2100 cm^-1
    nombre_onde, coef = absorptionCoefficient_Lorentz(SourceTables = 'CO2', 
                                         Diluent = {'air': 1.0},
                                         HITRAN_units = False)
    nombre_onde, transmittance = transmittanceSpectrum(nombre_onde, coef)
    nom_fichier = chemin_acces('Bases de données','CO2 Absorption HITRAN','csv')
    with open(nom_fichier, mode = 'w', newline = '') as fichier_csv:
        for lignes in range(0,len(nombre_onde)):
            fichier_csv.write(str(nombre_onde[lignes]) + ' ; ' + str(transmittance[lignes]) + '\n')
    return None

def chargement_données_HITRAN():
    """ Chargement données de la base donnée HITRAN en deux tableaux: 
    taux transmission CO2 (en %) et longueur d'onde (en m) """
    data = np.loadtxt(chemin_acces('Bases de données','CO2 Absorption HITRAN','csv'), delimiter = ';')
    taux_CO2 = []
    nombre_onde = 1E-2 * data[:,0]
    for i in range(0,np.size(data[:,0])):
        taux_CO2.append(data[i,1] * 100) 
    return 1 / nombre_onde, taux_CO2
