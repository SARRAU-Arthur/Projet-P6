# Importation de bibliothèques

import matplotlib.pyplot as plt # Bibliothèque graphique pour tracer des courbes fonctionelles
import numpy as np # Bibliothèque maths avec l'ensemble des opérateurs, outils statistiques et fonctions de référence
from scipy.interpolate import interp1d

from base_données import chargement_données
from modèle import évaluer_fonction_interpolation

# Constantes

h = 6.62607015e-34 # Constante de Plank (J.s)
c = 2.998E8 # Célérité de la lumière (m.s^-1)
kB = 1.380649E-23 # Constante de Boltzmann (J.K^-1)
R = 8.31446262 # Constante universelle des gaz parfaits (J.mol^−1.K^−1)
g = 9.80665 # Constante gravitationelle terrestre (m.s^-2)
T_0 = 2.8815E2 # Température 15 degrés Celcius (K)

# Fonctions et implémentation

# Loi de Plank, calcul de la luminance énergétique spectrale à une longueur d'onde et une température données (W.s.m^−2.sr^−1)
def planck_function(lambda_wavelength, T):
    term1 = (2 * h * c**2) / lambda_wavelength**5
    term2 = np.exp((h * c) / (lambda_wavelength * kB * T)) - 1
    return term1 / term2

# Formule du modèle de nivellement barométrique, calcul de la pression à une altitude donnée (T et M constants)
def pressure(z):
    P_0 = 1.01325E4 # Pression au niveau de la mer (Pa)
    T = T_0 # Température caractéristique (K)
    M = 2.8966E-3 # Masse molaire de l'air (g.mol^-1)
    H = (R * T) / (M * g) # Altitude caractéristique notée delta_h (m)
    return P_0 * np.exp(-z / H)

# Création d'un tableau rempli de T_0, mais la fonction n'est réutilisée nulle part dans le code... à approfondir 
def temperature_uniform(z):
    print(np.ones_like(z))
    return T_0 * np.ones_like(z)

# Modèle simple: Calcul de la fonction de température (K) selon l'altitude (m). Découpage en 2 domaines
def temperature_simple(z):
    z_trop = 1.1E4  # Hauteur de la troposphère (m)
    Gamma = -6.5E-3 # Gradient de température par unité de longueur (K.m^-1)
    T_trop = T_0 + Gamma * z_trop # Équation modèle affine (z<z_trop)
    return np.piecewise(z,
                        [z < z_trop, 
                         z >= z_trop],
                        [lambda z: T_0 + Gamma * z,
                         lambda z: T_trop])
                        # Fonctions lambda et piecewise permettent de définir une fonction par morceaux selon des conditions ajustables

# Modèle complexe: Calcul de la fonction de température (K) selon l'altitude (m). Découpage en 7 domaines
def temperature_US1976(z):
    z_km = z/1E3  # Altitude convertie (km)
    # Troposphère (0 à 11 km)
    T0 = 288.15
    z_trop = 11
    # Tropopause (11 à 20 km)
    T_tropopause = 216.65
    z_tropopause = 20
    # Stratosphère 1 (20 à 32 km)
    T_strat1 = T_tropopause
    z_strat1 = 32
    # Stratosphère 2 (32 à 47 km)
    T_strat2 = 228.65
    z_strat2 = 47
    # Stratopause (47 à 51 km)
    T_stratopause = 270.65
    z_stratopause = 51
    # Mesosphère 1 (51 à 71 km)
    T_meso1 = T_stratopause
    z_meso1 = 71
    # Mesosphère 2 (71 à ...)
    T_meso2 = 214.65
    return (z_km,
                        [z_km < z_trop,
                         (z_km >= z_trop) & (z_km < z_tropopause),
                         (z_km >= z_tropopause) & (z_km < z_strat1),
                         (z_km >= z_strat1) & (z_km < z_strat2),
                         (z_km >= z_strat2) & (z_km < z_stratopause),
                         (z_km >= z_stratopause) & (z_km < z_meso1),
                         z_km >= z_meso1],
                        [lambda z: T0 - 6.5 * z,
                         lambda z: T_tropopause,
                         lambda z: T_strat1 + 1 * (z - z_tropopause),
                         lambda z: T_strat2 + 2.8 * (z - z_strat1),
                         lambda z: T_stratopause,
                         lambda z: T_meso1 - 2.8 * (z - z_stratopause),
                         lambda z: T_meso2 - 2 * (z - z_meso1)]) 
                        # Modèles d'équation exclusivement affines, avec des coefficients variables selon la couhe atmosphérique
                        
# Choix du modèle de température
def temperature(z):
    return temperature_US1976(z)

# Calcul de la densité de l'air (kg.m^-3) à une altitude donnée (m)
def air_number_density(z):
    return pressure(z) / (kB * temperature(z))

# Calcul de section efficace du C02 (m^2) à la longueur d'onde (m) donnée
def cross_section_CO2(wavelength):
    LAMBDA_0 = 1.50E-5 # Bande centrale (m)
    exponent = -22.5 - 24 * np.abs((wavelength - LAMBDA_0) / LAMBDA_0)
    sigma = 10 ** exponent
    return sigma

def fonction_mathématique_interpolation (longueur_onde, k_abs):
    """ Retourne l'absorbance et la transmittance du CO2 en fonction de 
    la longueur d'onde (en m) sous forme d'un objet de classe <function> 
    qui est mathématiquement continue """
    
    def valeur_interpolation (longueur_onde, k_abs):
        """ Retourne la transmittance du CO2 à une longueur d'onde donnée, 
        sous la forme d'un objet de classe <interp1D> """
        return interp1d(longueur_onde, k_abs, kind = 'linear', 
                        bounds_error = False, fill_value = (100, 100))
    
    k_abs_fonction = lambda x: évaluer_fonction_interpolation(valeur_interpolation(longueur_onde, k_abs), x) \
                    if min(longueur_onde) <= x <= max(longueur_onde) \
                    else 1
    return k_abs_fonction

# RADIATIVE TRANSFER SIMULATION, All wavelengths are treated in parallel using vectorization
def simulate_radiative_transfer(CO2_fraction, z_max = 8E4, delta_z = 1E1, lambda_min = 1E-7, lambda_max = 1E-4, delta_lambda = 1E-8):
    z_range = np.arange(0, z_max, delta_z) # Tableau d'altitudes de 0 à z_max avec un écart constant de delta_z entre deux altitude successives
    lambda_range = np.arange(lambda_min, lambda_max, delta_lambda) # Idem pour la longueur d'onde, "arrange" crée des plages de données homogénement réparties
    upward_flux = np.zeros((len(z_range), len(lambda_range))) # Flux vers le haut
    optical_thickness = np.zeros((len(z_range), len(lambda_range))) # Initialisation tableaux à zéro
    earth_flux = np.pi * planck_function(lambda_range, temperature(0)) * delta_lambda # Calcul flux émis par la Terre selon la loi de Plank, P=sigma*T^4
    print(f"Total earth surface flux in wavelength range: {earth_flux.sum():.2f} W/m^2") # Affichage de la valeur trouvée
    flux_in = earth_flux
    for i, z in enumerate(z_range): # Instruction "enumerate" permet de boucler sur i et z simultanément, et évite un double boucle for imbriquée
        n_CO2 = air_number_density(z) * CO2_fraction # Densité CO2 en fonction de l'altitude dans l'atmosphère
        # kappa = cross_section_CO2(lambda_range) * n_CO2 # Calcul du coefficient d'absorption du C02 en fonction de sa section efficcace dépendant de la longueur d'onde
        fonction_k_abs = fonction_mathématique_interpolation(chargement_données())
        print('OK')
        kappa = fonction_k_abs(lambda_range) * n_CO2
        # Compute fluxes within the layer
        optical_thickness[i,:] = kappa * delta_z
        absorbed_flux = np.minimum(kappa * delta_z * flux_in , flux_in)
        emitted_flux = optical_thickness[i,:] * np.pi * planck_function(lambda_range, temperature(z)) * delta_lambda
        upward_flux[i, :] = flux_in - absorbed_flux + emitted_flux
        # The flux leaving the layer becomes the flux entering the next layer
        flux_in = upward_flux[i, :]
    print(f"Total outgoing flux at the top of the atmosphère: {upward_flux[-1,:].sum():.2f} W/m^2")
    return lambda_range, z_range, upward_flux, optical_thickness

# Programme principal

CO2_fraction = 2.8E-4 # Pourcentage de CO2 atmosphérique
lambda_range, z_range, upward_flux, optical_thickness = simulate_radiative_transfer(CO2_fraction) # Calcul paramètres correspondants
CO2_fraction *= 2 # Multiplication de la teneur actuelle en CO2
lambda_range, z_range,  upward_flux2, optical_thickness2 = simulate_radiative_transfer(CO2_fraction) # Idem

# Affichage sous forme de courbe 

range_factor=1E6
delta_lambda = lambda_range[1] - lambda_range[0]

plt.figure(figsize=(14, 9)) # Taille de la case du graphique, peu d'importance
plt.plot(range_factor * lambda_range, np.pi * planck_function(lambda_range, temperature(0))/range_factor,'--k') # Courbe de spectre de corps noir: T = Température surface Terre
plt.plot(range_factor * lambda_range, np.pi * planck_function(lambda_range, 216)/range_factor,'--k') # Courbe de spectre de corps noir: T = 220 K
plt.plot(range_factor * lambda_range, upward_flux[-1, :]/delta_lambda/range_factor,'-g')
plt.plot(range_factor * lambda_range, upward_flux2[-1, :]/delta_lambda/range_factor,'-r')
plt.fill_between(range_factor * lambda_range, upward_flux[-1, :]/delta_lambda/range_factor, upward_flux2[-1, :]/delta_lambda/range_factor, color='yellow', alpha=0.9)
plt.xlabel("Longueur d'onde (μm)") # Légende x
plt.ylabel("Luminance spectrale (W/m²/μm/sr)") # Légende y
plt.xlim(0,50) # Intervalle x
plt.ylim(0,30) # Intervalle y
plt.grid(True) # Affichage grille
plt.show()

plt.plot(temperature(np.arange(0, 100, 1)))
plt.xlabel("Altitude (km)")
plt.ylabel("Température (K)")
plt.xlim(0,100)
plt.ylim(-50,20)
plt.grid(True)
plt.show()