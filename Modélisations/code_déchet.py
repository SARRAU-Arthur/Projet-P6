def fonction_taux_transmission_CO2(longueur_onde, taux_CO2):
    fonction_par_morceaux = []
    for i in range(len(longueur_onde)-1):
        fonction_par_morceaux.append(lambda x, 
                                     borne_inf = longueur_onde[i + 1], 
                                     borne_sup = longueur_onde[i], 
                                     y = taux_CO2[i]:                                         
                                     0 if borne_inf <= x <= borne_sup else y)
        
    def fonction(x):
        for f in fonction_par_morceaux:
            if f(x) != 0:
                return f(x)
            else:
                pass
    
    return fonction

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