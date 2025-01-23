[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/oOQR1xPR)
# Librairie _Python filters_

## Membres du groupe :
- Maxime Magnenat
- Francisco Oliveira Barbosa
- Maxime Otero
- Sébastien Pfister
- David Vuillemier

## Description du projet 
_Python filters_ est un projet de 5 étudiants de l'école d'ingénieur d'Yverdon-les-Bainse (HEIG-VD).

Le projet est une librairie open-source d'aide à la conception de filtres électroniques par le calcul de valeurs de composants, nous abordons en détail les fonctionnalités plus bas.
Les filtres pris en charge par la librairie sont les suivants:
- Actifs (ordre 1 à 10)
  - Tchebychev
  - Bessel
  - Butterworth
- Passifs (ordre 1 et 2)

Pour chacun des types de filtres sus-mentionnés, la librairie prend en charge les passe haut et passe bas.
Les filtres actifs sont modélisés par des couplages de cellules de Sallen & Key [lien Wikipedia sur les cellules de Sallen & Key](https://en.wikipedia.org/wiki/Sallen%E2%80%93Key_topology)
Ces derniers sont dimmensionnés à l'aide du tableau de valeurs normalisées (wk/wr et facteur de qualité)

### Fonctionnalités
La conception de filtres requièrent de longs calculs chronophages et qui s'oublient facilement, dans ce sens, l'utilisateur de la librairie indique les paramètres connus et désirés, à savoir:
- Nature des composants (actif, passif)
- Type de filtre (Tchebychev, Bessel, Butterworth)
- Ordre
- Fréquence de coupure
- Composants connus (résistances ou condensateurs)

et en retour le code va calculer les valeurs de composants qui lui permettent d'avoir le filtre désiré.
Condensateurs en entrée  --->  valeurs des résistances en sortie
Résistances en entrée    --->  valeurs des condensateurs en sortie

### Structure
Notre bibliothèque se découpera sous la forme de fonctions.
La première qui calcule les composants prend en paramètre l'ordre du filtre, les pulsations voulues, ainsi qu'au moin une liste de composants que ce soit des résistances ou alors des condensateurs.
Pour la deuxième qui prend en paramètre l'ordre du filtre, les pulsations voulues ainsi que les composants qui composent ce filtre.
Il sera possible d'afficher deux graphiques qui le diagramme de Bode en amplitude ainsi qu'en phase

## Caractéristiques de la librairie
- Tout le code est en python
- Nous utilisons poetry afin de faciliter la gestion des dépendances
- Pytest est utilisé afin de s'assurer en tout temps du bon fonctionnement de la librairie et de la justesse des calculs, certains tests sont basés sur le support de cours _Electronique analogique_
- La qualité du code est contrôlée avec Black and Ruff, pour suivre PEP8
- Le fichier "requierements.txt" liste les librairies externes utilisées par _Python filters_


## Déploiement

La librairie a pour ambition de perdurer au délà du projet de classe, afin d'être open source et améliorable, en ajoutant d'autres fonctionnalités d'aide à la conception

## Exemple d'utilisation

```python
from filters.snk.bessel import lowpass, highpass

highpass_test = highpass()
r_vals = [1000,10000,1000,10000]
highpass_test.graphs(order =1, cutoff_freq=1000,r_vals=r_vals)

tf,stages = highpass_test.components(order=4, cutoff_freq=1000, r_vals=r_vals)
print("Fonction de transfert combinée :", tf)
for i, stage in enumerate(stages):
    print(f"Stage {i+1}:")
    print("Fonction de transfert :", stage['tf'])
    print("Paramètres :", stage['params'])
```

### Détail du code

Import de la classe qui nous intéresse :
```python
from filters.snk.bessel import lowpass, highpass
```
On indique les valeurs de résistances :
```python
r_vals = [1000,10000,1000,10000]
```
Le filtre voulu est un ordre 1 avec une fréquence de coupure de 1000 Hz
```python
highpass_test.graphs(order =1, cutoff_freq=1000,r_vals=r_vals)
```
Dans le code ci-dessous, on indique des valeurs de résistances

## Graphique
Voici un exemple de graphique obtenu à l'issue de l'exécution du code :
<img width="601" alt="image" src="https://github.com/user-attachments/assets/f150da63-2435-482a-9850-15398067f003" />


## Commandes importantes

### Exécution des tests

```python
pytest -v -s
```

### Préparation de la distribution

```python
poetry build
```

### Installation

```python
pip install filters-0.1.0-py3-none-any.whl
```

## Références
Support de cours _Electronique analogique_ du département TIN, rédigé par Monsieur Blaise Grandjean
