[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/oOQR1xPR)
# Mini-projet PythonGE2 
## Librairie "Pythonfilters"

### Membres du groupe :
- Maxime Magnenat
- Francisco Oliveira Barbosa
- Maxime Otero
- Sébastien Pfister
- David Vuillemier

### Description du projet 
Le but est de faire une librairie open source qui permet de calculer les différents paramètres (ω, R1 ,R2 ,... ,Rx ,c1 ,c2, ..., Cx) de filtres actifs et passifs.
Les filtres actifs seront traités avec des cellules à gain fixe (cellules de Sallen & Key) https://en.wikipedia.org/wiki/Sallen%E2%80%93Key_topology et pourront être de 3 types différents selon les besoins :
- Tchebychev
- Bessel
- Butterworth

Avec évidemment pour chacun la possibilité de faire des filtres :
- Passe-bas
- Passe-haut

### Fonctionnalités
Notre bibliothèque se découpera sous la forme de fonctions.
La première qui calcule les composants prend en paramètre l'ordre du filtre, les pulsations voulues, ainsi qu'au moin une liste de composants que ce soit des résistances ou alors des condensateurs.
Pour la deuxième qui prend en paramètre l'ordre du filtre, les pulsations voulues ainsi que les composants qui composent ce filtre.
Il sera possible d'afficher deux graphiques qui le diagramme de Bode en amplitude ainsi qu'en phase

### Détails
Les valeurs de qualité sont tiré du tableau normalisée jusqu'a l'ordre 10, donc l'ordre maximal de calcul sera 10.

### Déploiement

Nos objectifs sont de commencer cette librairie qui pourrait être développée par la suite par d'autres groupes pour étendre les fonctionnalités.

Nous devrons respecter les délais imposés par notre professeur et lui fournir un travail qui correspond "au mieux" au cahier des charges.

## Calculs en détails
### ButterWorth
#### Passe Bas
Fonction de transfert est donnée par : 
H(jw) = 1 / (1+jw(R1+R2) + C1*C2*R1*R2*(jw)^2)

Quand les condensateurs sont inconnus :
C2 = 1 / ((R1+R2)*W0*Q0)
C1 = (R1+R2)*Q0 / R1*R2*W0

Quand les résistances sont inconnues :
Une contrainte est imposé qui est : C1 >= 4*Q0^2*C2
Donc : R1+R2 = 1 / C2*W0*Q0
R1*R2 = 1 / C1*C2*W0^2
On en tire donc : R1 = (1 / C2*W0*Q0) - R2
et R2^2 - (1 / C2*W0*Q0)*R2 + (1 / C1*C2*W0^2) = 0
R2 = le résultat de ce trinôme

#### Passe Haut
Fonction de transfert est donnée par : 
H(jw) = C1*C2*R1*R2*(jw)^2 / (1+jw(R1+R2) + C1*C2*R1*R2*(jw)^2)

Quand les resistances sont inconnus :
C2 = 1 / ((R1+R2)*W0*Q0)
C1 = (R1+R2)*Q0 / R1*R2*W0

Quand les résistances sont inconnues :
Une contrainte est imposé qui est : C1 >= 4*Q0^2*C2
Donc : R1+R2 = 1 / C2*W0*Q0
R1*R2 = 1 / C1*C2*W0^2
On en tire donc : R1 = (1 / C2*W0*Q0) - R2
et R2^2 - (1 / C2*W0*Q0)*R2 + (1 / C1*C2*W0^2) = 0
R2 = le résultat de ce trinôme

## Installation

```bash
poetry install
...
```

## (Pour les étudiants, à supprimer une fois fait)

### Comment créer le module

1. Créer un nouveau répertoire avec le nom du module
2. Créer un fichier `__init__.py` vide
3. Créer un fichier `__main__.py` vide
4. Mettre à jour le fichier `README.md`
5. Créer un projet Poetry avec `poetry new`
6. Ajouter les fichiers à Git
7. Commit et push
