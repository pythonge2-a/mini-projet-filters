[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/oOQR1xPR)
# Mini-projet PythonGE2 
## Librairie "Pythonfilters"

### Membres du groupe :
- David Vuillemier
- Sébastien Pfister
- Maxime Magnenat
- Francisco Oliveira Barbosa
- Maxime Otero

### Description du projet 
Le but est de faire une librairie open source qui permet de calculer les différents paramètres (w, r1,r2,c1,c2) d'un filtre passif.
Les types de filtres abordés seront :
- Tchebychev
- Bessel
- Butterworth

Nous utliserons des cellules à gain fixe (Sallen & Key) qui permettront de réaliser des filtres à gains fixe de type passe-bas, passe-haut,
passe-bande et coupe-bande.

### Fonctionnalités
Notre bibliothèque se découpera sous la forme de fonctions prenant en paramètre le type de filtre, les pulsations voulues, etc....

Il sera possible pour chaque filtre calculé de déterminer le diagramme de bode (si voulu).

Une autre fonctionnalité sera d'afficher le schéma du circuit, avec les valeurs de résistances.

### Détails
Nous explorerons la simulation par le modèle Spice et verront si c'est nécessaire ou non dans notre projet.

On importera un fichier qui permettra de donner les pulsations et les facteurs de qualités des filtres à un certain ripple.

On s'arretêra à l'ordre 10.

Nous aimerions aussi traiter le cas des filtres RLC passifs standards.

### Déploiement

Nos objectifs sont de commencer cette librairie qui pourrait être développée par la suite par d'autres groupes pour étendre les fonctionnalités.

Nous devrons respecter les délais imposés par notre professeur et lui fournir un travail qui correspond "au mieux" au cahier des charges.

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
