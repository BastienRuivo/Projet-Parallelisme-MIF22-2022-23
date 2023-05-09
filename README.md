# Projet MIF22

## Compilation

Il est possible de compiler le projet avec les commandes

```makefile
make all # Compile tout les exécutables
make mpi # Compile uniquement les executables MPI
make gcc # Compile uniquement les executables monothread

make clean # Supprime tout les exécutables
```

## Exécution

Chaque fichier peut être exécuté avec la commande

```bash
./run_NAIVE # pour la version de base
./run_SEQ # pour la version séquentielle optimisé
./run_MPI # pour la version MPI
./run_OMP # pour la version OpenMP
```

Il existe également les commandes suivante pour lancer les tests sur les différentes versions

```bash
make run DIM=Dimension # Lance les tests sur les versions optimisés
make run_with_naive DIM=Dimension # Lance les tests sur les versions optimisés et sur la version de base
```

  
## Partie sequentielle

### Fonction Compute

- Tout d'abord, les if de la fonciton compute peuvent simplement être remplacés par l - 1 et m - 1 dans l'appel au tableau
  - Ceci retire tous les if et réduit la fonction compute a une ligne, en faisant gagner un temps considérable.
  - Cette fonction ne faisant qu'une ligne, je décide de la supprimer complètement, et de rapatrier le calcul dans kernel.

### Fonction Kernel

- Ici, on peut directement voir que la première boucle for ne fait que N fois le calcul, et donc on peut la supprimer et multiplier le résultat par N.
  - 6431 étant diviser par N, et le résultat multiplier par N, on peut donc retirer la division et la multiplication, et multiplier 6431 par 4, supprimant aussi cette multiplication, en la mettant en dur, et aussi mettre la constant directement en double pour retirer le cast inutile.
- Les boucle For du masque peuvent être remplacées en ajoutant += 9 fois succesivement le resultat par le masque, ce retrait fait énormément perdre en lisibilité mais fait gagner un temps considérable.
- La multiplication par 25724 (6431 * 4) est appliqué sur chaque partie du masque, on peut donc la factoriser, et appliquer une fois a l'emsemble du masque en x y.
- J'ai aussi imaginer retirer la deuxième boucle for, mais le calcul des indices est plus couteux que la boucle for, et donc ne fait pas gagner de temps.

### Résultat séquentielles

| Taille | Temps d'execution V0 | Temps d'execution V1 | Gain | Efficacité |
| ------ | ----------------- | -------------------------- | ---- | ---------- |
| 8      | 0.000000          | 0.000000                   | X    | X          |
| 256    | 0.039             | 0.000053                   |  1000   | 1000    |
| 4096   | 176               | 0.017474                   |  10072  | 10072   |

L'efficacité et le gain sont égaux car on compare deux versions séquentielles
La compilation inclu les optimisation gcc via -O3

On voit bien que sans surprise, les gains sont significative en supprimant le plus de calcul possible, et en factorisant les calculs.

Pour la suite, il n'y aura plus de comparaison avec la V0 car elle est trop lente pour être comparée.

## Partie parallele MPI

### Vérification

Vu qu'en monothread, le calcul est juste, j'ai changé le calcul de vérification pour comparé la matrice final rendue par le programme, et la même matrice calculé en monothread, au lieu de vérifier si une petite matrice fonctionne.

### Fonction Kernel MPI

- Ajout d'un argument taille de colonne, pour pouvoir avoir un nombre différent de ligne et de colonne dans la fonction.

### Main, et agencement des processus

- Plutôt que les processus gère X colonnes, comme décrit dans le sujet, j'ai préféré prendre pour chaque procéssus la ligne avant, nLigne, et aprés. nLigne est calculé en fonction du nombre de processus, et de la taille de la matrice.
  - Ceci permet de ne pas rendre la gestion de la mémoire très complexe, vu que chaque ligne est continue en mémoire, on envoit un bloc X continue de donnée, alors qu'en mode colonne, il y aurait des trous.
  - Celles si vont calculer avec la fonction KERNEL_MPI les x + 2 lignes qui leurs sont dédiés
  - Si le nombre de processus ne divise pas parfaitement la taille de la matrice, les processus R premier processus vont calculer une ligne de plus, pour que la matrice soit calculé en entier, le mieux est quand même de prendre un nombre de processus qui divise parfaitement la taille de la matrice. Il faut pour ça tenir compte du fait que la première et la dernière ligne ne sont pas utiliser pour le calcul, et donc que le nombre de processus doit diviser la taille de la matrice - 2.
- Le processus Main (le 0) va avoir le rôle de maitre et initialisé la matrices avant de dispatché les lignes aux autres processus.
  - Il va recevoir les lignes calculé par les autres processus, et les remettre dans la matrice.
  - Il se charge également de calculer le temps écoulé entre le début et la récupération des données calculés par les autres processus.
  - Il va enfin calculer la matrice de vérification, et la comparer avec la matrice calculé par les autres processus. 

### Résultat MPI

| Taille | Temps d'execution V1 | Temps d'execution V2 | Gain     | Efficacité |
| ------ | ----------------- | -------------------------- | ----  | ---------- |
| 16     | 0.000001          | 0.09                   | 1e-05     |  1e-06     |
| 4096   | 0.021             | 0.13                   | 0.16     | 0.016       |
| 16000  | 0.27             | 1.32                    | 0.20      | 0.023       |
| 18000  | 0.34             | 1.60                    | 0.21      | 0.024       |

On voit bien le surcout du au temps de communication entre les processur MPI, qui est très important et qui fait perdre beaucoup de temps, a un tel point qu'il est préférable pour les matrice de ne pas utiliser MPI, et de rester en monothread.

On voit également que le gain croit avec la taille de la matrice, car le surcout de communication du a MPI augmente moins vite que le temps de calcul, et donc le gain augmente, cependant le gain est faible, et il faudrait une matrice de taille bien plus grande pour que l'algo MPI soit plus rapide que l'algo monothread, a suppose que celui ci continue de croitre avec la taille de la matrice.

## Partie parallele OpenMP

Il n'y a pas grand chose que soit possible de faire pour optimiser la fonction kernel,
Tout d'abord, j'ai passé le masque de convolution et la in matrice en const, pour que openMP sache qu'elles sont readonly et puisse les optimiser en conséquence, ceci pour éviter tout histoire de section critique, et de synchronisation que le compilateur pourrait faire automatiquement.
Ensuite j'ai ajouté un pragma omp parallel for, pour que openMP puisse paralléliser la boucle for, et donc que chaque thread puisse calculer une ligne de la matrice.
J'ai également parallélisé la boucle de l'ajout du processus 0, etc... mais cela n'a pas d'impact sur le temps d'execution, car le temps de calcul est négligeable par rapport au temps de communication.
J'ai essayé de collapse la boucle for, mais cela fait perdre en temps d'execution.

### Résultat OpenMP

Pour tester openMP, j'ai limité MPI a un seul processus, et j'ai utilisé openMP pour paralléliser le calcul de la matrice.

| Taille | Temps d'execution V1 | Temps d'execution V3 | Gain     | Efficacité |
| ------ | ----------------- | -------------------------- | ----  | ---------- |
| 16     | 0.000001          | 0.09                   | 0.0029    |  0.00029  |
| 4096   | 0.017             | 0.07                   | 0.24      | 0.016     |
| 16000  | 0.27             | 1.16                    | 0.23     | 0.023     |

On voit bien que le gain est bien plus important que pour MPI, car openMP ne fait pas de communication entre les threads, et donc le surcout est bien plus faible, et donc la perte de temps du a la parallélisation est plus faible, cependant le gain est toujours inférieur a 1 pour des matrices de taille raisonnable, et donc il est préférable de rester en monothread pour ses matrices, même si on peut imagine que l'algo vaille le coup pour des matrices gigantesqte a plusieurs millions de ligne type machine learning.