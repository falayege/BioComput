Tout le code est présent dans le fichier get_snps.py.
Il y a une option à passer en ligne de commande lorsqu'on exécute le fichier.
Les différentes options sont:

    - "histograms k1 k2" où k1 et k2 sont des entiers.
    Cela permet d'afficher les histogrammes des occurences des k-mers dans le
    fichier "salmonella-enterica.reads.fna" pour k compris entre k1 et k2.
    Elle affiche sur la sortie standard la moyenne, la médiane, et la variance.

    - "poisson" qui permet d'afficher l'histogramme des lois de Poisson des 2
    variables aléatoires détaillées dans le rapport.

    - "comparison" qui écrit dans un fichier "mutated.fasta" tous les k-mers
    (sans erreur de lecture) du fichier "salmonella-enterica-variant.reads.fna"
    où nous avons trouvé des mutations. Le fichier "mutated.fasta" doit avoir 
    été créé au préalable.