# -*- coding: utf-8 -*-
import sys

if len(sys.argv) > 1:
    l = list()
    g = dict()
    #lecture des fichiers
    with open(sys.argv[1], 'r') as pipe:
        x = (pipe.readlines())
    #enregistrements
    for i in x:
        if not i in l:
            l.append(i)
            g[i] = 1
        else:
            g[i] += 1
    #affichage
    for gene, nb in sorted(g.items(), key=lambda x: x[1], reverse=True):
        if nb > 1:
            print(gene,nb)

else:
    print("Entrer la liste des gènes en paramètre")
    sys.exit(1)


