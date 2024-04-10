#!/usr/bin/env python3

import math, string



########################################
#         PARSING TOOLS
########################################



def PDB_parser(infile):
    """version adaptee du parser PDB qui est retourne
            - un dico PDB (meme format que d'habitude
            input: nom du fichier pdb
            output : dPDB (dico PDB)

            """
    # lecture du fichier PDB
    # -----------------------
    f = open(infile, "r")
    lines = f.readlines()
    f.close()

    # var ini
    # ---------
    dPDB = {}
    atomlist = []
    dPDB["chains"] = []

    # parcoure le PDB
    # -----------------
    for line in lines:
        if line[0:4] == "ATOM":

            # on recupere l'info de la chaine
            chain = line[21]

            # si la chaine n'existe pas, on cree la cle correspondante et on ajoute la chaine a la liste des chaines
            if not chain in dPDB["chains"]:
                dPDB["chains"].append(chain)  # ajout de "chain" a la liste des chaines
                dPDB[chain] = {}  # creation du sous-dico pour la chaine
                # on prepare la structure de donnees pour cette chaine
                dPDB[chain]["reslist"] = []

            # on recupere l'info du residu
            curres = "%s" % (line[22:26]).strip()

            # si le residu pour cette chaine "chain" n'existe pas, on cree la cle correspondante et on ajoute le res a la liste des res
            if not curres in dPDB[chain]["reslist"]:
                dPDB[chain]["reslist"].append(curres)
                dPDB[chain][curres] = {}
                # on prepare la structure de donnees pour ce residu
                dPDB[chain][curres]["atomlist"] = []
                # on recupere l'info du residu
                dPDB[chain][curres]["resname"] = line[17:20].strip()


            # on recupere les info pour l'atome de ce res de cette chaine (type atomique + coords x, y, z)
            atomtype = line[12:16].strip()
            dPDB[chain][curres]["atomlist"].append(atomtype)
            dPDB[chain][curres][atomtype] = {}
            dPDB[chain][curres][atomtype]["x"] = float(line[30:38])
            dPDB[chain][curres][atomtype]["y"] = float(line[38:46])
            dPDB[chain][curres][atomtype]["z"] = float(line[46:54])
            dPDB[chain][curres][atomtype]["id"] = line[6:11].strip()
            dPDB[chain][curres][atomtype]["bfactor"] = line[60:67].strip()


    return dPDB


def writePDB(dPDB, filout = "out.pdb", bfactor = False) :
    """purpose: according to the coordinates in dPDB, writes the corresponding PDB file.
       If bfactor = True, writes also the information corresponding to the key bfactor
       of each residue (one key per residue) in dPDB.
       input: a dico with the dPDB format
       output: PDB file.
    """

    fout = open(filout, "w")

    for chain in dPDB["chains"]:
        for res in dPDB[chain]["reslist"] :
            for atom in dPDB[chain][res]["atomlist"] :
                if bfactor :
                    fout.write("ATOM  %5s  %-4s%3s %s%4s    %8.3f%8.3f%8.3f  1.00%7.3f X X\n"%(
                        dPDB[chain][res][atom]["id"], atom, dPDB[chain][res]["resname"],chain, res,
                        dPDB[chain][res][atom]["x"], dPDB[chain][res][atom]["y"],
                        dPDB[chain][res][atom]["z"],float(dPDB[chain][res][atom]["bfactor"]) ))
                else:
                    fout.write("ATOM  %5s  %-4s%3s %s%4s    %8.3f%8.3f%8.3f  1.00  1.00 X X\n"%(
                        dPDB[chain][res][atom]["id"], atom, dPDB[chain][res]["resname"],
                        chain, res,dPDB[chain][res][atom]["x"], dPDB[chain][res][atom]["y"],
                        dPDB[chain][res][atom]["z"] ))
                    
    fout.close()


def initBfactor(dPDB, val = 0):
    """purpose: initiation of the bfactor key for each residue
       input: a dico with the dPDB format
    """

    for chain in dPDB["chains"]:
        for res in dPDB[chain]["reslist"]:
            for atom in dPDB[chain][res]["atomlist"] :
                dPDB[chain][res][atom]["bfactor"] = val


            

