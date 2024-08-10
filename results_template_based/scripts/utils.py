import os
import re
import sys
import colorsys
from pymol import cmd,stored
def load_all_in_dir(dir):
    """load all files in a folder
    load_all_in_dir dir
    """
    for i in os.listdir(dir):
        cmd.load(os.path.join(dir,i))
def load_pdb_in_dir(dir):
    """load all PDB files in a folder
    load_pdb_in_dir dir
    """
    for i in os.listdir(dir):
        if os.path.splitext(i)[1] == ".pdb":
            cmd.load(os.path.join(dir,i))

def load_top_rank_pdb_in_dir(dir, top_rank=10, consensus_file=None, group=""):
    """load top10 PDB files in a folder, a rank is supposed to be contained in the name
    load_top10_pdb_in_dir dir. If the consensus file name is given, complexes listed in
    the file will be excluded. If group  = <goup name>, all loaded complexes will be put
    together in a group named <group name>
    """
    def get_rank_from_the_name(pdb_name):
        res = re.search('top(\d+)',pdb_name)
        if res:
            rank = int(res.group(1))
        else:
            res = re.search('[\D](\d+)\.pdb',pdb_name)
            if res:
                rank = int(res.group(1))
            else:
                res = re.search('[\D](\d+).+\.pdb',pdb_name)
                if res:
                    rank = int(res.group(1))
        return rank

    if not os.path.exists(dir):
        return

    if consensus_file and consensus_file != 'None':
        top10_consensus = []
        f=open(consensus_file,'r')
        for line in f.readlines():
            top10_consensus.append(line.split()[1])
        f.close()
    else:
        top10_consensus = []

    top_rank = int(top_rank)
    lst_pdb = []
    for i in os.listdir(dir):
        if os.path.splitext(i)[1] == ".pdb":
            lst_pdb.append(i)
    sorted_lst_pdb =  sorted(lst_pdb,key=get_rank_from_the_name)
    to_group = []
    for i in sorted_lst_pdb[0:min(top_rank,len(sorted_lst_pdb))]:
        if i not in top10_consensus:
            to_group.append(os.path.splitext(i)[0])
            cmd.load(os.path.join(dir,i))

    if group:
        cmd.group(group," ".join(to_group))
        cmd.do("group {}, action=open".format(group))


def load_top_rank_pdb_in_consensus(consensus_file="InterEvDock_Top10Consensus.txt", top_rank=10, group=""):
    """load top10 PDB files in the consensus file, a rank is supposed to be contained in the name
    load_top_rank_pdb_in_consensus consensus_file
    """
    to_group = []
    corres = {"IES":"IES_models","IESh":"IESh_models","FRODOCK":"frodock_models","SPP":"SPP_models","SPPh":"SPPh_models","ISCh":"ISCh_models","ISC":"ISC_models"}
    f=open(consensus_file,'r')
    it = 0
    for line in f:
        if line.startswith("#"):
            continue
        if it >= int(top_rank):
            break
        best_cplx = line.split()[1]
        res = re.search("Complex_([\D]+)\d+.pdb",best_cplx)
        if res:
            dir = corres[res.group(1)]
        else:
            continue
        cmd.load(os.path.join(dir,best_cplx))
        to_group.append(os.path.splitext(best_cplx)[0])
        it += 1
    f.close()

    if group:
        cmd.group(group," ".join(to_group))

def load_top_rank_pdb_in_dir_with_pml(dir, top_rank='10', group="", pml_rootname=""):
    """load top10 PDB files in a folder, a rank is supposed to be contained in the name
    load_top10_pdb_in_dir dir.
    If group  = <goup name>, all loaded complexes will be put
    together in a group named <group name>
    For every pdb opened, if a pml is provided, it will be run just after loading of every pdb
    A priori no consensus is used for this script but it can be added if required
    """
    def get_rank_from_the_name(pdb_name):
        res = re.search('top(\d+)',pdb_name)
        if res:
            rank = int(res.group(1))
        else:
            res = re.search('[\D](\d+)\.pdb',pdb_name)
            if res:
                rank = int(res.group(1))
            else:
                res = re.search('[\D](\d+).+\.pdb',pdb_name)
                if res:
                    rank = int(res.group(1))
        return rank

    if not os.path.exists(dir):
        return

    top_rank = int(top_rank)
    lst_pdb = []
    lst_rank = []
    for i in os.listdir(dir):
        if os.path.splitext(i)[1] == ".pdb":
            lst_pdb.append(i)
            lst_rank.append(get_rank_from_the_name(i))
    sorted_lst_pdb = sorted(lst_pdb,key=get_rank_from_the_name)
    lst_rank.sort()
    if pml_rootname:
        list_pml = []
        for rank in lst_rank:
            pmlfile = os.path.splitext(pml_rootname)[0] + "{}.pml".format(rank)
            list_pml.append(pmlfile)
    print("INFO: One after the other, pdb files are opened and the pml analyzing the contacts respected will be run")
    print("{}".format(sorted_lst_pdb[0:min(top_rank,len(sorted_lst_pdb))]))
    print("{}".format(list_pml))
    for index,i in enumerate(sorted_lst_pdb[0:min(top_rank,len(sorted_lst_pdb))]):
        cmd.load(os.path.join(dir,i))
        if group:
            cmd.group(group,os.path.splitext(i)[0])
        if pml_rootname:
            cmd.do("run {}".format(list_pml[index]))
            if group:
                cmd.group(group,"c_" + os.path.splitext(i)[0])
    if group:
        cmd.do("group {}, action=open".format(group))



def select_top_5_residue_consensus(consensus_file="consensus_top5_residues.txt"):
    """select top5 residues of each chain based on residues listed in the consensus file
    """
    f=open(consensus_file,'r')
    best_res_on_chainA = []
    best_res_on_chainB = []
    for line in f.readlines():
        if line[:1] == "#":
            continue
        line = line.split()
        # detect format (i.e. Y45 or Y.45.A?)
        if '.' in line[1]:
            best_res_on_chainA.append("chain "+line[1].split('.')[2]+" and resi "+line[1].split('.')[1])
        else:
            best_res_on_chainA.append("chain A and resi "+line[1][1:])
        if '.' in line[2]:
            best_res_on_chainB.append("chain "+line[2].split('.')[2]+" and resi "+line[2].split('.')[1])
        else:
            best_res_on_chainB.append("chain B and resi "+line[2][1:])
    cmd.select("top5_res_inputA","(%s)"%" or ".join(best_res_on_chainA))
    cmd.select("top5_res_inputB","(%s)"%" or ".join(best_res_on_chainB))
    f.close()


def align_all(selection='sele'):
    """ align all objects onto a selection (def='sele') or align all other objects onto an object
    """
    obj_list = cmd.get_object_list(selection='(all)')
    for obj in obj_list:
        if obj != selection:
            cmd.align(obj,selection)

cmd.extend("load_all_in_dir", load_all_in_dir)
cmd.extend("load_pdb_in_dir", load_pdb_in_dir)
cmd.extend("load_top_rank_pdb_in_dir", load_top_rank_pdb_in_dir)
cmd.extend("load_top_rank_pdb_in_dir_with_pml", load_top_rank_pdb_in_dir_with_pml)
cmd.extend("load_top_rank_pdb_in_consensus", load_top_rank_pdb_in_consensus)
cmd.extend("select_top_5_residue_consensus", select_top_5_residue_consensus)
cmd.extend("align_all", align_all)
cmd.auto_arg[0]['align_all'] = [lambda: cmd.Shortcut(cmd.get_names("objects")+cmd.get_names("selections")),'object or selection name', '']

def color_b_R4S_consurf(selection, keep_selections=False):
    """
    """
    # These are constants
    min = 0.0
    max = 100.0
    n_colours = 10
    colours = [
    [1, 1, 0.58823529] ,
    [0.062745098,0.78431373,0.81960784],
    [0.54901961,1,1],
    [0.84313725,1,1],
    [0.91764706,1,1],
    [1,1,1],
    [0.98823529,0.92941176,0.95686275],
    [0.98039216,0.78823529,0.87058824],
    [0.94117647,0.49019608,0.67058824],
    [0.62745098,0.14509804,0.37647059] ]

    bin_size = ((max - min)) / n_colours
    # Loop through colour intervals
    for i in range(n_colours):

        lower = min + i * bin_size
        upper = lower + bin_size
        colour = colours[i]

        # Print out B-factor limits and the colour for this group
        print(lower, " - ", upper, " = ", colour)

        # Define a unique name for the atoms which fall into this group
        group = selection + "_group_" + str(i+1)

        # Compose a selection command which will select all atoms which are
        #       a) in the original selection, AND
        #       b) have B factor in range lower <= b < upper
        sel_string = selection + " & ! b < " + str(lower)

        if(i < n_colours - 1):
            sel_string += " & b < " + str(upper)
        else:
            sel_string += " & ! b > " + str(upper)

        # Select the atoms
        cmd.select(group, sel_string)

        # Create a new colour
        colour_name = "colour_" + str(i+1)
        cmd.set_color(colour_name, colour)

        # Colour them
        cmd.color(colour_name, group)

        # Delete selection
        if not keep_selections:
            cmd.delete(group)


    # Create new colour for insufficient sequences
    insuf_colour = [1, 1, 0.58823529]
    cmd.set_color("insufficient_colour", insuf_colour)

    # Colour atoms with B-factor of 10 using the new colour
    cmd.select("insufficient", selection + " & b = 10")
    cmd.color("insufficient_colour", "insufficient")

    cmd.delete("insufficient")
cmd.extend("color_b_R4S_consurf", color_b_R4S_consurf)
cmd.auto_arg[0]['color_b_R4S_consurf'] = [lambda: cmd.Shortcut(cmd.get_names('all')),'object name', '']

def colour_consurf(selection="all", keep_selections=False):
    """color a selection by Consurf coloring method
    """
    # These are constants
    min = 0.0
    max = 9.0
    n_colours = 10
    colours = [
    [1, 1, 0.58823529] ,
    [0.062745098,0.78431373,0.81960784],
    [0.54901961,1,1],
    [0.84313725,1,1],
    [0.91764706,1,1],
    [1,1,1],
    [0.98823529,0.92941176,0.95686275],
    [0.98039216,0.78823529,0.87058824],
    [0.94117647,0.49019608,0.67058824],
    [0.62745098,0.14509804,0.37647059] ]

    bin_size = ((max - min) + 1) / n_colours

    # Loop through colour intervals
    for i in range(n_colours):

        lower = min + i * bin_size
        upper = lower + bin_size
        colour = colours[i]

        # Print out B-factor limits and the colour for this group
        print(lower, " - ", upper, " = ", colour)

        # Define a unique name for the atoms which fall into this group
        group = selection + "_group_" + str(i+1)

        # Compose a selection command which will select all atoms which are
        #       a) in the original selection, AND
        #       b) have B factor in range lower <= b < upper
        sel_string = selection + " & ! b < " + str(lower)

        if(i < n_colours - 1):
            sel_string += " & b < " + str(upper)
        else:
            sel_string += " & ! b > " + str(upper)

        # Select the atoms
        cmd.select(group, sel_string)

        # Create a new colour
        colour_name = "colour_" + str(i+1)
        cmd.set_color(colour_name, colour)

        # Colour them
        cmd.color(colour_name, group)

        # Delete selection
        if not keep_selections:
            cmd.delete(group)


    # Create new colour for insufficient sequences
    insuf_colour = [1, 1, 0.58823529]
    cmd.set_color("insufficient_colour", insuf_colour)

    # Colour atoms with B-factor of 10 using the new colour
    cmd.select("insufficient", selection + " & b = 10")
    cmd.color("insufficient_colour", "insufficient")

    cmd.delete("insufficient")
cmd.extend("colour_consurf", colour_consurf)
cmd.auto_arg[0]['colour_consurf'] = [lambda: cmd.Shortcut(cmd.get_names('all')),'object name', '']

# usage analyze_complex object(or selection or all)
def analyze_complex(object='all', name=None, ch1=None, ch2=None):
    """show interface contact for 2 chains for a complex (or all complexes);
    you can select 2 chains and then do 'analyse_complex sele'; you can give a name"
    for your interface so that it will not be overwritten when you analyze the second interface"
    analyze_complex(object='all', name=None)"
    """
    def ana_int(o,ii,name=None,ch1=None, ch2=None):
        if ch1 == None and ch2 == None:
            chains = cmd.get_chains(o)
        else:
            ch1 = re.sub('_',' ',ch1)
            ch2 = re.sub('_',' ',ch2)
            chains = [str(ch1), str(ch2)]

        j= str(ii+1)
        if name != None:
            j= name
        print(str(o)+" chains: "+ str(chains))
        if len(chains)<=1:
            return
        elif len(chains)>=2:
            chain1 = "'"+chains[0]+"'"
            chain2 = "'"+chains[1]+"'"
            cmd.select("chain1","/%s//%s"%(o,chain1))
            cmd.select("chain2","/%s//%s"%(o,chain2))
            cmd.select("face1","/%s//%s"%(o,chain1))
            print("/%s//%s"%(o,chain1))
            cmd.select("face2","/%s//%s"%(o,chain2))
            cmd.select("face1","(byres (face1 around 4)) and /%s//%s"%(o,chain2),enable=1)
            cmd.select("face2","(byres (face2 around 4)) and /%s//%s"%(o,chain1),enable=1)
            cmd.show("sticks","face2")
            cmd.show("sticks","face1")
            cmd.select("interface_%s"%(j),"face1 or face2")
            cmd.enable("interface_%s"%(j))
            cmd.dist("interface_%s_polar_conts"%(j),"interface_%s"%(j),"interface_%s"%(j),quiet=1,mode=2,label=0,reset=1)
            cmd.enable("interface_%s_polar_conts"%(j))
            #cmd.spectrum("count",selection="chain1",byres=1)
            #cmd.spectrum("count",selection="chain2",byres=1)
            util.cnc("chain1")
            util.cnc("chain2")
            cmd.do("delete chain1 chain2 face1 face2 interface_%s"%(j))
    # end function ana_int
    if object=='all':
        cmd.hide("everything","(all)")
        cmd.show("cartoon"   ,"(all)")
        lobj = cmd.get_object_list('(all)')
        for ii,o in enumerate(lobj):\
                ana_int(o,ii,ch1=ch1,ch2=ch2)
    else:
        cmd.hide("everything",object)
        cmd.show("cartoon"   ,object)
        ana_int(object,9998, name,ch1=ch1,ch2=ch2)

cmd.extend("analyze_complex", analyze_complex)
cmd.auto_arg[0]['analyze_complex'] = [lambda: cmd.Shortcut(cmd.get_names('all')),'object name', '']

def save_mol(object='all', prefix=''):
    """save a molecule (or a selection, or use 'all' for all objects but not any selection)
    in the current dir. For example, prefix='pymol_' and the object is called 'xxx', the file name
    will be pymol_xxx.pdb.
    """
    prefix = re.sub("\'\"","",prefix)
    if object=='all':
        list_obj = cmd.get_object_list('(all)')
        for obj in list_obj:
            cmd.save("%s%s.pdb"%(prefix,obj), obj, 0)
    else:
        cmd.save("%s%s.pdb"%(prefix,object), object, 0)
cmd.extend("save_mol", save_mol)
cmd.auto_arg[0]['save_mol'] = [lambda: cmd.Shortcut(cmd.get_names('all')),'object name', '']

def interfacial_waters(my_object,chain1=None,chain2=None,dist='3.5'):
    """Select interfacial water molecules (which are around of each chain at a distance < dist (def=3.5A)).
    If there are more than 2 chains in this object, 1st and 2nd chains are selected by default.
    You can also name them explicitly. "interfacial_waters all[, ch1, ch2, dist]" will deal with all objects
    Normally a water has a residue name HOH but it chould also be TIP3` generated by Haddock
    """
    def get_interfacial_waters(my_object,chain1=None,chain2=None,dist='3.5'):
        chains = cmd.get_chains(my_object)
        if chain1 == None or chain2 == None:
            if len(chains)==2:
                chain1 = "'"+chains[0]+"'"
                chain2 = "'"+chains[1]+"'"
            elif len(chains)==1:
                print("Only one chain in this object. Quit.")
                return 1
            else:
                chain1 = "'"+chains[0]+"'"
                chain2 = "'"+chains[1]+"'"
                print("chains %s %s selected. You can name two chains by their IDs"% (chain1,chain2))
        text = "(byres ((/%s//%s and not /%s///HOH and not /%s///TIP3`) around %s))"%(my_object,chain1,my_object,my_object,dist)
        text += " and (byres ((/%s//%s and not /%s///HOH and not /%s///TIP3`) around %s))"%(my_object,chain2,my_object,my_object,dist)
        text += " and (/%s///HOH or /%s///TIP3`)"%(my_object,my_object)
        cmd.select("sele_interf_waters",text)
        cmd.do("show spheres, sele_interf_waters")
        print(text)
        print("Number of interfacial waters: %d"%cmd.count_atoms("sele_interf_waters and (name O or name OH2 or name OW)"))
        return 0
    if my_object is None:
        print("input your object")
        return 0
    elif my_object=='all':
        for one_obj in cmd.get_object_list('(all)'):
            print(one_obj)
            get_interfacial_waters(one_obj,chain1,chain2,dist)
    else:
        print("object: %s"%my_object)
        get_interfacial_waters(my_object,chain1,chain2,dist)


cmd.extend("interfacial_waters", interfacial_waters)
cmd.auto_arg[0]['interfacial_waters'] = [lambda: cmd.Shortcut(cmd.get_names('all')),'object name', '']

def waters_near_interface(my_object,chain1=None,chain2=None,dist='5.5'):
    """show spheres of all the waters which are of distance < dist(def=5.5) from interface of 2 chains (def=1st and 2nd chains)
    for your given object. This selection has renamed as sele_waters_near_interf.
    "waters_near_interface all[, ch1, ch2, dist]" will deal with all objects
    Normally a water has a residue name HOH but it chould also be TIP3` generated by Haddock
    """
    def get_waters_near_interface(my_object,chain1=None,chain2=None,dist='5.5'):
        chains = cmd.get_chains(my_object)
        if chain1 == None or chain2 == None:
            if len(chains)==2:
                chain1 = "'"+chains[0]+"'"
                chain2 = "'"+chains[1]+"'"
            elif len(chains)==1:
                print("Only one chain in this object. Quit.")
                return 1
            else:
                chain1 = "'"+chains[0]+"'"
                chain2 = "'"+chains[1]+"'"
                print("chains %s %s selected. You can name two chains by their IDs"% (chain1,chain2))
        # chain2 residues in contact with chain1
        text1 = "(byres ((/%s//%s and not /%s///HOH and not /%s///TIP3`) around %s))"%(my_object,chain1,my_object,my_object,dist)
        text1 += " and (/%s//%s) and not /%s///HOH and not /%s///TIP3`"%(my_object,chain2,my_object,my_object)
        print(text1)
        # chain1 residues in contact with chain2
        text2 = "(byres ((/%s//%s and not /%s///HOH and not /%s///TIP3`) around %s))"%(my_object,chain2,my_object,my_object,dist)
        text2 += " and (/%s//%s) and not /%s///HOH and not /%s///TIP3`"%(my_object,chain1,my_object,my_object)
        print(text2)
        # waters near chain1 and chain2 in interface
        text = "byres ( ((%s) or (%s)) around %s) and (/%s///HOH or /%s///TIP3`)"%(text1, text2, dist, my_object,my_object)
        cmd.select("sele_waters_near_interf",text)
        cmd.do("show spheres, sele_waters_near_interf")
        print(text)
        print("Number of waters near interface: %d"%cmd.count_atoms("sele_waters_near_interf and (name O or name OH2 or name OW)"))
        return 0
    if my_object is None:
        print("input your object")
        return 0
    elif my_object=='all':
        for one_obj in cmd.get_object_list('(all)'):
            print(one_obj)
            get_waters_near_interface(one_obj,chain1,chain2,dist)
    else:
        print("object: %s"%my_object)
        get_waters_near_interface(my_object,chain1,chain2,dist)


cmd.extend("waters_near_interface", waters_near_interface)
cmd.auto_arg[0]['waters_near_interface'] = [lambda: cmd.Shortcut(cmd.get_names('all')),'object name', '']

def remove_waters_far_from_interface(my_object,chain1=None,chain2=None,dist='5.5'):
    """Remove all waters which are of distance > dist(def=5.5) from interface of 2 chains (def=1st and 2nd chains)
    for your given object. "remove_waters_far_from_interface all[, ch1, ch2, dist]" will deal with all objects
    Normally a water has a residue name HOH but it chould also be TIP3` generated by Haddock
    """
    def sub_function_remove_waters_far_from_interface(my_object,chain1=None,chain2=None,dist='5.5'):
        chains = cmd.get_chains(my_object)
        if chain1 == None or chain2 == None:
            if len(chains)==2:
                chain1 = "'"+chains[0]+"'"
                chain2 = "'"+chains[1]+"'"
            elif len(chains)==1:
                print("Only one chain in this object. Quit.")
                return 1
            else:
                chain1 = "'"+chains[0]+"'"
                chain2 = "'"+chains[1]+"'"
                print("chains %s %s selected. You can name two chains by their IDs"% (chain1,chain2))
        # chain2 residues in contact with chain1
        text1 = "(byres ((/%s//%s and not /%s///HOH and not /%s///TIP3`) around %s))"%(my_object,chain1,my_object,my_object,dist)
        text1 += " and (/%s//%s) and not /%s///HOH and not /%s///TIP3`"%(my_object,chain2,my_object,my_object)
        # chain1 residues in contact with chain2
        text2 = "(byres ((/%s//%s and not /%s///HOH and not /%s///TIP3`) around %s))"%(my_object,chain2,my_object,my_object,dist)
        text2 += " and (/%s//%s) and not /%s///HOH and not /%s///TIP3`"%(my_object,chain1,my_object,my_object)
        # # waters near chain1 and chain2 in interface
        text_wat_near = "(byres ( ((%s) or (%s)) around %s) and (/%s///HOH or /%s///TIP3`))"%(text1, text2, dist, my_object,my_object)
        # waters not near chain1 and chain2 in interface
        text = "(not %s) and (/%s///HOH or /%s///TIP3`)"%(text_wat_near, my_object,my_object)
        print(text)
        cmd.select("sele_waters_not_near_interf",text)
        # print(text)
        print("Number of waters not near interface: %d"%cmd.count_atoms("sele_waters_not_near_interf and (name O or name OH2 or name OW)"))
        print("They have been removed")
        cmd.remove("sele_waters_not_near_interf")
        cmd.delete("sele_waters_not_near_interf")
        return 0
    if my_object is None:
        print("input your object")
        return 0
    elif my_object=='all':
        for one_obj in cmd.get_object_list('(all)'):
            print(one_obj)
            sub_function_remove_waters_far_from_interface(one_obj,chain1,chain2,dist)
    else:
        print("object: %s"%my_object)
        sub_function_remove_waters_far_from_interface(my_object,chain1,chain2,dist)


cmd.extend("remove_waters_far_from_interface", remove_waters_far_from_interface)
cmd.auto_arg[0]['remove_waters_far_from_interface'] = [lambda: cmd.Shortcut(cmd.get_names('all')),'object name', '']

def renumber_waters(my_object):
    """renumber all waters (resn HOH or TIP3) in an object and their chain ID becomes 'W' and segID becomes '',
    if the total number of waters are greater than 9999, the function terminates.
    Caution: make sure initial structure has no chain 'W' which contains other molecules than waters.
    """
    def sub_function_renumber_waters(my_object):
        possible_alternative_chainID = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
        tmp_water_chain_id = None
        for ch_id in possible_alternative_chainID:
            if ch_id not in cmd.get_chains(my_object):
                tmp_water_chain_id = ch_id
                break
        if tmp_water_chain_id == None:
            print("No temp chain Id possible for waters.")
            return
        text = "(/%s///HOH or /%s///TIP3`)"%(my_object,my_object)
        stored.list_residues = []
        cmd.iterate(text+ " and (name O or name OH2 or name OW)","stored.list_residues.append(\"/%s/\'%%s\'/\'%%s\'/%%s`%%s\"%%(segi,chain,resn,resi))"%my_object)
        print("%s water molecules found"%len(stored.list_residues))
        if len(stored.list_residues) > 9999:
            print("Total number of water molecules > 9999. Quit.")
            return
        new_index = 1
        for residue in stored.list_residues:
            cmd.select("sele_to_be_changed",residue)
            # change chain ID to tmp_water_chain_id so that this water is not redundant with another existing water
            cmd.alter("sele_to_be_changed","chain='%s'"%tmp_water_chain_id)
            cmd.alter("sele_to_be_changed","resi=str(%d)"%new_index)
            new_index +=1
        cmd.delete("sele_to_be_changed")
        # alter all water: chain = 'W' and segi = ''
        cmd.alter(text, "chain='W'")
        cmd.alter(text,"segi=''")
        print("Renumbering %s finished."%my_object)
    if my_object is None:
        print("input your object")
        return 0
    elif my_object=='all':
        for one_obj in cmd.get_object_list('(all)'):
            print("Renumbering waters in object %s"%one_obj)
            sub_function_renumber_waters(one_obj)
    else:
        print("Renumbering waters in object %s"%my_object)
        sub_function_renumber_waters(my_object)
cmd.extend("renumber_waters", renumber_waters)
cmd.auto_arg[0]['renumber_waters'] = [lambda: cmd.Shortcut(cmd.get_names('all')),'object name', '']

def remove_clashing_waters(my_object, dist=1.2):
    """remove all clashing waters having a neighbor within dist(def=1.2) for an object.
    """
    def sub_remove_clashing_waters(my_object, dist=1.2):
        list_waters_to_be_removed = []
        text = "(/%s///HOH or /%s///TIP3`)"%(my_object,my_object)
        stored.list_residues = []
        cmd.iterate(text,"stored.list_residues.append('/%s/%%s/%%s/%%s`%%s'%%(segi,chain,resn,resi))"%my_object)
        if len(stored.list_residues) > 9999:
            print("Total number of water molecules > 9999. Quit.")
            return
        for residue in stored.list_residues:
            cmd.select("sele_wat",residue + " and (name O or name OH2 or name OW)")
            cmd.select("sele_near","(sele_wat around %s) and %s"%(dist,my_object))
            near_number = cmd.count_atoms("sele_near")
            if near_number == 0:
                continue
            else:
                print("Clashing water found %s" % residue)
                list_waters_to_be_removed.append(residue)
        cmd.delete("sele_wat")
        cmd.delete("sele_near")
        for residue in list_waters_to_be_removed:
            cmd.remove(residue)
    if my_object is None:
        print("input your object")
        return 0
    elif my_object=='all':
        for one_obj in cmd.get_object_list('(all)'):
            print("Try to remove clashing waters in object %s"%one_obj)
            sub_remove_clashing_waters(one_obj, dist)
    else:
        print("Try to remove clashing waters in object %s"%my_object)
        sub_remove_clashing_waters(my_object, dist)
cmd.extend("remove_clashing_waters", remove_clashing_waters)
cmd.auto_arg[0]['remove_clashing_waters'] = [lambda: cmd.Shortcut(cmd.get_names('all')),'object name', '']

# main function called from within PyMOL
def color_b(selection='all',item='b',mode='hist',gradient='wyr',nbins=11,sat=1.,value=1.,minimum='',maximum=''):
    """

    AUTHOR

      Robert L. Campbell with enhancements from James Stroud

    USAGE

      color_b(selection='sel',item='b' or 'q'
        gradient='bgr' or 'rgb' or 'bwr' or 'rwb' or 'bmr' or 'rmb' or
        'rw' or 'wr' or 'gw' or 'wg' or 'bw' or wb' or 'gy' or 'yg' or
        'gray' or 'reversegray'
        mode='hist' or 'ramp', [minimum=''],[maximum=20.],
        nbins=11, sat=1.0, value=1.0)

        The "item" argument allows specifying 'b' or 'q' as the item to
        color on.  The "color_q" function is the same as "color_b item=q".

        This function allows coloring of a selection as a function of
        B-value or occupancy, following a gradient of colours.  The
        gradients can be:

        'bgr': blue -> green   -> red
        'rgb': red  -> green   -> blue
        'bwr': blue -> white   -> red
        'rwb': red  -> white   -> blue
        'bmr': blue -> magenta -> red
        'rmb': red  -> magenta -> blue
        'rw' : red -> white
        'wr' : white -> red
        'gw' : green -> white
        'wg' : white -> green
        'bw' : blue -> white
        'wb' : white -> blue
        'gy' : green -> yellow
        'yg' : yellow -> green
        'gray' : black -> white
        'reversegray' : white -> black

        ('rainbow' and 'reverserainbow' can be used as synonyms for
        'bgr' and 'rgb' respectively and 'grey' can be used as a synonym for 'gray').

        The division of B-value ranges can in either of two modes: 'hist' or
        'ramp'. 'hist' is like a histogram (equal-sized B-value increments
        leading to unequal numbers of atoms in each bin). 'ramp' as a ramp
        of B-value ranges with the ranges chosen to provide an equal number
        of atoms in each group.

        You can also specify the lower or upper limits of the data used to determine
        the color bins (minimum,maximum). e.g. color_b my_molecule, minimum=15., maximum=25.

        You can also specify the saturation and value (i.e. the "s" and "v"
        in the "HSV" color scheme) to be used for the gradient. The defaults
        are 1.0 for both "sat" and "value".

        In the case of the gray scale gradients, "sat" sets the minimum intensity
        (normally black) and "value" sets the maximum (normally white)

      usage:
        from within PyMOL do "run color_b.py" to load the function definition.
        Then you can use for example:

            color_b (all),mode=ramp,gradient=wyr,nbins=30,sat=1., value=1.

        to color chains A and B with the Blue-White-Red gradient in 30 colors of equal
        numbers of atoms in each color.
    """

    nbins=int(nbins)
    sat=float(sat)
    value=float(value)
# make sure sat and value are in the range 0-1.0
    sat = min(sat, 1.0)
    sat = max(sat, 0.0)
    value = min(value, 1.0)
    value = max(value, 0.0)

# make sure lowercase
    gradient.lower()
    mode.lower()

# Sanity checking
    if nbins == 1:
        print("\n     WARNING: You specified nbins=1, which doesn't make sense...resetting nbins=11\n")
        nbins=11

    if mode not in ('hist','ramp'):
        print("\n     WARNING: Unknown mode ",mode, "    ----->   Nothing done.\n")
        return
    elif gradient not in ('A', 'wc', 'wcb','bgr','rgb','rainbow','reverserainbow','bwr','rwb',
                          'bmr','rmb','rw','wr','gw','wg','bw','wb','gy','yg','gray','grey','reversegray','reversegrey','wyr'):
        print("\n     WARNING: Unknown gradient: ",gradient, "    ----->   Nothing done.\n")
        return

    print("MODE, GRADIENT, NBINS: %s, %s, %d"% (mode,gradient, nbins))

# get list of B-factors from selection
    m = cmd.get_model(selection)
    sel = []
    b_list = []

    if len(m.atom) == 0:
        print("Sorry, no atoms selected")

    else:
        if item == 'b':
            for i in range(len(m.atom)):
                b_list.append(m.atom[i].b)
        elif item == 'q':
            for i in range(len(m.atom)):
                b_list.append(m.atom[i].q)

        else:
            print("Not configured to work on item %s" % item)
            return

        max_b = max(b_list)
        min_b = min(b_list)
        print("Minimum and Maximum B-values: %f, %f"%(min_b, max_b))

        if mode == 'ramp':
            # color in bins of equal numbers of atoms
            b_list.sort()

            # subtract 0.1 from the lowest B in order to ensure that the single
            # atom with the lowest B value doesn't get omitted
#      b_list[0] = b_list[0] - 0.1

            bin_num = int(len(b_list)/nbins)
#      sel.append(selection + " and (b < " + str(b_list[bin_num]) + " or b = " + str(b_list[bin_num]) + ")")
            sel.append(selection + " and (%s < %4.4g" % (item,b_list[bin_num]) + " or %s = %4.4g" % (item,b_list[bin_num]) + ")")
            for j in range(1,nbins):
#        sel.append(selection + " and b > " + str(b_list[j*bin_num]))
                sel.append(selection + " and %s > %4.4g" % (item,b_list[j*bin_num]))
#         print("Color select: ",sel[j])

        elif mode == 'hist':

# check if minimum or maximum was specified and use the entered values
            if minimum != '':
                min_b = float(minimum)
            if maximum != '':
                max_b = float(maximum)
            # histogram:
            # color in bins of equal B-value ranges
            # subtract 0.1 from the lowest B in order to ensure that the single
            # atom with the lowest B value doesn't get omitted
            bin_width = (max_b - min_b)/nbins
            sel.append(selection + " and (%s < %4.4g" % (item,min_b + bin_width) + " or %s = %4.4g" % (item,min_b + bin_width) + ")")
            for j in range(1,nbins):
                sel.append(selection + " and %s > %4.4g" % (item,min_b + j*bin_width))
                #print("Color select: ",sel[j])

# call the function to create the gradient which returns a list of colours
        colours = make_gradient(sel,gradient,nbins,sat,value,item)

# do the colouring now
        for j in range(nbins):
            print("Color %s select: %s"%(colours[j],sel[j]))
            cmd.color(colours[j],sel[j])

def color_q(selection="all",mode="hist",gradient="wg",nbins=11,sat=1.,value=1.,minimum='',maximum=''):
    """

    USAGE

      color_q(selection,gradient,mode,nbins,sat,value) ='sel',
        gradient='bgr' or 'rgb' or 'bwr' or 'rwb' or 'bmr' or 'rmb'
        'rw' or 'wr','gw' or 'wg' or 'bw' or 'wb' or 'gy' or 'yg' or 'gray' or 'reversegray'
        mode='hist' or 'ramp', q0=0.,q1=1.0,
        nbins=11, sat=1.0, value=1.0)

        This function allows coloring of a selection as a function of
        occupancy.  See color_b for details.
    """
    item='q'
    color_b(selection,item,mode,gradient,nbins,sat,value,minimum,maximum)

# function for creating the gradient
def make_gradient(sel,gradient,nbins,sat,value,item):
    if gradient == 'bgr' or gradient == 'rainbow':
        col=[]
        coldesc=[]
        for j in range(nbins):
            # must append the str(sel[j]) to the color name so that it is unique
            # for the selection
            coldesc.append('col' + item + gradient + str(j))
            # coldesc.append('col' + str(sel[j]) + str(j))

            # create colors using hsv scale (fractional) starting at blue(.6666667)
            # through red(0.00000) in intervals of .6666667/(nbins -1) (the "nbins-1"
            # ensures that the last color is, in fact, red (0)
            # rewrote this to use the colorsys module to convert hsv to rgb
            hsv = (colorsys.TWO_THIRD - colorsys.TWO_THIRD * float(j) / (nbins-1), sat, value)
            #convert to rgb and append to color list
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])
            # print("defined as ", str(sel[j]))



    elif gradient == 'rgb' or gradient == 'reverserainbow':
        col=[]
        coldesc=[]
        for j in range(nbins):
            # must append the str(sel[j]) to the color name so that it is unique
            # for the selection
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))

            # create colors using hsv scale (fractional) starting at red(.00000)
            # through blue(0.66667) in intervals of .6666667/(nbins -1) (the "nbins-1"
            # ensures that the last color is, in fact, red (0)
            # rewrote this to use the colorsys module to convert hsv to rgb
            hsv = (colorsys.TWO_THIRD * float(j) / (nbins-1), sat, value)
            #convert to rgb and append to color list
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])

    elif gradient == 'bmr':
        col=[]
        coldesc=[]
        for j in range(nbins):
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))
            # create colors in a gradient from blue through magenta to red
            rgb = [min(1.0, float(j)*2/(nbins-1)), 0.0, min(1.0, float(nbins-j-1)*2/(nbins-1))]

            # convert rgb to hsv,  modify saturation and value and convert back to rgb
            hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
            hsv[1] = hsv[1]*sat
            hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])

    elif gradient == 'rmb':
        col=[]
        coldesc=[]
        for j in range(nbins):
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))
            # create colors in a gradient from red through magenta to blue
            rgb = [min(1.0, float(nbins-j-1)*2/(nbins-1)), 0.0, min(1.0, float(j)*2/(nbins-1))]

            # convert rgb to hsv,  modify saturation and value and convert back to rgb
            hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
            hsv[1] = hsv[1]*sat
            hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])

    elif gradient == 'rw':
        col=[]
        coldesc=[]
        for j in range(nbins):
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))
            # create colors in a gradient from red through white
            rgb = [1.0, float(j)/(nbins-1), float(j)/(nbins-1)]

            # convert rgb to hsv,  modify saturation and value and convert back to rgb
            hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
            hsv[1] = hsv[1]*sat
            hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])

    elif gradient == 'wr':
        col=[]
        coldesc=[]
        for j in range(nbins):
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))
            # create colors in a gradient from white through red
            rgb = [1.0, float(nbins-j-1)/(nbins-1), float(nbins-j-1)/(nbins-1)]

            # convert rgb to hsv,  modify saturation and value and convert back to rgb
            hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
            hsv[1] = hsv[1]*sat
            hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])

    elif gradient == 'gw':
        col=[]
        coldesc=[]
        for j in range(nbins):
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))
            # create colors in a gradient from green through white
            rgb = [float(j)/(nbins-1), 1.0, float(j)/(nbins-1)]

            # convert rgb to hsv,  modify saturation and value and convert back to rgb
            hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
            hsv[1] = hsv[1]*sat
            hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])

    elif gradient == 'wg':
        col=[]
        coldesc=[]
        for j in range(nbins):
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))
            # create colors in a gradient from white through green
            rgb = [float(nbins-j-1)/(nbins-1), 1.0, float(nbins-j-1)/(nbins-1)]

            # convert rgb to hsv,  modify saturation and value and convert back to rgb
            hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
            hsv[1] = hsv[1]*sat
            hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])

    elif gradient == 'bw':
        col=[]
        coldesc=[]
        for j in range(nbins):
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))
            # create colors in a gradient from blue through white
            rgb = [float(j)/(nbins-1), float(j)/(nbins-1), 1.0 ]

            # convert rgb to hsv,  modify saturation and value and convert back to rgb
            hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
            hsv[1] = hsv[1]*sat
            hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])

    elif gradient == 'wr':
        col=[]
        coldesc=[]
        for j in range(nbins):
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))
            # create colors in a gradient from white through blue
            rgb = [float(nbins-j-1)/(nbins-1), float(nbins-j-1)/(nbins-1), 1.0]

            # convert rgb to hsv,  modify saturation and value and convert back to rgb
            hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
            hsv[1] = hsv[1]*sat
            hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])

    elif gradient == 'gy':
        col=[]
        coldesc=[]
        for j in range(nbins):
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))
            # create colors in a gradient from green through yellow
            rgb = [float(j)/(nbins-1), 1.0, 0.]

            # convert rgb to hsv,  modify saturation and value and convert back to rgb
            hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
            hsv[1] = hsv[1]*sat
            hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])

    elif gradient == 'yg':
        col=[]
        coldesc=[]
        for j in range(nbins):
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))
            # create colors in a gradient from green through yellow
            rgb = [float(nbins-j-1)/(nbins-1), 1.0, 0.]

            # convert rgb to hsv,  modify saturation and value and convert back to rgb
            hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
            hsv[1] = hsv[1]*sat
            hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])

    elif gradient == 'bwr':
        col=[]
        coldesc=[]
        for j in range(nbins//2):
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))
            # create colors in a gradient from blue to white
            rgb = [min(1.0, float(j)*2/(nbins-1)), min(1.0,float(j)*2/(nbins-1)), min(1.0, float(nbins-j-1)*2/(nbins-1))]

            # convert rgb to hsv,  modify saturation and value and convert back to rgb
            hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
            hsv[1] = hsv[1]*sat
            hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])

        for j in range(nbins//2,nbins):
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))
            # create colors in a gradient from white to red
            rgb = [min(1.0, float(j)*2/(nbins-1)), min(1.0,float(nbins-j-1)*2/(nbins-1)), min(1.0, float(nbins-j-1)*2/(nbins-1))]

            # convert rgb to hsv,  modify saturation and value and convert back to rgb
            hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
            hsv[1] = hsv[1]*sat
            hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])

    elif gradient == 'wyr_rev':
        col=[]
        coldesc=[]
        for j in range(nbins//2):
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))
            # create colors in a gradient from red to yellow
            rgb = [1.,min(1.0, float(j)*2/(nbins-1)),0.]

            # convert rgb to hsv,  modify saturation and value and convert back to rgb
            hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
            hsv[1] = hsv[1]*sat
            hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])

        for j in range(nbins//2,nbins):
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))
            # create colors in a gradient from yellow to white
            rgb = [1.,1., min(1.0, 2./(nbins-1)*float(j)-1)]

            # convert rgb to hsv,  modify saturation and value and convert back to rgb
            hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
            hsv[1] = hsv[1]*sat
            hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])
    elif gradient == 'wyr':
        col=[]
        coldesc=[]
        for j in range(nbins//2):
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))
            # create colors in a gradient from blue to white
            rgb = [1.,1.,min(1.0, float((nbins-1)/2-j)*2/(nbins-1))]

            # convert rgb to hsv,  modify saturation and value and convert back to rgb
            hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
            hsv[1] = hsv[1]*sat
            hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])

        for j in range(nbins//2,nbins):
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))
            # create colors in a gradient from white to red
            rgb = [1.,min(1.0, float((nbins-1-j)*2./(nbins-1))),0.0]

            # convert rgb to hsv,  modify saturation and value and convert back to rgb
            hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
            hsv[1] = hsv[1]*sat
            hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])

    elif gradient == 'wcb':
        col=[]
        coldesc=[]
        for j in range(nbins//2):
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))
            # create colors in a gradient from blue to white
            rgb = [min(1.0, float((nbins-1)/2-j)*2/(nbins-1)),1.0,1.0]

            # convert rgb to hsv,  modify saturation and value and convert back to rgb
            hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
            hsv[1] = hsv[1]*sat
            hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])

        for j in range(nbins//2,nbins):
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))
            # create colors in a gradient from white to red
            rgb = [0.,min(1.0, float((nbins-1-j)*2./(nbins-1))),1.0]

            # convert rgb to hsv,  modify saturation and value and convert back to rgb
            hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
            hsv[1] = hsv[1]*sat
            hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])

    elif gradient == 'wc':
        col=[]
        coldesc=[]
        for j in range(nbins):
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))
            # create colors in a gradient from white to red
            rgb = [min(1.0, float((nbins-1)/2-j)*2/(nbins-1)),1.0,1.0]

            # convert rgb to hsv,  modify saturation and value and convert back to rgb
            hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
            hsv[1] = hsv[1]*sat
            hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])

    elif gradient == 'ryw':
        col=[]
        coldesc=[]
        for j in range(nbins//2):
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))
            # create colors in a gradient from blue to white
#      rgb = [1.,1.,min(1.0, float((nbins-1)/2-j)*2/(nbins-1))]
            rgb = [1.,min(1.0, float(j)*2/(nbins-1)),0.0]
            # convert rgb to hsv,  modify saturation and value and convert back to rgb
            hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
            hsv[1] = hsv[1]*sat
            hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])

        for j in range(nbins//2,nbins):
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))
            # create colors in a gradient from white to red
            rgb = [1.,1.,min(1.0, float(j)/(nbins-1))]

            # convert rgb to hsv,  modify saturation and value and convert back to rgb
            hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
            hsv[1] = hsv[1]*sat
            hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])

    elif gradient == 'rwb':
        col=[]
        coldesc=[]
        for j in range(nbins//2):
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))
            # create colors in a gradient from red to white
            rgb = [min(1.0, float(nbins-j-1)*2/(nbins-1)), min(1.0,float(j)*2/(nbins-1)), min(1.0, float(j)*2/(nbins-1))]

            # convert rgb to hsv,  modify saturation and value and convert back to rgb
            hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
            hsv[1] = hsv[1]*sat
            hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])

        for j in range(nbins//2,nbins):
            coldesc.append('col' + item + gradient + str(j))
            # coldesc.append('col' + str(sel[j]) + str(j))
            # create colors in a gradient from white to blue
            rgb = [min(1.0, float(nbins-j-1)*2/(nbins-1)), min(1.0,float(nbins-j-1)*2/(nbins-1)), min(1.0, float(j)*2/(nbins-1))]

            # convert rgb to hsv,  modify saturation and value and convert back to rgb
            hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
            hsv[1] = hsv[1]*sat
            hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])

    elif gradient == 'gray' or gradient == 'grey':
        col=[]
        coldesc=[]
        for j in range(nbins):
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))
            # create colors in a gradient of grays from "sat" to "value"

            hsv = [0, 0, sat + (value-sat)*float(j)/(nbins-1)]
#      hsv[1] = hsv[1]*sat
#      hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])

    elif gradient == 'reversegray' or gradient == 'reversegrey':
        col=[]
        coldesc=[]
        for j in range(nbins):
            # coldesc.append('col' + str(sel[j]) + str(j))
            coldesc.append('col' + item + gradient + str(j))
            # create colors in a gradient of grays from "sat" to "value"

            hsv = [0, 0, value - (value-sat)*float(j)/(nbins-1)]
#      hsv[1] = hsv[1]*sat
#      hsv[2] = hsv[2]*value
            rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color("col" + item + gradient + str(j),col[j])


# return the gradient as a list of colors named by their index (i.e. col0,col1,col2,col3,...)
    return coldesc

# allow calling without parentheses: color_hist_b [selection=], [mode= ],[gradient= ],[nbins= ]
cmd.extend("color_b",color_b)
cmd.auto_arg[0]['color_b'] = [lambda: cmd.Shortcut(cmd.get_names('all')),'object name', '']
cmd.extend("color_q",color_q)
cmd.auto_arg[0]['color_q'] = [lambda: cmd.Shortcut(cmd.get_names('all')),'object name', '']



def help_utils():
    """ all functions help messages in utils.py
    """
    print("# 1. load all files in a folder")
    print("load_all_in_dir(dir)")
    print()
    print("# 2. load all PDB files in a folder")
    print("load_pdb_in_dir(dir)")
    print()
    print("# 3. align all objects onto a selection (def='sele') or align all other objects onto an object")
    print("align_all(selection='sele')")
    print()
    print("# 4. color a selection by Consurf coloring method")
    print("colour_consurf(selection='all', keep_selections=0)")
    print()
    print("# 5. show interface contact for 2 chains for a complex (or all complexes);")
    print("#    you can select 2 chains and then do analyse_complex sele; you can give a name")
    print("#    for your interface so that it will not be overwritten when you analyze the second interface")
    print("analyze_complex(object='all', name=None)")
    print()
    print("# 6. save a molecule (or a selection, or use 'all' for all objects but not any selection)")
    print("#    in the current folder. The file name will be pymol_xxx.pdb.")
    print("save_mol(object='all')")
    print()

cmd.extend("help_utils", help_utils)
