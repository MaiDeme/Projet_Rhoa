cmd.do("run scripts/utils.py")

cmd.do("load_top_rank_pdb_in_dir_with_pml CMAPnb_models, 10, CoevMapNbCst, CMAPnb_models/VisualizeContacts_interhetero_Complex_CMAPnb")
cmd.do("refresh")
cmd.do("load_top_rank_pdb_in_dir_with_pml CMAPscore_models, 10, CoevMapScoreCst, CMAPscore_models/VisualizeContacts_interhetero_Complex_CMAPscore")
cmd.do("refresh")
cmd.do("disable C*")

cmd.do("hide everything")
cmd.do("refresh")
cmd.do("show cartoon, *")
cmd.do("refresh")
cmd.do("show dash, *")
util.color_chains("all")
cmd.do("refresh")

cmd.do("enable CoevMapNbCst")
cmd.do("zoom Complex_CMAPnb1")
cmd.do("enable Complex_CMAPnb1")
cmd.do("refresh")
cmd.do("enable c_*nb1")
cmd.do("enable c*_*_*_*")
cmd.do("refresh")

print("\n\nINFO: \nFor every PDB model 'xxx', predicted contacts satisfied in the model (<8A) are indicated as dashed lines which can be enabled/disabled as objects in the group below labeled 'c_xxx' below each complex name.\n\n Two groups of models are reported CoevMapNbCst and CoevMapScoreCst (10 models in each):\n  - In group CoevMapNbCst: Models were sorted with respect to the number of non redundant constraints/contacts they satisfied.\n  - In group CoevMapScoreCst: Models were sorted with respect to the sum of the probabilities of all the non redundant constraints/contacts which were satisfied.\n\nThe format of contact names is :\n c<index_of_contact>_<index_of_line_in_coevmap_of_the_most_probable_contact_in_group_of_redundant_contacts>_<individual_contact_probability>_<highest_probability_for_group_of_redundant_contacts_used_in_scoring>\n ")