import os
if not os.path.exists("scripts/_analysis_top30.pml"): print("WARNING: If you are a Linux user, please first do cd <folder_results> !!!");
if os.path.exists("scripts/_analysis_top30.pml"): cmd.do("run scripts/_analysis_top30.pml");
