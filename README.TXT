To make plots (e.g. MCS bias/resolution plots), use MCS_test.py in the MCS_Technote directory. Also in this directory is a script to rerun the MCS algorithm: ELOSS_redo-mcs_mcc7_over2GeV_3mrad_2.py.

The other directories have various xmls for creating an MCS ntuple over the grid. The main one to use is probably in ntupleOnGridNewVer, which takes a SAM defined file. The others are examples of creating samples over the grid from scratch.

For space charge studies, I made two separate samples and matched them event by event for comparison. The script to do this is eventRemover_noSC.py (or withSC). To make a new TTREE with MCS info from both files, use the SCE_diff_TTREE_maker.py. 

An example of the old code I used to make plots from summer 2016 is diff_analyze1.py. 