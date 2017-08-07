# Makes TTree of variables needed for fractional difference 
# plot of MCS with and without SCE
# (Putting noSC and withSC MCS p's into same ntuple)
#
# Note: We filter out by truth momentum since we don't have eventID
# in original TTrees
#
# Edited for use with Giuseppe's ntuples

from ROOT import TFile, TTree
from array import array

# Open separate noSC and withSC files
f = TFile("mooneySC_mcc8_1_July26/test_muminus_0.1-2.0GeV_isotropic_uboone_10kEvents_noSC_FILTERED.root","READ")
fwith = TFile("mooneySC_mcc8_1_July26/test_muminus_0.1-2.0GeV_isotropic_uboone_10kEvents_withSC_FILTERED.root","READ")

print "Opened files..."

#f.cd("trajmcsntupleCosmic")
#fwith.cd("trajmcsntupleCosmic")

#print "Changed directories to 'Cosmic.'"

# Get TTree from each file
t = f.Get('tree')
twith = fwith.Get('tree')

print "Got TTrees..."

# Create new .root file that will contain the result of this script
newF = TFile("mooneySC_mcc8_1_July26/test_muminus_0.1-2.0GeV_isotropic_uboone_10kEvents_SCE_comparison_TTREE.root","recreate")
# Make this file have a TTree
newT = TTree("MCS_SCE_comparison", "MCS SCE Comparison Tree")

# Add variables to TTree (in Python they must be arrays to mimic C++ types)
mcs_mom_noSC = array('d', [0.])
mcs_mom_withSC = array('d', [0.])
truth_mom_noSC = array('d', [0.])
truth_mom_withSC = array('d', [0.])
contained_noSC = array('i', [0])
contained_withSC = array('i', [0])
newT.Branch('mcs_mom_noSC', mcs_mom_noSC, 'mcs_mom_noSC/D')
newT.Branch('mcs_mom_withSC', mcs_mom_withSC, 'mcs_mom_withSC/D')
newT.Branch('truth_mom_noSC', truth_mom_noSC, 'truth_mom_noSC/D')
newT.Branch('truth_mom_withSC', truth_mom_withSC, 'truth_mom_withSC/D')
newT.Branch('contained_noSC', contained_noSC, 'contained_noSC/I')
newT.Branch('contained_withSC', contained_withSC, 'contained_withSC/I')

# Get total number of entries in TTrees
entries = t.GetEntries()
print "Original file's entries (noSC):", entries
entries_with = twith.GetEntries()
print "Other file's entries (withSC):", entries_with

newF.cd()

# Loop through entries in original file's TTree
# If event found in both samples -> add values to arrays
for e in xrange(entries):

    # Grab the entry in this TTree    
    t.GetEntry(e)

    # Loop through all entries of the second TTree and compare event-by-event
    for i in xrange(entries_with):
        # Grab entry from second TTree                                
        twith.GetEntry(i)
        # Compare the two. Are the true momenta the same?
        #if (t.simMom == twith.simMom):
        if (abs(t.simMom - twith.simMom) <= 0.00000000000000000001):
#            print 'Found match: entry %i with true mom = %f matches entry %i with true mom = %f'%(e,t.true_mom,i,twith.true_mom)
            mcs_mom_noSC[0] = t.trkMom_Mu
            mcs_mom_withSC[0] = twith.trkMom_Mu
            truth_mom_noSC[0] = t.simMom
            truth_mom_withSC[0] = twith.simMom
            contained_noSC[0] = t.trkIsContained
            contained_withSC[0] = twith.trkIsContained
            newT.Fill()
            break

newF.Write()
newF.Close()

print "Finished!"
