# Space Charge Studies (comparing events with and without SCE):
#
# This script creates a new .root file with the same TTree structure as the
# original, but with only events that we want to keep for comparison for withSC.
# Note: We filter out by truth momentum since we don't have eventID
#
# Edited for use with Giuseppe's ntuples

from ROOT import TFile, TTree

# Open original file
f = TFile("mooneySC_mcc8_1_July26/test_muminus_0.1-2.0GeV_isotropic_uboone_10kEvents_withSC.root","READ")
# Open other file we want to compare to
fno = TFile("mooneySC_mcc8_1_July26/test_muminus_0.1-2.0GeV_isotropic_uboone_10kEvents_noSC.root","READ")

print "Opened files..."

f.cd("trajmcsntupleCosmic")
fno.cd("trajmcsntupleCosmic")

print "Changed directories to 'Cosmic.'"

# Get TTree for each sample
t = f.Get('/trajmcsntupleCosmic/tree')
tno = fno.Get('/trajmcsntupleCosmic/tree')

print "Got TTrees..."

# Create new .root file that will contain the result of this script
newF = TFile("mooneySC_mcc8_1_July26/test_muminus_0.1-2.0GeV_isotropic_uboone_10kEvents_withSC_FILTERED_precision1.root","recreate")
newT = TTree("name_of_tree", "tree")
newT = t.CloneTree(0)

print "Cloned tree structure..."

# Get total number of entries in TTrees
entries = t.GetEntries()
print "Original file's entries (withSC):",entries

entries_no = tno.GetEntries()
print "Other file's entries (noSC):",entries_no

# Vectors to store info for new version of original file
true_mom_v = []
reco_mom_v = []

newF.cd()

# Loop through entries in original file's TTree
# If event found in both samples -> add values to arrays
for e in xrange(entries):

    # Grab the entry in this TTree    
    t.GetEntry(e)

    # Variable to keep track of whether event has been matched
    matched = False

    # Loop through all entries of the second TTree and compare event-by-event
    for i in xrange(entries_no):
        # Grab entry from second TTree                                
        tno.GetEntry(i)
        # Compare the two. Are the true momenta the same?
#        if (t.simMom == tno.simMom):
        if ((abs(t.simMom - tno.simMom) <= 0.00000000000000000001) and (t.simMom != -999.000000) and (tno.simMom != -999.000000)):
            print 'Found match: entry %i with true mom = %f matches entry %i with true mom = %f'%(e,t.simMom,i,tno.simMom)
            matched = True
            newT.Fill()
            break

    if (matched == True):

        true_mom_v.append(t.simMom)
        reco_mom_v.append(t.trkMom_Mu)

print 'A total of %i entries were matched over %i total entries in the TTree'%(len(true_mom_v),entries)

newF.Write()
newF.Close()
