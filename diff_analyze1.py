# Python script to create graphs of: 
# (1) MCS p w/o SCE vs. MCS p with SCE       
# (2) Entries vs. delta for each E range
# (3) Delta mean vs. truth E
# (4) Delta std. dev. vs. truth E

import ROOT
from ROOT import TFile, gDirectory, TCanvas, TH1D, TColor, TGraphErrors, TMath, TLatex
from ROOT import TH2D
from operator import truediv

# Open the files 
f = TFile("mooneySC_May22_10k/testmultiscattermomentum_SCE_comparison_TTREE_withTruth.root","READ")

# See what is inside the file
print "Now looking inside the file:"
f.ls()

# See what variables are stored inside the tree by looking at just one entry
print "Showing one entry of the tree (note the variables):"
f.MCS_SCE_comparison.Show(0)

# Get TTree
print "Gettin' TTree"
t = f.Get("MCS_SCE_comparison")

print "Got TTree..."

# Making graph (1):
# With TTree.Draw, first argument is the variable we are plotting, second
# argument is cuts we want to make, third argument is draw styles
print "This is the current graph of MC p (no SC) vs. MCS p (w/ SC)."
noSCVswSCMcsRecoColz = TCanvas("noSCVswSCMcsRecoColz", "noSCVswSCMcsRecoColz", 720, 152, 682, 505)
noSCVswSCMcsRecoColz.cd()
# Cut is placed to exclude failed events (-1)
corr = TH2D("corr","Muon MCS Momentum Correlation: with SC vs. no SC",100,0,2,100,0,2)
f.MCS_SCE_comparison.Draw("mcs_mom_withSC:mcs_mom_noSC>>corr", "mcs_mom_noSC>0 && mcs_mom_withSC>0", "COLZ")
#&& true_mom==t2.true_mom
corr.GetXaxis().SetTitle("MCS Reco Momentum, no SC [GeV]")
corr.GetYaxis().SetTitle("MCS Reco Momentum, with SC [GeV]")

input = raw_input("Press enter to continue...")

print "This is the difference plot:"
diffPlot = TCanvas("diffPlot","diffPlot",720,152,682,505)
diffPlot.cd()
corr2 = TH2D("corr2","Fractional Difference vs. MCS Reco Momentum, no SC",100,0,2,100,0,2)
f.MCS_SCE_comparison.Draw("(mcs_mom_withSC-mcs_mom_noSC)/mcs_mom_noSC:mcs_mom_noSC","mcs_mom_noSC>0 && mcs_mom_withSC>0","COLZ")
corr2.GetXaxis().SetTitle("MCS Reco Momentum, no SC [GeV]")
corr2.GetYaxis().SetTitle("(MCS p (w/ SC) - MCS p (no SC)) / MCS p (no SC)")

input = raw_input("Press enter to continue...")

print "This is a 1D histogram of fractonal difference (MCS and truth), no SC."
frNoSC = TCanvas("frNoSC","frNoSC",720,152,682,505)
frNoSC.cd()
fr1 = TH1D("fr1","Fractional Difference, no SC",100,0,2)
f.MCS_SCE_comparison.Draw("(mcs_mom_noSC-truth_mom_noSC)/truth_mom_noSC","mcs_mom_noSC>0","COLZ")

input = raw_input("Press enter to continue...")

print "This is a 1D histogram of fractonal difference (MCS and truth), with SC."
frWithSC = TCanvas("frWithSC","frWithSC",720,152,682,505)
frWithSC.cd()
fr2 = TH1D("fr2","Fractional Difference, with SC",100,0,2)
f.MCS_SCE_comparison.Draw("(mcs_mom_withSC-truth_mom_withSC)/truth_mom_withSC","mcs_mom_withSC>0","COLZ")

input = raw_input("Press enter to continue...")

###-----------------

# Function to create graphs (2), (3), (4):
def getHists(contained_only = False):
    
    # Define the desired energy ranges
    ERanges = [(0.3,0.5), (0.5,0.7), (0.7,0.9), (0.9,1.1), (1.1,1.3), (1.3,1.5), (1.5,1.7), (1.7,1.9), (1.9,2.1)]
    
    # Creating empty histograms
    hist_dictionary = {}
    for myE in ERanges:
        myName = "h_min%s_max%s" % (myE[0], myE[1])
        if contained_only: myName = "cont_h_min%s_max%s" % (myE[0], myE[1])
        myTitle = "Energies from %0.1f to %0.1f;(MCS p (withSC) - MCS p (no SC)) / MCS p (no SC);Entries" % (myE[0],  myE[1])
        if contained_only: myTitle = "Contained Tracks: Energies from %0.1f to %0.1f;(MCS p (withSC) - MCS p (no SC)) / MCS p (no SC);Entries" % (myE[0], myE[1])
        hist_dictionary[myE] = TH1D(myName, myTitle, 100, -4, 4)
        
    # Loop over energy ranges to create graphs of type (2)
    for myE in ERanges:
        myPlot = "(mcs_mom_withSC - mcs_mom_noSC) / mcs_mom_noSC >> h_min%0.1f_max%0.1f" % (myE[0], myE[1])
        if contained_only: myPlot = "(mcs_mom_withSC - mcs_mom_noSC) / mcs_mom_noSC >> cont_h_min%0.1f_max%0.1f" % (myE[0], myE[1])
        myCut = "mcs_mom_withSC>0 && mcs_mom_noSC>0 && mcs_mom_noSC > %f && mcs_mom_noSC < %f" % (myE[0], myE[1])
        #if contained_only: myCut += " && mu_contained == 1"
        f.MCS_SCE_comparison.Draw(myPlot, myCut)
        
    # Create lists for graphs (3), (4) that will be filled in the next loop
    xpoints = []
    filMeanVals, filSqrtHistEntries = [], []
    filStdVals, filStdErrs = [], []
    
    for myE in ERanges: 
        # Filter out energy ranges with 0 entries
        if hist_dictionary[myE].GetEntries() != 0:
            xpoints.append(myE[0] + ( myE[1] - myE[0] ) / 2.)
            filMeanVals.append(hist_dictionary[myE].GetMean())
            filSqrtHistEntries.append(TMath.Sqrt(hist_dictionary[myE].GetEntries()))
            filStdVals.append(hist_dictionary[myE].GetRMS())
            filStdErrs.append(hist_dictionary[myE].GetRMSError())
        
    # Calculate error bars for graph (3): (stdDev / sqrt(entries))
    yerrs = map(truediv, filStdVals, filSqrtHistEntries)

    # Define a constant width (x-error bar) for each marker
    xerrs = [( myE[1] - myE[0] ) / 2. for myE in ERanges]

    # Graph (3) and (4) setup
    meanGraph = TGraphErrors()
    stdGraph = TGraphErrors()

    meanGraphTitle = "Mean vs. Energy;MCS Muon Energy, no SC [GeV];#left[p_{MCSwithSC}-p_{MCSnoSC} / p_{MCSnoSC}#right]"
    if contained_only: meanGraphTitle = "Contained Tracks: Mean vs. Energy;MCS Muon Energy, no SC [GeV];#left[p_{MCSwithSC}-p_{MCSnoSC} / p_{MCSnoSC}#right]"
    stdGraphTitle = "Standard Deviation vs. Energy;MCS Muon Energy, no SC [GeV];Fractional Momentum Resolution"
    if contained_only: stdGraphTitle = "Contained Tracks: Standard Deviation vs. Energy;MCS Muon Energy, no SC [GeV];Fractional Momentum Resolution"

    meanGraph.SetTitle(meanGraphTitle)
    stdGraph.SetTitle(stdGraphTitle)

    # Fill in data for graph (3)
    for i in xrange(len(xpoints)):
        meanGraph.SetPoint(i, xpoints[i], filMeanVals[i])
        meanGraph.SetPointError(i, xerrs[i], yerrs[i])
    meanGraph.Draw("ALP")

    # Fill in data for graph (4)
    for i in xrange(len(xpoints)):
        stdGraph.SetPoint(i, xpoints[i], filStdVals[i])
        stdGraph.SetPointError(i, xerrs[i], filStdErrs[i])
    stdGraph.Draw("ALP")

    return hist_dictionary, meanGraph, stdGraph

myDict, meanGraph, stdGraph = getHists()
myDictCont, meanGraphCont, stdGraphCont = getHists(contained_only=True)

print "Now writing the graphs..."
    
fout = TFile("mooneySC_May22_10k/testmultiscattermomentum_SCE_comparison_SPREAD_2.root", "RECREATE")
fout.cd()

for myE, myhist in myDict.iteritems():
    myhist.Write()
meanGraph.Write()
stdGraph.Write()

#for myE, myhist in myDictCont.iteritems():
#    myhist.Write()
#meanGraphCont.Write()
#stdGraphCont.Write()

fout.Close()
