##############################################################################
# BEFORE RUNNING CODE BE SURE TO CHANGE TO COSMIC OR NU IN SHARED_CONTENT.PY
# AND EDIT FILE NAMES ACCORDINGLY BELOW!
##############################################################################

#get_ipython().magic(u'matplotlib inline')

from shared_content import *

# Do you want to write figures to png files? Careful, it will overwrite!
write_figures = True
figdir = 'Figures/'

filedir = 'anafiles/'
#'SingleMuonRecoTrack'
#'MCBNBSelectedRecoTrack'
#'DataBNBSelectedRecoTrack'
#'MCBNBRecoTrack'
#'MCBNBMCTrack'
#'SingleMuonMCTrack'
#'MCBNBMCTrackExiting'
#'MCBNBRecoTrackExiting'
anatype = 'MCBNBRecoTrack'  

seglen = 14

myfile = 'mcc7_bnb_nuOnly.root' # 3mrad, full E spectrum, mcc7 sample
#myfile = 'QUICKTEST_mcc7_bnb_nuOnly_over2GeV_1mrad.root'
#myfile = 'REMOVE_SEGS_mcc7_bnb_nuOnly_contained_3mrad_eloss0_removeLast2Segs.root'
#myfile = 'QUICKTEST_mcc7_bnb_nuOnly_over2GeV_3mrad.root'
#myfile = 'MCSBiasStudy_SingleMuonRecoTrack_anaout_14cmseg_muplus.root'
#myfile = 'MCSBiasStudy_%s_anaout_14cmseg_3res_bothscatters_nonrelfix.root' %( anatype)
#myfile = 'MCSBiasStudy_%s_anaout_%dcmseg_2res_bothscatters_nonrelfix.root' %( anatype, seglen )
#myfile = 'MCSBiasStudy_%s_anaout_14cmseg_3res_bothscatters_nonrelfix_realdedx1_highlandconstantMOMENTUMDEPENDENT.root'%( anatype )
print myfile
#get_ipython().system(u'ls -ltr $filedir | grep $myfile')

df = get_dfs(filedir + myfile)
print len(df)
print df.columns.values

# In[37]:

# myrun = 4
# mysub = 5093
# myevt = 101842
# print df.\
# query('full_MCS_energy<1.5 and true_E > 2.5')\
# .query('run == %d and subrun == %d and eventid == %d' % (myrun,mysub,myevt))\
# [['full_MCS_energy','true_E','run','subrun','eventid','full_range_energy']]#\
# #.head(n=1)
# print segdf.query('run == %d and subrun == %d and eventid == %d' % (myrun,mysub,myevt))\
# [['delta_theta_x','delta_theta_y','dthetayoverpredictedRMS_fromMCS','dthetaxoverpredictedRMS_fromMCS']]


# In[38]:

# Basic plot of muon energy and angle spectra ala Bruce's request
if anatype == 'MCBNBMCTrack':
    # Energy spectra
    myx = df['true_E'].values
    plt.figure(figsize=(10,6))
    plt.hist(myx,bins=np.linspace(0,2.2,100),alpha=0.5)
    plt.grid(True,'both')
    plt.title('True Muons Passing Selection Cuts',fontsize=16)
    plt.xlabel('True Total Energy [GeV]',fontsize=14)
    plt.ylabel('Events')
    plt.yscale('log')
    if write_figures == True:
        plotname = 'MCBNBMCTrack_EnergySpectrum.png'
        print "Saving figure %s"%plotname
        plt.tight_layout()
        plt.savefig(figdir + plotname)
    
    # Angle spectra
    myx = df['theta'].values
    plt.figure(figsize=(10,6))
    plt.hist(myx,bins=np.linspace(0,3.14159,100),alpha=0.5)
    plt.grid(True,'both')
    plt.title('True Muons Passing Selection Cuts',fontsize=16)
    plt.xlabel('True Theta Angle w.r.t. Beam [rad]',fontsize=14)
    plt.ylabel('Events')
    if write_figures == True:
        plotname = 'MCBNBMCTrack_AngleSpectrum.png'
        print "Saving figure %s"%plotname
        plt.tight_layout()
        plt.savefig(figdir + plotname)


# In[39]:

def basic_comparison_fig(xvar, yvar, plotname=None, extraquery = None, addtext = None, nbins=50, ybinmin = 0, ybinmax = 6, xbinmin = 0, xbinmax = 6):
    plt.figure(figsize=(10,6)) 
    myquery = 'trkMom_Mu < 99999999 and trkMom_P < 99999999'
#    myquery = 'mcs_mom_noSC < 99999999 and mcs_mom_withSC < 99999999'
    if extraquery is not None: myquery += ' and %s'%extraquery
    myx = df.query(myquery)[xvar].values
    myy = df.query(myquery)[yvar].values
    blah = plt.hist2d(myx,myy,bins=((np.linspace(xbinmin,xbinmax,nbins),np.linspace(ybinmin,ybinmax,nbins))),cmin=1)
    blah = plt.colorbar()
    blah = plt.grid(True)
    blha = plt.xlabel('%s'%titles[xvar],fontsize=16)
    blha = plt.ylabel('%s'%titles[yvar],fontsize=16)
    blha = plt.title('%s'%titles[anatype],fontsize=16)
    blha.set_y(1.04)
    blah = plt.plot([0,100],[0,100],'g--',linewidth=3)
    
    if addtext is not None:
        plt.text(plt.xlim()[1]*0.55, plt.ylim()[1]*0.1, addtext, fontsize=20)
        
    if write_figures:
        if plotname is not None: 
            print "Saving figure %s"%plotname
            plt.tight_layout()
            plt.savefig(figdir + plotname)
        else: print "YOU WANTED TO SAVE A PLOT BUT DIDN'T GIVE A PLOT NAME!"


# In[40]:

# MCS momentum vs range momentum for all samples

# Muon
#basic_comparison_fig(xvar='trkMom_RangeMu',yvar='trkMom_Mu',
#                         plotname='giuseppeNuTest_noperfcut_%s.png'%anatype,
#                        extraquery='trkIsContained and simID==13 and trkLength>100', nbins = 100)
#basic_comparison_fig(xvar='trkMom_RangeMu',yvar='trkMom_Mu',
#                         plotname='giuseppeMuplusCosmicTest_%s.png'%anatype,
#                        extraquery='trkIsContained and simID==13 and trkLength>100 and (trkLength>=(0.9*simLength) and trkLength<(1.1*simLength))', nbins = 100)

#basic_comparison_fig(xvar='mcs_mom_noSC',yvar='mcs_mom_withSC',
#                         plotname='SCE_withVsNo.png',
#                        extraquery='mcs_mom_noSC>0 and mcs_mom_withSC>0', nbins = 100)

#df["newvar"] = df.apply(lambda x : (x['var1']-x['var2']) / x['var1'], axis=1)

#basic_comparison_fig(xvar='mcs_mom_noSC',yvar='fracDiff',
#                         plotname='SCE_fracDiffVsNoSC.png',
#                        extraquery='mcs_mom_noSC>0 and mcs_mom_withSC>0', nbins = 100)

#basic_comparison_fig(xvar='simMom',yvar='trkMom_Mu',
#                         plotname='ELOSS_2_mcc7_exiting_nu_vsTruth.png',
#                        extraquery='trkIsContained==False and simID==13 and trkLength>100 and trkMom_Mu<7', nbins = 100, binmax=4)

basic_comparison_fig(xvar='simMom',yvar='trkMom_Mu',
                         plotname='Aug04_test_vsTruth.png',
                        extraquery='trkIsContained==False and simID==13 and trkLength>100 and trkMom_Mu<7', nbins = 100)

basic_comparison_fig(xvar='simMom',yvar='trkMom_Mu',
                         plotname='Aug04_test2_vsTruth.png',
                        extraquery='trkIsContained==False and simID==13 and trkLength>100 and trkMom_Mu<7 and (trkLength>=(0.9*simLength) and trkLength<(1.1*simLength))', nbins = 100)
"""
# Want to compare fractional difference (MCS wrt Truth) vs. trk Length
basic_comparison_fig(xvar='trkLength',yvar='fracDiff',
                         plotname='Comparison_mcc7_nu_slice_1.97_2.38_scatterPlot.png',
                        extraquery='trkIsContained==False and simID==13 and trkLength>100 and trkMom_Mu<7 and simMom < 2.38 and simMom > 1.97', nbins = 100, ybinmin = -1, ybinmax = 8, xbinmin = 0,xbinmax=1000)

basic_comparison_fig(xvar='trkLength',yvar='fracDiff',
                         plotname='Comparison_mcc7_nu_slice_2.38_2.78_scatterPlot.png',
                        extraquery='trkIsContained==False and simID==13 and trkLength>100 and trkMom_Mu<7 and simMom < 2.78 and simMom > 2.38', nbins = 100, ybinmin = -1, ybinmax = 8, xbinmin = 0,xbinmax=1000)

basic_comparison_fig(xvar='trkLength',yvar='fracDiff',
                         plotname='Comparison_mcc7_nu_slice_1.57_1.97_scatterPlot.png',
                        extraquery='trkIsContained==False and simID==13 and trkLength>100 and trkMom_Mu<7 and simMom < 1.97 and simMom > 1.57', nbins = 100, ybinmin = -1, ybinmax = 8, xbinmin = 0,xbinmax=1000)

# Comparing fracDiff vs. # hits in track
basic_comparison_fig(xvar='trkNHits',yvar='fracDiff',
                         plotname='Comparison_trkNHits_mcc7_nu_slice_1.97_2.38_scatterPlot.png',
                        extraquery='trkIsContained==False and simID==13 and trkLength>100 and trkMom_Mu<7 and simMom < 2.38 and simMom > 1.97', nbins = 100, ybinmin = -1, ybinmax = 8, xbinmin = 0,xbinmax=7000)

basic_comparison_fig(xvar='trkNHits',yvar='fracDiff',
                         plotname='Comparison_trkNHits_mcc7_nu_slice_2.38_2.78_scatterPlot.png',
                        extraquery='trkIsContained==False and simID==13 and trkLength>100 and trkMom_Mu<7 and simMom < 2.78 and simMom > 2.38', nbins = 100, ybinmin = -1, ybinmax = 8, xbinmin = 0,xbinmax=7000)

basic_comparison_fig(xvar='trkNHits',yvar='fracDiff',
                         plotname='Comparison_trkNHits_mcc7_nu_slice_1.57_1.97_scatterPlot.png',
                        extraquery='trkIsContained==False and simID==13 and trkLength>100 and trkMom_Mu<7 and simMom < 1.97 and simMom > 1.57', nbins = 100, ybinmin = -1, ybinmax = 8, xbinmin = 0,xbinmax=7000)

# Comparing fracDiff vs. Shortest Seg Length in Track

#def findMin(vectorOfTrkLengths):
#    min = 0
#    for i in vectorOfTrkLengths: 
#        if i < min:
#            min = i
#    return min

#minSegLen = findMin(trkSegRadLengths)

basic_comparison_fig(xvar='trkSegRadLengthsMin',yvar='fracDiff',
                         plotname='Comparison_segLength_mcc7_nu_slice_1.97_2.38_scatterPlot.png',
                        extraquery='trkIsContained==False and simID==13 and trkLength>100 and trkMom_Mu<7 and simMom < 2.38 and simMom > 1.97', nbins = 100, ybinmin = -1, ybinmax = 8, xbinmin = 0,xbinmax=1.1)

basic_comparison_fig(xvar='trkSegRadLengthsMin',yvar='fracDiff',
                         plotname='Comparison_segLength_mcc7_nu_slice_2.38_2.78_scatterPlot.png',
                        extraquery='trkIsContained==False and simID==13 and trkLength>100 and trkMom_Mu<7 and simMom < 2.78 and simMom > 2.38', nbins = 100, ybinmin = -1, ybinmax = 8, xbinmin = 0,xbinmax=1.1)

#basic_comparison_fig(xvar='trkSegRadLengthsMin',yvar='fracDiff',
#                         plotname='Comparison_segLength_mcc7_nu_slice_1.97_2.38_scatterPlot.png',
#                        extraquery='trkIsContained==False and simID==13 and trkLength>100 and trkMom_Mu<7 and simMom < 2. and simMom > 2.78', nbins = 100, ybinmin = -1, ybinmax = 8, xbinmin = 0,xbinmax=1.1)

basic_comparison_fig(xvar='trkSegRadLengthsAvg',yvar='fracDiff',
                         plotname='Comparison_AvgSegLength_mcc7_nu_slice_1.97_2.38_scatterPlot.png',
                        extraquery='trkIsContained==False and simID==13 and trkLength>100 and trkMom_Mu<7 and simMom < 2.38 and simMom > 1.97', nbins = 100, ybinmin = -1, ybinmax = 8, xbinmin = 0.99,xbinmax=1.1)

basic_comparison_fig(xvar='trkSegRadLengthsAvg',yvar='fracDiff',
                         plotname='Comparison_AvgSegLength_mcc7_nu_slice_2.38_2.78_scatterPlot.png',
                        extraquery='trkIsContained==False and simID==13 and trkLength>100 and trkMom_Mu<7 and simMom < 2.78 and simMom > 2.38', nbins = 100, ybinmin = -1, ybinmax = 8, xbinmin = 0.99,xbinmax=1.1)
"""
"""
# Proton
basic_comparison_fig(xvar='trkMom_RangeP',yvar='trkMom_P',
                         plotname='giuseppeNuTest_P_%s.png'%anatype,
#                         extraquery='trkid2==2212', nbins = 100)
                        extraquery='trkIsContained and (trkLength>=(0.9*simLength) and trkLength<(1.1*simLength)) and simID==13 and trkLength>100', nbins = 100)


# MCS Muon momentum vs true momentum for all samples (except data)
if 'Data' not in anatype:
#    basic_comparison_fig(xvar='simMom',yvar='trkMom_Mu',
#                         plotname='giuseppeNu_MCS_true_comparison_noperfcut_%s.png'%anatype,
#                        extraquery='trkIsContained and simID==13 and trkLength>100', nbins = 100)
    basic_comparison_fig(xvar='simMom',yvar='trkMom_Mu',
                         plotname='giuseppeMuplusCosmic_MCS_true_comparison_%s.png'%anatype,
                        extraquery='trkIsContained and simID==13 and trkLength>100 and (trkLength>=(0.9*simLength) and trkLength<(1.1*simLength))', nbins = 100)


# In[41]:

# Range energy vs True Energy for single muon MCTrack analysis only
#if anatype == 'SingleMuonMCTrack':
if anatype == 'MCBNBMCTrack':
    basic_comparison_fig(xvar='true_E',yvar='full_range_energy',
                         plotname='true_range_comparison_MCTracks.png',
                        nbins=100)


# MCS momentum vs True momentum for EXITING MCTRACK ANALYSIS ONLY
if 'Exiting' in anatype:
    
    basic_comparison_fig(xvar='true_momentum',yvar='full_MCS_momentum',
                         plotname='MCS_true_comparison_%s.png'%anatype,
                        nbins=100,binmin=0,binmax=4)
    
    
    extraqueries = ['full_length < 150.',
                   'full_length > 150. and full_length < 250.',
                   'full_length > 250. and full_length < 350.',
                   'full_length > 350.']
    for extraquery in extraqueries:
        basic_comparison_fig(xvar='true_momentum',yvar='full_MCS_momentum',
                            nbins=100,binmin=0,binmax=4,extraquery=extraquery)


# In[43]:

#Breakdown of proton, pion MIDs and true muon IDs
if anatype == 'MCBNBSelectedRecoTrack':

    for (pdg,name) in [ (211,'Pion MIDs'), (13,'True Muons'), (2212, 'Proton MIDs')]:
        
        #Testing this way
        extraquery = 'MCT_PDG == %d or MCT_PDG == -%d'%(pdg,pdg)
        basic_comparison_fig(xvar='full_range_momentum',yvar='full_MCS_momentum',
                         plotname='MCS_range_momentum_MCBNBSelectedRecoTrack_%d.png' % pdg,
                            extraquery = extraquery,
                            addtext = name)
    # Pie chart!

    plt.figure(figsize=(10,6))
    
    sizes = [ len(df.query('MCT_PDG == %d or MCT_PDG == -%d'%(13,13))),
              len(df.query('MCT_PDG == %d or MCT_PDG == -%d'%(211,211))),
              len(df.query('MCT_PDG == %d or MCT_PDG == -%d'%(2212,2212))) ]
    total = float(np.sum(sizes))
    labels = ['Muon (%0.1f%%)' % (100.*(sizes[0]/total)), 
              'Pion (%0.1f%%)' % (100.*(sizes[1]/total)), 
              'Proton (%0.1f%%)' % (100.*(sizes[2]/total))]
    colors = ['yellowgreen', 'lightskyblue', 'lightcoral']
 
    explode = (0.1, 0, 0)  # explode muon slice
 
    # Plot
    plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            shadow=True, startangle=140)
     
    plt.axis('equal')
    if write_figures:
        print "Saving figure",'MCBNBSelectedRecoTrack_MID_piechart.png'
        plt.savefig(figdir+'MCBNBSelectedRecoTrack_MID_piechart.png')
    


# In[44]:

#Read in handscan info if working with data
if anatype == 'DataBNBSelectedRecoTrack':
    filedir = 'handscan_results/'
    myfile = 'handscan_results_kaleko.csv'
    
    hsdf = pd.read_csv(filedir + myfile,index_col=False)
       
    #Column names with spaces are hard to deal with, this is easy
    cols = hsdf.columns
    cols = cols.map(lambda x: x.replace(' ', '_') if isinstance(x, (str, unicode)) else x)
    hsdf.columns = cols
    
    #If no optional comments were typed in it shows up as NaN
    #Let's just make a column of "optional comments present" (boolean)
    
    hsdf['comments_present'] = hsdf['optional_comments'].notnull()
    
    #Let's rename things... 1_good_0_bad is only 0 if the track is definitely bad
    hsdf['definitely_bad'] = hsdf['1_good_0_bad'] == 0
    
    #maybe bad is either definitely_bad, or if comments are present (more conservative)
    hsdf['maybe_bad'] = hsdf['definitely_bad'] | hsdf['comments_present']
    
    print "total length of hsdf is",len(hsdf)
    print "number of definitely bad is",len(hsdf.query('definitely_bad'))
    print "number of maybe bad is",len(hsdf.query('maybe_bad'))
    
    #let's rename run, subrun, eventid columns to match the MCS df
    #for now just make a new column because it's easier
    hsdf['run'] = hsdf['Run']
    hsdf['subrun'] = hsdf['Subrun']
    hsdf['eventid'] = hsdf['Event_ID']
    
    df = df.merge(hsdf, on=['run','subrun','eventid'])

    segdf = segdf.merge(hsdf, on=['run','subrun','eventid'])


# In[45]:

#Breakdown of events that were handscanned as good, and those handscanned as bad
if anatype == 'DataBNBSelectedRecoTrack':

    extraquery = 'not maybe_bad'
    basic_comparison_fig(xvar='full_range_momentum',yvar='full_MCS_momentum',
                         plotname='MCS_range_momentum_DataRecoTracks_goodhandscan.png',
                            extraquery = extraquery,
                        addtext = 'Good Scan Sample')
    
    extraquery = 'maybe_bad'
    basic_comparison_fig(xvar='full_range_momentum',yvar='full_MCS_momentum',
                         plotname='MCS_range_momentum_DataRecoTracks_badhandscan.png',
                            extraquery = extraquery,
                        addtext = 'Bad Scan Sample')
"""

# In[46]:

#def highlandE(rms):
#    loverX = float(seglen)/14.
#    return np.sqrt((13.6*13.6*loverX*np.square((1+0.038*np.log(loverX))))/np.square(rms))


# In[47]:

def deflection_plot(binmin,binmax,nbins,extraquery=None,plotname=None, addtext=None, figname = None, 
                    customtitle = None, highland_constant = None, cut2percent = False, use_true_p = False):
    
    fig = plt.figure(figsize=(10,6))
    
    poop = plt.grid(True)
    
    myquery = 'run < 99999999'
    if extraquery is not None: myquery += ' and %s'%extraquery
    print myquery
    mybins = np.linspace(binmin,binmax,nbins)
    
    myvar1 = 'dthetayoverpredictedRMS_fromMCS'
    myvar2 = 'dthetaxoverpredictedRMS_fromMCS'
#     myvar1 = 'dthetayovertruepredictedRMS_momentumdepenentconstant'
#     myvar2 = 'dthetaxovertruepredictedRMS_momentumdepenentconstant'
    
    #testing out changing the 13.6 in the highland formula...
    #want to use true particle momentum to compute highland RMS
    if use_true_p:
        myvar1 = 'dthetayovertruepredictedRMS'
        myvar2 = 'dthetaxovertruepredictedRMS'

    #2.3265 standard deviations from 0 contain 98.000% of scatters
    thisquery1 = '%s > %f and %s < %f'%(myvar1,binmin,myvar1,binmax)
    thisquery2 = '%s > %f and %s < %f'%(myvar2,binmin,myvar2,binmax)

    myvals = np.append(
        segdf.query(myquery).\
        query(thisquery1)[myvar1].values,
        segdf.query(myquery).\
        query(thisquery2)[myvar2].values
    )
    
    #testing out changing the 13.6 in the highland formula...
    # it's /= because the highland formula is in the denominator
    if highland_constant is not None:
        myvals /= highland_constant/13.6
    
    mystd = np.std(myvals)
    mymean = np.mean(myvals)

    datahist = plt.hist(myvals,bins=mybins,normed=True,
                             alpha=0.5,label='$\Delta\\theta/RMS$ Values')#,
                            #weights = myweights)
   
    datahist_nonorm = np.histogram(myvals,bins=mybins,normed=False)
    
    binvals = datahist[0]
    bincenters = [ datahist[1][x] + (datahist[1][x+1]-datahist[1][x])/2 for x in xrange(len(datahist[1])-1) ]
    
    tempbinvals = binvals
    tempbincenters = bincenters
    if cut2percent is True:
        #work your way out from both sides removing bins until you're left with 98% of the central area
        totalentries = float(np.sum(binvals))
        subset = 0.02 * totalentries
        outercount = 0
        while outercount < subset:
            outercount += tempbinvals[0]
            outercount += tempbinvals[-1]
            tempbinvals = np.delete(tempbinvals,0)
            tempbinvals = np.delete(tempbinvals,-1)
            tempbincenters = np.delete(tempbincenters,0)
            tempbincenters = np.delete(tempbincenters,-1)
        
        binvals = tempbinvals
        bincenters = tempbincenters
        
    
    # Fit a normal distribution
    gmod = Model(gaussian)
    #initial random guesses of 1, 1, 2
    result = gmod.fit(binvals, x=bincenters, amp=1, cen=1, wid=2)

    #print(result.fit_report())
    #print help(result)
    #plt.plot(bincenters, binvals,         'bo')
    #plt.plot(bincenters, result.init_fit, 'k--')
    print result.params
    plt.plot(bincenters, result.best_fit, 'g-',             label='Gaussian Fit: $\sigma$ = %0.4f' % result.params['wid'],            linewidth=4)
  
    #data_mu, data_std = norm.fit(myvals)
    
    #myhighland = highlandE(data_std)

    # Plot the PDF.
    #x = np.linspace(binmin,binmax,100)
    #p = norm.pdf(x, data_mu, data_std)
    #plt.plot(x, p, 'g', linewidth=4,label='Gaussian Fit to Data')
    
    
    plt.title('%s'%titles[anatype],fontsize=14)
   
    if customtitle is not None:
        plt.title(customtitle,fontsize=14)
        
    plt.xlabel('Delta Theta / RMS',fontsize=16)
    blah = plt.legend(loc='best')
    leg = plt.legend(loc=2)
    
    if addtext is not None:
        plt.text(plt.xlim()[1]*0.25, plt.ylim()[1]*0.85, addtext, fontsize=20)
        
    if write_figures:
        fullfigname = figdir + 'Highland_validation_%s.png' % anatype
        if figname is not None:
            fullfigname = figdir + figname
            
        print "\n\n WRITING FIGURE %s! \n\n"%fullfigname
        plt.tight_layout()
        if plotname is not None: plt.savefig(figdir + plotname)
        else: plt.savefig(fullfigname)
            
    return result.params['wid']


# In[48]:

# if 'MCTrack' in anatype:
#     segdf.hist('true_segment_p',bins=np.linspace(0,0.2,100))
#     plt.yscale('log')
#     segdf.hist('true_segment_E',bins=np.linspace(0,1,100))
#     segdf.hist('resid_dist',bins=np.linspace(0,100,100))


# In[49]:

# # determining the "13.6" constant in the highland formula...
# if 'MCBNBMCTrack' in anatype:
# #     myconstants = np.linspace(12.0,12.3,10)
#     #myconstants = np.linspace(10,13,10)
#     myconstants = np.linspace(11.0,11.3,10)
# #    myconstants = np.linspace(11.5,11.75,10)
#     mywidths = []
#     myextraquery = None
#     #myextraquery = 'true_segment_E < 0.5'
#     myextraquery = 'true_segment_p > 1'
#     for myconstant in myconstants:
#         mywidths.append(deflection_plot(-5,5,50,extraquery=myextraquery,highland_constant=myconstant,cut2percent=False,use_true_p=True))

#     plt.figure(figsize=(10,6))
#     plt.plot(myconstants,mywidths,'ro--')
#     plt.title('Calibrating the 13.6 MCS Constant',fontsize=16)
#     plt.xlabel('Constant Value',fontsize=16)
#     plt.ylabel('Gaussian Fit Width',fontsize=16)
#     plt.grid(True)


# In[50]:

# myx = [.25,.35,.45,.55,.65,.75,.85,.95,1.05,1.15,1.25]
# lower_xerr = [ .05 ] * len(myx)
# upper_xerr = [ .05 ] * len(myx)
# upper_xerr[-1] = 0.5
# asymmetric_error = [lower_xerr,upper_xerr]

# # ax1.errorbar(x, y, xerr=asymmetric_error, fmt='o')
# myc = [12.7,11.8,11.5,11.4,11.3,11.25,11.15,11.15,11.05,11.05,11.03]
# plt.figure(figsize=(10,6))
# plt.errorbar(myx,myc,xerr=asymmetric_error,fmt='bo',label='Fit Values from MCTracks')
# plt.title("Highland Constant: 14cm Segments, Momentum Dependence",fontsize=16)
# plt.xlabel('True Segment Momentum [GeV]',fontsize=16)
# plt.ylabel('Fit Highland Constant',fontsize=16)
# plt.ylim((10,14))
# plt.xlim((0,2))
# plt.grid(True)

# from scipy.optimize import curve_fit
# def func(x,a,c):
#     return (a/(x*x)) + c
# popt, pcov = curve_fit(func,myx,myc)

# curvex = np.linspace(0.2,2,100)
# curvey = func(curvex,popt[0],popt[1])
# plt.plot(curvex,curvey,'r',label='Fit: [ a/(p^2) + c ]; a = %0.4f, c = %0.4f'%(popt[0],popt[1]))
# plt.legend()


# In[ ]:

# In[51]:
"""
# For MCBNBSelectedRecoTrack only make highland validation plot for the true muons
myextraquery = None
if anatype == 'MCBNBSelectedRecoTrack':
    myextraquery = 'MCT_PDG == 13 or MCT_PDG == -13'
# For real data events, only make plot for those events handscanned as "good"
if anatype == 'DataBNBSelectedRecoTrack':
    myextraquery = 'not maybe_bad'
    
dummy = deflection_plot(-5,5,50,extraquery=myextraquery,use_true_p=False)
#dummy = deflection_plot(-5,5,50,extraquery=myextraquery,highland_constant=13.6,use_true_p=False)

if anatype == 'DataBNBSelectedRecoTrack':
    myextraquery = 'not maybe_bad'
    dummy = deflection_plot(-5,5,50,extraquery=myextraquery,addtext='Good Scan Sample',                            figname = 'Highland_validation_DataBNBSelectedRecoTrack_goodscan.png')
    myextraquery = 'maybe_bad'
    dummy = deflection_plot(-5,5,50,extraquery=myextraquery,addtext='Bad Scan Sample',                            figname = 'Highland_validation_DataBNBSelectedRecoTrack_badscan.png',
                           customtitle = 'Selected, Poorly Reconstructed Tracks from NumuCC Data')

"""
# In[52]:

#reco-true/true
def fractional_bias_resolution_plot(xvar = 'trkMom_RangeMu', xbins = np.linspace(0.35,1,10),
                   yvar = 'trkMom_Mu', slicevar = 'trkMom_RangeMu',
                         plot_bin_distributions = False, extraquery = None, slicetitlebase = None,
                        slicebins = np.linspace(-0.8,.8,20),
                        biasplotname = None, resplotname = None, biasmainfig_ylims = None,
                         resmainfig_ylims = None,
                        usegausfit = False,
                                   myaddtext = None, myresylabel = None,mybiasylabel = None,
                                   myrestitle = None, mybiastitle = None):

    binning = xbins
    binwidth = float(binning[1]-binning[0])
    bincenters = binning + (binwidth/2)
    myreses, mystds, myerrs_bias, myerrs_res = [], [], [], []

    for x in xrange(len(binning)-1):
        
        binmin = binning[x]
        binmax = binning[x+1]
        #print "binmin = %0.2f, binmax = %0.2f"% ( binmin, binmax )
        myquery = '%s > %f and %s < %f'%(slicevar,binmin,slicevar,binmax)
        if anatype == 'DataBNBSelectedRecoTrack': 
            myquery += ' and not maybe_bad'
        if extraquery is not None:
            myquery += ' and %s' % extraquery

        print "myquery:",myquery

        mydf = df.query(myquery)
        true = mydf[xvar].values
        reco = mydf[yvar].values
        mymean = ((reco-true)/true).mean()
        mystd = ((reco-true)/true).std()
        
        if plot_bin_distributions:
            plt.figure(figsize=(5,3))
            datahist = plt.hist((reco-true)/true,color='r',alpha=0.5,label='%d Entries'%len((reco-true)/true),                     bins=slicebins)
            temp = 'GeV'
            if 'length' in slicevar: temp = 'cm'
            titlestring = '$\\frac{%s - %s}{%s}$ for $%s$ in %0.2f $\\rightarrow$ %0.2f %s'%             (latextitles[yvar],latextitles[xvar],latextitles[xvar],latextitles[slicevar],binmin,binmax,temp)
            t =plt.title(titlestring,fontsize=16)
            #move the title up a bit
            t.set_y(1.04) 
            plt.grid(False)
            
            # Plot gaussian on each bin distribution
            if usegausfit:
                 
                slicebinvals = datahist[0]
                slicebincenters = [ datahist[1][x] + (datahist[1][x+1]-datahist[1][x])/2 for x in xrange(len(datahist[1])-1) ]
    
                # Fit a normal distribution
                gmod = Model(gaussian)
                #initial random guesses of 1, 1, 2
                result = gmod.fit(slicebinvals, x=slicebincenters, amp=1, cen=1, wid=2)

                thisx = np.linspace(np.min(slicebins),np.max(slicebins),100)
                thisy = gaussian(thisx, result.params['amp'], result.params['cen'], result.params['wid'])
                plt.plot(thisx, thisy, 'g-',                     label='Gaus Fit',                     linewidth=2)
                addtext = 'Fit: \n$\sigma$ = %0.2f, \n$\mu$ = %0.2f'%                (np.abs(result.params['wid']),result.params['cen'])
                plt.text(plt.xlim()[1]*0.35, plt.ylim()[1]*0.1, addtext, fontsize=14)
                
                # If use gaus fit, use the result of that instead of straight mean and RMS
                mymean = result.params['cen']
                #Somehow when there are like 2 data points you get a negative width?!
                mystd  = np.abs(result.params['wid'])
                # If the fit doesn't converge, 'wid' will be the initial guess of 2.0... don't use these points then
                # (this only happens when no entries in the sliced histogram)
                if int(result.params['wid']) == 2:
                    print "WARNING: FIT DIDN'T CONVERGE!"
                    mymean = ((reco-true)/true).mean()
                    mystd = ((reco-true)/true).std()
                    
             
                
                
            plt.xlabel('$\\frac{%s - %s}{%s}$'%(latextitles[yvar],latextitles[xvar],latextitles[xvar]),fontsize=16)
            plt.ylabel('Counts',fontsize=16)
            plt.xlim((np.min(slicebins),np.max(slicebins)))
            plt.legend(loc=1)
            if write_figures and slicetitlebase is not None:
                fullfigname = figdir + slicetitlebase + '_slice_%0.2f_%0.2f.png'%(binmin,binmax)
                print '\n\n WRITING A FIGURE!! %s\n\n'%fullfigname
                plt.tight_layout()
                plt.savefig(fullfigname)
        
        
        myerr_bias = mystd / np.sqrt( float(len(true)) )
        myerr_res = mystd / np.sqrt( float(2*len(true)) )
        myreses.append( mymean )
        mystds.append( mystd )
        myerrs_bias.append( myerr_bias )
        myerrs_res.append( myerr_res )
        
        
    #BIAS PLOT 
    plt.figure(figsize=(10,6))
    plt.errorbar(bincenters[:-1],myreses,yerr=myerrs_bias,xerr=binwidth/2,fmt='ro',label='Mean of Gaussian Fit, Errors = std/sqrt(N)')
    plt.ylabel('GausMean($\\frac{%s - %s}{%s}$)'%(latextitles[yvar],latextitles[xvar],latextitles[xvar]),fontsize=25)
    if mybiasylabel is not None:
        plt.ylabel(mybiasylabel,fontsize=15)
    plt.xlabel('%s'%titles[slicevar],fontsize=15)
    plt.grid(True)
    plt.legend(loc='best')
    t = plt.title('Fractional Bias: %s'%titles[anatype],fontsize=16)
    if mybiastitle is not None:
        plt.title(mybiastitle,fontsize=16)
    if mybiastitle is None: t.set_y(1.04)
    if biasmainfig_ylims is not None:
        blah = plt.ylim(biasmainfig_ylims)
        
    
    if myaddtext is not None:
        plt.text(plt.xlim()[1]*0.45, plt.ylim()[0] + (plt.ylim()[1]-plt.ylim()[0])*0.05, myaddtext, fontsize=30)
        
    if write_figures and biasplotname is not None:
        print " \n\n Writing the main bias figure!! %s\n\n" % (figdir+biasplotname)
        plt.tight_layout()
        plt.savefig(figdir + biasplotname)
        
    #RESOLUTION PLOT
    plt.figure(figsize=(10,6))
    plt.errorbar(bincenters[:-1],mystds,yerr=myerrs_res,xerr=binwidth/2,fmt='bo',label='Width of Gaussian Fit, Errors = std/sqrt(2N)')
    
    plt.ylabel('GausWidth($\\frac{%s - %s}{%s}$)'%(latextitles[yvar],latextitles[xvar],latextitles[xvar]),fontsize=25)
    if myresylabel is not None:
        plt.ylabel(myresylabel,fontsize=15)
    plt.xlabel('%s'%titles[slicevar],fontsize=15)
    plt.grid(True)
    t = plt.title('Momentum Resolution: %s' % titles[anatype],fontsize=16)
    if myrestitle is not None:
        plt.title(myrestitle,fontsize=16)
    plt.legend(loc='best')
    #move the title up a bit
    if myrestitle is None: t.set_y(1.05) 
                    
    if resmainfig_ylims is not None:
        blah = plt.ylim(resmainfig_ylims)
        
    if myaddtext is not None:
        plt.text(plt.xlim()[1]*0.45, plt.ylim()[0] + (plt.ylim()[1]-plt.ylim()[0])*0.05, myaddtext, fontsize=30)
        
    if write_figures and resplotname is not None:
        print " \n\n Writing the main resolution figure!! %s \n\n" % (figdir+resplotname)
        plt.tight_layout()
        plt.savefig(figdir + resplotname)


# In[53]:

print len(df)
#print "len(df.query('trkMom_Mu < 7.4'))"
#print len(df.query('trkMom_Mu < 7.4'))


# In[54]:

# Fractional bias plot for True vs Range energy for single MCTrack section only
if anatype == 'MCBNBMCTrack':
    fractional_bias_resolution_plot(xvar='true_E',yvar='full_range_energy',xbins=np.linspace(0.35,2,10),
                       plot_bin_distributions = True,
                       slicevar = 'true_E',
                       slicetitlebase = 'true_range_resolution_%s'%anatype,
                       slicebins = np.linspace(-0.1,0.1,20),
                       biasplotname = 'true_range_bias_%s.png'%anatype,
                       resplotname = 'true_range_resolution_%s.png'%anatype,
                       biasmainfig_ylims = (-.05,.05),
                       resmainfig_ylims = (0,.05),
                       usegausfit = True,
                       myaddtext = 'MicroBooNE Simulation',
                                   myresylabel = 'Fractional Resolution',
                                   mybiasylabel='Fractional Bias',
                                   mybiastitle='Range-Based Fractional Bias',
                                   myrestitle='Range-Based Momentum Resolution')


# In[ ]:




# In[55]:
# Fractional bias plot for MCS vs true momentum for EXITING MCTRACK SECTION
if 'Exiting' in anatype:
    fractional_bias_resolution_plot(xvar='true_momentum_inverse',yvar='full_MCS_momentum_inverse',xbins=np.linspace(0.35,4,10),
                       plot_bin_distributions = True,
                       slicevar = 'true_momentum',
                       slicetitlebase = 'MCS_true_exiting_resolution_%s'%anatype,
                       slicebins = np.linspace(-2,2,50),
                       biasplotname = 'MCS_true_exiting_bias_%s.png'%anatype,
                       resplotname = 'MCS_true_exiting_resolution_%s.png'%anatype,
                       biasmainfig_ylims = (-.1,.1),
                       resmainfig_ylims = (0,.40),
                       usegausfit = True,
                       extraquery = 'full_MCS_momentum < 7.4')
    
    mylengths = df['full_length']
    plt.figure(figsize=(10,6))
    blah = plt.hist(mylengths,bins=np.linspace(0,600,100),alpha=0.5)
    blah = plt.grid(True)
    blah = plt.title('Track Length Distribution: %s'%anatype,fontsize=16)
    plt.xlabel('Length of Track Contained [cm]',fontsize=16)
    plt.ylabel('Counts',fontsize=16)
   
    plotname = 'track_length_distribution_%s.png'%anatype
    if write_figures:
        print 'WRITING A FIGURE! %s' % (figdir + plotname)
        plt.tight_layout()
        plt.xlim((50,600))
        plt.savefig(figdir + plotname)
    


# In[ ]:
"""
# Fractional bias plot for MCS vs true momentum for everything except real data
if 'Data' not in anatype:
    fractional_bias_resolution_plot(xvar='true_momentum_inverse',yvar='full_MCS_momentum_inverse',xbins=np.linspace(0.35,2,10),
                       plot_bin_distributions = True,
                       slicevar = 'true_momentum',
                       slicetitlebase = 'MCS_true_resolution_%s'%anatype,
                       slicebins = np.linspace(-0.8,0.8,50),
                       biasplotname = 'MCS_true_bias_%s.png'%anatype,
                       resplotname = 'MCS_true_resolution_%s.png'%anatype,
                       biasmainfig_ylims = (-.1,.1),
                       resmainfig_ylims = (0,.2),
                       usegausfit = True,
                       extraquery = 'full_MCS_momentum < 7.4')

"""
# In[ ]:

# Fractional bias plot for MCS vs true momentum for EXITING MCTRACK SECTION
# if anatype == 'MCBNBRecoTrackExiting':
#     fractional_bias_resolution_plot(xvar='true_momentum_inverse',yvar='full_MCS_momentum_inverse',xbins=np.linspace(0.35,4,10),
#                        plot_bin_distributions = True,
#                        slicevar = 'true_momentum',
#                        slicetitlebase = 'MCS_true_exiting_resolution_%s'%anatype,
#                        slicebins = np.linspace(-0.8,0.8,20),
#                        biasplotname = 'MCS_true_exiting_bias_%s.png'%anatype,
#                        resplotname = 'MCS_true_exiting_resolution_%s.png'%anatype,
#                        biasmainfig_ylims = (-.5,.5),
#                        resmainfig_ylims = (0,.15),
#                        usegausfit = True)


# In[ ]:

# Fractional bias plot for MCS vs range momentum for EXITING MCTRACK SECTION
# NOTE THAT RANGE MOMENTUM REALLY MEANS "LENGTH OF TRACK CONTAINED" HERE!!
if anatype == 'MCBNBMCTrackExiting':
    fractional_bias_resolution_plot(xvar='true_momentum_inverse',yvar='full_MCS_momentum_inverse',xbins=np.linspace(0,1000,10),
                       plot_bin_distributions = True,
                       slicevar = 'full_length',
                       slicetitlebase = 'MCS_range_exiting_resolution_%s'%anatype,
                       slicebins = np.linspace(-0.8,0.8,20),
                       biasplotname = 'MCS_range_exiting_bias_%s.png'%anatype,
                       resplotname = 'MCS_range_exiting_resolution_%s.png'%anatype,
                       biasmainfig_ylims = (-0.2,0.2),
                       resmainfig_ylims = (0,0.2),
                       usegausfit = True)

    


# In[ ]:

# Fractional bias plot for MCS vs true momentum for EXITING MCTRACK SECTION
# This is broken up for different bins of length!
if anatype == 'MCBNBMCTrackExiting':
    
    extraqueries = ['full_length < 150.',
                   'full_length > 150. and full_length < 250.',
                   'full_length > 250. and full_length < 350.',
                   'full_length > 350. and full_length < 450.']
    
    counter = 0
    temptitles = [ '100_150', '150_250', '250_350', '350_450' ]
    
    for extraquery in extraqueries:
        fractional_bias_resolution_plot(xvar='true_momentum_inverse',yvar='full_MCS_momentum_inverse',xbins=np.linspace(0.35,4,10),
                           plot_bin_distributions = True,
                           slicevar = 'true_momentum',
                           slicetitlebase = 'MCS_true_exiting_resolution_%s'%anatype,
                           slicebins = np.linspace(-0.8,0.8,20),
                           biasplotname = 'MCS_true_exiting_bias_L%s_%s.png'%(temptitles[counter],anatype),
                           resplotname = 'MCS_true_exiting_resolution_L%s_%s.png'%(temptitles[counter],anatype),
                           biasmainfig_ylims = (-.5,.5),
                           resmainfig_ylims = (0,.3),
                           extraquery=extraquery,
                           usegausfit = True)
        counter += 1


# In[ ]:

# Nominal fractional bias plot for all samples
# For MCBNBSelectedRecoTrack only make plot for the true muons
myextraquery = 'trkIsContained and simID==13 and trkLength>100 and trkMom_Mu<7'
myextraquery2 = 'trkIsContained and (trkLength>=(0.9*simLength) and trkLength<(1.1*simLength)) and simID==13 and trkLength>100 and trkMom_Mu<7'
myextraqueryExiting = 'trkIsContained==False and simID==13 and trkLength>100 and trkMom_Mu<7'
myextraqueryExiting2 = 'trkIsContained==False and simID==13 and trkLength>100 and trkMom_Mu<7 and (trkLength>=(0.9*simLength) and trkLength<(1.1*simLength))'
myextraqueryExitingTest = 'trkIsContained==False and simID==13 and trkLength>100 and trkMom_Mu<7 and trkSegRadLengthsAvg < 1.02'
#myextraqueryProton = 'trkIsContained and (trkLength>=(0.9*simLength) and trkLength<(1.1*simLength)) and trkid2==2212 and trkLength>100'
if anatype == 'MCBNBSelectedRecoTrack':
    myextraquery = 'MCT_PDG == 13 or MCT_PDG == -13'
# For real data events, only make plot for those events handscanned as "good"
if anatype == 'DataBNBSelectedRecoTrack':
    myextraquery = 'not maybe_bad'
"""
fractional_bias_resolution_plot(xvar='mcs_mom_noSC_inverse',yvar='mcs_mom_withSC_inverse',
                     xbins=np.linspace(0.35,2,10),
                       plot_bin_distributions = True,
                       slicevar = 'mcs_mom_noSC',
                       extraquery = None,
                       slicetitlebase = 'SCE_resolution_bin_%s'%anatype,
                       slicebins = np.linspace(-0.6,0.6,20),
                       biasplotname = 'SCE_bias_%s.png'%anatype,
                       resplotname = 'SCE_resolution_%s.png'%anatype,
                       biasmainfig_ylims = (-.10,.10),
                       resmainfig_ylims = (0,.15),
                        usegausfit = True)

fractional_bias_resolution_plot(xvar='mcs_mom_noSC',yvar='mcs_mom_withSC',
                     xbins=np.linspace(0.35,2,10),
                       plot_bin_distributions = True,
                       slicevar = 'mcs_mom_noSC',
                       extraquery = None,
                       slicetitlebase = 'SCE_noInverse_resolution_bin_%s'%anatype,
                       slicebins = np.linspace(-0.6,0.6,20),
                       biasplotname = 'SCE_noInverse_bias_%s.png'%anatype,
                       resplotname = 'SCE_noInverse_resolution_%s.png'%anatype,
                       biasmainfig_ylims = (-.10,.10),
                       resmainfig_ylims = (0,.15),
                        usegausfit = True)


  
# Trying inverse momentum (truth)
# Muon
fractional_bias_resolution_plot(xvar='simMom_inverse',yvar='trkMom_Mu_inverse',
                     xbins=np.linspace(0.35,2,10),
                       plot_bin_distributions = True,
                       slicevar = 'simMom',
                       extraquery = myextraquery2,
                       slicetitlebase = 'giuseppeMuplusCosmic_MCS_truth_resolution_%s'%anatype,
                       slicebins = np.linspace(-0.6,0.6,20),
                       biasplotname = 'giuseppeMuplusCosmic_MCS_truth_bias_%s.png'%anatype,
                       resplotname = 'giuseppeMuplusCosmic_MCS_truth_resolution_%s.png'%anatype,
                       biasmainfig_ylims = (-.10,.10),
                       resmainfig_ylims = (0,.15),
                        usegausfit = True)
"""
# Trying inverse momentum (truth)
fractional_bias_resolution_plot(xvar='simMom_inverse',yvar='trkMom_Mu_inverse',
                       xbins=np.linspace(0.35,4,10),
                       plot_bin_distributions = True,
                       slicevar = 'simMom',
                       extraquery = myextraqueryExiting2,
                       slicetitlebase = 'Aug04_test2_resolution_%s'%anatype,
                       slicebins = np.linspace(-0.6,0.6,20),
                       biasplotname = 'Aug04_test2_bias_%s.png'%anatype,
                       resplotname = 'Aug04_test2_resolution_%s.png'%anatype,
                       biasmainfig_ylims = (-.10,.30),
                       resmainfig_ylims = (0,.20),
                        usegausfit = True)

"""
fractional_bias_resolution_plot(xvar='simMom_inverse',yvar='trkMom_Mu_inverse',
                       xbins=np.linspace(0.35,4,10),
                       plot_bin_distributions = True,
                       slicevar = 'simMom',
                       extraquery = myextraqueryContained,
                       slicetitlebase = '1mrad_under7_mcc8_1_nuOnly_giuseppeNu_contained_MCS_truth_resolution_%s'%anatype,
                       slicebins = np.linspace(-0.6,0.6,20),
                       biasplotname = '1mrad_under7_mcc8_1_nuOnly_giuseppeNu_contained_MCS_truth_bias_%s.png'%anatype,
                       resplotname = '1mrad_under7_mcc8_1_nuOnly_giuseppeNu_contained_MCS_truth_resolution_%s.png'%anatype,
                       biasmainfig_ylims = (-.50,.50),
                       resmainfig_ylims = (0,.4),
                        usegausfit = True)






# Proton
fractional_bias_resolution_plot(xvar='trkMom_RangeP_inverse',yvar='trkMom_P_inverse',
                     xbins=np.linspace(0.35,2,10),
                       plot_bin_distributions = True,
                       slicevar = 'trkMom_RangeP',
                       extraquery = myextraquery2,
                       slicetitlebase = 'giuseppeCosmic_P_MCS_range_resolution_%s'%anatype,
                       slicebins = np.linspace(-0.6,0.6,20),
                       biasplotname = 'giuseppeCosmic_P_MCS_range_bias_%s.png'%anatype,
                       resplotname = 'giuseppeCosmic_P_MCS_range_resolution_%s.png'%anatype,
                       biasmainfig_ylims = (-.10,.10),
                       resmainfig_ylims = (0,.15),
                        usegausfit = True)


# Trying inverse momentum (range)
# Muon
fractional_bias_resolution_plot(xvar='trkMom_RangeMu_inverse',yvar='trkMom_Mu_inverse',
                     xbins=np.linspace(0.35,2,10),
                       plot_bin_distributions = True,
                       slicevar = 'trkMom_RangeMu',
                       extraquery = myextraquery2,
                       slicetitlebase = 'giuseppeMuplusCosmic_MCS_range_resolution_%s'%anatype,
                       slicebins = np.linspace(-0.6,0.6,20),
                       biasplotname = 'giuseppeMuplusCosmic_MCS_range_bias_%s.png'%anatype,
                       resplotname = 'giuseppeMuplusCosmic_MCS_range_resolution_%s.png'%anatype,
                       biasmainfig_ylims = (-.10,.10),
                       resmainfig_ylims = (0,.15),
                        usegausfit = True)

# Proton
fractional_bias_resolution_plot(xvar='trkMom_RangeP_inverse',yvar='trkMom_P_inverse',
                     xbins=np.linspace(0.35,2,10),
                       plot_bin_distributions = True,
                       slicevar = 'trkMom_RangeP',
                       extraquery = myextraquery2,
                       slicetitlebase = 'giuseppeCosmic_P_MCS_range_resolution_%s'%anatype,
                       slicebins = np.linspace(-0.6,0.6,20),
                       biasplotname = 'giuseppeCosmic_P_MCS_range_bias_%s.png'%anatype,
                       resplotname = 'giuseppeCosmic_P_MCS_range_resolution_%s.png'%anatype,
                       biasmainfig_ylims = (-.10,.10),
                       resmainfig_ylims = (0,.15),
                        usegausfit = True)

# In[ ]:

##Here I fit the resolution to a realistic resolution equation
#def res_eqtn(E,a,b,c):
#    return np.sqrt(np.square(a/np.sqrt(E)) + np.square(b/E) + np.square(c))
#
##def res_eqtn(x,a):
##    return a/np.sqrt(x)
#
##def res_eqtn(E,a,b,c):
##    return np.sqrt(np.square(a/np.sqrt(1./E)) + np.square(b/(1./E)) + np.square(c))
#
#
##reco-true/true
#def resolution_plot(xvar = 'full_range_momentum', xbins = np.linspace(0.35,1,10),
#                   yvar = 'full_MCS_momentum', slicevar = 'full_range_momentum',
#                         extraquery = None,
#                        plotname = None, mainfig_ylims = None):
#    
#    #This includes a fit to a realistic resolution equation
#    
#    binning = xbins
#    binwidth = float(binning[1]-binning[0])
#    bincenters = binning + (binwidth/2)
#    myreses, mystds, myerrs = [], [], []
#    myN = []
#    for x in xrange(len(binning)-1):
#        
#        binmin = binning[x]
#        binmax = binning[x+1]
#        #print "binmin = %0.2f, binmax = %0.2f"% ( binmin, binmax )
#        myquery = '%s > %f and %s < %f'%(slicevar,binmin,slicevar,binmax)
#        if anatype == 'DataBNBSelectedRecoTrack':
#            myquery += ' and not maybe_bad'
#        if extraquery is not None:
#            myquery += ' and %s' % extraquery
#            
#        mydf = df.query(myquery)
#        true = mydf[xvar].values
#        reco = mydf[yvar].values
#        mymean = ((reco-true)/true).mean()
#        mystd = ((reco-true)/true).std()
#        #Error of standard deviation is sigma/sqrt(2N)
#        myerr = mystd / np.sqrt( 2 * float(len(true)) )
#        myN.append(len(true))
#        myreses.append( mymean )
#        mystds.append( mystd )
#        myerrs.append( myerr )
#            
#    plt.figure(figsize=(10,6))
#    plt.errorbar(bincenters[:-1],mystds,yerr=myerrs,fmt='o--',label='Std of Distribution, Errors = std/sqrt(2N)')
#    plt.ylabel('std($\\frac{%s - %s}{%s}$)'%(latextitles[yvar],latextitles[xvar],latextitles[xvar]),fontsize=25)
#    plt.xlabel('%s'%titles[xvar],fontsize=15)
#    plt.grid(True)
#    t = plt.title('Momentum Resolution: %s' % titles[anatype],fontsize=16)
#    #move the title up a bit
#    t.set_y(1.04) 
#                    
#    if mainfig_ylims is not None:
#        blah = plt.ylim(mainfig_ylims)
#    if write_figures and plotname is not None:
#        print " \n\n Writing the main figure!! %s \n\n" % (figdir+plotname)
#        plt.tight_layout()
#        plt.savefig(figdir + plotname)
#    
##     if draw_fit:
##         popt, pcov = curve_fit(res_eqtn, bincenters[:-1], mystds, sigma=myerrs)
#        
##         fitx = np.linspace(0.35,1,100)
##         fity = res_eqtn(fitx,popt[0],popt[1],popt[2])
##         labelstring = 'Resolution Fit: '
##         labelstring += '$[\\frac{a}{\sqrt{E}} \circ \\frac{b}{E} \circ c]$'
##         plt.plot(fitx,fity,'r-',label=labelstring)
#   
##         plt.legend(loc='best')
#    
##         return popt, pcov
##     else:
##         plt.legend(loc='best')
##         return (-1, -1)


# In[ ]:

## Fractional resolution plot for True vs Range energy for single MCTrack section only
#if anatype == 'SingleMuonMCTrack':
#    resolution_plot(xvar='true_E',yvar='full_range_energy',xbins=np.linspace(0.35,2,10),
#                       slicevar = 'full_range_energy',
#                       plotname = 'true_range_resolution_%s.png'%anatype,
#                       mainfig_ylims = (0,.05))
#   
## Nominal fractional bias plot for all samples
## For MCBNBSelectedRecoTrack only plot for the true muons
#myextraquery = None
#if anatype == 'MCBNBSelectedRecoTrack':
#    myextraquery = 'MCT_PDG == 13 or MCT_PDG == -13'
## For real data events, only make plot for those events handscanned as "good"
#if anatype == 'DataBNBSelectedRecoTrack':
#    myextraquery = 'not maybe_bad'
#    w
#resolution_plot(xvar='full_range_momentum',yvar='full_MCS_momentum',xbins=np.linspace(0.35,2,10),
#                       slicevar = 'full_range_momentum',
#                       extraquery = myextraquery,
#                       plotname = 'MCS_range_resolution_%s.png'%anatype,
#                       mainfig_ylims = (0,.25))
#


# In[ ]:




# In[ ]:




# In[ ]:

#segdf.columns.values


# In[ ]:

#segdf.hist('llbf',bins=np.linspace(0,1600,100))


# In[ ]:


myquery = 'llbf < 1500'
myx = segdf.query(myquery)['llbf'].values
myy = segdf.query(myquery)['full_MCS_E'].values
blah = plt.hist2d(myx,myy,bins=100,cmin=1)
blah = plt.colorbar()


# In[ ]:




# In[30]:

def beta(true_p_GEV):
    return np.sqrt( 1 - ((0.106*0.106)/(true_p_GEV * true_p_GEV + 0.106*0.106)) )

def func(x,a,c):
    return (a/(x*x)) + c

def highland_modified(p):
    a = 0.1049
    c = 11.0038
    return func(p,a,c)/p*beta(p)

def highland_normal(p):
    return 13.6 / p*beta(p)


myps = np.linspace(0.1,4,100)
myhl = highland_modified(myps)
myhl_normal = highland_normal(myps)

plt.figure(figsize=(10,6))
plt.plot(myps,myhl,label='14cm Segment Highland With Momentum-Dependent Constant')
plt.plot(myps,myhl_normal,label='14cm Segment Highland With 13.6 Constant')
plt.legend()
plt.grid(True)
plt.ylim((0,20))
plt.xlabel('Segment Momentum [GeV]',fontsize=16)
plt.ylabel('RMS Angular Scatter [mrad]',fontsize=16)
plt.title('Highland Formula Visualized in Two Forms',fontsize=16)


# In[ ]:

segdf.hist('n_traj_points',bins=np.linspace(0,100,102))


# In[ ]:

np.mean(segdf['n_traj_points'])


# In[36]:

highland_normal(4.55)

"""



