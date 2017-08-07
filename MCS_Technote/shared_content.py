import numpy as np
import matplotlib.pylab as plt
import pandas as pd
from scipy.stats import norm
from root_numpy import root2array
from scipy.optimize import curve_fit
from lmfit import  Model #better gaussian model


titles = { 
           'SingleMuonMCTrack'        : 'Fully Contained Single Muon MCTracks',
#           'SingleMuonRecoTrack'      : 'Fully Contained, Well Reconstructed Single Muon Tracks',
           'SingleMuonRecoTrack'      : 'Exiting Reconstructed Single Muon Tracks',
           'DataBNBSelectedRecoTrack' : 'Selected, Well Reconstructed Tracks from NumuCC Data',
           'MCBNBSelectedRecoTrack'   : 'Selected, Well Reconstructed Tracks from NumuCC Simulation',
#           'MCBNBRecoTrack'           : 'MC numuCC BNB Truth-Selected, Well Reconstructed Tracks',
#           'MCBNBRecoTrack'           : 'MC numuCC BNB Exiting Reco Tracks',
           'MCBNBRecoTrack'           : 'MC numuCC BNB Tracks',
           'MCBNBMCTrack'             : 'MC numuCC BNB Truth-Selected MCTracks',
           'MCBNBMCTrackExiting'      : 'MC numuCC BNB Truth-Selected EXITING MCTracks',
           'full_MCS_energy'          : 'MCS Total Energy [GeV]',
           'full_range_energy'        : 'Range-Based Total Energy [GeV]',
           'full_integrated_range_energy'        : 'Integrated Range-Based Total Energy [GeV]',
           'trkMom_Mu'        : 'MCS Momentum [GeV]',
           'trkMom_MuFwd'        : 'MCS Momentum (Forward) [GeV]',
#           'mcs_mom_noSC'        : 'MCS Momentum, no SC [GeV]',
#           'mcs_mom_withSC'        : 'MCS Momentum, with SC [GeV]',
           'trkMom_P'        : 'MCS Proton Momentum [GeV]',
           'trkMom_RangeMu'      : 'Range-Based Momentum [GeV]',
           'trkMom_RangeP'      : 'Range-Based Proton Momentum [GeV]',
           'full_integrated_range_momentum'      : 'Integrated Range-Based Momentum [GeV]',
           'trkMom_Mu_inverse'        : 'Inverse MCS Momentum [GeV^-1]',
           'trkMom_MuFwd_inverse'        : 'Inverse MCS Momentum (Forward) [GeV^-1]',
           'trkMom_P_inverse'        : 'Inverse MCS Proton Momentum [GeV^-1]',
           'trkMom_RangeMu_inverse'      : 'Inverse Range-Based Momentum [GeV^-1]',
           'trkMom_RangeP_inverse'      : 'Inverse Range-Based Proton Momentum [GeV^-1]',
           'true_E'                   : 'True Total Energy [GeV]',
           'simMom'                   : 'True Momentum [GeV]',
           'simMom_inverse'                   : 'Inverse True Momentum [GeV^-1]',
           'full_length'              : 'Length Of Track Contained [cm]', 
           'trkLength' : 'Length of Track [cm]',
           'trkNHits' : 'Number of Hits in Track',
           'trkSegRadLengthsMin' : 'Minimum Segment Rad Length',
           'trkSegRadLengthsAvg' : 'Average Segment Rad Length',
           'fracDiff' : 'Fractional Difference'
         }

latextitles = {
           'full_range_energy'   : 'E_{Range}',
           'trkMom_RangeMu' : 'p_{Range}',
           'trkMom_RangeP' : 'Proton p_{Range}',
           'full_integrated_range_energy'   : 'E_{Integrated Range}',
           'full_integrated_range_momentum' : 'p_{Integrated Range}',
           'full_MCS_energy'     : 'E_{MCS}',
           'trkMom_Mu'   : 'p_{MCS}',
           'trkMom_MuFwd'   : 'p_{MCS} (Forward)',
#           'mcs_mom_noSC'   : 'p_{MCS, no SC}',
#           'mcs_mom_noSC_inverse'   : 'p_{MCS, no SC}^{-1}',
#           'mcs_mom_withSC'   : 'p_{MCS, with SC}',
#           'mcs_mom_withSC_inverse'   : 'p_{MCS, with SC}^{-1}',
           'trkMom_P'   : 'Proton p_{MCS}',
           'trkMom_Mu_inverse'        : 'p_{MCS}^{-1}',
           'trkMom_MuFwd_inverse'        : 'p_{MCS}^{-1} (Forward)',
           'trkMom_P_inverse'        : 'Proton p_{MCS}^{-1}',
           'trkMom_RangeMu_inverse'      : 'p_{Range}^{-1}',
           'trkMom_RangeP_inverse'      : 'Proton p_{Range}^{-1}',
           'true_E'              : 'E_{True}',
           'simMom'              : 'p_{True}',
           'simMom_inverse'              : 'p_{True}^{-1}',
           'full_length'         : 'L_{Track}',
           'trkLength' : 'Len. of Track [cm]',
           'trkNHits' : 'Num. Hits in Track',
           'trkSegRadLengthsMin' : 'Min. Segment Rad Length',
           'trkSegRadLengthsAvg' : 'Avg. Segment Rad Length',
           'fracDiff' : 'Fractional Diff.'
         }

def gaussian(x, amp, cen, wid):
    "1-d gaussian: gaussian(x, amp, cen, wid)"
    return (amp/(np.sqrt(2*np.pi)*wid)) * np.exp(-(x-cen)**2 /(2*wid**2))

def get_dfs(myfile):
    #This df has track-by-track information (MCS energy, range energy, etc)
#    df = pd.DataFrame( root2array ( myfile, 'MCS_bias_tree' ) )
    df = pd.DataFrame( root2array ( myfile, '/trajmcsntupleNu/tree' ) )
#    df = pd.DataFrame( root2array ( myfile, 'MCS_SCE_comparison' ) )

#    df['full_MCS_momentum_inverse'] = 1./df['full_MCS_momentum']
#    df['full_range_momentum_inverse'] = 1./df['full_range_momentum']
    df['trkMom_Mu_inverse'] = 1./df['trkMom_Mu']
    df['trkMom_MuFwd_inverse'] = 1./df['trkMom_MuFwd']
#    df['mcs_mom_noSC_inverse'] = 1./df['mcs_mom_noSC']
#    df['mcs_mom_withSC_inverse'] = 1./df['mcs_mom_withSC']
    df['trkMom_RangeMu_inverse'] = 1./df['trkMom_RangeMu']
#    df['trkMom_P_inverse'] = 1./df['trkMom_P']
#    df['trkMom_RangeP_inverse'] = 1./df['trkMom_RangeP']

#    df['trkSegRadLengthsMin'] = df['trkSegRadLengths'].min(level=0)

    # List for min seg array
    min_seg_arr = []
    avg_seg_arr = []

    for min_seg in range(len(df['trkSegRadLengths'])):

        single_track_segments_row     = df['trkSegRadLengths'][min_seg]

        single_track_segments_row_arr = np.asarray(single_track_segments_row)

        #print single_track_segments_row_arr
        if len(single_track_segments_row_arr)>0: 
            min_seg_arr.append( single_track_segments_row_arr.min() )
            avg_seg_arr.append( np.average(single_track_segments_row_arr) )
            
        else: 
            min_seg_arr.append(-999.)
            avg_seg_arr.append(-999.)

    print "The length of the 'avg_seg_arr' variable = %d." %( len(avg_seg_arr) )
    print "The first entry in 'avg_seg_arr' = %f." %( avg_seg_arr[0] )

    df['trkSegRadLengthsMin'] = min_seg_arr        
    df['trkSegRadLengthsAvg'] = avg_seg_arr        

    print df['trkSegRadLengthsAvg'][0]

    if 'simMom' in df.columns.values:
        df['simMom_inverse'] = 1./df['simMom']

    df['fracDiff'] = (df['trkMom_Mu_inverse']-df['simMom_inverse']) / df['simMom_inverse']

#    print df.columns.values

    return df
"""
    #This df has segment-by-segment deviation (scattering angle, etc)
    segdf = pd.DataFrame(  root2array ( myfile, 'TMC_debug_tree' ) )
    segdf['dthetayoverpredictedRMS'] = segdf['delta_theta_y']/segdf['predicted_RMS']
    segdf['dthetayovertruepredictedRMS'] = segdf['delta_theta_y']/segdf['true_predicted_RMS']
    segdf['dthetayoverpredictedRMS_fromMCS'] = segdf['delta_theta_y']/segdf['predicted_RMS_fromMCS']
    
    segdf['dthetaxoverpredictedRMS'] = segdf['delta_theta_x']/segdf['predicted_RMS']
    segdf['dthetaxovertruepredictedRMS'] = segdf['delta_theta_x']/segdf['true_predicted_RMS']
    segdf['dthetaxoverpredictedRMS_fromMCS'] = segdf['delta_theta_x']/segdf['predicted_RMS_fromMCS']
 
    
    #Optional driver DF tree that holds some MCTrack informationOA
    driverdf = pd.DataFrame( root2array ( myfile, 'driver_tree' ) )
    
    #Merge it into the main df by run,subrun,eventid
    df = df.merge(driverdf, on=['run','subrun','eventid'])
    #also merge into segdf
    segdf = segdf.merge(driverdf, on=['run','subrun','eventid'])

    return df, segdf
"""
