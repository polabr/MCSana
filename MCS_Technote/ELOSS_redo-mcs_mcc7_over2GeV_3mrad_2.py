import ROOT
import math, sys
from array import array

pMin_ = 0.01
pMax_ = 7.50
pStep_ = 0.01
eLossMode_ = 2 #Change to 2 for BB, 1 for constant (MIP); 0 is Landau MPV
nElossSteps_ = 10
lar_radl = 14.0
angResol_ = 3.0

file = ROOT.TFile.Open("anafiles/mcc7_bnb_nuOnly.root")
dir = "trajmcsntupleNu"

def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step

def mass():
    return 0.105658367

def MomentumDependentConstant(p):
    #these are from https://arxiv.org/abs/1703.06187
    a = 0.1049
    c = 11.0038
    return (a/(p*p)) + c

def energyLossLandau(mass2, e2, x):
    #
    # eq. (33.11) in http://pdg.lbl.gov/2016/reviews/rpp2016-rev-passage-particles-matter.pdf (except density correction is ignored)
    #
    if x<=0.: return 0.
    Iinv2 = 1./(188.E-6*188.E-6)
    matConst = 1.4*18./40. #density*Z/A
    me = 0.511
    kappa = 0.307075
    j = 0.200
    #
    beta2 = (e2-mass2)/e2
    gamma2 = 1./(1.0 - beta2)
    epsilon = 0.5*kappa*x*matConst/beta2
    #
    return 0.001*epsilon*( math.log(2.*me*beta2*gamma2*epsilon*Iinv2) + j - beta2 )

def energyLossBetheBloch(mass, e2):
    # stolen, mostly, from GFMaterialEffects.
    Iinv = 1./188.E-6
    matConst = 1.4*18./40. #density*Z/A
    me = 0.511
    kappa = 0.307075
    #
    beta2 = (e2-mass*mass)/e2
    gamma2 = 1./(1.0 - beta2)
    massRatio = me/mass
    argument = (2.*me*gamma2*beta2*Iinv) * math.sqrt(1+2*math.sqrt(gamma2)*massRatio + massRatio*massRatio)
    #
    dedx = kappa*matConst/beta2
    #
    if math.fabs(mass)<1e-10 : return 0.0
    if argument <= math.exp(beta2):
        dedx = 0.
    else:
        dedx *= (math.log(argument)-beta2)*1.E-3 # Bethe-Bloch, converted to GeV/cm
        if dedx<0.: dedx = 0.
    return dedx

def GetE(initial_E, length_travelled, m):
    #
    step_size = length_travelled / nElossSteps_
    #
    current_E = initial_E
    m2 = m*m
    #
    for i in range (0, nElossSteps_):
        if eLossMode_==2:
            dedx = energyLossBetheBloch(m,current_E)
            current_E = current_E-(dedx * step_size)
        else:
            #MPV of Landau energy loss distribution
            current_E = current_E-energyLossLandau(m2,current_E*current_E,step_size)
        if current_E <= m :
            #std::cout<<"WARNING: current_E less than mu mass. it is "<<current_E<<std::endl;
            return 0.
    return current_E

def mcsLikelihood(p, theta0x, dthetaij, seg_nradl, cumLen, fwd, momDepConst):
    beg  = 0 if fwd else (len(dthetaij)-1)
    end  = len(dthetaij) if fwd else -1
    incr = +1 if fwd else -1
    #
    # print = false;//(p>1.999 && p<2.001);
    #
    m = mass()
    m2 = m*m
    Etot = math.sqrt(p*p + m2) #Initial energy
    Eij2 = 0.
    #
    fixedterm = 0.5 * math.log( 2.0 * math.pi )
    result = 0
    for i in range (beg, end, incr):
        if dthetaij[i]<0:
            #cout << "skip segment with too few points" << endl;
            continue
        #
        if eLossMode_==1 :
            # ELoss mode: MIP (constant)
            kcal = 0.002105
            Eij = Etot - kcal*cumLen[i] #energy at this segment
            Eij2 = Eij*Eij
        else:
            # Non constant energy loss distribution
            Eij = GetE(Etot,cumLen[i],m)
            Eij2 = Eij*Eij
        #
        if Eij2 <= m2 :
            result = sys.float_info.max
            break
        pij = math.sqrt(Eij2 - m2) #momentum at this segment
        beta = math.sqrt( 1. - ((m2)/(pij*pij + m2)) )
        tuned_HL_term1 = 11.0038 # https://arxiv.org/abs/1703.06187
        HL_term2 = 0.038
        mdc = (MomentumDependentConstant(pij) if momDepConst else tuned_HL_term1)
        #print seg_nradl[i]
        tH0 = ( mdc / (pij*beta) ) * ( 1.0 + HL_term2 * math.log( seg_nradl[i] ) ) * math.sqrt( seg_nradl[i] )
        rms = math.sqrt( 2.0*( tH0 * tH0 + theta0x * theta0x ) )
        if math.fabs(rms)<1e-10 :
            print " Error : RMS cannot be zero ! "
            return sys.float_info.max
        arg = dthetaij[i]/rms
        result = result + ( math.log( rms ) + 0.5 * arg * arg + fixedterm)
        #if (print && fwd==true) cout << "TrajectoryMCSFitter pij=" << pij << " dthetaij[i]=" << dthetaij[i] << " tH0=" << tH0 << " rms=" << rms << " prob=" << ( std::log( rms ) + 0.5 * arg * arg + fixedterm) << " const=" << (momDepConst ? MomentumDependentConstant(pij) : tuned_HL_term1) << " beta=" << beta << " red_length=" << seg_nradl[i] << endl;
    return result

def doLikelihoodScan(dtheta, seg_nradlengths, cumLen, fwdFit, momDepConst):
    best_idx  = -1
    best_logL = sys.float_info.max
    best_p    = -1.0
    vlogL = []
    for p_test in drange (pMin_, pMax_, pStep_):
        logL = mcsLikelihood(p_test, angResol_, dtheta, seg_nradlengths, cumLen, fwdFit, momDepConst)
        if logL < best_logL :
            best_p    = p_test
            best_logL = logL
            best_idx  = len(vlogL)
        vlogL.append(logL);
    #
    #uncertainty from left side scan
    lunc = -1.0
    if best_idx>0:
        for j in range(best_idx-1,-1,-1):
            dLL = vlogL[j]-vlogL[best_idx]
            if dLL<0.5:
                lunc = (best_idx-j)*pStep_
            else: break
    #uncertainty from right side scan
    runc = -1.0
    if best_idx<(len(vlogL)-1):
        for j in range(best_idx+1,len(vlogL)):
            dLL = vlogL[j]-vlogL[best_idx]
            if dLL<0.5:
                runc = (j-best_idx)*pStep_
            else: break;
    best_p_err = max(lunc,runc)
    return [best_p, best_p_err, best_logL]

tree = file.Get(dir+"/tree")

filemod = ROOT.TFile.Open("ELOSS_mcc7_bnb_nuOnly_over2GeV_3mrad_eloss2.root","RECREATE")
filemod.mkdir(dir)
filemod.cd(dir)
treemod = tree.CloneTree(0)

trkSegRadLengths = ROOT.vector('double')()
tree.SetBranchAddress("trkSegRadLengths", trkSegRadLengths)

trkScattAngles = ROOT.vector('double')()
tree.SetBranchAddress("trkScattAngles", trkScattAngles)

trkMom_MuFwd_tmp = array ( 'd', [0] )
treemod.SetBranchAddress( "trkMom_MuFwd", trkMom_MuFwd_tmp )
trkMomErr_MuFwd_tmp = array ( 'd', [0] )
treemod.SetBranchAddress( "trkMomErr_MuFwd", trkMomErr_MuFwd_tmp )
trkLL_MuFwd_tmp = array ( 'd', [0] )
treemod.SetBranchAddress( "trkLL_MuFwd", trkLL_MuFwd_tmp )
trkMom_MuBwd_tmp = array ( 'd', [0] )
treemod.SetBranchAddress( "trkMom_MuBwd", trkMom_MuBwd_tmp )
trkMomErr_MuBwd_tmp = array ( 'd', [0] )
treemod.SetBranchAddress( "trkMomErr_MuBwd", trkMomErr_MuBwd_tmp )
trkLL_MuBwd_tmp = array ( 'd', [0] )
treemod.SetBranchAddress( "trkLL_MuBwd", trkLL_MuBwd_tmp )
trkMom_Mu_tmp = array ( 'd', [0] )
treemod.SetBranchAddress( "trkMom_Mu", trkMom_Mu_tmp )
trkMomErr_Mu_tmp = array ( 'd', [0] )
treemod.SetBranchAddress( "trkMomErr_Mu", trkMomErr_Mu_tmp )
trkLL_Mu_tmp = array ( 'd', [0] )
treemod.SetBranchAddress( "trkLL_Mu", trkLL_Mu_tmp )

count = 0
total = tree.GetEntries()

for t in tree: 
    count = count + 1
    if (count % 10000 == 0): 
        print count, total
#if count > 1: break
    if tree.simMom < 2 or tree.trkLength < 100 or tree.simID!=13: 
        continue
    trkLL_MuFwd = tree.trkLL_MuFwd
    trkMom_MuFwd = tree.trkMom_MuFwd
    trkMomErr_MuFwd = tree.trkMomErr_MuFwd
    trkLL_MuBwd = tree.trkLL_MuBwd
    trkMom_MuBwd = tree.trkMom_MuBwd
    trkMomErr_MuBwd = tree.trkMomErr_MuBwd
    trkLL_Mu = tree.trkLL_Mu
    trkMom_Mu = tree.trkMom_Mu
    trkMomErr_Mu = tree.trkMomErr_Mu
    #
    cumseglens = []
    cumseglens.append(0)
    for i in range (0, len(trkSegRadLengths)):
        cumseglens.append(cumseglens[-1]+trkSegRadLengths[i]*lar_radl)
    cumLenFwd = []
    cumLenBwd = []
    for i in range (0, len(cumseglens)-2) :
        cumLenFwd.append( cumseglens[i] );
        cumLenBwd.append( cumseglens[-1]-cumseglens[i+2] );
    fwdresult = doLikelihoodScan(trkScattAngles, trkSegRadLengths, cumLenFwd, True, True)
    bwdresult = doLikelihoodScan(trkScattAngles, trkSegRadLengths, cumLenBwd, False, True)
    bestresult = fwdresult if fwdresult[2]<bwdresult[2] else bwdresult
    #
    #print trkMom_MuFwd, trkMomErr_MuFwd, trkLL_MuFwd, fwdresult
    #print trkMom_MuBwd, trkMomErr_MuBwd, trkLL_MuBwd, bwdresult
    #print trkMom_Mu, trkMomErr_Mu, trkLL_Mu, bestresult
    trkMom_MuFwd_tmp[0] = fwdresult[0]
    trkMomErr_MuFwd_tmp[0] = fwdresult[1]
    trkLL_MuFwd_tmp[0] = fwdresult[2]
    trkMom_MuBwd_tmp[0] = bwdresult[0]
    trkMomErr_MuBwd_tmp[0] = bwdresult[1]
    trkLL_MuBwd_tmp[0] = bwdresult[2]
    trkMom_Mu_tmp[0] = bestresult[0]
    trkMomErr_Mu_tmp[0] = bestresult[1]
    trkLL_Mu_tmp[0] = bestresult[2]
    treemod.Fill()

filemod.Write()
filemod.Close()
