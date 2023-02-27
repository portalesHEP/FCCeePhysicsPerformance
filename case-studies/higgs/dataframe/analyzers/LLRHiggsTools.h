
#ifndef  LLRHiggsTools_ANALYZERS_H
#define  LLRHiggsTools_ANALYZERS_H

#include <cmath>
#include <vector>

#include "TLorentzVector.h"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/ParticleIDData.h"
#include "MCParticle.h"

#include "ReconstructedParticle2MC.h"

#include "Math/Vector4D.h"
#include "ROOT/RVec.hxx"
#include "fastjet/JetDefinition.hh"
#include "TRandom3.h"

namespace LLRHiggsTools{

    ROOT::VecOps::RVec<fastjet::PseudoJet> get_compositeJets(ROOT::VecOps::RVec<fastjet::PseudoJet> in);

    /// build the resonance from 2 particles from an arbitrary list of input ReconstructedPartilces. Keep the closest to the mass given as input
    struct jetResonanceBuilder {
        float m_resonance_mass;
        jetResonanceBuilder(float arg_resonance_mass);
        ROOT::VecOps::RVec<fastjet::PseudoJet> operator()(ROOT::VecOps::RVec<fastjet::PseudoJet> in);
    };

    ROOT::VecOps::RVec<fastjet::PseudoJet> get_compositeSubJets(ROOT::VecOps::RVec<fastjet::PseudoJet> in);

    ROOT::VecOps::RVec<float> get_subjetpl(int index, ROOT::VecOps::RVec<fastjet::PseudoJet> in);
    ROOT::VecOps::RVec<float> get_subjetpt(int index, ROOT::VecOps::RVec<fastjet::PseudoJet> in);
    ROOT::VecOps::RVec<float> get_subjetE(int index, ROOT::VecOps::RVec<fastjet::PseudoJet> in);
    ROOT::VecOps::RVec<float> get_subjetEta(int index, ROOT::VecOps::RVec<fastjet::PseudoJet> in);
    ROOT::VecOps::RVec<float> get_subjetDeltaEta(ROOT::VecOps::RVec<fastjet::PseudoJet> in);
    ROOT::VecOps::RVec<float> get_subjetDeltaPhi(ROOT::VecOps::RVec<fastjet::PseudoJet> in);
    ROOT::VecOps::RVec<float> get_subjetAcoll(ROOT::VecOps::RVec<fastjet::PseudoJet> in);
    ROOT::VecOps::RVec<float> get_subjetAcollPaPb(ROOT::VecOps::RVec<fastjet::PseudoJet> in);
    ROOT::VecOps::RVec<float> get_jetNo(ROOT::VecOps::RVec<fastjet::PseudoJet> in);

    ROOT::VecOps::RVec<float> get_m_rec(float ecm, ROOT::VecOps::RVec<fastjet::PseudoJet> in);

    struct jetRecoilBuilder {
     jetRecoilBuilder(float arg_sqrts);
     float sqrts = 240.0;
     ROOT::VecOps::RVec<fastjet::PseudoJet> operator()(ROOT::VecOps::RVec<fastjet::PseudoJet> in) ;
    };

    ROOT::VecOps::RVec<float> get_ht(ROOT::VecOps::RVec<fastjet::PseudoJet> in);

    struct selMC_leg{
     selMC_leg( int idx );
     int m_idx;
     ROOT::VecOps::RVec<edm4hep::MCParticleData> operator() (ROOT::VecOps::RVec<int> list_of_indices,
                                                             ROOT::VecOps::RVec<edm4hep::MCParticleData> in) ;
    };


    /// return the pdg ID of the parent of a lepton (pre-FSR)
    int get_lepton_origin( int idx, ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> ind) ;
    int get_lepton_origin( edm4hep::MCParticleData p, ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> ind) ;

    ROOT::VecOps::RVec<int> get_leptons_origin( ROOT::VecOps::RVec<edm4hep::MCParticleData> particles, ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> ind) ;

    ROOT::VecOps::RVec<edm4hep::MCParticleData> sel_leptons_origin( ROOT::VecOps::RVec<edm4hep::MCParticleData> particles, ROOT::VecOps::RVec<int> flags_origin, int code ) ;



}
#endif
