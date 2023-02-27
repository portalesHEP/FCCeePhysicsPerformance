#include "LLRHiggsTools.h"
using namespace LLRHiggsTools;




ROOT::VecOps::RVec<fastjet::PseudoJet> LLRHiggsTools::get_compositeJets(ROOT::VecOps::RVec<fastjet::PseudoJet> in){ //MUSS schon z -> qq sein!!!
  ROOT::VecOps::RVec<fastjet::PseudoJet> result; //Z-> ij: Immer gerade anzahl jets
  for (size_t i = 0; i < in.size() - 1; i++) { //i++ erhöht i nach ausführung
    result.push_back(fastjet::join(in[i], in[i+1]));
    i += 1; 
  }
  return result;
}


LLRHiggsTools::jetResonanceBuilder::jetResonanceBuilder(float arg_resonance_mass) {m_resonance_mass = arg_resonance_mass;}
ROOT::VecOps::RVec<fastjet::PseudoJet> LLRHiggsTools::jetResonanceBuilder::operator()(ROOT::VecOps::RVec<fastjet::PseudoJet> in) {
  ROOT::VecOps::RVec<fastjet::PseudoJet> result;
  int n = in.size();
  if (n >1) {
    ROOT::VecOps::RVec<bool> v(n);
    std::fill(v.end() - 2, v.end(), true);
    do {
      fastjet::PseudoJet reso;
      fastjet::PseudoJet reso_lv; 
      for (int i = 0; i < n; ++i) {
          if (v[i]) {
            //fastjet::PseudoJet leg_lv = fastjet::PseudoJet(in[i].px(), in[i].py(), in[i].pz(), in[i].e());
            //leg_lv.SetXYZM(in[i].momentum.x, in[i].momentum.y, in[i].momentum.z, in[i].mass);
            //reso_lv += leg_lv;
            reso_lv = fastjet::join(in[i], reso_lv);
          }
      }
      //reso = fastjet::PseudoJet(reso_lv.px(), reso_lv.py(), reso_lv.pz(), reso_lv.e());
      reso = fastjet::join(reso_lv);
      result.emplace_back(reso);
    } while (std::next_permutation(v.begin(), v.end()));
  }
  if (result.size() > 1) {
    auto resonancesort = [&] (fastjet::PseudoJet i ,fastjet::PseudoJet j) { return (abs( m_resonance_mass -i.m())<abs(m_resonance_mass-j.m())); };
    std::sort(result.begin(), result.end(), resonancesort);
    ROOT::VecOps::RVec<fastjet::PseudoJet>::const_iterator first = result.begin();
    ROOT::VecOps::RVec<fastjet::PseudoJet>::const_iterator last = result.begin() + 1;
    ROOT::VecOps::RVec<fastjet::PseudoJet> onlyBestReso(first, last);
    return onlyBestReso;
  } else {
    return result;
  }
}

ROOT::VecOps::RVec<fastjet::PseudoJet> LLRHiggsTools::get_compositeSubJets(ROOT::VecOps::RVec<fastjet::PseudoJet> in){ //MUSS schon z -> qq sein!!!
  ROOT::VecOps::RVec<fastjet::PseudoJet> result; //Z-> ij: Immer gerade anzahl jets
  //std::cout << in.size();
  if (in.size() == 1) {
    //fastjet::PseudoJet tmp;
    result = (((in.at(0)).pieces()).at(0)).pieces(); //unpack the dijet
    //std::cout << result.size();
    return result;
  } else {
    return in;
  }
}

ROOT::VecOps::RVec<float> LLRHiggsTools::get_subjetpl(int index, ROOT::VecOps::RVec<fastjet::PseudoJet> in){
  ROOT::VecOps::RVec<float> result;
  if (in.size() == 2) {
    result.push_back(((in.at(index)).pz())); //nutze push_back weil wir einen vektor haben ROOT::VecOps::RVec
    return result;
  } else {
    for (auto & p: in) {
      //result.push_back(p.pt()); //nutze push_back weil wir einen vektor haben ROOT::VecOps::RVec
      result.push_back(-12);
  }
  return result;
  }
}

ROOT::VecOps::RVec<float> LLRHiggsTools::get_subjetpt(int index, ROOT::VecOps::RVec<fastjet::PseudoJet> in){
  ROOT::VecOps::RVec<float> result;
  if (in.size() == 2) {
    result.push_back(((in.at(index)).pt())); //nutze push_back weil wir einen vektor haben ROOT::VecOps::RVec
    return result;
  } else {
    for (auto & p: in) {
      //result.push_back(p.pt()); //nutze push_back weil wir einen vektor haben ROOT::VecOps::RVec
      result.push_back(-12);
  }
  return result;
  }
}

ROOT::VecOps::RVec<float> LLRHiggsTools::get_subjetE(int index, ROOT::VecOps::RVec<fastjet::PseudoJet> in){
  ROOT::VecOps::RVec<float> result;
  if (in.size() == 2) {
    result.push_back(((in.at(index)).E())); //nutze push_back weil wir einen vektor haben ROOT::VecOps::RVec
    return result;
  } else {
    for (auto & p: in) {
      //result.push_back(p.E()); //nutze push_back weil wir einen vektor haben ROOT::VecOps::RVec
      result.push_back(-12);
  }
  return result;
  }
}

ROOT::VecOps::RVec<float> LLRHiggsTools::get_subjetEta(int index, ROOT::VecOps::RVec<fastjet::PseudoJet> in){
  ROOT::VecOps::RVec<float> result;
  if (in.size() == 2) {
    result.push_back(((in.at(index)).eta())); //nutze push_back weil wir einen vektor haben ROOT::VecOps::RVec
    return result;
  } else {
    for (auto & p: in) {
      //result.push_back(p.eta()); //nutze push_back weil wir einen vektor haben ROOT::VecOps::RVec
      result.push_back(-120);
  }
  return result;
  }
}


ROOT::VecOps::RVec<float> LLRHiggsTools::get_subjetDeltaEta(ROOT::VecOps::RVec<fastjet::PseudoJet> in){
  ROOT::VecOps::RVec<float> result;
  if (in.size() == 2) {
    //std::cout << abs(((in.at(0)).theta() - (in.at(1)).theta()));
    //result.push_back(fabs(((in.at(0)).theta() - (in.at(1)).theta())));
    result.push_back(fabs(((in.at(0)).eta() - (in.at(1)).eta())));
    return result;
  } else {
    for (auto & p: in) {
      //result.push_back(p.theta()); //nutze push_back weil wir einen vektor haben ROOT::VecOps::RVec
      result.push_back(-12); //irgendeine zahl ausserhalb des betrachteten intervalls
  }
  return result;
  }
}

ROOT::VecOps::RVec<float> LLRHiggsTools::get_subjetDeltaPhi(ROOT::VecOps::RVec<fastjet::PseudoJet> in){
  ROOT::VecOps::RVec<float> result;
  if (in.size() == 2) {
    //std::cout << abs(((in.at(0)).theta() - (in.at(1)).theta()));
    //result.push_back(fabs(((in.at(0)).theta() - (in.at(1)).theta())));
    result.push_back(fabs(((in.at(0)).phi() - (in.at(1)).phi())));
    return result;
  } else {
    for (auto & p: in) {
      //result.push_back(p.theta()); //nutze push_back weil wir einen vektor haben ROOT::VecOps::RVec
      result.push_back(-12); //irgendeine zahl ausserhalb des betrachteten intervalls
  }
  return result;
  }
}


ROOT::VecOps::RVec<float> LLRHiggsTools::get_subjetAcoll(ROOT::VecOps::RVec<fastjet::PseudoJet> in){
  ROOT::VecOps::RVec<float> result;
  if (in.size() == 2) {
    //std::cout << abs(((in.at(0)).theta() - (in.at(1)).theta()));
    //result.push_back(fabs(((in.at(0)).theta() - (in.at(1)).theta())));
    result.push_back(M_PI - fabs(((in.at(0)).theta() + (in.at(1)).theta()))); //M_PI ist Pi
    return result;
  } else {
    for (auto & p: in) {
      //result.push_back(p.theta()); //nutze push_back weil wir einen vektor haben ROOT::VecOps::RVec
      result.push_back(-12); //irgendeine zahl ausserhalb des betrachteten intervalls
  }
  return result;
  }
}

ROOT::VecOps::RVec<float> LLRHiggsTools::get_subjetAcollPaPb(ROOT::VecOps::RVec<fastjet::PseudoJet> in){
  ROOT::VecOps::RVec<float> result;
  if (in.size() == 2) {
    //std::cout << abs(((in.at(0)).theta() - (in.at(1)).theta()));
    //result.push_back(fabs(((in.at(0)).theta() - (in.at(1)).theta())));
    float spPaPb;
    float absPa;
    float absPb;
    spPaPb = ((in.at(0)).px() * (in.at(1)).px()) + ((in.at(0)).py() * (in.at(1)).py()) + ((in.at(0)).pz() * (in.at(1)).pz());
    absPa = std::sqrt(((in.at(0)).px() * (in.at(0)).px()) + ((in.at(0)).py() * (in.at(0)).py()) + ((in.at(0)).pz() * (in.at(0)).pz()));
    absPb = std::sqrt(((in.at(1)).px() * (in.at(1)).px()) + ((in.at(1)).py() * (in.at(1)).py()) + ((in.at(1)).pz() * (in.at(1)).pz()));
    result.push_back(std::acos(spPaPb / (absPa * absPb)));
    return result;
  } else {
    for (auto & p: in) {
      //result.push_back(p.theta()); //nutze push_back weil wir einen vektor haben ROOT::VecOps::RVec
      result.push_back(-12); //irgendeine zahl ausserhalb des betrachteten intervalls
  }
  return result;
  }
}



ROOT::VecOps::RVec<float> LLRHiggsTools::get_jetNo(ROOT::VecOps::RVec<fastjet::PseudoJet> in){
  ROOT::VecOps::RVec<float> result;
  result.push_back(in.size());
  return result;
}







ROOT::VecOps::RVec<float> LLRHiggsTools::get_m_rec(float ecm, ROOT::VecOps::RVec<fastjet::PseudoJet> in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(sqrt((pow(ecm, 2)) - (2*p.e()*ecm) + (pow(p.m(), 2))));
  }
  return result;
}

LLRHiggsTools::jetRecoilBuilder::jetRecoilBuilder(float arg_sqrts)  {sqrts = arg_sqrts;}
ROOT::VecOps::RVec<fastjet::PseudoJet> LLRHiggsTools::jetRecoilBuilder::operator() (ROOT::VecOps::RVec<fastjet::PseudoJet> in) {
  ROOT::VecOps::RVec<fastjet::PseudoJet> result;
  auto recoil_p4 = fastjet::PseudoJet(0, 0, 0, sqrts);
  //auto recoil_p4 = TLorentzVector(0, 0, 0, m_sqrts);
  for (auto & v1: in) {
    fastjet::PseudoJet tv1 = fastjet::PseudoJet(v1.px(), v1.py(), v1.pz(), v1.e());
    recoil_p4 -= tv1;
  }
  //auto recoil_fcc = edm4hep::ReconstructedParticleData();
  //recoil_fcc.momentum.x = recoil_p4.Px();
  //recoil_fcc.momentum.y = recoil_p4.Py();
  //recoil_fcc.momentum.z = recoil_p4.Pz();
  //recoil_fcc.mass = recoil_p4.M();
  result.push_back(recoil_p4);
  return result;
};

ROOT::VecOps::RVec<float> LLRHiggsTools::get_ht(ROOT::VecOps::RVec<fastjet::PseudoJet> in){
  ROOT::VecOps::RVec<float> result;
  float temp = 0.0;
  //for (size_t i = 0; i < in.size(); ++i) {
  //  temp += in[i].pt()
  //}
  for (auto & p: in) {
    temp += p.pt();
  }
  result.push_back(temp); //nutze push_back weil wir einen vektor haben ROOT::VecOps::RVec
  return result;
}

// --------------------------------------------------------------------------------------------------
// Courtesy of the Flavor group

// To retrieve a given MC leg corresponding to the Bs decay

selMC_leg::selMC_leg( int idx ) {
	m_idx = idx;
};

   // I return a vector instead of a single particle :
   //   - such that the vector is empty when there is no such decay mode (instead
   //     of returning a dummy particle)
   //   - such that I can use the getMC_theta etc functions, which work with a
   //     ROOT::VecOps::RVec of particles, and not a single particle

ROOT::VecOps::RVec<edm4hep::MCParticleData> selMC_leg::operator() ( ROOT::VecOps::RVec<int> list_of_indices,  ROOT::VecOps::RVec<edm4hep::MCParticleData> in) {
	ROOT::VecOps::RVec<edm4hep::MCParticleData>  res;
	if ( list_of_indices.size() == 0) return res;
	if ( m_idx < list_of_indices.size() ) {
		res.push_back( FCCAnalyses::MCParticle::sel_byIndex( list_of_indices[m_idx], in ) );
		return res;
	 }
	else {
		std::cout << "   !!!  in selMC_leg:  idx = " << m_idx << " but size of list_of_indices = " << list_of_indices.size() << std::endl;
	}
	return res;
}


// --------------------------------------------------------------------------------------------------
// Courtesy of E. Perez
int LLRHiggsTools::get_lepton_origin( edm4hep::MCParticleData p, ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> ind)  {

 // std::cout  << std::endl << " enter in MCParticle::get_lepton_origin  PDG = " << p.PDG << std::endl;

 int pdg = std::abs( p.PDG ) ; // Nur ein teilchen ist ausgewählt
 if ( pdg != 11 && pdg != 13 && pdg  != 15 && pdg  != 12 &&  pdg != 14 && pdg != 16 ) return -1; // Wollen nur leptonen (charged + neutrinos) betrachten, alle anderen sind egal!

 int result  = 0;

 //std::cout << " p.parents_begin p.parents_end " << p.parents_begin <<  " "  << p.parents_end << std::endl;
    for (unsigned j = p.parents_begin; j != p.parents_end; ++j) { //Loope über die parents der leptonen
      int index = ind.at(j);
      int pdg_parent = in.at(index).PDG ;
      //std::cout  << " parent has pdg = " << in.at(index).PDG <<  "  status = " << in.at(index).generatorStatus << std::endl;

      if ( abs( pdg_parent ) == 23 || abs( pdg_parent ) == 24 ) {
        result = pdg_parent ;
        //std::cout <<  " ... Lepton is from W or Z ,  return code = " << result <<  std::endl;
        break; //Unterbreche den Loop, wenn parent ein W oder Z ist
      }

      if ( abs( pdg_parent ) == 22 ) {
        result = pdg_parent ;
        //std::cout <<  " ... Lepton is from a virtual photon ,  return code = " << result <<  std::endl;
        break;
      }

      if ( abs( pdg_parent ) == 15 ) {
         result = pdg_parent ;
         //std::cout <<  " ... Lepton is from a tau,  return code = " << result <<  std::endl;
         break;
      }

      if ( abs( pdg_parent ) == 11 ) {    // beam particle ?
			// beam particles should have generatorStatus = 4,
			// but that is not the case in files produced from Whizard + p6
        if ( in.at(index).generatorStatus == 4 || ind.at  ( in.at(index).parents_begin ) == 0 ) { //zuerst oder, .parents_begin = erstes teilchen im decay prozess?
           result = 0; //Ist egal, denn wir suchen am schluss nach -11 und betrachten wzp6 sample
           //std::cout <<  " ... Lepton is from the hard subprocess, return code = " << result <<  std::endl;
           break;
        }
      }

      if ( pdg == 11 && abs( pdg_parent ) == 13 ) {	// mu -> e
          result  = pdg_parent;
          //std::cout <<  " ... Electron from a muon decay, return code = " << result <<  std::endl;
          break;
      }

      if ( abs( pdg_parent ) == pdg  ) {
	//std::cout << " ... iterate ... " << std::endl;
	return get_lepton_origin( in.at(index),  in, ind  ); //Do loop if parent is no W, Z, mu, gamma or beam particle => has to be coming from hadrons
      }

      // This must come from a hadron decay
      result = pdg_parent;
      //std::cout <<  " ... Lepton from a hadron decay " << std::endl;

    }
 
 return result;

}


int LLRHiggsTools::get_lepton_origin( int index, ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> ind)  {
  if ( index < 0 || index >= in.size() ) return -1;
  edm4hep::MCParticleData p = in[index];
  return get_lepton_origin( p, in, ind );
}


ROOT::VecOps::RVec<int> LLRHiggsTools::get_leptons_origin( ROOT::VecOps::RVec<edm4hep::MCParticleData> particles, ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> ind)  {

  ROOT::VecOps::RVec<int> result;
  result.reserve(particles.size());
  for (size_t i = 0; i < particles.size(); ++i) {
    auto & p = particles[i]; //betrachte ein einzelnes mc teilchen 
    int origin = LLRHiggsTools::get_lepton_origin( p, in, ind );
    //std::cout << origin << "\n";
    result.push_back( origin );
  }

  return result;

}

// Courtesy of E. Perez
ROOT::VecOps::RVec<edm4hep::MCParticleData> LLRHiggsTools::sel_leptons_origin( ROOT::VecOps::RVec<edm4hep::MCParticleData> particles, ROOT::VecOps::RVec<int> flags_origin, int code ) {
  /// flags_origin and particles should have the same size. 

 ROOT::VecOps::RVec<edm4hep::MCParticleData>  result;
 for (size_t i = 0; i < particles.size(); ++i) {
    auto & p = particles[i];
    int origin = flags_origin[i] ;
    if ( origin == code ) result.push_back( p ) ;
 }
 return result;
}
