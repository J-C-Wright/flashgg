#ifndef FLASHgg_VBFTagTruth_h
#define FLASHgg_VBFTagTruth_h

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"

namespace flashgg {

    class VBFTagTruth : public TagTruthBase
    {

    public:

        VBFTagTruth();
        ~VBFTagTruth();
        //        VBFTagTruth(const VBFTagTruth &b);

        // Functions for dumping
        float pt_genJetMatchingToJ1() const { return ( hasClosestGenJetToLeadingJet() ? closestGenJetToLeadingJet()->pt() : -1. ); }
        float eta_genJetMatchingToJ1() const { return ( hasClosestGenJetToLeadingJet() ? closestGenJetToLeadingJet()->eta() : -999. );}
        float phi_genJetMatchingToJ1() const { return ( hasClosestGenJetToLeadingJet() ? closestGenJetToLeadingJet()->phi() : -999. );}
        float pt_genJetMatchingToJ2() const { return ( hasClosestGenJetToSubLeadingJet() ? closestGenJetToSubLeadingJet()->pt() : -1. );}
        float eta_genJetMatchingToJ2() const { return ( hasClosestGenJetToSubLeadingJet() ? closestGenJetToSubLeadingJet()->eta() : -999. );}
        float phi_genJetMatchingToJ2() const { return ( hasClosestGenJetToSubLeadingJet() ? closestGenJetToSubLeadingJet()->phi() : -999. );}
        float pt_genJetMatchingToJ3() const { return ( hasClosestGenJetToSubSubLeadingJet() ? closestGenJetToSubSubLeadingJet()->pt() : -1. );}
        float eta_genJetMatchingToJ3() const { return ( hasClosestGenJetToSubSubLeadingJet() ? closestGenJetToSubSubLeadingJet()->eta() : -999. );}
        float phi_genJetMatchingToJ3() const { return ( hasClosestGenJetToSubSubLeadingJet() ? closestGenJetToSubSubLeadingJet()->phi() : -999. );}
        float pt_genPartMatchingToJ1() const { return ( hasClosestParticleToLeadingJet() ? closestParticleToLeadingJet()->pt() : -1. );}
        float eta_genPartMatchingToJ1() const { return ( hasClosestParticleToLeadingJet() ? closestParticleToLeadingJet()->eta() : -999. );}
        float phi_genPartMatchingToJ1() const { return ( hasClosestParticleToLeadingJet() ? closestParticleToLeadingJet()->phi() : -999. );}
        float pt_genPartMatchingToJ2() const { return ( hasClosestParticleToSubLeadingJet() ? closestParticleToSubLeadingJet()->pt() : -1. );}
        float eta_genPartMatchingToJ2() const { return ( hasClosestParticleToSubLeadingJet() ? closestParticleToSubLeadingJet()->eta() : -999. );}
        float phi_genPartMatchingToJ2() const { return ( hasClosestParticleToSubLeadingJet() ? closestParticleToSubLeadingJet()->phi() : -999. );}
        float pt_genPartMatchingToJ3() const { return ( hasClosestParticleToSubSubLeadingJet() ? closestParticleToSubSubLeadingJet()->pt() : -1. );}
        float eta_genPartMatchingToJ3() const { return ( hasClosestParticleToSubSubLeadingJet() ? closestParticleToSubSubLeadingJet()->eta() : -999. );}
        float phi_genPartMatchingToJ3() const { return ( hasClosestParticleToSubSubLeadingJet() ? closestParticleToSubSubLeadingJet()->phi() : -999. );}
        float pt_genPartMatchingToPho1() const { return ( hasClosestParticleToLeadingPhoton() ? closestParticleToLeadingPhoton()->pt() : -1. );}
        float eta_genPartMatchingToPho1() const { return ( hasClosestParticleToLeadingPhoton() ? closestParticleToLeadingPhoton()->eta() : -999. );}
        float phi_genPartMatchingToPho1() const { return ( hasClosestParticleToLeadingPhoton() ? closestParticleToLeadingPhoton()->phi() : -999. );}
        float pt_genPartMatchingToPho2() const { return ( hasClosestParticleToSubLeadingPhoton() ? closestParticleToSubLeadingPhoton()->pt() : -1. );}
        float eta_genPartMatchingToPho2() const { return ( hasClosestParticleToSubLeadingPhoton() ? closestParticleToSubLeadingPhoton()->eta() : -999. );}
        float phi_genPartMatchingToPho2() const { return ( hasClosestParticleToSubLeadingPhoton() ? closestParticleToSubLeadingPhoton()->phi() : -999. );}
        float pt_P1() const { return ( hasLeadingParton() ? leadingParton()->pt() : -1. ); }
        float eta_P1() const { return ( hasLeadingParton() ? leadingParton()->eta() : -999. ); }
        float phi_P1() const { return ( hasLeadingParton() ? leadingParton()->phi() : -999. ); }
        float pt_P2() const { return ( hasSubLeadingParton() ? subLeadingParton()->pt() : -1. ); }
        float eta_P2() const { return ( hasSubLeadingParton() ? subLeadingParton()->eta() : -999. ); }
        float phi_P2() const { return ( hasSubLeadingParton() ? subLeadingParton()->phi() : -999. ); }
        float pt_P3() const { return ( hasSubSubLeadingParton() ? subSubLeadingParton()->pt() : -1. ); }
        float eta_P3() const { return ( hasSubSubLeadingParton() ? subSubLeadingParton()->eta() : -999. ); }
        float phi_P3() const { return ( hasSubSubLeadingParton() ? subSubLeadingParton()->phi() : -999. ); }

        //DeltaRs between Jet and truth
        float dR_genJetMatchingToJ1() const { return ( hasClosestGenJetToLeadingJet() ? deltaR(closestGenJetToLeadingJet()->eta(),closestGenJetToLeadingJet()->phi(),
                                                                                               ptOrderedFggJets()[0]->eta(),ptOrderedFggJets()[0]->phi()) : -1. );}
        float dR_genJetMatchingToJ2() const { return ( hasClosestGenJetToLeadingJet() ? deltaR(closestGenJetToSubLeadingJet()->eta(),closestGenJetToSubLeadingJet()->phi(),
                                                                                               ptOrderedFggJets()[1]->eta(),ptOrderedFggJets()[1]->phi()) : -1. );}
        float dR_genJetMatchingToJ3() const { return ( hasClosestGenJetToLeadingJet() ? deltaR(closestGenJetToSubSubLeadingJet()->eta(),closestGenJetToSubSubLeadingJet()->phi(),
                                                                                               ptOrderedFggJets()[2]->eta(),ptOrderedFggJets()[2]->phi()) : -1. );}
        float dR_particleMatchingToJ1() const { return ( hasClosestParticleToLeadingJet() ? deltaR(closestParticleToLeadingJet()->eta(),closestParticleToLeadingJet()->phi(),
                                                                                               ptOrderedFggJets()[0]->eta(),ptOrderedFggJets()[0]->phi()) : -1. );}
        float dR_particleMatchingToJ2() const { return ( hasClosestParticleToLeadingJet() ? deltaR(closestParticleToSubLeadingJet()->eta(),closestParticleToSubLeadingJet()->phi(),
                                                                                               ptOrderedFggJets()[1]->eta(),ptOrderedFggJets()[1]->phi()) : -1. );}
        float dR_particleMatchingToJ3() const { return ( hasClosestParticleToLeadingJet() ? deltaR(closestParticleToSubSubLeadingJet()->eta(),
                                                         closestParticleToSubSubLeadingJet()->phi(), ptOrderedFggJets()[2]->eta(),ptOrderedFggJets()[2]->phi()) : -1. );}
        float dR_partonMatchingToJ1() const { return ( hasClosestPartonToLeadingJet() ? deltaR(closestPartonToLeadingJet()->eta(),closestPartonToLeadingJet()->phi(),
                                                                                               ptOrderedFggJets()[0]->eta(),ptOrderedFggJets()[0]->phi()) : -1. );}
        float dR_partonMatchingToJ2() const { return ( hasClosestPartonToLeadingJet() ? deltaR(closestPartonToSubLeadingJet()->eta(),closestPartonToSubLeadingJet()->phi(),
                                                                                               ptOrderedFggJets()[1]->eta(),ptOrderedFggJets()[1]->phi()) : -1. );}
        float dR_partonMatchingToJ3() const { return ( hasClosestPartonToLeadingJet() ? deltaR(closestPartonToSubSubLeadingJet()->eta(),
                                                         closestPartonToSubSubLeadingJet()->phi(), ptOrderedFggJets()[2]->eta(),ptOrderedFggJets()[2]->phi()) : -1. );}
                  

    //MVA vars
        //Hemispheres
        //Flashgg jets
        int hemisphere_J1() const { if (numberOfFggJets() > 1) { return ( ptOrderedFggJets()[0]->eta() > 0 ? 1 : -1 ); }else{ return 0;}}
        int hemisphere_J2() const { if (numberOfFggJets() > 1) { return ( ptOrderedFggJets()[1]->eta() > 0 ? 1 : -1 ); }else{ return 0;}}
        int hemisphere_J3() const { if (numberOfFggJets() > 2) { return ( ptOrderedFggJets()[2]->eta() > 0 ? 1 : -1 ); }else{ return 0;}}
        //Truth partons
        int hemisphere_P1() const { if (hasLeadingParton()) {return (leadingParton()->eta() > 0 ? 1 : -1 ); }else{ return 0;}} 
        int hemisphere_P2() const { if (hasSubLeadingParton()) {return (subLeadingParton()->eta() > 0 ? 1 : -1 ); }else{ return 0;}} 
        int hemisphere_P3() const { if (hasSubSubLeadingParton()) {return (subSubLeadingParton()->eta() > 0 ? 1 : -1 ); }else{ return 0;}} 

        //Delta Rs
        float dR_J1J2_FggJet() const {if(numberOfFggJets() > 1) {return deltaR(ptOrderedFggJets()[0]->eta(),ptOrderedFggJets()[0]->phi(),
                                                                               ptOrderedFggJets()[1]->eta(),ptOrderedFggJets()[1]->phi()); }else{return -1.;}}
        float dR_J1J3_FggJet() const {if(numberOfFggJets() > 2) {return deltaR(ptOrderedFggJets()[0]->eta(),ptOrderedFggJets()[0]->phi(),
                                                                               ptOrderedFggJets()[2]->eta(),ptOrderedFggJets()[2]->phi()); }else{return -1.;}}
        float dR_J2J3_FggJet() const {if(numberOfFggJets() > 2) {return deltaR(ptOrderedFggJets()[1]->eta(),ptOrderedFggJets()[1]->phi(),
                                                                               ptOrderedFggJets()[2]->eta(),ptOrderedFggJets()[2]->phi()); }else{return -1.;}}
        float dR_J1J2_GenJet() const {if(hasClosestGenJetToLeadingJet() && hasClosestGenJetToSubLeadingJet()) {
                                                                 return deltaR(closestGenJetToLeadingJet()->eta(),closestGenJetToLeadingJet()->phi(),
                                                                               closestGenJetToSubLeadingJet()->eta(),closestGenJetToSubLeadingJet()->phi());
                                     }else{return -1.;}}
        float dR_J1J3_GenJet() const {if(hasClosestGenJetToLeadingJet() && hasClosestGenJetToSubSubLeadingJet()) {
                                                                 return deltaR(closestGenJetToLeadingJet()->eta(),closestGenJetToLeadingJet()->phi(),
                                                                               closestGenJetToSubSubLeadingJet()->eta(),closestGenJetToSubSubLeadingJet()->phi());
                                     }else{return -1.;}}
        float dR_J2J3_GenJet() const {if(hasClosestGenJetToSubLeadingJet() && hasClosestGenJetToSubSubLeadingJet()) {
                                                                 return deltaR(closestGenJetToSubLeadingJet()->eta(),closestGenJetToSubLeadingJet()->phi(),
                                                                               closestGenJetToSubSubLeadingJet()->eta(),closestGenJetToSubSubLeadingJet()->phi());
                                     }else{return -1.;}}
        float dR_J1J2_Particle() const {if(hasClosestParticleToLeadingJet() && hasClosestParticleToSubLeadingJet()) {
                                                                 return deltaR(closestParticleToLeadingJet()->eta(),closestParticleToLeadingJet()->phi(),
                                                                               closestParticleToSubLeadingJet()->eta(),closestParticleToSubLeadingJet()->phi());
                                     }else{return -1.;}}
        float dR_J1J3_Particle() const {if(hasClosestParticleToLeadingJet() && hasClosestParticleToSubSubLeadingJet()) {
                                                                 return deltaR(closestParticleToLeadingJet()->eta(),closestParticleToLeadingJet()->phi(),
                                                                               closestParticleToSubSubLeadingJet()->eta(),closestParticleToSubSubLeadingJet()->phi());
                                     }else{return -1.;}}
        float dR_J2J3_Particle() const {if(hasClosestParticleToSubLeadingJet() && hasClosestParticleToSubSubLeadingJet()) {
                                                                 return deltaR(closestParticleToSubLeadingJet()->eta(),closestParticleToSubLeadingJet()->phi(),
                                                                               closestParticleToSubSubLeadingJet()->eta(),closestParticleToSubSubLeadingJet()->phi());
                                     }else{return -1.;}}
        float dR_P1P2_Partons() const {if(hasLeadingParton() && hasSubLeadingParton()) { return deltaR(leadingParton()->eta(),leadingParton()->phi(),
                                                                                                       subLeadingParton()->eta(),subLeadingParton()->phi());
                                      }else{return -1.;}}
        float dR_P1P3_Partons() const {if(hasLeadingParton() && hasSubSubLeadingParton()) { return deltaR(leadingParton()->eta(),leadingParton()->phi(),
                                                                                                          subSubLeadingParton()->eta(),subSubLeadingParton()->phi());
                                      }else{return -1.;}}
        float dR_P2P3_Partons() const {if(hasSubLeadingParton() && hasSubSubLeadingParton()) { return deltaR(subLeadingParton()->eta(),subLeadingParton()->phi(),
                                                                                                             subSubLeadingParton()->eta(),subSubLeadingParton()->phi());
                                      }else{return -1.;}}
       
        //Invariant Masses 
        //(mjj)
        float mjj_J1J2_FggJet() const {if (numberOfFggJets() > 1) { return (ptOrderedFggJets()[0]->p4() + ptOrderedFggJets()[1]->p4()).mass(); }else{return -1.;}}
        float mjj_J1J3_FggJet() const {if (numberOfFggJets() > 2) { return (ptOrderedFggJets()[0]->p4() + ptOrderedFggJets()[2]->p4()).mass(); }else{return -1.;}}
        float mjj_J2J3_FggJet() const {if (numberOfFggJets() > 2) { return (ptOrderedFggJets()[1]->p4() + ptOrderedFggJets()[2]->p4()).mass(); }else{return -1.;}}
        float mjj_J1J2_GenJet() const {if (hasClosestGenJetToLeadingJet() && hasClosestGenJetToSubLeadingJet()) {
                                            return (closestGenJetToLeadingJet()->p4() + closestGenJetToSubLeadingJet()->p4()).mass(); }else{return -1.;}}
        float mjj_J1J3_GenJet() const {if (hasClosestGenJetToLeadingJet() && hasClosestGenJetToSubSubLeadingJet()) {
                                            return (closestGenJetToLeadingJet()->p4() + closestGenJetToSubSubLeadingJet()->p4()).mass(); }else{return -1.;}}
        float mjj_J2J3_GenJet() const {if (hasClosestGenJetToSubLeadingJet() && hasClosestGenJetToSubSubLeadingJet()) {
                                            return (closestGenJetToSubLeadingJet()->p4() + closestGenJetToSubSubLeadingJet()->p4()).mass(); }else{return -1.;}}
        float mjj_J1J2_GenParticle() const {if (hasClosestParticleToLeadingJet() && hasClosestParticleToSubLeadingJet()) {
                                            return (closestParticleToLeadingJet()->p4() + closestParticleToSubLeadingJet()->p4()).mass(); }else{return -1.;}}
        float mjj_J1J3_GenParticle() const {if (hasClosestParticleToLeadingJet() && hasClosestParticleToSubSubLeadingJet()) {
                                            return (closestParticleToLeadingJet()->p4() + closestParticleToSubSubLeadingJet()->p4()).mass(); }else{return -1.;}}
        float mjj_J2J3_GenParticle() const {if (hasClosestParticleToSubLeadingJet() && hasClosestParticleToSubSubLeadingJet()) {
                                            return (closestParticleToSubLeadingJet()->p4() + closestParticleToSubSubLeadingJet()->p4()).mass(); }else{return -1.;}}
        float mjj_P1P2_Partons() const {if (hasLeadingParton() && hasSubLeadingParton()) {
                                            return (leadingParton()->p4() + subLeadingParton()->p4()).mass();}else{return -1.;}}
        float mjj_P1P3_Partons() const {if (hasLeadingParton() && hasSubSubLeadingParton()) {
                                            return (leadingParton()->p4() + subSubLeadingParton()->p4()).mass();}else{return -1.;}}
        float mjj_P2P3_Partons() const {if (hasSubLeadingParton() && hasSubSubLeadingParton()) {
                                            return (subLeadingParton()->p4() + subSubLeadingParton()->p4()).mass();}else{return -1.;}}
        //(mjjj)
        float mjjj_FggJet() const {if (numberOfFggJets() > 2) {return (ptOrderedFggJets()[0]->p4() + ptOrderedFggJets()[1]->p4() + ptOrderedFggJets()[2]->p4()).mass(); 
                                    }else{return -1.;}}
        float mjjj_GenJet() const {if (hasClosestGenJetToLeadingJet() && hasClosestGenJetToSubLeadingJet() && hasClosestGenJetToSubSubLeadingJet()) {
                                            return (closestGenJetToLeadingJet()->p4() + closestGenJetToSubLeadingJet()->p4() + closestGenJetToSubSubLeadingJet()->p4()).mass(); 
                                    }else{return -1.;}}
        float mjjj_GenParticle() const {if (hasClosestParticleToLeadingJet() && hasClosestParticleToSubLeadingJet() && hasClosestParticleToSubSubLeadingJet()) {
                                            return (closestParticleToLeadingJet()->p4() + closestParticleToSubLeadingJet()->p4() + closestParticleToSubSubLeadingJet()->p4()).mass();
                                    }else{return -1.;}}
        float mjjj_Parton() const {if ( hasLeadingParton() && hasSubLeadingParton() && hasSubSubLeadingParton() ) {
                                            return (leadingParton()->p4() + subLeadingParton()->p4() + subSubLeadingParton()->p4()).mass();
                                    }else{return -1.;}}

        //dEtas
        float dEta_J1J2_FggJet() const {if (numberOfFggJets() > 1) { return fabs(ptOrderedFggJets()[0]->eta()-ptOrderedFggJets()[1]->eta()); }else{return -1.;}}
        float dEta_J1J3_FggJet() const {if (numberOfFggJets() > 2) { return fabs(ptOrderedFggJets()[0]->eta()-ptOrderedFggJets()[2]->eta()); }else{return -1.;}}
        float dEta_J2J3_FggJet() const {if (numberOfFggJets() > 2) { return fabs(ptOrderedFggJets()[1]->eta()-ptOrderedFggJets()[2]->eta()); }else{return -1.;}}
        float dEta_J1J2_GenJet() const {if (hasClosestGenJetToLeadingJet() && hasClosestGenJetToSubLeadingJet()) {
                                            return fabs(closestGenJetToLeadingJet()->eta() - closestGenJetToSubLeadingJet()->eta()); }else{return -1.;}}
        float dEta_J1J3_GenJet() const {if (hasClosestGenJetToLeadingJet() && hasClosestGenJetToSubSubLeadingJet()) {
                                            return fabs(closestGenJetToLeadingJet()->eta() - closestGenJetToSubSubLeadingJet()->eta()); }else{return -1.;}}
        float dEta_J2J3_GenJet() const {if (hasClosestGenJetToSubLeadingJet() && hasClosestGenJetToSubSubLeadingJet()) {
                                            return fabs(closestGenJetToSubLeadingJet()->eta() - closestGenJetToSubSubLeadingJet()->eta()); }else{return -1.;}}
        float dEta_J1J2_Particle() const {if (hasClosestParticleToLeadingJet() && hasClosestParticleToSubLeadingJet()) {
                                            return fabs(closestParticleToLeadingJet()->eta() - closestParticleToSubLeadingJet()->eta()); }else{return -1.;}}
        float dEta_J1J3_Particle() const {if (hasClosestParticleToLeadingJet() && hasClosestParticleToSubSubLeadingJet()) {
                                            return fabs(closestParticleToLeadingJet()->eta() - closestParticleToSubSubLeadingJet()->eta()); }else{return -1.;}}
        float dEta_J2J3_Particle() const {if (hasClosestParticleToSubLeadingJet() && hasClosestParticleToSubSubLeadingJet()) {
                                            return fabs(closestParticleToSubLeadingJet()->eta() - closestParticleToSubSubLeadingJet()->eta()); }else{return -1.;}}
        float dEta_P1P2_Partons() const {if (hasLeadingParton() && hasSubLeadingParton()) {
                                            return fabs(leadingParton()->eta()-subLeadingParton()->eta()); }else{return -1.;}}
        float dEta_P1P3_Partons() const {if (hasLeadingParton() && hasSubSubLeadingParton()) {
                                            return fabs(leadingParton()->eta() - subLeadingParton()->eta()); }else{return -1.;}}
        float dEta_P2P3_Partons() const {if (hasSubLeadingParton() && hasSubSubLeadingParton()) {
                                            return fabs(subLeadingParton()->eta() - subSubLeadingParton()->eta()); }else{return -1.;}}
        //Zeppenfelds
        //(Zepjj)
        float zepjj_J1J2_FggJet() const {if (numberOfFggJets() > 1) {
                                            return fabs(diPhoton()->eta() - 0.5*(ptOrderedFggJets()[0]->eta() + ptOrderedFggJets()[1]->eta())); }else{return -1.;}}
        float zepjj_J1J3_FggJet() const {if (numberOfFggJets() > 2) {
                                            return fabs(diPhoton()->eta() - 0.5*(ptOrderedFggJets()[0]->eta() + ptOrderedFggJets()[2]->eta())); }else{return -1.;}}
        float zepjj_J2J3_FggJet() const {if (numberOfFggJets() > 2) {
                                            return fabs(diPhoton()->eta() - 0.5*(ptOrderedFggJets()[1]->eta() + ptOrderedFggJets()[2]->eta())); }else{return -1.;}}
        float zepjj_J1J2_GenJet() const {if (hasClosestGenJetToLeadingJet() && hasClosestGenJetToSubLeadingJet()) {
                                            return fabs(diPhoton()->eta() - 0.5*(closestGenJetToLeadingJet()->eta() + closestGenJetToSubLeadingJet()->eta())); 
                                         }else{return -1.;}}
        float zepjj_J1J3_GenJet() const {if (hasClosestGenJetToLeadingJet() && hasClosestGenJetToSubSubLeadingJet()) {
                                            return fabs(diPhoton()->eta() - 0.5*(closestGenJetToLeadingJet()->eta() + closestGenJetToSubSubLeadingJet()->eta())); 
                                         }else{return -1.;}}
        float zepjj_J2J3_GenJet() const {if (hasClosestGenJetToSubLeadingJet() && hasClosestGenJetToSubSubLeadingJet()) {
                                            return fabs(diPhoton()->eta() - 0.5*(closestGenJetToSubLeadingJet()->eta() + closestGenJetToSubSubLeadingJet()->eta()));
                                         }else{return -1.;}}
        float zepjj_J1J2_Particle() const {if (hasClosestParticleToLeadingJet() && hasClosestParticleToSubLeadingJet()) {
                                            return fabs(diPhoton()->eta() - 0.5*(closestParticleToLeadingJet()->eta() + closestParticleToSubLeadingJet()->eta()));
                                           }else{return -1.;}}
        float zepjj_J1J3_Particle() const {if (hasClosestParticleToLeadingJet() && hasClosestParticleToSubSubLeadingJet()) {
                                            return fabs(diPhoton()->eta() - 0.5*(closestParticleToLeadingJet()->eta() + closestParticleToSubSubLeadingJet()->eta())); 
                                           }else{return -1.;}}
        float zepjj_J2J3_Particle() const {if (hasClosestParticleToSubLeadingJet() && hasClosestParticleToSubSubLeadingJet()) {
                                            return fabs(diPhoton()->eta() - 0.5*(closestParticleToSubLeadingJet()->eta() + closestParticleToSubSubLeadingJet()->eta())); 
                                           }else{return -1.;}}
        float zepjj_P1P2_Partons() const {if (hasLeadingParton() && hasSubLeadingParton()) {
                                            return fabs(diPhoton()->eta() - 0.5*(leadingParton()->eta() + subLeadingParton()->eta())); }else{return -1.;}}
        float zepjj_P1P3_Partons() const {if (hasLeadingParton() && hasSubSubLeadingParton()) {
                                            return fabs(diPhoton()->eta() - 0.5*(leadingParton()->eta() + subSubLeadingParton()->eta())); }else{return -1.;}}
        float zepjj_P2P3_Partons() const {if (hasSubLeadingParton() && hasSubSubLeadingParton()) {
                                            return fabs(diPhoton()->eta() - 0.5*(subLeadingParton()->eta() + subSubLeadingParton()->eta())); }else{return -1.;}}
        //(Zepjjj)
        float zepjjj_FggJet() const {if (numberOfFggJets() > 2) {
                                            return fabs(diPhoton()->eta() - (ptOrderedFggJets()[0]->eta() + ptOrderedFggJets()[1]->eta() + ptOrderedFggJets()[2]->eta())/3 ); 
                                    }else{return -1.;}}
        float zepjjj_GenJet() const {if (hasClosestGenJetToLeadingJet() && hasClosestGenJetToSubLeadingJet() && hasClosestGenJetToSubSubLeadingJet()) {
                                            return fabs(diPhoton()->eta() - (closestGenJetToLeadingJet()->eta() + closestGenJetToSubLeadingJet()->eta() 
                                                                                + closestGenJetToSubSubLeadingJet()->eta())/3); }else{return -1.;}}
        float zepjjj_Particle() const {if (hasClosestParticleToLeadingJet() && hasClosestParticleToSubLeadingJet() && hasClosestParticleToSubSubLeadingJet()) {
                                            return fabs(diPhoton()->eta() - (closestParticleToLeadingJet()->eta() + closestParticleToSubLeadingJet()->eta() 
                                                                                + closestParticleToSubSubLeadingJet()->eta())/3); }else{return -1.;}}
        float zepjjj_Partons() const {if (hasLeadingParton() && hasSubLeadingParton() && hasSubSubLeadingParton()) {
                                            return fabs(diPhoton()->eta() - (leadingParton()->eta() + subLeadingParton()->eta() + subSubLeadingParton()->eta())); 
                                      }else{return -1.;}}
        //dPhi
        //(jj)
        float dPhijj_J1J2_FggJet() const {if (numberOfFggJets() > 1) {
                                            return deltaPhi(diPhoton()->phi(),(ptOrderedFggJets()[0]->p4()+ptOrderedFggJets()[1]->p4()).phi()); }else{return -1.;}}       
        float dPhijj_J1J3_FggJet() const {if (numberOfFggJets() > 2) {
                                            return deltaPhi(diPhoton()->phi(),(ptOrderedFggJets()[0]->p4()+ptOrderedFggJets()[2]->p4()).phi()); }else{return -1.;}}       
        float dPhijj_J2J3_FggJet() const {if (numberOfFggJets() > 2) {
                                            return deltaPhi(diPhoton()->phi(),(ptOrderedFggJets()[2]->p4()+ptOrderedFggJets()[2]->p4()).phi()); }else{return -1.;}}       
        float dPhijj_J1J2_GenJet() const {if (hasClosestGenJetToLeadingJet() && hasClosestGenJetToSubLeadingJet()) {
                                            return deltaPhi(diPhoton()->phi(),(closestGenJetToLeadingJet()->p4() + closestGenJetToSubLeadingJet()->p4()).phi()); 
                                          }else{return -1.;}}
        float dPhijj_J1J3_GenJet() const {if (hasClosestGenJetToLeadingJet() && hasClosestGenJetToSubSubLeadingJet()) {
                                            return deltaPhi(diPhoton()->phi(),(closestGenJetToLeadingJet()->p4() + closestGenJetToSubSubLeadingJet()->p4()).phi());
                                          }else{return -1.;}}
        float dPhijj_J2J3_GenJet() const {if (hasClosestGenJetToSubLeadingJet() && hasClosestGenJetToSubSubLeadingJet()) {
                                            return deltaPhi(diPhoton()->phi(),(closestGenJetToSubLeadingJet()->p4() + closestGenJetToSubSubLeadingJet()->p4()).phi());
                                         }else{return -1.;}}
        float dPhijj_J1J2_Particle() const {if (hasClosestParticleToLeadingJet() && hasClosestParticleToSubLeadingJet()) {
                                            return deltaPhi(diPhoton()->phi(),(closestParticleToLeadingJet()->p4() + closestParticleToSubLeadingJet()->p4()).phi());
                                         }else{return -1.;}}
        float dPhijj_J1J3_Particle() const {if (hasClosestParticleToLeadingJet() && hasClosestParticleToSubSubLeadingJet()) {
                                            return deltaPhi(diPhoton()->phi(),(closestParticleToLeadingJet()->p4() + closestParticleToSubSubLeadingJet()->p4()).phi());
                                         }else{return -1.;}}
        float dPhijj_J2J3_Particle() const {if (hasClosestParticleToSubLeadingJet() && hasClosestParticleToSubSubLeadingJet()) {
                                            return deltaPhi(diPhoton()->phi(),(closestParticleToSubLeadingJet()->p4() + closestParticleToSubSubLeadingJet()->p4()).phi());
                                         }else{return -1.;}}
        float dPhijj_P1P2_Partons()  const {if (hasLeadingParton() && hasSubLeadingParton()) {
                                            return deltaPhi(diPhoton()->phi(),(leadingParton()->p4() + subLeadingParton()->p4()).phi()); }else{return -1.;}}
        float dPhijj_P1P3_Partons()  const {if (hasLeadingParton() && hasSubSubLeadingParton()) {
                                            return deltaPhi(diPhoton()->phi(),(leadingParton()->p4() + subSubLeadingParton()->p4()).phi()); }else{return -1.;}}
        float dPhijj_P2P3_Partons()  const {if (hasSubLeadingParton() && hasSubSubLeadingParton()) {
                                            return deltaPhi(diPhoton()->phi(),(subLeadingParton()->p4() + subSubLeadingParton()->p4()).phi()); }else{return -1.;}}
        //(jjj)
        float dPhijjj_FggJet() const {if (numberOfFggJets() > 2) {  
                                            return deltaPhi(diPhoton()->phi(),(ptOrderedFggJets()[0]->p4()+ptOrderedFggJets()[1]->p4()+ptOrderedFggJets()[2]->p4()).phi());
                                      }else{return -1.;}}       
        float dPhijjj_GenJet() const {if (hasClosestGenJetToLeadingJet() && hasClosestGenJetToSubLeadingJet() && hasClosestGenJetToSubSubLeadingJet()) {
                                            return deltaPhi(diPhoton()->phi(),(closestGenJetToLeadingJet()->p4() + closestGenJetToSubLeadingJet()->p4()
                                                                                +closestGenJetToSubSubLeadingJet()->p4()).phi()); }else{return -1.;}}
        float dPhijjj_Particle() const {if (hasClosestParticleToLeadingJet() && hasClosestParticleToSubLeadingJet() && hasClosestParticleToSubSubLeadingJet()) {
                                            return deltaPhi(diPhoton()->phi(),(closestParticleToLeadingJet()->p4() + closestParticleToSubLeadingJet()->p4()
                                                                                +closestParticleToSubSubLeadingJet()->p4()).phi()); }else{return -1.;}}
        float dPhijjj_Partons() const {if (hasLeadingParton() && hasSubLeadingParton() && hasSubSubLeadingParton()) {
                                            return deltaPhi(diPhoton()->phi(),(leadingParton()->p4()+subLeadingParton()->p4()+subSubLeadingParton()->p4()).phi());
                                       }else{return -1.;}}

        //Has (thing) methods
        bool hasClosestGenJetToLeadingJet() const { return closestGenJetToLeadingJet_.isNonnull(); }
        bool hasClosestGenJetToSubLeadingJet() const { return closestGenJetToSubLeadingJet_.isNonnull(); }
        bool hasClosestGenJetToSubSubLeadingJet() const { return closestGenJetToSubSubLeadingJet_.isNonnull(); }
        bool hasClosestParticleToLeadingJet() const { return closestParticleToLeadingJet_.isNonnull(); }
        bool hasClosestParticleToSubLeadingJet() const { return closestParticleToSubLeadingJet_.isNonnull(); }
        bool hasClosestParticleToSubSubLeadingJet() const { return closestParticleToSubSubLeadingJet_.isNonnull(); }
        bool hasClosestPartonToLeadingJet() const { return closestPartonToLeadingJet_.isNonnull(); }
        bool hasClosestPartonToSubLeadingJet() const { return closestPartonToSubLeadingJet_.isNonnull(); }
        bool hasClosestPartonToSubSubLeadingJet() const { return closestPartonToSubSubLeadingJet_.isNonnull(); }
        bool hasClosestParticleToLeadingPhoton() const { return closestParticleToLeadingPhoton_.isNonnull(); }
        bool hasClosestParticleToSubLeadingPhoton() const { return closestParticleToSubLeadingPhoton_.isNonnull(); }
        bool hasLeadingParton() const { return leadingParton_.isNonnull(); }
        bool hasSubLeadingParton() const { return subLeadingParton_.isNonnull(); }
        bool hasSubSubLeadingParton() const { return subSubLeadingParton_.isNonnull(); }
        bool hasLeadingGenJet() const { return leadingGenJet_.isNonnull(); }
        bool hasSubLeadingGenJet() const { return subLeadingGenJet_.isNonnull(); }
        bool hasSubSubLeadingGenJet() const { return subSubLeadingGenJet_.isNonnull(); }
        bool hasDijet() const {return numberOfFggJets() > 1;}
        bool hasTrijet() const {return numberOfFggJets() > 2;}

        const edm::Ptr<reco::GenJet> closestGenJetToLeadingJet() const { return closestGenJetToLeadingJet_; }
        const edm::Ptr<reco::GenJet> closestGenJetToSubLeadingJet() const { return closestGenJetToSubLeadingJet_; }
        const edm::Ptr<reco::GenJet> closestGenJetToSubSubLeadingJet() const { return closestGenJetToSubSubLeadingJet_; }
        const edm::Ptr<reco::GenParticle> closestParticleToLeadingJet() const { return closestParticleToLeadingJet_; }
        const edm::Ptr<reco::GenParticle> closestParticleToSubLeadingJet() const { return closestParticleToSubLeadingJet_; }
        const edm::Ptr<reco::GenParticle> closestParticleToSubSubLeadingJet() const { return closestParticleToSubSubLeadingJet_; }
        const edm::Ptr<reco::GenParticle> closestPartonToLeadingJet() const { return closestPartonToLeadingJet_; }
        const edm::Ptr<reco::GenParticle> closestPartonToSubLeadingJet() const { return closestPartonToSubLeadingJet_; }
        const edm::Ptr<reco::GenParticle> closestPartonToSubSubLeadingJet() const { return closestPartonToSubSubLeadingJet_; }
        const edm::Ptr<reco::GenParticle> closestParticleToLeadingPhoton() const { return closestParticleToLeadingPhoton_; }
        const edm::Ptr<reco::GenParticle> closestParticleToSubLeadingPhoton() const { return closestParticleToSubLeadingPhoton_; }

        const edm::Ptr<reco::GenParticle> leadingParton() const { return leadingParton_; }
        const edm::Ptr<reco::GenParticle> subLeadingParton() const { return subLeadingParton_; }
        const edm::Ptr<reco::GenParticle> subSubLeadingParton() const { return subSubLeadingParton_; }
        const edm::Ptr<reco::GenJet> leadingGenJet() const { return leadingGenJet_; }
        const edm::Ptr<reco::GenJet> subLeadingGenJet() const { return subLeadingGenJet_; }
        const edm::Ptr<reco::GenJet> subSubLeadingGenJet() const { return subSubLeadingGenJet_; }
        const edm::Ptr<flashgg::DiPhotonCandidate> diPhoton() const { return diPhoton_; }
        const edm::Ptr<flashgg::Jet> leadingJet() const {return ptOrderedFggJets()[0]; }
        const edm::Ptr<flashgg::Jet> subLeadingJet() const {return ptOrderedFggJets()[1]; }
        const edm::Ptr<flashgg::Jet> subSubLeadingJet() const {return ptOrderedFggJets()[2]; }

        void setClosestGenJetToLeadingJet( const edm::Ptr<reco::GenJet> &val ) { closestGenJetToLeadingJet_ = val; }
        void setClosestGenJetToSubLeadingJet( const edm::Ptr<reco::GenJet> &val ) { closestGenJetToSubLeadingJet_ = val; }
        void setClosestGenJetToSubSubLeadingJet( const edm::Ptr<reco::GenJet> &val ) { closestGenJetToSubSubLeadingJet_ = val; }
        void setClosestParticleToLeadingJet( const edm::Ptr<reco::GenParticle> &val ) { closestParticleToLeadingJet_ = val; }
        void setClosestParticleToSubLeadingJet( const edm::Ptr<reco::GenParticle> &val ) { closestParticleToSubLeadingJet_ = val; }
        void setClosestParticleToSubSubLeadingJet( const edm::Ptr<reco::GenParticle> &val ) { closestParticleToSubSubLeadingJet_ = val; }
        void setClosestPartonToLeadingJet( const edm::Ptr<reco::GenParticle> &val ) { closestPartonToLeadingJet_ = val; }
        void setClosestPartonToSubLeadingJet( const edm::Ptr<reco::GenParticle> &val ) { closestPartonToSubLeadingJet_ = val; }
        void setClosestPartonToSubSubLeadingJet( const edm::Ptr<reco::GenParticle> &val ) { closestPartonToSubSubLeadingJet_ = val; }
        void setClosestParticleToLeadingPhoton( const edm::Ptr<reco::GenParticle> &val ) { closestParticleToLeadingPhoton_ = val; }
        void setClosestParticleToSubLeadingPhoton( const edm::Ptr<reco::GenParticle> &val ) { closestParticleToSubLeadingPhoton_ = val; }
        void setLeadingParton( const edm::Ptr<reco::GenParticle> &val ) { leadingParton_ = val; }
        void setSubLeadingParton( const edm::Ptr<reco::GenParticle> &val ) { subLeadingParton_ = val; }
        void setSubSubLeadingParton( const edm::Ptr<reco::GenParticle> &val ) { subSubLeadingParton_ = val; }
        void setLeadingGenJet( const edm::Ptr<reco::GenJet> &val ) { leadingGenJet_ = val; }
        void setSubLeadingGenJet( const edm::Ptr<reco::GenJet> &val ) { subLeadingGenJet_ = val; }
        void setSubSubLeadingGenJet( const edm::Ptr<reco::GenJet> &val ) { subSubLeadingGenJet_ = val; }

        //Diphoton
        void setDiPhoton( const edm::Ptr<flashgg::DiPhotonCandidate> &val ) {diPhoton_ = val;}

        //Pt ordered collection methods
        void setPtOrderedPartons( const std::vector<edm::Ptr<reco::GenParticle>> &val ) { ptOrderedPartons_ = val; }
        void setPtOrderedGenJets( const std::vector<edm::Ptr<reco::GenJet>> &val ) { ptOrderedGenJets_ = val; }
        void setPtOrderedFggJets( const std::vector<edm::Ptr<flashgg::Jet>> &val ) { ptOrderedFggJets_ = val; }

        const std::vector<edm::Ptr<reco::GenParticle>> ptOrderedPartons() const {return ptOrderedPartons_;}
        const std::vector<edm::Ptr<reco::GenJet>> ptOrderedGenJets() const {return ptOrderedGenJets_;}
        const std::vector<edm::Ptr<flashgg::Jet>> ptOrderedFggJets() const {return ptOrderedFggJets_;}
        
        unsigned int numberOfPartons() const {return ptOrderedPartons_.size();} 
        unsigned int numberOfGenJets() const {return ptOrderedGenJets_.size();} 
        unsigned int numberOfFggJets() const {return ptOrderedFggJets_.size();} 

        VBFTagTruth *clone() const;

    private:
        edm::Ptr<reco::GenJet> closestGenJetToLeadingJet_;
        edm::Ptr<reco::GenJet> closestGenJetToSubLeadingJet_;
        edm::Ptr<reco::GenJet> closestGenJetToSubSubLeadingJet_;
        edm::Ptr<reco::GenParticle> closestParticleToLeadingJet_;
        edm::Ptr<reco::GenParticle> closestParticleToSubLeadingJet_;
        edm::Ptr<reco::GenParticle> closestParticleToSubSubLeadingJet_;
        edm::Ptr<reco::GenParticle> closestPartonToLeadingJet_;
        edm::Ptr<reco::GenParticle> closestPartonToSubLeadingJet_;
        edm::Ptr<reco::GenParticle> closestPartonToSubSubLeadingJet_;
        edm::Ptr<reco::GenParticle> closestParticleToLeadingPhoton_;
        edm::Ptr<reco::GenParticle> closestParticleToSubLeadingPhoton_;

        edm::Ptr<reco::GenParticle> leadingParton_;
        edm::Ptr<reco::GenParticle> subLeadingParton_;
        edm::Ptr<reco::GenParticle> subSubLeadingParton_;

        edm::Ptr<reco::GenJet> leadingGenJet_;
        edm::Ptr<reco::GenJet> subLeadingGenJet_;
        edm::Ptr<reco::GenJet> subSubLeadingGenJet_;

        edm::Ptr<flashgg::DiPhotonCandidate> diPhoton_;

        std::vector<edm::Ptr<reco::GenParticle>> ptOrderedPartons_;
        std::vector<edm::Ptr<reco::GenJet>> ptOrderedGenJets_;
        std::vector<edm::Ptr<flashgg::Jet>> ptOrderedFggJets_;
        
    };
}

#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
