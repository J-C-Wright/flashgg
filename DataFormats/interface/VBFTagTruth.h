#ifndef FLASHgg_VBFTagTruth_h
#define FLASHgg_VBFTagTruth_h

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "flashgg/DataFormats/interface/Jet.h"

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
        float pt_genPartMatchingToJ1() const { return ( hasClosestParticleToLeadingJet() ? closestParticleToLeadingJet()->pt() : -1. );}
        float eta_genPartMatchingToJ1() const { return ( hasClosestParticleToLeadingJet() ? closestParticleToLeadingJet()->eta() : -999. );}
        float phi_genPartMatchingToJ1() const { return ( hasClosestParticleToLeadingJet() ? closestParticleToLeadingJet()->phi() : -999. );}
        float pt_genPartMatchingToJ2() const { return ( hasClosestParticleToSubLeadingJet() ? closestParticleToSubLeadingJet()->pt() : -1. );}
        float eta_genPartMatchingToJ2() const { return ( hasClosestParticleToSubLeadingJet() ? closestParticleToSubLeadingJet()->eta() : -999. );}
        float phi_genPartMatchingToJ2() const { return ( hasClosestParticleToSubLeadingJet() ? closestParticleToSubLeadingJet()->phi() : -999. );}
        float pt_genPartMatchingToPho1() const { return ( hasClosestParticleToLeadingPhoton() ? closestParticleToLeadingPhoton()->pt() : -1. );}
        float eta_genPartMatchingToPho1() const { return ( hasClosestParticleToLeadingPhoton() ? closestParticleToLeadingPhoton()->eta() : -999. );}
        float phi_genPartMatchingToPho1() const { return ( hasClosestParticleToLeadingPhoton() ? closestParticleToLeadingPhoton()->phi() : -999. );}
        float pt_genPartMatchingToPho2() const { return ( hasClosestParticleToSubLeadingPhoton() ? closestParticleToSubLeadingPhoton()->pt() : -1. );}
        float eta_genPartMatchingToPho2() const { return ( hasClosestParticleToSubLeadingPhoton() ? closestParticleToSubLeadingPhoton()->eta() : -999. );}
        float phi_genPartMatchingToPho2() const { return ( hasClosestParticleToSubLeadingPhoton() ? closestParticleToSubLeadingPhoton()->phi() : -999. );}
        float pt_Q1() const { return ( hasLeadingParton() ? leadingParton()->pt() : -1. ); }
        float eta_Q1() const { return ( hasLeadingParton() ? leadingParton()->eta() : -999. ); }
        float phi_Q1() const { return ( hasLeadingParton() ? leadingParton()->phi() : -999. ); }
        float pt_Q2() const { return ( hasSubLeadingParton() ? subLeadingParton()->pt() : -1. ); }
        float eta_Q2() const { return ( hasSubLeadingParton() ? subLeadingParton()->eta() : -999. ); }
        float phi_Q2() const { return ( hasSubLeadingParton() ? subLeadingParton()->phi() : -999. ); }

        bool hasClosestGenJetToLeadingJet() const { return closestGenJetToLeadingJet_.isNonnull(); }
        bool hasClosestGenJetToSubLeadingJet() const { return closestGenJetToSubLeadingJet_.isNonnull(); }
        bool hasClosestGenJetToSubSubLeadingJet() const { return closestGenJetToSubSubLeadingJet_.isNonnull(); }
        bool hasClosestParticleToLeadingJet() const { return closestParticleToLeadingJet_.isNonnull(); }
        bool hasClosestParticleToSubLeadingJet() const { return closestParticleToSubLeadingJet_.isNonnull(); }
        bool hasClosestParticleToSubSubLeadingJet() const { return closestParticleToSubSubLeadingJet_.isNonnull(); }
        bool hasClosestParticleToLeadingPhoton() const { return closestParticleToLeadingPhoton_.isNonnull(); }
        bool hasClosestParticleToSubLeadingPhoton() const { return closestParticleToSubLeadingPhoton_.isNonnull(); }
        bool hasLeadingParton() const { return leadingParton_.isNonnull(); }
        bool hasSubLeadingParton() const { return subLeadingParton_.isNonnull(); }
        bool hasSubSubLeadingParton() const { return subSubLeadingParton_.isNonnull(); }
        bool hasLeadingGenJet() const { return leadingGenJet_.isNonnull(); }
        bool hasSubLeadingGenJet() const { return subLeadingGenJet_.isNonnull(); }
        bool hasSubSubLeadingGenJet() const { return subSubLeadingGenJet_.isNonnull(); }

        const edm::Ptr<reco::GenJet> closestGenJetToLeadingJet() const { return closestGenJetToLeadingJet_; }
        const edm::Ptr<reco::GenJet> closestGenJetToSubLeadingJet() const { return closestGenJetToSubLeadingJet_; }
        const edm::Ptr<reco::GenJet> closestGenJetToSubSubLeadingJet() const { return closestGenJetToSubSubLeadingJet_; }
        const edm::Ptr<reco::GenParticle> closestParticleToLeadingJet() const { return closestParticleToLeadingJet_; }
        const edm::Ptr<reco::GenParticle> closestParticleToSubLeadingJet() const { return closestParticleToSubLeadingJet_; }
        const edm::Ptr<reco::GenParticle> closestParticleToSubSubLeadingJet() const { return closestParticleToSubSubLeadingJet_; }
        const edm::Ptr<reco::GenParticle> closestParticleToLeadingPhoton() const { return closestParticleToLeadingPhoton_; }
        const edm::Ptr<reco::GenParticle> closestParticleToSubLeadingPhoton() const { return closestParticleToSubLeadingPhoton_; }

        const edm::Ptr<reco::GenParticle> leadingParton() const { return leadingParton_; }
        const edm::Ptr<reco::GenParticle> subLeadingParton() const { return subLeadingParton_; }
        const edm::Ptr<reco::GenParticle> subSubLeadingParton() const { return subSubLeadingParton_; }
        const edm::Ptr<reco::GenJet> leadingGenJet() const { return leadingGenJet_; }
        const edm::Ptr<reco::GenJet> subLeadingGenJet() const { return subLeadingGenJet_; }
        const edm::Ptr<reco::GenJet> subSubLeadingGenJet() const { return subSubLeadingGenJet_; }
        const edm::Ptr<flashgg::DiPhotonCandidate> diPhoton() const { return diPhoton_; }

        void setClosestGenJetToLeadingJet( const edm::Ptr<reco::GenJet> &val ) { closestGenJetToLeadingJet_ = val; }
        void setClosestGenJetToSubLeadingJet( const edm::Ptr<reco::GenJet> &val ) { closestGenJetToSubLeadingJet_ = val; }
        void setClosestGenJetToSubSubLeadingJet( const edm::Ptr<reco::GenJet> &val ) { closestGenJetToSubSubLeadingJet_ = val; }
        void setClosestParticleToLeadingJet( const edm::Ptr<reco::GenParticle> &val ) { closestParticleToLeadingJet_ = val; }
        void setClosestParticleToSubLeadingJet( const edm::Ptr<reco::GenParticle> &val ) { closestParticleToSubLeadingJet_ = val; }
        void setClosestParticleToSubSubLeadingJet( const edm::Ptr<reco::GenParticle> &val ) { closestParticleToSubSubLeadingJet_ = val; }
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
        
        unsigned int numberOfPartons() {return ptOrderedPartons_.size();} 
        unsigned int numberOfGenJets() {return ptOrderedGenJets_.size();} 
        unsigned int numberOfFggJets() {return ptOrderedFggJets_.size();} 

        VBFTagTruth *clone() const;

    private:
        edm::Ptr<reco::GenJet> closestGenJetToLeadingJet_;
        edm::Ptr<reco::GenJet> closestGenJetToSubLeadingJet_;
        edm::Ptr<reco::GenJet> closestGenJetToSubSubLeadingJet_;
        edm::Ptr<reco::GenParticle> closestParticleToLeadingJet_;
        edm::Ptr<reco::GenParticle> closestParticleToSubLeadingJet_;
        edm::Ptr<reco::GenParticle> closestParticleToSubSubLeadingJet_;
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
