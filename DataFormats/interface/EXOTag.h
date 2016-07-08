#ifndef FLASHgg_EXOTag_h
#define FLASHgg_EXOTag_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"

#include "flashgg/DataFormats/interface/DiPhotonTagBase.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "DataFormats/Math/interface/deltaPhi.h"

namespace flashgg {

    class EXOTag: public DiPhotonTagBase
    {
    public:
        EXOTag();
        ~EXOTag();

        EXOTag *clone() const { return ( new EXOTag( *this ) ); }

        EXOTag( edm::Ptr<DiPhotonCandidate> &diphoton, edm::Handle<edm::View<flashgg::Jet>> &jets, 
                edm::Handle<edm::View<flashgg::Electron>> &electrons, double rhoFixedGrid, unsigned eventNumber);

        const unsigned getEventNumber() const ;

        int   getDiphotonCutsPass() const ;

        float getDiphotonMass() const ;
        int   getDiphotonNConv() const ;
        float getDiphotonPullConv() const ;

        float getDiphotonLeadPt() const ; 
        float getDiphotonSubleadPt() const ; 
        float getDiphotonLeadEta() const ; 
        float getDiphotonSubleadEta() const ; 

        float getDiphotonLeadR9() const ;
        float getDiphotonSubleadR9() const ;
        float getDiphotonLeadEtaSC() const ;
        float getDiphotonSubleadEtaSC() const ;
        float getDiphotonLeadPhiSC() const ;
        float getDiphotonSubleadPhiSC() const ;
        int   getDiphotonCategory() const ;

        float getDiphotonLeadCHI() const ;
        float getDiphotonSubleadCHI() const ;
        float getDiphotonLeadEGPhoIso() const ;
        float getDiphotonSubleadEGPhoIso() const ;
        float getDiphotonLeadF5x5SigmaIetaIeta() const ;
        float getDiphotonSubleadF5x5SigmaIetaIeta() const ;
        float getDiphotonLeadFull5x5R9() const ;
        float getDiphotonSubleadFull5x5R9() const ;
        float getDiphotonLeadHadronicOverEM() const ;
        float getDiphotonSubleadHadronicOverEM() const ;

        int getDiphotonLeadIsSaturated() const ;
        int getDiphotonSubleadIsSaturated() const ;
        int getDiphotonLeadPassElectronVeto() const ;
        int getDiphotonSubleadPassElectronVeto() const ;

        //Jet variables
        int getJetMultiplicities_All() const ;
        int getJetMultiplicities_EGT20() const ;
        int getJetMultiplicities_EGT30() const ;
        int getJetMultiplicities_EGT40() const ;

        float getDijetLeadPt() const ; 
        float getDijetSubleadPt() const ; 
        float getDijetLeadEta() const ; 
        float getDijetSubleadEta() const ; 
        float getDijetMass() const ;
        float getDijetDeltaEta() const ;
        float getDijetZeppenfeld() const ;
        float getDijetDeltaPhi_jj() const ;
        float getDijetDeltaPhi_ggjj() const ;

        //Electron variables
        int getElectronMultiplicity_EGT35() const ;
        int getElectronMultiplicity_EGT75() const ;

        float getDielectronLeadPt() const ; 
        float getDielectronSubleadPt() const ; 
        float getDielectronLeadEta() const ; 
        float getDielectronSubleadEta() const ; 
        float getDielectronMass() const ;
        float getDielectronDeltaEta() const ;
        float getDielectronZeppenfeld() const ;
        float getDielectronDeltaPhi_ee() const ;
        float getDielectronDeltaPhi_ggee() const ;

    private:
        unsigned eventNumber_;

        edm::Ptr<DiPhotonCandidate> diphoton_;
        std::vector<edm::Ptr<flashgg::Jet>> jetsVector_;
        std::vector<edm::Ptr<flashgg::Electron>> electronsVector_;

        typedef std::pair<edm::Ptr<flashgg::Jet>,edm::Ptr<flashgg::Jet>> Dijet;
        typedef std::pair<edm::Ptr<flashgg::Electron>,edm::Ptr<flashgg::Electron>> Dielectron;

        bool hasDiphoton_;
        bool hasJets_;
        bool hasElectrons_;

        double rhoFixedGrid_;

        Dijet dijet_;
        Dielectron dielectron_;
        bool hasDijet_;
        bool hasDielectron_;

        void setDijet();
        void setDielectron();
        void setHasJets();
        void setHasElectrons();

        int countJetsOverPT(float ptCut) const ;
        int countElectronsOverPT(float ptCut) const ;

        bool passPhotonIDCuts(const flashgg::Photon* pho) const ;
        float correctIsoGam(const flashgg::Photon* pho) const ;
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

