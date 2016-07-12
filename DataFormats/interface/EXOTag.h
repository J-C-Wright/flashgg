#ifndef FLASHgg_EXOTag_h
#define FLASHgg_EXOTag_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"

#include "flashgg/DataFormats/interface/DiPhotonTagBase.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "DataFormats/Math/interface/deltaPhi.h"

namespace flashgg {

    class EXOTag: public DiPhotonTagBase
    {
    public:
        EXOTag();
        ~EXOTag();

        EXOTag *clone() const { return ( new EXOTag( *this ) ); }

        EXOTag( edm::Ptr<DiPhotonCandidate> &diphoton, edm::Handle<edm::View<flashgg::Jet>> &jets, 
                edm::Handle<edm::View<flashgg::Electron>> &electrons,  edm::Handle<edm::View<flashgg::Muon>> &muons, double rhoFixedGrid, unsigned eventNumber);

        const unsigned getEventNumber() const ;
        
        //Photon variables
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
        
        float getDiphotonCosThetaStar() const ;

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

        //Muon variables 
        int getMuonMultiplicity_EGT35() const ;
        int getMuonMultiplicity_EGT75() const ;

        float getDimuonLeadPt() const ;
        float getDimuonSubleadPt() const ;
        float getDimuonLeadEta() const ;
        float getDimuonSubleadEta() const ;
        float getDimuonMass() const ;
        float getDimuonDeltaEta() const ;
        float getDimuonZeppenfeld() const ;
        float getDimuonDeltaPhi_mm() const ;
        float getDimuonDeltaPhi_ggmm() const ;

    private:
        unsigned eventNumber_;

        edm::Ptr<DiPhotonCandidate> diphoton_;

        std::vector<edm::Ptr<flashgg::Jet>> jetsVector_; 
        std::vector<edm::Ptr<flashgg::Electron>> electronsVector_;
        std::vector<edm::Ptr<flashgg::Muon>> muonsVector_;

        typedef std::pair<edm::Ptr<flashgg::Jet>,edm::Ptr<flashgg::Jet>> Dijet;
        typedef std::pair<edm::Ptr<flashgg::Electron>,edm::Ptr<flashgg::Electron>> Dielectron;
        typedef std::pair<edm::Ptr<flashgg::Muon>,edm::Ptr<flashgg::Muon>> Dimuon;

        bool hasDiphoton_;
        bool hasJets_;
        bool hasElectrons_;
        bool hasMuons_;

        double rhoFixedGrid_;

        Dijet dijet_;
        Dielectron dielectron_;
        Dimuon dimuon_; 

        bool hasDijet_;
        bool hasDielectron_;
        bool hasDimuon_; 
        
        void setDijet();
        void setDielectron();
        void setDimuon(); 
    
        void setHasJets();
        void setHasElectrons();
        void setHasMuons(); 
        
        int countJetsOverPT(float ptCut) const ;
        int countElectronsOverPT(float ptCut) const ;
        int countMuonsOverPT(float ptCut) const ; 
            
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

