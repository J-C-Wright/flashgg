#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/EXOTag.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

using namespace flashgg;

EXOTag::EXOTag() {}

EXOTag::~EXOTag() {}

EXOTag::EXOTag( edm::Ptr<DiPhotonCandidate> &diphoton, edm::Handle<edm::View<flashgg::Jet>> &jets, 
                edm::Handle<edm::View<flashgg::Electron>> &electrons, double rhoFixedGrid, unsigned eventNumber):
            DiPhotonTagBase::DiPhotonTagBase()
{

    eventNumber_=eventNumber;

    diphoton_=diphoton;
    hasDiphoton_ = !diphoton.isNull();

    rhoFixedGrid_ = rhoFixedGrid;

    if (diphoton_.isNull()){
        hasJets_ = false;
        hasDijet_ = false;
        hasElectrons_ = false;
        hasDielectron_ = false;
    }else{
        
        for (unsigned i=0;i<jets->size();i++){
            edm::Ptr<flashgg::Jet> jet = jets->ptrAt(i);
            float dR_leadDP = deltaR(jet->eta(),jet->phi(),diphoton_->leadingPhoton()->eta(),diphoton_->leadingPhoton()->phi());
            float dR_subLeadDP = deltaR(jet->eta(),jet->phi(),diphoton_->subLeadingPhoton()->eta(),diphoton_->subLeadingPhoton()->phi());
            if (dR_leadDP < 0.5 || dR_subLeadDP < 0.5) {
                continue;
            }else{
                jetsVector_.push_back(jet);
            }
        }

        for (unsigned i=0;i<electrons->size();i++){
            edm::Ptr<flashgg::Electron> electron = electrons->ptrAt(i);
            float dR_leadDP = deltaR(electron->eta(),electron->phi(),diphoton_->leadingPhoton()->eta(),diphoton_->leadingPhoton()->phi());
            float dR_subLeadDP = deltaR(electron->eta(),electron->phi(),diphoton_->subLeadingPhoton()->eta(),diphoton_->subLeadingPhoton()->phi());
            if (dR_leadDP < 0.5 || dR_subLeadDP < 0.5) continue;
            electronsVector_.push_back(electron);
        }

        if (electronsVector_.size() > 0){
            setHasElectrons();
            setDielectron();
        }else{
            hasElectrons_ = false;
            hasDielectron_ = false;
        }
            
        if (jetsVector_.size() > 0){
            setHasJets();
            setDijet();
        }else{
            hasJets_ = false;
            hasDijet_ = false;
        }

    }
}

void EXOTag::setDijet(){

    float leadPt(0), subleadPt(0);
    int leadIndex(-1), subleadIndex(-1);

    for (unsigned i=0;i<jetsVector_.size();i++){
        edm::Ptr<flashgg::Jet> jet = jetsVector_[i];

        if (jet->pt() > leadPt){
            subleadPt = leadPt;
            subleadIndex = leadIndex;
            leadPt = jet->pt();
            leadIndex = i;
        }else if (jet->pt() > subleadPt){
            subleadPt = jet->pt();
            subleadIndex = i;
        }
    }

    if (leadIndex == -1 || subleadIndex == -1){
        hasDijet_=false;
    }else{
        hasDijet_=true;
        dijet_.first = jetsVector_[leadIndex];
        dijet_.second = jetsVector_[subleadIndex];
    }

}

void EXOTag::setDielectron(){

    float leadPt(0), subleadPt(0);
    int leadIndex(-1), subleadIndex(-1);

    for (unsigned i=0;i<electronsVector_.size();i++){
        edm::Ptr<flashgg::Electron> electron = electronsVector_[i];

        if (electron->pt() > leadPt){
            subleadPt = leadPt;
            subleadIndex = leadIndex;
            leadPt = electron->pt();
            leadIndex = i;
        }else if (electron->pt() > subleadPt){
            subleadPt = electron->pt();
            subleadIndex = i;
        }
    }

    if (leadIndex == -1 || subleadIndex == -1){
        hasDielectron_=false;
    }else{
        hasDielectron_=true;
        dielectron_.first = electronsVector_[leadIndex];
        dielectron_.second = electronsVector_[subleadIndex];
    }

}

void EXOTag::setHasJets(){
    if (jetsVector_.size() == 0){
        hasJets_ = false;
    }else{
        hasJets_ = true;
    }
}

void EXOTag::setHasElectrons(){
    if (electronsVector_.size() == 0){
        hasElectrons_ = false;
    }else{
        hasElectrons_ = true;
    }
}

const unsigned EXOTag::getEventNumber() const { return eventNumber_; }

float EXOTag::getDiphotonLeadPt() const { return hasDiphoton_ ? diphoton_->leadingPhoton()->pt() : -999.; }

float EXOTag::getDiphotonSubleadPt() const { return hasDiphoton_ ? diphoton_->subLeadingPhoton()->pt() : -999.; }

float EXOTag::getDiphotonLeadEta() const { return hasDiphoton_ ? diphoton_->leadingPhoton()->eta() : -999.; }

float EXOTag::getDiphotonSubleadEta() const { return hasDiphoton_ ? diphoton_->subLeadingPhoton()->eta() : -999.; }

float EXOTag::getDiphotonMass() const { return hasDiphoton_ ? diphoton_->mass() : -999.; }

float EXOTag::getDiphotonLeadR9() const { return hasDiphoton_ ? diphoton_->leadingPhoton()->r9() : -999.; }

float EXOTag::getDiphotonSubleadR9() const { return hasDiphoton_ ? diphoton_->subLeadingPhoton()->r9() : -999.; }

float EXOTag::getDiphotonLeadEtaSC() const { return hasDiphoton_ ? diphoton_->leadingPhoton()->superCluster()->eta() : -999.; }

float EXOTag::getDiphotonSubleadEtaSC() const { return hasDiphoton_ ? diphoton_->subLeadingPhoton()->superCluster()->eta() : -999.; }

float EXOTag::getDiphotonLeadPhiSC() const { return hasDiphoton_ ? diphoton_->leadingPhoton()->superCluster()->phi() : -999.; }

float EXOTag::getDiphotonSubleadPhiSC() const { return hasDiphoton_ ? diphoton_->subLeadingPhoton()->superCluster()->phi() : -999.; }

int EXOTag::getDiphotonCategory() const {

    if (!hasDiphoton_){
        return -999;
    }else{
        
        float boundaryEB(1.4442);
        float boundaryEELo(1.566), boundaryEEHi(2.5);
        int leadCat(-1), subLeadCat(-1);

        if (fabs(diphoton_->leadingPhoton()->eta()) < boundaryEB){
            leadCat = 0;
        }else if (fabs(diphoton_->leadingPhoton()->eta()) > boundaryEELo
                    && fabs(diphoton_->subLeadingPhoton()->eta()) < boundaryEEHi){
            leadCat = 1;
        }
        if (fabs(diphoton_->subLeadingPhoton()->eta()) < boundaryEB){
            subLeadCat = 0;
        }else if (fabs(diphoton_->subLeadingPhoton()->eta()) > boundaryEELo
                    && fabs(diphoton_->subLeadingPhoton()->eta()) < boundaryEEHi){
            subLeadCat = 1;
        }
        if (leadCat == 0 && subLeadCat == 0){
            return 0;
        }else if ((leadCat == 0 && subLeadCat == 1) || (leadCat == 1 && subLeadCat == 0)){
            return 1;
        }else{
            return -999;
        }
    }
}


int   EXOTag::getDiphotonCutsPass() const {

    if (!hasDiphoton_){
        return -999.;
    }else{

        float boundaryEB(1.4442);
        float boundaryEELo(1.566), boundaryEEHi(2.5);

        if (!passPhotonIDCuts(diphoton_->leadingPhoton())) return 0;
        if (!passPhotonIDCuts(diphoton_->subLeadingPhoton())) return 0;
        
        int cat = getDiphotonCategory();

        //Preselection
        float ptCut(75), mggCutEBEB(230), mggCutEBEE(320);
        float diphotonLeadPt = getDiphotonLeadPt();
        float diphotonSubLeadPt = getDiphotonSubleadPt();
        float diphotonLeadEtaSC = getDiphotonLeadEtaSC();
        float diphotonSubLeadEtaSC = getDiphotonSubleadEtaSC();

        bool passesPtCut = diphotonLeadPt > ptCut && diphotonSubLeadPt > ptCut;
        bool passesLeadEtaSCCut = fabs(diphotonLeadEtaSC) < boundaryEEHi && !(fabs(diphotonLeadEtaSC) > boundaryEB && fabs(diphotonLeadEtaSC) < boundaryEELo);
        bool passesSubLeadEtaSCCut = fabs(diphotonSubLeadEtaSC) < boundaryEEHi && !(fabs(diphotonSubLeadEtaSC) > boundaryEB && fabs(diphotonSubLeadEtaSC) < boundaryEELo);
        bool passesOneBarrelEtaSCCut = fabs(diphotonLeadEtaSC) < boundaryEB || fabs(diphotonSubLeadEtaSC) < boundaryEB;
        bool passesMassCut(false);

        if (cat == 0){
            passesMassCut = diphoton_->mass() > mggCutEBEB;
        }else if (cat == 1){
            passesMassCut = diphoton_->mass() > mggCutEBEE;
        }

        if (passesPtCut && passesLeadEtaSCCut && passesSubLeadEtaSCCut && passesOneBarrelEtaSCCut && passesMassCut){
            return 1;
        }else{
            return 0;
        }
    }
}



bool EXOTag::passPhotonIDCuts(const flashgg::Photon* pho) const {
    float eta = pho->superCluster()->eta();
    int saturated = int(pho->checkStatusFlag(flashgg::Photon::rechitSummaryFlags_t::kSaturated));
    int weird = int(pho->checkStatusFlag(flashgg::Photon::rechitSummaryFlags_t::kWeird));
    float isoCh = pho->egChargedHadronIso();

    float hoe = pho->hadTowOverEm() ;
    float sieie = pho->full5x5_sigmaIetaIeta(); 
    int eleVeto = pho->passElectronVeto();
    bool pass=0;

    float correctedIsoGam = correctIsoGam(pho);
    if (eleVeto){
        if (fabs(eta) < 1.4442 && !saturated){
             if (isoCh <5 && correctedIsoGam < 2.75 && hoe <0.05 && sieie<0.0105){
                 pass=1;
             }
        }
        if (fabs(eta) < 1.4442 && saturated && !weird){
             if (isoCh <5 && correctedIsoGam < 2.75 && hoe <0.05 && sieie<0.0112){
                 pass=1;
             }
        }
        if (fabs(eta) > 1.566 && !saturated){
             if (isoCh <5 && correctedIsoGam < 2.0 && hoe <0.05 && sieie<0.028){
                 pass=1;
             }
        }
        if (fabs(eta) < 1.566 && saturated && !weird){
             if (isoCh <5 && correctedIsoGam < 2.0 && hoe <0.05 && sieie<0.03){
                 pass=1;
             }
        }
    }
    return pass;
}



float EXOTag::correctIsoGam(const flashgg::Photon* pho) const {
    float eta = pho->superCluster()->eta();
    float isoGam =  pho->egPhotonIso();
    float pt = pho->pt();
    float alpha =2.5;
    float A =-1.;
    float kappa =-1.;
    if (fabs(eta) < 0.9 ){ A=0.17 ; kappa =4.5e-3 ;}
    if (fabs(eta) >= 0.9 && fabs(eta)<1.4442 ){ A=0.14 ; kappa =4.5e-3;}
    if (fabs(eta) >= 1.566 && fabs(eta)<2.0 ){ A=0.11 ; kappa =3e-3;}
    if (fabs(eta) >= 2.0 && fabs(eta)<2.2 ){ A=0.14 ; kappa =3e-3;}
    if (fabs(eta) >= 2.2 && fabs(eta)<2.5 ){ A=0.22 ; kappa =3e-3;}
                                                    
    float corrIsoGam = alpha + isoGam - rhoFixedGrid_*A  - kappa *pt;

    return corrIsoGam;
}






float EXOTag::getDiphotonLeadCHI() const { return hasDiphoton_ ? diphoton_->leadingPhoton()->egChargedHadronIso() : -999.; }

float EXOTag::getDiphotonSubleadCHI() const { return hasDiphoton_ ? diphoton_->subLeadingPhoton()->egChargedHadronIso() : -999.; }

float EXOTag::getDiphotonLeadEGPhoIso() const { return hasDiphoton_ ? diphoton_->leadingPhoton()->egPhotonIso() : -999.; }

float EXOTag::getDiphotonSubleadEGPhoIso() const { return hasDiphoton_ ? diphoton_->subLeadingPhoton()->egPhotonIso() : -999.; }

float EXOTag::getDiphotonLeadF5x5SigmaIetaIeta() const { return hasDiphoton_ ? diphoton_->leadingPhoton()->full5x5_sigmaIetaIeta() : -999.; }

float EXOTag::getDiphotonSubleadF5x5SigmaIetaIeta() const { return hasDiphoton_ ? diphoton_->subLeadingPhoton()->full5x5_sigmaIetaIeta() : -999.; }

float EXOTag::getDiphotonLeadFull5x5R9() const { return hasDiphoton_ ? diphoton_->leadingPhoton()->full5x5_r9() : -999.; }

float EXOTag::getDiphotonSubleadFull5x5R9() const { return hasDiphoton_ ? diphoton_->subLeadingPhoton()->full5x5_r9() : -999.; }

float EXOTag::getDiphotonLeadHadronicOverEM() const { return hasDiphoton_ ? diphoton_->leadingPhoton()->hadTowOverEm() : -999.; }

float EXOTag::getDiphotonSubleadHadronicOverEM() const { return hasDiphoton_ ? diphoton_->subLeadingPhoton()->hadTowOverEm() : -999.; }

int EXOTag::getDiphotonLeadIsSaturated() const {
    if (!hasDiphoton_){
        return -999.;
    }else{
        return (int)(diphoton_->leadingPhoton()->checkStatusFlag(flashgg::Photon::rechitSummaryFlags_t::kSaturated) 
                                              && !diphoton_->leadingPhoton()->checkStatusFlag(flashgg::Photon::rechitSummaryFlags_t::kWeird));
    }
}

int EXOTag::getDiphotonSubleadIsSaturated() const {
    if (!hasDiphoton_){
        return -999.;
    }else{
        return (int)(diphoton_->subLeadingPhoton()->checkStatusFlag(flashgg::Photon::rechitSummaryFlags_t::kSaturated) 
                                              && !diphoton_->subLeadingPhoton()->checkStatusFlag(flashgg::Photon::rechitSummaryFlags_t::kWeird));
    }
}

int EXOTag::getDiphotonLeadPassElectronVeto() const { return hasDiphoton_ ? diphoton_->leadingPhoton()->passElectronVeto() : -999.; }

int EXOTag::getDiphotonSubleadPassElectronVeto() const { return hasDiphoton_ ? diphoton_->subLeadingPhoton()->passElectronVeto() : -999.; }

int EXOTag::getDiphotonNConv() const { return hasDiphoton_ ? diphoton_->nConv() : -999; }

float EXOTag::getDiphotonPullConv() const { return hasDiphoton_ ? diphoton_->pullConv() : -999; }

int EXOTag::countJetsOverPT(float ptCut) const {
    unsigned count(0);
    for (unsigned i=0;i<jetsVector_.size();i++){
        if (jetsVector_[i]->pt() > ptCut) count++;
    }
    
    return count;
}

int EXOTag::countElectronsOverPT(float ptCut) const {
    unsigned count(0);
    for (unsigned i=0;i<electronsVector_.size();i++){
        if (electronsVector_[i]->pt()) count++;
    }
    
    return count;
}
        
int EXOTag::getJetMultiplicities_All() const {return hasJets_ ? countJetsOverPT(0.0) : 0;}
int EXOTag::getJetMultiplicities_EGT20() const {return hasJets_ ? countJetsOverPT(20) : 0;}
int EXOTag::getJetMultiplicities_EGT30() const {return hasJets_ ? countJetsOverPT(30) : 0;}
int EXOTag::getJetMultiplicities_EGT40() const {return hasJets_ ? countJetsOverPT(40) : 0;}

float EXOTag::getDijetLeadPt() const { return hasDijet_ ? dijet_.first->pt() : -999.; }
float EXOTag::getDijetSubleadPt() const { return hasDijet_ ? dijet_.second->pt() : -999.; }
float EXOTag::getDijetLeadEta() const { return hasDijet_ ? dijet_.first->eta() : -999.; }
float EXOTag::getDijetSubleadEta() const { return hasDijet_ ? dijet_.second->eta() : -999.; }
float EXOTag::getDijetMass() const { return hasDijet_ ? (dijet_.first->p4()+dijet_.second->p4()).M() : -999.; }
float EXOTag::getDijetDeltaEta() const { return hasDijet_ ? fabs(dijet_.first->eta()-dijet_.second->eta()) : -999.; }
float EXOTag::getDijetZeppenfeld() const { return hasDijet_ ? fabs(diphoton_->eta() - 0.5*(dijet_.first->eta()-dijet_.second->eta())) : -999.; }
float EXOTag::getDijetDeltaPhi_jj() const { return hasDijet_ ? fabs(deltaPhi(dijet_.first->phi() , dijet_.second->phi())) : -999.; } 
float EXOTag::getDijetDeltaPhi_ggjj() const { return hasDijet_ ? fabs(deltaPhi(diphoton_->phi() , (dijet_.first->p4()+dijet_.second->p4()).Phi())) : -999; }

int EXOTag::getElectronMultiplicity_EGT35() const {return hasElectrons_ ? countElectronsOverPT(35) : 0;}
int EXOTag::getElectronMultiplicity_EGT75() const {return hasElectrons_ ? countElectronsOverPT(75) : 0;}

float EXOTag::getDielectronLeadPt() const { return hasDielectron_ ? dielectron_.first->pt() : -999.; }
float EXOTag::getDielectronSubleadPt() const { return hasDielectron_ ? dielectron_.second->pt() : -999.; }
float EXOTag::getDielectronLeadEta() const { return hasDielectron_ ? dielectron_.first->eta() : -999.; }
float EXOTag::getDielectronSubleadEta() const { return hasDielectron_ ? dielectron_.second->eta() : -999.; }
float EXOTag::getDielectronMass() const { return hasDielectron_ ? (dielectron_.first->p4()+dielectron_.second->p4()).M() : -999.; }
float EXOTag::getDielectronDeltaEta() const { return hasDielectron_ ? fabs(dielectron_.first->eta()-dielectron_.second->eta()) : -999.; }
float EXOTag::getDielectronZeppenfeld() const { return hasDielectron_ ? fabs(diphoton_->eta() - 0.5*(dielectron_.first->eta()-dielectron_.second->eta())) : -999.; }
float EXOTag::getDielectronDeltaPhi_ee() const { return hasDielectron_ ? fabs(deltaPhi(dielectron_.first->phi() , dielectron_.second->phi())) : -999.; } 
float EXOTag::getDielectronDeltaPhi_ggee() const { return hasDielectron_ ? fabs(deltaPhi(diphoton_->phi() , (dielectron_.first->p4()+dielectron_.second->p4()).Phi())) : -999; }


// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

