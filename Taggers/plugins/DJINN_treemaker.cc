// DJINNTreeMaker.cc by J. Wright
//
// * Produces dataset trees for DJINN training/evaluation
//
// Adapted from TagTestAnalyzer code by S. Zenz, which was
// Adapted from globelikeTreeMakerWithTagSorter code by L. D. Corpe, which was
// Adapted from the flashggCommissioning tree maker code  by C. Favaro et al.

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include "flashgg/DataFormats/interface/VBFTag.h"
#include "flashgg/DataFormats/interface/UntaggedTag.h"
#include "flashgg/DataFormats/interface/DiPhotonTagBase.h"
#include "flashgg/DataFormats/interface/TTHHadronicTag.h"
#include "flashgg/DataFormats/interface/TTHLeptonicTag.h"
#include "flashgg/DataFormats/interface/VHTightTag.h"
#include "flashgg/DataFormats/interface/VHEtTag.h"
#include "flashgg/DataFormats/interface/VHLooseTag.h"
#include "flashgg/DataFormats/interface/VHHadronicTag.h"
#include "flashgg/DataFormats/interface/VBFTagTruth.h"
#include "flashgg/DataFormats/interface/ZPlusJetTag.h"

#include "TTree.h"

#include "flashgg/DataFormats/interface/Jet.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TMVA/Reader.h"

using namespace std;
using namespace edm;

// ******************************************************************************************
// ******************************************************************************************

//Define the structs for filling the trees
struct JetStruct {
    std::vector<int> pdgIds;
    vector<float> constituents;
};

struct SystWeightStruct {

    float UnmatchedPUWeightUp01sigma;
    float MvaLinearSystUp01sigma;
    float LooseMvaSFUp01sigma;
    float PreselSFUp01sigma;
    float electronVetoSFUp01sigma;
    float TriggerWeightUp01sigma;
    float FracRVWeightUp01sigma;
    float FracRVNvtxWeightUp01sigma;
    float ElectronWeightUp01sigma;
    float MuonWeightUp01sigma;
    float MuonMiniIsoWeightUp01sigma;
    float JetBTagCutWeightUp01sigma;
    float JetBTagReshapeWeightUp01sigma;

    float UnmatchedPUWeightDown01sigma;
    float MvaLinearSystDown01sigma;
    float LooseMvaSFDown01sigma;
    float PreselSFDown01sigma;
    float electronVetoSFDown01sigma;
    float TriggerWeightDown01sigma;
    float FracRVWeightDown01sigma;
    float FracRVNvtxWeightDown01sigma;
    float ElectronWeightDown01sigma;
    float MuonWeightDown01sigma;
    float MuonMiniIsoWeightDown01sigma;
    float JetBTagCutWeightDown01sigma;
    float JetBTagReshapeWeightDown01sigma;

    const TString SystWeightString = TString( "UnmatchedPUWeightUp01sigma/F:"
                                              "MvaLinearSystUp01sigma/F:"
                                              "LooseMvaSFUp01sigma/F:"
                                              "PreselSFUp01sigma/F:"
                                              "electronVetoSFUp01sigma/F:"
                                              "TriggerWeightUp01sigma/F:"
                                              "FracRVWeightUp01sigma/F:"
                                              "FracRVNvtxWeightUp01sigma/F:"
                                              "ElectronWeightUp01sigma/F:"
                                              "MuonWeightUp01sigma/F:"
                                              "MuonMiniIsoWeightUp01sigma/F:"
                                              "JetBTagCutWeightUp01sigma/F:"
                                              "JetBTagReshapeWeightUp01sigma/F:"

                                              "UnmatchedPUWeightDown01sigma/F:"
                                              "MvaLinearSystDown01sigma/F:"
                                              "LooseMvaSFDown01sigma/F:"
                                              "PreselSFDown01sigma/F:"
                                              "electronVetoSFDown01sigma/F:"
                                              "TriggerWeightDown01sigma/F:"
                                              "FracRVWeightDown01sigma/F:"
                                              "FracRVNvtxWeightDown01sigma/F:"
                                              "ElectronWeightDown01sigma/F:"
                                              "MuonWeightDown01sigma/F:"
                                              "MuonMiniIsoWeightDown01sigma/F:"
                                              "JetBTagCutWeightDown01sigma/F:"
                                              "JetBTagReshapeWeightDown01sigma/F" );

};

struct EventStruct {

    float weight;

    float lead_pt_m;
    float sublead_pt_m;
    float total_pt_m;
    float mass_gg;
    float mass_jj;
    float abs_dEta;
    float centrality;
    float dPhi_jj;
    float dPhi_ggjj;
    float dPhi_ggjj_trunc;
    float minDR;

    float leadPhoton_eta;
    float subleadPhoton_eta;
    float cos_dPhi_gg;
    float leadPhotonID;
    float subleadPhotonID;
    float sigma_rv;
    float sigma_wv;
    float prob_vtx;

    float diphoton_score;

    float leadJetPt;
    float leadJetEta;
    float leadJetPhi;

    float subleadJetPt;
    float subleadJetEta;
    float subleadJetPhi;

    float leadJetRMS;
    float subleadJetRMS;

    float leadJetPUJID;
    float subleadJetPUJID;

    float dijet_BDT_score;
    float combined_BDT_score;

    float dZ;
    float HTXSstage0cat;

    const TString eventVariableString = TString("weight/F:lead_pt_m/F:sublead_pt_m/F:total_pt_m/F:mass_gg/F:"
                                                "mass_jj/F:abs_dEta/F:centrality/F:dPhi_jj/F:dPhi_ggjj/F:dPhi_ggjj_trunc/F:minDR/F:"
                                                "leadPhoton_eta/F:subleadPhoton_eta/F:cos_dPhi_gg/F:leadPhotonID/F:subleadPhotonID/F:"
                                                "sigma_rv/F:sigma_wv/F:prob_vtx/F:diphoton_score/F:"
                                                "leadJetPt/F:leadJetEta/F:leadJetPhi/F:"
                                                "subleadJetPt/F:subleadJetEta/F:subleadJetPhi/F:"
                                                "leadJetRMS/F:subleadJetRMS/F:"
                                                "leadJetPUJID/F:subleadJetPUJID/F:"
                                                "dijet_BDT_score/F:combined_BDT_score/F:"
                                                "dZ/F:HTXSstage0cat/F"
                                                );
};


vector<int> partonMatchPdgIds(edm::Ptr<flashgg::Jet> jet, Handle<View<reco::GenParticle>> genParticles){

    vector<int> pdgIds;
    for (unsigned i=0;i<genParticles->size();i++){
        edm::Ptr<reco::GenParticle> particle = genParticles->ptrAt(i);
        if (particle->isHardProcess() &&
            (particle->pdgId() == 21 || abs(particle->pdgId()<7)) &&
            particle->pt() > 0.0){

            float dr = deltaR(jet->eta(),jet->phi(),particle->eta(),particle->phi());
            if (dr <= 0.4){
                pdgIds.push_back(particle->pdgId());
            }

        }
    }

    return pdgIds;
};
// **********************************************************************

namespace flashgg {

    class DJINNTreeMaker : public edm::EDAnalyzer
    {
    public:
        explicit DJINNTreeMaker( const edm::ParameterSet & );
        ~DJINNTreeMaker();

        static void fillDescriptions( edm::ConfigurationDescriptions &descriptions );


    private:

        edm::Service<TFileService> fs_;

        TTree *tree_;

        virtual void beginJob() override;
        virtual void analyze( const edm::Event &, const edm::EventSetup & ) override;
        virtual void endJob() override;

        edm::EDGetTokenT<edm::View<flashgg::DiPhotonTagBase> > TagSorterToken_;

        EDGetTokenT<View<DiPhotonCandidate>> diphotonToken_;
        EDGetTokenT<View<DiPhotonMVAResult>> mvaResultToken_;
        std::vector<edm::InputTag> inputTagJets_;
        std::vector<edm::EDGetTokenT<View<flashgg::Jet>>> tokenJets_;
        EDGetTokenT<GenEventInfoProduct> genInfoToken_;
        edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puInfoToken_;
        EDGetTokenT<View<reco::GenParticle>> genPartToken_;

        double lumiWeight_;
        double xs_;

        bool expectMultiples_;

        typedef std::vector<edm::Handle<edm::View<flashgg::Jet>>> JetCollectionVector;

        //Jet rejection parameters
        double _rmsforwardCut;
        string _JetIDLevel;
        std::vector<double> _pujid_wp_pt_bin_1;
        std::vector<double> _pujid_wp_pt_bin_2;
        std::vector<double> _pujid_wp_pt_bin_3;

        bool _useJetID;

        //Pileup weight parameters
        std::vector<double> _dataPu;
        std::vector<double> _mcPu;
        std::vector<double> _puBins;
        edm::InputTag _puInfo;
        bool _puReWeight;
        edm::InputTag _rho;
        bool _useTruePu;
        edm::InputTag _vertexes;

        bool _isData;
        bool _getPu;

        std::vector<double> _puWeights;

        //Stuff for outuput
        SystWeightStruct systInfo_;
        EventStruct eventInfo_;

        JetStruct leadJetInfo_;
        JetStruct subleadJetInfo_;

        vector<int> leadPdgIds_;
        vector<int> subleadPdgIds_;

        //Old BDTs for comparison studies
        unique_ptr<TMVA::Reader> dijet_BDT_;
        unique_ptr<TMVA::Reader> combined_BDT_;
        FileInPath dijet_BDT_XML_;
        FileInPath combined_BDT_XML_;
        string     BDTMethod_;
    };


// ******************************************************************************************
// ******************************************************************************************

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
    DJINNTreeMaker::DJINNTreeMaker( const edm::ParameterSet &iConfig ):
        TagSorterToken_( consumes<edm::View<flashgg::DiPhotonTagBase> >( iConfig.getParameter<InputTag> ( "TagSorter" ) ) ),
        diphotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
        mvaResultToken_( consumes<View<flashgg::DiPhotonMVAResult> >( iConfig.getParameter<InputTag> ( "MVAResultTag" ) ) ),
        inputTagJets_ ( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
        genInfoToken_( consumes<GenEventInfoProduct>( iConfig.getParameter<InputTag>( "genInfoTag" ) ) ),
        puInfoToken_(consumes<std::vector<PileupSummaryInfo> >( iConfig.getParameter<InputTag> ( "pileupInfo" ) ) ),
        genPartToken_( consumes<View<reco::GenParticle>> ( iConfig.getParameter<InputTag>("GenParticleTag"))),
        lumiWeight_( iConfig.getParameter<double>( "lumiWeight" ) ),
        xs_( iConfig.getParameter<double>( "xs" ) ),
        expectMultiples_( iConfig.getUntrackedParameter<bool>( "ExpectMultiples", false) ),
        _rmsforwardCut( iConfig.getParameter<double> ( "rmsforwardCut") ),
        _JetIDLevel   ( iConfig.getParameter<string> ( "JetIDLevel"   ) ),
        _pujid_wp_pt_bin_1  ( iConfig.getParameter<std::vector<double> > ( "pujidWpPtBin1" ) ),
        _pujid_wp_pt_bin_2  ( iConfig.getParameter<std::vector<double> > ( "pujidWpPtBin2" ) ),
        _pujid_wp_pt_bin_3  ( iConfig.getParameter<std::vector<double> > ( "pujidWpPtBin3" ) ),
        _useJetID( iConfig.getParameter<bool>( "useJetID" ) ),
        _dataPu  ( iConfig.getParameter<std::vector<double> > ( "dataPu" ) ),
        _mcPu    ( iConfig.getParameter<std::vector<double> > ( "mcPu" ) ),
        _puBins  ( iConfig.getParameter<std::vector<double> > ( "puBins" ) ),
        _puInfo  ( iConfig.getParameter<edm::InputTag > ( "puInfo" ) ),
        _puReWeight( iConfig.getParameter<bool>( "puReWeight" ) ),
        _rho     ( iConfig.getParameter<edm::InputTag > ( "rho" ) ),
        _useTruePu( iConfig.getParameter<bool>( "useTruePu" ) ),
        _vertexes( iConfig.getParameter<edm::InputTag > ( "vertexes" ) ),
        _isData( iConfig.getParameter<bool>( "isData" ) ),
        _getPu( iConfig.exists("puInfo") ),
        dijet_BDT_XML_ ( iConfig.getParameter<edm::FileInPath> ( "dijet_BDT_XML" ) ),
        combined_BDT_XML_ ( iConfig.getParameter<edm::FileInPath> ( "combined_BDT_XML" ) ),
        BDTMethod_ ( iConfig.getParameter<string> ( "BDTMethod" ) )
    {
        //Prep jet collection
        for (unsigned i = 0 ; i < inputTagJets_.size() ; i++) {
            auto token = consumes<View<flashgg::Jet> >(inputTagJets_[i]);
            tokenJets_.push_back(token);
        }

        //Prep tree
        tree_ = fs_->make<TTree>("JetData","Jet images and jet variables");

        tree_->Branch("eventVars",&eventInfo_.weight,eventInfo_.eventVariableString);
        tree_->Branch("systWeightVars",&systInfo_.UnmatchedPUWeightUp01sigma,systInfo_.SystWeightString);
        tree_->Branch("leadConstituents",&leadJetInfo_.constituents);
        tree_->Branch("subleadConstituents",&subleadJetInfo_.constituents);
        tree_->Branch("leadPdgIds",&leadPdgIds_);
        tree_->Branch("subleadPdgIds",&subleadPdgIds_);

        //Pileup weights
        _puWeights = _dataPu;

        auto scl  = std::accumulate(_mcPu.begin(),_mcPu.end(),0.) / std::accumulate(_puWeights.begin(),_puWeights.end(),0.); // rescale input distribs to unit ara
        for( size_t ib = 0; ib<_puWeights.size(); ++ib ) { _puWeights[ib] *= scl / _mcPu[ib]; }

        //Legacy BDTs Setup
        //Dijet
        dijet_BDT_.reset(new TMVA::Reader( "!Color:Silent" ) );

        dijet_BDT_->AddVariable( "dijet_LeadJPt"          , &eventInfo_.leadJetPt );
        dijet_BDT_->AddVariable( "dijet_SubJPt"           , &eventInfo_.subleadJetPt );
        dijet_BDT_->AddVariable( "dijet_abs_dEta"         , &eventInfo_.abs_dEta );
        dijet_BDT_->AddVariable( "dijet_Mjj"              , &eventInfo_.mass_jj );
        dijet_BDT_->AddVariable( "dijet_centrality_gg"    , &eventInfo_.centrality );
        dijet_BDT_->AddVariable( "dijet_dipho_dphi_trunc" , &eventInfo_.dPhi_ggjj_trunc );
        dijet_BDT_->AddVariable( "dijet_dphi"             , &eventInfo_.dPhi_jj );
        dijet_BDT_->AddVariable( "dijet_minDRJetPho"      , &eventInfo_.minDR );
        dijet_BDT_->AddVariable( "leadPho_PToM"           , &eventInfo_.lead_pt_m );
        dijet_BDT_->AddVariable( "sublPho_PToM"           , &eventInfo_.sublead_pt_m );

        dijet_BDT_->BookMVA( BDTMethod_.c_str() , dijet_BDT_XML_.fullPath() );

        //Combined
        combined_BDT_.reset(new TMVA::Reader( "!Color:Silent" ) );

        combined_BDT_->AddVariable( "dipho_mva",  &eventInfo_.diphoton_score );
        combined_BDT_->AddVariable( "dijet_mva",  &eventInfo_.dijet_BDT_score );
        combined_BDT_->AddVariable( "dipho_PToM", &eventInfo_.total_pt_m );

        combined_BDT_->BookMVA( "BDT", combined_BDT_XML_.fullPath() );
        
    }

    DJINNTreeMaker::~DJINNTreeMaker()
    {

    }

    void
    DJINNTreeMaker::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup )
    {

        // ********************************************************************************
        // access edm objects

        Handle<edm::View<flashgg::DiPhotonTagBase> > TagSorter;
        iEvent.getByToken( TagSorterToken_, TagSorter );

        Handle<View<flashgg::DiPhotonMVAResult> > mvaResults;
        iEvent.getByToken( mvaResultToken_, mvaResults );

        Handle<View<flashgg::DiPhotonCandidate> > diphotons;
        iEvent.getByToken( diphotonToken_, diphotons );

        Handle<View<reco::GenParticle>> genParticles;
        iEvent.getByToken(genPartToken_,genParticles);

        JetCollectionVector jetCollection( inputTagJets_.size() );
        for( size_t j = 0; j < inputTagJets_.size(); ++j ) {
            iEvent.getByToken( tokenJets_[j], jetCollection[j] );
        }

        edm::Handle<GenEventInfoProduct> genInfo;
        iEvent.getByToken(genInfoToken_, genInfo);

        edm::Handle<std::vector<PileupSummaryInfo> > puInfo;
        iEvent.getByToken(puInfoToken_, puInfo);

        //Scale (XS weight * lumi weight)
        //std::cout << "[DJINN DEBUG] lumiWeight " << lumiWeight_ << std::endl;
        float scale = 1.0;
        if (!_isData){
            scale = xs_*lumiWeight_;
        }
        //std::cout << "Scale " << scale << std::endl;

        //Generator weight
        float genWeight = 1.0;
        if (!_isData){
            genWeight = genInfo->weight();
            //std::cout << "genWeight " << genWeight << std::endl;
        }

        //Pileup weight
        float puWeight = 1.0;
        if (!_isData && _getPu){
            
            double truePu = 0.0, obsPu = 0.0;
            for (auto & frame : *puInfo){
                if (frame.getBunchCrossing() == 0){
                    truePu = frame.getTrueNumInteractions();
                    obsPu = frame.getPU_NumInteractions();
                    break;
                }
            }

            float npu = (_useTruePu ? truePu : obsPu);
            if (_puReWeight){
                if (npu <= _puBins.front() || npu >= _puBins.back()){
                    puWeight = 0.0;
                }else{
                    int ibin = (std::lower_bound(_puBins.begin(),
                                                 _puBins.end(),
                                                 npu) - _puBins.begin()) - 1;
                    puWeight = _puWeights[ibin];
                }
            }
            
        }
        //std::cout << "puWeight " << puWeight << std::endl;

        float event_weight = scale*genWeight*puWeight;

        edm::Ptr<flashgg::DiPhotonCandidate> diphoton;
        edm::Handle<edm::View<flashgg::Jet>> jets;
        flashgg::DiPhotonMVAResult mvares;

        unsigned count = 0;
        for ( auto tag = TagSorter.product()->begin(); tag != TagSorter.product()->end(); tag++) {

            count++;

            //Get systematics weights
            systInfo_.UnmatchedPUWeightUp01sigma = tag->weight("UnmatchedPUWeightUp01sigma");
            systInfo_.MvaLinearSystUp01sigma = tag->weight("MvaLinearSystUp01sigma");
            systInfo_.LooseMvaSFUp01sigma = tag->weight("LooseMvaSFUp01sigma");
            systInfo_.PreselSFUp01sigma = tag->weight("PreselSFUp01sigma");
            systInfo_.electronVetoSFUp01sigma = tag->weight("electronVetoSFUp01sigma");
            systInfo_.TriggerWeightUp01sigma = tag->weight("TriggerWeightUp01sigma");
            systInfo_.FracRVWeightUp01sigma = tag->weight("FracRVWeightUp01sigma");
            systInfo_.FracRVNvtxWeightUp01sigma = tag->weight("FracRVNvtxWeightUp01sigma");
            systInfo_.ElectronWeightUp01sigma = tag->weight("ElectronWeightUp01sigma");
            systInfo_.MuonWeightUp01sigma = tag->weight("MuonWeightUp01sigma");
            systInfo_.MuonMiniIsoWeightUp01sigma = tag->weight("MuonMiniIsoWeightUp01sigma");
            systInfo_.JetBTagCutWeightUp01sigma = tag->weight("JetBTagCutWeightUp01sigma");
            systInfo_.JetBTagReshapeWeightUp01sigma = tag->weight("JetBTagReshapeWeightUp01sigma");

            systInfo_.UnmatchedPUWeightDown01sigma = tag->weight("UnmatchedPUWeightDown01sigma");
            systInfo_.MvaLinearSystDown01sigma = tag->weight("MvaLinearSystDown01sigma");
            systInfo_.LooseMvaSFDown01sigma = tag->weight("LooseMvaSFDown01sigma");
            systInfo_.PreselSFDown01sigma = tag->weight("PreselSFDown01sigma");
            systInfo_.electronVetoSFDown01sigma = tag->weight("electronVetoSFDown01sigma");
            systInfo_.TriggerWeightDown01sigma = tag->weight("TriggerWeightDown01sigma");
            systInfo_.FracRVWeightDown01sigma = tag->weight("FracRVWeightDown01sigma");
            systInfo_.FracRVNvtxWeightDown01sigma = tag->weight("FracRVNvtxWeightDown01sigma");
            systInfo_.ElectronWeightDown01sigma = tag->weight("ElectronWeightDown01sigma");
            systInfo_.MuonWeightDown01sigma = tag->weight("MuonWeightDown01sigma");
            systInfo_.MuonMiniIsoWeightDown01sigma = tag->weight("MuonMiniIsoWeightDown01sigma");
            systInfo_.JetBTagCutWeightDown01sigma = tag->weight("JetBTagCutWeightDown01sigma");
            systInfo_.JetBTagReshapeWeightDown01sigma = tag->weight("JetBTagReshapeWeightDown01sigma");



            // ********************************************************************************
            // get the objects
            diphoton = tag->diPhoton();
            jets = jetCollection[diphoton->jetCollectionIndex()];
            mvares = tag->diPhotonMVA();


            //----Select jets
            std::pair<int,int> dijet_indices(-999,-999);
            std::pair<float,float> dijet_pts(0.0,0.0);

            for (unsigned i=0;i<jets->size();i++){

                //The Jet
                edm::Ptr<flashgg::Jet> jet = jets->ptrAt(i);

                

                //Photon removal
                float dr_leadPhoton = deltaR(jet->eta(),jet->phi(),
                                             diphoton->leadingPhoton()->eta(),
                                             diphoton->leadingPhoton()->phi());
                float dr_subleadPhoton = deltaR(jet->eta(),jet->phi(),
                                                diphoton->subLeadingPhoton()->eta(),
                                                diphoton->subLeadingPhoton()->phi());
                if (dr_leadPhoton <= 0.4 || dr_subleadPhoton <= 0.4){
                    continue;
                }

                //----Jet ID/PUJID/RMS Cuts
                //PUJID
                std::vector<std::pair<double,double>> eta_cuts_(4);
                eta_cuts_[0] = std::make_pair (0    ,2.50 );
                eta_cuts_[1] = std::make_pair (2.50 ,2.75 );
                eta_cuts_[2] = std::make_pair (2.75 ,3.00 );
                eta_cuts_[3] = std::make_pair (3.00 ,10);

                if ( (!_pujid_wp_pt_bin_1.empty())  &&
                     (!_pujid_wp_pt_bin_2.empty())  &&
                     (!_pujid_wp_pt_bin_3.empty())  ){
                    bool pass=false;
                    for (UInt_t eta_bin=0; eta_bin < _pujid_wp_pt_bin_1.size(); eta_bin++ ){
                        if ( fabs( jet->eta() ) >  eta_cuts_[eta_bin].first &&
                             fabs( jet->eta() ) <= eta_cuts_[eta_bin].second){
                            if ( jet->pt() >  20 &&
                                 jet->pt() <= 30 && jet->puJetIdMVA() > _pujid_wp_pt_bin_1[eta_bin] )
                                pass=true;
                            if ( jet->pt() >  30 &&
                                 jet->pt() <= 50 && jet->puJetIdMVA() > _pujid_wp_pt_bin_2[eta_bin] )
                                pass=true;
                            if ( jet->pt() >  50 &&
                                 jet->pt() <= 100&& jet->puJetIdMVA() > _pujid_wp_pt_bin_3[eta_bin] )
                                pass=true;
                            if (jet->pt() > 100) pass = true;
                        }
                    }
                    if (!pass) continue;
                }

                //RMS
                if( fabs( jet->eta() ) > 2.5 && jet->rms() > _rmsforwardCut ){ 
                    continue; 
                }

                //JetID
                if( _useJetID ){
                    if( _JetIDLevel == "Loose" && !jet->passesJetID( flashgg::Loose ) ) continue;
                    if( _JetIDLevel == "Tight" && !jet->passesJetID( flashgg::Tight ) ) continue;
                }
               
                // Abs eta cut at 4.7
                if( fabs( jet->eta() ) > 4.7 ) { continue; }

                //----Selection
                //Conditionals for dijet selection
                if (jets->ptrAt(i)->pt() > dijet_pts.first){
                    dijet_indices.second = dijet_indices.first;
                    dijet_pts.second = dijet_pts.first;
                    dijet_indices.first = i;
                    dijet_pts.first = jets->ptrAt(i)->pt();
                }else if (jets->ptrAt(i)->pt() > dijet_pts.second){
                    dijet_indices.second = i;
                    dijet_pts.second = jets->ptrAt(i)->pt();
                }

            }

            //Apply loose preselection to keep the size down...
            bool loose_ps_pass = false;
            if (dijet_indices.first >= 0 && dijet_indices.second >= 0){

                edm::Ptr<flashgg::Jet> leadJet = jets->ptrAt(dijet_indices.first);
                edm::Ptr<flashgg::Jet> subleadJet = jets->ptrAt(dijet_indices.second);

                bool pass_mjj_cut = (leadJet->p4()+subleadJet->p4()).mass() > 100.0;
                bool pass_leadPt_cut = leadJet->pt() > 30.0;
                bool pass_subleadPt_cut = subleadJet->pt() > 20.0;
                bool pass_leadEta_cut = fabs(leadJet->eta()) < 4.7;
                bool pass_subleadEta_cut = fabs(subleadJet->eta()) < 4.7;
                bool pass_lead_ptom = mvares.leadptom > 1.0/4.0;
                bool pass_sublead_ptom = mvares.subleadptom > 1.0/5.0;
                bool pass_mgg_cut = (diphoton->mass() > 100 && diphoton->mass() < 180);
                bool pass_lead_photonID = mvares.leadmva > -0.2;
                bool pass_sublead_photonID = mvares.subleadmva > -0.2;

                loose_ps_pass = pass_mjj_cut && pass_leadPt_cut && 
                                pass_subleadPt_cut && pass_leadEta_cut &&
                                pass_subleadEta_cut && pass_lead_ptom &&
                                pass_sublead_ptom && pass_mgg_cut &&
                                pass_lead_photonID && pass_sublead_photonID;
            }


            // ********************************************************************************
            // do stuff with selected objects

            //If there's no valid dijet we're not interested...
            if (loose_ps_pass){

                edm::Ptr<flashgg::Jet> leadJet = jets->ptrAt(dijet_indices.first);
                edm::Ptr<flashgg::Jet> subleadJet = jets->ptrAt(dijet_indices.second);

                reco::Candidate::LorentzVector p1 = diphoton->leadingPhoton()->p4();
                reco::Candidate::LorentzVector p2 = diphoton->subLeadingPhoton()->p4();
                reco::Candidate::LorentzVector j1 = leadJet->p4();
                reco::Candidate::LorentzVector j2 = subleadJet->p4();

                eventInfo_.weight = event_weight;

                eventInfo_.lead_pt_m = mvares.leadptom;
                eventInfo_.sublead_pt_m = mvares.subleadptom;
                eventInfo_.total_pt_m = diphoton->pt()/diphoton->mass();
                eventInfo_.mass_gg = diphoton->mass();
                eventInfo_.mass_jj = (j1+j2).mass();
                eventInfo_.abs_dEta = fabs(j1.eta()-j2.eta());
                eventInfo_.centrality = exp((-4.0/pow(j1.eta()-j2.eta(),2))*pow((p1+p2).eta() - 0.5*(j1.eta()+j2.eta()),2));
                eventInfo_.dPhi_jj = fabs(deltaPhi(j1.phi(),j2.phi()));
                eventInfo_.dPhi_ggjj = fabs(deltaPhi( (p1+p2).phi(), (j1+j2).phi() ));
                eventInfo_.dPhi_ggjj_trunc = std::min((float) eventInfo_.dPhi_ggjj, (float) 2.9416);
                
                float dr_p1_j1 = deltaR(p1.eta(),p1.phi(),j1.eta(),j1.phi());
                float dr_p2_j1 = deltaR(p2.eta(),p2.phi(),j1.eta(),j1.phi());
                float dr_p1_j2 = deltaR(p1.eta(),p1.phi(),j2.eta(),j2.phi());
                float dr_p2_j2 = deltaR(p2.eta(),p2.phi(),j2.eta(),j2.phi());
                eventInfo_.minDR = min(min(dr_p1_j1,dr_p2_j1),min(dr_p1_j2,dr_p2_j2));

                eventInfo_.leadPhoton_eta = mvares.leadeta;
                eventInfo_.subleadPhoton_eta = mvares.subleadeta;
                eventInfo_.cos_dPhi_gg = mvares.CosPhi;
                eventInfo_.leadPhotonID = mvares.leadmva;
                eventInfo_.subleadPhotonID = mvares.subleadmva;
                eventInfo_.sigma_rv = mvares.sigmarv;
                eventInfo_.sigma_wv = mvares.sigmawv;
                eventInfo_.prob_vtx = mvares.vtxprob;
                eventInfo_.diphoton_score = mvares.result;

                eventInfo_.leadJetPt = j1.pt();
                eventInfo_.leadJetEta = j1.eta();
                eventInfo_.leadJetPhi = j1.phi();
                eventInfo_.subleadJetPt = j2.pt();
                eventInfo_.subleadJetEta = j2.eta();
                eventInfo_.subleadJetPhi = j2.phi();

                eventInfo_.leadJetRMS = leadJet->rms();
                eventInfo_.subleadJetRMS = subleadJet->rms();
                eventInfo_.leadJetPUJID = leadJet->puJetIdMVA();
                eventInfo_.subleadJetPUJID = subleadJet->puJetIdMVA();

                //Jet structured Info
                leadJetInfo_.constituents = leadJet->getConstituentInfo();
                subleadJetInfo_.constituents = subleadJet->getConstituentInfo();

                leadPdgIds_ = partonMatchPdgIds(leadJet,genParticles); 
                subleadPdgIds_ = partonMatchPdgIds(subleadJet,genParticles); 

                //Include legacy BDT values for comparison studies
                eventInfo_.dijet_BDT_score = dijet_BDT_->EvaluateMVA( BDTMethod_.c_str() );
                eventInfo_.combined_BDT_score = combined_BDT_->EvaluateMVA( "BDT" );

                eventInfo_.dZ = tag->tagTruth()->genPV().z()-tag->diPhoton()->vtx()->z();
                eventInfo_.HTXSstage0cat = tag->tagTruth()->HTXSstage0cat();


            }else{

                //Placeholder jet-based values

                eventInfo_.weight = event_weight;

                eventInfo_.lead_pt_m = mvares.leadptom;
                eventInfo_.sublead_pt_m = mvares.subleadptom;
                eventInfo_.total_pt_m = diphoton->pt()/diphoton->mass();
                eventInfo_.mass_gg = diphoton->mass();

                eventInfo_.mass_jj =  -999.0;
                eventInfo_.abs_dEta =  -999.0;
                eventInfo_.centrality =  -999.0;
                eventInfo_.dPhi_jj =  -999.0;
                eventInfo_.dPhi_ggjj =  -999.0;
                eventInfo_.dPhi_ggjj_trunc =  -999.0;
                eventInfo_.minDR =  -999.0;

                eventInfo_.leadPhoton_eta = mvares.leadeta;
                eventInfo_.subleadPhoton_eta = mvares.subleadeta;
                eventInfo_.cos_dPhi_gg = mvares.CosPhi;
                eventInfo_.leadPhotonID = mvares.leadmva;
                eventInfo_.subleadPhotonID = mvares.subleadmva;
                eventInfo_.sigma_rv = mvares.sigmarv;
                eventInfo_.sigma_wv = mvares.sigmawv;
                eventInfo_.prob_vtx = mvares.vtxprob;
                eventInfo_.diphoton_score = mvares.result;

                eventInfo_.leadJetPt =  -999.0;
                eventInfo_.leadJetEta =  -999.0;
                eventInfo_.leadJetPhi =  -999.0;
                eventInfo_.subleadJetPt =  -999.0;
                eventInfo_.subleadJetEta =  -999.0;
                eventInfo_.subleadJetPhi =  -999.0;

                eventInfo_.leadJetRMS =  -999.0;
                eventInfo_.subleadJetRMS =  -999.0;
                eventInfo_.leadJetPUJID =  -999.0;
                eventInfo_.subleadJetPUJID =  -999.0;

                //Jet structured Info
                vector<float> dummy_f(1,-999.0);
                vector<int> dummy_i(1,-999);

                leadJetInfo_.constituents = dummy_f;
                subleadJetInfo_.constituents = dummy_f;

                leadPdgIds_ = dummy_i;
                subleadPdgIds_ = dummy_i;

                eventInfo_.dijet_BDT_score = -999.0;
                eventInfo_.combined_BDT_score = -999.0;

                eventInfo_.dZ = tag->tagTruth()->genPV().z()-tag->diPhoton()->vtx()->z();
                eventInfo_.HTXSstage0cat = tag->tagTruth()->HTXSstage0cat();

            }

            tree_->Fill();
        }

             
    } // analyse

    void
    DJINNTreeMaker::beginJob()
    {
    }

    void
    DJINNTreeMaker::endJob()
    {
    }

    void
    DJINNTreeMaker::fillDescriptions( edm::ConfigurationDescriptions &descriptions )
    {
        //The following says we do not know what parameters are allowed so do no validation
        // Please change this to state exactly what you do use, even if it is no parameters
        edm::ParameterSetDescription desc;
        desc.setUnknown();
        descriptions.addDefault( desc );
    }

} // namespace flashgg

typedef flashgg::DJINNTreeMaker FlashggDJINNTreeMaker;
DEFINE_FWK_MODULE( FlashggDJINNTreeMaker );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

