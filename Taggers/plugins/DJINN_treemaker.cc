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

#include "DNN/Tensorflow/interface/Graph.h"
#include "DNN/Tensorflow/interface/Tensor.h"

#include "flashgg/Taggers/interface/JetImage.h"

using namespace std;
using namespace edm;

// ******************************************************************************************
// ******************************************************************************************

//Define the structs for filling the trees
struct JetStruct {
    std::vector<int> pdgIds;
    vector<float> constituents;
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

    const TString eventVariableString = TString("weight/F:lead_pt_m/F:sublead_pt_m/F:total_pt_m/F:mass_gg/F:"
                                                "mass_jj/F:abs_dEta/F:centrality/F:dPhi_jj/F:dPhi_ggjj/F:minDR/F:"
                                                "leadPhoton_eta/F:subleadPhoton_eta/F:cos_dPhi_gg/F:leadPhotonID/F:subleadPhotonID/F:"
                                                "sigma_rv/F:sigma_wv/F:prob_vtx/F:diphoton_score/F:"
                                                "leadJetPt/F:leadJetEta/F:leadJetPhi/F:"
                                                "subleadJetPt/F:subleadJetEta/F:subleadJetPhi/F"
                                                );
};

//Dummy maker functions
JetStruct make_dummy_jet_struct(){

    JetStruct jetInfo;
    std::vector<float> dummy_vals(1,-999.);
    jetInfo.constituents = dummy_vals;

    return jetInfo;
};

EventStruct make_dummy_event_struct(){

    EventStruct eventInfo;

    eventInfo.weight = -999.;

    eventInfo.lead_pt_m = -999.;
    eventInfo.sublead_pt_m = -999.;
    eventInfo.total_pt_m = -999.;
    eventInfo.mass_gg = -999.;
    eventInfo.mass_jj = -999.;
    eventInfo.abs_dEta = -999.;
    eventInfo.centrality = -999.;
    eventInfo.dPhi_jj = -999.;
    eventInfo.dPhi_ggjj = -999.;
    eventInfo.minDR = -999.;

    eventInfo.leadPhoton_eta = -999.;
    eventInfo.subleadPhoton_eta = -999.;
    eventInfo.cos_dPhi_gg = -999.;
    eventInfo.leadPhotonID = -999.;
    eventInfo.subleadPhotonID = -999.;
    eventInfo.sigma_rv = -999.;
    eventInfo.sigma_wv = -999.;
    eventInfo.prob_vtx = -999.;

    eventInfo.diphoton_score = -999.;


    return eventInfo;
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

        EventStruct eventInfo_;
        JetStruct leadJetInfo_;
        JetStruct subleadJetInfo_;

        vector<int> leadPdgIds_;
        vector<int> subleadPdgIds_;

        dnn::tf::Graph* g_;
        dnn::tf::Tensor* x_im_;
        dnn::tf::Tensor* x_ef_;

        dnn::tf::Tensor* kp_conv_;
        dnn::tf::Tensor* kp_hidd_;

        dnn::tf::Tensor* y_;

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
        expectMultiples_( iConfig.getUntrackedParameter<bool>( "ExpectMultiples", false) )
    {
        for (unsigned i = 0 ; i < inputTagJets_.size() ; i++) {
            auto token = consumes<View<flashgg::Jet> >(inputTagJets_[i]);
            tokenJets_.push_back(token);
        }

        tree_ = fs_->make<TTree>("JetData","Jet images and jet variables");

        tree_->Branch("eventVars",&eventInfo_.weight,eventInfo_.eventVariableString);
        tree_->Branch("leadConstituents",&leadJetInfo_.constituents);
        tree_->Branch("subleadConstituents",&subleadJetInfo_.constituents);
        tree_->Branch("leadPdgIds",&leadPdgIds_);
        tree_->Branch("subleadPdgIds",&subleadPdgIds_);
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

        float genWeight;
        genWeight = genInfo->weight();

        float puWeight = 1.0;

        edm::Ptr<flashgg::DiPhotonCandidate> diphoton;
        edm::Handle<edm::View<flashgg::Jet>> jets;
        flashgg::DiPhotonMVAResult mvares;

        unsigned count = 0;
        for ( auto tag = TagSorter.product()->begin(); tag != TagSorter.product()->end(); tag++) {
            count++;
            // ********************************************************************************
            // get the objects
            diphoton = tag->diPhoton();
            jets = jetCollection[diphoton->jetCollectionIndex()];
            mvares = tag->diPhotonMVA();

            float event_weight = genWeight*lumiWeight_*xs_*puWeight;

            std::pair<int,int> dijet_indices(-999,-999);
            std::pair<float,float> dijet_pts(0.0,0.0);

            for (unsigned i=0;i<jets->size();i++){

                edm::Ptr<flashgg::Jet> jet = jets->ptrAt(i);
                float dr_leadPhoton = deltaR(jet->eta(),jet->phi(),
                                             diphoton->leadingPhoton()->eta(),
                                             diphoton->leadingPhoton()->phi());
                float dr_subleadPhoton = deltaR(jet->eta(),jet->phi(),
                                                diphoton->subLeadingPhoton()->eta(),
                                                diphoton->subLeadingPhoton()->phi());
                
                if (dr_leadPhoton <= 0.4 || dr_subleadPhoton <= 0.4){
                    continue;
                }

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

            // ********************************************************************************
            // do stuff with objects

            //If there's no valid dijet we're not interested...
            if (dijet_indices.first >= 0 && dijet_indices.second >= 0){

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

                //Jet Struct Info
                leadJetInfo_.constituents = leadJet->getConstituentInfo();

                subleadJetInfo_.constituents = subleadJet->getConstituentInfo();


                leadPdgIds_ = partonMatchPdgIds(leadJet,genParticles); 
                subleadPdgIds_ = partonMatchPdgIds(subleadJet,genParticles); 

            
                //Tensorflow DJINN test stuff...

                std::cout << "Making images..." << std::endl;

                std::vector<std::vector<float>> lead_reshaped;
                constituent_vector_from_flat(leadJetInfo_.constituents,lead_reshaped);
                std::vector<std::vector<std::vector<float>>> lead_image = blank_image(24,3);
                image_from_vector(lead_image,lead_reshaped,24,0.5,0);
                print_image(lead_image);     
                std::cout << "Lead image done" << std::endl;

                std::vector<std::vector<float>> sublead_reshaped;
                constituent_vector_from_flat(subleadJetInfo_.constituents,sublead_reshaped);
                std::vector<std::vector<std::vector<float>>> sublead_image = blank_image(24,3);
                image_from_vector(sublead_image,sublead_reshaped,24,0.5,0);
                print_image(sublead_image);

                std::cout << "Sublead image done" << std::endl;

                //Switch dropout off...
                kp_conv_->setValue<float>(0,1.0);
                kp_hidd_->setValue<float>(0,1.0);

                for (unsigned i=0;i<24;i++){
                    for (unsigned j=0;j<24;j++){
                        for (unsigned k=0;k<3;k++){
                            x_im_->setValue<float>(0,i,j,k,lead_image[i][j][k]);
                            x_im_->setValue<float>(0,i,j,k+3,sublead_image[i][j][k]);
                        }
                    }
                }

//                ['eventInfo_.lead_pt_m', 'eventInfo_.sublead_pt_m', 'eventInfo_.total_pt_m', 
//                 'eventInfo_.mass_jj', 'eventInfo_.abs_dEta', 'eventInfo_.centrality', 'eventInfo_.dPhi_jj', 
//                 'eventInfo_.dPhi_ggjj', 'eventInfo_.minDR', 'leadJetVars_pt', 'subleadJetVars_pt', 
//                 'leadJetVars_eta', 'subleadJetVars_eta']

                x_ef_->setValue<float>(0,0,eventInfo_.lead_pt_m);
                x_ef_->setValue<float>(0,1,eventInfo_.sublead_pt_m);
                x_ef_->setValue<float>(0,2,eventInfo_.total_pt_m);
                x_ef_->setValue<float>(0,3,eventInfo_.mass_jj);
                x_ef_->setValue<float>(0,4,eventInfo_.abs_dEta);
                x_ef_->setValue<float>(0,5,eventInfo_.centrality);
                x_ef_->setValue<float>(0,6,eventInfo_.dPhi_jj);
                x_ef_->setValue<float>(0,7,eventInfo_.dPhi_ggjj);
                x_ef_->setValue<float>(0,8,eventInfo_.minDR);
                x_ef_->setValue<float>(0,9,leadJet->pt());
                x_ef_->setValue<float>(0,10,subleadJet->pt());
                x_ef_->setValue<float>(0,11,leadJet->pt());
                x_ef_->setValue<float>(0,12,subleadJet->pt());

//                float in_value = 0.5;
//                for (unsigned i=0;i<13;i++){
//                    x_ef_->setValue<float>(0,i,in_value);
//                }

                std::cout << "Evaluating graph..." << std::endl;
                g_->eval();

                float out_value_1 = y_->getValue<float>(0,0);
                float out_value_2 = y_->getValue<float>(0,1);

                std::cout << "The output values: " << std::endl;
                std::cout << out_value_1 << " " << out_value_2 << std::endl;

                tree_->Fill();

            }        
        }
             
    } // analyze

    void
    DJINNTreeMaker::beginJob()
    {

        g_ = new dnn::tf::Graph("/home/hep/jw3914/Work/flashgg_tensorflow/CMSSW_8_0_28/src/flashgg/models/Model_dummy");

        dnn::tf::Shape x_im_Shape[] = {1,24,24,6};
        dnn::tf::Shape x_ef_Shape[] = {1,13};
        dnn::tf::Shape kp_Shape[] = {1};

        x_im_ = g_->defineInput(new dnn::tf::Tensor("im_in:0", 4, x_im_Shape));
        x_ef_ = g_->defineInput(new dnn::tf::Tensor("ef_in:0", 2, x_ef_Shape));

        kp_conv_ = g_->defineInput(new dnn::tf::Tensor("kp_conv:0", 1, kp_Shape));
        kp_hidd_ = g_->defineInput(new dnn::tf::Tensor("kp_hidd:0", 1, kp_Shape));

        y_ = g_->defineOutput(new dnn::tf::Tensor("y_prob:0"));

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

