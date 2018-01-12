// TagTestAnalyzer.cc by S. Zenz
//
// * Tests getting tags out of the event, sorting them, casting them
// * Dumps debugging output to the screen
// * Useful for quick tests of code changes, and should be kept up-to-date as tags are added/changed
// * Should NOT be included in productions
//
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

#include "DNN/Tensorflow/interface/Graph.h"
#include "DNN/Tensorflow/interface/Tensor.h"

using namespace std;
using namespace edm;

// **********************************************************************

namespace flashgg {

    class TagTestAnalyzer : public edm::EDAnalyzer
    {
    public:
        explicit TagTestAnalyzer( const edm::ParameterSet & );
        ~TagTestAnalyzer();

        static void fillDescriptions( edm::ConfigurationDescriptions &descriptions );


    private:

        edm::Service<TFileService> fs_;

        virtual void beginJob() override;
        virtual void analyze( const edm::Event &, const edm::EventSetup & ) override;
        virtual void endJob() override;

        edm::EDGetTokenT<edm::View<flashgg::DiPhotonTagBase> > TagSorterToken_;
        bool expectMultiples_;

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
    TagTestAnalyzer::TagTestAnalyzer( const edm::ParameterSet &iConfig ):
        TagSorterToken_( consumes<edm::View<flashgg::DiPhotonTagBase> >( iConfig.getParameter<InputTag> ( "TagSorter" ) ) ),
        expectMultiples_( iConfig.getUntrackedParameter<bool>( "ExpectMultiples", false) )
    {


    }

    TagTestAnalyzer::~TagTestAnalyzer()
    {

    }

    void
    TagTestAnalyzer::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup )
    {

        //Switch dropout off...
        kp_conv_->setValue<float>(0,1.0);
        kp_hidd_->setValue<float>(0,1.0);


        float in_value = 0.5;

        for (unsigned i=0;i<24;i++){
            for (unsigned j=0;j<24;j++){
                for (unsigned k=0;k<6;k++){

                    x_im_->setValue<float>(0,i,j,k,in_value);
                }
            }
        }

        for (unsigned i=0;i<13;i++){
            x_ef_->setValue<float>(0,i,in_value);
        }

        std::cout << "Evaluating graph..." << std::endl;
        g_->eval();

        float out_value_1 = y_->getValue<float>(0,0);
        float out_value_2 = y_->getValue<float>(0,1);

        std::cout << "The output values: " << std::endl;
        std::cout << out_value_1 << " " << out_value_2 << std::endl;


    } 

    void
    TagTestAnalyzer::beginJob()
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
    TagTestAnalyzer::endJob()
    {
    }

    void
    TagTestAnalyzer::fillDescriptions( edm::ConfigurationDescriptions &descriptions )
    {
        //The following says we do not know what parameters are allowed so do no validation
        // Please change this to state exactly what you do use, even if it is no parameters
        edm::ParameterSetDescription desc;
        desc.setUnknown();
        descriptions.addDefault( desc );
    }

} // namespace flashgg

typedef flashgg::TagTestAnalyzer FlashggTagTestAnalyzer;
DEFINE_FWK_MODULE( FlashggTagTestAnalyzer );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

