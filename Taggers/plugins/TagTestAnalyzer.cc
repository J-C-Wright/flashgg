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
#include "flashgg/DataFormats/interface/DiPhotonUntaggedCategory.h"
#include "flashgg/DataFormats/interface/DiPhotonTagBase.h"
#include "flashgg/DataFormats/interface/TTHHadronicTag.h"
#include "flashgg/DataFormats/interface/TTHLeptonicTag.h"
#include "flashgg/DataFormats/interface/VHTightTag.h"
#include "flashgg/DataFormats/interface/VHLooseTag.h"
#include "flashgg/DataFormats/interface/VHHadronicTag.h"

#include "TVector.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH2F.h"

class Particle {
  public:
  float pt;
  float eta;
  float phi;
  int   pdgid;
  int   status;
  int   number;
  float rapidity;  
  bool  match;
  int   quark;

  Particle(void){

	pt=0;
	eta=0;
	phi=0;
	pdgid=0;
	status=0;
	number=0;
	rapidity=0;
	match = true;
	
  }

  float deltaR(Particle quark) {

	float error(0);
	error += pow((eta-quark.eta),2);
	error += pow((phi-quark.phi),2);

	return sqrt(error);
  
  }

  void print(void) {

    std::cout << std::setw(12) << pt;
    std::cout << std::setw(12) << eta;
    std::cout << std::setw(12) << phi;
    std::cout << std::setw(12) << status;
    std::cout << std::setw(12) << pdgid;
    std::cout << std::setw(12) << rapidity;
    std::cout << std::endl;

  }

};

using namespace std;
using namespace edm;

// **********************************************************************

namespace flashgg {

	class TagTestAnalyzer : public edm::EDAnalyzer {
		public:
			explicit TagTestAnalyzer(const edm::ParameterSet&);
			~TagTestAnalyzer();

			static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

		private:

	  		edm::Service<TFileService> fs_;

	  		virtual void beginJob() override;
	  		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
			virtual void endJob() override;

			edm::EDGetTokenT<edm::OwnVector<flashgg::DiPhotonTagBase> > TagSorterToken_;

			int eventCount;
			int misIdJets_1;
			int misIdJets_2;

			edm::InputTag src_;
			TFile * outputFile_;
			TH1F * h_m_untagged_0;
			TH1F * h_m_untagged_1;
			TH1F * h_m_untagged_2;
			TH1F * h_m_untagged_3;
			TH1F * h_m_untagged_4;
			TH1F * h_m_vbf_0;
			TH1F * h_m_vbf_1;
			TH1F * h_m_vbf_2;
			TH1F * h_m_tthhad;
			TH1F * h_m_tthlep;
			TH1F * h_m_vhtight;
			TH1F * h_m_vhloose;
			TH1F * h_m_vhhad;
			
			TH1F * deltaR_VBF;
			TH1F * misIdJet_pt;
			TH1F * misIdJet_eta;
			TH1F * misIdJet_phi;
			TH2F * etaVsEta;	
			TH1F * trueJetEta;


			EDGetTokenT< edm::View<reco::GenParticle> > genPartToken_;
			Int_t eventNumber;
			
			TTree * eventTree;
			Particle particle;

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
	TagTestAnalyzer::TagTestAnalyzer(const edm::ParameterSet& iConfig):
		TagSorterToken_(consumes<edm::OwnVector<flashgg::DiPhotonTagBase> >(iConfig.getUntrackedParameter<InputTag> ("TagSorter", InputTag("flashggTagSorter")))),
		genPartToken_ (consumes<View<reco::GenParticle> >(iConfig.getUntrackedParameter<InputTag> ("GenParticleTag", InputTag("prunedGenParticles"))))
	{
		eventNumber = 0;
		misIdJets_1 = 0;
		misIdJets_2 = 0;
	}

	TagTestAnalyzer::~TagTestAnalyzer()
	{
		eventNumber =0;
	}

	void
	TagTestAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
				
		eventCount++;
		vector<int> jetQuarks;

		Handle<View<reco::GenParticle> > genParticles;
  		iEvent.getByToken(genPartToken_,genParticles);
		std::cout << setw(12) << "Pt" << setw(12) << "Eta" << setw(12) << "Phi"; 
		std::cout << setw(12) << "Status" << setw(12) << "PDG id" << setw(12) << "Rapidity\n";

  		for( unsigned int genLoop =0 ; genLoop < genParticles->size(); genLoop++){

    		  particle.pt = genParticles->ptrAt(genLoop)->pt();
    		  particle.eta = genParticles->ptrAt(genLoop)->eta();
    		  particle.phi = genParticles->ptrAt(genLoop)->phi();
    		  particle.status = genParticles->ptrAt(genLoop)->status();
    		  particle.pdgid = genParticles->ptrAt(genLoop)->pdgId();
		  particle.number = eventCount;
		  particle.number = genParticles->ptrAt(genLoop)->rapidity();
 
		  //Get the status = 3 quarks

		  int partStatus = genParticles->ptrAt(genLoop)->status();
		  int partID = genParticles->ptrAt(genLoop)->pdgId();

		  if (partStatus == 3) {
			if (abs(partID) <= 6 && abs(partID) >= 1) {
				jetQuarks.push_back(genLoop);				
		  		std::cout << setw(12) << genParticles->ptrAt(genLoop)->pt();
	    			std::cout << setw(12) << genParticles->ptrAt(genLoop)->eta();
	    			std::cout << setw(12) << genParticles->ptrAt(genLoop)->phi();
    			  	std::cout << setw(12) << genParticles->ptrAt(genLoop)->status();
    			  	std::cout << setw(12) << genParticles->ptrAt(genLoop)->pdgId();
			 	std::cout << setw(12) << genParticles->ptrAt(genLoop)->rapidity();
			  	std::cout << std::endl;
			}
		  }

		  //Fill TTree
		  eventTree->Fill();

		}
		  //Take the last two and create particle objects
		  Particle quark1,quark2;
		if (jetQuarks.size() > 1) {

		  quark1.pt = genParticles->ptrAt(jetQuarks.at(jetQuarks.size()-1))->pt();
		  quark1.eta = genParticles->ptrAt(jetQuarks.at(jetQuarks.size()-1))->eta();
		  quark1.phi = genParticles->ptrAt(jetQuarks.at(jetQuarks.size()-1))->phi();
		  quark1.rapidity = genParticles->ptrAt(jetQuarks.at(jetQuarks.size()-1))->rapidity();

		  quark2.pt = genParticles->ptrAt(jetQuarks.at(jetQuarks.size()-2))->pt();
		  quark2.eta = genParticles->ptrAt(jetQuarks.at(jetQuarks.size()-2))->eta();
		  quark2.phi = genParticles->ptrAt(jetQuarks.at(jetQuarks.size()-2))->phi();
		  quark2.rapidity = genParticles->ptrAt(jetQuarks.at(jetQuarks.size()-2))->rapidity();
		

		}
		trueJetEta->Fill(quark1.eta);
 		trueJetEta->Fill(quark2.eta);

		Handle<edm::OwnVector<flashgg::DiPhotonTagBase> > TagSorter;
		iEvent.getByToken(TagSorterToken_,TagSorter);

		if (TagSorter.product()->size() > 0 ) 
		{
			const flashgg::DiPhotonTagBase* chosenTag = &*(TagSorter.product()->begin());

			const	DiPhotonUntaggedCategory *untagged = dynamic_cast<const DiPhotonUntaggedCategory*>(chosenTag);
			if(untagged != NULL) {
				int category = untagged->categoryNumber();
				
				if (category == 0) {

					std::cout << "Tagged as untagged - category " << category << std::endl;
					std::cout << "Mass  " << untagged->diPhoton()->mass() << " GeV" << std::endl;
					h_m_untagged_0->Fill(untagged->diPhoton()->mass()); 

				}else if (category == 1) {

                                        std::cout << "Tagged as untagged - category " << category << std::endl;
                                        std::cout << "Mass  " << untagged->diPhoton()->mass() << " GeV" << std::endl;
                                        h_m_untagged_1->Fill(untagged->diPhoton()->mass());

                                }else if (category == 2) {

                                        std::cout << "Tagged as untagged - category " << category << std::endl;
                                        std::cout << "Mass  " << untagged->diPhoton()->mass() << " GeV" << std::endl;
                                        h_m_untagged_2->Fill(untagged->diPhoton()->mass());

                                }else if (category == 3) {

                                        std::cout << "Tagged as untagged - category " << category << std::endl;
                                        std::cout << "Mass  " << untagged->diPhoton()->mass() << " GeV" << std::endl;
                                        h_m_untagged_3->Fill(untagged->diPhoton()->mass());

                                }else if (category == 4) {

                                        std::cout << "Tagged as untagged - category " << category << std::endl;
                                        std::cout << "Mass  " << untagged->diPhoton()->mass() << " GeV" << std::endl;
                                        h_m_untagged_4->Fill(untagged->diPhoton()->mass());

                                }

			}

			const	VBFTag *vbftag = dynamic_cast<const VBFTag*>(chosenTag);
			if(vbftag != NULL) {
				
				Particle leadJet, subLeadJet;

				int category = vbftag->categoryNumber();
				std::cout << "Tagged as VBF - category " << category << std::endl;
				std::cout << setw(12) << "Pt" << setw(12) << "Eta" << setw(12) << "Phi" << setw(12) << "Rapidity" << std::endl;

				if (category == 0) {

					h_m_vbf_0->Fill(vbftag->diPhoton()->mass()); 
									
				}else if (category == 1) {

                                        h_m_vbf_1->Fill(vbftag->diPhoton()->mass());
						
                                }else if (category == 2) {

                                        h_m_vbf_2->Fill(vbftag->diPhoton()->mass());
						
				}

				leadJet.pt = vbftag->leadingJet().pt();
				leadJet.eta = vbftag->leadingJet().eta();
				leadJet.phi = vbftag->leadingJet().phi();
				leadJet.rapidity = vbftag->leadingJet().rapidity();
				subLeadJet.pt = vbftag->subLeadingJet().pt();
				subLeadJet.eta = vbftag->subLeadingJet().eta();
				subLeadJet.phi = vbftag->subLeadingJet().phi();
				subLeadJet.rapidity = vbftag->subLeadingJet().rapidity();

				std::cout << "Jets: " << std::endl;
				leadJet.print();
				subLeadJet.print();
				std::cout << "Quarks: " << std::endl;
				quark1.print();
				quark2.print();
				std::cout << "delta R: " << std::endl;
				std::cout << setw(12) << leadJet.deltaR(quark1) << setw(12) << leadJet.deltaR(quark2) << std::endl; 
				std::cout << setw(12) << subLeadJet.deltaR(quark1) << setw(12) << subLeadJet.deltaR(quark2) << std::endl; 

			//Jet-quark allocation
				int numMissIdJets(0);
			//Leading jet
				//Associate jet to nearest quark
				if (leadJet.deltaR(quark1) < leadJet.deltaR(quark2)) {
					leadJet.quark = 1;
					if (leadJet.deltaR(quark1) > 0.2) {
						leadJet.match=false;
						numMissIdJets++;
					}
				} else {
					leadJet.quark = 2;
					if (leadJet.deltaR(quark2) > 0.2) {
						leadJet.match=false;
						numMissIdJets++;
					}
				}
			//Subleading jet
				//Associate jet to nearest quark
				if (subLeadJet.deltaR(quark1) < subLeadJet.deltaR(quark2)) {
					subLeadJet.quark = 1;
					if (subLeadJet.deltaR(quark1) > 0.2) {
						subLeadJet.match=false;
						numMissIdJets++;
					}
				} else {
					subLeadJet.quark = 2;
					if (subLeadJet.deltaR(quark2) > 0.2) {
						subLeadJet.match=false;
						numMissIdJets++;
					}
				}

				//Single wrong jet case
				if (numMissIdJets == 1) {

					misIdJets_1++;

					if (!leadJet.match) {

						if (subLeadJet.quark==1) {
							etaVsEta->Fill(quark2.eta,leadJet.eta);
						}else{
							etaVsEta->Fill(quark1.eta,leadJet.eta);
						}

						misIdJet_pt->Fill(leadJet.pt);
						misIdJet_eta->Fill(leadJet.eta);
						misIdJet_phi->Fill(leadJet.phi);

					} else if (!subLeadJet.match) {

						if (leadJet.quark==1) {
							etaVsEta->Fill(quark2.eta,subLeadJet.eta);
						}else{
							etaVsEta->Fill(quark1.eta,subLeadJet.eta);
						}

						misIdJet_pt->Fill(subLeadJet.pt);
						misIdJet_eta->Fill(subLeadJet.eta);
						misIdJet_phi->Fill(subLeadJet.phi);

					}
				}

				//Double wrong jet case
				if (numMissIdJets == 2) {
				
					misIdJets_2++;
					//If they're allocated to the same quark
					if (leadJet.quark == subLeadJet.quark && leadJet.quark==1) {
						if (leadJet.deltaR(quark1) > subLeadJet.deltaR(quark1)) {
							etaVsEta->Fill(quark2.eta,leadJet.eta);
							etaVsEta->Fill(quark1.eta,subLeadJet.eta);
						} else {
							etaVsEta->Fill(quark1.eta,leadJet.eta);
							etaVsEta->Fill(quark2.eta,subLeadJet.eta);
						}
					}else if (leadJet.quark == subLeadJet.quark && leadJet.quark==2) {
                                                if (leadJet.deltaR(quark2) > subLeadJet.deltaR(quark2)) {
                                                        etaVsEta->Fill(quark1.eta,leadJet.eta);
                                                        etaVsEta->Fill(quark2.eta,subLeadJet.eta);
                                                } else {
                                                        etaVsEta->Fill(quark2.eta,leadJet.eta);
                                                        etaVsEta->Fill(quark1.eta,subLeadJet.eta);
                                                }				
					}
					//Different quarks
					if (leadJet.quark==1 && leadJet.quark != subLeadJet.quark) {
						etaVsEta->Fill(quark1.eta,leadJet.eta);
						etaVsEta->Fill(quark2.eta,subLeadJet.eta);
					}else if (leadJet.quark==2 && leadJet.quark != subLeadJet.quark) {
						etaVsEta->Fill(quark2.eta,leadJet.eta);
						etaVsEta->Fill(quark1.eta,subLeadJet.eta);
					}	
				}

			}

	                const   TTHHadronicTag *tthhadronictag = dynamic_cast<const TTHHadronicTag*>(chosenTag);
        	        if(tthhadronictag != NULL) {
				std::cout << "Mass  " << tthhadronictag->diPhoton()->mass() << " GeV"  << std::endl;
				std::cout << "Tagged as TTH Hadronic" << std::endl;
	                	h_m_tthhad->Fill(tthhadronictag->diPhoton()->mass());
			}

	                const   TTHLeptonicTag *tthleptonictag = dynamic_cast<const TTHLeptonicTag*>(chosenTag);
        	        if(tthleptonictag != NULL) {
				std::cout << "Mass  " << tthleptonictag->diPhoton()->mass() << " GeV" << std::endl;
				std::cout << "Tagged as TTH Leptonic" << std::endl;
				h_m_tthlep->Fill(tthleptonictag->diPhoton()->mass());
			}

			const   VHTightTag *vhtighttag = dynamic_cast<const VHTightTag*>(chosenTag);
			if(vhtighttag != NULL) {
				h_m_vhtight->Fill(vhtighttag->diPhoton()->mass());
				std::cout << "Mass  " << vhtighttag->diPhoton()->mass() << " GeV" << std::endl;
				std::cout << "Tagged as VH Tight" << std::endl;
			}

                	const   VHLooseTag *vhloosetag = dynamic_cast<const VHLooseTag*>(chosenTag);
			if(vhloosetag != NULL) {
				h_m_vhloose->Fill(vhloosetag->diPhoton()->mass());
				std::cout << "Mass  " << vhloosetag->diPhoton()->mass() << " GeV" << std::endl;
				std::cout << "Tagged as VH Loose" << std::endl;
			}

			const   VHHadronicTag *vhhadronictag = dynamic_cast<const VHHadronicTag*>(chosenTag);
			if(vhhadronictag != NULL) {
				h_m_vhhad->Fill(vhhadronictag->diPhoton()->mass());
				std::cout << "Mass  " << vhhadronictag->diPhoton()->mass() << " GeV" << std::endl;
				std::cout << "Tagged as VH Hadronic" << std::endl;
			}

			// IMPORTANT: All future Tags must be added in the way of untagged and vbftag.	

			if (untagged == NULL && vbftag == NULL && tthhadronictag == NULL && tthleptonictag == NULL && vhtighttag == NULL && vhloosetag == NULL && vhhadronictag==NULL) {
				std::cout << "[FAILED TO CONVERT TAG] with SumPt " << chosenTag->sumPt() << std::endl;
			}

		} else { //case where TagSorter[0] doesn't exist
			std::cout << "[NO TAG]" <<std::endl;
		}
		std::cout << std::endl;	
	
	} // analyze

	void 
	TagTestAnalyzer::beginJob()
	{
		outputFile_ = new TFile("TagTest.root","RECREATE");
	
		eventCount = 0;
	
		h_m_untagged_0 = new TH1F("untagged 0","Untagged 0 m_{H}",50,100,150);
		h_m_untagged_1 = new TH1F("untagged 1","Untagged 1 m_{H}",50,100,150);
		h_m_untagged_2 = new TH1F("untagged 2","Untagged 2 m_{H}",50,100,150);
		h_m_untagged_3 = new TH1F("untagged 3","Untagged 3 m_{H}",50,100,150);
		h_m_untagged_4 = new TH1F("untagged 4","Untagged 4 m_{H}",50,100,150);
		
		h_m_vbf_0 = new TH1F("vbf 0","VBF 0 m_{H}",50,100,150);
		h_m_vbf_1 = new TH1F("vbf 1","VBF 1 m_{H}",50,100,150);
		h_m_vbf_2 = new TH1F("vbf 2","VBF 2 m_{H}",50,100,150);

		h_m_tthhad = new TH1F("tthhad","ttH Hadronic m_{H}",50,100,150);
		h_m_tthlep = new TH1F("tthlep","ttH Leptonic m_{H}",50,100,150);

		h_m_vhtight = new TH1F("vhtight","VH Tight m_{H}",50,100,150);
		h_m_vhloose = new TH1F("vhloose","VH Loose m_{H}",50,100,150);
		h_m_vhhad = new TH1F("vhhad","VH Hadronic m_{H}",50,100,150);

		deltaR_VBF = new TH1F("jetErrorVBF","VBF Jet E",50,0,2);
		misIdJet_pt = new TH1F("MisIdedJetPt","Misid'ed jet pt",50,0,150);
		misIdJet_eta = new TH1F("MisIdedJetEta","Misid'ed jet eta",50,-6,-6);
		misIdJet_phi = new TH1F("MisIdedJetPhi","Misid'ed jet phi",50,-4,4);
		etaVsEta = new TH2F("trueEtaVsMeasured","True eta vs measured eta",50,-6,6,50,-6,6);
		trueJetEta = new TH1F("trueJetEta","True jet eta",50,-6,6);

		eventTree = new TTree("EventTree","Event Tree");
		eventTree->Branch("pt"     ,&particle.pt      ,"pt/F" );
  		eventTree->Branch("eta"    ,&particle.eta     ,"eta/F");
  		eventTree->Branch("phi"    ,&particle.phi     ,"phi/F");
  		eventTree->Branch("status" ,&particle.status  ,"status/I" );
  		eventTree->Branch("pdgid"  ,&particle.pdgid   ,"pdgid/I");
		eventTree->Branch("rapidity",&particle.rapidity, "rapidity/F"); 
		eventTree->Branch("number" ,&particle.number  ,"number/I");
		
	}

	void
	TagTestAnalyzer::endJob()
	{

		std::cout << "Closing file" << std::endl;

		outputFile_->cd();

		h_m_untagged_0->Write();	
		h_m_untagged_1->Write();
		h_m_untagged_2->Write();
		h_m_untagged_3->Write();
		h_m_untagged_4->Write();

		h_m_vbf_0->Write();
		h_m_vbf_1->Write();
		h_m_vbf_2->Write();

		h_m_tthhad->Write();
		h_m_tthlep->Write();

		h_m_vhtight->Write();
		h_m_vhloose->Write();
		h_m_vhhad->Write();
		
		deltaR_VBF->Write();
		misIdJet_pt->Write();
		misIdJet_eta->Write();
		misIdJet_phi->Write();
		etaVsEta->Write();
		trueJetEta->Write();

		eventTree->Write();
		eventTree->Print();

		TVectorF properties(1);
		properties[0] = eventCount;
		properties.Write("properties");

		outputFile_->Close();

		std::cout << "Number of events: " << properties[0] << std::endl;
		std::cout << "Number of wrong jets\n One wrong: " << misIdJets_1 << " two wrong: " << misIdJets_2 << std::endl; 
	}

	void
	TagTestAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
		//The following says we do not know what parameters are allowed so do no validation
				
		// Please change this to state exactly what you do use, even if it is no parameters
		edm::ParameterSetDescription desc;
		desc.setUnknown();
		descriptions.addDefault(desc);
	}

} // namespace flashgg

typedef flashgg::TagTestAnalyzer FlashggTagTestAnalyzer;
DEFINE_FWK_MODULE(FlashggTagTestAnalyzer);
