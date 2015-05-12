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

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/Ptr.h"

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
  bool  charm;
  int   mother0ID;

  Particle(void){

	pt=0;
	eta=0;
	phi=0;
	pdgid=0;
	status=0;
	number=0;
	rapidity=0;
	match = true;
	charm = false;
 	mother0ID = 0;

  }

  float deltaR(Particle quark) {

	float deltaR(0);
	deltaR += pow((eta-quark.eta),2);
	deltaR += pow((phi-quark.phi),2);

	return sqrt(deltaR);
  
  }

  void print(void) {

    std::cout << std::setw(12) << pt;
    std::cout << std::setw(12) << eta;
    std::cout << std::setw(12) << phi;
    std::cout << std::setw(12) << status;
    std::cout << std::setw(12) << pdgid;
    std::cout << std::setw(12) << rapidity;
    std::cout << std::setw(12) << mother0ID;
    std::cout << std::endl;

  }

  void charmTag(void) {

    TString pdgidString;
    pdgidString.Form("%d",pdgid);
    pdgidString = pdgidString(0,2);

    if (pdgidString.Contains("4") && !pdgidString.Contains("14")) {
      charm = true;
    } else {
      charm = false;
    }

  }



};

using namespace std;
using namespace edm;

// **********************************************************************

namespace flashgg {

	struct VBFTruth {

		edm::Ptr<reco::GenJet> closestGenJetToLeadingJet;
		edm::Ptr<reco::GenJet> closestGenJetToSubLeadingJet;
		edm::Ptr<reco::GenJet> closestGenJetToLeadingPhoton;
		edm::Ptr<reco::GenJet> closestGenJetToSubLeadingPhoton;

		edm::Ptr<reco::GenParticle> closestParticleToLeadingJet;
		edm::Ptr<reco::GenParticle> closestParticleToSubLeadingJet;
		edm::Ptr<reco::GenParticle> closestParticleToLeadingPhoton;
		edm::Ptr<reco::GenParticle> closestParticleToSubLeadingPhoton;

		edm::Ptr<reco::GenParticle> leadingQuark;
		edm::Ptr<reco::GenParticle> subLeadingQuark;

	};

	bool isCharmed(int pdgId) {

		bool charm;
		TString pdgIdString;

		pdgIdString.Form("%d",pdgId);
		pdgIdString = pdgIdString(0,2);

		if (pdgIdString.Contains("4") && !pdgIdString.Contains("14")) {
			charm = true;
		} else {
			charm = false;
		}

		return charm;

	}

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
			int charmCount;
			int VBFCount;

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
			
			TH1F * misIdJet_pt;
			TH1F * misIdJet_eta;
			TH1F * misIdJet_phi;
			TH2F * etaVsEta;	
			TH1F * misIdJet_dEta;

			EDGetTokenT< edm::View<reco::GenParticle> > genPartToken_;
			EDGetTokenT< edm::View<reco::GenJet> > genJetToken_;
			EDGetTokenT< edm::View<reco::Vertex> > vertexToken_;
			EDGetTokenT< edm::View<pat::PackedGenParticle> > pgParticleToken_;

			Int_t eventNumber;
		
			//Particle variables	
			bool charm_p;
			float pt_p, eta_p, phi_p, rapidity_p;
			int   status_p, pdgId_p, nMothers_p, nDaughters_p;
			//Jet variables
			bool charm_j, quarkMatch;
			float pt_j, eta_j, phi_j, pdgId_j, rapidity_j;

			//Counts
			int oppositeSign, oneCharmJet, twoCharmJet, charmedVBF;
			
			TTree * eventTree;
			TTree * jetTree;
			TTree * VBFTree;
			TTree * quarkTree;
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
		genPartToken_ (consumes<View<reco::GenParticle> >(iConfig.getUntrackedParameter<InputTag> ("GenParticleTag", InputTag("prunedGenParticles")))),
		genJetToken_ (consumes<View<reco::GenJet> >(iConfig.getUntrackedParameter<InputTag> ("GenJetTag", InputTag("slimmedGenJets")))),
		vertexToken_ (consumes<View<reco::Vertex> >(iConfig.getUntrackedParameter<InputTag> ("VertexTag", InputTag("offlineSlimmedPrimaryVertices")))),
		pgParticleToken_ (consumes<edm::View<pat::PackedGenParticle> >(iConfig.getUntrackedParameter<InputTag>("PackedGenParticleTag", InputTag("packedGenParticles"))))
	{
		eventNumber = 0;
	}

	TagTestAnalyzer::~TagTestAnalyzer()
	{
	}

	void
	TagTestAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
	
		eventCount++;

		//GenParticles
		Handle<View<reco::GenParticle> > genParticles;
  		iEvent.getByToken(genPartToken_,genParticles);
		//GenJets	
		Handle<View<reco::GenJet> > genJets;
		iEvent.getByToken(genJetToken_,genJets);
		//Vertices
		Handle<View<reco::Vertex> > primaryVertices;
		iEvent.getByToken(vertexToken_,primaryVertices);
		//PackedGenParticles
		Handle<View<pat::PackedGenParticle> > pgParticles;
		iEvent.getByToken(pgParticleToken_,pgParticles);		
	
		std::cout << "TESTING\nVertices:" << std::endl;
		for (unsigned int vertLoop = 0; vertLoop < primaryVertices->size(); vertLoop++) {

			edm::Ptr<reco::Vertex> vertex = primaryVertices->ptrAt(vertLoop);

			std::cout << setw(12) << vertex->x();
			std::cout << setw(12) << vertex->y();
			std::cout << setw(12) << vertex->z();
			std::cout << std::endl;

		}

		std::cout << "TESTING\nPackedGenParticles:" << std::endl;
		for (unsigned int pgLoop = 0; pgLoop < pgParticles->size(); pgLoop++) {

			edm::Ptr<pat::PackedGenParticle> pgp = pgParticles->ptrAt(pgLoop);

			std::cout << setw(12) << pgp->vertex().x();
			std::cout << setw(12) << pgp->vertex().y();
			std::cout << setw(12) << pgp->vertex().z();
			std::cout << setw(12) << pgp->pt();
			std::cout << setw(12) << pgp->eta();
			std::cout << setw(12) << pgp->phi();
			std::cout << setw(12) << pgp->pdgId();
			std::cout << setw(12) << pgp->status();
			std::cout << std::endl;

		}

		std::cout << "GenParticles:" << std::endl;
		for (unsigned int partLoop = 0; partLoop < genParticles->size(); partLoop++) {

			edm::Ptr<reco::GenParticle> part = genParticles->ptrAt(partLoop);		


			std::cout << setw(12) << part->vertex().x();
			std::cout << setw(12) << part->vertex().y();
			std::cout << setw(12) << part->vertex().z();
			std::cout << std::endl;

			pt_p = part->pt();			
			eta_p = part->eta();			
			phi_p = part->phi();			
			status_p = part->status();			
			pdgId_p = part->pdgId();			
			rapidity_p = part->rapidity();			
			nMothers_p = part->numberOfMothers();			
			nDaughters_p = part->numberOfDaughters();			
			charm_p = isCharmed(pdgId_p);
		
			if (nMothers_p>0) {eventTree->Fill();}

		}

		for (unsigned int jetLoop = 0; jetLoop < genJets->size(); jetLoop++) {

			edm::Ptr<reco::GenJet> jet = genJets->ptrAt(jetLoop);		

			pt_j  = jet->pt();			
			eta_j = jet->eta();			
			phi_j = jet->phi();			
			rapidity_j = jet->rapidity();			

			jetTree->Fill();
		}


		//Tagging
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
				
				int category = vbftag->categoryNumber();
				std::cout << "Tagged as VBF - category " << category << std::endl;

				if (category == 0) {

					h_m_vbf_0->Fill(vbftag->diPhoton()->mass()); 
									
				}else if (category == 1) {

                                        h_m_vbf_1->Fill(vbftag->diPhoton()->mass());
						
                                }else if (category == 2) {

                                        h_m_vbf_2->Fill(vbftag->diPhoton()->mass());
						
				}

			//VBF Truth
				float dr_leadjet=999., dr_subLeadjet=999.;
				float dr_leadPhoton=999., dr_subLeadPhoton=999.;
				unsigned int index_leadJet(0), index_subLeadJet(0), index_leadPhoton(0), index_subLeadPhoton(0);
				float dr;
				VBFTruth truth;

			//Particles
				//Checking VBF jets and VBF diphoton candidates against GenParticles 
				for (unsigned int partLoop = 0; partLoop < genParticles->size(); partLoop++) {

					edm::Ptr<reco::GenParticle> part = genParticles->ptrAt(partLoop);		

					if (part->status() == 3) {

						dr = deltaR(vbftag->leadingJet().eta(),vbftag->leadingJet().phi(),part->eta(),part->phi());
						if (dr < dr_leadjet) {
                        				dr_leadjet = dr;
							index_leadJet = partLoop;
                    				}

						dr = deltaR(vbftag->subLeadingJet().eta(),vbftag->subLeadingJet().phi(),part->eta(),part->phi());
						if (dr < dr_subLeadjet) {
                        				dr_subLeadjet = dr;
							index_subLeadJet = partLoop;
                    				}

						dr = deltaR(vbftag->diPhoton()->leadingPhoton()->eta(),vbftag->diPhoton()->leadingPhoton()->phi(),part->eta(),part->phi());
						if (dr < dr_leadPhoton) {
							dr_leadPhoton = dr;
							index_leadPhoton = partLoop;
						}

						dr = deltaR(vbftag->diPhoton()->subLeadingPhoton()->eta(),vbftag->diPhoton()->subLeadingPhoton()->phi(),part->eta(),part->phi());
						if (dr < dr_subLeadPhoton) {
							dr_subLeadPhoton = dr;
							index_subLeadPhoton = partLoop;
						}

					}

				}

				//Store the GenParticle candidates
				truth.closestParticleToLeadingJet = genParticles->ptrAt(index_leadJet);
				truth.closestParticleToSubLeadingJet =genParticles->ptrAt(index_subLeadJet);
				truth.closestParticleToLeadingPhoton =genParticles->ptrAt(index_leadPhoton);
				truth.closestParticleToSubLeadingPhoton =genParticles->ptrAt(index_subLeadPhoton);

			//Jets
				//Checking VBF jets and VBF diphoton candidates against GenJets
				for (unsigned int jetLoop = 0; jetLoop < genJets->size(); jetLoop++) {

					edm::Ptr<reco::GenJet> jet = genJets->ptrAt(jetLoop);
					dr = deltaR(vbftag->leadingJet().eta(),vbftag->leadingJet().phi(),jet->eta(),jet->phi());
					if (dr < dr_leadjet) {
                        			dr_leadjet = dr;
						index_leadJet = jetLoop;
                    			}
					dr = deltaR(vbftag->subLeadingJet().eta(),vbftag->subLeadingJet().phi(),jet->eta(),jet->phi());
					if (dr < dr_subLeadjet) {
                        			dr_subLeadjet = dr;
						index_subLeadJet = jetLoop;
                    			}
					dr = deltaR(vbftag->diPhoton()->leadingPhoton()->eta(),vbftag->diPhoton()->leadingPhoton()->phi(),jet->eta(),jet->phi());
					if (dr < dr_leadPhoton) {
						dr_leadPhoton = dr;
						index_leadPhoton = jetLoop;
					}
					dr = deltaR(vbftag->diPhoton()->subLeadingPhoton()->eta(),vbftag->diPhoton()->subLeadingPhoton()->phi(),jet->eta(),jet->phi());
					if (dr < dr_subLeadPhoton) {
						dr_subLeadPhoton = dr;
						index_subLeadPhoton = jetLoop;
					}
					

				}

				//Store the GenJet candidates
				truth.closestGenJetToLeadingJet = genJets->ptrAt(index_leadJet);
				truth.closestGenJetToSubLeadingJet =genJets->ptrAt(index_subLeadJet);
				truth.closestGenJetToLeadingPhoton =genJets->ptrAt(index_leadPhoton);
				truth.closestGenJetToSubLeadingPhoton =genJets->ptrAt(index_subLeadPhoton);

			//Leading quarks	
				int index_leadQ(0), index_subLeadQ(0);
				float pt_leadQ(0), pt_subLeadQ(0);

				for (unsigned int partLoop = 0; partLoop < genParticles->size(); partLoop++) {
										
					edm::Ptr<reco::GenParticle> part = genParticles->ptrAt(partLoop);		
				
					if (part->status() == 3 && part->numberOfMothers()>1) {
						if (abs(part->pdgId()) <= 5) {
							if (part->pt() > pt_leadQ) {
								index_subLeadQ = index_leadQ;
								pt_subLeadQ = pt_leadQ;
								index_leadQ = partLoop;
								pt_leadQ = part->pt();
							} else if (part->pt() > pt_subLeadQ) {
								index_subLeadQ = partLoop;
								pt_subLeadQ = part->pt();
							}
						}
					}
				}

				//Store the quark candidates
				truth.leadingQuark = genParticles->ptrAt(index_leadQ);
				truth.subLeadingQuark = genParticles->ptrAt(index_subLeadQ);
				
				//Add to quark tree
				pt_p = truth.leadingQuark->pt();			
				eta_p = truth.leadingQuark->eta();			
				phi_p = truth.leadingQuark->phi();			
				status_p = truth.leadingQuark->status();			
				pdgId_p = truth.leadingQuark->pdgId();			
				rapidity_p = truth.leadingQuark->rapidity();			
				nMothers_p = truth.leadingQuark->numberOfMothers();			
				nDaughters_p = truth.leadingQuark->numberOfDaughters();			
				charm_p = isCharmed(pdgId_p);
				quarkTree->Fill();

				pt_p = truth.subLeadingQuark->pt();			
				eta_p = truth.subLeadingQuark->eta();			
				phi_p = truth.subLeadingQuark->phi();			
				status_p = truth.subLeadingQuark->status();			
				pdgId_p = truth.subLeadingQuark->pdgId();			
				rapidity_p = truth.subLeadingQuark->rapidity();			
				nMothers_p = truth.subLeadingQuark->numberOfMothers();			
				nDaughters_p = truth.subLeadingQuark->numberOfDaughters();			
				charm_p = isCharmed(pdgId_p);
				quarkTree->Fill();

			//Quark/VBF Jet matching
				std::cout << "VBF Jets" << std::endl;
				std::cout << setw(12) << "Pt" << setw(12) << "Eta" << setw(12) << "Phi" << std::endl;
				std::cout << setw(12) << vbftag->leadingJet().pt() << setw(12) << vbftag->leadingJet().eta();
				std::cout << setw(12) << vbftag->leadingJet().phi() << std::endl;
				std::cout << setw(12) << vbftag->subLeadingJet().pt() << setw(12) << vbftag->subLeadingJet().eta();
				std::cout << setw(12) << vbftag->subLeadingJet().phi() << std::endl;

				std::cout << "VBF Jet/Quark matching\nQuarks:" << std::endl;
				std::cout << setw(12) << "DeltaR" << setw(12) << "Pt" << setw(12) << "Eta" << setw(12) << "Phi" << std::endl;
				float dr1, dr2;
				bool leadVBFMatch=false, subLeadVBFMatch=false;

				//lead jet
				dr1 = deltaR(vbftag->leadingJet().eta(),vbftag->leadingJet().phi(),truth.leadingQuark->eta(),truth.leadingQuark->phi());	
				dr2 = deltaR(vbftag->leadingJet().eta(),vbftag->leadingJet().phi(),truth.subLeadingQuark->eta(),truth.subLeadingQuark->phi());
				charm_j = false;		
	
				if (dr1<dr2) {

					std::cout << setw(12) << dr1;
					std::cout << setw(12) << truth.leadingQuark->pt();
					std::cout << setw(12) << truth.leadingQuark->eta();
					std::cout << setw(12) << truth.leadingQuark->phi();

					if (isCharmed(truth.leadingQuark->pdgId())) {
						std::cout << setw(12) << truth.leadingQuark->pdgId() << " (Charmed)" << std::endl;
					}else{
						std::cout << setw(12) << truth.leadingQuark->pdgId() << std::endl;
					}
					if (dr1 > 0.2) {
						misIdJet_pt->Fill(vbftag->leadingJet().pt());
						misIdJet_eta->Fill(vbftag->leadingJet().eta());
						misIdJet_dEta->Fill(fabs(vbftag->leadingJet().eta()-truth.leadingQuark->eta()));
					}
					if (dr1 < 0.2) {leadVBFMatch = true;charm_j=isCharmed(truth.leadingQuark->pdgId());}

				} else {

					std::cout << setw(12) << dr2;
					std::cout << setw(12) << truth.subLeadingQuark->pt();
					std::cout << setw(12) << truth.subLeadingQuark->eta();
					std::cout << setw(12) << truth.subLeadingQuark->phi();

					if (isCharmed(truth.subLeadingQuark->pdgId())) {
						std::cout << setw(12) << truth.subLeadingQuark->pdgId() << " (Charmed)" << std::endl;
					}else{
						std::cout << setw(12) << truth.subLeadingQuark->pdgId() << std::endl;
					}
					if (dr2 > 0.2) {
						misIdJet_pt->Fill(vbftag->leadingJet().pt());
						misIdJet_eta->Fill(vbftag->leadingJet().eta());
						misIdJet_dEta->Fill(fabs(vbftag->leadingJet().eta()-truth.subLeadingQuark->eta()));
					}
					if (dr2 < 0.2) {leadVBFMatch = true;charm_j=isCharmed(truth.subLeadingQuark->pdgId());}

				}
				//Add jets to the vbf tree				
				pt_j  = vbftag->leadingJet().pt();			
				eta_j  = vbftag->leadingJet().eta();			
				phi_j  = vbftag->leadingJet().phi();			
				rapidity_j  = vbftag->leadingJet().rapidity();			
				quarkMatch = leadVBFMatch;
				VBFTree->Fill();

				//sublead jet
				dr1 = deltaR(vbftag->subLeadingJet().eta(),vbftag->subLeadingJet().phi(),truth.leadingQuark->eta(),truth.leadingQuark->phi());
				dr2 = deltaR(vbftag->subLeadingJet().eta(),vbftag->subLeadingJet().phi(),truth.subLeadingQuark->eta(),truth.subLeadingQuark->phi());
				charm_j = false;

				if (dr1<dr2) {

					std::cout << setw(12) << dr1;
					std::cout << setw(12) << truth.leadingQuark->pt();
					std::cout << setw(12) << truth.leadingQuark->eta();
					std::cout << setw(12) << truth.leadingQuark->phi();

					if (isCharmed(truth.leadingQuark->pdgId())) {
						std::cout << setw(12) << truth.leadingQuark->pdgId() << " (Charmed)" << std::endl;
					}else{
						std::cout << setw(12) << truth.leadingQuark->pdgId() << std::endl;
					}
					if (dr1 > 0.2) {
						misIdJet_pt->Fill(vbftag->subLeadingJet().pt());
						misIdJet_eta->Fill(vbftag->subLeadingJet().eta());
						misIdJet_dEta->Fill(fabs(vbftag->subLeadingJet().eta()-truth.leadingQuark->eta()));
					}
					if (dr1 < 0.2) {subLeadVBFMatch = true;charm_j=isCharmed(truth.leadingQuark->pdgId());}
						
				} else {

					std::cout << setw(12) << dr2;
					std::cout << setw(12) << truth.subLeadingQuark->pt();
					std::cout << setw(12) << truth.subLeadingQuark->eta();
					std::cout << setw(12) << truth.subLeadingQuark->phi();

					if (isCharmed(truth.subLeadingQuark->pdgId())) {
						std::cout << setw(12) << truth.subLeadingQuark->pdgId() << " (Charmed)" << std::endl;
					}else{
						std::cout << setw(12) << truth.subLeadingQuark->pdgId() << std::endl;
					}
					if (dr2 > 0.2) {
						misIdJet_pt->Fill(vbftag->subLeadingJet().pt());
						misIdJet_eta->Fill(vbftag->subLeadingJet().eta());
						misIdJet_dEta->Fill(fabs(vbftag->subLeadingJet().eta()-truth.subLeadingQuark->eta()));
					}
					if (dr2 < 0.2) {subLeadVBFMatch = true;charm_j=isCharmed(truth.subLeadingQuark->pdgId());}

				}

				pt_j  = vbftag->subLeadingJet().pt();			
				eta_j  = vbftag->subLeadingJet().eta();			
				phi_j  = vbftag->subLeadingJet().phi();			
				rapidity_j  = vbftag->subLeadingJet().rapidity();			
				quarkMatch = subLeadVBFMatch;
				VBFTree->Fill();
				
				std::cout << "Other Jets: " << std::endl;
				std::cout << setw(12) << "DeltaR LQ" << setw(12) << "DeltaR SLQ" << setw(12) << "Pt" << setw(12) << "Eta" << setw(12) << "Phi" << std::endl;
				for (unsigned int jetLoop = 0; jetLoop < genJets->size(); jetLoop++) {

					edm::Ptr<reco::GenJet> jet = genJets->ptrAt(jetLoop);
					dr = deltaR(jet->eta(),jet->phi(),truth.leadingQuark->eta(),truth.leadingQuark->phi());
					std::cout << setw(12) << dr;
					dr = deltaR(jet->eta(),jet->phi(),truth.subLeadingQuark->eta(),truth.subLeadingQuark->phi());
					std::cout << setw(12) << dr;
					std::cout << setw(12) << jet->pt();
					std::cout << setw(12) << jet->eta();
					std::cout << setw(12) << jet->phi() << std::endl;
				
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
		std::cout << "DEBUG1" << std::endl;
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

		misIdJet_pt = new TH1F("MisIdedJetPt","Misid'ed jet pt",50,0,150);
		misIdJet_eta = new TH1F("MisIdedJetEta","Misid'ed jet eta",50,-6,-6);
		misIdJet_phi = new TH1F("MisIdedJetPhi","Misid'ed jet phi",50,-4,4);
		etaVsEta = new TH2F("trueEtaVsMeasured","True eta vs measured eta",50,-6,6,50,-6,6);
		misIdJet_dEta = new TH1F("MisIdedJetdEta","Misid'ed jet dEta",50,-6,-6);

		eventTree = new TTree("EventTree","Event Tree");
		eventTree->Branch("pt"     ,&pt_p      ,"pt/F" );
  		eventTree->Branch("eta"    ,&eta_p     ,"eta/F");
  		eventTree->Branch("phi"    ,&phi_p     ,"phi/F");
  		eventTree->Branch("status" ,&status_p  ,"status/I" );
  		eventTree->Branch("pdgid"  ,&pdgId_p   ,"pdgid/I");
		eventTree->Branch("rapidity",&rapidity_p, "rapidity/F"); 
		eventTree->Branch("event" ,&eventCount,  "event/I");
		eventTree->Branch("charm", &charm_p, "charm/O");
		eventTree->Branch("numMothers" ,&nMothers_p,  "numMothers/I");
		eventTree->Branch("numDaughters" ,&nDaughters_p,  "numDaughters/I");

		jetTree = new TTree("JetTree","Jet Tree");
		jetTree->Branch("pt", &pt_j, "pt/F");
		jetTree->Branch("eta", &eta_j, "eta/F");
		jetTree->Branch("phi", &phi_j, "phi/F");
		jetTree->Branch("event", &eventCount, "event/I");

		VBFTree = new TTree("VBFTree","VBF Jets Tree");
		VBFTree->Branch("pt", &pt_j, "pt/F");
		VBFTree->Branch("eta", &eta_j, "eta/F");
		VBFTree->Branch("phi", &phi_j, "phi/F");
		VBFTree->Branch("charm", &charm_j, "charm/O");
		VBFTree->Branch("event", &eventCount, "event/I");
		VBFTree->Branch("charm", &charm_j, "charm/O");
		VBFTree->Branch("quarkMatch", &quarkMatch, "quarkMatch/O");

		quarkTree = new TTree("quarkTree","quark Jets Tree");
		quarkTree->Branch("pt", &pt_p, "pt/F");
		quarkTree->Branch("eta", &eta_p, "eta/F");
		quarkTree->Branch("phi", &phi_p, "phi/F");
		quarkTree->Branch("charm", &charm_p, "charm/O");
		quarkTree->Branch("event", &eventCount, "event/I");
		quarkTree->Branch("pdgId", &pdgId_p, "pdgId/I");
		quarkTree->Branch("numMothers" ,&nMothers_p,  "numMothers/I");
		quarkTree->Branch("numDaughters" ,&nDaughters_p,  "numDaughters/I");
	}

	void
	TagTestAnalyzer::endJob()
	{

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
		
		misIdJet_pt->Write();
		misIdJet_eta->Write();
		misIdJet_phi->Write();
		etaVsEta->Write();
		misIdJet_dEta->Write();

		eventTree->Write();
		eventTree->Print();

		jetTree->Write();
		jetTree->Print();

		VBFTree->Write();
		VBFTree->Print();

		quarkTree->Write();
		quarkTree->Print();

		TVectorF properties(1);
		properties[0] = eventCount;
		properties.Write("properties");

		outputFile_->Close();

		std::cout << "Total number of events: " << properties[0] << std::endl;
		
		std::cout << "Number of wrong VFB jets\n One wrong: " << misIdJets_1 << " two wrong: " << misIdJets_2 << std::endl; 
		
		std::cout << "Number tagged as VBF: " << VBFCount << std::endl;
		std::cout << "Number of VBF with a charm jet: " << charmCount << std::endl;

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
