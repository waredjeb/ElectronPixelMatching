// -*- C++ -*-
//
// Package:    SeedFromTrack/SeedFromTrackAnalyzer
// Class:      SeedFromTrackAnalyzer
//
/**\class SeedFromTrackAnalyzer SeedFromTrackAnalyzer.cc SeedFromTrack/SeedFromTrackAnalyzer/plugins/SeedFromTrackAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Wahid Redjeb
//         Created:  Wed, 29 Jul 2020 09:31:23 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetType.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include <stdexcept>
// #include "SeedFromTrack/SeedFromTrackAnalyzer/src/ElectronMatch.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"

#include <iostream>
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include <cmath>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// // This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class ElectronMatchSeedNew : public edm::EDAnalyzer {
    public:
        explicit ElectronMatchSeedNew(const edm::ParameterSet& iConfig);
        virtual ~ElectronMatchSeedNew(){};

    private:
        //----edm control---
        virtual void beginJob();
        virtual void beginRun(edm::Run const&, edm::EventSetup const&);
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob();
        virtual void endRun(edm::Run const&, edm::EventSetup const&);
        virtual void clear();
        


        // check on lumi
        std::vector<edm::LuminosityBlockRange> lumisTemp ;
        const edm::EDGetTokenT<reco::GenParticleCollection> gensToken_;              //GenParticles
        const edm::EDGetTokenT<reco::ElectronCollection> electronCollection;  //GenParticles
        //basic infos on the event
        unsigned long long int event_;
        int       run_;
        int       lumi_;
        double DeltaR_;

        //Important objects;


        //booleans
        bool verbose_;

        //Tree
        std::string t;
        TTree* tree_;
        bool seed_;

        //Kinematics
        std::vector<double>* eR_pt = new std::vector<double>;
        std::vector<double>* eR_eta = new std::vector<double>;
        std::vector<double>* eR_phi = new std::vector<double>;
        std::vector<double>* eG_pt = new std::vector<double>;
        std::vector<double>* eG_eta = new std::vector<double>;
        std::vector<double>* eG_phi = new std::vector<double>;
        std::vector<double>* eDouble = new std::vector<double>;
        std::vector<double>* eTrip = new std::vector<double>;
        std::vector<int>* eventPassed = new std::vector<int>;
        std::vector<int>* eventIdNotPassed = new std::vector<int>;
        std::vector<int>* eventIdPassed = new std::vector<int>;
        std::vector<int>* eMatched = new std::vector<int>;
        std::vector<int>* electronReco = new std::vector<int>;
        std::vector<std::vector<int>>* eMatchedBPIX = new std::vector<std::vector<int>>;
        std::vector<std::vector<int>>* eMatchedFPIX = new std::vector<std::vector<int>>;
        std::vector<double>* eSim_pt = new std::vector<double>;
        std::vector<double>* eSim_eta = new std::vector<double>;
        std::vector<double>* eSim_phi = new std::vector<double>;
        // double totMatchedNew = 0;
        // double totMatchedNew = 0;

    

};

using namespace reco;

ElectronMatchSeedNew::ElectronMatchSeedNew(const edm::ParameterSet& iConfig): 
     gensToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genP"))),
      electronCollection(
          consumes<reco::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
      DeltaR_(iConfig.getParameter<double>("DeltaR")),
      verbose_(iConfig.getParameter<bool>("verbose")),
      t(iConfig.getParameter<std::string>("tree")),
      seed_(iConfig.getParameter<bool>("seed"))

{   
};

void ElectronMatchSeedNew::beginJob()
{   
    edm::Service<TFileService> fs;
    tree_ = fs -> make<TTree>(t.c_str(), t.c_str());

    tree_->Branch("eR_pt", &eR_pt);
    tree_->Branch("eR_eta", &eR_eta);
    tree_->Branch("eR_phi", &eR_phi);
    tree_->Branch("eG_pt", &eG_pt);
    tree_->Branch("eG_eta", &eG_eta);
    tree_->Branch("eG_phi", &eG_phi);
    tree_->Branch("eventPassed", &eventPassed);
    tree_->Branch("eventIdPassed", &eventIdPassed);
    tree_->Branch("eventIdNotPassed", &eventIdNotPassed);
    tree_->Branch("eDouble", &eDouble);
    tree_->Branch("eTrip", &eTrip);
    tree_->Branch("eMatched", &eMatched);
    tree_->Branch("electronReco", &electronReco);
    tree_->Branch("BPIX", &eMatchedBPIX);
    tree_->Branch("FPIX", &eMatchedFPIX);
    tree_->Branch("eSim_pt", &eSim_pt); 
    tree_->Branch("eSim_eta", &eSim_eta); 
    tree_->Branch("eSim_phi", &eSim_phi); 
    return;
}


void ElectronMatchSeedNew::clear(){

    eR_pt->clear();
    eR_eta->clear();
    eR_phi->clear();
    eG_pt->clear();
    eG_eta->clear();
    eG_phi->clear();  
    eMatched->clear();
    eventPassed->clear();
    eventIdPassed->clear();
    eventIdNotPassed->clear();
    eDouble->clear();
    eTrip->clear();
    eMatchedBPIX->clear();
    eMatchedFPIX->clear();
    electronReco->clear();
    eSim_pt->clear();
    eSim_eta->clear();
    eSim_phi->clear(); 

    return;
}

reco::GenParticle get_lastcopy_prefsrSN(reco::GenParticle part) {
    auto daughters = part.daughterRefVector();
    if (daughters.size() == 1 && daughters.at(0)->pdgId() == part.pdgId()) {
        return get_lastcopy_prefsrSN(*(daughters.at(0)));
    } else {
        return part;
    }
}

reco::GenParticle get_lastcopySN(reco::GenParticle part) {
    auto daughters = part.daughterRefVector();
    for (size_t p = 0; p != daughters.size(); ++p) {
        if (daughters.at(p)->pdgId() == part.pdgId()) {
        return get_lastcopySN(*(daughters.at(p)));
        }
    }
    return part;
    }

reco::GenParticleCollection get_genpartsSN(reco::GenParticleCollection genparts,
                                            int pid = 11,
                                            bool antipart = true,
                                            int status = 2) {
    std::vector<reco::GenParticle> selected;
    for (size_t i = 0; i != genparts.size(); ++i) {
        auto part = genparts.at(i);
        auto pdg_id = part.pdgId();
        if (pdg_id == pid || (antipart && abs(pdg_id) == abs(pid))) {
        if(part.mother()->pdgId() == 23){
            if (part.isHardProcess()) {
            if (status == 1) {
                selected.push_back(part);
            } else if (status == 2) {
                selected.push_back(get_lastcopy_prefsrSN(part));
            } else if (status == 3) {
                selected.push_back(get_lastcopySN(part));
            } else {
                throw std::runtime_error("error status not implemented");
            }
            }
        }
        } else {
        if (part.numberOfMothers() == 0) {
            selected.push_back(part);
        }
        }
    }
    return selected;
    }

std::pair<reco::GenParticle*, int> match_to_genSN(double eta_reco,
                                                    double phi_reco,
                                                    reco::GenParticleCollection genParticles,
                                                    int pid = 11,
                                                    bool antipart = true,
                                                    double max_dr = 0.1) {
    reco::GenParticle* best_match = nullptr;
    double best_dr2 = max_dr * max_dr;
    auto selected_parts = get_genpartsSN(genParticles);
    // std::cout << "Number selected particles: " << selected_parts.size() << std::endl;
    for (size_t i = 0; i != selected_parts.size(); ++i) {
        auto s_part = selected_parts.at(i);
        auto dr2 = std::pow((eta_reco - s_part.eta()), 2) + std::pow((phi_reco - s_part.phi()), 2);
        if (dr2 < best_dr2) {
        best_match = &(selected_parts.at(i));
        best_dr2 = dr2;
        }
    }
    // return best_match
    return std::make_pair(best_match, selected_parts.size());
    }

reco::GenParticle* matchGenSN(reco::Electron electron, const reco::GenParticleCollection& genParticle ){
    auto eta = electron.eta();
    auto phi = electron.phi();
    auto matchedEle = match_to_genSN(eta,phi, genParticle).first;
    return matchedEle;
  }


void ElectronMatchSeedNew::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{   
    using namespace edm;
    using namespace reco;

    event_ = iEvent.id().event();
    run_   = iEvent.id().run();
    lumi_  = iEvent.luminosityBlock();


    clear();

    edm::Handle<reco::GenParticleCollection> genParticleH;
    edm::Handle<reco::ElectronCollection> gsfElectron;
    
    iEvent.getByToken(gensToken_, genParticleH);
    const reco::GenParticleCollection& genParticle = (*genParticleH.product());
    iEvent.getByToken(electronCollection, gsfElectron);

    edm::ESHandle<TrackerGeometry> pDD;
    iSetup.get<TrackerDigiGeometryRecord>().get(pDD);
    // const TrackerGeometry &hpDD = *pDD;
    edm::ESHandle<TrackerTopology> httopo;
    iSetup.get<TrackerTopologyRcd>().get(httopo);
    const TrackerTopology &ttopo = *httopo;
    


    if(gsfElectron.isValid()){
      eventIdPassed->push_back(event_);
      
      typename reco::ElectronCollection::const_iterator gsfIter(gsfElectron->begin());

      for(; gsfIter != gsfElectron->end(); ++gsfIter){
          eventPassed->push_back(0); //true_e
          electronReco->push_back(0); //true
        auto seed = gsfIter->gsfTrack()->seedRef();
        if(verbose_){
            std::cout << "New Electron Reco" << std::endl;
        }

        auto matchedEle = matchGenSN(*gsfIter, genParticle);
        
        if(matchedEle != nullptr){
            if(verbose_){
                std::cout << "Matched New electron" << std::endl;
            }
          // eLegMatched->push_back(0);
          eMatched->push_back(0);//true
          eG_pt->push_back(matchedEle->pt());
          eG_phi->push_back(matchedEle->phi());
          eG_eta->push_back(matchedEle->eta());
          eR_pt->push_back(gsfIter->pt());
          eR_phi->push_back(gsfIter->phi());
          eR_eta->push_back(gsfIter->eta());
          auto gsfTrack = gsfIter->gsfTrack();
          auto seed = gsfTrack->seedRef();
          auto rhits = seed->recHits();
        //   pixelSubDet(*gsfTrack, hpDD, ttopo );
          std::vector<int>* bpixlayer = new std::vector<int>;
          std::vector<int>* fpixlayer = new std::vector<int>;

            for(auto const& rhit: rhits){
                //   auto rhit = *iH;
                  if((rhit.isValid() && rhit.det() != nullptr)){
                    //   std::cout << "Rhit valid" << std::endl;
                    DetId det = rhit.geographicalId();
                    int subdet = rhit.det()->geographicalId().subdetId();
                    auto trkSubDet = pDD->geomDetSubDetector(subdet);
                    if(GeomDetEnumerators::isTrackerPixel(trkSubDet)){
                        if(subdet == PixelSubdetector::PixelBarrel){
                            // std::cout << "Barrel Layer: " << ttopo.pxbLayer(det) << std::endl;
                            bpixlayer->push_back(ttopo.pxbLayer(det));
                        }
                        if(subdet == PixelSubdetector::PixelEndcap){
                            // std::cout << "Endcap disk: " << ttopo.pxfDisk(det) << std::endl;
                            fpixlayer->push_back(ttopo.pxfDisk(det));
                        }
                    }
                  }
            }
            eMatchedBPIX->push_back(*bpixlayer);
            eMatchedFPIX->push_back(*fpixlayer);
            bpixlayer->clear();
            fpixlayer->clear();

          if(seed_){
            if(seed->nHits() == 2){
                eDouble->push_back(1);
                eTrip->push_back(0);
            }
            else{
                eDouble->push_back(0);
                eTrip->push_back(1);
            }
            }
        }
        else{
        eMatched->push_back(1);//false
        eG_pt->push_back(-999);
        eG_phi->push_back(-999);
        eG_eta->push_back(-999);
        eR_pt->push_back(-999);
        eR_phi->push_back(-999);
        eR_eta->push_back(-999);
        eDouble->push_back(0);
        eTrip->push_back(0);
        std::vector<int> fail= {-1};
        eMatchedBPIX->push_back(fail);
        eMatchedFPIX->push_back(fail);
        }
      }
    }
    else{
      eventPassed->push_back(1);//false
      eventIdNotPassed->push_back(event_);
      
    }

    tree_->Fill();


};

//------------------Starting Event-----------------------------------------

void ElectronMatchSeedNew::beginRun(edm::Run const& iRun, edm::EventSetup const& iEvent)
{

};

void ElectronMatchSeedNew::endJob(){

};

void ElectronMatchSeedNew::endRun(edm::Run const& iRun, edm::EventSetup const& iEvent){

};

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(ElectronMatchSeedNew);

