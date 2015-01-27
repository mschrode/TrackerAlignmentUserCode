// -*- C++ -*-
//
// Package:    VertexAnalysis
// Class:      VertexAnalysis
// 
/**\class VertexAnalysis VertexAnalysis.cc Alignment/VertexAnalysis/src/VertexAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Schroeder,32 3-B20,+41227677557,
//         Created:  Mo 12. Jan 16:46:27 CET 2015
// $Id$
//
//


// system include files
#include <memory>
#include <string>

// user include files
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/RefToBase.h" 
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "TrackerAlignmentUserCode/VertexAnalysis/interface/VertexAnalysis.h"

#include "TTree.h"


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
VertexAnalysis::VertexAnalysis(const edm::ParameterSet& iConfig) 
  : maxNTracks_(iConfig.getParameter<int>("MaxNTracks")),
    maxNVertices_(iConfig.getParameter<int>("MaxNVertices")),
    maxNGenVertices_(10),
    maxNVerticesSplitTracksRefit_(2) {
  treeName_ = iConfig.getParameter<std::string>("TreeName");
  tracksTag_ = iConfig.getParameter<edm::InputTag>("TrackCollection");
  vertexTag_ = iConfig.getParameter<edm::InputTag>("VertexCollection");
  genParticlesTag_ = iConfig.getParameter<edm::InputTag>("GenParticlesCollection");
}


VertexAnalysis::~VertexAnalysis() {
  delete [] tracksPt_;
  delete [] tracksEta_;
  delete [] tracksPhi_;
  delete [] tracksPx_;
  delete [] tracksPy_;
  delete [] tracksPz_;
  delete [] tracksNChi2_;
  delete [] tracksChi2_;
  delete [] tracksNdof_;
  delete [] tracksQualityIsLoose_;
  delete [] tracksQualityIsTight_;
  delete [] tracksQualityIsHighPurity_;

  delete [] vtxX_;
  delete [] vtxY_;
  delete [] vtxZ_;
  delete [] vtxXErr_;
  delete [] vtxYErr_;
  delete [] vtxZErr_;
  delete [] vtxNChi2_;
  delete [] vtxChi2_;
  delete [] vtxNdof_;
  delete [] vtxNTracks_;

  delete [] genVtxX_;
  delete [] genVtxY_;
  delete [] genVtxZ_;
  delete [] genPDGId_;

  delete [] vtxSplitTracksRefitX_;
  delete [] vtxSplitTracksRefitY_;
  delete [] vtxSplitTracksRefitZ_;
  delete [] vtxSplitTracksRefitXErr_;
  delete [] vtxSplitTracksRefitYErr_;
  delete [] vtxSplitTracksRefitZErr_;
  delete [] vtxSplitTracksRefitNChi2_;
  delete [] vtxSplitTracksRefitChi2_;
  delete [] vtxSplitTracksRefitNdof_;
  delete [] vtxSplitTracksRefitNTracks_;
}


//
// member functions
//

// ------------ method called for each event  ------------
void
VertexAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // --- reset tree variables to default values -----------------------------
  setTreeVariablesToDefault();

  // --- store event information --------------------------------------------
  run_   = iEvent.id().run();
  event_ = iEvent.id().event();
  lumi_  = iEvent.luminosityBlock();


  // --- store track information --------------------------------------------

  // get collection of reconstructed tracks from event
  edm::Handle<reco::TrackCollection> tracksHandle;
  iEvent.getByLabel(tracksTag_,tracksHandle);

  if( tracksHandle.isValid() ) {
    const reco::TrackCollection* tracks = tracksHandle.product();

    // store information of tracks
    for(size_t itk = 0; itk < tracks->size() && itk < maxNTracks_; ++itk) {
      tracksPt_[tracksN_] = tracks->at(itk).pt();
      tracksEta_[tracksN_] = tracks->at(itk).eta();
      tracksPhi_[tracksN_] = tracks->at(itk).phi();
      tracksPx_[tracksN_] = tracks->at(itk).px();
      tracksPy_[tracksN_] = tracks->at(itk).py();
      tracksPz_[tracksN_] = tracks->at(itk).pz();
      tracksNChi2_[tracksN_] = tracks->at(itk).normalizedChi2();
      tracksChi2_[tracksN_] = tracks->at(itk).chi2();
      tracksNdof_[tracksN_] = tracks->at(itk).ndof();
      tracksQualityIsLoose_[tracksN_] = tracks->at(itk).quality(reco::TrackBase::loose);
      tracksQualityIsTight_[tracksN_] = tracks->at(itk).quality(reco::TrackBase::tight);
      tracksQualityIsHighPurity_[tracksN_] = tracks->at(itk).quality(reco::TrackBase::highPurity);
      ++tracksN_;
    }
  }


  // --- store vertex information -------------------------------------------

  // get collection of reconstructed vertices from event
  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByLabel(vertexTag_,vertexHandle);

  if( vertexHandle.isValid() ) {
    const reco::VertexCollection* vertices = vertexHandle.product();

    // store index of leading selected primary vertex
    int pvIdx = -1;

    // store information of selected vertices
    for(size_t ivtx = 0; ivtx < vertices->size() && ivtx < maxNVertices_; ++ivtx) {
      reco::Vertex theVtx = vertices->at(ivtx);
      
      // apply quality cuts
      if( !theVtx.isValid() ) continue;
      if( theVtx.ndof() < 4 ) continue;
      if( std::abs(theVtx.z()) > 24 ) continue;
      if( theVtx.position().Rho() > 2. ) continue;

      if( vtxN_ == 0 ) pvIdx = static_cast<int>(ivtx);
      
      // store vertex information
      vtxX_[vtxN_] = theVtx.x();
      vtxY_[vtxN_] = theVtx.y();
      vtxZ_[vtxN_] = theVtx.z();
      vtxXErr_[vtxN_] = theVtx.xError();
      vtxYErr_[vtxN_] = theVtx.yError();
      vtxZErr_[vtxN_] = theVtx.zError();
      vtxNChi2_[vtxN_] = theVtx.normalizedChi2();
      vtxChi2_[vtxN_] = theVtx.chi2();
      vtxNdof_[vtxN_] = theVtx.ndof();
      vtxNTracks_[vtxN_] = theVtx.nTracks();

      vtxN_++;
    }


    // --- split-track vertex -------------------------------------------------
    // Re-fit leading primary vertex from a randomly split track collection to
    // get resolution in a databased way.

    if( pvIdx > -1 ) {		// use leading good PV

      // Need transient tracks to re-fit vertex (they contain trajectory information)
      // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideTransientTracks
      std::vector<reco::TransientTrack> tracksFromPVSet1;
      std::vector<reco::TransientTrack> tracksFromPVSet2;
      edm::ESHandle<TransientTrackBuilder> transientTrackBuilder;
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transientTrackBuilder);

      // Loop over all tracks that were used to fit the leading good PV
      // and split into two sets of tracks
      unsigned int idx = 0;
      const reco::Vertex* pv = &(vertices->at(pvIdx));
      for(reco::Vertex::trackRef_iterator tv = pv->tracks_begin();
    	  tv != pv->tracks_end(); ++tv, ++idx) {
    	const reco::Track& track = *(tv->get()); 
    	if( idx % 2 == 0 ) {	// split into tracks with odd and event index
    	  tracksFromPVSet1.push_back( transientTrackBuilder->build(track) ); 
	} else {
    	  tracksFromPVSet2.push_back( transientTrackBuilder->build(track) ); 
    	}
      }

      // std::cout << "  " << tracksFrom1stPV1.size() << " new tracks:" << std::endl;
      // for(std::vector<reco::TransientTrack>::const_iterator ti = tracksFrom1stPV1.begin(); ti != tracksFrom1stPV1.end(); ++ti) {
      // 	const reco::Track& track = ti->track();
      // 	std::cout << "  " << track.outerPx() << std::endl;
      // }
      
      // fit the vertex with only the selected tracks
      AdaptiveVertexFitter avf;
      avf.setWeightThreshold(0.1); //weight per track. allow almost every fit, else --> exception
      for(unsigned trackSetIdx = 0; trackSetIdx < maxNVerticesSplitTracksRefit_; ++trackSetIdx) {
	TransientVertex transVtx;
	bool fitOK = true;
	try {
	  if( trackSetIdx == 0 ) transVtx = avf.vertex(tracksFromPVSet1);
	  else                   transVtx = avf.vertex(tracksFromPVSet2);
	} catch(...) {
	  fitOK = false;
	}
	if( fitOK ) {
	  const reco::Vertex pvReFit(transVtx);
	  //std::cout << "  Re-Fitted PV: x = " << pvReFit.x() << std::endl;
	  
	  vtxSplitTracksRefitX_[vtxSplitTracksRefitN_] = pvReFit.x();
	  vtxSplitTracksRefitY_[vtxSplitTracksRefitN_] = pvReFit.y();
	  vtxSplitTracksRefitZ_[vtxSplitTracksRefitN_] = pvReFit.z();
	  vtxSplitTracksRefitXErr_[vtxSplitTracksRefitN_] = pvReFit.xError();
	  vtxSplitTracksRefitYErr_[vtxSplitTracksRefitN_] = pvReFit.yError();
	  vtxSplitTracksRefitZErr_[vtxSplitTracksRefitN_] = pvReFit.zError();
	  vtxSplitTracksRefitNChi2_[vtxSplitTracksRefitN_] = pvReFit.normalizedChi2();
	  vtxSplitTracksRefitChi2_[vtxSplitTracksRefitN_] = pvReFit.chi2();
	  vtxSplitTracksRefitNdof_[vtxSplitTracksRefitN_] = pvReFit.ndof();
	  vtxSplitTracksRefitNTracks_[vtxSplitTracksRefitN_] = pvReFit.nTracks();
	  
	  vtxSplitTracksRefitN_++;
	}
      }
    }

  } // end if vertexHandle is valid


  // store vertices of generator-level status-3 (before showering) particles
  edm::Handle<reco::GenParticleCollection> genParticlesHandle;
  iEvent.getByLabel(genParticlesTag_,genParticlesHandle); 
  if( genParticlesHandle.isValid() ) {
    const reco::GenParticleCollection* genParticles = genParticlesHandle.product();

    for(unsigned int ipart = 0; ipart < genParticles->size() && genVtxN_ < maxNGenVertices_; ipart ++){
      const int status = genParticles->at(ipart).status();
      const int pdgId = genParticles->at(ipart).pdgId();
      if( status == 3 && std::abs(pdgId) != 2212) {
	// status-3 particles: before showering
	// PDG ID 2212: proton. They are simulated at (0,0,0)
	genVtxX_[genVtxN_] = genParticles->at(ipart).vx();
	genVtxY_[genVtxN_] = genParticles->at(ipart).vy();
	genVtxZ_[genVtxN_] = genParticles->at(ipart).vz();
	genPDGId_[genVtxN_] = pdgId;
	genVtxN_++;
      }
    }
  }


  tree_->Fill();
}


void
VertexAnalysis::setTreeVariablesToDefault() {
  run_   = 0;
  event_ = 0;
  lumi_  = 0;

  tracksN_ = 0;
  for(unsigned int i = 0; i < maxNTracks_; ++i) {
    tracksPt_[i] = -99999.;
    tracksEta_[i] = -99999.;
    tracksPhi_[i] = -99999.;
    tracksPx_[i] = -99999.;
    tracksPy_[i] = -99999.;
    tracksPz_[i] = -99999.;
    tracksNChi2_[i] = -99999.;
    tracksChi2_[i] = -99999.;
    tracksNdof_[i] = -99999.;
    tracksQualityIsLoose_[i] = false;
    tracksQualityIsTight_[i] = false;
    tracksQualityIsHighPurity_[i] = false;
  }

  vtxN_ = 0;
  for(unsigned int i = 0; i < maxNVertices_; ++i) {
    vtxX_[i] = -99999.;
    vtxY_[i] = -99999.;
    vtxZ_[i] = -99999.;
    vtxXErr_[i] = -99999.;
    vtxYErr_[i] = -99999.;
    vtxZErr_[i] = -99999.;
    vtxNChi2_[i] = -99999.;
    vtxChi2_[i] = -99999.;
    vtxNdof_[i] = -99999.;
    vtxNTracks_[i] = 0;
  }

  genVtxN_ = 0;
  for(unsigned int i = 0; i < maxNGenVertices_; ++i) {
    genVtxX_[i] = -99999.;
    genVtxY_[i] = -99999.;
    genVtxZ_[i] = -99999.;
    genPDGId_[i] = -99999;
  }

  vtxSplitTracksRefitN_ = 0;
  for(unsigned int i = 0; i < maxNVerticesSplitTracksRefit_; ++i) {
    vtxSplitTracksRefitX_[i] = -99999.;
    vtxSplitTracksRefitY_[i] = -99999.;
    vtxSplitTracksRefitZ_[i] = -99999.;
    vtxSplitTracksRefitXErr_[i] = -99999.;
    vtxSplitTracksRefitYErr_[i] = -99999.;
    vtxSplitTracksRefitZErr_[i] = -99999.;
    vtxSplitTracksRefitNChi2_[i] = -99999.;
    vtxSplitTracksRefitChi2_[i] = -99999.;
    vtxSplitTracksRefitNdof_[i] = -99999.;
    vtxSplitTracksRefitNTracks_[i] = 0;
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
VertexAnalysis::beginJob() {
  tracksPt_ = new double[maxNTracks_];
  tracksEta_ = new double[maxNTracks_];
  tracksPhi_ = new double[maxNTracks_];
  tracksPx_ = new double[maxNTracks_];
  tracksPy_ = new double[maxNTracks_];
  tracksPz_ = new double[maxNTracks_];
  tracksNChi2_ = new double[maxNTracks_];
  tracksChi2_ = new double[maxNTracks_];
  tracksNdof_ = new double[maxNTracks_];
  tracksQualityIsLoose_ = new bool[maxNTracks_];
  tracksQualityIsTight_ = new bool[maxNTracks_];
  tracksQualityIsHighPurity_ = new bool[maxNTracks_];

  vtxX_ = new double[maxNVertices_];
  vtxY_ = new double[maxNVertices_];
  vtxZ_ = new double[maxNVertices_];
  vtxXErr_ = new double[maxNVertices_];
  vtxYErr_ = new double[maxNVertices_];
  vtxZErr_ = new double[maxNVertices_];
  vtxNChi2_ = new double[maxNVertices_];
  vtxChi2_ = new double[maxNVertices_];
  vtxNdof_ = new double[maxNVertices_];
  vtxNTracks_ = new unsigned int[maxNVertices_];

  genVtxX_ = new double[maxNGenVertices_];
  genVtxY_ = new double[maxNGenVertices_];
  genVtxZ_ = new double[maxNGenVertices_];
  genPDGId_ = new int[maxNGenVertices_];

  vtxSplitTracksRefitX_ = new double[maxNVerticesSplitTracksRefit_];
  vtxSplitTracksRefitY_ = new double[maxNVerticesSplitTracksRefit_];
  vtxSplitTracksRefitZ_ = new double[maxNVerticesSplitTracksRefit_];
  vtxSplitTracksRefitXErr_ = new double[maxNVerticesSplitTracksRefit_];
  vtxSplitTracksRefitYErr_ = new double[maxNVerticesSplitTracksRefit_];
  vtxSplitTracksRefitZErr_ = new double[maxNVerticesSplitTracksRefit_];
  vtxSplitTracksRefitNChi2_ = new double[maxNVerticesSplitTracksRefit_];
  vtxSplitTracksRefitChi2_ = new double[maxNVerticesSplitTracksRefit_];
  vtxSplitTracksRefitNdof_ = new double[maxNVerticesSplitTracksRefit_];
  vtxSplitTracksRefitNTracks_ = new unsigned int[maxNVerticesSplitTracksRefit_];

  setTreeVariablesToDefault();

  // set up output file and TTree
  edm::Service<TFileService> fs;
  if( !fs ) {
    throw edm::Exception(edm::errors::Configuration, "TFile Service is not registered in cfg file");
  }
  tree_ = fs->make<TTree>(treeName_.c_str(),treeName_.c_str());
  tree_->SetAutoSave(10000000000);
  tree_->SetAutoFlush(1000000);

  tree_->Branch("Run",&run_,"Run/i");
  tree_->Branch("Lumi",&lumi_,"Lumi/i");
  tree_->Branch("Event",&event_,"Event/i");

  tree_->Branch("TracksN",&tracksN_,"TracksN/i");
  tree_->Branch("TracksPt",tracksPt_,"TracksPt[TracksN]/D");
  tree_->Branch("TracksEta",tracksEta_,"TracksEta[TracksN]/D");
  tree_->Branch("TracksPhi",tracksPhi_,"TracksPhi[TracksN]/D");
  tree_->Branch("TracksPx",tracksPx_,"TracksPx[TracksN]/D");
  tree_->Branch("TracksPy",tracksPy_,"TracksPy[TracksN]/D");
  tree_->Branch("TracksPz",tracksPz_,"TracksPz[TracksN]/D");
  tree_->Branch("TracksNChi2",tracksNChi2_,"TracksNChi2[TracksN]/D");
  tree_->Branch("TracksChi2",tracksChi2_,"TracksChi2[TracksN]/D");
  tree_->Branch("TracksNdof",tracksNdof_,"TracksNdof[TracksN]/D");
  tree_->Branch("TracksQualityIsLoose",tracksQualityIsLoose_,"TracksQualityIsLoose[TracksN]/O");
  tree_->Branch("TracksQualityIsTight",tracksQualityIsTight_,"TracksQualityIsTight[TracksN]/O");
  tree_->Branch("TracksQualityIsHighPurity",tracksQualityIsHighPurity_,"TracksQualityIsHighPurity[TracksN]/O");

  tree_->Branch("VtxN",&vtxN_,"VtxN/i");
  tree_->Branch("VtxX",vtxX_,"VtxX[VtxN]/D");
  tree_->Branch("VtxY",vtxY_,"VtxY[VtxN]/D");
  tree_->Branch("VtxZ",vtxZ_,"VtxZ[VtxN]/D");
  tree_->Branch("VtxXErr",vtxXErr_,"VtxXErr[VtxN]/D");
  tree_->Branch("VtxYErr",vtxYErr_,"VtxYErr[VtxN]/D");
  tree_->Branch("VtxZErr",vtxZErr_,"VtxZErr[VtxN]/D");
  tree_->Branch("VtxNChi2",vtxNChi2_,"VtxNChi2[VtxN]/D");
  tree_->Branch("VtxChi2",vtxChi2_,"VtxChi2[VtxN]/D");
  tree_->Branch("VtxNdof",vtxNdof_,"VtxNdof[VtxN]/D");
  tree_->Branch("VtxNTracks",vtxNTracks_,"VtxNTracks[VtxN]/i");

  tree_->Branch("GenVtxN",&genVtxN_,"GenVtxN/i");
  tree_->Branch("GenVtxX",genVtxX_,"GenVtxX[GenVtxN]/D");
  tree_->Branch("GenVtxY",genVtxY_,"GenVtxY[GenVtxN]/D");
  tree_->Branch("GenVtxZ",genVtxZ_,"GenVtxZ[GenVtxN]/D");
  tree_->Branch("GenPDGId",genPDGId_,"GenPDGId[GenVtxN]/I");

  tree_->Branch("VtxSplitTracksRefitN",&vtxSplitTracksRefitN_,"VtxSplitTracksRefitN/i");
  tree_->Branch("VtxSplitTracksRefitX",vtxSplitTracksRefitX_,"VtxSplitTracksRefitX[VtxSplitTracksRefitN]/D");
  tree_->Branch("VtxSplitTracksRefitY",vtxSplitTracksRefitY_,"VtxSplitTracksRefitY[VtxSplitTracksRefitN]/D");
  tree_->Branch("VtxSplitTracksRefitZ",vtxSplitTracksRefitZ_,"VtxSplitTracksRefitZ[VtxSplitTracksRefitN]/D");
  tree_->Branch("VtxSplitTracksRefitXErr",vtxSplitTracksRefitXErr_,"VtxSplitTracksRefitXErr[VtxSplitTracksRefitN]/D");
  tree_->Branch("VtxSplitTracksRefitYErr",vtxSplitTracksRefitYErr_,"VtxSplitTracksRefitYErr[VtxSplitTracksRefitN]/D");
  tree_->Branch("VtxSplitTracksRefitZErr",vtxSplitTracksRefitZErr_,"VtxSplitTracksRefitZErr[VtxSplitTracksRefitN]/D");
  tree_->Branch("VtxSplitTracksRefitNChi2",vtxSplitTracksRefitNChi2_,"VtxSplitTracksRefitNChi2[VtxSplitTracksRefitN]/D");
  tree_->Branch("VtxSplitTracksRefitChi2",vtxSplitTracksRefitChi2_,"VtxSplitTracksRefitChi2[VtxSplitTracksRefitN]/D");
  tree_->Branch("VtxSplitTracksRefitNdof",vtxSplitTracksRefitNdof_,"VtxSplitTracksRefitNdof[VtxSplitTracksRefitN]/D");
  tree_->Branch("VtxSplitTracksRefitNTracks",vtxSplitTracksRefitNTracks_,"VtxSplitTracksRefitNTracks[VtxSplitTracksRefitN]/i");
}


// ------------ method called once each job just after ending the event loop  ------------
void 
VertexAnalysis::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
VertexAnalysis::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
VertexAnalysis::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
VertexAnalysis::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
VertexAnalysis::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
VertexAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(VertexAnalysis);
