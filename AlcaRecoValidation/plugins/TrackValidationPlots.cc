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
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/RefToBase.h" 
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TH1.h"
#include "TH1D.h"



class TrackValidationPlots : public edm::EDAnalyzer {
public:
  explicit TrackValidationPlots(const edm::ParameterSet&);
  ~TrackValidationPlots() {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  // ----------member data ---------------------------
  
  // input tags
  edm::InputTag tracksTag_;
  edm::InputTag vertexTag_;

  TH1* hTrackN_;
  TH1* hTrackPt_;
  TH1* hTrackEta_;
  TH1* hTrackPhi_;
  TH1* hTrackChi2Ndof_;
  TH1* hLeadingTrackPt_;
  TH1* hLeadingTrackEta_;
  TH1* hLeadingTrackPhi_;
  TH1* hLeadingTrackChi2Ndof_;
};



//
// constructors and destructor
//
TrackValidationPlots::TrackValidationPlots(const edm::ParameterSet& iConfig) {
  tracksTag_ = iConfig.getParameter<edm::InputTag>("TrackCollection");
  vertexTag_ = iConfig.getParameter<edm::InputTag>("VertexCollection");
}


//
// member functions
//

// ------------ method called for each event  ------------
void
TrackValidationPlots::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // --- store track information --------------------------------------------

  // get collection of reconstructed tracks from event
  edm::Handle<reco::TrackCollection> tracksHandle;
  iEvent.getByLabel(tracksTag_,tracksHandle);

  if( tracksHandle.isValid() ) {
    const reco::TrackCollection* tracks = tracksHandle.product();

    // store information of tracks
    size_t itk = 0;
    for(; itk < tracks->size(); ++itk) {
      if( itk == 0 ) {
	hLeadingTrackPt_->Fill(tracks->at(itk).pt());
	hLeadingTrackEta_->Fill(tracks->at(itk).eta());
	hLeadingTrackPhi_->Fill(tracks->at(itk).phi());
	hLeadingTrackPhi_->Fill(tracks->at(itk).phi());
	hLeadingTrackChi2Ndof_->Fill(tracks->at(itk).normalizedChi2());
      }
      hTrackPt_->Fill(tracks->at(itk).pt());
      hTrackEta_->Fill(tracks->at(itk).eta());
      hTrackPhi_->Fill(tracks->at(itk).phi());
      hTrackPhi_->Fill(tracks->at(itk).phi());
      hTrackChi2Ndof_->Fill(tracks->at(itk).normalizedChi2());
    }
    hTrackN_->Fill(itk);
  }


  // --- store vertex information -------------------------------------------

  // // get collection of reconstructed vertices from event
  // edm::Handle<reco::VertexCollection> vertexHandle;
  // iEvent.getByLabel(vertexTag_,vertexHandle);

  // if( vertexHandle.isValid() ) {
  //   const reco::VertexCollection* vertices = vertexHandle.product();
  // }

}


// ------------ method called once each job just before starting event loop  ------------
void 
TrackValidationPlots::beginJob() {
  // set up output file and TTree
  edm::Service<TFileService> fs;
  if( !fs ) {
    throw edm::Exception(edm::errors::Configuration, "TFile Service is not registered in cfg file");
  }
  hTrackN_ = fs->make<TH1D>("TrackN",";N(tracks);N(events)",100,0,100);
  hTrackPt_ = fs->make<TH1D>("TrackPt",";p_{T}(track) [GeV];N(tracks)",100,0,200);
  hTrackEta_ = fs->make<TH1D>("TrackEta",";#eta(track);N(tracks)",60,-3,3);
  hTrackPhi_ = fs->make<TH1D>("TrackPhi",";#phi(track);N(tracks)",70,-3.5,3.5);
  hTrackChi2Ndof_ = fs->make<TH1D>("TrackChi2Ndof",";#chi^{2}/ndof(track);N(tracks)",100,0.,10.);
  hLeadingTrackPt_ = fs->make<TH1D>("LeadingTrackPt",";p_{T}(leading track) [GeV];N(events)",100,0,200);
  hLeadingTrackEta_ = fs->make<TH1D>("LeadingTrackEta",";#eta(leading track);N(events)",60,-3,3);
  hLeadingTrackPhi_ = fs->make<TH1D>("LeadingTrackPhi",";#phi(leading track);N(events)",70,-3.5,3.5);
  hLeadingTrackChi2Ndof_ = fs->make<TH1D>("LeadingTrackChi2Ndof",";#chi^{2}/ndof(track);N(events)",100,0.,10.);
}


// ------------ method called once each job just after ending the event loop  ------------
void 
TrackValidationPlots::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
TrackValidationPlots::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TrackValidationPlots::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TrackValidationPlots::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TrackValidationPlots::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrackValidationPlots::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackValidationPlots);
