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
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "TTree.h"

//
// class declaration
//

class VertexAnalysis : public edm::EDAnalyzer {
public:
  explicit VertexAnalysis(const edm::ParameterSet&);
  ~VertexAnalysis();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  void setTreeVariablesToDefault();

  // ----------member data ---------------------------

  // constants
  const unsigned int maxNTracks_;
  const unsigned int maxNVertices_;
  const unsigned int maxNGenVertices_;
  const unsigned int maxNVerticesSplitTracksRefit_;
  
  // input tags
  edm::InputTag tracksTag_;
  edm::InputTag vertexTag_;
  edm::InputTag genParticlesTag_;

  // tree variables
  TTree* tree_;
  std::string treeName_;

  unsigned int run_;
  unsigned int event_;
  unsigned int lumi_;

  // tracks (from event)
  unsigned int tracksN_;
  double* tracksPt_;
  double* tracksEta_;
  double* tracksPhi_;
  double* tracksPx_;
  double* tracksPy_;
  double* tracksPz_;
  double* tracksNChi2_;
  double* tracksChi2_;
  double* tracksNdof_;
  bool* tracksQualityIsLoose_;
  bool* tracksQualityIsTight_;
  bool* tracksQualityIsHighPurity_;
  
  // primary vertices (from event)
  unsigned int vtxN_;
  double* vtxX_;
  double* vtxY_;
  double* vtxZ_;
  double* vtxXErr_;
  double* vtxYErr_;
  double* vtxZErr_;
  double* vtxNChi2_;
  double* vtxChi2_;
  double* vtxNdof_;
  unsigned int* vtxNTracks_;

  unsigned int genVtxN_;
  double* genVtxX_;
  double* genVtxY_;
  double* genVtxZ_;
  int* genPDGId_;

  // re-fitted primary vertices from split track sample
  // of leading-PV associated tracks
  unsigned int vtxSplitTracksRefitN_;
  double* vtxSplitTracksRefitX_;
  double* vtxSplitTracksRefitY_;
  double* vtxSplitTracksRefitZ_;
  double* vtxSplitTracksRefitXErr_;
  double* vtxSplitTracksRefitYErr_;
  double* vtxSplitTracksRefitZErr_;
  double* vtxSplitTracksRefitNChi2_;
  double* vtxSplitTracksRefitChi2_;
  double* vtxSplitTracksRefitNdof_;
  unsigned int* vtxSplitTracksRefitNTracks_;
};
