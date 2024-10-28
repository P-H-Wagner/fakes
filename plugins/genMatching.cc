#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
// for vertex fitting (both global and sequential)
#include "TrackingTools/Records/interface/TransientTrackRecord.h" 
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h" 
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h" 
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h" 
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h" 
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h" 
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h" 
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h" 
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h" 
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicConstraint.h" 
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "RecoVertex/KinematicFitPrimitives/interface/Matrices.h" 
#include "DataFormats/GeometryVector/interface/GlobalPoint.h" 

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // for the tracks!
#include "FWCore/Framework/interface/MakerMacros.h"
#include <limits>
#include <algorithm>
#include <cmath>
#include "KinVtxFitter.h" //--> not needed now 
#include "helper.h" // helper functions
// include "PxPyPzMVector.h" // to new :(
#include "TLorentzVector.h" // use this instead 
#include "TVector3.h" // for boost vector
// for gen matching
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h" 

// B field
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"
// for 3DPoint
#include "DataFormats/Math/interface/Point3D.h"
// for the cov matrix correction
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h" 


// counters
int nEventsGen = 0;
int nMuonsGen  = 0;
int nTracksGen = 0;
int nRecoCandidates = 0;
int nGenMatched = 0;

int muSel1CounterGen = 0;
int muSel2CounterGen = 0;
int k1Sel1CounterGen = 0;
int k1Sel2CounterGen = 0;
int k2Sel1CounterGen = 0;
int k2Sel2CounterGen = 0;
int piSel1CounterGen = 0;
int piSel2CounterGen = 0;
int nKKPiMuGen = 0;
int nFoundPhi  = 0;
int nFoundDs   = 0;
int nFoundB    = 0;
int nBMassCut   = 0;

class genMatching : public edm::global::EDProducer<> {

public:

  //define collections which dont exist by default  
  typedef std::vector<reco::GenParticle> GenParticleCollection;
  typedef std::vector<pat::PackedGenParticle> PackedGenParticleCollection;
  //constructor
  explicit genMatching(const edm::ParameterSet&);
  //destructor
  ~genMatching() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  virtual void endJob() override;

  int  getPVIdx(const reco::VertexCollection*,const reco::TransientTrack&) const;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}


  //must be constant as it takes constant arguments!! otherwise compiler rises an error
  //because it thinks it may modify the input!!
  reco::TransientTrack getTransientTrack(const reco::Track track) const {    
      reco::TransientTrack transientTrack(track, paramField);

      return transientTrack;
    }


private:
  
  //Bfield
  OAEParametrizedMagneticField *paramField = new OAEParametrizedMagneticField("3_8T");

  //cuts 
  const StringCutObjectSelector<pat::PackedCandidate> hadSelection_; // cut on hadrons
  //const StringCutObjectSelector<pat::PackedGenParticle> hadSelectionGen_; // cut on gen hadrons
  const StringCutObjectSelector<reco::GenParticle> hadSelectionGen_; // cut on gen hadrons for test with pruned only!

  const double minMuPt_;
  const double maxMuEta_;
  const double maxdRHadMuon_;
  const double mindRHadMuon_;
  const double maxdzDiffHadMuon_; 
  const double phiMassAllowance_;
  const double dsMassAllowance_;
  const double drMatchGen_;
  const double ptMatchGen_;
  const double maxBsMass_;
  const double piMass_;
  const double kMass_;
  const double phiMass_;
  const double dsMass_;
  const double dsStarMass_;
  const double muMass_;
  const double tauMass_;
  const double bsMass_;
  const double isoCone_;
  //tokens to access data later
  //edm::Input tag can not be directly initialized inside the construcor! Why did it work fro Trigger.cc??
  //anyway ... 

  //for the muons

  // vertices
  const edm::InputTag primaryVtxTag;
  const edm::EDGetTokenT<reco::VertexCollection> primaryVtx_;

  const edm::InputTag bsTag;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> bs_;

  //gen for gen-matching
  const edm::InputTag prunedGenTag; //pruned is a compressed packed format
  const edm::EDGetTokenT<reco::GenParticleCollection> prunedGen_;

  const edm::InputTag packedGenTag; //packed contains much more info->most likely not needed!
  const edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGen_;

};

//define the constructor
genMatching::genMatching(const edm::ParameterSet& iConfig):
    // f.e. hadSelection_ = cfg.getPatameter...

    hadSelection_(iConfig.getParameter<std::string>("hadSelection")),
    hadSelectionGen_(iConfig.getParameter<std::string>("hadSelectionGen")),
    minMuPt_(iConfig.getParameter<double>("minMuPt")),
    maxMuEta_(iConfig.getParameter<double>("maxMuEta")),
    maxdRHadMuon_(iConfig.getParameter<double>("maxdRHadMuon")),
    mindRHadMuon_(iConfig.getParameter<double>("mindRHadMuon")),
    maxdzDiffHadMuon_(iConfig.getParameter<double>("maxdzDiffHadMuon")),
    phiMassAllowance_(iConfig.getParameter<double>("phiMassAllowance")),
    dsMassAllowance_(iConfig.getParameter<double>("dsMassAllowance")),
    drMatchGen_(iConfig.getParameter<double>("drMatchGen")),
    ptMatchGen_(iConfig.getParameter<double>("ptMatchGen")),
    maxBsMass_(iConfig.getParameter<double>("maxBsMass")),

    piMass_(iConfig.getParameter<double>("piMass")),
    kMass_(iConfig.getParameter<double>("kMass")),
    phiMass_(iConfig.getParameter<double>("phiMass")),
    dsMass_(iConfig.getParameter<double>("dsMass")),
    dsStarMass_(iConfig.getParameter<double>("dsStarMass")),
    muMass_(iConfig.getParameter<double>("muMass")),
    tauMass_(iConfig.getParameter<double>("tauMass")),
    bsMass_(iConfig.getParameter<double>("bsMass")),
    isoCone_(iConfig.getParameter<double>("isoCone")),

    primaryVtxTag(iConfig.getParameter<edm::InputTag>("pvCand")),
    primaryVtx_(consumes<reco::VertexCollection>(primaryVtxTag)),

    bsTag(iConfig.getParameter<edm::InputTag>("bs")),
    bs_(consumes<pat::CompositeCandidateCollection>(bsTag)), 

    prunedGenTag(iConfig.getParameter<edm::InputTag>("prunedCand")),
    prunedGen_(consumes<reco::GenParticleCollection>(prunedGenTag)),
    packedGenTag(iConfig.getParameter<edm::InputTag>("packedCand")),
    packedGen_(consumes<pat::PackedGenParticleCollection>(packedGenTag)){
       // output collection
       produces<pat::CompositeCandidateCollection>("gen");
       //produces<pat::CompositeCandidateCollection>("gen");
       //produces<TransientTrackCollection>("kkTransientTracks");
    }

//check const keywords 

// this starts the event loop
void genMatching::produce(edm::StreamID, edm::Event &iEvent, const edm::EventSetup &iSetup) const {

  //std::cout << "I entered here" << std::endl;
  //input
  edm::Handle<reco::VertexCollection> primaryVtx;
  iEvent.getByToken(primaryVtx_,primaryVtx);

  edm::Handle<pat::CompositeCandidateCollection> bsColl;
  iEvent.getByToken(bs_,bsColl);
 
  edm::Handle<reco::GenParticleCollection> prunedGen;
  iEvent.getByToken(prunedGen_,prunedGen);

  edm::Handle<pat::PackedGenParticleCollection> packedGen;
  iEvent.getByToken(packedGen_,packedGen);

  // to save 
  std::unique_ptr<pat::CompositeCandidateCollection> ret_value(new pat::CompositeCandidateCollection());
  //std::unique_ptr<pat::CompositeCandidateCollection> ret_value_gen(new pat::CompositeCandidateCollection());
  //std::unique_ptr<TransientTrackCollection> kkpi_ttrack(new TransientTrackCollection);

  nEventsGen++;

  //////////////////////////////////////////////////////
  // Match the bsCand (and all its final states) with //
  // a gen particle. We loop over all candidates      //
  // per event!                                       //
  //////////////////////////////////////////////////////

  pat::CompositeCandidate gen; 
  for(size_t bsIdx = 0; bsIdx < bsColl->size(); ++bsIdx){
  

    nRecoCandidates++;
 
    pat::CompositeCandidate gen; 
    edm::Ptr<pat::CompositeCandidate> bsPtr(bsColl, bsIdx);

    //get kinematic info of final states
    auto muBs = bsPtr->userCand("mu");
    auto k1Bs = bsPtr->userCand("k1");
    auto k2Bs = bsPtr->userCand("k2");
    auto piBs = bsPtr->userCand("pi");
 

    // reco cand
    //std::cout << "mu pt: " << muBs->pt() << std::endl;
    //std::cout << "k1 pt: " << k1Bs->pt() << std::endl;


    int sigId = -9999;
    int bId = 0;
    int genMatchSuccess = 0;

    //count the number of gen matches we find, ideally only 1
    int nGenMatches = 0;

    ////////////////////////////////////////////////////
    // find the gen-matched muon                      //
    ////////////////////////////////////////////////////
    
    //std::cout << "reco mu pt: " << muBs->pt() << std::endl;
    //std::cout << "reco mu eta: " << muBs->eta() << std::endl;
    //std::cout << "reco mu phi: " << muBs->phi() << std::endl;

    /*
    for(size_t muIdxGen = 0; muIdxGen < prunedGen->size(); ++muIdxGen){

      nMuonsGen++;

      //define a pointer to the gen muon    
      edm::Ptr<reco::GenParticle> muPtrGen(prunedGen, muIdxGen);

      //if ( fabs(muPtrGen->pdgId()) == 13 ){
      //std::cout << "mu pt: " << muPtrGen->pt() << std::endl;
      //std::cout << "mu pdg id: " << muPtrGen->pdgId() << std::endl;
      //std::cout << "mother id: " << muPtrGen->mother(0)->pdgId() << std::endl;
      //std::cout << "charge: " << muPtrGen->charge() << std::endl;
      //}

     if (reco::deltaR(*muBs,*muPtrGen) < 0.1){

     //std::cout << "Gen particle in cone with ID: " << muPtrGen->pdgId() << std::endl;
     //std::cout << "Gen particle in cone with pt: " << muPtrGen->pt() << std::endl;
     //std::cout << "Gen particle in cone with eta: " << muPtrGen->eta() << std::endl;
     //std::cout << "Gen particle in cone with phi: " << muPtrGen->phi() << std::endl;
     //std::cout << "Gen particle in cone with charge: " << muPtrGen->charge() << std::endl;

     }

    }
    */

    for(size_t muIdxGen = 0; muIdxGen < prunedGen->size(); ++muIdxGen){

      nMuonsGen++;

      //define a pointer to the gen muon    
      edm::Ptr<reco::GenParticle> muPtrGen(prunedGen, muIdxGen);

      //std::cout << "mu pt: " << muPtrGen->pt() << std::endl;
      //std::cout << "pdg id: " << muPtrGen->pdgId() << std::endl;
      //std::cout << "mother id: " << muPtrGen->mother(0)->pdgId() << std::endl;
      //std::cout << "charge: " << muPtrGen->charge() << std::endl;

      //select only useful gen muons -> check this selection!
      if( (fabs(muPtrGen->pdgId()) != 13) || muPtrGen->pt() < minMuPt_ || fabs(muPtrGen->eta()) > maxMuEta_ || (muBs->charge() * muPtrGen->charge() < 0)) continue; 

      muSel1CounterGen++;

      //now check the dR of the reco muon wrt to the gen Muon 
      float drMuonMatch = reco::deltaR(*muBs,*muPtrGen);
      float ptMuonMatch = fabs(muBs->pt() - muPtrGen->pt());

      //TODO:define as variable
      if (drMuonMatch > drMatchGen_) continue;
      if (ptMuonMatch > ptMatchGen_) continue;
      muSel2CounterGen++;


      //std::cout << "yes!!" <<  std::endl;
      //std::cout << muPtrGen->pdgId() <<  std::endl;
      //std::cout << muPtrGen->pt() << std::endl;
      //std::cout << "found a gen matched muon" << std::endl;

      ////////////////////////////////////////////////
      // find gen matched k1                        //
      ////////////////////////////////////////////////

      for(size_t k1IdxGen = 0; k1IdxGen < prunedGen->size(); ++k1IdxGen){
           
        nTracksGen++;
        //define a pointer to the gen kaon    
        edm::Ptr<reco::GenParticle> k1PtrGen(prunedGen, k1IdxGen);

        //select only useful kaons -> check this selection!
        if((fabs(k1PtrGen->pdgId()) != 321) || !hadSelectionGen_(*k1PtrGen) || (k1Bs->charge() * k1PtrGen->charge() < 0)) continue; 
        k1Sel1CounterGen++;

        //std::cout << k1PtrGen->pt() << std::endl;
        //std::cout << "passed selection for k1!" << std::endl;
        //now check the dR of the reco muon wrt to the gen Muon 
        float drK1Match = reco::deltaR(*k1Bs,*k1PtrGen);

        //TODO:define as variable
        if(drK1Match > drMatchGen_) continue;
        k1Sel2CounterGen++;
        //std::cout << "found a gen matched k1!" << std::endl;

        ////////////////////////////////////////////////
        // find gen matched k2                        //
        ////////////////////////////////////////////////
     
        for(size_t k2IdxGen = 0; k2IdxGen < prunedGen->size(); ++k2IdxGen){
     
             //avoid picking the same gen particle as for k1
             if(k2IdxGen == k1IdxGen) continue; 

             //define a pointer to the gen kaon    
             edm::Ptr<reco::GenParticle> k2PtrGen(prunedGen, k2IdxGen);
     
             //select only useful kaons -> check this selection!
             if((fabs(k2PtrGen->pdgId()) != 321) || !hadSelectionGen_(*k2PtrGen) || (k2Bs->charge() * k2PtrGen->charge() < 0 )) continue; 
             //std::cout << "reco  pt " << k2Bs->pt() << "reco eta " << k2Bs->eta() << std::endl;
             //std::cout << "gen pt " << k2PtrGen->pt() << "gen eta " << k2PtrGen->eta() << std::endl;
             k2Sel1CounterGen++;
             //std::cout << "passed selection for k2" << std::endl;
             //now check the dR of the reco muon wrt to the gen Muon 
             float drK2Match = reco::deltaR(*k2Bs,*k2PtrGen);
     
             //TODO:define as variable
             if(drK2Match > drMatchGen_ ) continue;
             k2Sel2CounterGen++;

             //std::cout << "found a gen matched k2!\n" << std::endl;

             ////////////////////////////////////////////////
             // find gen matched pion                      //
             ////////////////////////////////////////////////

             for(size_t piIdxGen = 0; piIdxGen < prunedGen->size(); ++piIdxGen){
      
               //avoid picking the same gen particle as for k1 or k2
               if((piIdxGen == k1IdxGen) || (piIdxGen == k2IdxGen)) continue; 
                       
               //define a pointer to the gen kaon    
               edm::Ptr<reco::GenParticle> piPtrGen(prunedGen, piIdxGen);
     
               //select only useful kaons -> check this selection!
               if((fabs(piPtrGen->pdgId()) != 211) || !hadSelectionGen_(*piPtrGen) || (piBs->charge() * piPtrGen->charge() < 0 )) continue; 
               piSel1CounterGen++;
     
               //std::cout << "gen pt " << piPtrGen->pt() << "gen eta " << piPtrGen->eta() << std::endl;
               //now check the dR of the reco muon wrt to the gen Muon 
               float drPiMatch = reco::deltaR(*piBs,*piPtrGen);
               //std::cout << drPiMatch << std::endl; 
               //TODO:define as variable
               if(drPiMatch > drMatchGen_) continue;
               //std::cout << "found a gen matched pion!" << std::endl;
               piSel2CounterGen++;


               //////////////////////////////////////////////////
               // Find resonances at gen level                 //
               //////////////////////////////////////////////////
               nKKPiMuGen++;     
 
               //Should we pick the best gen match (in terms of dR) only? -> No, like this is better 

               const reco::Candidate* k1Reco = k1PtrGen.get(); 
               const reco::Candidate* k2Reco = k2PtrGen.get(); 
               const reco::Candidate* piReco = piPtrGen.get(); 
               const reco::Candidate* muReco = muPtrGen.get(); 

               // searching for phi resonance 
               auto phiFromK1 = getAncestor(k1Reco,333);
               auto phiFromK2 = getAncestor(k2Reco,333);
               if( (phiFromK1 != phiFromK2) || (phiFromK1 == nullptr) || (phiFromK2 == nullptr)) continue; 
               nFoundPhi++;               
     
               // searching for ds resonance 
               auto dsFromPhi = getAncestor(phiFromK1,431);
               auto dsFromPi  = getAncestor(piReco,431);
               if( (dsFromPhi != dsFromPi) || (dsFromPhi == nullptr) || (dsFromPi == nullptr)) continue; 
               nFoundDs++;               

               // we dont know what b mother we have
               int bMotherId = 0;

               // first search for b baryon (isAncestor also checks neg. Ids)
               for(int bIdx = 5000; bIdx < 6000; bIdx++){
                 if (isAncestor(dsFromPhi, bIdx)){
                   bMotherId = bIdx;
                   break;
                 }
               }

               // Then search for possible B meson coming from B-baryon
               for(int bIdx = 500; bIdx < 600; bIdx++){
                 if (isAncestor(dsFromPhi, bIdx)){
                   bMotherId = bIdx;
                   break;
                 }
               }
              
               if (bMotherId == 0) break; // no b mother found

               // Even if the mu is not required to come brom the b mother directly (would be signal case)
               // if it comes from another D meson (double charm background case), we still want
               // that this D meson is coming from the b mother. So the muon should share
               // the same ancestor as the Ds.
               auto bsFromDs = getAncestor(dsFromPhi,bMotherId);
               auto bsFromMu = getAncestor(muReco,   bMotherId);

               if( (bsFromDs != bsFromMu) || (bsFromDs == nullptr) || (bsFromMu == nullptr)) continue; 
   
               nFoundB++;
               
               if (bsFromMu->mass() > maxBsMass_) continue;
               nBMassCut++;

               //remove oscillations
               auto bsFromMuWOOsc = removeOscillations(bsFromMu);

               nGenMatches++;
               genMatchSuccess = 1;
               nGenMatched++;

               gen.addUserFloat("mu_gen_px"      ,muPtrGen->px());
               gen.addUserFloat("mu_gen_py"      ,muPtrGen->py());
               gen.addUserFloat("mu_gen_pz"      ,muPtrGen->pz());
               gen.addUserFloat("mu_gen_pt"      ,muPtrGen->pt());
               gen.addUserFloat("mu_gen_eta"     ,muPtrGen->eta());
               gen.addUserFloat("mu_gen_phi"     ,muPtrGen->phi());
               gen.addUserFloat("mu_gen_m"    ,muPtrGen->mass());
               gen.addUserFloat("mu_gen_charge"  ,muPtrGen->charge());
               gen.addUserInt(  "mu_gen_pdgid"   ,muPtrGen->pdgId());
     
               gen.addUserFloat("k1_gen_px"      ,k1PtrGen->px());
               gen.addUserFloat("k1_gen_py"      ,k1PtrGen->py());
               gen.addUserFloat("k1_gen_pz"      ,k1PtrGen->pz());
               gen.addUserFloat("k1_gen_pt"      ,k1PtrGen->pt());
               gen.addUserFloat("k1_gen_eta"     ,k1PtrGen->eta());
               gen.addUserFloat("k1_gen_phi"     ,k1PtrGen->phi());
               gen.addUserFloat("k1_gen_m"    ,k1PtrGen->mass());
               gen.addUserFloat("k1_gen_charge"  ,k1PtrGen->charge());
               gen.addUserInt(  "k1_gen_pdgid"   ,k1PtrGen->pdgId());
     
               gen.addUserFloat("k2_gen_px"      ,k2PtrGen->px());
               gen.addUserFloat("k2_gen_py"      ,k2PtrGen->py());
               gen.addUserFloat("k2_gen_pz"      ,k2PtrGen->pz());
               gen.addUserFloat("k2_gen_pt"      ,k2PtrGen->pt());
               gen.addUserFloat("k2_gen_eta"     ,k2PtrGen->eta());
               gen.addUserFloat("k2_gen_phi"     ,k2PtrGen->phi());
               gen.addUserFloat("k2_gen_m"    ,k2PtrGen->mass());
               gen.addUserFloat("k2_gen_charge"  ,k2PtrGen->charge());
               gen.addUserInt(  "k2_gen_pdgid"   ,k2PtrGen->pdgId());
     
               gen.addUserFloat("pi_gen_px"      ,piPtrGen->px());
               gen.addUserFloat("pi_gen_py"      ,piPtrGen->py());
               gen.addUserFloat("pi_gen_pz"      ,piPtrGen->pz());
               gen.addUserFloat("pi_gen_pt"      ,piPtrGen->pt());
               gen.addUserFloat("pi_gen_eta"     ,piPtrGen->eta());
               gen.addUserFloat("pi_gen_phi"     ,piPtrGen->phi());
               gen.addUserFloat("pi_gen_m"    ,piPtrGen->mass());
               gen.addUserFloat("pi_gen_charge"  ,piPtrGen->charge());
               gen.addUserInt(  "pi_gen_pdgid"   ,piPtrGen->pdgId());

               ///////////////////////////////////////////////////////// 
               // now find the channel ID, we have the following scheme:
              
               // SIGNAL:
               // Bs -> Ds mu   0
               // Bs -> Ds tau  1
               //
               // Bs -> Ds* mu  5
               // Bs -> Ds* tau 6
               //
               // HB: (the mother base is defined in step1 for each mother)
               //
               // B mother -> Ds  +  D                              mother_base + 0
               // B mother -> Ds  +  D*                             mother_base + 1
               // B mother -> Ds  +  (something else producing mu)  mother_base + 2

               // B mother -> Ds* +  D                              mother_base + 5
               // B mother -> Ds* +  D*                             mother_base + 6
               // B mother -> Ds* +  (something else producing mu)  mother_base + 7

               // ground state charmed mesons:  
               std::vector<int> dMesons{
                 411, // "           : yes 1869
                 421, // "           : no  1864
                 431, // "           : yes 1968
                 4122 // "           : yes 4122
               };

               // excited charmed mesons:  
               std::vector<int> dMesonsExc{
                 10411, // charged/mass: yes 2400 
                 10421, // "           : no  2400
                 413,   // "           : yes 2010
                 423,   // "           : no  2007
                 10413, // "           : yes 2420
                 10423, // "           : no  2420
                 20413, // "           : yes ????
                 20423, // "           : no  2430
                 415,   // "           : yes 2460
                 425,   // "           : no  2460
                 10431, // "           : yes 2317
                 433,   // "           : yes todo
                 10433, // "           : yes 2536
                 20433, // "           : yes 2460 (resp. 2457)
                 435    // "           : yes 2573
               };


               //bId = bMotherId;

               // Step1: the b Mother fixes the 10s
               switch(abs(bMotherId)){
                 case 521:  sigId = 100;  break;  // B+
                 case 511:  sigId = 200;  break;  // B0
                 case 531:  sigId = 300;  break;  // Bs
                 case 5122: sigId = 400;  break;  // LambdaB
                 default:   sigId = 500;          // anything else
               }

               //bool isNotDoubleCharm = false;
                
               //std::cout << "candidate nr: " << nRecoCandidates << std::endl;
               int dsID = getDsID(piPtrGen); // get charmed strange ID
               //std::cout << "ds Id is: " << dsID << std::endl;
               int dID  = getSecondCharmID(muPtrGen); // get charmed ID
               //std::cout << "d Id is: " << dID << std::endl;
               bool isTau = isAncestor(muPtrGen, 15); 

               switch(dsID){
                 case 431:   sigId += 0;  break; // Ds+
                 case 433:   sigId += 10; break; // Ds+*
                 case 10431: sigId += 20; break; // Ds+(2317)*
                 case 20433: sigId += 30; break; // Ds+(2457)*
                 default:    sigId += 40; break; // anything else
               }

               // Signal candidates enter here
               if ((dID == 0) && ((dsID == 431) ||(dsID == 433))&& (abs(bMotherId) == 531)) {
                 //std::cout << "i enter the signal tag!" << std::endl;
                 sigId -= 300; // signal should live in 0
                 if (isTau) sigId += 1; 
               }

               // special case of B+ -> K nu mu / B+ -> K tau nu 
               else if((dID == 0) && (dsID == 0) && (abs(bMotherId) == 521)){
               
                 if (!isTau) sigId += 7;
                 if (isTau)  sigId  += 8;

               }

               else {
                 switch(dID){
                   case 411:   sigId += 0;  break; // D+
                   case 421:   sigId += 1; break;  // D0
                   case 431:   sigId += 2; break;  // Ds
                   case 413:   sigId += 3; break;  // D+*
                   case 423:   sigId += 4; break;  // D0*
                   case 433:   sigId += 5; break;  // Ds+*
                   case 4122:   sigId += 6; break;  // Lambda c

                   default:    sigId += 9; break; // anything else
                 }
               }


               goto fakes;
               //////////////////////////////////////////////////

             }//close gen matching pi loop 
            //break;
            }//close gen matching k2 loop 
          //break;
          } //close gen matching k1 loop
        //break;
        } //close gen matching mu loop

        //if (genMatchSuccess == 0) continue; 
 
        fakes:
          //gen.addUserInt("b_mother_id",bId);
          gen.addUserInt("gen_match_success",genMatchSuccess);
 
          if (genMatchSuccess == 0){
          //Didnt find a gen match! :( Now lets look for fakes!


            for(size_t muIdxGen = 0; muIdxGen < prunedGen->size(); ++muIdxGen){
        
              nMuonsGen++;
        
              //define a pointer to the gen muon    
              edm::Ptr<reco::GenParticle> muPtrGen(prunedGen, muIdxGen);
        
              //std::cout << "mu pt: " << muPtrGen->pt() << std::endl;
              //std::cout << "pdg id: " << muPtrGen->pdgId() << std::endl;
              //std::cout << "mother id: " << muPtrGen->mother(0)->pdgId() << std::endl;
              //std::cout << "charge: " << muPtrGen->charge() << std::endl;
        
              //select only FAKES!! 
              if( (fabs(muPtrGen->pdgId()) == 13) || muPtrGen->pt() < minMuPt_ || fabs(muPtrGen->eta()) > maxMuEta_ || (muBs->charge() * muPtrGen->charge() < 0)) continue; 
        
              muSel1CounterGen++;
        
              //now check the dR of the reco muon wrt to the gen Muon 
              float drMuonMatch = reco::deltaR(*muBs,*muPtrGen);
        
              //TODO:define as variable
              if(drMuonMatch > drMatchGen_) continue;
              muSel2CounterGen++;
        
        
              //std::cout << "yes!!" <<  std::endl;
              //std::cout << muPtrGen->pdgId() <<  std::endl;
              //std::cout << muPtrGen->pt() << std::endl;
              //std::cout << "found a gen matched muon" << std::endl;
        
              ////////////////////////////////////////////////
              // find gen matched k1                        //
              ////////////////////////////////////////////////
        
              for(size_t k1IdxGen = 0; k1IdxGen < prunedGen->size(); ++k1IdxGen){
                   
                nTracksGen++;
                //define a pointer to the gen kaon    
                edm::Ptr<reco::GenParticle> k1PtrGen(prunedGen, k1IdxGen);
        
                //select only useful kaons -> check this selection!
                if((fabs(k1PtrGen->pdgId()) != 321) || !hadSelectionGen_(*k1PtrGen) || (k1Bs->charge() * k1PtrGen->charge() < 0)) continue; 
                k1Sel1CounterGen++;
        
                //std::cout << k1PtrGen->pt() << std::endl;
                //std::cout << "passed selection for k1!" << std::endl;
                //now check the dR of the reco muon wrt to the gen Muon 
                float drK1Match = reco::deltaR(*k1Bs,*k1PtrGen);
        
                //TODO:define as variable
                if(drK1Match > drMatchGen_) continue;
                k1Sel2CounterGen++;
                //std::cout << "found a gen matched k1!" << std::endl;
        
                ////////////////////////////////////////////////
                // find gen matched k2                        //
                ////////////////////////////////////////////////
             
                for(size_t k2IdxGen = 0; k2IdxGen < prunedGen->size(); ++k2IdxGen){
             
                     //avoid picking the same gen particle as for k1
                     if(k2IdxGen == k1IdxGen) continue; 
        
                     //define a pointer to the gen kaon    
                     edm::Ptr<reco::GenParticle> k2PtrGen(prunedGen, k2IdxGen);
             
                     //select only useful kaons -> check this selection!
                     if((fabs(k2PtrGen->pdgId()) != 321) || !hadSelectionGen_(*k2PtrGen) || (k2Bs->charge() * k2PtrGen->charge() < 0 )) continue; 
                     //std::cout << "reco  pt " << k2Bs->pt() << "reco eta " << k2Bs->eta() << std::endl;
                     //std::cout << "gen pt " << k2PtrGen->pt() << "gen eta " << k2PtrGen->eta() << std::endl;
                     k2Sel1CounterGen++;
                     //std::cout << "passed selection for k2" << std::endl;
                     //now check the dR of the reco muon wrt to the gen Muon 
                     float drK2Match = reco::deltaR(*k2Bs,*k2PtrGen);
             
                     //TODO:define as variable
                     if(drK2Match > drMatchGen_ ) continue;
                     k2Sel2CounterGen++;
        
                     //std::cout << "found a gen matched k2!\n" << std::endl;
        
                     ////////////////////////////////////////////////
                     // find gen matched pion                      //
                     ////////////////////////////////////////////////
        
                     for(size_t piIdxGen = 0; piIdxGen < prunedGen->size(); ++piIdxGen){
              
                       //avoid picking the same gen particle as for k1 or k2
                       if((piIdxGen == k1IdxGen) || (piIdxGen == k2IdxGen)) continue; 
                               
                       //define a pointer to the gen kaon    
                       edm::Ptr<reco::GenParticle> piPtrGen(prunedGen, piIdxGen);
             
                       //select only useful kaons -> check this selection!
                       if((fabs(piPtrGen->pdgId()) != 211) || !hadSelectionGen_(*piPtrGen) || (piBs->charge() * piPtrGen->charge() < 0 )) continue; 
                       piSel1CounterGen++;
             
                       //std::cout << "gen pt " << piPtrGen->pt() << "gen eta " << piPtrGen->eta() << std::endl;
                       //now check the dR of the reco muon wrt to the gen Muon 
                       float drPiMatch = reco::deltaR(*piBs,*piPtrGen);
                       //std::cout << drPiMatch << std::endl; 
                       //TODO:define as variable
                       if(drPiMatch > drMatchGen_) continue;
                       //std::cout << "found a gen matched pion!" << std::endl;
                       piSel2CounterGen++;


                       //////////////////////////////////////////
                       // If we're here we found a true fake!! //
                       //////////////////////////////////////////
                       std::cout << "we have a fake!!" << std::endl;
                       sigId = -10;
                       genMatchSuccess = 1;
        
                       gen.addUserFloat("mu_gen_px"      ,muPtrGen->px());
                       gen.addUserFloat("mu_gen_py"      ,muPtrGen->py());
                       gen.addUserFloat("mu_gen_pz"      ,muPtrGen->pz());
                       gen.addUserFloat("mu_gen_pt"      ,muPtrGen->pt());
                       gen.addUserFloat("mu_gen_eta"     ,muPtrGen->eta());
                       gen.addUserFloat("mu_gen_phi"     ,muPtrGen->phi());
                       gen.addUserFloat("mu_gen_m"    ,muPtrGen->mass());
                       gen.addUserFloat("mu_gen_charge"  ,muPtrGen->charge());
                       gen.addUserInt(  "mu_gen_pdgid"   ,muPtrGen->pdgId());
             
                       gen.addUserFloat("k1_gen_px"      ,k1PtrGen->px());
                       gen.addUserFloat("k1_gen_py"      ,k1PtrGen->py());
                       gen.addUserFloat("k1_gen_pz"      ,k1PtrGen->pz());
                       gen.addUserFloat("k1_gen_pt"      ,k1PtrGen->pt());
                       gen.addUserFloat("k1_gen_eta"     ,k1PtrGen->eta());
                       gen.addUserFloat("k1_gen_phi"     ,k1PtrGen->phi());
                       gen.addUserFloat("k1_gen_m"    ,k1PtrGen->mass());
                       gen.addUserFloat("k1_gen_charge"  ,k1PtrGen->charge());
                       gen.addUserInt(  "k1_gen_pdgid"   ,k1PtrGen->pdgId());
             
                       gen.addUserFloat("k2_gen_px"      ,k2PtrGen->px());
                       gen.addUserFloat("k2_gen_py"      ,k2PtrGen->py());
                       gen.addUserFloat("k2_gen_pz"      ,k2PtrGen->pz());
                       gen.addUserFloat("k2_gen_pt"      ,k2PtrGen->pt());
                       gen.addUserFloat("k2_gen_eta"     ,k2PtrGen->eta());
                       gen.addUserFloat("k2_gen_phi"     ,k2PtrGen->phi());
                       gen.addUserFloat("k2_gen_m"    ,k2PtrGen->mass());
                       gen.addUserFloat("k2_gen_charge"  ,k2PtrGen->charge());
                       gen.addUserInt(  "k2_gen_pdgid"   ,k2PtrGen->pdgId());
             
                       gen.addUserFloat("pi_gen_px"      ,piPtrGen->px());
                       gen.addUserFloat("pi_gen_py"      ,piPtrGen->py());
                       gen.addUserFloat("pi_gen_pz"      ,piPtrGen->pz());
                       gen.addUserFloat("pi_gen_pt"      ,piPtrGen->pt());
                       gen.addUserFloat("pi_gen_eta"     ,piPtrGen->eta());
                       gen.addUserFloat("pi_gen_phi"     ,piPtrGen->phi());
                       gen.addUserFloat("pi_gen_m"    ,piPtrGen->mass());
                       gen.addUserFloat("pi_gen_charge"  ,piPtrGen->charge());
                       gen.addUserInt(  "pi_gen_pdgid"   ,piPtrGen->pdgId());

                       goto nans;          
                       } // close pi loop
                     } // close k2 loop
                   } // close k1 loop
                 } // close mu loop
               } // close if condition 


          nans: 
            gen.addUserInt("sig",sigId);
            gen.addUserInt("gen_match_success",genMatchSuccess);

            if (genMatchSuccess == 0){
   
              gen.addUserFloat("mu_gen_px"      ,std::nan(""));
              gen.addUserFloat("mu_gen_py"      ,std::nan(""));
              gen.addUserFloat("mu_gen_pz"      ,std::nan(""));
              gen.addUserFloat("mu_gen_pt"      ,std::nan(""));
              gen.addUserFloat("mu_gen_eta"     ,std::nan(""));
              gen.addUserFloat("mu_gen_phi"     ,std::nan(""));
              gen.addUserFloat("mu_gen_m"    ,std::nan(""));
              gen.addUserFloat("mu_gen_charge"  ,std::nan(""));
              gen.addUserInt(  "mu_gen_pdgid"   ,-9999);
    
              gen.addUserFloat("k1_gen_px"      ,std::nan(""));
              gen.addUserFloat("k1_gen_py"      ,std::nan(""));
              gen.addUserFloat("k1_gen_pz"      ,std::nan(""));
              gen.addUserFloat("k1_gen_pt"      ,std::nan(""));
              gen.addUserFloat("k1_gen_eta"     ,std::nan(""));
              gen.addUserFloat("k1_gen_phi"     ,std::nan(""));
              gen.addUserFloat("k1_gen_m"    ,std::nan(""));
              gen.addUserFloat("k1_gen_charge"  ,std::nan(""));
              gen.addUserInt(  "k1_gen_pdgid"   ,-9999);
    
              gen.addUserFloat("k2_gen_px"      ,std::nan(""));
              gen.addUserFloat("k2_gen_py"      ,std::nan(""));
              gen.addUserFloat("k2_gen_pz"      ,std::nan(""));
              gen.addUserFloat("k2_gen_pt"      ,std::nan(""));
              gen.addUserFloat("k2_gen_eta"     ,std::nan(""));
              gen.addUserFloat("k2_gen_phi"     ,std::nan(""));
              gen.addUserFloat("k2_gen_m"    ,std::nan(""));
              gen.addUserFloat("k2_gen_charge"  ,std::nan(""));
              gen.addUserInt(  "k2_gen_pdgid"   ,-9999);
    
              gen.addUserFloat("pi_gen_px"      ,std::nan(""));
              gen.addUserFloat("pi_gen_py"      ,std::nan(""));
              gen.addUserFloat("pi_gen_pz"      ,std::nan(""));
              gen.addUserFloat("pi_gen_pt"      ,std::nan(""));
              gen.addUserFloat("pi_gen_eta"     ,std::nan(""));
              gen.addUserFloat("pi_gen_phi"     ,std::nan(""));
              gen.addUserFloat("pi_gen_m"    ,std::nan(""));
              gen.addUserFloat("pi_gen_charge"  ,std::nan(""));
              gen.addUserInt(  "pi_gen_pdgid"   ,-9999);
   
            }
  
          /////////////////////// END OF VARIABLE DEFINITION //////////////////////
  
          //append candidate at the end of our return value :)
          //ret_value can be a vector!!
          ret_value->emplace_back(gen);
          //ret_value_gen->emplace_back(gen);
  
    } //closing loop over Bs
    iEvent.put(std::move(ret_value), "gen");
}//closing event loop


void genMatching::endJob(){
// Printouts:
std::cout << "\n--------- GEN MATCHING MODULE ----------\n" << std::endl;
std::cout << "#Events in file                                           : " << nEventsGen  << std::endl;
std::cout << "#Gen Muons in file                                        : " << nMuonsGen   << std::endl;
std::cout << "#Gen Tracks in file                                       : " << nTracksGen  << std::endl;
std::cout << "#Reco candidates found                                    : " << nRecoCandidates << std::endl;
std::cout << "#Gen matched candidates                                   : " << nGenMatched << std::endl;

std::cout << "#Gen Kaon 1 which passed the hadronic selection           : " << k1Sel1CounterGen << std::endl;
std::cout << "#Gen Kaon 1 which passed the dR matching                  : " << k1Sel2CounterGen << std::endl;
std::cout << "#Gen Kaon 2 which passed the hadronic selection           : " << k2Sel1CounterGen << std::endl;
std::cout << "#Gen Kaon 2 which passed the dR matching                  : " << k2Sel2CounterGen << std::endl;
std::cout << "#Gen Pions  which passed the hadronic selection           : " << piSel1CounterGen << std::endl;
std::cout << "#Gen Pions  which passed the dR matching                  : " << piSel2CounterGen << std::endl;

std::cout << "\n#KKPiMu Gen combinations:" << nKKPiMuGen << std::endl;
std::cout << "#KKPiMu Gen combinations for which we found a Phi         : " << nFoundPhi  << std::endl;
std::cout << "#KKPiMu Gen combinations for which we found a Ds          : " << nFoundDs   << std::endl;
std::cout << "#KKPiMu Gen combinations for which we found a B-mom       : " << nFoundB    << std::endl;
std::cout << "#KKPiMu Gen combinations for which the B-mom < B mass cut : " << nBMassCut << std::endl;
}


DEFINE_FWK_MODULE(genMatching);
