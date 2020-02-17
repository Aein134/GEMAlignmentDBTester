// -*- C++ -*-
//
// Package:    MyDB/MyGEMRcdMaker
// Class:      MyGEMRcdMaker
//
/**\class MyGEMRcdMaker MyGEMRcdMaker.cc MyDB/MyGEMRcdMaker/plugins/MyGEMRcdMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ian James Watson
//         Created:  Wed, 13 Feb 2019 14:29:19 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "FWCore/Utilities/interface/InputTag.h"
 #include "DataFormats/TrackReco/interface/Track.h"
 #include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "CondFormats/AlignmentRecord/interface/GEMAlignmentRcd.h"
#include "CondFormats/AlignmentRecord/interface/GEMAlignmentErrorRcd.h"
#include "CondFormats/AlignmentRecord/interface/GEMAlignmentErrorExtendedRcd.h"

#include "CondFormats/Alignment/interface/Alignments.h"
#include "CondFormats/Alignment/interface/AlignmentErrors.h"
#include "CondFormats/Alignment/interface/AlignmentErrorsExtended.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"

#include "CLHEP/Vector/EulerAngles.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Geometry/Transform3D.h"
#include "CondFormats/Alignment/interface/AlignTransform.h"

 #include "DataFormats/DetId/interface/DetId.h"
 #include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
 #include "DataFormats/MuonDetId/interface/GEMDetId.h"

#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "Alignment/CommonAlignment/interface/Utilities.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class AlignmentReader : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit AlignmentReader(const edm::ParameterSet&);
      ~AlignmentReader();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      float xShift,yShift,zShift,rotX,rotY,rotZ;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
AlignmentReader::AlignmentReader(const edm::ParameterSet& p)
{
  xShift = p.getParameter<double>("xShift");
  yShift = p.getParameter<double>("yShift");
  zShift = p.getParameter<double>("zShift");
  rotX = p.getParameter<double>("rotX");
  rotY = p.getParameter<double>("rotY");
  rotZ = p.getParameter<double>("rotZ");
}


AlignmentReader::~AlignmentReader()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

using CLHEP::Hep3Vector;
using CLHEP::HepRotation;
// ------------ method called for each event  ------------
void
AlignmentReader::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Get the GEM Geometry
  edm::ESHandle<GEMGeometry> gemGeo;
  iSetup.get<MuonGeometryRecord>().get(gemGeo);

  Alignments* MyGEMAlignment = new Alignments();
  // AlignmentErrors* MyGEMAlignmentError = new AlignmentErrors();
  AlignmentErrorsExtended* MyGEMAlignmentErrorExtended = new AlignmentErrorsExtended();

  // Form the data here

  // Needs to match with what GEMGeometry hands out, which is:
  //
  // 1 detid per eta partition
  // 1 det id for chamber with roll=0
  // 1 det id for the superchamber with layer=0, roll=0
  //

  int detNum; int endcap;
  std::string line, cell, DetNum, dx, dy, dz, dphix, dphiy, dphiz;
  std::ifstream maptype("gemAl.csv");
  while(std::getline(maptype, line)){
    std::cout << line << std::endl;
    std::stringstream ssline(line);

    getline(ssline, DetNum, ',');
    getline(ssline, dx, ',');
    getline(ssline, dy, ',');
    getline(ssline, dz, ',');
    getline(ssline, dphix, ',');
    getline(ssline, dphiy, ',');
    getline(ssline, dphiz, ',');
    //DetNum >> detNum; dx >> xShift; dy >> yShift; dz >> zShift; dphix >> rotX; dphiy >> rotY; dphiz >> rotZ;
    detNum = (float)atof(DetNum.c_str());
    xShift = (float)atof(dx.c_str());
    yShift = (float)atof(dy.c_str());
    zShift = (float)atof(dz.c_str());
    rotX = (float)atof(dphix.c_str());
    rotY = (float)atof(dphiy.c_str());
    rotZ = (float)atof(dphiz.c_str());
    std::cout << detNum << xShift << yShift << zShift << rotX << rotY << rotZ << std::endl;
    if (detNum >= 0){endcap = 1;}
    else {endcap = -1;}
    std::cout << "endcap is " << endcap << std::endl;
    GEMDetId id = GEMDetId(endcap, 1, abs(detNum/100), 0, abs(detNum%100), 0);
    auto sch = gemGeo->superChamber(id);
    std::cout << "sch->id() " << sch->id() << std::endl;
    for (auto ch : sch->chambers()){
      std::cout << "ch->id() " << ch->id() << std::endl;
      const BoundPlane& chSurface(ch->surface());
      for (auto roll : ch->etaPartitions()){
        std::cout << "roll->id() " << roll->id() << std::endl;
        auto center = roll->surface().toGlobal(LocalPoint(xShift, yShift, zShift));
        auto rot = roll->surface().rotation();
        auto hrot = HepRotation(Hep3Vector(rot.xx(), rot.xy(), rot.xz()).unit(),
                                Hep3Vector(rot.yx(), rot.yy(), rot.yz()).unit(),
                                Hep3Vector(rot.zx(), rot.zy(), rot.zz()).unit());
        hrot.rotateX(rotX); hrot.rotateY(rotY); hrot.rotateZ(rotZ);
        auto euler = hrot.inverse().eulerAngles();
        const BoundPlane& rollSurface(roll->surface());
        LocalPoint  roll_lCentre( 0., 0., 0. );
        GlobalPoint roll_gCentre(rollSurface.toGlobal(roll_lCentre));
        LocalPoint ch_lCentre(chSurface.toLocal(roll_gCentre));
        auto l_hrot = HepRotation(Hep3Vector(1, 0, 0),
                                  Hep3Vector(0, 1, 0),
                                  Hep3Vector(0, 0, 1));
        l_hrot.rotateX(rotX); l_hrot.rotateY(rotY); l_hrot.rotateZ(rotZ);
        auto oldP = Hep3Vector(0, ch_lCentre.y(), ch_lCentre.z());
        auto newP = l_hrot*oldP;
        auto oldGP(chSurface.toGlobal(LocalPoint(oldP.x(), oldP.y(), oldP.z())));
        auto newGP(chSurface.toGlobal(LocalPoint(newP.x(), newP.y(), newP.z())));
        auto dP = newGP - oldGP;
        MyGEMAlignment->m_align.push_back(AlignTransform(AlignTransform::Translation(center.x()+dP.x(), center.y()+dP.y(), center.z()+dP.z()), euler, roll->id()));
        MyGEMAlignmentErrorExtended->m_alignError.push_back(AlignTransformErrorExtended(AlignTransformErrorExtended::SymMatrix(6), roll->id()));
      }
      auto center = ch->surface().toGlobal(LocalPoint(xShift, yShift, zShift));
      auto rot = ch->surface().rotation();
      auto hrot = HepRotation(Hep3Vector(rot.xx(), rot.xy(), rot.xz()).unit(),
                              Hep3Vector(rot.yx(), rot.yy(), rot.yz()).unit(),
                              Hep3Vector(rot.zx(), rot.zy(), rot.zz()).unit());
      hrot.rotateX(rotX); hrot.rotateY(rotY); hrot.rotateZ(rotZ);
      auto euler = hrot.inverse().eulerAngles();
      MyGEMAlignment->m_align.push_back(AlignTransform(AlignTransform::Translation(center.x(), center.y(), center.z()), euler, ch->id()));
      MyGEMAlignmentErrorExtended->m_alignError.push_back(AlignTransformErrorExtended(AlignTransformErrorExtended::SymMatrix(6), ch->id()));
    }
    auto center = sch->surface().toGlobal(LocalPoint(xShift, yShift, zShift));
    auto rot = sch->surface().rotation();
    auto hrot = HepRotation(Hep3Vector(rot.xx(), rot.xy(), rot.xz()).unit(),
                            Hep3Vector(rot.yx(), rot.yy(), rot.yz()).unit(),
                            Hep3Vector(rot.zx(), rot.zy(), rot.zz()).unit());
    hrot.rotateX(rotX); hrot.rotateY(rotY); hrot.rotateZ(rotZ);
    auto euler = hrot.inverse().eulerAngles();
    MyGEMAlignment->m_align.push_back(AlignTransform(AlignTransform::Translation(center.x(), center.y(), center.z()), euler, sch->id()));
    MyGEMAlignmentErrorExtended->m_alignError.push_back(AlignTransformErrorExtended(AlignTransformErrorExtended::SymMatrix(6), sch->id()));

  }


  // GeometryAligner expects ordering by raw ID
  std::sort(MyGEMAlignment->m_align.begin(), MyGEMAlignment->m_align.end(), [](AlignTransform a, AlignTransform b){return a.rawId() < b.rawId();});
  std::sort(MyGEMAlignmentErrorExtended->m_alignError.begin(), MyGEMAlignmentErrorExtended->m_alignError.end(), [](auto a, auto b){return a.rawId() < b.rawId();});
  edm::Service<cond::service::PoolDBOutputService> poolDbService;
  if( poolDbService.isAvailable() ) {
      poolDbService->writeOne( MyGEMAlignment, poolDbService->currentTime(),
                                               "GEMAlignmentRcd"  );
      // poolDbService->writeOne( MyGEMAlignmentError, poolDbService->currentTime(),
      //                                          "GEMAlignmentErrorRcd"  );
      poolDbService->writeOne( MyGEMAlignmentErrorExtended, poolDbService->currentTime(),
                                               "GEMAlignmentErrorExtendedRcd"  );
  }
  else
      throw std::runtime_error("PoolDBService required.");
}


// ------------ method called once each job just before starting event loop  ------------
void
AlignmentReader::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
AlignmentReader::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AlignmentReader::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AlignmentReader);
