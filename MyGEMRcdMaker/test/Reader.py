import FWCore.ParameterSet.Config as cms

#process = cms.Process("MyGEMAlignmentRcdWriter")
#process.load('Configuration.Geometry.GeometryExtended2021_cff')

from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('SIM2RAW',Run3)


process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.GEMGeometryESProducer = cms.ESProducer("GEMGeometryESModule",
                                               DDDetector = cms.ESInputTag('',''),
                                               applyAlignment = cms.bool(False),
                                               alignmentsLabel = cms.string(''),
                                               attribute = cms.string('MuStructure'),
                                               value = cms.string('MuonEndCapGEM'),
                                               useDDD = cms.bool(True),
                                               useDD4hep = cms.untracked.bool(False)                                        
                                               )
process.DDSpecParRegistryESProducer = cms.ESProducer("DDSpecParRegistryESProducer",
                                                     appendToDataLabel = cms.string('MUON')
                                                     )
process.MuonNumberingESProducer = cms.ESProducer("MuonNumberingESProducer",
                                                 label = cms.string('MUON'),
                                                 key = cms.string('MuonCommonNumbering')
                                                 )

process.load("CondCore.CondDB.CondDB_cfi")
process.OutDB = process.CondDB.clone()
process.OutDB.connect = 'sqlite_file:GEMAl.db'

process.source = cms.Source("EmptyIOVSource",
    timetype = cms.string('runnumber'),
    firstValue = cms.uint64(1),
    lastValue = cms.uint64(1),
    interval = cms.uint64(1)
)

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    process.OutDB,
    timetype = cms.untracked.string('runnumber'),
    toPut = cms.VPSet(cms.PSet(
        record = cms.string('GEMAlignmentRcd'),
        tag = cms.string('GEMAlignment_test')
    ),
    cms.PSet(
        record = cms.string('GEMAlignmentErrorExtendedRcd'),
        tag = cms.string('GEMAlignmentErrorExtended_test')
    ))
)

process.gem_maker = cms.EDAnalyzer("AlignmentReader",
    loggingOn= cms.untracked.bool(True),
    SinceAppendMode=cms.bool(True),
    Source=cms.PSet(
        IOVRun=cms.untracked.uint32(1)
    ),
    xShift = cms.double(0),		#Shift (cm)
    yShift = cms.double(0),
    zShift = cms.double(0),
    rotX = cms.double(0),		#Shift (Radians)
    rotY = cms.double(0),
    rotZ = cms.double(0)
)

process.path = cms.Path(process.gem_maker)
