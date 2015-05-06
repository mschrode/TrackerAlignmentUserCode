import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import copy, sys, os

process = cms.Process("Misaligner")

###################################################################
# Setup 'standard' options
###################################################################
options = VarParsing.VarParsing()
options.register('myGT',
                 "START53_V24::All", # default value
                 VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                 VarParsing.VarParsing.varType.string, # string, int, or float
                 "Global Tag to use (START53_V24::All is default)")

options.register('myScenario',
                 "MisalignmentScenario100Mu", # default value
                 VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                 VarParsing.VarParsing.varType.string, # string, int, or float
                 "scenario to apply")

options.parseArguments()

###################################################################
# Message logger service
###################################################################
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cout = cms.untracked.PSet(
    threshold = cms.untracked.string('INFO'),
    default = cms.untracked.PSet(
        limit = cms.untracked.int32(10000000)
    )
)
# replace MessageLogger.debugModules = { "*" }
# service = Tracer {}

###################################################################
# Ideal geometry producer and standard includes
###################################################################
process.load("Configuration.Geometry.GeometryIdeal_cff")

###################################################################
# Just state the Global Tag (and pick some run)
###################################################################
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = options.myGT #'MCRUN2_72_V0::All' # take your favourite
print "Using global tag: %s" % process.GlobalTag.globaltag._value

###################################################################
# This uses the object from the tag and applies the misalignment scenario on top of that object
###################################################################
process.load("Alignment.CommonAlignmentProducer.AlignmentProducer_cff")
process.AlignmentProducer.doMisalignmentScenario=True
process.AlignmentProducer.applyDbAlignment=True
process.AlignmentProducer.checkDbAlignmentValidity=False #otherwise error thrown for IOV dependent GTs
from TrackerAlignmentUserCode.VertexAnalysis.MisalignmentScenarios_cff import *

isMatched = False
print "Using scenario:",options.myScenario

for objname,oid in globals().items():
    #print objname
    if (str(objname) == str(options.myScenario)):
        isMatched = True
        print "Using scenario:",objname
        process.AlignmentProducer.MisalignmentScenario = oid
   
if isMatched is not True:
    print "----- Begin Fatal Exception -----------------------------------------------"
    print "Unrecognized",options.myScenario,"misalignment scenario !!!!"
    print "Aborting cmsRun now, please check your input"
    print "----- End Fatal Exception -------------------------------------------------"
    os._exit(1)
  
process.AlignmentProducer.saveToDB=True
process.AlignmentProducer.saveApeToDB=True

###################################################################
# Misalignment example scenario producer
# >>>>>>> THIS WORKS ONLY IF YOU WANT TO START FROM IDEAL
###################################################################
#process.load("Alignment.TrackerAlignment.MisalignedTracker_cfi")
#process.MisalignedTracker.saveToDbase = True # to store to DB
#import Alignment.TrackerAlignment.Scenarios_cff as _Scenarios
#process.MisalignedTracker.scenario = _Scenarios.TrackerCSA14Scenario
#process.MisalignedTracker.scenario = _Scenarios.Tracker10pbScenario
#process.MisalignedTracker.scenario = _Scenarios.SurveyLASOnlyScenario
#process.MisalignedTracker.scenario = _Scenarios.SurveyLASCosmicsScenario
#process.MisalignedTracker.scenario = _Scenarios.TrackerCRAFTScenario

###################################################################
# The module
###################################################################
process.prod = cms.EDAnalyzer("TestAnalyzer",
    fileName = cms.untracked.string('misaligned.root')
)

###################################################################
# Output name
###################################################################
outputfilename = None
scenariolabel = str(options.myScenario)

###################################################################
# Source
###################################################################
process.source = cms.Source("EmptySource",
                            #firstRun = cms.untracked.uint32(210659) # choose your run (for data only!)
                            )

outputfilename = "geometry_"+str(scenariolabel)+"_from"+process.GlobalTag.globaltag._value.replace('::All','')+"_fromRandomTool.db"

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1) )

###################################################################
# Database output service
###################################################################

import CondCore.DBCommon.CondDBSetup_cfi
process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    CondCore.DBCommon.CondDBSetup_cfi.CondDBSetup,
    # Writing to oracle needs the following shell variable setting (in zsh):
    # export CORAL_AUTH_PATH=/afs/cern.ch/cms/DB/conddb
    # connect = cms.string('oracle://cms_orcoff_prep/CMS_COND_ALIGNMENT'),  # preparation/develop. DB
    timetype = cms.untracked.string('runnumber'),
    connect = cms.string('sqlite_file:'+outputfilename),                                      
    toPut = cms.VPSet(
        cms.PSet(
            record = cms.string('TrackerAlignmentRcd'),
            tag = cms.string('Alignments')
        ), 
       cms.PSet(
           record = cms.string('TrackerAlignmentErrorRcd'),
           tag = cms.string('AlignmentErrors')
       )
    )
)
process.PoolDBOutputService.DBParameters.messageLevel = 2

process.p1 = cms.Path(process.prod)



