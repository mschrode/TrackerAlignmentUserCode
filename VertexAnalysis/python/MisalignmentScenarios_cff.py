import FWCore.ParameterSet.Config as cms

# -----------------------------------------------------------------------
# General settings common to all scenarios
MisalignmentScenarioSettings = cms.PSet(
    setRotations = cms.bool(True),
    setTranslations = cms.bool(True),
    seed = cms.int32(1234567),
    distribution = cms.string('gaussian'),
    setError = cms.bool(True)
)


MisalignmentScenario100Mu = cms.PSet(
    MisalignmentScenarioSettings,
    scale = cms.double(0.01),# shifts in 100mum

    TPBHalfBarrels = cms.PSet(
        DetUnits = cms.PSet(
            dZlocal = cms.double(1),
            dYlocal = cms.double(1),
            dXlocal = cms.double(1),
        ),
    ),

    TIBHalfBarrels = cms.PSet(
        DetUnits = cms.PSet(
            dZlocal = cms.double(1),
            dYlocal = cms.double(1),
            dXlocal = cms.double(1),
        ),
    ),

    TOBHalfBarrels = cms.PSet(
        DetUnits = cms.PSet(
            dZlocal = cms.double(1),
            dYlocal = cms.double(1),
            dXlocal = cms.double(1),
        ),
    ),

    TPEEndcaps = cms.PSet(
        DetUnits = cms.PSet(
            dZlocal = cms.double(1),
            dYlocal = cms.double(1),
            dXlocal = cms.double(1),
        ),
    ),

    TIDEndcaps = cms.PSet(
        DetUnits = cms.PSet(
            dZlocal = cms.double(1),
            dYlocal = cms.double(1),
            dXlocal = cms.double(1),
        ),
    ),

    TECEndcaps = cms.PSet(
        DetUnits = cms.PSet(
            dZlocal = cms.double(1),
            dYlocal = cms.double(1),
            dXlocal = cms.double(1),
        ),
    ),
)


MisalignmentScenario200Mu = cms.PSet(
    MisalignmentScenarioSettings,
    scale = cms.double(0.01),# shifts in 100mum

    TPBHalfBarrels = cms.PSet(
        DetUnits = cms.PSet(
            dZlocal = cms.double(2),
            dYlocal = cms.double(2),
            dXlocal = cms.double(2),
        ),
    ),

    TIBHalfBarrels = cms.PSet(
        DetUnits = cms.PSet(
            dZlocal = cms.double(2),
            dYlocal = cms.double(2),
            dXlocal = cms.double(2),
        ),
    ),

    TOBHalfBarrels = cms.PSet(
        DetUnits = cms.PSet(
            dZlocal = cms.double(2),
            dYlocal = cms.double(2),
            dXlocal = cms.double(2),
        ),
    ),

    TPEEndcaps = cms.PSet(
        DetUnits = cms.PSet(
            dZlocal = cms.double(2),
            dYlocal = cms.double(2),
            dXlocal = cms.double(2),
        ),
    ),

    TIDEndcaps = cms.PSet(
        DetUnits = cms.PSet(
            dZlocal = cms.double(2),
            dYlocal = cms.double(2),
            dXlocal = cms.double(2),
        ),
    ),

    TECEndcaps = cms.PSet(
        DetUnits = cms.PSet(
            dZlocal = cms.double(2),
            dYlocal = cms.double(2),
            dXlocal = cms.double(2),
        ),
    ),
)


MisalignmentScenario300Mu = cms.PSet(
    MisalignmentScenarioSettings,
    scale = cms.double(0.01),# shifts in 100mum

    TPBHalfBarrels = cms.PSet(
        DetUnits = cms.PSet(
            dZlocal = cms.double(3),
            dYlocal = cms.double(3),
            dXlocal = cms.double(3),
        ),
    ),

    TIBHalfBarrels = cms.PSet(
        DetUnits = cms.PSet(
            dZlocal = cms.double(3),
            dYlocal = cms.double(3),
            dXlocal = cms.double(3),
        ),
    ),

    TOBHalfBarrels = cms.PSet(
        DetUnits = cms.PSet(
            dZlocal = cms.double(3),
            dYlocal = cms.double(3),
            dXlocal = cms.double(3),
        ),
    ),

    TPEEndcaps = cms.PSet(
        DetUnits = cms.PSet(
            dZlocal = cms.double(3),
            dYlocal = cms.double(3),
            dXlocal = cms.double(3),
        ),
    ),

    TIDEndcaps = cms.PSet(
        DetUnits = cms.PSet(
            dZlocal = cms.double(3),
            dYlocal = cms.double(3),
            dXlocal = cms.double(3),
        ),
    ),

    TECEndcaps = cms.PSet(
        DetUnits = cms.PSet(
            dZlocal = cms.double(3),
            dYlocal = cms.double(3),
            dXlocal = cms.double(3),
        ),
    ),
)



MisalignmentScenarioBPIX100Mu = cms.PSet(
    MisalignmentScenarioSettings,
    scale = cms.double(0.01),# shifts in 100mum

    TPBHalfBarrels = cms.PSet(
        DetUnits = cms.PSet(
            dZlocal = cms.double(1),
            dYlocal = cms.double(1),
            dXlocal = cms.double(1),
        ),
    ),
)


