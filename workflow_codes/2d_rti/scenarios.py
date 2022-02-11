#!/usr/bin/env python

import numpy as np

base_dic              = {}
train_points          = {}
test_points           = {}

sample_mesh_fractions = {}
leverage_scores_betas = {}

rom_energies          = {}
rom_basis_sets        = {}

'''
GalerkinFull
GalerkinCollocation
GalerkinGappy
GalerkinMaskedCollocation
GalerkinMaskedGappy
'''
rom_algos             = {}

'''
SplitVars or DontSplitVars
'''
rom_basis_styles      = {}

interpolation_sets    = {}
polyfit_sets          = {}

################################
################################

'''
------------------------------------------------------
SCENARIO 1
------------------------------------------------------
'''
base_dic[1] = {
  'general' : {
    'meshDir': "tbd",
  },

  'fom' : {
    'finalTime': 1.5,
    'inviscidFluxReconstruction' : "Weno3",
    'odeScheme': "RungeKutta4",
    'dt' : 0.00025,
    'stateSamplingFreq' : 5,
    'velocitySamplingFreq' : 5
  },

  'rom' : {
    'finalTime': 1.5,
    'inviscidFluxReconstruction' : "Weno3",
    'odeScheme': "RungeKutta4",
    'dt' : 0.00025,
    'stateSamplingFreq' : 5,
    # the following params are here but left tbd
    # because these are changed by driver scripts
    'numModes'            : 'tbd',
    'podFile'             : 'tbd',
    'algoName'            : 'tbd',
    'romInitialStateFile' : 'tbd',
    'referenceStateFile'  : 'tbd'
  },

  'physicalCoefficients' : {
    'amplitude' : "tbd",
    'gamma' : 5./3.
  }
}

train_points[1] = {
  0: 0.025,
  1: 0.050
}

test_points[1]  = {
  0: 0.025,
  1: 0.040
}

sample_mesh_fractions[1] = [0.025, 0.05] #, 0.05, 0.10]
leverage_scores_betas[1] = [0.05, 0.5, 0.95] #, 0.95, 0.99]

rom_algos[1]        = ["GalerkinGappy"]
rom_basis_styles[1] = ["DontSplitVars"]
rom_energies[1]     = [99.99, 99.9995, 99.99999, \
                       99.999995, 99.9999999] #, 99.999999]
rom_basis_sets[1]   = { 0: [0,1] }




'''
------------------------------------------------------
SCENARIO 2
------------------------------------------------------
'''
base_dic[2] = {
  'general' : {
    'meshDir': "tbd",
  },

  'fom' : {
    'finalTime': 1.2,
    'inviscidFluxReconstruction' : "Weno3",
    'odeScheme': "RungeKutta4",
    'dt' : 0.00025,
    'stateSamplingFreq' : 5,
    'velocitySamplingFreq' : 5
  },

  'rom' : {
    'finalTime': 1.2,
    'inviscidFluxReconstruction' : "Weno3",
    'odeScheme': "RungeKutta4",
    'dt' : 0.00025,
    'stateSamplingFreq' : 5,
    # the following params are here but left tbd
    # because these are changed by driver scripts
    'numModes'            : 'tbd',
    'podFile'             : 'tbd',
    'algoName'            : 'tbd',
    'romInitialStateFile' : 'tbd',
    'referenceStateFile'  : 'tbd'
  },

  'physicalCoefficients' : {
    'amplitude' : 0.025,
    'gamma' : "tbd"
  }
}

train_points[2] = {
  0: 1.45,
  1: 1.5,
  2: 1.66
}

test_points[2]  = {
  0: 1.5,
  1: 1.55
}

sample_mesh_fractions[2] = [0.01, 0.025, 0.05, 0.10]
leverage_scores_betas[2] = [0.05, 0.5, 0.95, 0.99]

rom_algos[2]        = ["GalerkinGappy"]
rom_basis_styles[2] = ["DontSplitVars"]
rom_energies[2]     = [99.999995, 99.999999]
rom_basis_sets[2]   = {
  0: [0,2],
  #1: [0,1,2]
}


'''
male list of valid scenarios, figure out from keys in dic
'''
valid_scenarios_ids = list(base_dic.keys())
