#!/usr/bin/env python

import sys, os, time, shutil, subprocess, yaml, re
import numpy as np
from extract_from_mesh_directory import *


def find_mesh_from_run_directory(dirPath):
  with open(dirPath+'/input.yaml') as file:
    ifile = yaml.load(file, Loader=yaml.FullLoader)
  return ifile['general']['meshDir']


def find_dimensionality_from_run_directory(dirPath):
  meshDir = find_mesh_from_run_directory(dirPath)
  return find_dimensionality_from_info_file(meshDir)


def find_num_cells_from_run_directory(dirPath):
  meshDir = find_mesh_from_run_directory(dirPath)
  dim = find_dimensionality_from_run_directory(dirPath)
  if dim == 1:
    nx = find_num_cells_from_info_file(meshDir, 'nx')
    return [nx]
  elif dim==2:
    nx = find_num_cells_from_info_file(meshDir, 'nx')
    ny = find_num_cells_from_info_file(meshDir, 'ny')
    return [nx,ny]
  elif dim==3:
    nx = find_num_cells_from_info_file(meshDir, 'nx')
    ny = find_num_cells_from_info_file(meshDir, 'ny')
    nz = find_num_cells_from_info_file(meshDir, 'nz')
    return [nx,ny,nz]


def find_total_dofs_sample_mesh_from_run_directory(dirPath, numDofsPerCell):
  meshDir = find_mesh_from_run_directory(dirPath)
  smCount = find_sample_mesh_count_from_info_file(meshDir)
  return int(smCount * numDofsPerCell)


def find_total_dofs_from_run_directory(dirPath, numDofsPerCell):
  meshDir = find_mesh_from_run_directory(dirPath)
  dim = find_dimensionality_from_run_directory(dirPath)
  N = find_num_cells_from_run_directory(dirPath)
  if dim == 1:
    return np.int64(N[0]*numDofsPerCell)
  elif dim==2:
    return np.int64(N[0]*N[1]*numDofsPerCell)
  elif dim==3:
    return np.int64(N[0]*N[1]*N[2]*numDofsPerCell)


def find_pod_modes_file_from_run_directory(dirPath):
  with open(dirPath+'/input.yaml') as file:
    ifile = yaml.load(file, Loader=yaml.FullLoader)
  return ifile['rom']['podFile']


def find_num_modes_from_run_directory(dirPath):
  with open(dirPath+'/input.yaml') as file:
    ifile = yaml.load(file, Loader=yaml.FullLoader)
  return int(ifile['rom']['numModes'])


def find_reference_state_file_from_run_directory(dirPath):
  with open(dirPath+'/input.yaml') as file:
    ifile = yaml.load(file, Loader=yaml.FullLoader)
  return ifile['rom']['referenceStateFile']


def find_final_time_from_fom_run_directory(dirPath):
  with open(dirPath+'/input.yaml') as file:
    ifile = yaml.load(file, Loader=yaml.FullLoader)
  return float(ifile['fom']['finalTime'])


def find_final_time_from_rom_run_directory(dirPath):
  with open(dirPath+'/input.yaml') as file:
    ifile = yaml.load(file, Loader=yaml.FullLoader)
  return float(ifile['rom']['finalTime'])


#
# def find_time_step_size_from_run_directory(dirPath):
#   with open(dirPath+'/input.yaml') as file:
#     ifile = yaml.load(file, Loader=yaml.FullLoader)
#   return float(ifile['general']['dt'])

#
# def find_sampling_freq_from_run_directory(dirPath):
#   with open(dirPath+'/input.yaml') as file:
#     ifile = yaml.load(file, Loader=yaml.FullLoader)
#   return int(ifile['general']['sf'])

#
# def find_reaction_from_run_directory(dirPath):
#   with open(dirPath+'/input.yaml') as file:
#     ifile = yaml.load(file, Loader=yaml.FullLoader)
#   return float(ifile['parameters']['reaction'])

#
# def find_diffusion_from_run_directory(dirPath):
#   with open(dirPath+'/input.yaml') as file:
#     ifile = yaml.load(file, Loader=yaml.FullLoader)
#   return float(ifile['parameters']['diffusion'])
