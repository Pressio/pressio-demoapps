
import matplotlib.pyplot as plt
import sys, os, time
import numpy as np
from numpy import linspace, meshgrid
from matplotlib import cm
import collections
import random

from mesh_utils import *

class OneDimMesh:
  def __init__(self, Nx, dx, xL, xR, stencilSize, enablePeriodicBc):
    self.stSize_ = stencilSize
    self.numNeighbors_ = int((stencilSize-1))
    self.Nx_ = Nx
    self.numCells_ = Nx
    self.dx_ = dx
    self.xBounds_ = [xL, xR]
    self.x_ = np.zeros(self.numCells_)
    self.y_ = np.zeros(self.numCells_)
    self.gids_ = np.zeros(self.numCells_)
    self.G_ = {}
    self.createFullGrid()

    # these tags are used when we deal with nonperiodic BC
    self.cellOutsideWest_  = -1
    self.cellOutsideEast_  = -1
    self.buildGraph(enablePeriodicBc)

  def getCoordinates(self):
    return [self.x_, self.y_]

  def getGIDs(self):
    return self.gids_

  # return the grap (dictionary) repr of the connectivity
  def getGraph(self):
    return self.G_

  # return the sparse matrix repr of the connectivity
  def getSparseMatrixRepr(self):
    return self.spMat_

  def gigjToGlobalID(self, gi, gj):
    return np.int64( gj*self.Nx_ + gi )

  def createFullGrid(self):
    ox = self.xBounds_[0] + 0.5*self.dx_
    for i in range(self.numCells_):
      self.x_[i] = ox + i * self.dx_
      self.gids_[i] = i

  def buildGraph(self, enablePeriodicBc):
    if self.numNeighbors_ == 2:
      self._buildGraphStencil3(enablePeriodicBc)
    elif self.numNeighbors_ == 4:
      self._buildGraphStencil5(enablePeriodicBc)
    elif self.numNeighbors_ == 6:
      self._buildGraphStencil7(enablePeriodicBc)

  def _buildGraphStencil3(self, pBc):
    # neighbors are listed as west, north, east, south
    # since we have PBC, ensure these are met
    for iPt in range(self.numCells_):
      tmpList = np.zeros(2, dtype=np.int64)

      # west neighbor
      # if iPt==0, we are on left BD
      if iPt==0:
        tmpList[0]=iPt+self.Nx_-1 if pBc else self.cellOutsideWest_
      if iPt>0:
        tmpList[0]=iPt-1

      # east neighbor
      # if iPt==self.Nx_-1, we are on Right BD
      if iPt==self.Nx_-1:
        tmpList[1]=iPt-self.Nx_+1 if pBc else self.cellOutsideEast_
      if iPt<self.Nx_-1:
        tmpList[1]=iPt+1

      # store currrent neighboring list
      self.G_[iPt] = tmpList

  # def _buildGraphStencil5(self, pBc):
  #   for iPt in range(self.numCells_):
  #     # convert to (i,j)
  #     [gi, gj] = self.globalIDToGiGj(iPt)

  #     # temporary list holding neighbors
  #     tmpList = np.zeros(8, dtype=np.int64)

  #     # west neighbors
  #     # if gi==0, we are on left BD
  #     if gi==0:
  #       tmpList[0]=iPt+self.Nx_-1
  #       tmpList[4]=iPt+self.Nx_-2
  #     # if gi==1, we are one inside from left BD
  #     elif gi==1:
  #       tmpList[0]=iPt-1
  #       tmpList[4]=iPt+self.Nx_-2
  #     elif gi>1:
  #       tmpList[0]=iPt-1
  #       tmpList[4]=iPt-2

  #     # east neighbor
  #     # if gi==self.Nx_-1, we are on Right BD
  #     if gi==self.Nx_-1:
  #       tmpList[2]=iPt-self.Nx_+1
  #       tmpList[6]=tmpList[2]+1
  #     if gi==self.Nx_-2:
  #       tmpList[2]=iPt+1
  #       tmpList[6]=iPt-self.Nx_+2
  #     if gi<self.Nx_-2:
  #       tmpList[2]=iPt+1
  #       tmpList[6]=iPt+2

  #     # store currrent neighboring list
  #     self.G_[iPt] = tmpList

  def _buildGraphStencil7(self, pBc):
    for iPt in range(self.numCells_):
      tmpList = np.zeros(6, dtype=np.int64)

      # west neighbors
      # if iPt==0, we are on left BD
      if iPt==0:
        tmpList[0]=self.Nx_-1 if pBc else self.cellOutsideWest_
        tmpList[2]=self.Nx_-2 if pBc else self.cellOutsideWest_
        tmpList[4]=self.Nx_-3 if pBc else self.cellOutsideWest_
      elif iPt==1:
        tmpList[0]=0
        tmpList[2]=self.Nx_-1 if pBc else self.cellOutsideWest_
        tmpList[4]=self.Nx_-2 if pBc else self.cellOutsideWest_
      elif iPt==2:
        tmpList[0]=iPt-1
        tmpList[2]=iPt-2
        tmpList[4]=self.Nx_-1 if pBc else self.cellOutsideWest_
      elif iPt>2:
        tmpList[0]=iPt-1
        tmpList[2]=iPt-2
        tmpList[4]=iPt-3

      # east neighbor
      # if iPt==self.Nx_-1, we are on Right BD
      if iPt==self.Nx_-1:
        tmpList[1]=0 if pBc else self.cellOutsideEast_
        tmpList[3]=1 if pBc else self.cellOutsideEast_
        tmpList[5]=2 if pBc else self.cellOutsideEast_
      if iPt==self.Nx_-2:
        tmpList[1]=self.Nx_-1
        tmpList[3]=0 if pBc else self.cellOutsideEast_
        tmpList[5]=1 if pBc else self.cellOutsideEast_
      if iPt==self.Nx_-3:
        tmpList[1]=self.Nx_-2
        tmpList[3]=self.Nx_-1
        tmpList[5]=0 if pBc else self.cellOutsideEast_
      if iPt<self.Nx_-3:
        tmpList[1]=iPt+1
        tmpList[3]=iPt+2
        tmpList[5]=iPt+3

      # store currrent neighboring list
      self.G_[iPt] = tmpList
