
import numpy as np
import scipy.sparse as sp
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee

# natural row ordering
#
#  ...
#  10 11 12 13 14
#  5  6  7  8  9
#  0  1  2  3  4
#

class NatOrdMeshRow2d:
  def __init__(self, Nx, Ny, dx, dy, \
               xL, xR, yL, yR,
               stencilSize, enablePeriodicBc):
    self.stSize_ = stencilSize
    self.numNeighbors_ = int((stencilSize-1))
    self.Nx_ = Nx
    self.Ny_ = Ny
    self.numCells_ = Nx * Ny
    self.dx_ = dx
    self.dy_ = dy
    self.xBounds_ = [xL, xR]
    self.yBounds_ = [yL, yR]
    self.x_ = np.zeros(self.numCells_)
    self.y_ = np.zeros(self.numCells_)
    self.gids_ = np.zeros(self.numCells_)
    self.G_ = {}
    self.createFullGrid()
    self.buildGraph(enablePeriodicBc)
    #self.spMat_ = convertGraphDicToSparseMatrix(self.G_)

  def getTotCells(self):
    return self.Nx_ * self.Ny_

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

  # lower-left corner is where i,j=(0,0)
  def globalIDToGiGj(self,ID):
    return [ID % self.Nx_, np.int64(ID/self.Nx_)]

  def gigjToGlobalID(self, gi, gj):
    return np.int64( gj*self.Nx_ + gi )

  def createFullGrid(self):
    ox = self.xBounds_[0] + 0.5*self.dx_
    oy = self.yBounds_[0] + 0.5*self.dy_
    for i in range(self.numCells_):
      [gi, gj] = self.globalIDToGiGj(i)
      self.x_[i] = ox + gi * self.dx_
      self.y_[i] = oy + gj * self.dy_
      self.gids_[i] = i

  def buildGraph(self, enablePeriodicBc):
    if self.numNeighbors_ == 2:
      self._buildGraphStencil3(enablePeriodicBc)
    elif self.numNeighbors_ == 4:
      self._buildGraphStencil5(enablePeriodicBc)
    elif self.numNeighbors_ == 6:
      self._buildGraphStencil7(enablePeriodicBc)

  def _buildGraphStencil3(self, pBc):
    # neighbors are listed as left, front, right, bottom
    for iPt in range(self.numCells_):
      [gi, gj] = self.globalIDToGiGj(iPt)
      tmpList = np.zeros(4, dtype=np.int64)

      # left
      if gi==0:
        tmpList[0]=iPt+self.Nx_-1 if pBc else -1
      if gi>0:
        tmpList[0]=iPt-1

      # front
      if gj==self.Ny_-1:
        tmpList[1]=iPt-self.Nx_*(self.Ny_-1) if pBc else -1
      if gj<self.Ny_-1:
        tmpList[1]=iPt+self.Nx_

      # right
      if gi==self.Nx_-1:
        tmpList[2]=iPt-self.Nx_+1 if pBc else -1
      if gi<self.Nx_-1:
        tmpList[2]=iPt+1

      # back
      if gj==0:
        tmpList[3]=iPt+self.Nx_*(self.Ny_-1) if pBc else -1
      if gj>0:
        tmpList[3]=iPt-self.Nx_

      # store currrent neighboring list
      self.G_[iPt] = tmpList

  def _buildGraphStencil5(self, pBc):
    '''
    left0 front0 right0 back0 left1 front1 right1 back1
      0     1      2      3     4    5       6      7
    '''

    for iPt in range(self.numCells_):
      [gi, gj] = self.globalIDToGiGj(iPt)

      # temporary list holding neighbors
      tmpList = np.zeros(8, dtype=np.int64)

      # left
      if gi==0:
        tmpList[0]=self.gigjToGlobalID(self.Nx_-1, gj) if pBc else -1
        tmpList[4]=self.gigjToGlobalID(self.Nx_-2, gj) if pBc else -1
      elif gi==1:
        tmpList[0]=self.gigjToGlobalID(0,          gj)
        tmpList[4]=self.gigjToGlobalID(self.Nx_-1, gj) if pBc else -1
      elif gi>1:
        tmpList[0]=self.gigjToGlobalID(gi-1,       gj)
        tmpList[4]=self.gigjToGlobalID(gi-2,       gj)

      # front
      if gj==self.Ny_-1:
        tmpList[1]=self.gigjToGlobalID(gi, 0) if pBc else -1
        tmpList[5]=self.gigjToGlobalID(gi, 1) if pBc else -1
      if gj==self.Ny_-2:
        tmpList[1]=self.gigjToGlobalID(gi, self.Ny_-1)
        tmpList[5]=self.gigjToGlobalID(gi, 0) if pBc else -1
      if gj<self.Ny_-2:
        tmpList[1]=self.gigjToGlobalID(gi, gj+1)
        tmpList[5]=self.gigjToGlobalID(gi, gj+2)

      # right
      if gi==self.Nx_-1:
        tmpList[2]=self.gigjToGlobalID(0, gj) if pBc else -1
        tmpList[6]=self.gigjToGlobalID(1, gj) if pBc else -1
      if gi==self.Nx_-2:
        tmpList[2]=self.gigjToGlobalID(self.Nx_-1, gj)
        tmpList[6]=self.gigjToGlobalID(0,          gj) if pBc else -1
      if gi<self.Nx_-2:
        tmpList[2]=self.gigjToGlobalID(gi+1, gj)
        tmpList[6]=self.gigjToGlobalID(gi+2, gj)

      # back
      if gj==0:
        tmpList[3]=self.gigjToGlobalID(gi, self.Ny_-1) if pBc else -1
        tmpList[7]=self.gigjToGlobalID(gi, self.Ny_-2) if pBc else -1
      if gj==1:
        tmpList[3]=self.gigjToGlobalID(gi, 0)
        tmpList[7]=self.gigjToGlobalID(gi, self.Ny_-1) if pBc else -1
      if gj>1:
        tmpList[3]=self.gigjToGlobalID(gi,  gj-1)
        tmpList[7]=self.gigjToGlobalID(gi,  gj-2)

      # store currrent neighboring list
      self.G_[iPt] = tmpList

  def _buildGraphStencil7(self, pBc):
    '''
    neighbors are listed as:
    l0,f0,r0,b0, l1,f1,r1,b1, l2,f2,r2,b2

    0,1,2 indicate the closeness. So 0 is for the closest points,
    1 is for one layer out, and 2 is for 2 layers out

    left0 front0 right0 back0 left1 front1 right1 back1 left2 front2 right2 back2
      0     1      2      3     4    5       6      7     8     9      10     11
    '''
    for iPt in range(self.numCells_):
      [gi, gj] = self.globalIDToGiGj(iPt)
      tmpList = np.zeros(12, dtype=np.int64)

      # left
      if gi==0:
        tmpList[0]=self.gigjToGlobalID(self.Nx_-1, gj) if pBc else -1
        tmpList[4]=self.gigjToGlobalID(self.Nx_-2, gj) if pBc else -1
        tmpList[8]=self.gigjToGlobalID(self.Nx_-3, gj) if pBc else -1
      elif gi==1:
        tmpList[0]=self.gigjToGlobalID(0, gj)
        tmpList[4]=self.gigjToGlobalID(self.Nx_-1, gj) if pBc else -1
        tmpList[8]=self.gigjToGlobalID(self.Nx_-2, gj) if pBc else -1
      elif gi==2:
        tmpList[0]=self.gigjToGlobalID(gi-1, gj)
        tmpList[4]=self.gigjToGlobalID(gi-2, gj)
        tmpList[8]=self.gigjToGlobalID(self.Nx_-1, gj) if pBc else -1
      elif gi>2:
        tmpList[0]=self.gigjToGlobalID(gi-1, gj)
        tmpList[4]=self.gigjToGlobalID(gi-2, gj)
        tmpList[8]=self.gigjToGlobalID(gi-3, gj)

      # front
      if gj==self.Ny_-1:
        tmpList[1]=self.gigjToGlobalID(gi, 0) if pBc else -1
        tmpList[5]=self.gigjToGlobalID(gi, 1) if pBc else -1
        tmpList[9]=self.gigjToGlobalID(gi, 2) if pBc else -1
      if gj==self.Ny_-2:
        tmpList[1]=self.gigjToGlobalID(gi, self.Ny_-1)
        tmpList[5]=self.gigjToGlobalID(gi, 0) if pBc else -1
        tmpList[9]=self.gigjToGlobalID(gi, 1) if pBc else -1
      if gj==self.Ny_-3:
        tmpList[1]=self.gigjToGlobalID(gi, self.Ny_-2)
        tmpList[5]=self.gigjToGlobalID(gi, self.Ny_-1)
        tmpList[9]=self.gigjToGlobalID(gi, 0) if pBc else -1
      if gj<self.Ny_-3:
        tmpList[1]=self.gigjToGlobalID(gi, gj+1)
        tmpList[5]=self.gigjToGlobalID(gi, gj+2)
        tmpList[9]=self.gigjToGlobalID(gi, gj+3)

      # right
      if gi==self.Nx_-1:
        tmpList[2]=self.gigjToGlobalID(0, gj) if pBc else -1
        tmpList[6]=self.gigjToGlobalID(1, gj) if pBc else -1
        tmpList[10]=self.gigjToGlobalID(2, gj) if pBc else -1
      if gi==self.Nx_-2:
        tmpList[2]=self.gigjToGlobalID(self.Nx_-1, gj)
        tmpList[6]=self.gigjToGlobalID(0, gj) if pBc else -1
        tmpList[10]=self.gigjToGlobalID(1, gj) if pBc else -1
      if gi==self.Nx_-3:
        tmpList[2]=self.gigjToGlobalID(self.Nx_-2, gj)
        tmpList[6]=self.gigjToGlobalID(self.Nx_-1, gj)
        tmpList[10]=self.gigjToGlobalID(0, gj) if pBc else -1
      if gi<self.Nx_-3:
        tmpList[2]=self.gigjToGlobalID(gi+1, gj)
        tmpList[6]=self.gigjToGlobalID(gi+2, gj)
        tmpList[10]=self.gigjToGlobalID(gi+3, gj)

      # back
      if gj==0:
        tmpList[3]=self.gigjToGlobalID(gi, self.Ny_-1) if pBc else -1
        tmpList[7]=self.gigjToGlobalID(gi, self.Ny_-2) if pBc else -1
        tmpList[11]=self.gigjToGlobalID(gi, self.Ny_-3) if pBc else -1
      if gj==1:
        tmpList[3]=self.gigjToGlobalID(gi, 0)
        tmpList[7]=self.gigjToGlobalID(gi, self.Ny_-1) if pBc else -1
        tmpList[11]=self.gigjToGlobalID(gi, self.Ny_-2) if pBc else -1
      if gj==2:
        tmpList[3]=self.gigjToGlobalID(gi, 1)
        tmpList[7]=self.gigjToGlobalID(gi, 0)
        tmpList[11]=self.gigjToGlobalID(gi, self.Ny_-1) if pBc else -1
      if gj>2:
        tmpList[3]=self.gigjToGlobalID(gi,  gj-1)
        tmpList[7]=self.gigjToGlobalID(gi,  gj-2)
        tmpList[11]=self.gigjToGlobalID(gi, gj-3)

      # store currrent neighboring list
      self.G_[iPt] = tmpList



# # natural col ordering
# #
# #  ...
# #  2 5 8 11 14 ...
# #  1 4 7 10 13 ...
# #  0 3 6 9  12 ...
# #

# class NatOrdMeshCol:
#   def __init__(self, Nx, Ny, dx, dy, \
#                xL, xR, yL, yR, stencilSize):
#     self.stSize_ = stencilSize
#     self.numNeighbors_ = int((stencilSize-1)/2)
#     self.Nx_ = Nx
#     self.Ny_ = Ny
#     self.numCells_ = Nx * Ny
#     self.dx_ = dx
#     self.dy_ = dy
#     self.xBounds_ = [xL, xR]
#     self.yBounds_ = [yL, yR]
#     self.x_ = np.zeros(self.numCells_)
#     self.y_ = np.zeros(self.numCells_)
#     self.gids_ = np.zeros(self.numCells_)
#     self.G_ = {}

#     self.createFullGrid()
#     self.buildGraph()
#     self.spMat_ = convertGraphDicToSparseMatrix(self.G_)

#   def getXY(self):
#     return [self.x_, self.y_]

#   def getGIDs(self):
#     return self.gids_

#   # return the grap (dictionary) repr of the connectivity
#   def getGraph(self):
#     return self.G_

#   # return the sparse matrix repr of the connectivity
#   def getSparseMatrixRepr(self):
#     return self.spMat_

#   # lower-left corner is where i,j=(0,0)
#   def globalIDToGiGj(self,ID):
#     return [np.int64(ID/self.Ny_), ID % self.Ny_]

#   # lower-left corner is where origin is
#   def createFullGrid(self):
#     ox = self.xBounds_[0] + 0.5*self.dx_
#     oy = self.yBounds_[0] + 0.5*self.dy_
#     for i in range(self.numCells_):
#       [gi, gj] = self.globalIDToGiGj(i)
#       self.x_[i] = ox + gi * self.dx_
#       self.y_[i] = oy + gj * self.dy_
#       self.gids_[i] = i

#   def buildGraph(self):
#     # neighbors are listed as west, north, east, south
#     # since we have PBC, ensure these are met
#     for iPt in range(self.numCells_):
#       # convert from globa enumeration to (i,j)
#       [gi, gj] = self.globalIDToGiGj(iPt)

#       # temporary list holding neighbors
#       tmpList = np.zeros(4, dtype=np.int64)

#       # west neighbor
#       # if gi==0, we are on left BD
#       if gi==0: tmpList[0]=iPt+self.Ny_*(self.Nx_-1)
#       if gi>0 : tmpList[0]=iPt-self.Ny_

#       # # north naighbor
#       # # if gj==self.Ny_-1, we are on TOP BD
#       if gj==self.Ny_-1: tmpList[1]=iPt-self.Ny_+1
#       if gj<self.Ny_-1 : tmpList[1]=iPt+1

#       # east neighbor
#       # if gi==self.Nx_-1, we are on Right BD
#       if gi==self.Nx_-1: tmpList[2]=iPt-self.Ny_*(self.Nx_-1)
#       if gi<self.Nx_-1 : tmpList[2]=iPt+self.Ny_

#       # south naighbor
#       # if gj==0, we are on bottom BD
#       if gj==0: tmpList[3]=iPt+self.Ny_-1
#       if gj>0 : tmpList[3]=iPt-1

#       # store currrent neighboring list
#       self.G_[iPt] = tmpList
