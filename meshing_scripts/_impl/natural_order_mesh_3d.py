
import numpy as np

class NatOrdMesh3d:
  def __init__(self, Nx, Ny, Nz, \
               dx, dy, dz, \
               xL, xR, yL, yR, zL, zR, \
               stencilSize, \
               enablePeriodicBc):
    self.stSize_ = stencilSize
    self.numNeighbors_ = int((stencilSize-1))
    self.Nx_ = Nx
    self.Ny_ = Ny
    self.Nz_ = Nz
    self.numCells_ = Nx * Ny * Nz
    self.dx_ = dx
    self.dy_ = dy
    self.dz_ = dz
    self.xBounds_ = [xL, xR]
    self.yBounds_ = [yL, yR]
    self.zBounds_ = [zL, zR]
    self.x_ = np.zeros(self.numCells_)
    self.y_ = np.zeros(self.numCells_)
    self.z_ = np.zeros(self.numCells_)
    self.gids_ = []
    self.G_ = {}
    self.createFullGrid()
    self.buildGraph(enablePeriodicBc)

  def getTotCells(self):
    return self.Nx_ * self.Ny_ * self.Nz_

  def getCoordinates(self):
    return [self.x_, self.y_, self.z_]

  def getGIDs(self):
    return self.gids_

  # return the grap (dictionary) repr of the connectivity
  def getGraph(self):
    return self.G_

  def globalIDToGiGjGk(self,ID):
    gi = np.int64(ID % (self.Nx_*self.Ny_) % self.Nx_)
    gj = np.int64((ID % (self.Nx_*self.Ny_)) / self.Nx_)
    gk = np.int64(ID / (self.Nx_*self.Ny_))
    return [gi,gj,gk]

  def gigjgkToGlobalID(self, gi, gj, gk):
    return np.int64( gk*self.Nx_*self.Ny_ + gj*self.Nx_ + gi )

  def createFullGrid(self):
    ox = self.xBounds_[0] + 0.5*self.dx_
    oy = self.yBounds_[0] + 0.5*self.dy_
    oz = self.zBounds_[0] + 0.5*self.dz_
    for ptID in range(self.numCells_):
      [gi, gj, gk] = self.globalIDToGiGjGk(ptID)
      self.x_[ptID] = ox + gi * self.dx_
      self.y_[ptID] = oy + gj * self.dy_
      self.z_[ptID] = oz + gk * self.dz_
      self.gids_.append(ptID)
    assert(len(self.gids_)==self.numCells_)
    self.gids_ = np.array(self.gids_)

  def buildGraph(self, enablePeriodicBc):
    if self.numNeighbors_ == 2:
      self._buildGraphStencil3(enablePeriodicBc)
    elif self.numNeighbors_ == 4:
      self._buildGraphStencil5(enablePeriodicBc)

  def _buildGraphStencil3(self, pBc):
    for iPt in range(self.numCells_):
      [gi, gj, gk] = self.globalIDToGiGjGk(iPt)
      tmpList = np.zeros(6, dtype=np.int64)

      '''
      (z)  (y)
      k   j
      |  /
      | /
      |/
      ----- i (x)

      left,right: along i (x)
      back,front: along j (y)
      bottom,top: along k (z)
      '''

      # left neighbor
      if gi==0:
        tmpList[0]=self.gigjgkToGlobalID(self.Nx_-1, gj, gk) if pBc else -1
      if gi>0:
        tmpList[0]=self.gigjgkToGlobalID(gi-1, gj, gk)

      # front
      if gj==self.Ny_-1:
        tmpList[1]=self.gigjgkToGlobalID(gi, 0, gk) if pBc else -1
      if gj<self.Ny_-1:
        tmpList[1]=self.gigjgkToGlobalID(gi, gj+1, gk)

      # right
      if gi==self.Nx_-1:
        tmpList[2]=self.gigjgkToGlobalID(0, gj, gk) if pBc else -1
      if gi<self.Nx_-1:
        tmpList[2]=self.gigjgkToGlobalID(gi+1, gj, gk)

      # back
      if gj==0:
        tmpList[3]=self.gigjgkToGlobalID(gi, self.Ny_-1, gk) if pBc else -1
      if gj>0:
        tmpList[3]=self.gigjgkToGlobalID(gi, gj-1, gk)

      # bottom
      if gk==0:
        tmpList[4]=self.gigjgkToGlobalID(gi, gj, self.Nz_-1) if pBc else -1
      if gk>0:
        tmpList[4]=self.gigjgkToGlobalID(gi, gj, gk-1)

      # top
      if gk==self.Nz_-1:
        tmpList[5]=self.gigjgkToGlobalID(gi, gj, 0) if pBc else -1
      if gk<self.Nz_-1:
        tmpList[5]=self.gigjgkToGlobalID(gi, gj, gk+1)

      # store currrent neighboring list
      self.G_[iPt] = tmpList


  def _buildGraphStencil5(self, pBc):
    for iPt in range(self.numCells_):
      [gi, gj, gk] = self.globalIDToGiGjGk(iPt)
      tmpList = np.zeros(12, dtype=np.int64)

      '''
      (z)  (y)
      k   j
      |  /
      | /
      |/
      ----- i (x)

      left,right: along i (x)
      back,front: along j (y)
      bottom,top: along k (z)

      l0 f0 r0 ba0 bot0 top0 l1 f1 r1 ba1 bot1 top1
      0   1  2  3   4    5    6  7  8  9   10   11
      '''

      # left neighbor
      if gi==0:
        tmpList[0]=self.gigjgkToGlobalID(self.Nx_-1, gj,  gk) if pBc else -1
        tmpList[6]=self.gigjgkToGlobalID(self.Nx_-2, gj,  gk) if pBc else -1
      if gi==1:
        tmpList[0]=self.gigjgkToGlobalID(0,          gj,  gk)
        tmpList[6]=self.gigjgkToGlobalID(self.Nx_-1, gj,  gk) if pBc else -1
      if gi>1:
        tmpList[0]=self.gigjgkToGlobalID(gi-1,       gj,  gk)
        tmpList[6]=self.gigjgkToGlobalID(gi-2,       gj,  gk)

      # front
      if gj==self.Ny_-1:
        tmpList[1]=self.gigjgkToGlobalID(gi,          0,  gk) if pBc else -1
        tmpList[7]=self.gigjgkToGlobalID(gi,          1,  gk) if pBc else -1
      if gj==self.Ny_-2:
        tmpList[1]=self.gigjgkToGlobalID(gi, self.Ny_-1,  gk)
        tmpList[7]=self.gigjgkToGlobalID(gi,          0,  gk) if pBc else -1
      if gj<self.Ny_-2:
        tmpList[1]=self.gigjgkToGlobalID(gi,       gj+1,  gk)
        tmpList[7]=self.gigjgkToGlobalID(gi,       gj+2,  gk)

      # right
      if gi==self.Nx_-1:
        tmpList[2]=self.gigjgkToGlobalID(0,          gj,  gk) if pBc else -1
        tmpList[8]=self.gigjgkToGlobalID(1,          gj,  gk) if pBc else -1
      if gi==self.Nx_-2:
        tmpList[2]=self.gigjgkToGlobalID(self.Nx_-1, gj,  gk)
        tmpList[8]=self.gigjgkToGlobalID(0,          gj,  gk) if pBc else -1
      if gi<self.Nx_-2:
        tmpList[2]=self.gigjgkToGlobalID(gi+1,       gj,  gk)
        tmpList[8]=self.gigjgkToGlobalID(gi+2,       gj,  gk)

      # back
      if gj==0:
        tmpList[3]=self.gigjgkToGlobalID(gi, self.Ny_-1,  gk) if pBc else -1
        tmpList[9]=self.gigjgkToGlobalID(gi, self.Ny_-2,  gk) if pBc else -1
      if gj==1:
        tmpList[3]=self.gigjgkToGlobalID(gi,          0,  gk)
        tmpList[9]=self.gigjgkToGlobalID(gi, self.Ny_-1,  gk) if pBc else -1
      if gj>1:
        tmpList[3]=self.gigjgkToGlobalID(gi,       gj-1,  gk)
        tmpList[9]=self.gigjgkToGlobalID(gi,       gj-2,  gk)

      # bottom
      if gk==0:
        tmpList[4] =self.gigjgkToGlobalID(gi,        gj, self.Nz_-1) if pBc else -1
        tmpList[10]=self.gigjgkToGlobalID(gi,        gj, self.Nz_-2) if pBc else -1
      if gk==1:
        tmpList[4] =self.gigjgkToGlobalID(gi,        gj, 0,        )
        tmpList[10]=self.gigjgkToGlobalID(gi,        gj, self.Nz_-1) if pBc else -1
      if gk>1:
        tmpList[4] =self.gigjgkToGlobalID(gi,        gj, gk-1)
        tmpList[10]=self.gigjgkToGlobalID(gi,        gj, gk-2)

      # top
      if gk==self.Nz_-1:
        tmpList[5] =self.gigjgkToGlobalID(gi,        gj, 0) if pBc else -1
        tmpList[11]=self.gigjgkToGlobalID(gi,        gj, 1) if pBc else -1
      if gk==self.Nz_-2:
        tmpList[5] =self.gigjgkToGlobalID(gi,        gj, self.Nz_-1)
        tmpList[11]=self.gigjgkToGlobalID(gi,        gj, 0         ) if pBc else -1
      if gk<self.Nz_-2:
        tmpList[5] =self.gigjgkToGlobalID(gi,        gj, gk+1)
        tmpList[11]=self.gigjgkToGlobalID(gi,        gj, gk+2)

      # store currrent neighboring list
      self.G_[iPt] = tmpList


  # def _buildGraphStencil7(self, pBc):
  #   for iPt in range(self.numCells_):
  #     [gi, gj, gk] = self.globalIDToGiGjGk(iPt)
  #     tmpList = np.zeros(18, dtype=np.int64)

  #     '''
  #     (z)  (y)
  #     k   j
  #     |  /
  #     | /
  #     |/
  #     ----- i (x)

  #     left,right: along i (x)
  #     back,front: along j (y)
  #     bottom,top: along k (z)

  #     For 3d, graph is:
  #     l0 f0 r0 ba0 bot0 top0 l1 f1 r1 ba1 bot1 top1 l2 f2 r2 ba2 bot2 top2
  #     0   1  2  3   4    5    6  7  8  9   10   11  12 13 14  15  16   17
  #     '''

  #     # left
  #     if gi==0:
  #       tmpList[0] =self.gigjgkToGlobalID(self.Nx_-1, gj, gk) if pBc else -1
  #       tmpList[6] =self.gigjgkToGlobalID(self.Nx_-2, gj, gk) if pBc else -1
  #       tmpList[12]=self.gigjgkToGlobalID(self.Nx_-3, gj, gk) if pBc else -1
  #     elif gi==1:
  #       tmpList[0] =self.gigjgkToGlobalID(0,          gj, gk)
  #       tmpList[6] =self.gigjgkToGlobalID(self.Nx_-1, gj, gk) if pBc else -1
  #       tmpList[12]=self.gigjgkToGlobalID(self.Nx_-2, gj, gk) if pBc else -1
  #     elif gi==2:
  #       tmpList[0] =self.gigjgkToGlobalID(gi-1,       gj, gk)
  #       tmpList[6] =self.gigjgkToGlobalID(gi-2,       gj, gk)
  #       tmpList[12]=self.gigjgkToGlobalID(self.Nx_-1, gj, gk) if pBc else -1
  #     elif gi>2:
  #       tmpList[0] =self.gigjgkToGlobalID(gi-1,       gj, gk)
  #       tmpList[6] =self.gigjgkToGlobalID(gi-2,       gj, gk)
  #       tmpList[12]=self.gigjgkToGlobalID(gi-3,       gj, gk)

  #     # front
  #     if gj==self.Ny_-1:
  #       tmpList[1] =self.gigjgkToGlobalID(gi, 0,          gk) if pBc else -1
  #       tmpList[7] =self.gigjgkToGlobalID(gi, 1,          gk) if pBc else -1
  #       tmpList[13]=self.gigjgkToGlobalID(gi, 2,          gk) if pBc else -1
  #     if gj==self.Ny_-2:
  #       tmpList[1] =self.gigjgkToGlobalID(gi, self.Ny_-1, gk)
  #       tmpList[7] =self.gigjgkToGlobalID(gi, 0,          gk) if pBc else -1
  #       tmpList[13]=self.gigjgkToGlobalID(gi, 1,          gk) if pBc else -1
  #     if gj==self.Ny_-3:
  #       tmpList[1] =self.gigjgkToGlobalID(gi, self.Ny_-2, gk)
  #       tmpList[7] =self.gigjgkToGlobalID(gi, self.Ny_-1, gk)
  #       tmpList[13]=self.gigjgkToGlobalID(gi, 0,          gk) if pBc else -1
  #     if gj<self.Ny_-3:
  #       tmpList[1] =self.gigjgkToGlobalID(gi, gj+1,       gk)
  #       tmpList[7] =self.gigjgkToGlobalID(gi, gj+2,       gk)
  #       tmpList[13]=self.gigjgkToGlobalID(gi, gj+3,       gk)

  #     # right
  #     if gi==self.Nx_-1:
  #       tmpList[2] =self.gigjgkToGlobalID(0,          gj, gk) if pBc else -1
  #       tmpList[8] =self.gigjgkToGlobalID(1,          gj, gk) if pBc else -1
  #       tmpList[14]=self.gigjgkToGlobalID(2,          gj, gk) if pBc else -1
  #     if gi==self.Nx_-2:
  #       tmpList[2] =self.gigjgkToGlobalID(self.Nx_-1, gj, gk)
  #       tmpList[8] =self.gigjgkToGlobalID(0,          gj, gk) if pBc else -1
  #       tmpList[14]=self.gigjgkToGlobalID(1,          gj, gk) if pBc else -1
  #     if gi==self.Nx_-3:
  #       tmpList[2] =self.gigjgkToGlobalID(self.Nx_-2, gj, gk)
  #       tmpList[8] =self.gigjgkToGlobalID(self.Nx_-1, gj, gk)
  #       tmpList[14]=self.gigjgkToGlobalID(0,          gj, gk) if pBc else -1
  #     if gi<self.Nx_-3:
  #       tmpList[2] =self.gigjgkToGlobalID(gi+1,       gj, gk)
  #       tmpList[8] =self.gigjgkToGlobalID(gi+2,       gj, gk)
  #       tmpList[14]=self.gigjgkToGlobalID(gi+3,       gj, gk)

  #     # back
  #     if gj==0:
  #       tmpList[3] =self.gigjgkToGlobalID(gi,  self.Ny_-1, gk) if pBc else -1
  #       tmpList[9] =self.gigjgkToGlobalID(gi,  self.Ny_-2, gk) if pBc else -1
  #       tmpList[15]=self.gigjgkToGlobalID(gi,  self.Ny_-3, gk) if pBc else -1
  #     if gj==1:
  #       tmpList[3] =self.gigjgkToGlobalID(gi,  0,          gk)
  #       tmpList[9] =self.gigjgkToGlobalID(gi,  self.Ny_-1, gk) if pBc else -1
  #       tmpList[15]=self.gigjgkToGlobalID(gi,  self.Ny_-2, gk) if pBc else -1
  #     if gj==2:
  #       tmpList[3] =self.gigjgkToGlobalID(gi,  1,          gk)
  #       tmpList[9] =self.gigjgkToGlobalID(gi,  0,          gk)
  #       tmpList[15]=self.gigjgkToGlobalID(gi,  self.Ny_-1, gk) if pBc else -1
  #     if gj>2:
  #       tmpList[3] =self.gigjgkToGlobalID(gi,  gj-1,       gk)
  #       tmpList[9] =self.gigjgkToGlobalID(gi,  gj-2,       gk)
  #       tmpList[15]=self.gigjgkToGlobalID(gi,  gj-3,       gk)

  #     # bottom
  #     if gk==0:
  #       tmpList[4] =self.gigjgkToGlobalID(gi,  gj,         self.Nz_-1) if pBc else -1
  #       tmpList[10]=self.gigjgkToGlobalID(gi,  gj,         self.Nz_-2) if pBc else -1
  #       tmpList[16]=self.gigjgkToGlobalID(gi,  gj,         self.Nz_-3) if pBc else -1
  #     if gk==1:
  #       tmpList[4] =self.gigjgkToGlobalID(gi,  gj,         0)
  #       tmpList[10]=self.gigjgkToGlobalID(gi,  gj,         self.Nz_-1) if pBc else -1
  #       tmpList[16]=self.gigjgkToGlobalID(gi,  gj,         self.Nz_-2) if pBc else -1
  #     if gk==2:
  #       tmpList[4] =self.gigjgkToGlobalID(gi,  gj,         1)
  #       tmpList[10]=self.gigjgkToGlobalID(gi,  gj,         0)
  #       tmpList[16]=self.gigjgkToGlobalID(gi,  gj,         self.Nz_-1) if pBc else -1
  #     if gk>2:
  #       tmpList[4] =self.gigjgkToGlobalID(gi,  gj,         gk-1)
  #       tmpList[10]=self.gigjgkToGlobalID(gi,  gj,         gk-2)
  #       tmpList[16]=self.gigjgkToGlobalID(gi,  gj,         gk-3)

  #     # top
  #     if gk==self.Nz_-1:
  #       tmpList[5] =self.gigjgkToGlobalID(gi, gj, 0) if pBc else -1
  #       tmpList[11]=self.gigjgkToGlobalID(gi, gj, 1) if pBc else -1
  #       tmpList[17]=self.gigjgkToGlobalID(gi, gj, 2) if pBc else -1
  #     if gk==self.Nz_-2:
  #       tmpList[5] =self.gigjgkToGlobalID(gi, gj, self.Nz_-1)
  #       tmpList[11]=self.gigjgkToGlobalID(gi, gj, 0) if pBc else -1
  #       tmpList[17]=self.gigjgkToGlobalID(gi, gj, 1) if pBc else -1
  #     if gk==self.Nz_-3:
  #       tmpList[5] =self.gigjgkToGlobalID(gi, gj, self.Nz_-2)
  #       tmpList[11]=self.gigjgkToGlobalID(gi, gj, self.Nz_-1)
  #       tmpList[17]=self.gigjgkToGlobalID(gi, gj, 0) if pBc else -1
  #     if gk<self.Nz_-3:
  #       tmpList[5] =self.gigjgkToGlobalID(gi, gj, gk+1)
  #       tmpList[11]=self.gigjgkToGlobalID(gi, gj, gk+2)
  #       tmpList[17]=self.gigjgkToGlobalID(gi, gj, gk+3)

  #     # store currrent neighboring list
  #     self.G_[iPt] = tmpList
