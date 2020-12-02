#!/usr/bin/python3

import datetime
import numpy
import MITgcmutils as mitgcm

class mitgcm:
  def __init__(self,directory,deltaTclock,iterate,dateTimeStart=[2000,1,1,0,0,0])

    self.iterateState = iterate
    self.dateTimeStart = datetime.datetime(dateTimeStart[0],
                                           dateTimeStart[1],
                                           dateTimeStart[2],
                                           dateTimeStart[3],
                                           dateTimeStart[4],
                                           dateTimeStart[5])
    self.iterate = iterate
    self.dateTime = self.simStart
    self.deltaTclock = deltaTclock
    self.directory = directory

    # Load Grid information
    self.xc = mitgcm.mds.rdmds(directory+'/XC')
    self.yc = mitgcm.mds.rdmds(directory+'/YC')
    self.xg = mitgcm.mds.rdmds(directory+'/XG')
    self.yg = mitgcm.mds.rdmds(directory+'/XG')
    self.rf = mitgcm.mds.rdmds(directory+'/RF')
    self.rc = mitgcm.mds.rdmds(directory+'/RC')
    self.dxf = mitgcm.mds.rdmds(directory+'/DXF')
    self.dyf = mitgcm.mds.rdmds(directory+'/DYF')
    self.dxg = mitgcm.mds.rdmds(directory+'/DXG')
    self.dyg = mitgcm.mds.rdmds(directory+'/DYG')
    self.hfacc = mitgcm.mds.rdmds(directory+'/hFacC')
    self.hfacS = mitgcm.mds.rdmds(directory+'/hFacS')
    self.hfacW = mitgcm.mds.rdmds(directory+'/hFacW')

    self.nx, self.ny, self.nz = np.shape(self.hFacC)
    print( "%d %d %d" % (self.nx, self.ny, self.nz) )
    
    # Load the initial iterate
    self.stateload()

  
  def stateload(self, iterate=0):
  
    self.iterate = iterate
    self.u = mitgcm.mds.rdmds(self.directory+'/U',iterate)
    self.v = mitgcm.mds.rdmds(self.directory+'/V',iterate)
    self.w = mitgcm.mds.rdmds(self.directory+'/W',iterate)
    self.temperature = mitgcm.mds.rdmds(self.directory+'/T',iterate)
    self.salinity = mitgcm.mds.rdmds(self.directory+'/S',iterate)
    self.eta = mitgcm.mds.rdmds(self.directory+'/Eta',iterate)
    self.pressure = mitgcm.mds.rdmds(self.directory+'/PH',iterate)

  #END stateload

  def setDateTimeByIterate(self,iterate):

    self.iterate = iterate
    tSec = float(self.iterate-self.iterateStart)*float(self.deltaTclock)
    self.dateTime = self.dateTimeStart + datetime.timedelta(seconds=tSec)

  #END setDateTimeByIterate
	
  def getIterateList(self):

    sFiles = glob.glob(self.directory+'/S.*.data')
    iterates = []
    for s in sFiles:
      iterates.append(int(s.split('.')[1]))

    return iterates

  #END getIterateList

  def zBoxModel(self, zsplit)
  """Splits the ocean into upper ocean [0,zsplit] and lower ocean (zsplit,bottom] and computes state averages for upper and lower ocean"""

    # Find the depth index for the upper ocean
    ku = np.zeros(2)
    kl = np.zeros(2)
    for i in range(self.nz):
      if self.rc[i] <= zsplit :
        ksplit = i
        break
    kl[0] = ksplit
    ku[1] = ksplit+1
    kl[1] = self.nz

    temp = np.zeros(2)
    salt = np.zeros(2)
    vol = np.zeros(2)

    for ibox in range(2)
      for k in range(ku[ibox],kl[ibox]):
        for j in range(ny):
          for i in range(nx):
            vol[ibox] += self.dxc[j,i]*self.dyc[j,i]*self.drf[k]*self.hfacc[k,j,i]

      for k in range(ku[ibox],kl[ibox]):
        for j in range(ny):
          for i in range(nx):
            temp[ibox] += self.temperature[k,j,i]*self.dxc[j,i]*dy[j,i]*self.drf[k]*self.hfacc[k,j,i]
            salt[ibox] += self.salinity[k,j,i]*self.dxc[j,i]*dyc[j,i]*self.drf[k]*self.hfacc[k,j,i]

      temp[ibox] = np.asscalar(temp[ibox]/vol)
      salt[ibox] = np.asscalar(salt[ibox]/vol)

    return {'temperature':temp,'salinity':salt,'volume':vol,'kl':kl,'ku',ku}

  #END zBoxModel
