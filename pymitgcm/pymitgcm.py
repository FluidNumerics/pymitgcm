#!/usr/bin/python3

import subprocess
import shlex
import datetime
import numpy as np
import MITgcmutils as mitgcm
import glob
from scipy import interpolate
from scipy.signal import convolve2d
from pyevtk.hl import gridToVTK, pointsToVTK

class pymitgcm:
  def __init__(self,directory='null',deltaTclock=1.0,iterate=0,dateTimeStart=[2000,1,1,0,0,0],loadState=True,loadGrid=True):

    self.deltaTclock = deltaTclock
    self.dateTimeStart = datetime.datetime(dateTimeStart[0],
                                           dateTimeStart[1],
                                           dateTimeStart[2],
                                           dateTimeStart[3],
                                           dateTimeStart[4],
                                           dateTimeStart[5])
    # Initialize fields
    self.iterate = iterate
    self.setDateTimeByIterate(iterate)
    self.directory = directory
    self.nz = 0
    self.ny = 0
    self.nx = 0
    self.xc = np.zeros(1) 
    self.yc = np.zeros(1)
    self.xg = np.zeros(1)
    self.yg = np.zeros(1)
    self.rf = np.zeros(1)
    self.rc = np.zeros(1)
    self.dxc = np.zeros(1)
    self.dyc = np.zeros(1)
    self.dxg = np.zeros(1)
    self.dyg = np.zeros(1)
    self.drf = np.zeros(1)
    self.drc = np.zeros(1)
    self.hfacc = np.zeros(1)
    self.hfacs = np.zeros(1)
    self.hfacw = np.zeros(1)
    self.u = np.zeros(1) 
    self.v = np.zeros(1)
    self.w = np.zeros(1)
    self.temperature = np.zeros(1)
    self.salinity = np.zeros(1)
    self.eta = np.zeros(1)
    self.atmfield = np.zeros(1)
#    self.pressure = mitgcm.mds.rdmds(self.directory+'/PH',iterate)

    if loadGrid :
      print('Loading Grid')
      self.gridload( )
      self.nz, self.ny, self.nx = np.shape(self.hfacc)
      print( "%d %d %d" % (self.nx, self.ny, self.nz) )

    
    # Load the initial iterate
    if loadState :
      self.stateload(iterate)
      self.nz, self.ny, self.nx = np.shape(self.temperature)
  
  def gridload(self):
      self.xc = mitgcm.mds.rdmds(self.directory+'/XC')
      self.yc = mitgcm.mds.rdmds(self.directory+'/YC')
      self.xg = mitgcm.mds.rdmds(self.directory+'/XG')
      self.yg = mitgcm.mds.rdmds(self.directory+'/YG')
      self.rf = mitgcm.mds.rdmds(self.directory+'/RF')
      self.rc = mitgcm.mds.rdmds(self.directory+'/RC')
      self.dxc = mitgcm.mds.rdmds(self.directory+'/DXC')
      self.dyc = mitgcm.mds.rdmds(self.directory+'/DYC')
      self.dxg = mitgcm.mds.rdmds(self.directory+'/DXG')
      self.dyg = mitgcm.mds.rdmds(self.directory+'/DYG')
      self.drf = mitgcm.mds.rdmds(self.directory+'/DRF')
      self.drc = mitgcm.mds.rdmds(self.directory+'/DRC')
      self.hfacc = mitgcm.mds.rdmds(self.directory+'/hFacC')
      self.hfacs = mitgcm.mds.rdmds(self.directory+'/hFacS')
      self.hfacw = mitgcm.mds.rdmds(self.directory+'/hFacW')

  def stateload(self, iterate=0):
  
    print('Loading Model State at iterate = %d'%iterate)
    self.setDateTimeByIterate(iterate)
    self.u = mitgcm.mds.rdmds(self.directory+'/U',iterate)
    self.v = mitgcm.mds.rdmds(self.directory+'/V',iterate)
    self.w = mitgcm.mds.rdmds(self.directory+'/W',iterate)
    self.temperature = mitgcm.mds.rdmds(self.directory+'/T',iterate)
    self.salinity = mitgcm.mds.rdmds(self.directory+'/S',iterate)
    self.eta = mitgcm.mds.rdmds(self.directory+'/Eta',iterate)
#    self.pressure = mitgcm.mds.rdmds(self.directory+'/PH',iterate)

  #END stateload

  def print_grid_size(self):

    print( ' MITgcm Grid Size ' )
    print( '==================' )
    print( 'Nx = '+str(self.nx) )
    print( 'Ny = '+str(self.ny) )
    print( 'Nz = '+str(self.nz) )

  #END print_grid_size

  def print_field_statistics(self):

    print('    MITgcm State Statistics    ')
    print('===============================')
    print('Min(Temp) : %d, %.3f' % (self.iterate, np.min(self.temperature)))
    print('Max(Temp) : %d, %.3f' % (self.iterate, np.max(self.temperature)))
    print('Min(Salt) : %d, %.3f' % (self.iterate, np.min(self.salinity)))
    print('Max(Salt) : %d, %.3f' % (self.iterate, np.max(self.salinity)))
    print('Min(U)    : %d, %.3f' % (self.iterate, np.min(self.u)))
    print('Max(U)    : %d, %.3f' % (self.iterate, np.max(self.u)))
    print('Min(V)    : %d, %.3f' % (self.iterate, np.min(self.v)))
    print('Max(V)    : %d, %.3f' % (self.iterate, np.max(self.v)))
    print('Min(W)    : %d, %.3f' % (self.iterate, np.min(self.w)))
    print('Max(W)    : %d, %.3f' % (self.iterate, np.max(self.w)))
    print('Min(Eta)  : %d, %.3f' % (self.iterate, np.min(self.eta)))
    print('Max(Eta)  : %d, %.3f' % (self.iterate, np.max(self.eta)))
    print('===============================')

  #END print_field_statistics

  def write_state_binary(self,path):

    ext = str(self.iterate).zfill(10)+'.bin'

    a = np.float32(self.temperature)
    output_file = open(path+'/temperature.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(self.salinity)
    output_file = open(path+'/salinity.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(self.u)
    output_file = open(path+'/u.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(self.v)
    output_file = open(path+'/v.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(self.w)
    output_file = open(path+'/w.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(self.eta)
    output_file = open(path+'/eta.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()


  def write_boundary_state_binary(self,path):

    ext = str(self.iterate).zfill(10)+'.bin'

    # south
    a = np.float32(np.squeeze(self.temperature[:,0,:]))
    output_file = open(path+'/temperature.south.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(np.squeeze(self.salinity[:,0,:]))
    output_file = open(path+'/salinity.south.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(np.squeeze(self.u[:,0,:]))
    output_file = open(path+'/u.south.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(np.squeeze(self.v[:,0,:]))
    output_file = open(path+'/v.south.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(np.squeeze(self.w[:,0,:]))
    output_file = open(path+'/w.south.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(np.squeeze(self.eta[0,:]))
    output_file = open(path+'/eta.south.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    # north
    a = np.float32(np.squeeze(self.temperature[:,-1,:]))
    output_file = open(path+'/temperature.north.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(np.squeeze(self.salinity[:,-1,:]))
    output_file = open(path+'/salinity.north.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(np.squeeze(self.u[:,-1,:]))
    output_file = open(path+'/u.north.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(np.squeeze(self.v[:,-1,:]))
    output_file = open(path+'/v.north.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(np.squeeze(self.w[:,-1,:]))
    output_file = open(path+'/w.north.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(np.squeeze(self.eta[-1,:]))
    output_file = open(path+'/eta.north.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    # west
    a = np.float32(np.squeeze(self.temperature[:,:,0]))
    output_file = open(path+'/temperature.west.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(np.squeeze(self.salinity[:,:,0]))
    output_file = open(path+'/salinity.west.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(np.squeeze(self.u[:,:,0]))
    output_file = open(path+'/u.west.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(np.squeeze(self.v[:,:,0]))
    output_file = open(path+'/v.west.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(np.squeeze(self.w[:,:,0]))
    output_file = open(path+'/w.west.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(np.squeeze(self.eta[:,0]))
    output_file = open(path+'/eta.west.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    # east
    a = np.float32(np.squeeze(self.temperature[:,:,-1]))
    output_file = open(path+'/temperature.east.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(np.squeeze(self.salinity[:,:,-1]))
    output_file = open(path+'/salinity.east.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(np.squeeze(self.u[:,:,-1]))
    output_file = open(path+'/u.east.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(np.squeeze(self.v[:,:,-1]))
    output_file = open(path+'/v.east.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(np.squeeze(self.w[:,:,-1]))
    output_file = open(path+'/w.east.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()

    a = np.float32(np.squeeze(self.eta[:,-1]))
    output_file = open(path+'/eta.east.'+ext, 'wb')
    a.tofile(output_file)
    output_file.close()


  def subsample(self,north,south,east,west):

    newself = pymitgcm(loadState=False,loadGrid=False)
    newself.deltaTclock = self.deltaTclock
    newself.dateTimeStart = self.dateTimeStart 
    newself.directory = self.directory
    newself.iterate = self.iterate

    # Find the indices nearest to the subsample domain boundaries
    x = np.squeeze( self.xc[0,:])
    y = np.squeeze( self.yc[:,0])
 
    i1 = np.argmin( np.abs(x-west) )
    i2 = np.argmin( np.abs(x-east) )
    j1 = np.argmin( np.abs(y-south) )
    j2 = np.argmin( np.abs(y-north) )
    
    newself.nx = (i2-i1)
    newself.ny = (j2-j1)
    newself.nz = self.nz

    newself.temperature = self.temperature[:,j1:j2,i1:i2]
    newself.salinity = self.salinity[:,j1:j2,i1:i2]
    newself.u = self.u[:,j1:j2,i1:i2]
    newself.v = self.v[:,j1:j2,i1:i2]
    newself.w = self.w[:,j1:j2,i1:i2]
    newself.eta = self.eta[j1:j2,i1:i2]

    if np.shape(self.atmfield) == np.shape(self.xc) :
      newself.atmfield = self.atmfield[j1:j2,i1:i2]

    return newself
    

  def refine(self,factor):

    newself = pymitgcm(loadState=False,loadGrid=False)
    newself.deltaTclock = self.deltaTclock
    newself.dateTimeStart = self.dateTimeStart 
    newself.directory = self.directory
    newself.iterate = self.iterate
    newself.nx = self.nx*factor
    newself.ny = self.ny*factor
    newself.nz = self.nz

    # Create a computational grid  [0,1]x[0,1] with nx*ny points (source grid)
    xi1 = np.linspace(0.0,1.0,num=self.nx)
    xi2 = np.linspace(0.0,1.0,num=self.ny)

    # Create a computational grid  [0,1]x[0,1] with (factor*nx)*(factor*ny) points (target grid)
    xi1p = np.linspace(0.0,1.0,num=newself.nx)
    xi2p = np.linspace(0.0,1.0,num=newself.ny)

    # Eta (Free surface height)
    var = np.squeeze(self.eta)
    mask = self.hfacc[0,:,:]*0
    mask[ self.hfacc[0,:,:] == 0.0 ] = 1.0

    for j in range(self.ny):
      for i in range(self.nx-1,0,-1):
        if mask[j,i] == 1.0: 
          var[j,i] = vlast
        else:
          vlast = var[j,i]
    f = interpolate.interp2d(xi1,xi2,var,kind='linear')
    newself.eta = f(xi1p,xi2p)
 
    newself.temperature = np.zeros([newself.nz,newself.ny,newself.nx])
    newself.salinity = np.zeros([newself.nz,newself.ny,newself.nx])
    newself.u = np.zeros([newself.nz,newself.ny,newself.nx])
    newself.v = np.zeros([newself.nz,newself.ny,newself.nx])
    newself.w = np.zeros([newself.nz,newself.ny,newself.nx])
        
    f = interpolate.interp2d(xi1,xi2,self.xc,kind='linear')
    newself.xc = f(xi1p,xi2p)

    f = interpolate.interp2d(xi1,xi2,self.yc,kind='linear')
    newself.yc = f(xi1p,xi2p)

    if np.shape(self.atmfield) == np.shape(self.xc) :
      f = interpolate.interp2d(xi1,xi2,self.atmfield,kind='linear')
      newself.atmfield = f(xi1p,xi2p)
    
    for k in range(newself.nz):

      var = np.squeeze(self.temperature[k,:,:])
      mask = self.hfacc[k,:,:]*0
      mask[ self.hfacc[k,:,:] == 0.0 ] = 1.0

      for j in range(self.ny):
        for i in range(self.nx-1,0,-1):
          if mask[j,i] == 1.0: 
            var[j,i] = tlast
          else:
            tlast = var[j,i]
        
      f = interpolate.interp2d(xi1,xi2,var,kind='linear')
      newself.temperature[k,:,:] = f(xi1p,xi2p)

      var = np.squeeze(self.salinity[k,:,:])
      mask = self.hfacc[k,:,:]*0
      mask[ self.hfacc[k,:,:] == 0.0 ] = 1.0

      for j in range(self.ny):
        for i in range(self.nx-1,0,-1):
          if mask[j,i] == 1.0: 
            var[j,i] = tlast
          else:
            tlast = var[j,i]

      f = interpolate.interp2d(xi1,xi2,var,kind='linear')
      newself.salinity[k,:,:] = f(xi1p,xi2p)

      var = np.squeeze(self.u[k,:,:])
      mask = self.hfacw[k,:,:]*0
      mask[ self.hfacw[k,:,:] == 0.0 ] = 1.0

      for j in range(self.ny):
        for i in range(self.nx-1,0,-1):
          if mask[j,i] == 1.0: 
            var[j,i] = tlast
          else:
            tlast = var[j,i]

      f = interpolate.interp2d(xi1,xi2,var,kind='linear')
      newself.u[k,:,:] = f(xi1p,xi2p)

      var = np.squeeze(self.v[k,:,:])
      mask = self.hfacs[k,:,:]*0
      mask[ self.hfacs[k,:,:] == 0.0 ] = 1.0

      for j in range(self.ny):
        for i in range(self.nx-1,0,-1):
          if mask[j,i] == 1.0: 
            var[j,i] = tlast
          else:
            tlast = var[j,i]

      f = interpolate.interp2d(xi1,xi2,var,kind='linear')
      newself.v[k,:,:] = f(xi1p,xi2p)

      var = np.squeeze(self.w[k,:,:])
      mask = self.hfacc[k,:,:]*0
      mask[ self.hfacc[k,:,:] == 0.0 ] = 1.0

      for j in range(self.ny):
        for i in range(self.nx-1,0,-1):
          if mask[j,i] == 1.0: 
            var[j,i] = tlast
          else:
            tlast = var[j,i]

      f = interpolate.interp2d(xi1,xi2,var,kind='linear')
      newself.w[k,:,:] = f(xi1p,xi2p)
    
    return newself

  def setDateTimeByIterate(self,iterate):

    self.iterate = iterate
    tSec = float(self.iterate)*float(self.deltaTclock)
    self.dateTime = self.dateTimeStart + datetime.timedelta(seconds=tSec)

  #END setDateTimeByIterate
	
  def getIterateList(self):

    sFiles = glob.glob(self.directory+'/S.*.data')
    sFiles.sort()
    iterates = []
    for s in sFiles:
      iterates.append(int(s.split('.')[1]))

    return iterates

  #END getIterateList

  def monstats(self):

    stats = {'time_tsnumber':[],
             'time_secondsf':[],
             'dynstat_eta_max':[],
             'dynstat_eta_min':[],
             'dynstat_eta_mean':[],
             'dynstat_eta_sd':[],
             'dynstat_eta_del2':[],
             'dynstat_uvel_max':[],
             'dynstat_uvel_min':[],
             'dynstat_uvel_mean':[],
             'dynstat_uvel_sd':[],
             'dynstat_uvel_del2':[],
             'dynstat_vvel_max':[],
             'dynstat_vvel_min':[],
             'dynstat_vvel_mean':[],
             'dynstat_vvel_sd':[],
             'dynstat_vvel_del2':[],
             'dynstat_wvel_max':[],
             'dynstat_wvel_min':[],
             'dynstat_wvel_mean':[],
             'dynstat_wvel_sd':[],
             'dynstat_wvel_del2':[],
             'dynstat_theta_max':[],
             'dynstat_theta_min':[],
             'dynstat_theta_mean':[],
             'dynstat_theta_sd':[],
             'dynstat_theta_del2':[],
             'dynstat_salt_max':[],
             'dynstat_salt_min':[],
             'dynstat_salt_mean':[],
             'dynstat_salt_sd':[],
             'dynstat_salt_del2':[],
             'trAdv_CFL_u_max':[],
             'trAdv_CFL_v_max':[],
             'trAdv_CFL_w_max':[],
             'advcfl_uvel_max':[],
             'advcfl_vvel_max':[],
             'advcfl_wvel_max':[],
             'advcfl_W_hf_max':[],
             'pe_b_mean':[],
             'ke_max':[],
             'ke_mean':[],
             'ke_vol':[],
             'vort_r_min':[],
             'vort_r_max':[],
             'vort_a_mean':[],
             'vort_a_sd':[],
             'vort_p_mean':[],
             'vort_p_sd':[],
             'surfExpan_theta_mean':[],
             'surfExpan_salt_mean':[]}

    with open(self.directory+'/STDOUT.0000', 'r') as fp:
      for line in fp:
        for var in stats.keys():
          if var in line:
            stats[var].append(np.float(line.split('=')[-1].lstrip().rstrip()))

    stats['time_secondsf'] = np.unique(np.array(stats['time_secondsf']))

    output_payload = {'monitor_stats':[]}
    for var in stats.keys():
      for k in range(len(stats['time_secondsf'])):
        simDateTime = self.dateTimeStart + datetime.timedelta(seconds=stats['time_secondsf'][k])
        pl = {'simulation_datetime': simDateTime,
              'metric_name':var,
              'metric_value':stats[var][k]}
        output_payload['monitor_stats'].append(pl)

    return output_payload

  #END monStats

  def maxAdvFrequency(self):

    ucfl = np.zeros((self.nz,self.ny,self.nx))
    vcfl = np.zeros((self.nz,self.ny,self.nx))
    wcfl = np.zeros((self.nz,self.ny,self.nx))
    for k in range(self.nz):
      for j in range(self.ny):
        for i in range(self.nx):
          ucfl[k,j,i] = abs(self.u[k,j,i]/self.dxc[j,i])
          vcfl[k,j,i] = abs(self.v[k,j,i]/self.dyc[j,i])
          wcfl[k,j,i] = abs(self.w[k,j,i]/self.drf[k])

    return {'max_UAdvFreq':np.asscalar(ucfl.max()),
            'max_VAdvFreq':np.asscalar(vcfl.max()),
            'max_WAdvFreq':np.asscalar(wcfl.max())}

  #END maxAdvFrequency 

  def globalStats(self):

    return {'max_temperature': np.asscalar(self.temperature.max()),
            'max_salinity': np.asscalar(self.salinity.max()),
            'min_temperature': np.asscalar(self.temperature.min()),
            'min_salinity': np.asscalar(self.salinity.min())}

  #END globalStats
          

  def zBoxModel(self, zsplit):
#"""Splits the ocean into upper ocean [0,zsplit] and lower ocean (zsplit,bottom] and computes state averages for upper and lower ocean"""

    # Find the depth index for the upper ocean
    ku = np.zeros(2,dtype=np.int8)
    kl = np.zeros(2,dtype=np.int8)
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

    for ibox in range(2):
      for k in range(ku[ibox],kl[ibox]):
        for j in range(self.ny):
          for i in range(self.nx):
            vol[ibox] += self.dxc[j,i]*self.dyc[j,i]*self.drf[k]*self.hfacc[k,j,i]
            temp[ibox] += self.temperature[k,j,i]*self.dxc[j,i]*self.dyc[j,i]*self.drf[k]*self.hfacc[k,j,i]
            salt[ibox] += self.salinity[k,j,i]*self.dxc[j,i]*self.dyc[j,i]*self.drf[k]*self.hfacc[k,j,i]

      temp[ibox] = np.asscalar(temp[ibox]/vol[ibox])
      salt[ibox] = np.asscalar(salt[ibox]/vol[ibox])

    return {'zm_temperature_upper':temp[0],'zm_temperature_lower':temp[1],'zm_salinity_upper':salt[0],'zm_salinity_lower':salt[1],'zm_volume_upper':vol[0],'zm_volume_lower':vol[1]}

  #END zBoxModel

  def writeToVTK(self, outdir):

    fillValue=1000.0
    subprocess.check_call(shlex.split('mkdir -p {}'.format(outdir)))

    ncells = self.nx*self.ny*self.nz
    npoints = (self.nx+1)*(self.ny+1)*(self.nz+1) 

    x = np.zeros((self.nx + 1))
    y = np.zeros((self.ny + 1))
    z = np.zeros((self.nz + 1))

    # Append the east-most xg point xg
    for i in range(self.nx):
      x[i] = self.xg[0,i]
    x[self.nx] = self.xg[0,self.nx-1] + self.xg[0,self.nx-1] - self.xg[0,self.nx-2] 

    for j in range(self.ny):
      y[j] = self.yg[j,0]
    y[self.ny] = self.yg[self.ny-1,0] + self.yg[self.ny-1,0] - self.yg[self.ny-2,0] 

    for k in range(self.nz+1):
      z[k] = self.rf[k]/1000.0

    # Create the topographic mask vtk output
    wetdryMask = np.transpose(self.hfacc, (2, 1, 0))
    wetdryMask[wetdryMask<1.0] = 0.0

    gridToVTK(outdir+"/mask",x,y,z, cellData = {"wetdry_mask": wetdryMask})

    temperature = np.transpose(self.temperature,(2,1,0))
    salinity = np.transpose(self.salinity,(2,1,0))
    u = np.transpose(self.u,(2,1,0))
    v = np.transpose(self.v,(2,1,0))
    w = np.transpose(self.w,(2,1,0))
    eta = np.transpose(self.eta,(1,0))

    temperature[wetdryMask==0.0] = fillValue
    salinity[wetdryMask==0.0] = fillValue
    u[wetdryMask==0.0] = fillValue
    v[wetdryMask==0.0] = fillValue
    w[wetdryMask==0.0] = fillValue
    eta[wetdryMask[:,:,0]==0.0] = fillValue
    eta = eta[:,:,np.newaxis]

    vtkFile = 'state.%s'%str(self.iterate).zfill(10)
    gridToVTK(outdir+"/"+vtkFile,x,y,z, cellData = {"temperature": temperature,
                                                    "salinity": salinity,
                                                    "u": u,
                                                    "v": v,
                                                    "w": w})

    vtkFile = 'eta.%s'%str(self.iterate).zfill(10)
    gridToVTK(outdir+"/"+vtkFile,x,y,z[0:1],cellData = {"eta":eta})
