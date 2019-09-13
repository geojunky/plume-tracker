#!/usr/bin/python
'''
 * =====================================================================================
 *
 * Copyright (C) 2014 Rakib Hassan (rakib.hassan@sydney.edu.au)
 *
 * This program is free software; you can redistribute it and/or modify it under 
 * the terms of the GNU General Public License as published by the Free Software 
 * Foundation; either version 2 of the License, or (at your m_option) any later 
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple 
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * ===================================================================================== 
'''

import sys
from optparse import OptionParser
import h5py
import math
from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy import spatial
from scipy.interpolate import griddata
from collections import defaultdict

import numpy as np
import os
import glob

from matplotlib.tri import Triangulation, UniformTriRefiner, LinearTriInterpolator
import matplotlib.cm as cm
import matplotlib.patches as mpatches
from scipy.interpolate import LinearNDInterpolator
from pyrr import quaternion, matrix44, vector3, vector4

import Pycluster as pc

#spatial
from matplotlib.collections import PatchCollection
from descartes import PolygonPatch
import fiona
from shapely.geometry import Polygon, MultiPolygon, shape
#========================================================================

msg=\
"""

Plots results from cluster-analysis from plumeTrackFast

Usage 1: plotPlumeAscentRate.py -b <base file-name> -d <reconstruction-directory> -n <neighbour file-name> -r <radial-velocity minimum> -s <start-age> -o <output file-name>

"""

RADIUS = 6371e3
THERM_DIFF = 1e-6
MIN_PLUME_DIST = 500e3
MAX_PLUME_OFFSET = 1000e3
PLUME_PLOT_RADIUS = 400e3
VSCALE = 1e-4 # convert from m/Myr to cm/yr
DENSITY = 4e3 #kg/m3
THERMAL_EXPANSIVITY = 2.875e-5
DELTA_T = 3500

SYMBOL_SCALE = 20

def rtp2xyz(r, theta, phi):
    rst = r * math.sin(theta)
    xout = [0]*3
    xout[0] = rst * math.cos(phi)	    # x
    xout[1] = rst * math.sin(phi) 	    # y
    xout[2] = r * math.cos(theta)       # z

    return xout
#end

def xyz2rtp(x, y, z):
    tmp1 = x*x + y*y
    tmp2 = tmp1 + z*z
    rout = [0]*3
    rout[0] = math.sqrt(tmp2)
    rout[1] = math.atan2(math.sqrt(tmp1),z)
    rout[2] = math.atan2(y,x)
    return rout
#end

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.asin(math.sqrt(a))

    dist = RADIUS * c
    return dist
#end func

class Plume:
    def rotateAboutOrigin(self, xyz, t, p):
        rp = p - math.pi/2
        rt = -t

        q1 = quaternion.create_from_z_rotation(rp)
        q2 = quaternion.create_from_x_rotation(rt)
        v = vector3.create(xyz[0], xyz[1], xyz[2])

        v = quaternion.apply_to_vector(q1, v)
        v = quaternion.apply_to_vector(q2, v)

        return v
    #end func

    # data must have the following columns: [lon, lat, temperature, vr, concentration]
    def __init__(self, options, data, meanTemp, suppFigName):
        self.options = options
        self.data = data
        self.meanTemp = meanTemp
        self.suppFigName = suppFigName

        self.mlon = np.mean(self.data[:,0], axis=0)
        self.mlat = np.mean(self.data[:,1], axis=0)

        self.nodes = []
        for i in range(np.shape(self.data)[0]):
            self.nodes.append(rtp2xyz(RADIUS, self.data[i,1]*math.pi/180, self.data[i,0]*math.pi/180))
        #end for
        self.nodes = np.array(self.nodes)

        mx = np.mean(self.nodes[:,0], axis=0)
        my = np.mean(self.nodes[:,1], axis=0)
        mz = np.mean(self.nodes[:,2], axis=0)

        rtp = xyz2rtp(mx, my, mz)
        if(rtp[2]<0): rtp[2] = rtp[2] + 2*math.pi

        # rotate nodes to the pole
        for i in range(np.shape(self.data)[0]):
            res = self.rotateAboutOrigin(self.nodes[i,:], rtp[1], rtp[2])
            self.nodes[i,:] = np.squeeze(res) / 1e3 # scale to Km
        #end for

        self.t = Triangulation(self.nodes[:,0], self.nodes[:,1])
        self.interpTempObject = LinearNDInterpolator(self.nodes[:,:2], self.data[:,2])
        self.interpVrObject = LinearNDInterpolator(self.nodes[:,:2], self.data[:,3])
        self.interpConcObject = LinearNDInterpolator(self.nodes[:,:2], self.data[:,4])

        self.conduitFound = False
        self.LocateConduits()
    #end func

    def LocateConduits(self):
        plt.figure()
        plt.gca().set_aspect('equal')

        plt.xlim(-PLUME_PLOT_RADIUS/1e3, PLUME_PLOT_RADIUS/1e3)
        plt.ylim(-PLUME_PLOT_RADIUS/1e3, PLUME_PLOT_RADIUS/1e3)

        localMeanVr = np.mean(self.data[:,3]*VSCALE)
        localMinVr  = np.min(self.data[:,3]*VSCALE)
        localMaxVr  = np.max(self.data[:,3]*VSCALE)
        localArgMaxVr  = np.argmax(self.data[:,3]*VSCALE)

        localVrFast=[]
        for i in range(np.shape(self.nodes)[0]):
            if(self.data[i,3]>localMeanVr): localVrFast.append(self.data[i,3])
        #end for
        localVrFast = np.array(localVrFast)
        localVrFastStd  = np.std(localVrFast*VSCALE)

        levels = [localMaxVr-2*localVrFastStd]
        cmap = cm.get_cmap(name='jet', lut=None)

        plt.tricontourf(self.t, self.data[:,3]*VSCALE-localMeanVr, cmap=cmap)
        plt.colorbar()

        cs = plt.tricontour(self.t, self.data[:,3]*VSCALE-localMeanVr, levels=levels, linewidths=[3], colors=['k'], linestyles=['dashed'])

        if(len(cs.collections)):
            for collection in cs.collections:
                self.plumeConduitNodes = []
                broken = False
                for path in collection.get_paths():
                    polygon = path.to_polygons()
                    vertices = np.squeeze(np.array(polygon))

                    if(path.contains_point([self.nodes[localArgMaxVr,0],self.nodes[localArgMaxVr,1]])):

                        #patch = mpatches.PathPatch(path, edgecolor='g', fill=False, lw=3, linestyle='dashed')
                        #plt.gca().add_patch(patch)

                        # Add plumeConduit nodes
                        for v in vertices: self.plumeConduitNodes.append(v)
                        interiorPointsX = np.random.uniform(low=-PLUME_PLOT_RADIUS/1e3, high=PLUME_PLOT_RADIUS/1e3, size=1000)
                        interiorPointsY = np.random.uniform(low=-PLUME_PLOT_RADIUS/1e3, high=PLUME_PLOT_RADIUS/1e3, size=1000)

                        for i in range(np.shape(interiorPointsX)[0]):
                            if(path.contains_point([interiorPointsX[i], interiorPointsY[i]])):
                                self.plumeConduitNodes.append([interiorPointsX[i], interiorPointsY[i]])
                            #end if
                        #end for
                        self.plumeConduitNodes = np.array(self.plumeConduitNodes)
                        self.pct = Triangulation(self.plumeConduitNodes[:,0], self.plumeConduitNodes[:,1])
                        plt.triplot(self.pct, lw=0.4, color='green')

                        #print vertices
                        self.conduitFound = True
                        broken = True
                        break
                    #end if

                #end for
                if(broken): break
            #end for
        #end if

        plt.clabel(cs, inline=1, fontsize=8)
        plt.gca().set_xlabel('x / km')
        plt.gca().set_ylabel('y / km')
        plt.gca().set_title('$V_r / (cm/yr)$')

        if(self.conduitFound==False):
            # Naive approach where conduit radii is assumed fixed. This acts as a fallback
            # usually when conduit velocities are too high at the time of eruption
            maxDist = ((PLUME_PLOT_RADIUS/2)**2/1e6)
            self.plumeConduitNodes = []
            for i in range(np.shape(self.nodes)[0]):
                if((self.nodes[i,0]**2 + self.nodes[i,1]**2) < maxDist):
                    self.plumeConduitNodes.append([self.nodes[i,0], self.nodes[i,1]])
                #end if
            #end for
            self.plumeConduitNodes = np.array(self.plumeConduitNodes)
            self.pct = Triangulation(self.plumeConduitNodes[:,0], self.plumeConduitNodes[:,1])
            plt.triplot(self.pct, lw=0.4, color='green')
            self.conduitFound = True

            print '\tWarning: failed to locate conduit - resorting to fallback option..'
        #end if

        plt.savefig(self.suppFigName)
        plt.close()
    #end func

    def ComputeBuoyancyFlux(self):
        if(self.conduitFound == False): return {'flux':0, 'area':0, 'avgConduitTemp':0, 'avgConduitConcentration':0}
                
        flux = 0
        totalArea = 0
        avgConduitTemp = 0
        avgConduitConcentration = 0
        for tri in self.pct.triangles:
            mat = np.matrix( [self.plumeConduitNodes[tri[0],0:2] - self.plumeConduitNodes[tri[2],0:2], \
                              self.plumeConduitNodes[tri[1],0:2] - self.plumeConduitNodes[tri[2],0:2]] )

            area = math.fabs(0.5*np.linalg.det(mat))

            totalArea = totalArea + area
            center = self.plumeConduitNodes[tri,0:2].mean(axis=0)


            flux = flux + area * \
                          DENSITY * THERMAL_EXPANSIVITY * (self.interpTempObject(center[0], center[1]) - self.meanTemp) * DELTA_T *\
                          self.interpVrObject(center[0], center[1]) * \
                          THERM_DIFF/RADIUS / 1e3 # Mg/s

            avgConduitTemp += self.interpTempObject(center[0], center[1])
            avgConduitConcentration += self.interpConcObject(center[0], center[1])
        #end for

        if(len(self.pct.triangles)>0): avgConduitTemp /= len(self.pct.triangles)
        if(len(self.pct.triangles)>0): avgConduitConcentration /= len(self.pct.triangles)

        # Note: spatial coordinates were converted to kms in __init__
        # Areas therefore need to be scaled by 1e6
        flux *= 1e6
        return {'flux':flux, 'area':totalArea, 'avgConduitTemp':avgConduitTemp,
                'avgConduitConcentration':avgConduitConcentration} #area:km2
    #end func
#end class

#================================ Local functions
def processPlumes(options):
    eps=1e-5

    #=============================== Local functions
    def checkConduitContinuity(nnodes):
        maxVrShallow = -1e20
        maxVrDeep = -1e20
        selectedNodeShallow = -1
        selectedNodeDeep = -1
        for n in nnodes:
            if((cdata[n,5]>maxVrShallow) and (cdata[n,9]==0)):
                maxVrShallow = cdata[n,5]
                selectedNodeShallow = n
            #end if

            if((cdata[n,11]>maxVrDeep) and (cdata[n,15]==0)):
                maxVrDeep = cdata[n,11]
                selectedNodeDeep = n
            #end if
        #end for

        # Check for continuity of identified plume-region at depth
        if(cdata[selectedNodeShallow, 9] == cdata[selectedNodeDeep, 15]):
            return selectedNodeShallow, selectedNodeDeep, True,
        else:
            return selectedNodeShallow, selectedNodeDeep, False
        #end if
    #end func

    def isPlumeEruption(plumeIds, theta, phi):
        mindist = 1e20
        for p in plumeIds:
            dist = haversine(cdata[p,2], cdata[p,1], phi, theta)

            if(dist < mindist): mindist = dist
        #end for
        if(mindist > MIN_PLUME_DIST): return True
        else: return False
    #end func

    def shouldAdd(plist, pid):
        nxyz = rtp2xyz(cdata[pid,0]*RADIUS, \
                       (90-cdata[pid,1])*math.pi/180, \
                       cdata[pid,2]*math.pi/180)

        for epid in plist:
            exyz = rtp2xyz(cdata[epid,0]*RADIUS, \
                           (90-cdata[epid,1])*math.pi/180, \
                           cdata[epid,2]*math.pi/180)
            dist = (nxyz[0] - exyz[0])**2 + \
                   (nxyz[1] - exyz[1])**2 + \
                   (nxyz[2] - exyz[2])**2
            dist = math.sqrt(dist)
            if(dist < MIN_PLUME_DIST): return 0
        #end for

        return 1
    #end func

    def getNeighbours(nbFileName):
        nf = open(nbFileName)
        neighbours = defaultdict(list)

        for line in nf:
            vals = line.split(' ')
            for i in range(len(vals)):
                if(i>0):neighbours[int(vals[0])].append(int(vals[i]))
            #end for
        #end for
        nf.close()

        return neighbours
    #end func

    def getNeighbouringNodes(tree, xyz, radius):
        l = tree.query_ball_point(xyz, radius)
        return l
    #end func

    #=============================== Begin processing
    def key_func(x): return int(x.split('.')[-2])
    fileNames = glob.glob('%s*.txt'%(options.baseName))
    fileNames = sorted(fileNames, key=key_func)

    minNablaV = 1e30
    maxNablaV = -1e30
    cdataList = []
    for fileName in fileNames:

        ''' Data-format=========================================
        r, theta, phi, time, magGradVr, v_r, T, Eta, data cid magGradVrVD, v_rVD, TVD, EtaVD, dataVD cidVD tracComp
        0   1      2    3        4       5   6   7     8   9       10        11   12    13      14     15     16  */
        ====================================================='''

        # Load clusteredData
        cdata = np.loadtxt(fileName)
        cdataList.append(cdata)

        if(np.min(cdata[:,4]) < minNablaV): minNablaV = np.min(cdata[:,4])
        if(np.max(cdata[:,4]) > maxNablaV): maxNablaV = np.max(cdata[:,4])
    #end for

    nodes=[]
    tree = None
    plumeIdsList = [] # List of identified plumes for each time-slice
    outFile = open(options.outFileName, 'w')
    for counter,fileName in enumerate(fileNames):

        print 'Processing: %s..'%(fileName)

        currTime = int(fileName.split('.')[-2])
        dirName = os.path.dirname(options.outFileName)
        if(len(dirName)==0):dirName = './'
        baseName = os.path.basename(fileName)
        suppDirName = '%s/supp/'%(dirName)

        # create supplimentary directory if it doesn't exist
        if (not os.path.exists(suppDirName)):
            os.makedirs(suppDirName)
        #end if

        cdata = cdataList[counter]
        if(counter==0):
            for i in range(np.shape(cdata)[0]):
                xyz=rtp2xyz(RADIUS,math.pi/2-cdata[i,1]*math.pi/180,cdata[i,2]*math.pi/180)
                nodes.append(xyz)
            #end for
            tree = spatial.cKDTree(nodes, 16)
        #end if

        dataOut = np.zeros((361,181))
        xOut = np.zeros((361,181))
        yOut = np.zeros((361,181))
        for i in range(361):
            for j in range(181):
                theta = float(j)
                phi = float(i)

                xOut[i,j] = phi*math.pi/180
                yOut[i,j] = theta*math.pi/180
                xyz = rtp2xyz(RADIUS, yOut[i,j], xOut[i,j])

                (d, l) = tree.query(xyz)
                dataOut[i,j] = cdata[l, 4] * THERM_DIFF / RADIUS #magGradVrDepth

                #dataOut[i,j] = cdata[l, 9] #cidDepth
                #dataOut[i,j] = cdata[l, 15] #cidValidationDepth
            #end for
        #end for

        #=======================================Compute fluxes
        numNodes = np.shape(cdata)[0]
        plumeIds=[]
        meanTemp = np.mean(cdata[:,6], axis=0)

        meanBgTemp = 0.
        hitCount = 0.
        for i in range(numNodes):
            if(cdata[i,6] > meanTemp):
                meanBgTemp += cdata[i,6]
                hitCount += 1.
            #end if
        #end for
        meanBgTemp /= hitCount

        print 'mean: %f, meanBg: %f'%(meanTemp, meanBgTemp)

        for i in range(numNodes):
            if((cdata[i, 9]==0)):

                nn = np.array(getNeighbouringNodes(tree, nodes[i], MAX_PLUME_OFFSET))
                selectedNodeShallow, selectedNodeDeep, conduitIsContinuous = checkConduitContinuity(nn)

                if(conduitIsContinuous == False): continue


                # Check if plume-conduit meets the minimum radial velocity criteria
                if(cdata[selectedNodeShallow, 5]*VSCALE < options.radialVelocityMin): continue

                if(shouldAdd(plumeIds, selectedNodeShallow)):
                    plumeIds.append(selectedNodeShallow)
                #end if
            #end if
        #end for

        # Save plume-ids list
        plumeIdsList.append(plumeIds)

        flux = np.zeros(len(plumeIds))
        area = np.zeros(len(plumeIds))
        avgConduitTemp = np.zeros(len(plumeIds))
        avgConduitConcentration = np.zeros(len(plumeIds))

        for i,pid in enumerate(plumeIds):
            #print "plume %d: T(%d) = %f"%(i, pid, cdata[pid,6])

            suppFigName = '%s/%s.plume%d.png'%(suppDirName, baseName, i)
            nn = np.array(getNeighbouringNodes(tree, nodes[pid], PLUME_PLOT_RADIUS))
            p = Plume(options, np.array((cdata[nn,2], 90-cdata[nn,1], cdata[nn,6], cdata[nn, 5], cdata[nn, 16])).T, 
                      meanBgTemp, suppFigName)
            result = p.ComputeBuoyancyFlux()
            flux[i] = result['flux']
            area[i] = result['area']
            avgConduitTemp[i] = result['avgConduitTemp']
            avgConduitConcentration[i] = result['avgConduitConcentration']

            # Report plume=======
            mlatp = 90-p.mlat
            mlonp = p.mlon
            if(mlonp > 180): mlonp = mlonp-360
            print '\tplume(%d) loc: (%f %f)'%(i, mlatp, mlonp)

        #end for
        #================================== Done computing fluxes


        #================================== Prepare  data and axes
        dataOutFlipped = np.zeros(np.shape(dataOut))
        dataOutFlipped[0:181:,:] = dataOut[180:361,:]
        dataOutFlipped[181:361:,:] = dataOut[0:180,:]

        fig, gax = plt.subplots(nrows=1, ncols=1)
        fig.set_size_inches(11.69, 11.69)
        m = Basemap(ax=gax, projection='moll',lon_0=0,resolution='c')


        #=======================================Plot flux
        pxOut,pyOut = m(xOut*180/math.pi-180, 90-yOut*180/math.pi)
        cmap = cm.get_cmap(name='jet', lut=None)
        cbInfo = m.pcolormesh(pxOut, pyOut, dataOutFlipped, cmap=cmap)
        m.colorbar(cbInfo, "bottom", size="5%", pad="15%")

        legList=[]
        labels = []
        for i,pid in enumerate(plumeIds):
            #print "plume %d: T(%d) = %f"%(i, pid, cdata[pid,6])

            #coords = circle(m, cdata[pid,2], cdata[pid,1], 50*flux[i])
            #casa = SPolygon(coords)
            #cpatch = PolygonPatch(casa, fc='magenta', ec='#333333', lw=0, alpha=.25, zorder=2)
            #gax.add_patch(cpatch)

            if(flux[i]>0):
                plumeEruption = 0
                if(counter):
                    plumeEruption = int(isPlumeEruption(plumeIdsList[counter-1], cdata[pid,1], cdata[pid,2]))
                #end if

                #(time, theta, phi, T, Vr, Flux, Area, Eruption, meanBgTemp, avgConduitTemp, concentration)
                outFile.write('%f %f %f %f %f %f %f %f %f %f %f\n'%( currTime, cdata[pid,1], cdata[pid,2],
                                                      cdata[pid, 6], cdata[pid, 5], flux[i], area[i], plumeEruption,
                                                      meanBgTemp, avgConduitTemp[i], avgConduitConcentration[i] ))


                px,py = m(cdata[pid,2], cdata[pid,1])
                #m.scatter(px,py,flux[i]*SYMBOL_SCALE, marker='o',color='green', edgecolors='none', alpha=0.5)
                m.scatter(px,py,flux[i]*SYMBOL_SCALE, marker='o',color='none', edgecolors='white', lw=0.5, alpha=1)

                leg = m.scatter([],[],s=0, edgecolors='black', color='none')
                legList.append(leg)
                labels.append('(%d) : %0.1f'%(i, flux[i]))

                phil = cdata[pid,2]
                thetal = cdata[pid,1]
                pxl,pyl = m(phil, thetal)
                gax.text(pxl,pyl,'%d'%(i),color='white', alpha=0.5)

                #gax.text(px,py,'%.1f Mg/s'%(flux[i]),color='white')
                #gax.text(px,py,'%d: %.0f km2'%(i, area[i]),color='white')
                #gax.text(px,py,'%.3f'%(cdata[pid, 6]),color='white') # temperature
            else:
                print 'Negative flux: %f detected for plume %d..'%(flux[i], i)
            #end if
        #end for
        #=======================================Done plotting flux

        #==================================================Plot reconstruction
        #m.drawparallels(np.arange(-90.,210.,30.))
        #m.drawmeridians(np.arange(-180.,240.,60.))
        age = options.startAge-int(fileName.split('.')[-2])
        if(age<0):age = 0
        if 1:
            sfName = '%s/reconstructed_%d.00Ma.shp'%(options.reconsPath,age)
            sf = fiona.open(sfName)

            features = []
            for feature in sf:
                #print feature['geometry']['type']
                features.append(shape(feature['geometry']))
            #end for

            patches = []
            for feature in features:
                if(isinstance(feature, Polygon)):
                    polygon = feature
                    x, y = polygon.exterior.coords.xy
                    px, py = m(x, y)
                    ppolygon = Polygon(zip(px,py))
                    patches.append(PolygonPatch(ppolygon, fill=False, edgecolor=(0.5, 0.5, 0.5), alpha=0.5))
                else:
                    multiPolugon = feature

                    for polygon in multiPolugon:
                        x, y = polygon.exterior.coords.xy
                        px, py = m(x, y)
                        ppolygon = Polygon(zip(px,py))
                        patches.append(PolygonPatch(ppolygon, fill=False, edgecolor=(0.5, 0.5, 0.5), alpha=0.5))
                    #end for
            #end for
            gax.add_collection(PatchCollection(patches, match_original=True))
        #end if
        m.drawmapboundary()
        #=======================================Done plotting reconstruction

        leg = plt.legend(legList, labels, ncol=6, frameon=False, fontsize=12,
                          handlelength=1,
                          loc = 'lower center', borderaxespad = -13,
                          handletextpad=0.2, title='Flux / (Mg/s)', scatterpoints = 1,
                          labelspacing=0.2)

        gax.set_title('Plumes Identified at %d Ma'%(age))
        gax.set_xlabel('$\\Vert\\nabla V_r \\Vert [s^{-1}]$')
        
        plt.savefig('%s/%s.png'%(dirName, os.path.basename(fileName)), dpi=300)
        plt.close()
    #end for
    outFile.close()
#end function

def main():
    parser = OptionParser(msg)

    parser.add_option("-b", "--base-name", dest="baseName",
                      help="File name for clustered data", type="string")
    parser.add_option("-d", "--directory", dest="reconsPath",
                      help="Path to reconstruction files", type="string")
    parser.add_option("-s", "--start-age", dest="startAge",
                      help="Start age", type="int")    
    parser.add_option("-n", "--neighbours", dest="nbFileName",
                      help="File containing node-neighbours from spherical triangulation", type="string")
    parser.add_option("-r", "--radial-velocity-minimum (cm/yr)", dest="radialVelocityMin",
                      help="Minimum radial velocity that a plume-conduit must meet to be considered as such", type="float")
    parser.add_option("-o", "--output-file", dest="outFileName",
                      help="Output-file name", type="string")

    (options, args) = parser.parse_args()

    if (options.baseName is None): parser.error('File-name not provided')
    if (options.startAge is None): parser.error('Start age not provided')
    if (options.reconsPath is None): parser.error('Reconstruction-path not provided')
    if (options.radialVelocityMin is None): parser.error('Minimum radial velocity not provided')
    if (options.nbFileName is None): parser.error('Neighbour-file not provided')
    if (options.outFileName is None): parser.error('Output-file name not provided')

    processPlumes(options)
#end function

if __name__=="__main__":
    main()
#end if

