/*
 * =====================================================================================
 *
 * Copyright (C) 2014 Rakib Hassan (rakib.hassan@sydney.edu.au)
 *
 * This program is free software; you can redistribute it and/or modify it under 
 * the terms of the GNU General Public License as published by the Free Software 
 * Foundation; either version 2 of the License, or (at your option) any later 
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
 */
/*
 * =====================================================================================
 *
 *       Filename:  TimeSlices.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  20/07/14 14:02:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Rakib Hassan (rakib.hassan@sydney.edu.au)         
 *
 * =====================================================================================
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <TimeSlices.hh>
#include <TextFileReader.hh>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif
    #include <f2c.h>
    int intrc0_(integer *n, doublereal *plat, doublereal *plon, doublereal *x, 
        doublereal *y, doublereal *z__, doublereal *w, integer *list, integer *lptr, integer *
        lend, integer *ist, doublereal *pw, integer *ier);

    int intrc1_(integer *n, doublereal *plat, doublereal *plon, doublereal *x, 
        doublereal *y, doublereal *z__, doublereal *f, integer *list, integer *lptr, integer *
        lend, integer *iflgs, doublereal *sigma, integer *iflgg, doublereal *grad, 
        integer *ist, doublereal *fp, integer *ier);

    int trmesh_(integer *n, doublereal *x, doublereal *y, doublereal *z__, integer 
        *list, integer *lptr, integer *lend, integer *lnew, integer *near__, 
        integer *next, doublereal *dist, integer *ier);

    int trans_(integer *n, doublereal *rlat, doublereal *rlon, doublereal *x, doublereal 
        *y, doublereal *z__);

    int gradg_(integer *n, doublereal *x, doublereal *y, doublereal *z__, doublereal *f, 
        integer *list, integer *lptr, integer *lend, integer *iflgs, doublereal *
        sigma, integer *nit, doublereal *dgmax, doublereal *grad, integer *ier);

    int getsig_(integer *n, doublereal *x, doublereal *y, doublereal *z__, doublereal *
        h__, integer *list, integer *lptr, integer *lend, doublereal *grad, doublereal *
        tol, doublereal *sigma, doublereal *dsmax, integer *ier);
    
    void kcluster (int nclusters, int ngenes, int ndata, double** data,
    int** mask, double weight[], int transpose, int npass, char method, char dist,
    int clusterid[], double* error, int* ifound);
#ifdef __cplusplus
}
#endif

double LARGE_FLOAT = 1e30;

namespace src { namespace mesh {
using namespace std;
using src::io::TextFileReader;

    void TimeSlices::vxyz2vrtp(vector<float> &u1, vector<float> &u2, vector<float> &u3, 
                               vector<float> &x1, vector<float> &x2, vector<float> &x3,
                               vector<float> &vr, vector<float> &vtheta, vector<float> &vphi)
    {
        int dataNRows = u1.size();

        for (int i=0; i<dataNRows; i++)
        {
            double x = double(x1[i]);
            double y = double(x2[i]);
            double z = double(x3[i]);

            double r = sqrt(x*x + y*y + z*z);
            double eps = 2.220446049250313e-16;
            double xy = max(eps * r, sqrt(x*x + y*y));

            double T11 = x / r;
            double T21 = x * z / (xy * r);
            double T31 = -y / xy;
            double T12 = y / r;
            double T22 = y * z / (xy * r);
            double T32 = x / xy;
            double T13 = z / r;
            double T23 = -xy / r;
            double T33 = 0.;

            vr[i]      = T11*float(u1[i]) + T12*float(u2[i]) + T13*float(u3[i]);
            vtheta[i]  = T21*float(u1[i]) + T22*float(u2[i]) + T23*float(u3[i]);
            vphi[i]    = T31*float(u1[i]) + T32*float(u2[i]) + T33*float(u3[i]);
        }
    }
    
    void TimeSlices::vrtp2vxyz(vector<float> &vr, vector<float> &vtheta, vector<float> &vphi, 
                               vector<float> &r, vector<float> &theta, vector<float> &phi, 
                               vector<float> &vx, vector<float> &vy, vector<float> &vz)
    {
        int dataNRows = r.size();

        for (int i=0; i<dataNRows; i++)
        {
            float sinTheta = sin(theta[i]);
            float cosTheta = cos(theta[i]);
            
            float sinPhi = sin(phi[i]);
            float cosPhi = cos(phi[i]);

            float vr_i      = vr[i];
            float vtheta_i  = vtheta[i];
            float vphi_i    = vphi[i];

            vx[i] = sinTheta*cosPhi*vr_i + cosTheta*cosPhi*vtheta_i - sinPhi*vphi_i;
            vy[i] = sinTheta*sinPhi*vr_i + cosTheta*sinPhi*vtheta_i + cosPhi*vphi_i;
            vz[i] = cosTheta*vr_i - sinTheta*vtheta_i;
        }
    }

    TimeSlices::TimeSlices(string fBaseName, vector<float> times, vector<int> &timeSteps, 
                           float velocityScale, Geometry *g, bool processOpt)
    :m_fileBaseName(fBaseName),
    m_times(times),
    m_timeSteps(timeSteps),
    m_velocityScale(velocityScale),
    m_geometry(g),
    m_processOpt(processOpt)
    {
        /* Sort time and timeStep vector - although they should be sorted to begin with */
        sort(m_timeSteps.begin(), m_timeSteps.end()); 
        sort(m_times.begin(), m_times.end());         
        
        m_minTime = *min_element(m_times.begin(), m_times.end());
        m_maxTime = *max_element(m_times.begin(), m_times.end());
        
        /*-----------------------------------------------------------------------------
         * Resize storage vectors
         *-----------------------------------------------------------------------------*/
        m_timeStepDataList.resize(m_times.size());
        int dataNRows = m_geometry->m_r.size();
        for(unsigned int i=0; i<m_times.size(); i++)
        {
            m_timeStepDataList[i].m_vr.resize(dataNRows);
            m_timeStepDataList[i].m_vtheta.resize(dataNRows);
            m_timeStepDataList[i].m_vphi.resize(dataNRows);
            m_timeStepDataList[i].m_vx.resize(dataNRows);
            m_timeStepDataList[i].m_vy.resize(dataNRows);
            m_timeStepDataList[i].m_vz.resize(dataNRows);
            m_timeStepDataList[i].m_temperature.resize(dataNRows);
            m_timeStepDataList[i].m_viscosity.resize(dataNRows);
        }
        
        /* Need to check if optional files exist */
        m_hasOpt = 0;
        m_optNCol = 0;

        loadData();
        if(m_processOpt) loadOpt();
    }
    
    TimeSlices::~TimeSlices()
    {
        if(m_hasOpt)
        {
            for(unsigned int i=0; i<m_times.size(); i++)
            {
                for(unsigned int j=0; j<m_timeStepDataList[i].m_opt.size(); j++)
                {
                    delete m_timeStepDataList[i].m_opt[j];
                }
            }
        }
    }
    

    /*
     *--------------------------------------------------------------------------------------
     *       Class:  TimeSlices
     *      Method:  TimeSlices :: getClosestTimeIndex
     * Description:  Get index of time-slice closest to a given time instance
     *--------------------------------------------------------------------------------------
     */
    float TimeSlices::getClosestTimeIndex(float time, int *idx)
    {
        float minDiff = 1e32;
        *idx = -1;
        for(unsigned int i=0; i<m_times.size(); i++)
        {
            float diff = fabs(m_times[i] - time);

            if(diff < minDiff)
            {
                minDiff = diff;
                *idx = i;
            }
        }

        return minDiff;
    }

    /*
     *--------------------------------------------------------------------------------------
     *       Class:  TimeSlices
     *      Method:  TimeSlices :: getTimeExtentIndices
     * Description:  Get time-slice indices that 'time' falls in between
     *--------------------------------------------------------------------------------------
     */
    bool TimeSlices::getTimeExtentIndices(float time, int *lowIdx, int *highIdx)
    {
        *lowIdx = -1;
        *highIdx = -1;

        for(unsigned int i=0; i<m_times.size(); i++)
        {
            if(m_times[i]<=time)
            {
                *lowIdx = i;
            }
            else
            {
                break;
            }
        }
        
        *highIdx = *lowIdx;
        for(unsigned int i=0; i<m_times.size(); i++)
        {
            if(m_times[i]>time)
            {
                *highIdx = i;
                break;
            }
        }
        
        if((*lowIdx==-1) || (*highIdx==-1)) 
            return false; /* Error */
        else 
            return true;
    }    
    
    /*
     *--------------------------------------------------------------------------------------
     *       Class:  TimeSlices
     *      Method:  TimeSlices :: getPlumeMetric
     * Description:  Scans model-space at 'depth (m)' using a threshold dT and returns a 
     *               list of indices.
     *--------------------------------------------------------------------------------------
     */
    void TimeSlices::getPlumeMetric(float time, float depth, float cutOffPercentile,
                            int *closestTimeIdx, double *closestR, vector<double> *gradVr, 
                            vector<double> *vr, vector<double> *temperature, vector<double> *eta,
                            ClusterResult *cr, int flavourId, vector<double> *composition)
    {
        getGradVr(time, depth, gradVr, closestTimeIdx, closestR, vr, temperature, eta);
        
        /*
        vector<float> Vrr(m_geometry->m_shellNumNodes);
        vector<float> Vrt(m_geometry->m_shellNumNodes);
        vector<float> Vrp(m_geometry->m_shellNumNodes);
        {
            int shellNumNodes = m_geometry->m_shellNumNodes;
            
            vector<float> Vrx(shellNumNodes);
            vector<float> Vry(shellNumNodes);
            vector<float> Vrz(shellNumNodes);
            
            for(int i=0; i<shellNumNodes; i++)
            {
                Vrx[i] = (*gradVr)[i*3+0];
                Vry[i] = (*gradVr)[i*3+1];
                Vrz[i] = (*gradVr)[i*3+2];
            }
            
            vector<float> x(m_geometry->m_uxs.begin(), m_geometry->m_uxs.end());
            vector<float> y(m_geometry->m_uys.begin(), m_geometry->m_uys.end());
            vector<float> z(m_geometry->m_uzs.begin(), m_geometry->m_uzs.end());
            
            vxyz2vrtp(Vrx, Vry, Vrz, x, y, z, Vrr, Vrt, Vrp);
        }*/
        
        /*-----------------------------------------------------------------------------
         * Perform k-means clustering on vr_x
         *-----------------------------------------------------------------------------*/
        {   
            int shellNumNodes = m_geometry->m_shellNumNodes;
            double weight[1] = {1};
            double error;
            int ifound;
            double **data = new double*[shellNumNodes];
            int **mask = new int*[shellNumNodes];
            data[0] = new double[shellNumNodes];
            mask[0] = new int[shellNumNodes];

            cr->m_cid = vector<int>(shellNumNodes);
            cr->m_mask = vector<int>(shellNumNodes);
            cr->m_clusterInput = vector<float>(shellNumNodes);
            
            /*-----------------------------------------------------------------------------
             * Compute magnitude of gradVr
             *-----------------------------------------------------------------------------*/
            vector<double> magGradVr(m_geometry->m_shellNumNodes);
            for(int i=0; i<shellNumNodes; i++)
            {
                magGradVr[i] = sqrt( pow((*gradVr)[i*3+0], 2.) + 
                                     pow((*gradVr)[i*3+1], 2.) + 
                                     pow((*gradVr)[i*3+2], 2.) );
            }
            vector<double> magGradVrSorted(magGradVr);
            sort(magGradVrSorted.begin(), magGradVrSorted.end(), greater<double>());
            
            int dataCount = 0;
            for(int i=0; i<shellNumNodes; i++)
            {
                data[i] = data[0]+i;
                mask[i] = mask[0]+i;
                
                if( (magGradVr[i] > magGradVrSorted[shellNumNodes*cutOffPercentile/100.]) &&
                    ((*vr)[i]>0) )
                {
                    dataCount++;
                    data[i][0] = magGradVr[i];
                    mask[i][0] = 1;
                }
                else
                {
                    data[i][0] = 0;
                    mask[i][0] = 0;
                }

                /* save cluster-data and mask */
                cr->m_clusterInput[i] = data[i][0];
                cr->m_mask[i] = mask[i][0];
            }
            
            if(dataCount==0)
            {
                cout << "\tNo data found. Skipping cluster analysis.." << endl;
            }
            else
            {
                cout << "\tPerforming cluster analysis at " << depth/1e3 << " km depth.." << endl;
                
                kcluster( NCLUSTER,  shellNumNodes, 1, data, mask, weight, 0, 10, 'a', 'b', cr->m_cid.data(), &error, &ifound);
            }
            
            /* Debug output */
            /*{
                char fname[1024] = {};

                sprintf(fname, "/tmp/cluster%d.txt", m_callCount++);
                FILE *f = fopen(fname, "w");
                for(int i=0; i<m_geometry->m_shellNumNodes; i++)
                {
                    fprintf(f, "%f %f %d %d\n", m_geometry->m_ulons[i], m_geometry->m_ulats[i], mask[i][0], m_cid[i]);
                }
                fclose(f);
            }*/

            /*cout << "Neighbours" << endl;
            for(unsigned int i=0; i<m_geometry->m_ulons.size(); i++)
            {
                cout << i;

                for( set<int>::iterator nit = m_geometry->m_unodeNeighbours[i]->begin(); 
                     nit != m_geometry->m_unodeNeighbours[i]->end(); nit++)
                {
                    cout << " " << *nit;
                }

                cout << endl;
            }*/
            
            if(dataCount)
            {
                /*-----------------------------------------------------------------------------
                 * The clustering algorithm labels hotter and colder regions arbitrarily, which 
                 * necessitates identifying which cluster-id we should propagate for locating 
                 * blob-centres
                 *-----------------------------------------------------------------------------*/
                
                int cid0 = 0; /* define cluster-ids */
                int cid1 = 1;
                map<int, int> cidHitCount;
                
                cidHitCount.insert(pair<int, int>(cid0, 0));
                cidHitCount.insert(pair<int, int>(cid1, 0));
                for(int i=0; i<shellNumNodes; i++)
                {
                    if(mask[i][0])
                    {
                        if(cr->m_cid[i]==cid0) cidHitCount[cid0]++;
                        if(cr->m_cid[i]==cid1) cidHitCount[cid1]++;
                    }
                    else
                    {
                        cr->m_cid[i] = -1;
                    }
                }

                printf("\tCluster population-count %d(%d), %d(%d)\n", cid0, cidHitCount[cid0], cid1, cidHitCount[cid1]);
                
                /*-----------------------------------------------------------------------------
                 * Choose cluster-id based on population-count. The hot-rounded blobs we seek, 
                 * which represent plume birth-places, should almost always belong to the cluster 
                 * with a smaller population.
                 *-----------------------------------------------------------------------------*/
                int propClusterId = cid0;
                if(cidHitCount[cid0]>cidHitCount[cid1]) propClusterId = cid1;

                /*-----------------------------------------------------------------------------
                 * Switch cluster ids as per population-count - such that 0 is always the
                 * cluster-id for the smaller population count.
                 *-----------------------------------------------------------------------------*/
                if(propClusterId!=0)
                {
                    for(int i=0; i<shellNumNodes; i++)
                    {
                        if(mask[i][0]) cr->m_cid[i] = !(cr->m_cid[i]);
                    }
                }
            }

            /* Release allocated memory */
            delete data[0];
            delete mask[0];
            delete data;
            delete mask;
        }

        /*-----------------------------------------------------------------------------
         * Tracer composition related 
         *-----------------------------------------------------------------------------*/
        int rid = -1;
        double dummy;
        vector<int>    gids(m_geometry->m_shellNumNodes);
        if(m_hasOpt)
        {
            float r = (m_geometry->m_rScale - depth)/m_geometry->m_rScale;
            rid     = m_geometry->getClosestShellRadiusIndex(r, &dummy);

            int count = 0;
            for(set<GeoLocation, GeoCompareStruct>::iterator sit = m_geometry->m_shells[rid].begin();
                sit != m_geometry->m_shells[rid].end(); sit++, count++) 
            {
                gids[count] = sit->globalId;
            }        

            composition->resize(m_geometry->m_shellNumNodes);
            for(int i=0; i<m_geometry->m_shellNumNodes; i++)
            {
                int oid            = m_geometry->m_order[i];
                
                (*composition)[i]  = (*(m_timeStepDataList[*closestTimeIdx].m_opt[flavourId-1]))[gids[oid]];
            }
        }
    }

    /*
     *--------------------------------------------------------------------------------------
     *       Class:  TimeSlices
     *      Method:  TimeSlices :: computeGradient
     * Description:  Computes spatial gradient of a given field.
     *--------------------------------------------------------------------------------------
     */
    void TimeSlices::computeGradient(vector<double> *field, vector<double> *result)
    {
        result->clear();
        result->resize(m_geometry->m_shellNumNodes * 3);
        /*-----------------------------------------------------------------------------
         * Compute gradient
         *-----------------------------------------------------------------------------*/
        {
            vector<double> sigma(m_geometry->m_shellNumNodes*6);
            int iflgs = 0; /* uniform tension */
            int ier = 0; /* error reporter */
            int maxit = 10;
            int nit = maxit;
            double dgmax = 1e-7;
            double dgmx = dgmax;
            
            vector<double> sxs(m_geometry->m_shellNumNodes);
            vector<double> sys(m_geometry->m_shellNumNodes);
            vector<double> szs(m_geometry->m_shellNumNodes);
            vector<double> slats(m_geometry->m_shellNumNodes);
            vector<double> slons(m_geometry->m_shellNumNodes);
            
            for(int i=0; i<m_geometry->m_shellNumNodes; i++) 
            {
                slons[i] = m_geometry->m_slons[i] + drand48()/1e4; /* Add white noise */
                slats[i] = m_geometry->m_slats[i] + drand48()/1e4; /* Add white noise */
            }
    
            trans_( &(m_geometry->m_shellNumNodes), 
                    slats.data(), 
                    slons.data(), 
                    sxs.data(), 
                    sys.data(), 
                    szs.data() );
    
            gradg_( &(m_geometry->m_shellNumNodes), 
                    sxs.data(), 
                    sys.data(), 
                    szs.data(), 
                    field->data(), 
                    m_geometry->m_list.data(), 
                    m_geometry->m_lptr.data(), 
                    m_geometry->m_lend.data(),                 
                    &iflgs, 
                    sigma.data(), 
                    &maxit, 
                    &dgmx, 
                    result->data(), 
                    &ier );      
            
            fprintf (stderr, "%s: GRADG :  tolerance = %g max change = %g  maxit = %d no. iterations = %d ier = %d\n",
                        "SphericalSplineInterpolator", dgmax, dgmx, maxit, nit, ier);         
            if(ier<0) 
            {
                cout << "Error in gradg_" << endl;
                exit(EXIT_FAILURE);
            }
        }
    }

    /*
     *--------------------------------------------------------------------------------------
     *       Class:  TimeSlices
     *      Method:  TimeSlices :: getGradVr
     * Description:  
     *--------------------------------------------------------------------------------------
     */
    void TimeSlices::getGradVr(float time, float depth, vector<double> *result, int *closestTimeIdx,
                              double *closestR, vector<double> *vr, vector<double> *temperature, 
                              vector<double> *eta)
    {
        float r = (m_geometry->m_rScale - depth)/m_geometry->m_rScale;
        int rid = m_geometry->getClosestShellRadiusIndex(r, closestR);
        int lowIdx, highIdx;
        bool retVal;

        /*-----------------------------------------------------------------------------
         * Reset output arrays 
         *-----------------------------------------------------------------------------*/
        vr->clear(); vr->resize(m_geometry->m_shellNumNodes);
        temperature->clear(); temperature->resize(m_geometry->m_shellNumNodes);
        eta->clear(); eta->resize(m_geometry->m_shellNumNodes);

        /*-----------------------------------------------------------------------------
         * Find closest timeSlice index
         *-----------------------------------------------------------------------------*/
        int timeSliceIdx = -1;
        
        retVal = getTimeExtentIndices(time, &lowIdx, &highIdx);
        assert(retVal);

        if(fabs(time-m_times[lowIdx]) < fabs(time-m_times[highIdx])) timeSliceIdx = lowIdx;
        else timeSliceIdx = highIdx;
        
        
        /*-----------------------------------------------------------------------------
         * Gather temperature field 
         *-----------------------------------------------------------------------------*/
        vector<double> field(m_geometry->m_shellNumNodes);
        vector<int>    gids(m_geometry->m_shellNumNodes);

        int count = 0;
        for(set<GeoLocation, GeoCompareStruct>::iterator sit = m_geometry->m_shells[rid].begin();
            sit != m_geometry->m_shells[rid].end(); sit++, count++) 
        {
            gids[count] = sit->globalId;
        }
        
        /*-----------------------------------------------------------------------------
         * Get ordered field
         *-----------------------------------------------------------------------------*/
        for(int i=0; i<m_geometry->m_shellNumNodes; i++)
        {
            int oid            = m_geometry->m_order[i];
            
            field[i]           = m_timeStepDataList[timeSliceIdx].m_vr[gids[oid]];
            (*vr)[i]           = m_timeStepDataList[timeSliceIdx].m_vr[gids[oid]];
            (*temperature)[i]  = m_timeStepDataList[timeSliceIdx].m_temperature[gids[oid]];
            (*eta)[i]          = m_timeStepDataList[timeSliceIdx].m_viscosity[gids[oid]];
        }
        
        computeGradient(&field, result);
        
        *closestR       *= m_geometry->m_rScale;
        *closestTimeIdx = timeSliceIdx;
    }

    /*
     *--------------------------------------------------------------------------------------
     *       Class:  TimeSlices
     *      Method:  TimeSlices :: getShellVelocities
     * Description:  Returns velocities at given locations from a citcom shell. Note that 
     *               colats->[-pi/2, pi/2] and lons->[0, 2pi]
     *--------------------------------------------------------------------------------------
     */
    void TimeSlices::getShellVelocities(int timeIndex, int shellIndex, 
                                                   float lateralRadius,
                                                   vector<float> &colats, 
                                                   vector<float> &lons, 
                                                   vector<float> *vr,
                                                   vector<float> *vt,
                                                   vector<float> *vp)
    {
        /* Sanity check */
        assert(colats.size() == lons.size());
        
        /* Resize result-arrays */
        vr->clear(); vr->resize(colats.size());
        vt->clear(); vt->resize(colats.size());
        vp->clear(); vp->resize(colats.size());

        vector<int> gids(m_geometry->m_shellNumNodes);
        /*-----------------------------------------------------------------------------
         * Get global IDs for the requested shell
         *-----------------------------------------------------------------------------*/
        int count = 0;
        for(set<GeoLocation, GeoCompareStruct>::iterator sit = m_geometry->m_shells[shellIndex].begin();
            sit != m_geometry->m_shells[shellIndex].end(); sit++, count++) 
        {
            gids[count] = sit->globalId;
        }


        /*-----------------------------------------------------------------------------
         * Nearest node interpolation for the time being 
         * TODO: implement IDW interpolation
         *-----------------------------------------------------------------------------*/
        float tp[2] = {0};
        for(unsigned int i=0; i<colats.size(); i++)
        {
            vector<float> distances;
            vector<int> indices;
            
            tp[0] = colats[i]; // convert to colat
            tp[1] = lons[i];
            
            m_geometry->queryShellPoint(tp, lateralRadius, &distances, &indices);

            int idx = gids[indices[0]]; // finally, get gidx after undoing shuffle
            (*vr)[i] = m_timeStepDataList[timeIndex].m_vr[idx];
            (*vt)[i] = m_timeStepDataList[timeIndex].m_vtheta[idx];
            (*vp)[i] = m_timeStepDataList[timeIndex].m_vphi[idx];
        }
    }

    /*
     *--------------------------------------------------------------------------------------
     *       Class:  TimeSlices
     *      Method:  TimeSlices :: loadOpt
     * Description:  Setup storage for opt-data
     *--------------------------------------------------------------------------------------
     */
    void TimeSlices::loadOpt()
    {
        struct Splitter
        {
            static vector<string> Split(string str, char delimiter) 
            {
                vector<string> internal;
                stringstream ss(str); // Turn the string into a stream.
                string tok;
              
                while(getline(ss, tok, delimiter)) 
                {
                    internal.push_back(tok);
                }
              
                return internal;
            }
        };
        
        /*-----------------------------------------------------------------------------
         * Check if they exist and get number of columns
         *-----------------------------------------------------------------------------*/
        {
            char fName[1024] = {};
            sprintf(fName, "%s.opt%02d.%d", m_fileBaseName.c_str(), 0, m_timeSteps[0]);
            ifstream f(fName);

            if(f.good())
            {
                cout << "\nFound opt files..\n";
                m_hasOpt = 1;
                
                string line;
                
                getline(f,line); // get header
                getline(f,line); // get first line

                vector<string> tokens = Splitter::Split(line, ' ');
                m_optNCol = int(tokens.size());
                
                int dataNRows = m_geometry->m_r.size();
                
                for(unsigned int i=0; i<m_times.size(); i++)
                {
                    for(int j=0; j<m_optNCol; j++)
                    {
                        m_timeStepDataList[i].m_opt.push_back(new vector<float>);
                        m_timeStepDataList[i].m_opt[j]->resize(dataNRows);
                    }
                }
                f.close();
            }
        }
        
        if(m_hasOpt==0) 
        {
            cout << "\nNo opt files detected..\n";
            return;
        }
        /*-----------------------------------------------------------------------------
         * Load opt data
         *-----------------------------------------------------------------------------*/
        int length = m_times.size();
        for (int itime=0; itime < length; itime++)
        {
            /*-----------------------------------------------------------------------------
             * Status
             *-----------------------------------------------------------------------------*/
            {
                cout << "\tLoading opt-data for time-step " << m_timeSteps[itime]
                     << " (time: " << m_times[itime] << " Myr) .." << endl;
            }

            int counter = 0;
            for(int icap=0; icap<NCAP; icap++)
            {
                char fName[1024] = {};

                sprintf(fName, "%s.opt%02d.%d", m_fileBaseName.c_str(), icap, m_timeSteps[itime]);

                //cout << fName << endl;
                TextFileReader tr(fName, 1);
                
                vector<int> c;
                for(int i=0; i<m_optNCol; i++) c.push_back(1); // get all columns
                vector< vector<float>* > res;
                
                /* Read data */
                tr.read(c, &res);
                
                int dataNRows = res[0]->size();
                
                /* Note: (vx,vy,vz) <- (vr, vt, vp) */
                for(int i=0; i<dataNRows; i++)
                {
                    for(int j=0; j<m_optNCol; j++)
                    {
                        (*(m_timeStepDataList[itime].m_opt[j]))[counter] = (*res[j])[i]; //vr

                    }
                    counter++;
                }
                
                for (unsigned int i=0; i<c.size(); i++) delete res[i];
            }        
        }
    }

    /*
     *--------------------------------------------------------------------------------------
     *       Class:  TimeSlices
     *      Method:  TimeSlices :: loadData
     * Description:  Loads time-slice data
     *--------------------------------------------------------------------------------------
     */
    void TimeSlices::loadData()
    {
        int length = m_times.size();
        
        for (int itime=0; itime < length; itime++)
        {
            /*-----------------------------------------------------------------------------
             * Status
             *-----------------------------------------------------------------------------*/
            {
                cout << "\tLoading data for time-step " << m_timeSteps[itime]
                     << " (time: " << m_times[itime] << " Myr) .." << endl;
            }

            int counter = 0;
            map<float, float> shellNodeCount;
            for(int icap=0; icap<NCAP; icap++)
            {
                char fName[1024] = {};

                sprintf(fName, "%s.cap%02d.%d", m_fileBaseName.c_str(), icap, m_timeSteps[itime]);

                //cout << fName << endl;
                TextFileReader tr(fName, 1);
                
                int colsArray[] = {0,0,0,1,1,1,1,1};
                vector<int> c(8);
                vector< vector<float>* > res;
                
                copy(colsArray, colsArray+8, c.begin());

                /* Read data */
                tr.read(c, &res);
                
                int dataNRows = res[3]->size();
                
                /* Note: (vx,vy,vz) <- (vr, vt, vp) */
                for(int i=0; i<dataNRows; i++)
                {
                    m_timeStepDataList[itime].m_vr[counter]             = (*res[5])[i]; //vr
                    m_timeStepDataList[itime].m_vtheta[counter]         = (*res[3])[i]; //vt
                    m_timeStepDataList[itime].m_vphi[counter]           = (*res[4])[i]; //vp
                    
                    /* Note: vx,vy,vz are transformed further down */
                    m_timeStepDataList[itime].m_vx[counter]             = (*res[5])[i]; //vr
                    m_timeStepDataList[itime].m_vy[counter]             = (*res[3])[i]; //vt
                    m_timeStepDataList[itime].m_vz[counter]             = (*res[4])[i]; //vp
                    m_timeStepDataList[itime].m_temperature[counter]    = (*res[6])[i];
                    m_timeStepDataList[itime].m_viscosity[counter]      = (*res[7])[i];

                    /* Average-related */
                    float r = m_geometry->m_r[i];
                    m_timeStepDataList[itime].m_vrShellAverage[r] += m_timeStepDataList[itime].m_vr[counter];
                    m_timeStepDataList[itime].m_temperatureShellAverage[r] += m_timeStepDataList[itime].m_temperature[counter];
                    m_timeStepDataList[itime].m_viscosityShellAverage[r] += m_timeStepDataList[itime].m_viscosity[counter];
                    shellNodeCount[r] += 1;

                    counter++;
                }
                
                for (unsigned int i=0; i<c.size(); i++) delete res[i];
            }
            
            /* Compute averages for each shell */
            map<float, float>::iterator vrAveIt = m_timeStepDataList[itime].m_vrShellAverage.begin();
            map<float, float>::iterator tAveIt = m_timeStepDataList[itime].m_temperatureShellAverage.begin();
            map<float, float>::iterator vAveIt = m_timeStepDataList[itime].m_viscosityShellAverage.begin();
            map<float, float>::iterator countIt = shellNodeCount.begin();

            for( ; (vrAveIt != m_timeStepDataList[itime].m_vrShellAverage.end()) &&
                   (tAveIt != m_timeStepDataList[itime].m_temperatureShellAverage.end()) &&
                   (vAveIt != m_timeStepDataList[itime].m_viscosityShellAverage.end()) &&
                   (countIt != shellNodeCount.end());
                   vrAveIt++, tAveIt++, vAveIt++, countIt++ )
            {
                m_timeStepDataList[itime].m_vrShellAverage[tAveIt->first] /= countIt->second;
                m_timeStepDataList[itime].m_temperatureShellAverage[tAveIt->first] /= countIt->second;
                m_timeStepDataList[itime].m_viscosityShellAverage[vAveIt->first] /= countIt->second;
            }

            /*-----------------------------------------------------------------------------
             * Compute m_temperatureShellAverageHot for each shell, which simply excludes
             * material colder than m_temperatureShellAverage at each depth-level. This is 
             * useful for tracking active upwellings better.
             *-----------------------------------------------------------------------------*/
            map<float, float> shellNodeCountHot;
            for(int inode=0; inode<counter; inode++)
            {
                float r = m_geometry->m_r[inode];
                if(m_timeStepDataList[itime].m_temperature[inode] >= m_timeStepDataList[itime].m_temperatureShellAverage[r])
                {
                    m_timeStepDataList[itime].m_temperatureShellAverageHot[r] += m_timeStepDataList[itime].m_temperature[inode];
                    shellNodeCountHot[r] += 1;
                }
            }
            
            tAveIt = m_timeStepDataList[itime].m_temperatureShellAverageHot.begin();
            countIt = shellNodeCountHot.begin();
            for( ; (tAveIt != m_timeStepDataList[itime].m_temperatureShellAverageHot.end()) &&
                   (countIt != shellNodeCountHot.end());
                   tAveIt++, countIt++ )
            {
                m_timeStepDataList[itime].m_temperatureShellAverageHot[tAveIt->first] /= countIt->second;
                
                /*cout << "average: " << m_timeStepDataList[itime].m_temperatureShellAverage[tAveIt->first] << 
                        " hot-average: " << m_timeStepDataList[itime].m_temperatureShellAverageHot[tAveIt->first] << endl;*/
            }

            /*cout << "tAve-map size " << m_timeStepDataList[itime].m_temperatureShellAverage.size() << endl;
            cout << "vAve-map size " << m_timeStepDataList[itime].m_viscosityShellAverage.size() << endl;
            cout << "count-map size " << shellNodeCount.size() << endl;*/
        }

        /* Note input and output vectors are the same here */
        for(int itime=0; itime<length; itime++)
        {
            vrtp2vxyz( m_timeStepDataList[itime].m_vx, 
                       m_timeStepDataList[itime].m_vy, 
                       m_timeStepDataList[itime].m_vz,
                       
                       m_geometry->m_r, 
                       m_geometry->m_theta,
                       m_geometry->m_phi,
                       
                       m_timeStepDataList[itime].m_vx, 
                       m_timeStepDataList[itime].m_vy, 
                       m_timeStepDataList[itime].m_vz );
        }
        
        /* applying velocity-scale */
        {
            for(int itime=0; itime<length; itime++)
            {
                int dataNRows = m_timeStepDataList[itime].m_vx.size();
                for (int i=0; i<dataNRows; i++)
                {
                    m_timeStepDataList[itime].m_vx[i] *= m_velocityScale;
                    m_timeStepDataList[itime].m_vy[i] *= m_velocityScale;
                    m_timeStepDataList[itime].m_vz[i] *= m_velocityScale;
                }
            }
        }
    }
}}
