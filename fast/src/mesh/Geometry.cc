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
 *       Filename:  Geometry.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  19/07/14 14:02:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Rakib Hassan (rakib.hassan@sydney.edu.au)         
 *
 * =====================================================================================
 */

#include <math.h>
#include <iostream>
#include <Geometry.hh>
#include <TextFileReader.hh>

#ifdef __cplusplus
extern "C" {
#endif

#include <f2c.h>
int trmesh_(integer *n, doublereal *x, doublereal *y, doublereal *z__, integer 
    *list, integer *lptr, integer *lend, integer *lnew, integer *near__, 
    integer *next, doublereal *dist, integer *ier);

int trlist_(integer *n, integer *list, integer *lptr, 
    integer *lend, integer *nrow, integer *nt, integer *ltri, integer *
    ier);

int trans_(integer *n, doublereal *rlat, doublereal *rlon, doublereal *x, doublereal 
    *y, doublereal *z__);

#ifdef __cplusplus
}
#endif

namespace src { namespace mesh {

using namespace std;
using src::io::TextFileReader;
using src::mesh::KdTree;
 
    void Geometry::rtp2xyz()
    {
        unsigned int length = m_r.size();

        for(unsigned int i=0; i<length; i++)
        {
            double r        = m_r[i];
            double theta    = m_theta[i];
            double phi      = m_phi[i];

            double rst = r * sin(theta);
            
            m_x[i] = rst * cos(phi);
            m_y[i] = rst * sin(phi); 
            m_z[i] = r * cos(theta);  
        }
    }

    void Geometry::xyz2rtp(float *xyz, float *rtp)
    {
        float x = xyz[0];
        float y = xyz[1];
        float z = xyz[2];

        float tmp1 = x*x + y*y;
        float tmp2 = tmp1 + z*z;
        
        rtp[0] = sqrt(tmp2);     
        rtp[1] = atan2(sqrt(tmp1),z); 
        rtp[2] = atan2(y,x);  
    }

    void Geometry::rtp2xyz(float *rtp, float *xyz)
    {
        double r        = rtp[0];
        double theta    = rtp[1];
        double phi      = rtp[2];

        double rst = r * sin(theta);
        
        xyz[0] = rst * cos(phi);      // x 
        xyz[1] = rst * sin(phi);      // y 
        xyz[2] = r * cos(theta);      // z
    }
    
    float Geometry::getDistance(int nid1, int nid2)
    {
        float dx = m_x[nid1] - m_x[nid2];
        float dy = m_y[nid1] - m_y[nid2];
        float dz = m_z[nid1] - m_z[nid2];

        return sqrt(dx*dx + dy*dy + dz*dz);
    }

    /*
     *--------------------------------------------------------------------------------------
     *       Class:  Geometry
     *      Method:  Geometry :: Geometry
     * Description:  Constructor
     *--------------------------------------------------------------------------------------
     */
    Geometry::Geometry(string fBaseName, float rScale, bool initializeTriangulation)
    :m_initializeTriangulation(initializeTriangulation),
    m_fileBaseName(fBaseName),
    m_rScale(rScale)
    {
        /*-----------------------------------------------------------------------------
         * Instantiate KdTrees
         *-----------------------------------------------------------------------------*/
        m_tree      = new KdTree();
        m_shellTree = new KdTree();

        for(int icap=0; icap<NCAP; icap++)
        {
            char fName[1024] = {};

            sprintf(fName, "%s.cap%02d.0", m_fileBaseName.c_str(), icap);

            //cout << fName << endl;
            TextFileReader tr(fName, 1);
            
            int colsArray[] = {1,1,1,0,0,0,0,0};
            vector<int> c(8);
            vector< vector<float>* > res;
            
            copy(colsArray, colsArray+8, c.begin());

            /* Read coordinates */
            tr.read(c, &res);
            
            
            /* Save rtp triplets */
            int length = res[0]->size();            
            for (int i=0; i<length; i++)
            {
                float phi = (*res[1])[i];

                if(phi > PI) phi -= 2*PI;
                
                m_r.push_back      ((*res[2])[i] * m_rScale); // applying r-scale
                m_theta.push_back  ((*res[0])[i]);
                m_phi.push_back    (phi);
                
                m_rNormalized.push_back((*res[2])[i]); // normalized r

                m_ur.insert((*res[2])[i]); // unique normalized r
            }
            
            for (unsigned int i=0; i<c.size(); i++) delete res[i];
        }   
        
        /*-----------------------------------------------------------------------------
         * Initialize indices for unique r-values
         *-----------------------------------------------------------------------------*/
        {
            int ucount = 0;
            for(set<float>::iterator sit = m_ur.begin(); sit != m_ur.end(); sit++, ucount++)
            {
                m_uri[*sit] = ucount;
                m_url.push_back(*sit);
            }
        }

        /* Compute and store xyz triplets */
        int length = m_r.size();
        m_x.resize(length);
        m_y.resize(length);
        m_z.resize(length);
        rtp2xyz();
        
        /* Populate tree */
        {
            float coord[3];
            int counter = 0;
            for(int i=0; i<length; i++) 
            {
                coord[0] = m_x[i];
                coord[1] = m_y[i];
                coord[2] = m_z[i];
                
                //cout << "adding: " << coord[0] << " " << coord[1] << " " << coord[2] << endl;

                m_tree->Add(coord, counter);

                /*-----------------------------------------------------------------------------
                 * Initialize shell-related structures
                 *-----------------------------------------------------------------------------*/
                pair<set<GeoLocation, GeoCompareStruct>::iterator, bool> ret;
                GeoLocation ugl;
                float phi = m_phi[i];
                if(phi < 0) phi += 2*PI; /* phi must be in [0, 360] */
                float theta = PI/2 - m_theta[i]; /* converting from co-lat */
                
                ugl.lon = phi;
                ugl.lat = theta;
                ugl.globalId = counter;
                
                ret = m_shells[m_uri[m_rNormalized[counter]]].insert(ugl);

                counter++;
            }        
        }
        m_shellNumNodes = m_shells[0].size();

        m_minR = *min_element(m_r.begin(), m_r.end());
        m_maxR = *max_element(m_r.begin(), m_r.end());

        /*-----------------------------------------------------------------------------
         * Initialize triangulation
         *-----------------------------------------------------------------------------*/
        if(m_initializeTriangulation) 
        {
            initTriangulation();

            /*-----------------------------------------------------------------------------
             * Initialize shell-tree
             *-----------------------------------------------------------------------------*/
            float coord[3];
            int counter = 0;
            for(set<GeoLocation, GeoCompareStruct>::iterator sit = m_shells[0].begin();
                sit != m_shells[0].end(); sit++, counter++)
            {
                int gid = sit->globalId;

                coord[0] = m_x[gid];
                coord[1] = m_y[gid];
                coord[2] = m_z[gid];
                
                m_shellTree->Add(coord, counter);
            }        
        }
    }
    
    int Geometry::getClosestShellRadiusIndex(float r, double *closestR)
    {
        int result = -1;
        float minDist = float(INT_MAX); 
        /*-----------------------------------------------------------------------------
         * r must in normalized form 
         *-----------------------------------------------------------------------------*/
        int ucount = 0;
        for(set<float>::iterator sit = m_ur.begin(); sit != m_ur.end(); sit++, ucount++)
        {
            if(fabs(*sit - r) < minDist)
            {
                minDist   = fabs(*sit - r);
                result    = ucount;
                *closestR = *sit;
            }
        }

        return result;
    }

    void Geometry::queryBallPoint(float *pos, float r, vector<float> *distance, vector<int> *id)
    {
        if(!m_tree) return;

        m_tree->QueryBallPoint(pos, r, distance, id);
    }
    
    /*
     *--------------------------------------------------------------------------------------
     *       Class:  Geometry
     *      Method:  Geometry :: queryShellPoint
     * Description:  Returns closest points to tp(theta, phi) within a distance 
     *               of lateralRadius.
     *--------------------------------------------------------------------------------------
     */
    void Geometry::queryShellPoint(float *tp, float lateralRadius, vector<float> *distance, vector<int> *id)
    {
        float rtp[3] = {0};
        float xyz[3] = {0};
        if(!m_shellTree) return;
        
        rtp[0] = m_url[0]*m_rScale; // The shell-tree is populated using the CMB shell
        rtp[1] = tp[0];
        rtp[2] = tp[1];
        
        rtp2xyz(rtp, xyz);

        m_shellTree->QueryBallPoint(xyz, lateralRadius, distance, id);
    }

    void Geometry::initTriangulation()
    {
        cout << "\tComputing triangulation.. ";
        
        for(int i=0; i<m_shellNumNodes; i++)
        {
            m_order.push_back(i);
        }

        m_slons.resize(m_shellNumNodes);
        m_slats.resize(m_shellNumNodes);
        m_sxs.resize(m_shellNumNodes);
        m_sys.resize(m_shellNumNodes);
        m_szs.resize(m_shellNumNodes);
        m_nodeNeighbours.resize(m_shellNumNodes);

        /*-----------------------------------------------------------------------------
         * Shuffle data to ensure the first 3 points are not colinear
         *-----------------------------------------------------------------------------*/
        random_shuffle(m_order.begin(), m_order.end(), myrandom);
        
        vector<float>lonTemp(m_shellNumNodes);
        vector<float>latTemp(m_shellNumNodes);
        int sitCounter = 0;
        for(set<GeoLocation, GeoCompareStruct>::iterator sit = m_shells[0].begin();
            sit != m_shells[0].end(); sit++, sitCounter++) 
        {
            lonTemp[sitCounter] = sit->lon;
            latTemp[sitCounter] = sit->lat;
        }
        
        for(int i=0; i<m_shellNumNodes; i++) 
        {
            m_slons[i] = lonTemp[m_order[i]];
            m_slats[i] = latTemp[m_order[i]];
        }
        
        /*-----------------------------------------------------------------------------
         *  Convert nodes to cartesian.
         *-----------------------------------------------------------------------------*/
        trans_( &m_shellNumNodes, 
                m_slats.data(), 
                m_slons.data(), 
                m_sxs.data(), 
                m_sys.data(), 
                m_szs.data() );
        
        /*-----------------------------------------------------------------------------
         * Create triangulation 
         *-----------------------------------------------------------------------------*/
        m_list.resize(m_shellNumNodes*6);
        m_lptr.resize(m_shellNumNodes*6);
        m_lend.resize(m_shellNumNodes);
        m_near.resize(m_shellNumNodes);
        m_next.resize(m_shellNumNodes);
        m_dist.resize(m_shellNumNodes);
         
        int lnew  = 0;
        int error = 0;
        trmesh_( &m_shellNumNodes, 
                 m_sxs.data(), 
                 m_sys.data(), 
                 m_szs.data(), 
                 m_list.data(), 
                 m_lptr.data(), 
                 m_lend.data(), 
                 &lnew, 
                 m_near.data(), 
                 m_next.data(), 
                 m_dist.data(), 
                 &error );
        
        if(error!=0)
        {
            cout << "Error encountered in triangulation routine: " << error << endl;
            exit(EXIT_FAILURE);
        }
        cout << "\tDone.. " << endl;

        /*-----------------------------------------------------------------------------
         * Get triangle-indices and neighbours
         *-----------------------------------------------------------------------------*/
        vector<int> triangleIndicesNeighbours(12*m_shellNumNodes);
        for(int i=0; i<m_shellNumNodes; i++)
        {
            m_nodeNeighbours[i] = new set<int>;
        }
        
        
        cout << "\tComputing node and triangle neighbours.. ";
        int nrow = 6; /* We just want triangle indices and triangle neighbours */
        trlist_( &m_shellNumNodes, 
                 m_list.data(), 
                 m_lptr.data(), 
                 m_lend.data(), 
                 &nrow, 
                 &m_numTriangles, 
                 triangleIndicesNeighbours.data(), 
                 &error );
        
        if(error<0)
        {
            cout << "Error encountered while finding triangle-indices: " << error << endl;
            exit(EXIT_FAILURE);
        }
        
        /*-----------------------------------------------------------------------------
         * Store triangle-indices and find node-neighbours. 
         * NOTE: Fortran indices start from 1.
         *-----------------------------------------------------------------------------*/
        m_triangleIndices.resize(m_numTriangles*3);
        for(int i=0; i<m_numTriangles; i++)
        {
            int idx1 = triangleIndicesNeighbours[i*6]-1;
            int idx2 = triangleIndicesNeighbours[i*6+1]-1;
            int idx3 = triangleIndicesNeighbours[i*6+2]-1;

            m_triangleIndices[i*3]   = idx1;
            m_triangleIndices[i*3+1] = idx2;
            m_triangleIndices[i*3+2] = idx3;

            m_nodeNeighbours[idx1]->insert(idx2);
            m_nodeNeighbours[idx1]->insert(idx3);
            
            m_nodeNeighbours[idx2]->insert(idx1);
            m_nodeNeighbours[idx2]->insert(idx3);
            
            m_nodeNeighbours[idx3]->insert(idx1);
            m_nodeNeighbours[idx3]->insert(idx2);
        }
        
        /*for(int i=0; i<m_shellNumNodes; i++)
        {
            cout << "\t\t node:" << i << endl;
            for( set<int>::iterator nit=m_nodeNeighbours[i]->begin(); 
                 nit!=m_nodeNeighbours[i]->end(); nit++ )
            {
                float dx = m_xs[i]-m_xs[*nit];
                float dy = m_ys[i]-m_ys[*nit];
                float dz = m_zs[i]-m_zs[*nit];

                cout << "\t\t\t neighbour: " << *nit << " dist(" << sqrt(dx*dx+dy*dy+dz*dz) << ")" << endl;
            }
        }*/

        cout << "\tDone.. " << endl;
    }

    Geometry::~Geometry()
    {
        //m_tree->Print();
        delete m_tree;
        delete m_shellTree;

        if(m_initializeTriangulation)
        {
            for(int i=0; i<m_shellNumNodes; i++) delete m_nodeNeighbours[i];
        }
    }
}}
