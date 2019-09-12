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
 *       Filename:  Geometry.hh
 *
 *    Description:  Implements output routines
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
#ifndef SRC_MESH_GEOMETRY_HH
#define SRC_MESH_GEOMETRY_HH

#include <KdTree.hh>
#include <algorithm>
#include <set>
#include <map>

namespace src { namespace mesh {
using namespace std;
    const double PI = 3.141592653589793238463;
    const int NCAP = 12;

    class TimeSlices;
    class ChemicalPile;
    
    typedef struct GeoLocation_t
    {
        float lon;
        float lat;
        int   globalId;
    }GeoLocation;
    struct GeoCompareStruct
    {
        bool operator() (const GeoLocation &a, const GeoLocation &b)
        {
            if( (a.lon < b.lon) ||
                ((a.lon == b.lon) && (a.lat < b.lat)) )
            {
                return true;
            }

            return false;
        }
    };
    
    class Geometry
    {
        public:


        friend class TimeSlices;
        friend class ChemicalPile;
        Geometry(string fBaseName, float rScale, bool initializeTriangulation);
        ~Geometry();
        static void rtp2xyz(float *rtp, float *xyz);
        static void xyz2rtp(float *xyz, float *rtp);

        void queryBallPoint(float *pos, float r, vector<float> *distance, vector<int> *id);
        void queryShellPoint(float *tp, float lateralRadius, vector<float> *distance, vector<int> *id);
        void initTriangulation();
        float getDistance(int nid1, int nid2);
        int getClosestShellRadiusIndex(float r, double *closestR); 

        private:
        bool m_initializeTriangulation;
        void rtp2xyz();
        string m_fileBaseName;
        KdTree *m_tree;
        KdTree *m_shellTree; // For radial shells
        
        static int myrandom (int i){ return rand()%i; }

        public: 
        vector<float>  m_r;
        set<float>     m_ur; // Unique r-values
        map<float,int> m_uri; // Unique r-value idx
        vector<float>  m_url; // Unique r-value list
        vector<float>  m_theta;
        vector<float>  m_phi;
        vector<float>  m_rNormalized;
        
        vector<float>  m_x;
        vector<float>  m_y;
        vector<float>  m_z;

        float m_minR;
        float m_maxR;
        float m_rScale;

        /*-----------------------------------------------------------------------------
         * Radial Shells from CMB (id 0) to Surface (id n)
         *-----------------------------------------------------------------------------*/
        map< int, set<GeoLocation, GeoCompareStruct> > m_shells;

        /*-----------------------------------------------------------------------------
         * Triangulation related
         *-----------------------------------------------------------------------------*/
        vector<double> m_slons; // Shell lons; shuffled and ranges between [0, 2pi]
        vector<double> m_slats; // Shell lats; shuffled and ranges between [-pi/2, pi/2]

        vector<double> m_sxs;
        vector<double> m_sys;
        vector<double> m_szs;
        
        int m_shellNumNodes;

        vector<int> m_list;
        vector<int> m_lptr;
        vector<int> m_lend;
        vector<int> m_near;
        vector<int> m_next;
        vector<double> m_dist;        
        int m_numTriangles;
        
        vector< set<int>* > m_nodeNeighbours;
        vector<int> m_triangleIndices;
        
        vector<int> m_order;
    };
}}
#endif

