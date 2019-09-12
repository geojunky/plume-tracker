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
 *       Filename:  TimeSlices.hh
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
#ifndef SRC_MESH_TIMESLICES_HH
#define SRC_MESH_TIMESLICES_HH

#include <Geometry.hh>
#include <map>
namespace src { namespace solver{
    class VelocityExtractor;
}}

namespace src { namespace mesh {
using namespace std;
    
    class ClusterResult
    {
        public:
        vector<int> m_cid;
        vector<int> m_mask;
        vector<float> m_clusterInput;
    };

    class TimeSlices
    {
        public:
        
        friend class src::solver::VelocityExtractor;

        /* Define constants */
        static const int NCLUSTER = 2;

        class TimeStepData
        {
            public:
            vector<float> m_vx;
            vector<float> m_vy;
            vector<float> m_vz;
            vector<float> m_temperature;
            vector<float> m_viscosity;
            vector<float> m_vr;
            vector<float> m_vtheta;
            vector<float> m_vphi;
            map<float, float> m_vrShellAverage;
            map<float, float> m_temperatureShellAverage;
            map<float, float> m_temperatureShellAverageHot; /* average excluding cold material */
            map<float, float> m_viscosityShellAverage;
            vector< vector<float>* > m_opt;
        };

        TimeSlices(string fBaseName, vector<float> times, vector<int> &timeSteps, float velocityScale, Geometry *g,
                   bool processOpt);
        ~TimeSlices();
        
        bool getTimeExtentIndices(float time, int *lowIdx, int *highIdx);
        float getClosestTimeIndex(float time, int *idx);
        
        void getPlumeMetric(float time, float depth, float cutOffPercentile,
                            int *closestTimeIdx, double *closestR, vector<double> *gradVr, 
                            vector<double> *vr, vector<double> *temperature, vector<double> *eta,
                            ClusterResult *cr, int flavourId, vector<double> *composition);

        /* Gradient related */
        void getGradVr(float time, float depth, vector<double> *result, int *closestTimeIdx, double *closestR,
                       vector<double> *vr, vector<double> *temperature, vector<double> *eta);
        
        void getShellVelocities(int timeIndex, int shellIndex, 
                                                   float lateralRadius,
                                                   vector<float> &colats, 
                                                   vector<float> &lons, 
                                                   vector<float> *vr,
                                                   vector<float> *vt,
                                                   vector<float> *vp);

        private:
        void loadData();
        void loadOpt();
        void vrtp2vxyz(vector<float> &vr, vector<float> &vtheta, vector<float> &vphi, 
                       vector<float> &r, vector<float> &theta, vector<float> &phi, 
                       vector<float> &vx, vector<float> &vy, vector<float> &vz);
        
        void vxyz2vrtp(vector<float> &u1, vector<float> &u2, vector<float> &u3, 
                               vector<float> &x1, vector<float> &x2, vector<float> &x3,
                               vector<float> &vr, vector<float> &vtheta, vector<float> &vphi);

        void computeGradient(vector<double> *field, vector<double> *result);
        
        string                  m_fileBaseName;

        public:
        vector<float>           m_times;
        vector<int>             m_timeSteps;
        vector<TimeStepData>    m_timeStepDataList;
        float                   m_velocityScale;
        Geometry                *m_geometry;

        float                   m_minTime;
        float                   m_maxTime;
        bool                    m_processOpt;
        int                     m_hasOpt;
        int                     m_optNCol;
    };
}}
#endif

