/*
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
 */
/*
 * =====================================================================================
 *
 *       Filename:  plumeTrackFast.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  21/07/14 16:22:56
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Rakib Hassan (rakib.hassan@sydney.edu.au)          
 *
 * =====================================================================================
 */
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <exception>
#include <stdexcept>
#include <fstream>
#include <sstream>

#include <math.h>

#include <Allocator.hh>
#include <TextFileReader.hh>
#include <KdTree.hh>
#include <Geometry.hh>
#include <TimeSlices.hh>

#include <utils.hh>

#include <ezOptionParser.hh>
#include <exception>

using namespace std;
using src::mem::PoolAllocator;
using src::io::TextFileReader;
using src::mesh::KdTree;
using src::mesh::Geometry;
using src::mesh::TimeSlices;
using namespace src::mesh;
using namespace ez;

/*
 * =====================================================================================
 *        Class:  ParameterSet
 *  Description:  Class for harvesting parameters from the command-line arguments 
 *                provided.
 * =====================================================================================
 */
class ParameterSet
{
    public:
    ParameterSet(int argc, char **argv)
    {
        m_opt.overview = "plumeTrackFast features:";
        m_opt.syntax = "./plumeTrackFast [REQUIRED OPTIONS] ";
        m_opt.example = "./plumeTrackFast --data-file <gpm19> --data-dir <pathToCapFiles> --time-file <modelTimes.txt> --velocity-scale 4.953209857165279 --radius 6371e3 --start-time <0> --stop-time <230> --output-file-basename <gpm19.plumes> --start-depth 350e3 --validation-depth 1500e3 --cutoff-percentile 5 --tracer-flavour 5 \n\n";
        m_opt.footer = "plumeTrackFast v(0.1), Copyright (C) 2014 Rakib Hassan (rakib.hassan@sydney.edu.au) \nThis program is free and without warranty.\n";

        m_opt.add(
            "", // Default.
            0, // Required?
            0, // Number of args expected.
            0, // Delimiter if expecting multiple args.
            "Display usage instructions.", // Help description.
            "-h",     // Flag token. 
            "-help",  // Flag token.
            "--help", // Flag token.
            "--usage" // Flag token.
        );
        
        m_opt.add(
            "", // Default.
            1, // Required?
            1, // Number of args expected.
            0, // Delimiter if expecting multiple args.
            "Model name.", // Help description.
            "--data-file"     // Flag token. 
        );
        
        m_opt.add(
            "", // Default.
            1, // Required?
            1, // Number of args expected.
            0, // Delimiter if expecting multiple args.
            "Path to available cap-files.", // Help description.
            "--data-dir"     // Flag token. 
        );


        m_opt.add(
            "", // Default.
            1, // Required?
            1, // Number of args expected.
            0, // Delimiter if expecting multiple args.
            "A two-column text file listing model-times (Myr) and model-time-steps corresponding to available cap-files.", // Help description.
            "--time-file"     // Flag token. 
        );
        
        m_opt.add(
            "", // Default.
            1, // Required?
            1, // Number of args expected.
            0, // Delimiter if expecting multiple args.
            "Velocity-scale (m/Myr)", // Help description.
            "--velocity-scale"     // Flag token. 
        );
        
        m_opt.add(
            "", // Default.
            1, // Required?
            1, // Number of args expected.
            0, // Delimiter if expecting multiple args.
            "Radius (m)", // Help description.
            "--radius"     // Flag token. 
        );
        
        m_opt.add(
            "", // Default.
            1, // Required?
            1, // Number of args expected.
            0, // Delimiter if expecting multiple args.
            "Time (Myr) from which the model-space is scanned for plumes", // Help description.
            "--start-time"     // Flag token. 
        );
        
        m_opt.add(
            "", // Default.
            1, // Required?
            1, // Number of args expected.
            0, // Delimiter if expecting multiple args.
            "Time (Myr) from which the model-space is no longer scanned for plumes", // Help description.
            "--stop-time"     // Flag token. 
        );
        
        m_opt.add(
            "", // Default.
            1, // Required?
            1, // Number of args expected.
            0, // Delimiter if expecting multiple args.
            "Depth (m) at which the model-space is scanned for plumes", // Help description.
            "--start-depth"     // Flag token. 
        );
        
        m_opt.add(
            "", // Default.
            1, // Required?
            1, // Number of args expected.
            0, // Delimiter if expecting multiple args.
            "Depth (m) at which the model-space is scanned to validate that plume-conduits found at depth '--start-depth' are indeed plumes", // Help description.
            "--validation-depth"     // Flag token. 
        );
        
        m_opt.add(
            "", // Default.
            1, // Required?
            1, // Number of args expected.
            0, // Delimiter if expecting multiple args.
            "Cut-off percentile below which field values are ignored for cluster analysis.", // Help description.
            "--cutoff-percentile"     // Flag token. 
        );   

        m_opt.add(
            "", // Default.
            1, // Required?
            1, // Number of args expected.
            0, // Delimiter if expecting multiple args.
            "Flavour (1 based) of deep tracers associated with LLSVPs. Set to -1 to ignore tracers.", // Help description.
            "--tracer-flavour"     // Flag token. 
        );        

        m_opt.add(
            "", // Default.
            1, // Required?
            1, // Number of args expected.
            0, // Delimiter if expecting multiple args.
            "Base-name for output files", // Help description.
            "--output-file-basename"     // Flag token. 
        );

        m_opt.parse(argc, const_cast<const char**>(argv));
        if (m_opt.isSet("-h")) 
        {
            Usage();
            exit(EXIT_SUCCESS);
        }

        std::vector<std::string> badOptions;
        if(!m_opt.gotRequired(badOptions)) {
            for(unsigned int i=0; i < badOptions.size(); ++i)
                std::cerr << "ERROR: Missing required option " << badOptions[i] << ".\n\n";
            Usage();
            exit(EXIT_SUCCESS);
        }

        /*-----------------------------------------------------------------------------
         * Extract parameters and set default values
         *-----------------------------------------------------------------------------*/
        m_cutOffPercentile = 5;

        try
        {
            m_opt.get("--data-file")->getString(m_dataFile);
            m_opt.get("--data-dir")->getString(m_dataDir);
            m_opt.get("--time-file")->getString(m_timeFile);
            m_opt.get("--velocity-scale")->getFloat(m_velocityScale);
            m_opt.get("--radius")->getFloat(m_radius);
            m_opt.get("--start-time")->getFloat(m_startTime);
            m_opt.get("--stop-time")->getFloat(m_stopTime);
            m_opt.get("--start-depth")->getFloat(m_startDepth);
            m_opt.get("--validation-depth")->getFloat(m_validationDepth);
            m_opt.get("--cutoff-percentile")->getFloat(m_cutOffPercentile);
            m_opt.get("--tracer-flavour")->getInt(m_tracerFlavour);
            m_opt.get("--output-file-basename")->getString(m_outputFileBaseName);
        }
        catch(exception &e)
        {
            cout << "Exception encountered while processing arguments: " << e.what() << endl;
            exit(EXIT_FAILURE);
        }

        /*-----------------------------------------------------------------------------
         * Derived parameters 
         *-----------------------------------------------------------------------------*/
        m_capFileBaseNameWithPath = m_dataDir + "/" + m_dataFile;
        
        try
        {
            /* Read times and time-steps */
            {
                char buffer[1024] = {0};
                FILE *f = fopen(m_timeFile.c_str(), "r");

                while(fgets(buffer, sizeof(buffer), f)!=NULL)
                {
                    float t = 0;
                    int ts = 0;
                    if(sscanf(buffer, "%f %d", &t, &ts)) 
                    {
                        m_times.push_back(t);
                        m_timeSteps.push_back(ts);
                    }
                }

                fclose(f);
            }                        
        }
        catch(exception &e)
        {
            cout << "Invalid text-files: " << e.what() << endl;
        }
    }

    ~ParameterSet(){}
    
    void Usage()
    {
        string usage;
        m_opt.getUsage(usage);
        cout << usage;
    }

    ezOptionParser  m_opt;

    /* parameters */
    string          m_dataFile;
    string          m_dataDir;
    string          m_timeFile;
    float           m_velocityScale;
    float           m_radius;
    float           m_startTime;
    float           m_stopTime;
    float           m_startDepth;
    float           m_validationDepth;
    float           m_cutOffPercentile;
    int             m_tracerFlavour;
    string          m_outputFileBaseName;

    /* derived parameters */
    struct RTP { float triplet[3]; };
    string          m_capFileBaseNameWithPath;
    vector<float>   m_times;
    vector<int>     m_timeSteps;
};

/*-----------------------------------------------------------------------------
 * Main 
 *-----------------------------------------------------------------------------*/
int main(int argc, char **argv)
{
    /*-----------------------------------------------------------------------------
     * Save arguments 
     *-----------------------------------------------------------------------------*/
    saveArguments(argc, argv);
    
    /*-----------------------------------------------------------------------------
     * Output banner 
     *-----------------------------------------------------------------------------*/
    cout << endl << endl;
    cout << "************************************************************************" << endl;
    cout << "* plumeTrackFast (v 0.1)                                               *" << endl;
    cout << "************************************************************************" << endl << endl;
    
    ParameterSet ps(argc, argv);

    /*-----------------------------------------------------------------------------
     * Report input parameters 
     *-----------------------------------------------------------------------------*/
    cout << "\nInput parameters:" << endl;
    cout <<   "-----------------" << endl << endl;
    
    cout << left << setw(50) << "\tModel name" << ": " << ps.m_dataFile << endl;
    cout << left << setw(50) << "\tCap-file path" << ": " << ps.m_dataDir << endl;
    cout << left << setw(50) << "\tTime-file" << ": " << ps.m_timeFile << endl;
    cout << left << setw(50) << "\tVelocity-scale" << ": " << ps.m_velocityScale << endl;
    cout << left << setw(50) << "\tEarth radius" << ": " << ps.m_radius << endl;
    cout << left << setw(50) << "\tScan start-time" << ": " << ps.m_startTime << endl;
    cout << left << setw(50) << "\tScan stop-time" << ": " << ps.m_stopTime << endl;
    cout << left << setw(50) << "\tScan start-depth" << ": " << ps.m_startDepth << endl;
    cout << left << setw(50) << "\tScan validation-depth" << ": " << ps.m_validationDepth << endl;
    cout << left << setw(50) << "\tCut-off percentile" << ": " << ps.m_cutOffPercentile << endl;
    cout << left << setw(50) << "\tFlavour of tracer to track" << ": " << ps.m_tracerFlavour << endl;
    cout << left << setw(50) << "\tOutput-file base-name" << ": " << ps.m_outputFileBaseName << endl << endl;
    
    /*-----------------------------------------------------------------------------
     * Instantiate objects to be used by multiple threads simultaneously
     *-----------------------------------------------------------------------------*/
    cout << "\nLoading geometry .." << endl;
    Geometry g(ps.m_capFileBaseNameWithPath, ps.m_radius, true);
    
    cout << "\nLoading time-slice data .." << endl;
    TimeSlices ts(ps.m_capFileBaseNameWithPath, ps.m_times, ps.m_timeSteps, ps.m_velocityScale, &g, true);
     
    cout << "\nDetecting plume.." << endl << endl;
    for(unsigned int it=0; it<ts.m_times.size(); it++)
    {
        vector<double> gradVrDepth;
        vector<double> vrDepth;
        vector<double> temperatureDepth;
        vector<double> etaDepth;
        vector<double> compositionDepth;
        double closestRDepth;
        int closestTimeIdxDepth;
        ClusterResult crDepth;
        
        vector<double> gradVrValidationDepth;
        vector<double> vrValidationDepth;
        vector<double> temperatureValidationDepth;
        vector<double> etaValidationDepth;
        vector<double> compositionValidationDepth;
        double closestRValidationDepth;
        int closestTimeIdxValidationDepth;
        ClusterResult crValidationDepth;
        
        
        cout << "\tProcessing time-slice " << ts.m_times[it] << " Myr.." << endl;

        ts.getPlumeMetric(ts.m_times[it], ps.m_startDepth, ps.m_cutOffPercentile,
                          &closestTimeIdxDepth, &closestRDepth, &gradVrDepth, &vrDepth, 
                          &temperatureDepth, &etaDepth, &crDepth, ps.m_tracerFlavour, 
                          &compositionDepth);
        
        ts.getPlumeMetric(ts.m_times[it], ps.m_validationDepth, ps.m_cutOffPercentile,
                          &closestTimeIdxValidationDepth, &closestRValidationDepth, 
                          &gradVrValidationDepth, &vrValidationDepth, 
                          &temperatureValidationDepth, &etaValidationDepth, &crValidationDepth,
                          ps.m_tracerFlavour, &compositionValidationDepth);

        /*-----------------------------------------------------------------------------
         * Output files and cleanup 
         *-----------------------------------------------------------------------------*/
        ostringstream oss;
        ofstream ofs;
        string outputFileName;
        
        int ti = (int)round(ts.m_times[it]);
        oss << ps.m_outputFileBaseName << ".clustered." << ti << ".txt";
        outputFileName = oss.str();
        
        ofs.open(outputFileName.c_str(), ofstream::out);
        {
            for(int i=0; i<g.m_shellNumNodes; i++)
            {
                double magGradVrDepth = sqrt( pow(gradVrDepth[i*3+0], 2.) +
                                              pow(gradVrDepth[i*3+1], 2.) +
                                              pow(gradVrDepth[i*3+2], 2.) );
                
                double magGradVrValidationDepth = sqrt( pow(gradVrValidationDepth[i*3+0], 2.) +
                                                        pow(gradVrValidationDepth[i*3+1], 2.) +
                                                        pow(gradVrValidationDepth[i*3+2], 2.) );
                                                

                double composition = 0.;
                if(ts.m_hasOpt) composition = compositionDepth[i];

                char line[1048] = {0};
                /* r, theta, phi, time, magGradVr, v_r, T, Eta, data cid magGradVrVD, v_rVD, TVD, EtaVD, dataVD cidVD tracComp
                 * 0   1      2    3        4       5   6   7     8   9       10        11   12    13      14     15     16   */

                sprintf(line, "% 08.6e % 08.6e % 08.6e % 08.6e % 08.6e % 08.6e % 08.6e % 08.6e % 08.6e % 08.6e % 08.6e % 08.6e % 08.6e % 08.6e % 08.6e % 08.6e % 08.6e\n",
                        closestRDepth/ps.m_radius,
                        g.m_slats[i]*180./PI, 
                        g.m_slons[i]*180./PI, 
                        ts.m_times[closestTimeIdxDepth],
                        magGradVrDepth,
                        vrDepth[i],
                        temperatureDepth[i],
                        etaDepth[i],
                        (double)crDepth.m_clusterInput[i],
                        (double)crDepth.m_cid[i], 
                        
                        magGradVrValidationDepth,
                        vrValidationDepth[i],
                        temperatureValidationDepth[i],
                        etaValidationDepth[i],
                        (double)crValidationDepth.m_clusterInput[i],
                        (double)crValidationDepth.m_cid[i],
                        composition);
                ofs << line;
            }
        }
        ofs.close();
    }    

    return EXIT_SUCCESS;
}

