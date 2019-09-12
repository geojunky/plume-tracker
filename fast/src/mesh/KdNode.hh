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
 *       Filename:  KdNode.hh
 *
 *    Description:  Implements output routines
 *
 *        Version:  1.0
 *        Created:  11/06/14 14:02:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Rakib Hassan (rakib.hassan@sydney.edu.au)         
 *
 * =====================================================================================
 */
#ifndef SRC_MESH_KDNODE_HH
#define SRC_MESH_KDNODE_HH

#include <KdItem.hh>
#include <vector>

namespace src { namespace mesh {
using namespace std;

    class KdTree;

    /*
     * =====================================================================================
     *        Class:  KdNode
     *  Description:  KdNode
     * =====================================================================================
     */
    class KdNode
    {
        /*-----------------------------------------------------------------------------
         * Public interface
         *-----------------------------------------------------------------------------*/
        public:
        KdNode( KdTree *tree );

        KdNode( KdNode *node );

        void Add( KdItem *m );

        void Split( KdItem *m );

        void Expand( float *newCoord );

        void Print();

        void DeleteHelper();

        bool Intersects( float *up0, float *low0, float *up1, float *low1 );

        void Range( float *upper, float *lower, vector<KdItem*> *result );
        
        bool Contains( float *upper, float *lower, float *pnt );

        ~KdNode();

        private:
        
        /*-----------------------------------------------------------------------------
         * Private internals
         *-----------------------------------------------------------------------------*/
        void RangeHelper( float *upper, float *lower, vector<KdItem*> *result );
        
        struct KdItemSort
        {
            bool operator()( const KdItem *a, const KdItem *b )
            {
                if(DIM==2)
                {
                    return ( (a->m_coord[0] < b->m_coord[0]) || 
                             ((a->m_coord[0] == b->m_coord[0]) && (a->m_coord[1] < b->m_coord[1])) );
                }
                else if(DIM==3)
                {
                    return ( (a->m_coord[0] < b->m_coord[0]) || 
                             ((a->m_coord[0] == b->m_coord[0]) && (a->m_coord[1] < b->m_coord[1])) || 
                             ((a->m_coord[0] == b->m_coord[0]) && (a->m_coord[1] == b->m_coord[1]) && (a->m_coord[2] < b->m_coord[2])) );
                }

                return false;
            }
        }m_itemSort;

        KdTree *m_owner;
        int m_currDim;
        int m_currSize;

        vector< KdItem* > *m_bucket;

        KdNode *m_left;
        KdNode *m_right;

        float *m_upper;
        float *m_lower;

        float m_slice;
    };
}}
#endif

