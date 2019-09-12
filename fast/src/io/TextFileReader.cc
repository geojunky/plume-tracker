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
 *       Filename:  TextFileReader.cc
 *
 *    Description:  TextFileReader
 *
 *        Version:  1.0
 *        Created:  19/07/14 15:56:58
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Rakib Hassan (rakib.hassan@sydney.edu.au)          
 *
 * =====================================================================================
 */
#include <TextFileReader.hh>
#include <pthread.h>

namespace src { namespace io {
   
    TextFileReader::TextFileReader(string fileName, int headerLineCount)
    :m_fileName(fileName),
    m_headerLineCount(headerLineCount)
    {
        mapFile();
    }

    TextFileReader::~TextFileReader()
    {

        if(munmap(m_addr, m_fileLength)==-1)
        {
            handleError("munmap");
        }

        if(close(m_fd)==-1)
        {
            handleError("close");
        }
    }
    
    void TextFileReader::handleError(const char* msg) 
    {
        cout << "Error in :" << msg << endl; 
        exit(0);
    }   

    void TextFileReader::mapFile()
    {
        m_fd = open(m_fileName.c_str(), O_RDONLY);
        if (m_fd == -1)
        {
            handleError("open");
        }   

        posix_fadvise(m_fd, 0, 0, 1);  // FDADVICE_SEQUENTIAL

        // obtain file size
        struct stat sb;
        if (fstat(m_fd, &sb) == -1)
        {
            handleError("fstat");
        }

        m_fileLength = sb.st_size;
        m_addr = (char*) mmap(NULL, m_fileLength, PROT_READ, MAP_PRIVATE, m_fd, 0u);
        if(m_addr == MAP_FAILED)
        {
            handleError("mmap");
        }
    }

    void TextFileReader::read(vector<int> columnMap, vector< vector<float>* > *result)
    {
        int columnMapSize = columnMap.size();
        vector<char> buff(1024);
        /*-----------------------------------------------------------------------------
         * Allocate vectors 
         *-----------------------------------------------------------------------------*/
        result->clear();
        for (vector<int>::iterator it=columnMap.begin(); it!= columnMap.end(); it++)
        {
            result->push_back(new vector<float>);
        }


        int headerLineCount = 0;
        size_t numLines = 0;

        char *f = m_addr;
        char *l = f + m_fileLength;
        while (f && f!=l)
        {
            char *start = f;
            char *end = f;
            if ((f = (char*)(memchr(f, '\n', l-f))))
            {
                end = f;
                numLines++;
                f++;
            }
            
            if(headerLineCount<m_headerLineCount)
            {
                headerLineCount++;
                continue;
            }

            /*-----------------------------------------------------------------------------
             * Parse values 
             *-----------------------------------------------------------------------------*/
            {
                /* copy line */
                copy(start, end+1, buff.begin());
                vector< vector<float>* > &resultRef = *result;
                
                int tokCount = 0;
                char *tok = strtok(buff.data(), " ");
                while(tok != NULL)
                {
                    if(tokCount<columnMapSize)
                    {
                        if(columnMap[tokCount])
                        {
                            resultRef[tokCount]->push_back(atof(tok));
                        }
                    }

                    tok = strtok(NULL, " ");
                    tokCount++;
                }

                if(tokCount != columnMapSize) handleError("Input column-count != column-count on file");
            }
        }
        
        //std::cout << "numLines = " << numLines << "\n";
    }
}}

