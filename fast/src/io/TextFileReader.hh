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
 *       Filename:  TextFileReader.hh
 *
 *    Description:  Fast text-file reader
 *
 *        Version:  1.0
 *        Created:  19/07/14 15:55:57
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Rakib Hassan (rakib.hassan@sydney.edu.au)          
 *
 * =====================================================================================
 */
#ifndef PROPAT_TEXT_FILE_READER
#define PROPAT_TEXT_FILE_READER

#include <iostream>
#include <vector>

#include <algorithm>
#include <iostream>

// for mmap:
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstring>

using namespace std;
namespace src { namespace io {

    class TextFileReader
    {
        public:
        TextFileReader(string fileName, int headerLineCount);
        ~TextFileReader();

        void read(vector<int> columnMap, vector< vector<float>* > *result);
        
        private:
        
        void mapFile();
        void handleError(const char* msg) ;

        string m_fileName;
        int m_headerLineCount;

        int m_fd;
        char *m_addr;
        size_t m_fileLength;
    };
}}
#endif
