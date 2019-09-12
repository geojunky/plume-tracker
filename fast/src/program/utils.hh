/*
 * =====================================================================================
 *
 *       Filename:  utils.hh
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  23/02/15 11:47:37
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef UTILS_HH
#define UTILS_HH

#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include <time.h>

using namespace std;

/*-----------------------------------------------------------------------------
 * Get current date and time
 *-----------------------------------------------------------------------------*/
const std::string currentDateTime() 
{
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
}

/*-----------------------------------------------------------------------------
 * Save commandline arguments to a file for reference 
 *-----------------------------------------------------------------------------*/
void saveArguments(int argc, char **argv)
{
    char cwd[4096] = {0};
    ostringstream oss;
    
    if (getcwd(cwd, sizeof(cwd)) != NULL){}
    else
    {
        perror("getcwd() error");
    }
    

    oss << cwd << "/" << ".run_history.txt";
    
    FILE *f = fopen(oss.str().c_str(), "a+");
    ostringstream line;

    line << currentDateTime() << " : ";

    for(int i=0; i<argc; i++)
    {
        line << argv[i] << " ";
    }
    
    fprintf(f, "%s\n\n", line.str().c_str());
    fclose(f);
}

#endif
