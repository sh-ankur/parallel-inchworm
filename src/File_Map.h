//
//  File_Map.h
//  
//
//  Created by Ankur Sharma on 10/11/13.
//
//

#ifndef ____File_Map__
#define ____File_Map__

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <omp.h>
#include <vector>

using namespace std;

class File_Map {
public:
    File_Map(string filename, int parts, int offset);
    int partitions;
    vector <char *> pmap;
    long get_limit(int index);
    bool create_map();
    
private:
    string filename;
    long pages_per_partition;
    long filesize;
    long pagesize;
    int fd, offset;
};

#endif /* defined(____File_Map__) */
