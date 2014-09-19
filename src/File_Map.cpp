//
//  File_Map.cpp
//  
//
//  Created by Ankur Sharma on 10/11/13.
//
//

#include "File_Map.h"

File_Map :: File_Map(string filename, int parts, int offset){
    
    this->filename = filename;
    this->partitions = parts;
    this->offset = offset;
    this->pmap.resize(this->partitions);
    this->fd = open(filename.c_str(), O_RDONLY);
    
    if(this->fd == -1){
        perror("File Not Open:");
        exit(1);
    }
    
    struct stat filestat;
    
    if(fstat(this->fd,&filestat)){
        perror("Stat:");
        exit(1);
    }
    
    this->filesize = filestat.st_size;
    
    this->pagesize = getpagesize();
    
    if(this->partitions > 0)
        this->pages_per_partition = (this->filesize / this->pagesize) / this->partitions;
    
    create_map();	
    cerr << "File mapped to memory" << endl;
    
}

long File_Map :: get_limit(int index){
 
    if (index < (this->partitions - 1)) {
        return ((this->pages_per_partition * this->pagesize) + this->offset);
    } else {
        return (this->filesize - (index * this->pages_per_partition * this->pagesize) - 1);
    }
    
    
}

bool File_Map :: create_map(){
    
    omp_set_num_threads(this->partitions);
    long limit;
#pragma omp parallel for
    for(int i = 0 ; i < this->partitions ; i++){

	limit = this->get_limit(i);

        this->pmap[i] = (char *)mmap(0, limit, PROT_READ, MAP_SHARED, this->fd, i * pages_per_partition * this->pagesize);
        if(this->pmap[i] == MAP_FAILED){
            perror("MMAP FAILED:");
            exit(1);
        }
        
    }
}


