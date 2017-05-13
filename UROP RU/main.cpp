//
//  main.cpp
//  RU_IE_Project
//
//  Created by María B Jónsdóttir on 8/8/16.
//  Copyright © 2016 María B Jónsdóttir. All rights reserved.
//

#define SEQAN_HAS_ZLIB 1
#define SEQAN_HAS_BZIP2 1

#include <seqan/bam_io.h>
#include <zlib.h>
#include <bzlib.h>
#include "bamread.hpp"

using namespace std;

int main(int argc, char **argv)
{
    if(argc != 6)
    {
        cout << "Missing arguments!" << endl;
        return 1;
    }
    
    int intervalLength = atoi(argv[1]);
    int mappingQuality = atoi(argv[2]);
    int unmappedMates = atoi(argv[3]);
    string velvet = argv[4];
    string path = argv[5];

    bamread bam = bamread(path, intervalLength, mappingQuality, unmappedMates, velvet);
    bam.readAllFile();
    
    return 0;
}