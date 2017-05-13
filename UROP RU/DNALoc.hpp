//
//  DNALoc.hpp
//  RU_IE_Project
//
//  Created by María B Jónsdóttir on 9/1/16.
//  Copyright © 2016 María B Jónsdóttir. All rights reserved.
//

#ifndef DNALoc_hpp
#define DNALoc_hpp

#define SEQAN_HAS_ZLIB 1

#include <stdio.h>
#include <seqan/bam_io.h>
#include <zlib.h>
#include <vector>
#include <iostream>
#include <utility>
#include <map>
#include <deque>
#include <unordered_map>
#include "interval.hpp"
#include "mates.hpp"

using namespace std;

class DNALoc
{
    
public:
    DNALoc(int noOfMates, int lengthOfInterval, string velv);
    vector<interval> intervals;
    deque<seqan::BamAlignmentRecord> records;
    vector<mates> listOfMates;
    
    bool checkIfInterval();
    void addUnMapped(seqan::BamAlignmentRecord unmapped, seqan::BamAlignmentRecord mapped);
    void writePairsToFile(string filename);
    void writeIntervalsToFile();

private:
    int countUnMapped;
    int k;
    int lengthOfLoc;
    string velvet;
    
    interval setIntervalObj();
    bool checkOverLimit();
    int getUnmappedCount(int current, int begin);
    void findTotalLeftRightMates(int &right, int &left);
};

#endif /* DNALoc_hpp */


