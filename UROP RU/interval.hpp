//
//  interval.hpp
//  RU_IE_Project
//
//  Created by María B Jónsdóttir on 2/22/17.
//  Copyright © 2017 María B Jónsdóttir. All rights reserved.
//

#ifndef interval_hpp
#define interval_hpp

#include <stdio.h>
#include <vector>
#include <stdio.h>
#include <seqan/bam_io.h>
#include <vector>

using namespace std;

class interval
{
public:
    interval();
    int chromosome;
    int beginPos;
    int endPos;
    int noOfReadsWithLeftUnmappedMates;
    int noOfReadsWithRightUnmappedMates;
    int noOfIntersections;
    size_t noOfReads;
    size_t noOfReadsWithUnmappedMates;

    void setBeginPos(int bp);
    void setEndPos(int ep);
    void setNoOfReads(size_t reads);
    void setNoOfReadsUnmappedMates(size_t noR);
    void setNoOfReadsLeftUnmapped(int noR);
    void setNoOfReadsRightUnmapped(int noR);
    void setChromosome(int chr);
    void pushRecordsInInterval(vector<seqan::BamAlignmentRecord> readsRecords);
    friend ostream& operator<<(ostream& os, const interval& inter);
private:
    
};

#endif /* interval_hpp */
