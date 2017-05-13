//
//  bamread.hpp
//  RU_IE_Project
//
//  Created by María B Jónsdóttir on 8/8/16.
//  Copyright © 2016 María B Jónsdóttir. All rights reserved.
//

#ifndef bamread_hpp
#define bamread_hpp

#define SEQAN_HAS_ZLIB 1

#include <stdio.h>
#include <seqan/bam_io.h>
#include <seqan/store.h>
#include <zlib.h>
#include <vector>
#include "DNALoc.hpp" 
#include <array>
#include <fstream>
#include <iomanip>
#include <string>
#include <unordered_map>
#include <utility>
#include <map>
#include "mates.hpp"

using namespace std;

class bamread
{
public:
    bamread(string path, int intervalLen, int mapQual, int noOfmates, string velv);
    void readAllFile();
    void readFileHeader();
    void readFileRecord();
private:
    DNALoc inter;
    string pathToFile;
    int mapQuality;
    int intervalLength;
    int intervals;
    vector<seqan::BamAlignmentRecord> unmapped;
    vector<int> position;
    unordered_map<string, seqan::BamAlignmentRecord> mapped;

    string getKey(seqan::BamAlignmentRecord record);
    void printIntervalsWithIntersections(seqan::BamAlignmentRecord record);
    void printLastIntervalsWithIntersections();
    void eraser(seqan::BamAlignmentRecord record, size_t& begPos);
    void unmappedOrMapped(seqan::BamAlignmentRecord record);
    bool checkIfSoftClipped(seqan::BamAlignmentRecord record);
    void findMate(int i);
    void pairMatesAndCheckIfInterval();
    void readAll();
    void listOfUnmapped(seqan::BamAlignmentRecord record);
    void listOfMapped(seqan::BamAlignmentRecord record);
    void writeToFile(vector<seqan::BamAlignmentRecord> &listName, string filename);
    void storeRecords(seqan::BamAlignmentRecord record);
    void eraseUnmapped(size_t recordEndPos);
    void eraseMates(size_t recordEndPos);
    void eraseRecordsAndMapped(size_t recordEndPos, size_t& begPos, seqan::BamAlignmentRecord record);
    void finishChromosome();
};

#endif /* bamread_hpp */

