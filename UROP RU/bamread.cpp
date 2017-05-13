//
//  bamread.cpp
//  RU_IE_Project
//
//  Created by María B Jónsdóttir on 8/8/16.
//  Copyright © 2016 María B Jónsdóttir. All rights reserved.
//

#include "bamread.hpp"

int LENGTH_READ = 151;

bamread::bamread(string path, int intervalLen, int mapQual, int noOfmates, string velv) : inter(5, 600, "all") {
    pathToFile = path;
    mapQuality = mapQual;
    intervalLength = intervalLen;
    intervals = 0;
    inter = DNALoc(noOfmates, intervalLen, velv);
}

void bamread::readAllFile() {
    cout << "Here comes the output of the bam file from path " << pathToFile << endl << endl;
    readFileRecord();
    return;
}

                /********************* READ FROM BAM FILE *********************/

void bamread::readAll() {
    readFileHeader();
    readFileRecord();
    return;
}

void bamread::readFileHeader() {
    seqan::BamFileIn bamFile; //from seqan.de
    if (!open(bamFile, seqan::toCString(pathToFile))) { //from seqan.de
        std::cout << "Error! Could not open bam file" << std::endl;
        return;
    }
    // Open output file - from seqan.de
    seqan::BamFileOut bamFileOut(context(bamFile), cout, seqan::Sam());
    
    //from seqan.de
    try {
        // Copy header.
        seqan::BamHeader header;
        readHeader(header, bamFile);
        //writeHeader(bamFileOut, header);
        // Copy records.
    }
    
    catch (seqan::Exception const & e) {
        std::cout << "ERROR: " << e.what() << endl;
        return;
    }
    //from seqan.de
    return;
};

void bamread::readFileRecord() {
    int currentChromosome = 0;
    seqan::BamFileIn bamFile; //from seqan.de
    if (!open(bamFile, seqan::toCString(pathToFile))) { //from seqan.de
        cout << "Error! Could not open bam file" << endl;
        return;
    }
    // Open output file - from seqan.de
    seqan::BamFileOut bamFileOut(context(bamFile), std::cout, seqan::Sam());
    
    //from seqan.de
    try {
        seqan::BamHeader header;
        seqan::BamAlignmentRecord record;
        readHeader(header, bamFile);
        
        size_t begPos = 0; //begin position of the current interval

        clock_t start;
        double duration;
        start = clock();
        cout << "Start reading BAM file" << endl;
        
        while(!atEnd(bamFile)) {
            readRecord(record, bamFile);
            //writeRecord(bamFileOut, record);

            //check if record is in a new chromosome. If so, finish checking out the previous one and clean up datasets before the new one.
            if(record.rID + 1 != currentChromosome)
            {
                finishChromosome();
                currentChromosome = record.rID + 1;
            }
            
            if(inter.records.size() == 0) {
                begPos = record.beginPos;
            }
            else {
                //if current record does not fit in within the current interval, we start pairing mates for that interval
                //and clean up records from memory we don't have to use anymore
                if(record.beginPos + LENGTH_READ > begPos + intervalLength) {
                    pairMatesAndCheckIfInterval();
                    eraser(record, begPos);
                }
            }
            
            printIntervalsWithIntersections(record); //printing intervals we are sure that won't have any more intersection given the current record
            storeRecords(record);
            unmappedOrMapped(record);
        }
        
        //finish up the last interval although it might be shorter in bp and printing out the intervals left to be printed
        pairMatesAndCheckIfInterval();
        printLastIntervalsWithIntersections();
        
        duration = (clock() - start) / (double) CLOCKS_PER_SEC;
        cout << "\nTime reading BAM file: " << duration << " sec\nIntervals found: " << intervals << endl;
    }
    catch (seqan::Exception const & e) {
        cout << "ERROR: " << e.what() << endl;
        return;
    }
    //from seqan.de
    return;
}

void bamread::printIntervalsWithIntersections(seqan::BamAlignmentRecord record) {
    //printing out all intervals that could not intersect with an interval containing the current record
    while(inter.intervals.size() > 0 && inter.intervals.at(0).endPos < (record.beginPos + LENGTH_READ - intervalLength)) {
        int l = 0;
        //iterating through the position vector to check for intersections of an interval that will be printed
        while(l < position.size() && position.at(l) < inter.intervals.at(0).endPos) {
            if(position.at(l) < inter.intervals.at(0).endPos && position.at(l) + intervalLength > inter.intervals.at(0).beginPos) {
                inter.intervals.at(0).noOfIntersections++;
            }
            l++;
        }
    
        cout << inter.intervals.at(0) << endl;
        inter.intervals.erase(inter.intervals.begin());
        
        while(position.size() > 0 && inter.intervals.size() > 0 && position.at(0) + intervalLength < inter.intervals.at(0).beginPos) {
            position.erase(position.begin());
        }
    }
    
    return;
}

void bamread::printLastIntervalsWithIntersections() {
    //print out all intervals left for chromosome
    while(inter.intervals.size() > 0) {
        int l = 0;
        while(l < position.size() && position.at(l) < inter.intervals.at(0).endPos) {
            if(position.at(l) < inter.intervals.at(0).endPos && position.at(l) + intervalLength > inter.intervals.at(0).beginPos) {
                inter.intervals.at(0).noOfIntersections++;
            }
            l++;
        }
        
        cout << inter.intervals.at(0) << endl;
        inter.intervals.erase(inter.intervals.begin());
        //remove from position vector, positions that are no longer needed for intersection information for intervals left or to come
        while(position.size() > 0 && inter.intervals.size() > 0 && position.at(0) + intervalLength < inter.intervals.at(0).beginPos) {
            position.erase(position.begin());
        }
    }
    
    return;
}

bool bamread::checkIfSoftClipped(seqan::BamAlignmentRecord record) {
    bool softClipped = false;
    int i = 0;
   
    while(i < length(record.cigar)) {
        if(record.cigar[i].operation == 'S') {
            softClipped = true;
        }
        i++;
    }
    return softClipped;
}

void bamread::eraser(seqan::BamAlignmentRecord record, size_t& begPos) {
    size_t recordEndPos = record.beginPos + LENGTH_READ;
    
    //erasing records and mapped records we don't use anymore
    eraseRecordsAndMapped(recordEndPos, begPos, record);
    
    //erasing mates we don't use anymore
    eraseMates(recordEndPos);
    
    //erase unmapped we don't use anyore
    eraseUnmapped(recordEndPos);
    return;
}

void bamread::finishChromosome() {
    pairMatesAndCheckIfInterval();
    printLastIntervalsWithIntersections();
    unmapped.clear();
    mapped.clear();
    inter.records.clear();
    position.clear();
    return;
}
                    /********************* END READ FROM BAM FILE *********************/
                    /********************* INSERT RECORDS TO LISTS *********************/

void bamread::unmappedOrMapped(seqan::BamAlignmentRecord record) {
    //check if alignment record is unmapped
    if (hasFlagUnmapped(record)) {
        listOfUnmapped(record);
    }
    //if alignment record is mapped we check if it soft clipped or not
    else {
        //if alignment record is not soft clipped and has enough mapping quality, we add it to the mapped vector
        if(!checkIfSoftClipped(record) && static_cast<int>(record.mapQ) >= mapQuality) {
            listOfMapped(record);
        }
    }
    return;
}

void bamread::storeRecords(seqan::BamAlignmentRecord record) {
    inter.records.push_back(record);
    return;
}

void bamread::listOfUnmapped(seqan::BamAlignmentRecord record) {
    unmapped.push_back(record);
    return;
}

void bamread::listOfMapped(seqan::BamAlignmentRecord record) {
    string key = getKey(record);
    mapped[key] = record;
    return;
}

string bamread::getKey(seqan::BamAlignmentRecord record) {
    string key = to_string(record.beginPos) + " ";
    key += toCString(record.qName);
    return key;
}
                    /********************* END INSERT RECORDS TO LIST *********************/
                    /********************* FIND MATES - UNMAPPED/MAPPED *********************/

void bamread::pairMatesAndCheckIfInterval() {
    inter.listOfMates.clear();
    
    for(int i = 0; i < unmapped.size(); i++) {
        findMate(i);
    }
   
    //inter.writePairToFile("mates");
    if(inter.checkIfInterval()) {
        //Store the begin position of the leftmost record within interval
        position.push_back(inter.intervals.at(inter.intervals.size() - 1).beginPos);
        intervals++;
    }
    
    return;
}

void bamread::findMate(int i) {
    seqan::BamAlignmentRecord record = unmapped.at(i);
    string unmappedKey = getKey(record);
    
    unordered_map<string, seqan::BamAlignmentRecord>::iterator itMapped;
    itMapped = mapped.find(unmappedKey);
    if(itMapped != mapped.end()) {
        inter.addUnMapped(unmapped.at(i), itMapped->second);
    }
    return;
}
                    /********************* END FINDING MATES - UNMAPPED/MAPPED *********************/
                    /********************* ERASE RECORDS THAT WE DON'T USE ANYMORE *********************/

void bamread::eraseUnmapped(size_t recordEndPos) {
    int i = 0;
    if(unmapped.size() > i) {
        size_t unmPos = unmapped.at(i).beginPos;
        //delete all unmapped records that are not within the next interval being looked at
        while(recordEndPos > unmPos + intervalLength) {
            i++;
            if(unmapped.size() > i) {
                unmPos = unmapped.at(i).beginPos;
            }
            else {
                unmPos = recordEndPos;
            }
        }
    }
    if(unmapped.size() > 0) {
        unmapped.erase(unmapped.begin(), unmapped.begin() + i);
    }
    
    return;
}

void bamread::eraseMates(size_t recordEndPos) {
    int i = 0;
    if(inter.listOfMates.size() > i) {
        size_t matePos = inter.listOfMates.at(i).unmapped.beginPos;
        //delete all mates that are not within the next interval being looked at
        while(recordEndPos > matePos + intervalLength) {
            i++;
            if(inter.listOfMates.size() > i) {
                matePos = inter.listOfMates.at(i).unmapped.beginPos;
            }
            else {
                matePos = recordEndPos;
            }
        }
    }
    if(inter.listOfMates.size() > 0) {
        inter.listOfMates.erase(inter.listOfMates.begin(), inter.listOfMates.begin() + 1);
    }
    
    return;
}

void bamread::eraseRecordsAndMapped(size_t recordEndPos, size_t& begPos, seqan::BamAlignmentRecord record) {
    //delete all records and mapped records that are not within the next interval being looked at
    int i = 0;
    while(recordEndPos > begPos + intervalLength) {
        if(inter.records.size() > i) {
            //delete the mapped record of each record if one exists
            if(mapped.size() > 0) {
                string key = getKey(inter.records.at(i));
                
                unordered_map<string, seqan::BamAlignmentRecord>::iterator itMapped;
                itMapped = mapped.find(key);
                if(itMapped != mapped.end()) {
                    mapped.erase(itMapped);
                }
            }
            
            i++;
            if(inter.records.size() > i) {
                begPos = inter.records.at(i).beginPos;
            }
            else {
                begPos = record.beginPos;
            }
        }
    }
    if(inter.records.size() > 0) {
        inter.records.erase(inter.records.begin(), inter.records.begin() + i);
    }
    
    return;
}