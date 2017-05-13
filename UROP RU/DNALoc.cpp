//
//  DNALoc.cpp
//  RU_IE_Project
//
//  Created by María B Jónsdóttir on 9/1/16.
//  Copyright © 2016 María B Jónsdóttir. All rights reserved.
//

#include "DNALoc.hpp"

DNALoc::DNALoc(int noOfMates, int lengthOfInterval, string velv) {
    k = noOfMates;
    countUnMapped = 0;
    lengthOfLoc = lengthOfInterval;
    velvet = velv;
}

void DNALoc::addUnMapped(seqan::BamAlignmentRecord unmapped, seqan::BamAlignmentRecord mapped) {
    mates mate = mates(unmapped, mapped);
    listOfMates.push_back(mate);
    return;
}

                /********************* FINDING INTERVALS USING TERMS *********************/

bool DNALoc::checkIfInterval() {
    if(listOfMates.size() >= k) {
        interval inter = setIntervalObj();
        intervals.push_back(inter);
        return true;
    }
    return false;
}

int DNALoc::getUnmappedCount(int current, int begin) {
    return current - begin;
}

bool DNALoc::checkOverLimit() {
    if(countUnMapped < k) {
        return false;
    }
    else {
        return true;
    }
    return false;
}

                    /********************* END OF FINDING INTERVALS USING TERMS *********************/
                    /********************* FINDING AND SETTING INTERVAL FEATURES *********************/

interval DNALoc::setIntervalObj() {
    interval inter = interval();
    int left = 0;
    int right = 0;
    
    //Determing the total no. of unmapped left mates and unmapped right mates
    findTotalLeftRightMates(right, left);
    
    inter.setBeginPos(records.at(0).beginPos);
    inter.setEndPos(records.at(0).beginPos + lengthOfLoc);
    inter.setNoOfReadsUnmappedMates(listOfMates.size());
    inter.setNoOfReadsLeftUnmapped(left);
    inter.setNoOfReadsRightUnmapped(right);
    inter.setNoOfReads(records.size());
    inter.setChromosome(listOfMates.at(0).mapped.rID + 1);
    
    if(velvet == "unmapped") {
        //print all unmapped reads
    }
    else if(velvet == "all")
    {
        //print all reads in interval
    }
    return inter;
}

void DNALoc::findTotalLeftRightMates(int &right, int &left) {
    for(int m = 0; m < listOfMates.size(); m++) {
        //0x0040 This fragment is the first one in its template
        if(hasFlagFirst(listOfMates.at(m).unmapped)) {
            left++;
        }
        //0x0080 This fragment is the second in its template
        else {
            right++;
        }
    }
    return;
}

                    /********************* END OF FINDING AND SETTING INTERVAL FEATURES *********************/
                    /********************* WRITING/PRINTING OUT INTERVALS AND MATES *********************/

void DNALoc::writeIntervalsToFile() {
    ofstream myfile;
    myfile.open ("./Intervals.txt");
    myfile << intervals.size() << " Intervals found\n\n";
    for(int i = 0; i < intervals.size(); i++) {
        myfile << i+1 << ".\nBeginPos:\t\t\t\t" << intervals.at(i).beginPos << endl << "EndPos:\t\t\t\t\t" << intervals.at(i).endPos << endl << "Reads within interval:\t\t\t" << intervals.at(i).noOfReads << endl << "Reads with unmapped mates:\t\t" << intervals.at(i).noOfReadsWithUnmappedMates << endl << "Intervals with left unmapped mates:\t" << intervals.at(i).noOfReadsWithLeftUnmappedMates << endl << "Reads with right unmapped mates:\t" << intervals.at(i).noOfReadsWithRightUnmappedMates << endl << "No. of intersections:\t\t\t" << intervals.at(i).noOfIntersections << endl << endl;
    }
    return;
}

void DNALoc::writePairsToFile(string filename) {
    ofstream myfile;
    string filePath = "./" + filename + ".txt";
    myfile.open (filePath);
    
    if(myfile.is_open()) {
        myfile << "\tUNMAPPED\t\t\t\t\t\t\t\t\t\tMAPPED\n";
        myfile << "\tQNAME\t\t\tPOS\t\tCIGAR\t\tRNEXT\tPNEXT\t\t\tQNAME\t\t\tPOS\t\tCIGAR\t\tRNEXT\tPNEXT\n";
        for(int i = 0; i < listOfMates.size(); i++) {
            myfile << i << "\t";
            seqan::BamAlignmentRecord record;
            for(int l = 0; l < 2; l++) {
                if(l == 0) {
                    record = listOfMates.at(i).unmapped;
                }
                else {
                    record = listOfMates.at(i).mapped;
                }
                string cigarString = "";
                if(record.cigar.data_begin == NULL) {
                    cigarString = "*";
                }
                else {
                    int count = 0;
                    int m = 0;
                    while(count < length(record.cigar)) {
                        count += record.cigar[m].count;
                        cigarString += to_string(record.cigar[m].count) + record.cigar[m].operation;
                        m++;
                    }
                    
                }
                myfile << setw(20) <<left << toCString(record.qName);;
                myfile << setw(0) << "\t" << setw(8) << left << record.beginPos << setw(0) << "\t" << setw(10) << left << cigarString << "\t" << to_string(record.rNextId) << "\t" << setw(10) << record.pNext << "\t\t";
                
            }
            myfile << "\n";
        }
        myfile.close();
    }
    else {
        cout << "unable to open file";
    }
    return;
}

                    /********************* END OF WRITING/PRINTING OUT INTERVALS AND MATES *********************/
