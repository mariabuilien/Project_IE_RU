//
//  interval.cpp
//  RU_IE_Project
//
//  Created by María B Jónsdóttir on 2/22/17.
//  Copyright © 2017 María B Jónsdóttir. All rights reserved.
//

#include "interval.hpp"

interval::interval() {
    chromosome = 0;
    beginPos = 0;
    endPos = 0;
    noOfReadsWithLeftUnmappedMates = 0;
    noOfReadsWithRightUnmappedMates = 0;
    noOfIntersections = -1;
    noOfReads = 0;
    noOfReadsWithUnmappedMates = 0;
}

void interval::setBeginPos(int bp) {
    beginPos = bp;
    return;
}

void interval::setEndPos(int ep) {
    endPos = ep;
    return;
}

void interval::setNoOfReads(size_t reads) {
    noOfReads = reads;
    return;
}

void interval::setNoOfReadsUnmappedMates(size_t noR) {
    noOfReadsWithUnmappedMates = noR;
    return;
}

void interval::setNoOfReadsLeftUnmapped(int noR) {
    noOfReadsWithLeftUnmappedMates = noR;
    return;
}

void interval::setNoOfReadsRightUnmapped(int noR) {
    noOfReadsWithRightUnmappedMates = noR;
    return;
}

void interval::setChromosome(int chr) {
    chromosome = chr;
    return;
}

ostream& operator<<(ostream& os, const interval& inter) {
    os << inter.chromosome << " " << inter.beginPos << " " << inter.endPos << " " << inter.noOfReads << " " << inter.noOfReadsWithUnmappedMates << " " << inter.noOfReadsWithLeftUnmappedMates << " " <<inter.noOfReadsWithRightUnmappedMates << " " << inter.noOfIntersections;
    return os;
}
