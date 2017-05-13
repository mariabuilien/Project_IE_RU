//
//  mates.hpp
//  RU_IE_Project
//
//  Created by María B Jónsdóttir on 3/5/17.
//  Copyright © 2017 María B Jónsdóttir. All rights reserved.
//

#ifndef mates_hpp
#define mates_hpp

#include <stdio.h>
#include <seqan/bam_io.h>

using namespace std;

class mates
{
public:
    mates(seqan::BamAlignmentRecord unmapped, seqan::BamAlignmentRecord mapped);
    seqan::BamAlignmentRecord unmapped;
    seqan::BamAlignmentRecord mapped;
private:
};
#endif /* mates_hpp */
