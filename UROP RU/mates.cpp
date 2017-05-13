//
//  mates.cpp
//  RU_IE_Project
//
//  Created by María B Jónsdóttir on 3/5/17.
//  Copyright © 2017 María B Jónsdóttir. All rights reserved.
//

#include "mates.hpp"

mates::mates(seqan::BamAlignmentRecord unm, seqan::BamAlignmentRecord map)
{
    unmapped = unm;
    mapped = map;
}