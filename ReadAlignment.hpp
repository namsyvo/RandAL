/*
 *  Copyright (c) 2012 Nam Sy Vo.
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1. Redistributions of source code must retain the above Copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above Copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   3. Neither the name of the authors nor the names of its contributors
 *      may be used to endorse or promote products derived from this
 *      software without specific prior written permission.
 */

#ifndef _READ_ALIGNMENT_HPP_
#define _READ_ALIGNMENT_HPP_

#include "fmIndex.hpp"

using namespace std;

class ReadAlignment {

    vector<uint8_t> Ref; // reference string
    size_t Ref_Size; // size of reference string
    
    uint16_t Dist_Threshold;
    uint16_t Attempt_Number;
    uint16_t Repeat_Number;
    uint16_t Dist_Mode;
    uint16_t Min_Window_Size;

    fmIndex fm_backward,fm_forward;
    wat_array::WatArray wa_backward, wa_forward;
    vector<uint32_t> cf_backward, cf_forward;
    size_t wa_forward_len, wa_backward_len;
    
    void printVector(vector<uint8_t> v, ostream& os);
    void printVector(vector<uint8_t> v);
    
    inline size_t getIndexInOriginal(size_t index);
    inline uint16_t min3(uint16_t a, uint16_t b, uint16_t c);
    
    //uint16_t Hamming_dist(const vector<uint8_t> &s, const vector<uint8_t> &t);
    inline uint16_t Edit_dist(const vector<uint8_t> &s, const vector<uint8_t> &t);
    
    /*
    int getMatchHamming(const vector<uint8_t> &left_read, const vector<uint8_t> &right_read,
    size_t leftIndex, size_t rightIndex, uint16_t leftStep, uint16_t rightStep,
    vector<size_t> &matchIndexes, vector<uint16_t> &matchDistances);
    */
    inline int fm_backward_search(const vector<uint8_t> &read, size_t &b_sp, size_t &b_ep, uint16_t &b_step);
    inline int fm_forward_search(const vector<uint8_t> &read, size_t &f_sp, size_t &f_ep, uint16_t &f_step);

    int getMatchRandom(const vector<uint8_t> &left_read, const vector<uint8_t> &right_read,
                    size_t left_sp, size_t left_ep, size_t right_sp, size_t right_ep, uint16_t leftStep, uint16_t rightStep,
                    vector<size_t> &matchIndexes, vector<uint16_t> &matchDistances);
    
    int getMatch(const vector<uint8_t> &read, size_t left_sp, size_t left_ep, size_t right_sp, size_t right_ep, 
                        uint16_t leftStep, uint16_t rightStep, vector<size_t> &matchIndexes, vector<uint16_t> &matchDistances);
	
    inline int LF_mapping(size_t leftIndex, size_t rightIndex, uint16_t left_mate_flank_size, uint16_t right_mate_flank_size,
                    vector<uint8_t> &left_mate_flank, vector<uint8_t> &right_mate_flank);
    
public:
    int Init(int mode, int rep, int th, int att, int wds, const char *backwardIndex, const char *forwardIndex, const char *refFile);
    int AlignRead(string read, vector<size_t> &matchIndexes, vector<uint16_t> &matchDistances);
};


#endif // _READ_ALIGNMENT_H_
