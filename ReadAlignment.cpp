/*
 *   Copyright (c) 2012 Memphis-CS Genome Assembly Group
 *   Copyright (c) 2010 Daisuke Okanohara
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

#include "ReadAlignment.hpp"

using namespace std;

/*
 * Initialize FM-index data structure and other information
 *
 */
int ReadAlignment::Init(int mode, int rep, int th, int att, int minws, const char *backwardIndex, const char *forwardIndex, const char *refFile) {

    Dist_Mode = mode;
    Repeat_Number = rep;
    Dist_Threshold = th;
    Attempt_Number = att;
    Min_Window_Size = minws;
     
    //Load FM-Index for Ref from backwardIndex
    ifstream is_b(backwardIndex);
    fm_backward.load(is_b);
    //Load FM-Index for reverse Ref from forwardIndex
    ifstream is_f(forwardIndex);
    fm_forward.load(is_f);
    //Load Reference sequence
    fm_backward.read(refFile, Ref);
    Ref_Size = Ref.size();

    wa_backward = fm_backward.getWa();
    cf_backward = fm_backward.getCF();
    wa_forward = fm_forward.getWa();
    cf_forward = fm_forward.getCF();

    wa_backward_len = wa_backward.length();
    wa_forward_len = wa_forward.length();;

    return 0;
}

/*
 * Get matched index of the query on the original reference from the matched index on reverse reference (from forward search)
 *
 */
inline size_t ReadAlignment::getIndexInOriginal(size_t index){
    return (Ref_Size - 2) - index;
}

/*
 * Get minimum of 3 integer numbers
 *
 */
inline uint16_t ReadAlignment::min3(uint16_t a, uint16_t b, uint16_t c) {
    uint16_t min_ab = (a < b) ? a : b;
	uint16_t min_abc = (min_ab < c) ? min_ab : c;
	return min_abc;
}

/*
uint16_t ReadAlignment::Hamming_dist(const vector<uint8_t> &s, const vector<uint8_t> &t) {
    uint16_t i = 0, H = 0;
    for (i = 0 ; i < s.size(); i++) {
        if (s.at(i) != t.at(i))
            H += 1;
        if(H > Dist_Threshold) {
            break;
        }
    }
    return H;
}
*/

/*
 * Calculating Edit distance with Bound algorithm
 *
 */
inline uint16_t ReadAlignment::Edit_dist(const vector<uint8_t> &s, const vector<uint8_t> &t) {
    uint16_t s_len = s.size();
    uint16_t t_len = t.size();
    uint16_t m  = max(s_len, t_len);
    uint16_t a, b;
    
    vector<uint16_t> M0(m + 1), M1(m + 1);
    
    for (b = 0; b <= t_len; ++b)
        M0[b] = b;
    
    for (a = 1; a <= s_len; ++a)
    {
        M1[0] = a;
        for (b = 1; b <= t_len; ++b)
        {
            uint16_t x = M0[b] + 1;
            uint16_t y = M1[b-1] + 1;
            uint16_t z = M0[b-1] + (s[a-1] == t[b-1] ? 0 : 1);
            M1[b] = min3(x, y, z);
        }
    	// Break 
    	bool flag = true;
    	for (size_t c = 1; c <= t_len; ++c) {
    	    if (M1[c] <= Dist_Threshold)
    	    {
                flag = false;
                break;
    	    }
		}
        if (flag) {
    	   return (Dist_Threshold + 1);
		}
		// Swap two arrays
        swap(M0, M1);
    }
    return M0[t_len];
}

/*
 * Function to align a read
 * Input: read (a string)
 * Output: match positions (a vector) & distance between read and mates (a vector)
 */

int ReadAlignment::AlignRead(string read, vector<size_t> &matchIndexes, vector<uint16_t> &matchDistances) {

    //Start to align read with deterministic search
    string q = read;

    vector< vector<uint8_t> > qsL, qsR;
    fm_backward.buildQstring(q, qsL); //build query for backward search
    fm_forward.buildQstring(q, qsR); //build query for forward search (on reversed reference)
    
    vector<uint8_t> left_read  = qsL[0]; //query for backward search (from the end), encoded from original read
    vector<uint8_t> right_read = qsR[0]; //query for forward search (on reversed reference) (from the beginning)    

    size_t left_sp, left_ep, //sp, ep for backward search from the end (left -> right)
            right_sp, right_ep; //sp, ep for forward search from the beinning (right -> left)
    uint16_t leftStep, rightStep; // length of perfect match from the end and the beginning of read

    fm_backward_search(left_read, left_sp, left_ep, leftStep);
    fm_forward_search(right_read, right_sp, right_ep, rightStep);

    //Computing distance between remain part of query and all extended part of matches
    if (Dist_Mode == 3) {
        int matched = getMatch(left_read, left_sp, left_ep, right_sp, right_ep, leftStep, rightStep, matchIndexes, matchDistances);
        if (matched == 1) {
            return 1;
        }
    }
    //Search with randomized positions
    uint16_t pos;
    for (uint16_t att = 0 ; att < Attempt_Number ; att++) {
        pos = -1;
        while((pos < Repeat_Number) || (pos > read.size() - Repeat_Number - 1)) {
            pos = rand() % ((uint16_t)(read.size()));
        }

        string q = read;
        string qL = q.substr(0, pos);;
        string qR = q.substr(pos, q.size() - pos);

        vector< vector<uint8_t> > qsL, qsR, qsR_ori;
        fm_backward.buildQstring(qL, qsL);
        fm_forward.buildQstring(qR, qsR);
        fm_backward.buildQstring(qR, qsR_ori); // take encoded form for original query (not reversed)

        vector<uint8_t> left_read, right_read, right_read_ori;
        left_read  = qsL[0];
        right_read = qsR[0];
        right_read_ori = qsR_ori[0];

        size_t left_sp, left_ep, right_sp, right_ep; //sp, ep for left and right part (from the position to search) of partial match
        uint16_t leftStep, rightStep; // distance from the position to search to left-most and right-most of partial match

        fm_backward_search(left_read, left_sp, left_ep, leftStep);
        fm_forward_search(right_read, right_sp, right_ep, rightStep);
        if (Dist_Mode == 3) {
            int matched = getMatchRandom(left_read, right_read_ori, left_sp, left_ep, right_sp, right_ep, leftStep, rightStep, matchIndexes, matchDistances);
            if (matched == 1) {
                return 1;
            }
        }
    }
    return 0;
}

/*
 * FM-index-based backward search (search from end of the query)
 * Input: a query as a vector
 * Output: start pointer (b_sp) & end pointer (b_ep) (for fm-index) for maximum perfect match and length of match (b-step)
 */

inline int ReadAlignment::fm_backward_search(const vector<uint8_t> &read, size_t &b_sp, size_t &b_ep, uint16_t &b_step) {

    uint8_t c;    
    size_t sp, ep; // start, end pointer for fm index
    size_t prev_sp, prev_ep; //previous start, end pointer for fmindex
    uint16_t leftStep;
    int16_t idx;

    //Search from the end of read (backward search, from right -> left
    sp = 1;
    ep = wa_backward_len;
    prev_sp = 0;
    prev_ep = 0;
    leftStep = 0;
    idx = read.size() - 1;
    while ((sp < ep) && (idx > -1)) {
        c = read[idx];
        idx--;
        prev_sp = sp; //store previous sp and ep
        prev_ep = ep;
        sp = cf_backward[c] + wa_backward.Rank(c, sp);
        ep = cf_backward[c] + wa_backward.Rank(c, ep);
        leftStep++;
    }
    if (sp < ep) { // found match(es)
        //resL = make_pair(sp, ep - 1);
        b_sp = sp;
        b_ep = ep;
    }
    else { // found only local match(es)
        b_sp = prev_sp; //take previous sp, ep
        b_ep = prev_ep;
        leftStep--; // take previous step length
    }
    b_step = leftStep;
    return 0;
}

/*
 * FM-index-based forward search (search from beginning of the query)
 * Input:   a query as a vector
 * Output:  start pointer (f_sp) & end pointer (f_ep) (for fm-index) for maximum perfect match and length of match (f_step)
 *          on the reverse reference
 */

inline int ReadAlignment::fm_forward_search(const vector<uint8_t> &read, size_t &f_sp, size_t &f_ep, uint16_t &f_step) {
    //Search from the beginning of read (forward search, from left -> right
    uint8_t c;    
    size_t sp, ep; // start, end pointer for fm index
    size_t prev_sp, prev_ep; //previous start, end pointer for fmindex
    uint16_t rightStep;
    int16_t idx;

    sp = 1;
    ep = wa_forward_len;
    prev_sp = 0;
    prev_ep = 0;
    rightStep = 0;
    idx = 0;
    while ((sp < ep) && (idx < read.size())) {
        c = read[idx];
        idx++;
        prev_sp = sp; //store previous sp and ep
        prev_ep = ep;
        sp = cf_forward[c] + wa_forward.Rank(c, sp);
        ep = cf_forward[c] + wa_forward.Rank(c, ep);
        rightStep++;
    }
    if (sp < ep) { // found match(es)
        f_sp = sp;
        f_ep = ep;
    }
    else { // found only local match(es)
        //resR = make_pair(prev_sp, prev_ep - 1);
        f_sp = prev_sp; //take previous sp, ep
        f_ep = prev_ep;
        rightStep--; // take previous step length
    }
    f_step = rightStep;
    return 0;
}

/*
 * Deterministic search strategy: search from the beginning (forward) and the end (backward) of query
 *
 */
int ReadAlignment::getMatch(const vector<uint8_t> &read, size_t left_sp, size_t left_ep, size_t right_sp, size_t right_ep,
                    uint16_t leftStep, uint16_t rightStep, vector<size_t> &matchIndexes, vector<uint16_t> &matchDistances) {

    int match_count = 0;
    int left_count = 0;
    int right_count = 0;

    int16_t li, ri;
    size_t ls, rs;
    size_t leftIndex, rightIndex;

    size_t matchPos; // position for mate (first index on reference sequence)
    uint16_t left_matchDist, right_matchDist, matchDist; // distance for mate

    vector<uint8_t> left_read_flank, right_read_flank; // part of read from the left-most index of read
    vector<uint8_t> left_mate_flank, right_mate_flank; // part of mate from the left-most/right-most index of perfect mate
    vector<uint8_t> read_middle; // part of read from the left-most/right-most index of read
    vector<uint8_t> mate_middle; // part of mate from the left-most/right-most index of perfect mate

    uint16_t left_read_flank_size, right_read_flank_size; //
    uint16_t left_mate_flank_size, right_mate_flank_size; //
    uint16_t read_middle_size; //
    size_t mate_middle_size; //

    //
    if ( leftStep + rightStep < read.size() ) {

      read_middle_size = read.size() - leftStep - rightStep;
      vector<uint8_t> mate_middle_2; // part of mate from the left-most/right-most index of perfect mate
      bool match_flag = false;
      for(ls = left_sp ; ls < left_ep ; ls++) {
        leftIndex = fm_backward.locate(ls);

        for(rs = right_sp ; rs < right_ep ; rs++) {
            rightIndex = getIndexInOriginal(fm_forward.locate(rs));

            //Distance between left and right mate ~ distance between left and right read +/- Dist_Threshold
            mate_middle_size = leftIndex - 1 - rightIndex;
            if ((mate_middle_size <= read_middle_size + Dist_Threshold) && (mate_middle_size >= read_middle_size - Dist_Threshold)) {

                /*
                //Test: skip computing Edit distance
                matchDist = 0;
                matchPos = rightIndex + 1 - rightStep;
                matchIndexes.push_back(matchPos);
                matchDistances.push_back(matchDist);
                match_flag = true;
                */
                int check = LF_mapping(leftIndex, rightIndex, mate_middle_size, mate_middle_size, mate_middle_2, mate_middle);
                if (check != -1) {
                    read_middle.clear();
                    for (ri = rightStep ; ri < rightStep + read_middle_size ; ri++ ) {
                        //cout << (long unsigned int)right_read_ori.at(ri) << "-";
                        read_middle.push_back(read.at(ri));
                    }
                    
                    matchDist = Edit_dist(read_middle, mate_middle);
                    matchPos = rightIndex + 1 - rightStep;
                    if (matchDist <= Dist_Threshold) {
                        matchIndexes.push_back(matchPos);
                        matchDistances.push_back(matchDist);
                        match_flag = true;
		                match_count++;
                    }
                }
            }
        }
      }

      if (match_flag) {
        return 1;
      }
    }

    left_mate_flank_size = read.size() - leftStep;
    left_read_flank_size = read.size() - leftStep;
    bool left_match_flag = false;
    //Skip if match length is too small
    if (leftStep >= Min_Window_Size) {

        for(ls = left_sp ; ls < left_ep ; ls++) {
            leftIndex = fm_backward.locate(ls);

            int check = LF_mapping(leftIndex, leftIndex + leftStep - 1, left_mate_flank_size, 0, left_mate_flank, right_mate_flank);
            if (check != -1) {
                left_read_flank.clear();
                //for (li = 0 ; li < left_read_flank_size ; li++ ) {
                for (li = left_read_flank_size - 1 ; li > -1 ; li--) {
                    //cout << (long unsigned int)right_read_ori.at(ri) << "-";
                    left_read_flank.push_back(read.at(li));
                }
                
                matchPos = leftIndex - left_mate_flank_size;
                left_matchDist = Edit_dist(left_read_flank, left_mate_flank);
                if (left_matchDist <= Dist_Threshold) {
                    matchIndexes.push_back(matchPos);
                    matchDistances.push_back(left_matchDist);
                    left_match_flag = true;
		            left_count++;
                }
            }
        }
    }

    right_mate_flank_size = read.size() - rightStep;
    right_read_flank_size = read.size() - rightStep;
    bool right_match_flag = false;
    //Skip if match length is too small
    if (rightStep >= Min_Window_Size) {

        for(size_t rs = right_sp ; rs < right_ep ; rs++) {
            size_t rightIndex = getIndexInOriginal(fm_forward.locate(rs));

            int check = LF_mapping(rightIndex - rightStep + 1, rightIndex, 0, right_mate_flank_size, left_mate_flank, right_mate_flank);
            if (check != -1) {
                right_read_flank.clear();
                for (ri = rightStep ; ri < rightStep + right_read_flank_size ; ri++ ) {
                    right_read_flank.push_back(read.at(ri));
                }

                matchPos = rightIndex + 1 - rightStep;
                right_matchDist = Edit_dist(right_read_flank, right_mate_flank);
                if (right_matchDist <= Dist_Threshold) {
                    matchIndexes.push_back(matchPos);
                    matchDistances.push_back(right_matchDist);
                    right_match_flag = true;
					right_count++;
                }
            }
        }
    }

    if (left_match_flag||right_match_flag) {
        return 1;
    }
    else {
        return 0;
    }

}

/*
 * Random search strategy: search backward and forward from a random position
 *
 */
int ReadAlignment::getMatchRandom(const vector<uint8_t> &left_read, const vector<uint8_t> &right_read,
                                   size_t left_sp, size_t left_ep, size_t right_sp, size_t right_ep, uint16_t leftStep, uint16_t rightStep,
                                   vector<size_t> &matchIndexes, vector<uint16_t> &matchDistances) {
    int match_count = 0;
    int left_count = 0;
    int right_count = 0;

    int16_t li, ri;
    size_t ls, rs;
    size_t leftIndex, rightIndex;

    size_t matchPos; // position for mate (first index on reference)
    uint16_t leftDist, rightDist, matchDist; // distance for mate
    vector<uint8_t> left_read_flank, right_read_flank; // part of read from the left-most/right-most index of read
    vector<uint8_t> left_mate_flank, right_mate_flank; // part of mate from the left-most/right-most index of perfect mate
    uint16_t left_read_flank_size, right_read_flank_size; //
    uint16_t left_mate_flank_size, right_mate_flank_size; //

    //Computing distance between remain part of query and all extended part of matches
    left_read_flank_size = left_read.size() - leftStep;
    left_mate_flank_size = left_read.size() - leftStep;
    right_read_flank_size = right_read.size() - rightStep;    
    right_mate_flank_size = right_read.size() - rightStep;    
    bool match_flag = false;
    for(ls = left_sp ; ls < left_ep ; ls++) {
        leftIndex = fm_backward.locate(ls);

        for(rs = right_sp ; rs < right_ep ; rs++) {
            rightIndex = getIndexInOriginal(fm_forward.locate(rs));

            if ((rightIndex + 1 - rightStep) - (leftIndex + leftStep) == 0) { // left and right mate make a contiguous block
/*
                //match_count++;
                //Test: skip computing Edit distance
                matchPos = leftIndex - left_mate_flank_size;
                matchDist = 0;
                matchIndexes.push_back(matchPos);
                matchDistances.push_back(matchDist);
                match_flag = true;
*/
                int check = LF_mapping(leftIndex, rightIndex, left_mate_flank_size, right_mate_flank_size, left_mate_flank, right_mate_flank);
                if (check != -1) {
                    left_read_flank.clear();
                    for (li = left_read_flank_size - 1 ; li > -1 ; li--) {
                        //cout << (long unsigned int)left_read.at(li) << "-";
                        left_read_flank.push_back(left_read.at(li));
                    }
                    leftDist = Edit_dist(left_read_flank, left_mate_flank);

                    right_read_flank.clear();
                    for (ri = rightStep ; ri < rightStep + right_read_flank_size ; ri++ ) {
                        //cout << (long unsigned int)right_read.at(ri) << "-";
                        right_read_flank.push_back(right_read.at(ri));
                    }
                    rightDist = Edit_dist(right_read_flank, right_mate_flank);
                    
                    matchPos = leftIndex - left_mate_flank_size;
                    matchDist = leftDist + rightDist;
                    if (matchDist <= Dist_Threshold) {
                        matchIndexes.push_back(matchPos);
                        matchDistances.push_back(matchDist);
                        match_flag = true;
						match_count++;
                    }
                }
            }
        }
    }

    if (match_flag) {
        return 1;
    }

    left_read_flank_size = left_read.size() - leftStep;
    left_mate_flank_size = left_read.size() - leftStep;
    right_read_flank_size = right_read.size();
    right_mate_flank_size = right_read.size();
    bool left_match_flag = false;
    //Skip if match length is too small
    if (leftStep >= Min_Window_Size) {

        // Take mate with left alignment
        for(ls = left_sp ; ls < left_ep ; ls++) {
            leftIndex = fm_backward.locate(ls);

            //match_count++;

            int checkL = LF_mapping(leftIndex, leftIndex + leftStep - 1,
                                    left_mate_flank_size, right_mate_flank_size, left_mate_flank, right_mate_flank);
            if (checkL != -1) {
                left_read_flank.clear();
                for (li = left_read_flank_size - 1 ; li > -1 ; li--) {
                    //cout << (long unsigned int)left_read.at(li) << "-";
                    left_read_flank.push_back(left_read.at(li));
                }
                leftDist = Edit_dist(left_read_flank, left_mate_flank);
                
                right_read_flank.clear();
                for (ri = 0 ; ri < right_read_flank_size ; ri++ ) {
                    right_read_flank.push_back(right_read.at(ri));
                }
                rightDist = Edit_dist(right_read_flank, right_mate_flank);

                matchPos = leftIndex - left_mate_flank_size;
                matchDist = leftDist + rightDist;
                if (matchDist <= Dist_Threshold) {
                        matchIndexes.push_back(matchPos);
                        matchDistances.push_back(matchDist);
                        left_match_flag = true;
						left_count++;
                }
            }
        }
    }

    left_read_flank_size = left_read.size();
    left_mate_flank_size = left_read.size();
    right_read_flank_size = right_read.size() - rightStep;
    right_mate_flank_size = right_read.size() - rightStep;
    bool right_match_flag = false;
    //Skip if match length is too small
    if (rightStep >= Min_Window_Size) {
        // Take mate with right alignment
        for(rs = right_sp ; rs < right_ep ; rs++) {
            rightIndex = getIndexInOriginal(fm_forward.locate(rs));
            
            //match_count++;
    	    int checkR = LF_mapping(rightIndex - rightStep + 1, rightIndex,
                                    left_mate_flank_size, right_mate_flank_size, left_mate_flank, right_mate_flank);
    	    if (checkR != -1) {
                left_read_flank.clear();
    			for (li = left_read_flank_size - 1 ; li > -1 ; li--) {
    			    //cout << (long unsigned int)left_read.at(li) << "-";
    			    left_read_flank.push_back(left_read.at(li));
    			}
    			leftDist = Edit_dist(left_read_flank, left_mate_flank);

                right_read_flank.clear();
    			//for (ri = 0 ; ri < right_mate_flank_size; ri++ ) {
    			for (ri = rightStep ; ri < rightStep + right_read_flank_size ; ri++ ) {
                    //cout << (long unsigned int)right_read_ori.at(ri) << "-";
    			    right_read_flank.push_back(right_read.at(ri));
    			}
    			rightDist = Edit_dist(right_read_flank, right_mate_flank);

    			matchPos = rightIndex + 1 - rightStep - left_mate_flank_size;
                matchDist = leftDist + rightDist;
    			if (matchDist <= Dist_Threshold) {
    			    matchIndexes.push_back(matchPos);
    			    matchDistances.push_back(matchDist);
                    right_match_flag = true;
					right_count++;
    			}
            }
        }
    }

    if (left_match_flag||right_match_flag) {
        return 1;
    }
    else {
        return 0;
    }

}


/*
 * Take the mate of the query on the reference
 *
 */
inline int ReadAlignment::LF_mapping(size_t leftIndex, size_t rightIndex, uint16_t left_mate_flank_size, uint16_t right_mate_flank_size,
                                 vector<uint8_t> &left_mate_flank, vector<uint8_t> &right_mate_flank) {

    //Make sure we don't go to the index < 0
    int64_t left_check = leftIndex - left_mate_flank_size;
    if (left_check < 0) {
        return -1;
	}
    left_mate_flank.clear();
    for (size_t li = leftIndex - 1 ; li > leftIndex - 1 - left_mate_flank_size ; li--) {
        left_mate_flank.push_back(Ref.at(li));
    }
    //Make sure we don't go to the index > Ref_Size - 1
    uint64_t right_check = rightIndex + right_mate_flank_size;
    if (right_check > Ref_Size - 1) {
        return -1;
	}
    right_mate_flank.clear();
    for (size_t ri = rightIndex + 1 ; ri < rightIndex + 1 + right_mate_flank_size ; ri++) {
        right_mate_flank.push_back(Ref.at(ri));
    }
    return 0;
}