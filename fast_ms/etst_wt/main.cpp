//
//  test_wt.cpp
//  fast_ms
//
//  Created by denas on 11/4/16.
//  Copyright © 2016 denas. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/bit_vectors.hpp>

extern "C" {
#include "dbwt.h"
}

using namespace std;
using namespace sdsl;


int main(int argc, char **argv){
    string s {"babbabbbabaababbabaabaaaaabbbaabbbaaaaaaaaaabababbbaabbbbbabbbbabaabaaabbabaaabbbaabbaabbabbbababaaaabbbbbaababbaaabbbbabbbbaabaabbbaaabbbaabaaaaabbabbbaabaaaababbbabbbabaaabbaaabbaaabbbaabaaabbbbabbbbbbbbaabbababaabaabbbaaababbaabaaaaaabbbabbabaaababaabaabbaabaaaabbbabaabbaababbabbbbaabbabbabbaaaaabababbbababaabaabbbaabaaabbabbabbabbaaabbbbbbbaaabbabaaaabbbabababbbbbbbaabbaabbbaaaaaaaaaabaabbabaababbbabbaaaaabbbbaaabbbbbbbbaabbbbaabbbaaaabbbbabaaaaaabbbbbaababbbbaaabaaabaaababbbaaaaabaaababbaaaabbaabbabababababaabaaabaaaabbababbabaaabaababbababbaabbaaaaabababbbbaababbbaabababbbbabbaaaaababaaa"};
    unsigned int last = 0;
    unsigned char text[s.size() + 1];

    for(int i=0; i<s.size(); i++)
        text[i] = s[i];
    text[s.size()] = '\0';

    unsigned char *bwt = dbwt_bwt((unsigned char *)text, s.size(), &last, 0);
    bwt[last] = '#';
    for(int i=0; i<s.size(); i++)
        cout << bwt[i];
    cout << endl;
    return 0;
}

