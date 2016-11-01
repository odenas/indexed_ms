//
//  main.c
//  djamal_bwt
//
//  Created by denas on 10/30/16.
//  Copyright © 2016 denas. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dbwt.h"

int main(int argc, const char * argv[]) {
    unsigned int last;
    char *s = "babcabcbbabababcbcbbcbcbabcbacabcbbcbcbabbcba";

    unsigned char *bwt = dbwt_bwt((unsigned char*) s, (uint) strlen(s), &last, 0);
    if(bwt){
        for(int i=0; i<strlen(s) + 1; i++)
            printf("BWT[%d] = %c\n", i, (i == last ? '#' : bwt[i]));
    }
    return 0;
}
