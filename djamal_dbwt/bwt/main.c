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


int find_space (char *line){
    int i = 0;
    while(line[i] != ' ')
        i++;
    return i;
}


int main(int argc, const char * argv[]) {
    size_t line_length = 0, t_length = 0, s_length = 0;
    char *line, *t, *s;
    unsigned int last;

    while(getline(&line, &line_length, stdin) > 0){
        line_length = strlen(line);
        line[line_length - 1] = '\0'; // remove '\n'
        t_length = find_space(line);
        s = &line[t_length + 1];
        s_length = strlen(s);
        line[t_length] = '\0'; t = line;

        unsigned char *bwt = dbwt_bwt((unsigned char*) s, s_length, &last, 0);
        printf("%s %s ", t, s);
        if(bwt){
            for(int i=0; i<strlen(s) + 1; i++)
                printf("%c", (i == last ? '#' : bwt[i]));
            printf("\n");
        }
    }
    return 0;
}
