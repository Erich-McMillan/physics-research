//
//  main.c
//  Test
//
//  Created by erich on 5/19/15.
//  Copyright (c) 2015 erich. All rights reserved.
//

#include <stdio.h>

void test( int rand[] )
{
    for( int i = 1; i < 10; i++ )
    {
        rand[i] = i;
    }
}

int main(int argc, const char * argv[])
{
    int rand[10];
    test( rand );
    for( int i = 1; i < 10; i++ )
    {
        int f = rand[i];
        printf( "%d\n", f );
    }
}

