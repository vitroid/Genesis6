/*openMosix-related functions*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* get openMosix node id from /proc filesystem */
int getnodeid_()
{
    FILE* file;
    char  filename[100];
    int   node;

    sprintf( filename, "/proc/%d/where", getpid() );
    //fputs( filename, stderr );
    if ( NULL != ( file = fopen( filename, "r" ) ) ){
        fscanf( file, "%d", &node );
        fclose( file );
        return node;
    }
    return -1;
}
