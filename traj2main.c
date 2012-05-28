#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>

void usage( int argc, char * argv[] )
{
    printf( "usage: %s [-wn]\n", argv[0] );
    printf( "\t-n|--ngph\tOutput HBNetwork topology in @NGPH format.\n" );
    printf( "\t-w|--wgph\tOutput HBNetwork topology in @WGPH format.\n" );
    printf( "\t-m|--mdvw\tOutput molecular cofigurations in @MDVW format.\n" );
    printf( "\t-r|--rdf\tOutput Radial Distribution Function.\n" );
    printf( "\t-p|--pair\tOutput pain interactions.\n" );
    exit(1);
}

int main( int argc, char *argv[] )
{
    int c;
    int digid_optind = 0;
    int outputtype = 0; /* 1: WGPH */
    int stream = 6;
    
    while (1){
        int this_option_optind = optind ? optind : 1;
        int option_index = 0;
        static struct option long_options[] = {
            {"wgph", 0, 0, 'w'},
            {"ngph", 0, 0, 'n'},
            {"mdvw", 0, 0, 'm'},
            {"rdf", 0, 0, 'r'},
            {"pair", 0, 0, 'p'},
            {"near", 0, 0, 'l'},
            {"threebody", 0, 0, '3'},
            {0, 0, 0, 0}
        };

        c = getopt_long( argc, argv, "nwmrpl3", long_options, &option_index );

        if ( c == -1 )
            break;
        
        switch (c) {

        case 'w':
            outputtype = 1; /* WGPH */
            break;

        case 'n':
            outputtype = 2; /* NGPH */
            break;

        case 'm':
            outputtype = 3; /* MDVW */
            break;

        case 'r':
            outputtype = 4; /* RDF  */
            break;

        case 'p':
            outputtype = 5; /* PAIR */
            break;

        case 'l':
            outputtype = 6; /* NEAR */
            break;
        case '3':
            outputtype = 7; /* NEAR */
            break;

        }
    }
    
    if ( optind < argc ){
        usage ( argc, argv );
    }

    if ( outputtype == 0 )
        usage ( argc, argv );
    
    /* call f90 */
    traj2_( &stream, &outputtype );
    exit(0);
}

