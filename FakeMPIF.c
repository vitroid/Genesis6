#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "FakeMPI.h"

static sFakeMPI* fakempi;

/*
 *fortran interfaces
 */

void fakempif_new_( int *nprocs, int* ierr )
{
    fakempi = fakempi_new( *nprocs );
    *ierr = ( fakempi == NULL );
    //fakempi_verbose = stderr;
    resettimer();
}

void fakempif_finalize_( int* ierr )
{
    fakempi_finalize( fakempi );
    fakempi = NULL;
    *ierr = 0;
}

void fakempif_gather_( int* length, void* msg, char* array, int* ierr )
{
    *ierr = fakempi_gather( fakempi, *length, msg, array );
}

void fakempif_scatter_( int* length, void* msg, char* array, int* ierr )
{
    *ierr = fakempi_scatter( fakempi, *length, msg, array );
}

void fakempif_bcast_( int* length, void* msg, int* ierr )
{
    *ierr = fakempi_bcast( fakempi, *length, msg );
}

void fakempif_comm_size_( int* size, int* ierr )
{
    *ierr = 0;
    *size = fakempi->nprocs;
}

void fakempif_comm_rank_( int* rank, int* ierr )
{
    *ierr = 0;
    *rank = fakempi->myrank;
}

void fakempif_timer_delta_( int* delta )
{
    *delta = deltatime();
}
