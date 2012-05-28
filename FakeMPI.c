#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
/*for socket*/
#include <sys/types.h>
#include <sys/socket.h>
/*gethostbyname*/
#include <netdb.h>

#include "FakeMPI.h"
/*
 *MPI/f90の相性が非常に悪く、cとのリンクもままならないので、mpiの一部機能を自作する。
 *平成17年2月9日(水)
 *ノード間の通信はややこしいので、すべて親ノードに集約する。
 */

/*
 *エラーメッセージ出力のための大域変数。stderrにすれば表示される。
 */
FILE* fakempi_verbose = NULL;

int accept_child( int sock )
{
    while( 1 ){
        int msgsock = accept(sock, (struct sockaddr *)0, (int *)0);
        if (msgsock == -1) {
            //perror("accept");
        }
        else {
            fprintf( fakempi_verbose, "Accepted\n" );
            return msgsock;
        }
    }
}




int get_connected( int port )
{
    struct hostent *hp;
    struct sockaddr_in server;
    int sock = socket(AF_INET, SOCK_STREAM, 0);
    if (sock < 0) { 
	perror("opening stream socket");
	exit(1);
    }

    hp = gethostbyname( "localhost" );
    if ( hp == NULL ) {
	fprintf(fakempi_verbose, "unknown host\n" );
	exit(2);
    }

    server.sin_family = AF_INET;
    server.sin_port = htons( port );
    //memset(&server.sin_addr,0,sizeof(server.sin_addr));
    bcopy((char *)hp->h_addr, 
          &server.sin_addr,
          hp->h_length);

    while  (connect(sock, (struct sockaddr *)&server, 
		sizeof server) < 0) {
	perror("connecting stream socket");
        sleep(1);
    }
    return sock;
}



int lookup_unused_port( int port )
{
    int sock;
    int i;
    struct sockaddr_in client;

    if((sock=socket(AF_INET,SOCK_STREAM,0))<0){
        perror("server: socket");
        exit(1);
    }
    i=1;
    setsockopt(sock,SOL_SOCKET,SO_REUSEADDR,&i,sizeof(i));
    client.sin_family=AF_INET;
    while( 1 ){
        client.sin_port=htons( port );
        memset(&client.sin_addr,0,sizeof(client.sin_addr));
        if(0==bind(sock,(struct sockaddr *)&client,sizeof(client))){
            break;
        }
        port += 10;
    }
    close( sock );
    return port;
}



int listen_socket( int port )
{
    int sock;
    int i;
    struct sockaddr_in client;

    if((sock=socket(AF_INET,SOCK_STREAM,0))<0){
        perror("server: socket");
        exit(1);
    }
    i=1;
    setsockopt(sock,SOL_SOCKET,SO_REUSEADDR,&i,sizeof(i));
    client.sin_family=AF_INET;
    client.sin_port=htons( port );
    memset(&client.sin_addr,0,sizeof(client.sin_addr));
    if(bind(sock,(struct sockaddr *)&client,sizeof(client))){
        perror("server: bind");
        exit(1);
    }
    if(listen(sock,5)){
        perror("server: listen");
        exit(1);
    }
    return sock;
}


void fakempi_connect( sFakeMPI* fakempi )
{
    /*
     *親ノードなら
     */
    if ( fakempi->myrank == 0 ){
        /*
         *socketをnprocs-1だけ開く。ポートが重複しないようにするには、あらかじめそのポートが使えるかどうかを調べておいた方がよかろう。
         */
        int i;
        for( i=1; i<fakempi->nprocs; i++ ){
            fakempi->socks[i] = listen_socket( i + fakempi->portbase );
        }
        fprintf( fakempi_verbose, "Parent is listening.\n" );
        for( i=1; i<fakempi->nprocs; i++ ){
            fakempi->msgsocks[i] = accept_child( fakempi->socks[i] );
        }
    }
    else{
        /*
         *親のソケットに接続する。
         */
        sleep( 0 );
        fakempi->mysock = get_connected( fakempi->myrank + fakempi->portbase );
        fprintf( fakempi_verbose, "Child %d is get connected.( socket %d ) \n", fakempi->myrank, fakempi->mysock );
    }
}



/*
 *オブジェクトの初期化とプロセスの分岐
 */
sFakeMPI* fakempi_new( int nprocs )
{
    sFakeMPI* fakempi;
    int i, pid;
    fakempi = malloc( sizeof(sFakeMPI) );
    fakempi->nprocs = nprocs;
    fakempi->myrank = 0;
    /*
     *redirect it to fakempi_verbose if you want to get verbose messages.
     */
    if ( fakempi_verbose == NULL ){
        fakempi_verbose = fopen( "/dev/null", "w" );
    }
    fakempi->portbase = lookup_unused_port( 9800 ) - 1;
    fprintf( fakempi_verbose, "portbase: %d\n", fakempi->portbase );
    fakempi->pids   = malloc( sizeof(int) * nprocs );
    fakempi->socks  = malloc( sizeof(int) * nprocs );
    fakempi->msgsocks = malloc( sizeof(int) * nprocs );
    for( i=1; i<nprocs; i++ ){
        pid = fork();

        /*
         *まず、親プロセスで、子プロセスのPIDのリストを作る。
         *できればあとで子プロセスに渡す。
         *子プロセスは、自分が何番目の子であるかをこの時に知る。
         */
        if ( pid == 0 ){
            /*
             * i'm a child.
             */
            fakempi->myrank = i;
            break;
        }
        else{
            /*
             *i'm a parent.
             */
            fakempi->pids[i] = pid;
        }
    }
    fakempi->mypid = getpid();
    fprintf( fakempi_verbose, "%d/%d: %d\n", fakempi->myrank, fakempi->nprocs, fakempi->mypid );
    if ( fakempi->myrank == 0 ){
        fakempi->pids[0] = fakempi->mypid;
    }
    fakempi_connect( fakempi );
    return fakempi;
}




int receive_data( int msgsock, char* buf, int length )
{
    int rval;
    int count=0;
    bzero(buf, length);
    do {
        if ((rval = read(msgsock, buf, length)) < 0) {
            perror("reading stream message");
        }
        count += rval;
        buf += rval;
    } while ( count < length );
    return count;
    
}


void send_data( int sock, char* msg, int length )
{
    if ( write(sock, msg, length) < 0) {
	perror("writing on stream socket");
    }
}


       

/*
 *バイト列を親ノードに集約する。
 *endianのことが気になる。
 *同時通信は難しいが、できるかな。
 */
int fakempi_gather( sFakeMPI* fakempi, int length, void* msg, char* array )
{
    if ( fakempi->myrank == 0 ){
        int i;
        /*
         *親のデータはそのまま配列に。ただし、msgとarray[0]が同じアドレスを差している可能性があることに注意
         */
        memmove( (void*)array, msg, length );
        /*
         *親ノードは子の問いかけを待つ。同時に待つことはできるのか？
         */
        for( i=1; i<fakempi->nprocs; i++ ){
            int size;
            size = receive_data( fakempi->msgsocks[i], array + i*length, length );
            fprintf( fakempi_verbose, "%d>%d bytes received. \n", fakempi->myrank, size );
        }
    }
    else{
        /*
         *子ノードは自分のデータを送るのみ
         */
        send_data( fakempi->mysock, msg, length );
        fprintf( fakempi_verbose, "%d>Sent. \n", fakempi->myrank );
    }
}



/*
 *親ノードの配列の各要素を子にばらまく。あえて引数順はそのままに。
 */
int fakempi_scatter( sFakeMPI* fakempi, int length, void* msg, char* array )
{
    if ( fakempi->myrank == 0 ){
        int i;
        memmove( msg, (void*)array, length );
        for( i=1; i<fakempi->nprocs; i++ ){
            send_data( fakempi->msgsocks[i], array + i*length, length );
            fprintf( fakempi_verbose, "%d>Sent. \n", fakempi->myrank );
        }
    }
    else{
        int size;
        size = receive_data( fakempi->mysock, msg, length );
        fprintf( fakempi_verbose, "%d>%d bytes received. \n", fakempi->myrank, size );
    }
}




/*
 *親ノードのデータをすべての子にコピーする。
 */
int fakempi_bcast( sFakeMPI* fakempi, int length, void* msg )
{
    if ( fakempi->myrank == 0 ){
        int i;
        fprintf( fakempi_verbose, "Parent starts sending...\n" );
        for( i=1; i<fakempi->nprocs; i++ ){
            fprintf( fakempi_verbose, "Parent sends. %d\n", i );
            send_data( fakempi->msgsocks[i], msg, length );
            fprintf( fakempi_verbose, "%d>Sent. \n", fakempi->myrank );
        }
    }
    else{
        int size;
        size = receive_data( fakempi->mysock, msg, length );
        fprintf( fakempi_verbose, "%d>%d bytes received. \n", fakempi->myrank, size );
    }
}



/*
 *wait children die.
 */
void fakempi_finalize( sFakeMPI* fakempi )
{
    /*
     *wait children die.
     */
    if ( fakempi->myrank == 0 ){
        int i;
        for( i=1; i<fakempi->nprocs; i++ ){
            wait( fakempi->pids[i] );
        }
    }
}





//#define TESTCASE
#undef TESTCASE
#ifdef TESTCASE
/*
 *test case
 */
#define N 10
int main( int argc, char* argv[] )
{
    sFakeMPI* fakempi;
    int* childpid;
    int i, retrieve;
    childpid = (int*)malloc( sizeof( int ) * N );
    fakempi = fakempi_new( N );
    fprintf( fakempi_verbose, "fakempi_new done %d/%d\n", fakempi->myrank, fakempi->nprocs );
    fakempi_gather( fakempi, sizeof( int ), &fakempi->mypid, (char*)childpid );
    fprintf( fakempi_verbose, "fakempi_gather done %d/%d\n", fakempi->myrank, fakempi->nprocs );
    fakempi_scatter( fakempi, sizeof( int ), &retrieve, (char*)childpid );
    fprintf( fakempi_verbose, "fakempi_scatter done %d/%d\n", fakempi->myrank, fakempi->nprocs );
    for( i=0; i<N; i++ ){
        printf( "%d>PID %d: %d -- %d\n", fakempi->myrank, i, childpid[i], fakempi->pids[i] );
    }
    fakempi_bcast( fakempi, sizeof( int ) * N, (char*)childpid );
    fprintf( fakempi_verbose, "fakempi_bcast done %d/%d\n", fakempi->myrank, fakempi->nprocs );
    for( i=0; i<N; i++ ){
        printf( "%d>PID %d: %d -- %d\n", fakempi->myrank, i, childpid[i], fakempi->pids[i] );
    }
    sleep( 10 );
    fakempi_finalize( fakempi );
    exit(0);
}

            
#endif
