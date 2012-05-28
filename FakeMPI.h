#ifndef FAKEMPI_H
#define FAKEMPI_H
#include <stdio.h>

typedef struct {
    int nprocs;
    int myrank;
    int mypid;
    int portbase;
    /*for the parent.*/
    int* pids;
    int* socks;
    int* msgsocks;
    /*for children.*/
    int mysock;
} sFakeMPI;

/*
 *エラーメッセージ出力のための大域変数。stderrにすれば表示される。
 */
extern FILE* fakempi_verbose;
/*
 *オブジェクトの初期化とプロセスの分岐
 */
sFakeMPI* fakempi_new( int nprocs );
/*
 *バイト列を親ノードに集約する。
 */
int fakempi_gather( sFakeMPI* fakempi, int length, void* msg, char* array );
/*
 *親ノードの配列の各要素を子にばらまく。
 */
int fakempi_scatter( sFakeMPI* fakempi, int length, void* msg, char* array );
/*
 *親ノードのデータをすべての子にコピーする。
 */
int fakempi_bcast( sFakeMPI* fakempi, int length, void* msg );
/*
 *子プロセスを殺す。
 */
void fakempi_finalize( sFakeMPI* fakempi );

#endif
