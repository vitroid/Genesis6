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
 *���顼��å��������ϤΤ��������ѿ���stderr�ˤ����ɽ������롣
 */
extern FILE* fakempi_verbose;
/*
 *���֥������Ȥν�����ȥץ�����ʬ��
 */
sFakeMPI* fakempi_new( int nprocs );
/*
 *�Х������ƥΡ��ɤ˽��󤹤롣
 */
int fakempi_gather( sFakeMPI* fakempi, int length, void* msg, char* array );
/*
 *�ƥΡ��ɤ�����γ����Ǥ�ҤˤФ�ޤ���
 */
int fakempi_scatter( sFakeMPI* fakempi, int length, void* msg, char* array );
/*
 *�ƥΡ��ɤΥǡ����򤹤٤ƤλҤ˥��ԡ����롣
 */
int fakempi_bcast( sFakeMPI* fakempi, int length, void* msg );
/*
 *�ҥץ����򻦤���
 */
void fakempi_finalize( sFakeMPI* fakempi );

#endif
