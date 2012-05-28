module graph2
  !
  !@NGPH類似のデータ型式
  !
  use common_module
  implicit none
  integer, private, parameter :: MAXNEI=100
  type pGraph
     integer :: nnode
     integer :: nnei( MAXMOL )
     integer :: nei( MAXNEI, MAXMOL )
  end type pGraph

  !
  !同じものを別の表現にしたもの。
  !
  type pPair
     integer :: nbond
     integer, pointer :: a(:)
     integer, pointer :: b(:)
  end type pPair

  interface new
     module procedure graph_new
     module procedure pair_new
  end interface

  interface done
     module procedure graph_done
  end interface

contains

  subroutine graph2pair( graph, pair )
    type(pGraph), pointer :: graph
    type(pPair),  pointer :: pair

    integer :: nbond, i,j
    nbond = 0
    do i=1, graph%nnode
       nbond = nbond + graph%nnei(i)
    enddo
    call new( pair, nbond )
    nbond = 0
    do i=1, graph%nnode
       do j=1, graph%nnei(i)
          nbond = nbond + 1
          pair%a(nbond) = i
          pair%b(nbond) = graph%nei(j,i)
       enddo
    enddo
  end subroutine graph2pair

  subroutine graph_new( graph, size )
    type(pGraph), pointer :: graph
    integer, optional, intent(IN) :: size
    integer :: s
    s = 0
    allocate( graph )
    if ( present( size ) ) then
       s = size
    endif
    graph%nnode = s
    graph%nnei(1:s) = 0
  end subroutine graph_new

  subroutine pair_new( pair, size )
    type(pPair), pointer :: pair
    integer, optional, intent(IN) :: size
    integer :: s
    s = 0
    allocate( pair )
    if ( present( size ) ) then
       s = size
       allocate( pair%a(s) )
       allocate( pair%b(s) )
    endif
    pair%nbond = s
  end subroutine pair_new

  subroutine graph_done( graph )
    type(pGraph), pointer :: graph
    deallocate( graph )
  end subroutine graph_done

  subroutine graph_link( graph, i,j )
    type(pGraph), pointer :: graph
    integer, intent(IN)   :: i,j
    graph%nnei(i) = graph%nnei(i) + 1
    graph%nei( graph%nnei(i), i) = j
  end subroutine graph_link

  subroutine graph_perfect( graph, size )
    use error_module
    type(pGraph), pointer :: graph
    integer, intent(IN) :: size

    integer :: i,j
    if ( MAXNEI < size ) then
       call die( 0, "TOO LARGE GRAPH." )
    endif
    graph%nnode = size
    graph%nnei(:) = 0
    do i = 1, size
       do j = i+1, size
          call graph_link( graph, i, j )
          call graph_link( graph, j, i )
       enddo
    enddo
  end subroutine graph_perfect

  subroutine graph_read( graph, file )
    type(pGraph), pointer :: graph
    integer, intent(IN)   :: file

    integer :: i,j
    read( file,* ) graph%nnode
    graph%nnei(:) = 0
    do
       read( file, * ) i,j
       if ( i < 0 ) return
       i = i + 1
       j = j + 1
       call graph_link( graph, i, j )
       call graph_link( graph, j, i )
    enddo
  end subroutine graph_read

  !
  !グラフから連結な対を抽出する。各対が同一のノードを共有しないように選ぶ。
  !最大マッチであることを要請すると、組み合わせが固定的になってしまう可能性があるので、
  !もっと緩く探索する。
  subroutine graph_matching( graph, rand, subgraph )
    use random_module
    type(pRandom), pointer :: rand
    type(pGraph), pointer :: graph
    type(pGraph), pointer :: subgraph ! must not be allocated

    integer :: ii,i,j,k,n
    integer :: nnei, nei( graph%nnode ), node( graph%nnode )
    integer :: mark( graph%nnode )

    n = graph%nnode
    call new( subgraph, n )
    mark(:) = 0
    do i=1, n
       node(i) = i
    enddo
    call shuffle( rand, n, node )
    write(6,*) (node(i),i=1,n)
    do ii=1, n
       i = node(ii)
       if ( mark(i) == 0 ) then
          nnei = graph%nnei( i )
          if ( 0 < nnei ) then
             nei(1:nnei) = graph%nei( 1:nnei, i )
             call shuffle( rand, nnei, nei )
             inner:do j=1, nnei
                k = nei( j )
                if ( mark(k) == 0 ) then
                   mark(k) = 1
                   mark(i) = 1
                   call graph_link( subgraph, i, k )
                   call graph_link( subgraph, k, i )
                   exit inner
                endif
             enddo inner
          endif
       endif
    enddo
  end subroutine graph_matching

  subroutine Graph_show( graph, file )
    type(pGraph), pointer :: graph
    integer, intent(IN) :: file

    integer :: i,j
    do i=1, graph%nnode
       do j=1, graph%nnei(i)
          write(file,*) i,"-->",graph%nei(j,i)
       enddo
    enddo
  end subroutine Graph_show


  subroutine Pair_show( pair, file )
    type(pPair), pointer :: pair
    integer, intent(IN) :: file
    !local
    integer :: k
    do k=1, pair%nbond
       write(file,*) pair%a(k),"-->",pair%b(k)
    enddo
  end subroutine Pair_show
end module graph2

