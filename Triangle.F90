! -*- f90 -*-
!三体角(O-O-O)でMCを行うための算程群
!
module triangle
  use set_module
  use neighbors_module
  use distrib_module
  use umbrella
  implicit none

  type sTriangle
#ifdef DEBUG
     type(sHistogram), pointer :: lasthist
#ifdef DEBUG2
     type(sHistogram), pointer :: tmphist
     type(sNeighbors) :: nei2copy(MAXMOL)
#endif
#endif
     type(sHistogram), pointer :: hist, newhist
     type(sSet) :: set
     type(sNeighbors) :: nei(MAXMOL)  ! こちらは最近接4分子だけ。
     type(sNeighbors) :: nei2(MAXMOL) ! こちらは5A以内の近傍にある分子。
     !
     !temporary neighbor list
     !
     type(sNeighbors) :: newnei(99)   !newnei2と同じ大きさが必要。
     integer          :: newnei2idx(99)
     type(sNeighbors) :: newnei2(99)
     integer          :: nmodified
     logical          :: logarithmic
  end type sTriangle


  type sTrianglePtr
     type(sTriangle), pointer :: p
  end type sTrianglePtr

  interface new
     module procedure sTriangle_new
  end interface

contains

  subroutine sTriangle_new( tri, islog )
    type( sTriangle ), pointer :: tri
    logical, intent(IN) :: islog
    allocate( tri )
    tri%logarithmic = islog
    nullify( tri%hist )
    nullify( tri%newhist )
  end subroutine sTriangle_new

  subroutine sTriangle_prepareall( tri, nmol, ww, com, refhist, op )
    use interaction_module
    use vector_module
    type( sTriangle ) :: tri
    integer, intent(IN) :: nmol
    type( sInteraction ), intent(IN) :: ww
    type(vector3), intent(IN) :: com( * )
    type( sHistogram ), pointer :: refhist
    real(kind=8), intent(OUT) :: op

    integer :: i,j,k
    real(kind=8) :: costh
    !
    !ある分子の近傍4分子の成す角度の分布を求める。
    !1成分のみ。
    !
    call new( tri%hist, 21, 0.1d0, -1d0 )
    !
    !5A以内にいる隣接分子の表。1 MC stepの間維持する。
    !近傍分子があまり遠くにいってしまわないことが前提。
    !
#ifdef DEBUG2
    if ( loop .ne. 1 ) then
       do i=1, nmol
          tri%nei2copy(i) = tri%nei2(i)
       enddo
    endif
#endif
    call neighbors_all2( nmol, tri%nei2, ww, com, 0, ORDER_BY_DISTANCE, 5d0 )
    !
    !nei2を元に、各分子に最近接な4分子の表を生成する。必要に応じて生成する。
    !
    do i=1, nmol
       call neighbor1_compress2( tri%nei(i), tri%nei2(i), 4, ORDER_BY_DISTANCE )
    enddo
    do i=1, nmol
       do j=1, tri%nei(i)%nnear
          do k=j+1, tri%nei(i)%nnear
             !
             ! O-O-O角の余弦を計算する。
             !
             costh = tri%nei(i)%dx(j)*tri%nei(i)%dx(k) + tri%nei(i)%dy(j)*tri%nei(i)%dy(k) + tri%nei(i)%dz(j)*tri%nei(i)%dz(k)
             costh = costh / sqrt ( tri%nei(i)%dd(j) * tri%nei(i)%dd(k) )
             call histogram_accumulate( tri%hist, costh, 1d0 )
          enddo
       enddo
    enddo
    op = kullback_div( tri%hist, refhist )
    !call histogram_show( tri%hist, STDOUT )
    !write(6,*) op, log(op)
    if ( tri%logarithmic ) op = log(op)
    !stop
#ifdef DEBUG
    !
    !前回の積算ヒストグラムと、ここで再計算したヒストグラムが一致するかどうか
    !tmphistで完全に一致している(accept時のhistogramの更新は問題ない)のに、
    !ここで不一致が生じるとすれば、nei2が大きく変化しすぎている可能性がある。
    !
    if ( associated( tri%lasthist ) ) then
       diff = kullback_div( tri%hist, tri%lasthist )
       if ( diff /= 0d0 ) then
          write(STDERR,*) diff, "DIV"
          call histogram_show( tri%hist, STDERR )
          call histogram_show( tri%lasthist, STDERR )
          call die( 0, "Histogram differs." )
       endif
       call done( tri%lasthist )
    endif
#endif /*DEBUG*/
  end subroutine sTriangle_prepareall
    

  subroutine sTriangle_Differentiate( tri, target, deltax, deltay, deltaz, refhist, op )
    type( sTriangle ) :: tri
    integer, intent(IN) :: target
    real(kind=8), intent(IN) :: deltax,deltay,deltaz
    type( sHistogram ), pointer :: refhist
    real(kind=8), intent(OUT) :: op
    
    integer :: i, j, k, ii
    logical :: verbose
    real(kind=8) :: costh
    !
    !まず、5A以内の近傍分子リストを更新する。tri%newnei2には、移動した分の情報だけが返される。
    !
    call neighbors_displacement( tri%nei2, target, deltax, deltay, deltaz, tri%nmodified, tri%newnei2, tri%newnei2idx )
    !
    !target周辺の分子も全部表を作りなおす。
    !(実際には無駄な計算も多数あるだろう。nei2の距離閾値を適切に選べれば、計算量が最小化できる。)
    !
    do i=1, tri%nmodified
       call neighbor1_compress2( tri%newnei(i), tri%newnei2(i), 4, ORDER_BY_DISTANCE )
    enddo
    !
    !tri%nmodifiedの値は、targetが動くことによって、変化した近傍リストの要素数。勘違いでなければ、targetの隣接数+1(target自身も含まれるから)に等しくなると思う。
    !
#ifdef DEBUG
    if ( tri%nei2(target)%nnear+1 /= tri%nmodified ) then
       write(STDERR,*) tri%nei2(target)%nnear+1, tri%nmodified, "should be the same."
    endif
#endif
    call new( tri%set, MAXMOL )
    !
    !移動前の近傍群
    !
    call join( tri%set, target )
    do i=1,tri%nei(target)%nnear
       call join( tri%set, tri%nei(target)%near(i) )
    enddo
    do ii=1, tri%nei2(target)%nnear
       !
       !tgtの近傍5Aにある分子iを選ぶ
       !
       i = tri%nei2(target)%near(ii)
       !
       !iの近傍にtgtはあるはず。
       !
       k = 0
       do j=1, tri%nei(i)%nnear
          if ( tri%nei(i)%near(j) == target ) then
             call join( tri%set, i )
          endif
       enddo
    enddo
    !
    !移動後の近傍群
    !
    do i=1,tri%newnei(1)%nnear
       call join( tri%set, tri%newnei(1)%near(i) )
    enddo
    do i=2, tri%nmodified
       k = 0
       do j=1, tri%newnei(i)%nnear
          if ( tri%newnei(i)%near(j) == target ) then
             call join( tri%set, tri%newnei2idx(i) )
          endif
       enddo
    enddo
    !write(STDERR,*) "SET", (set%set(i), i=1, set%nmark )
    !
    !ある分子が移動すると、どのO-O-O角が変化するか
    !場合によっては、移動することで最近接4分子が代わる可能性がある。
    !
    !histogramを複製する
    !
    call histogram_dup( tri%newhist, tri%hist )
#ifdef DEBUG2
    !call histogram_dup( tmphist, hist )
#endif
    !
    !ヒストグラムから一旦差し引く
    !
    !verbose = loop == 1 .and. trial == 2
    !verbose = .true.
    verbose = .false.
    if ( verbose ) call histogram_show( tri%hist, STDERR )
    call partial_accumulate( tri%newhist, tri%set, tri%nei, tri%nei2, -1d0, verbose )
    !
    !まず、中心分子の、最近接4分子の表を更新する。tri%newnei2が正しく更新されているので、それをもとにすればよい。
    !
    !call neighbor1_compress( tri%newnei(1), tri%newnei2(1), 4 )
    if ( verbose ) then
       write(STDERR,*) target, ">", (tri%newnei(1)%near(i), i=1, tri%newnei(1)%nnear)
    endif
    !
    !ヒストグラムに新たに加算する。
    !
    call partial_accumulate_indexed( tri%newhist, tri%set, tri%newnei, tri%newnei2, +1d0, tri%nmodified, tri%newnei2idx, verbose )
    
#ifdef DEBUG
    if ( verbose ) then
       do i=1, nmol
          do j=1, tri%nei(i)%nnear
             do k=j+1, tri%nei(i)%nnear
                !
                ! O-O-O角の余弦を計算する。
                !
                costh = tri%nei(i)%dx(j)*tri%nei(i)%dx(k) + tri%nei(i)%dy(j)*tri%nei(i)%dy(k) + tri%nei(i)%dz(j)*tri%nei(i)%dz(k)
                costh = costh / sqrt ( tri%nei(i)%dd(j) * tri%nei(i)%dd(k) )
                !if ( 0.55d0 < costh .and. costh < 0.65d0 ) then
                !   write(STDERR,*) "Large:", nei(i)%near(j), i, nei(i)%near(k), nint( ( costh - hist%minvalue ) / hist%binwidth ) + 1
                !endif
                !if ( i .eq. target .or. nei(i)%near(j) .eq. target .or. nei(i)%near(k) .eq.target ) then
                !   write(STDERR,*) "Before:", nei(i)%near(j), i, nei(i)%near(k), nint( ( costh - hist%minvalue ) / hist%binwidth ) + 1
                !endif
             enddo
          enddo
       enddo
    endif
#endif
    
    !
    !ここまで読んだ限りではバグらしいものはない。
    !
    op = kullback_div( tri%newhist, refhist )
    if ( tri%logarithmic ) op = log(op)
  end subroutine sTriangle_Differentiate
  
  subroutine sTriangle_Integrate( tri )
    type( sTriangle ) :: tri
    
    integer :: i,j
    !
    !隣接表(nei2)を更新する。
    !
    do i=1, tri%nmodified
       j = tri%newnei2idx(i)
       tri%nei2(j) = tri%newnei2(i)
    enddo
    !
    !隣接表(nei)を更新する。
    !
    do i=1, tri%nmodified
       j = tri%newnei2idx(i)
       tri%nei(j) = tri%newnei(i)
    enddo
    !
    !ヒストグラムを更新する。
    !tri%newhistのポインタをhistに移し、旧histの内容を廃棄したいのだが、安全な方法がよくわからないので、まわりくどい書き方をする。
    !
    call done( tri%hist )
    call histogram_dup( tri%hist, tri%newhist )
    call done( tri%newhist )
    
#ifdef DEBUG2
    !
    !1trial後の積算ヒストグラムを、再計算して一致するかどうかを確認
    !
    if ( .true. ) then
       call new( tri%tmphist, 21, 0.1d0, -1d0 )
       do i=1, sys%mol(1)%nmol
          do j=1, tri%nei(i)%nnear
             do k=j+1, tri%nei(i)%nnear
                !
                ! O-O-O角の余弦を計算する。
                !
                costh = tri%nei(i)%dx(j)*tri%nei(i)%dx(k) + tri%nei(i)%dy(j)*tri%nei(i)%dy(k) + tri%nei(i)%dz(j)*tri%nei(i)%dz(k)
                costh = costh / sqrt ( tri%nei(i)%dd(j) * tri%nei(i)%dd(k) )
                call histogram_accumulate( tri%tmphist, costh, 1d0 )
             enddo
          enddo
       enddo
       if ( associated( tri%tmphist ) ) then
          diff  = kullback_div( tri%hist, tri%tmphist )
          if ( diff /= 0d0 ) then
             call histogram_show( tri%hist, STDERR )
             call histogram_show( tri%tmphist, STDERR )
             write(STDERR,*) diff, "DIV", loop, trial
             call die( 0, "TMP Histogram differs." )
          endif
       endif
       call done( tri%tmphist )
    endif
#endif /*DEBUG2*/
  end subroutine sTriangle_Integrate

  subroutine sTriangle_Reject( tri )
    type( sTriangle ) :: tri
    call done( tri%newhist )
  end subroutine sTriangle_Reject

  subroutine sTriangle_EndAStep( tri )
    type( sTriangle ) :: tri
#ifdef DEBUG
    call histogram_dup( tri%lasthist, tri%hist )
#endif
    call done( tri%hist )
  end subroutine sTriangle_EndAStep

end module triangle
