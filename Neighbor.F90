! -*- f90 -*-
module neighbors_module
  implicit none
  integer, parameter, private :: MAXNEAR=150
  integer, parameter :: ORDER_NONE = 0, ORDER_BY_DISTANCE = 1, ORDER_BY_LABEL = 2
  type sNeighbors
     integer :: nnear
     integer :: near( MAXNEAR+1 ) !Stackを使う
     integer :: order
     real( kind=8 ) :: dx( MAXNEAR+1 )
     real( kind=8 ) :: dy( MAXNEAR+1 )
     real( kind=8 ) :: dz( MAXNEAR+1 )
     real( kind=8 ) :: dd( MAXNEAR+1 )
  end type sNeighbors
  type sAdjacency
     real(kind=8) :: dx,dy,dz,dd
  end type sAdjacency
     
contains
  subroutine insert( nei, label, dx, dy, dz, dd, limit )
    use error_module
    type( sNeighbors ), intent(INOUT) :: nei
    integer, intent(IN)              :: label, limit
    real(kind=8), intent(IN)         :: dx, dy, dz
    integer :: i, j, max
    real(kind=8), intent(IN) :: dd
    !write(6,*)"ORDER", nei%order
    if ( ORDER_NONE == nei%order ) then
       !
       !if maximum number of the neighbor site is not specified
       ! == no sorting == append to the last
       !
       nei%nnear = nei%nnear + 1
#ifdef MEMUSAGE
       if ( MAXNEAR <= nei%nnear ) then
          write(STDERR,*) nei%nnear, MAXNEAR
          call die( 0, "NEIGHBORS_INSERT" )
       endif
#endif
       nei%dx(nei%nnear) = dx
       nei%dy(nei%nnear) = dy
       nei%dz(nei%nnear) = dz
       nei%dd(nei%nnear) = dd
       nei%near(nei%nnear)     = label
    else
       if ( ORDER_BY_DISTANCE == nei%order ) then
          do i=1, nei%nnear
             if ( dd < nei%dd(i) ) exit
          enddo
       else if ( ORDER_BY_LABEL == nei%order ) then
          do i=1, nei%nnear
             if ( label < nei%near(i) ) exit
          enddo
       endif
       !write(6,*) "===", i, limit
       !
       ! i番目に挿入する。
       !
       if ( i <= limit .or. limit == 0 ) then
          !write(6,*) (nei%near(j), nei%distance(j), j=1,nei%nnear )
          !write(6,*) ">>>", dd,i
          !
          !maxは挿入後のデータ数。ただしlimitを超えてはいけない。
          !
          max = nei%nnear + 1
          if ( limit /=0 .and.limit < max ) then
             max = limit
          else
             nei%nnear = max
          endif
          do j=max, i+1, -1
             nei%dx(j)   = nei%dx(j-1)
             nei%dy(j)   = nei%dy(j-1)
             nei%dz(j)   = nei%dz(j-1)
             nei%dd(j)   = nei%dd(j-1)
             nei%near(j) = nei%near(j-1)
          enddo
          nei%dx(i) = dx
          nei%dy(i) = dy
          nei%dz(i) = dz
          nei%dd(i) = dd
          nei%near(i)     = label
          !write(6,*) (nei%near(j), nei%distance(j), j=1,nei%nnear )
       endif
    endif
  end subroutine insert

  !
  !pair list --> neighbor list
  !
  subroutine neighbors_all2( nsite, nei, iv, coord, maxnei, order, radius, am )
    ! ivには、あらかじめ短いカットオフで求めた対表が与えられているものとする。
    ! radiusが与えられた場合はサイト間距離でさらに絞りこむ。
    !maxnei=0の場合は、無条件にすべて表に追加する。
    use common_module
    use interaction_module
    use vector_module
    type( sInteraction ),intent(IN) :: iv
    type( Vector3 ), intent(IN)     :: coord(*)
    type( sNeighbors ) :: nei(*)
    integer, intent(IN) :: nsite, maxnei, order
    real(kind=8), intent(IN), optional :: radius
    !type(sAdjacency), intent(INOUT), optional :: am(nsite,nsite)
    type(sAdjacency), intent(INOUT), optional :: am(:,:)
    integer :: i,j,k
    real(kind=8) :: dx,dy,dz,dd
    
    !write(6,*) "maxnei", maxnei
    !write(6,*) "maxnei", maxnei, nei%need_sort
    do i=1, nsite
       nei(i)%order = order
       nei(i)%nnear = 0
    enddo
    !write(6,*) "NPAIR",iv%npair
    do k=1,iv%Npair
       i = iv%pair_i(k)
       j = iv%pair_j(k)
       ! j-->i
       dx = coord(i)%vec(1) - coord(j)%vec(1)
       dy = coord(i)%vec(2) - coord(j)%vec(2)
       dz = coord(i)%vec(3) - coord(j)%vec(3)
       dx = dx - iv%ox(k)
       dy = dy - iv%oy(k)
       dz = dz - iv%oz(k)
       dd = dx**2 + dy**2 + dz**2
       !if ( i == 2 .and.j == 32 ) then
       !write(6,*) "pair", i,j, dx,dy,dz,dd
       !endif
       if ( present( radius ) ) then
          if ( radius**2 < dd ) cycle
       endif
       !
       !相対位置ベクトルは本当はivが保持しておいてほしいな。
       !
       call insert( nei(i), j, -dx, -dy, -dz, dd, maxnei )
       !call sneighbor1_write(nei(i),6)
       call insert( nei(j), i, dx, dy, dz, dd, maxnei )
       !if ( i.eq.10 .or. j.eq.10 ) then
       !   write(6,*) i,j
       !endif
       !write(6,*) j,nei(j)%nnear,"@"
       if ( present(am) ) then
          ! i-->j
          am(i,j)%dx = -dx
          am(i,j)%dy = -dy
          am(i,j)%dz = -dz
          am(i,j)%dd =  dd
          ! j-->i
          am(j,i)%dx = +dx
          am(j,i)%dy = +dy
          am(j,i)%dz = +dz
          am(j,i)%dd =  dd
       endif
    enddo
    do i=1,nsite
    !   write(6,*) "*", i, (nei(i)%near(j), nei(i)%distance(j), j=1,nei(i)%nnear )
       !write(6,*) "*", i, nei(i)%nnear, maxnei
    enddo
  end subroutine neighbors_all2


  subroutine neighbors_adjmatrix( nei, nmol, am )
    use common_module
    type( sNeighbors ) :: nei(*)
    integer, intent(IN) :: nmol
    type(sAdjacency), intent(INOUT) :: am(:,:)
    integer :: i,j,k
    !
    !q6s%hbondの内容は上書きする(rejectした場合は復旧必要)
    !
    am(:,:)%dx = 0d0
    am(:,:)%dy = 0d0
    am(:,:)%dz = 0d0
    am(:,:)%dd = 0d0
    do k=1, nmol
       do i=1, nei( k )%nnear
          j = nei( k )%near(i)
          ! k-->j
          am( k, j )%dx = nei( k )%dx(i)
          am( k, j )%dy = nei( k )%dy(i)
          am( k, j )%dz = nei( k )%dz(i)
          am( k, j )%dd = nei( k )%dd(i)
          ! j-->k
          am( j, k )%dx =-nei( k )%dx(i)
          am( j, k )%dy =-nei( k )%dy(i)
          am( j, k )%dz =-nei( k )%dz(i)
          am( j, k )%dd = nei( k )%dd(i)
       enddo
    enddo
  end subroutine neighbors_adjmatrix

  !
  !neighbor list-->よりコンパクトなneighbot list
  !
  subroutine neighbor1_compress2( newnei, nei, maxnei, order, radius )
    ! radiusが与えられた場合はサイト間距離でさらに絞りこむ。
    ! maxnei=0の場合は、無条件にすべて表に追加する。
    type(sNeighbors), intent(IN)  :: nei
    type(sNeighbors), intent(OUT) :: newnei
    integer, intent(IN) :: maxnei, order
    real(kind=8), intent(IN), optional :: radius

    integer :: i
    
    newnei%order = order
    newnei%nnear = 0
    !write(6,*) iv%npair
    do i=1, nei%nnear
       if ( present( radius ) ) then
          if ( radius**2 < nei%dd(i) ) cycle
       endif
       !write(6,*) "!!"
       call insert( newnei, nei%near(i), nei%dx(i), nei%dy(i), nei%dz(i), nei%dd(i), maxnei )
    enddo
  end subroutine neighbor1_compress2


  !
  !neiに登録されているサイトをソートしなおす。
  !とりあえずバブルソート
  !
  subroutine neighbors_bubblesort( nei )
    type( sNeighbors ), intent(INOUT) :: nei
    integer :: i,j
    integer :: n
    real(kind=8) :: x,y,z,d
    if ( nei%order == ORDER_BY_DISTANCE ) then
       do i=1, nei%nnear
          do j=1, i-1
             if ( nei%dd(i) < nei%dd(j) ) then
                n = nei%near(i)
                x = nei%dx(i)
                y = nei%dy(i)
                z = nei%dz(i)
                d = nei%dd(i)
                nei%near(i) = nei%near(j)
                nei%dx(i)   = nei%dx(j)
                nei%dy(i)   = nei%dy(j)
                nei%dz(i)   = nei%dz(j)
                nei%dd(i)   = nei%dd(j)
                nei%near(j) = n
                nei%dx(j)   = x
                nei%dy(j)   = y
                nei%dz(j)   = z
                nei%dd(j)   = d
             endif
          enddo
       enddo
    else if ( nei%order == ORDER_BY_LABEL ) then
       do i=1, nei%nnear
          do j=1, i-1
             if ( nei%near(i) < nei%near(j) ) then
                n = nei%near(i)
                x = nei%dx(i)
                y = nei%dy(i)
                z = nei%dz(i)
                d = nei%dd(i)
                nei%near(i) = nei%near(j)
                nei%dx(i)   = nei%dx(j)
                nei%dy(i)   = nei%dy(j)
                nei%dz(i)   = nei%dz(j)
                nei%dd(i)   = nei%dd(j)
                nei%near(j) = n
                nei%dx(j)   = x
                nei%dy(j)   = y
                nei%dz(j)   = z
                nei%dd(j)   = d
             endif
          enddo
       enddo
       endif
  end subroutine neighbors_bubblesort


  !
  !neiに登録されているサイトをソートしなおす。
  !insert sort
  !
  subroutine neighbors_sort( nei )
    type( sNeighbors ), intent(INOUT) :: nei
    integer :: i,j,k
    integer :: n
    real(kind=8) :: x,y,z,d
    if ( nei%order == ORDER_BY_DISTANCE ) then
       do i=1, nei%nnear
          do j=1, i-1
             if ( nei%dd(i) < nei%dd(j) ) then
                n = nei%near(i)
                x = nei%dx(i)
                y = nei%dy(i)
                z = nei%dz(i)
                d = nei%dd(i)
                do k=i, j+1, -1
                   nei%near(k) = nei%near(k-1)
                   nei%dx(k)   = nei%dx(k-1)
                   nei%dy(k)   = nei%dy(k-1)
                   nei%dz(k)   = nei%dz(k-1)
                   nei%dd(k)   = nei%dd(k-1)
                enddo
                nei%near(j) = n
                nei%dx(j)   = x
                nei%dy(j)   = y
                nei%dz(j)   = z
                nei%dd(j)   = d
             endif
          enddo
       enddo
    else if ( nei%order == ORDER_BY_LABEL ) then
       do i=1, nei%nnear
          do j=1, i-1
             if ( nei%near(i) < nei%near(j) ) then
                n = nei%near(i)
                x = nei%dx(i)
                y = nei%dy(i)
                z = nei%dz(i)
                d = nei%dd(i)
                do k=i, j+1, -1
                   nei%near(k) = nei%near(k-1)
                   nei%dx(k)   = nei%dx(k-1)
                   nei%dy(k)   = nei%dy(k-1)
                   nei%dz(k)   = nei%dz(k-1)
                   nei%dd(k)   = nei%dd(k-1)
                enddo
                nei%near(j) = n
                nei%dx(j)   = x
                nei%dy(j)   = y
                nei%dz(j)   = z
                nei%dd(j)   = d
             endif
          enddo
       enddo
       endif
  end subroutine neighbors_sort
    
  !
  !こっちの方が意外に遅かった。
  !
  subroutine neighbors_sort2( nei )
    type( sNeighbors ), intent(INOUT) :: nei
    integer :: i,j,k,n
    integer :: moved, last
    real(kind=8) :: x,y,z,d
    type( sNeighbors ) :: temp
    if ( nei%order == ORDER_BY_DISTANCE ) then
       do i=1, nei%nnear
          do j=1, i-1
             if ( nei%dd(i) < nei%dd(j) ) then
                n = nei%near(i)
                x = nei%dx(i)
                y = nei%dy(i)
                z = nei%dz(i)
                d = nei%dd(i)
                do k=i, j+1, -1
                   nei%near(k) = nei%near(k-1)
                   nei%dx(k)   = nei%dx(k-1)
                   nei%dy(k)   = nei%dy(k-1)
                   nei%dz(k)   = nei%dz(k-1)
                   nei%dd(k)   = nei%dd(k-1)
                enddo
                nei%near(j) = n
                nei%dx(j)   = x
                nei%dy(j)   = y
                nei%dz(j)   = z
                nei%dd(j)   = d
             endif
          enddo
       enddo
    else if ( nei%order == ORDER_BY_LABEL ) then
       n = nei%nnear
       do moved=1, n
          j = 1
          last = n - ( moved - 1 )
          do i=2, last
             if ( nei%near(i) < nei%near(j) ) then
                j = i
             endif
          enddo
          temp%near(moved) = nei%near(j)
          temp%dx(moved)   = nei%dx(j)
          temp%dy(moved)   = nei%dy(j)
          temp%dz(moved)   = nei%dz(j)
          temp%dd(moved)   = nei%dd(j)
          nei%near(j)      = nei%near(last)
          nei%dx(j)        = nei%dx(last)
          nei%dy(j)        = nei%dy(last)
          nei%dz(j)        = nei%dz(last)
          nei%dd(j)        = nei%dd(last)
       enddo
       nei%near(1:n) = temp%near(1:n)
       nei%dx(1:n)   = temp%dx(1:n)
       nei%dy(1:n)   = temp%dy(1:n)
       nei%dz(1:n)   = temp%dz(1:n)
       nei%dd(1:n)   = temp%dd(1:n)
    endif
  end subroutine neighbors_sort2
    
  !
  !target分子がdeltaxyzだけ変位することで、隣接分子表が更新される。
  !newneiには、更新された分子の情報だけが圧縮されて入る。(時には不便な仕様)
  !
  subroutine neighbors_displacement( nei, target, deltax, deltay, deltaz, nmodified, newnei, newneiidx )
    type(sNeighbors), intent(IN)  :: nei(*)
    integer, intent(IN)           :: target
    real(kind=8), intent(IN)      :: deltax,deltay,deltaz
    integer, intent(OUT)          :: nmodified
    type(sNeighbors), intent(OUT) :: newnei(*)
    integer, intent(OUT)          :: newneiidx(*)

    integer :: i,j,k
    logical :: isModified
    !
    !中心分子が移動した時に、近傍分子までの相対距離が変化する。
    !それらを計算しなおす。移動量が小さく、bookkeeping margin以下である
    !(から、隣接分子表を作りなおす必要はない)ことを仮定する。
    !
    !newnei(1)には中心分子を指定する。
    !neiとnewneiの添字の対応表はnewneiidxに記録する。
    !
    newnei(1) = nei( target )
    newneiidx(1) = target

    !
    !中心分子が移動したので、すべての距離が変化する。
    !
    do i=1, newnei(1)%nnear
       newnei(1)%dx(i) = newnei(1)%dx(i) - deltax
       newnei(1)%dy(i) = newnei(1)%dy(i) - deltay
       newnei(1)%dz(i) = newnei(1)%dz(i) - deltaz
       newnei(1)%dd(i) = newnei(1)%dx(i)**2 + newnei(1)%dy(i)**2 + newnei(1)%dz(i)**2
    enddo
    !
    !つぎに、targetの周辺分子について
    !周辺分子から見た、targetの位置を補正。
    !
    nmodified = 2
    do j=1, newnei( 1 )%nnear
       k = newnei( 1 )%near(j)
       !write(6,*) j,newnei( 1 )%nnear !, nei(target)%nnear
       newnei( nmodified )    = nei( k )
       newneiidx( nmodified ) = k
       ismodified = .false.
       i = sNeighbor1_Lookup( newnei( nmodified ), target )
       if ( 0 < i ) then
          newnei( nmodified )%dx(i) = newnei( nmodified )%dx(i) + deltax
          newnei( nmodified )%dy(i) = newnei( nmodified )%dy(i) + deltay
          newnei( nmodified )%dz(i) = newnei( nmodified )%dz(i) + deltaz
          newnei( nmodified )%dd(i) = newnei( nmodified )%dx(i)**2 + newnei( nmodified )%dy(i)**2 + newnei( nmodified )%dz(i)**2
          isModified = .true.
          nModified = nModified + 1
       endif
    enddo
    nModified = nModified - 1
    !write(6,*) nmodified, "nmodified", newnei( 1 )%nnear, nei(target)%nnear
    do i=1, nModified
       if ( newnei( i )%order /= ORDER_NONE ) call neighbors_sort( newnei( i ) )
    enddo
  end subroutine neighbors_displacement

  subroutine sneighbor1_write( nei, file )
    type(sNeighbors) :: nei
    integer :: file

    integer :: i
    write(file,*) "nnear=", nei%nnear
    write(file,*) "order=", nei%order
    write(file,*) (nei%near(i), i=1, nei%nnear)
    do i=1,nei%nnear
       write(file,*) nei%dx(i)
    enddo
  end subroutine sneighbor1_write

  subroutine sneighbor1_compare( n1, n2, file )
    type(sNeighbors) :: n1, n2
    integer :: file

    integer :: i
    if ( n1%nnear /= n2%nnear ) then
       write(file,*) n1%nnear, n2%nnear, " differs(nnear)."
    endif
    do i=1, n1%nnear
       if ( n1%near(i) /= n2%near(i) ) then
          write(file,*) i, n1%near(i), n2%near(i), " differs(near(i))."
       endif
    enddo
    do i=1,n1%nnear
       if ( n1%near(i) /= n2%near(i) ) then
          write(file,*) i, n1%dd(i),n2%dd(i),abs((n1%dd(i)-n2%dd(i))/n2%dd(i)), "% differs(dd(i))."
       endif
    enddo
  end subroutine sneighbor1_compare

  function sneighbor1_lookup( nei, label )
    integer :: sneighbor1_lookup
    integer, intent(IN) :: label
    type(sNeighbors) :: nei

    integer :: i
    do i=1, nei%nnear
       if ( nei%near(i) == label ) then
          sneighbor1_lookup = i
          return
       endif
    enddo
    sneighbor1_lookup = -1
  end function sneighbor1_lookup


end module neighbors_module
