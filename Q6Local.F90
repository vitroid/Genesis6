! -*- f90 -*-
#define DEBUG
#undef DEBUG5
!
!平成17年3月10日(木)新しいアルゴリズムではキャッシュする必要はほとんどない。
!キャッシュしない方が2/3程度に計算時間を短縮できる。
!
#undef USE_Q6CACHE
!
!最終デバッグ平成17年3月9日(水)時のメッセージを表示させるなら
!
#undef DEBUG20050309


!Q6でUmbrella MCを行うための算程群
!ヒストグラムを使わず、期待値をそのままOPとする。
!
!5A以内の、すべての水素結合(距離3A以下)の方位からQ6を計算する。
!これまで(平成17年2月23日(水))使っていたQ6Tは、切り出し半径6Aとしていたが、
!計算時間やMC特有の都合(変分計算のためのマージンが必要になる)で短縮した。
!
module q6local2
  !use set_module
  use neighbors_module
  !use distrib_module
  !use umbrella
  use common_module
  use set_module
#ifdef USE_Q6CACHE
#else
  use q6_module
#endif
  implicit none

  integer, parameter :: MAXBOND=10
  type sq6system
     type(sNeighbors) :: nei6(MAXMOL) ! 大元の隣接分子表
     type(sNeighbors) :: nei3(MAXMOL) ! HB距離。nei6から生成。
     type(sNeighbors) :: nei5(MAXMOL) ! Q6のレンジ。nei6から生成。
     real(kind=8)     :: q6( MAXMOL ), sumq6

     !for 2
     complex(kind=8)  :: qlmsum( -6:6, MAXMOL )  ! foreach mol
     complex(kind=8)  :: qlm( -6:6, MAXBOND, MAXMOL ) ! foreach bond
     type(sSet)       :: neis( MAXMOL )          ! ある分子から5A以内の分子
     type(sSet)       :: neib( MAXBOND, MAXMOL )      ! ある結合から5A以内の分子
     integer          :: nbond( MAXMOL )         ! ある分子から5A以内の結合の数
  end type sq6system

  type sQ6diff
     type(sNeighbors) :: newnei3(99)
     type(sNeighbors) :: newnei5(99)
     type(sNeighbors) :: newnei6(99)
     integer          :: newnei6idx(99)
     integer          :: nei6modified
     !for 2
     real(kind=8)     :: q6( MAXMOL ), sumq6
     complex(kind=8)  :: qlm( -6:6, MAXBOND ) ! foreach bond
     type(sSet)       :: neis            ! foreach mol
     type(sSet)       :: plus            ! foreach mol
     type(sSet)       :: minus           ! foreach mol
     type(sSet)       :: neib( MAXBOND )      ! foreach bond
  end type sQ6diff

  type sq6systemPtr
     type(sq6system), pointer :: p
  end type sq6systemPtr

  type sq6diffPtr
     type(sq6diff), pointer :: p
  end type sq6diffPtr

  interface new
     module procedure sq6system_new
  end interface

  interface done
     module procedure sq6system_done
  end interface

contains

  subroutine sq6system_new( tri )
    type( sq6system ), pointer :: tri
    allocate( tri )
    !
    !for qlmcache 10000で十分。1000だと足りない。大きすぎるとまた遅くなる。
    !10000から50000に増やすと、isequalの呼出しが増えるのでかえって遅くなるかも。
    !10000から5000に減らしてもやはりisequalの呼出しは増える。が、migrateできる場合も。
#ifdef USE_Q6CACHE
    call qlmcache_begin(5000)
#endif
  end subroutine sq6system_new

  subroutine sq6system_done( tri )
    type( sq6system ), pointer :: tri
    deallocate( tri )
  end subroutine sq6system_done


  subroutine sq6system_prepareall( q6s, nmol, ww, com, op )
    use interaction_module
    use vector_module
#ifdef USE_Q6CACHE
    use q6cache_module
#else
    use q6_module
#endif
    type(sq6system), pointer :: q6s
    integer, intent(IN) :: nmol
    type( sInteraction ), intent(IN) :: ww
    type(vector3), intent(IN) :: com( * )
    real(kind=8), intent(OUT), optional :: op

    complex(kind=8) :: dummy
    integer :: i,j,k,jj,kk,m
    real(kind=8) :: costh, newq6
    real(kind=8) :: pp, px, py, pz
    !
    !ある分子の近傍分子を記録しておく。
    !
    !call neighbors_all( nmol, q6s%nei3, ww, com, 0, 3d0, q6s%hbond )
    call neighbors_all2( nmol, q6s%nei6, ww, com, 0, ORDER_BY_LABEL, 6d0 )
    !
    !nei6を元に、必要に応じて生成する。
    !
    do i=1, nmol
       call neighbor1_compress2( q6s%nei3(i), q6s%nei6(i), 0, ORDER_BY_LABEL, 3d0 )
       call neighbor1_compress2( q6s%nei5(i), q6s%nei6(i), 0, ORDER_BY_LABEL, 5d0 )
#ifdef DEBUG20050309
       write(6,*)"NNEAR", q6s%nei6(i)%nnear,q6s%nei5(i)%nnear,q6s%nei3(i)%nnear
#endif
    enddo
    !
    !すべての分子について
    !
    do i=1,nmol
       !
       !その分子の近傍分子(自分を含む)の集合を作る。
       !
       call new( q6s%neis(i), nmol, q6s%nei5(i)%nnear, q6s%nei5(i)%near, .false. )
       call join( q6s%neis(i), i )
    enddo
    do i=1,nmol
       !
       !分子iのすべての結合について
       !
       do jj=1, q6s%nei3(i)%nnear
          j = q6s%nei3(i)%near(jj)
          !
          !ij結合を近傍に含む分子の集合を作る。
          !
          call intersection( q6s%neis(i), q6s%neis(j), q6s%neib( jj, i ) )
          !
          !ij結合のqlm値をあらかじめ計算しておく。
          !
          pp = 1d0 / sqrt( q6s%nei3(i)%dd(jj) )
          px = q6s%nei3(i)%dx(jj) * pp
          py = q6s%nei3(i)%dy(jj) * pp
          pz = q6s%nei3(i)%dz(jj) * pp
          do m=-6,6
#ifdef USE_Q6CACHE
             call qlm_cache( 6, m, px, py, pz, q6s%qlm( m, jj, i ) )
#else
             q6s%qlm( m, jj, i )  = qlm( 6, m, px, py, pz)
#endif
          enddo
          !write(6,*)q6s%qlm( -6, jj, i )
          !if ( i == 1 ) then
          !   write(STDERR,*) i,j,px,py,pz
          !endif
          !q6s%qlm( m, jj, i ) = qlm( 6, m, px, py, pz )
       enddo
    enddo

    !
    !qlmsum、分子ごとの、qlm値の合計を計算する。
    !
    q6s%qlmsum(:,:) = 0d0
    q6s%nbond(1:nmol) = 0
    !
    !実際には、ループは結合ごとに回す。
    !
    do i=1,nmol
       !
       !分子iのすべての結合について
       !
       do jj=1, q6s%nei3(i)%nnear
          j = q6s%nei3(i)%near(jj)
          !
          !ij結合を近傍に含む分子について
          !
          do kk=1, q6s%neib( jj, i )%size
             k = q6s%neib( jj, i )%set(kk)
             q6s%nbond(k) = q6s%nbond(k) + 1
#ifdef DEBUG20050309
             if ( k==1 .and. i<j ) write(STDOUT,*) "NEIB",i,j !,q6s%qlm( -6, jj, i )
#endif
             do m=-6,6
                q6s%qlmsum(m,k) = q6s%qlmsum(m,k) + q6s%qlm( m, jj, i )
             enddo
          enddo
       enddo
    enddo
    !
    !Average all q6
    !
    q6s%sumq6 = 0d0
    do i=1, nmol
       newq6 = 0d0
       do m=-6,6
          newq6 = newq6 + abs(q6s%qlmsum(m,i))**2
       enddo
       newq6 = sqrt( newq6 * 4d0 * pi / ( 2d0 * 6 + 1d0 ) ) / q6s%nbond(i)
#ifdef DEBUG20050309
       write(6,*) "Q6:", q6s%nbond(i), newq6
#endif
       q6s%q6(i) = newq6
       q6s%sumq6 = q6s%sumq6 + q6s%q6(i)
    enddo
#ifdef USE_Q6CACHE
    call dict_stat()
#endif
    if ( present( op ) ) then
    !write(6,*) "q6:", q6s%sumq6 / nmol
       op = q6s%sumq6 / nmol
    endif
    !
    !ここまではQ6Local.F90と完全に同一動作を確認 平成17年3月7日(月)
    !
  end subroutine sq6system_prepareall
    
  !
  !dirtyなベクトルのみをリストし、差分計算を試みる。平成17年3月2日(水)
  !
  subroutine sq6system_Differentiate( q6s, target, deltax, deltay, deltaz, nmol, q6d, op )
#ifdef USE_Q6CACHE
    use q6cache_module
#else
    use q6_module
#endif
    type( sq6system ) :: q6s
    type( sq6diff )   :: q6d
    integer, intent(IN) :: target
    real(kind=8), intent(IN) :: deltax,deltay,deltaz
    integer, intent(IN) :: nmol
    real(kind=8), intent(OUT) :: op
    
    integer :: i, j, k, ii, jj, kk, m
    real(kind=8) :: px, py, pz, pp, newq6

#ifdef DEBUG20050309
    call sq6system_write( q6s, nmol,STDOUT )
#endif
    !
    !移動前に、targetに水素結合していた相手は
    !
#ifdef DEBUG20050309
    write(STDERR,*) "OLD NEAR", target, ":", ( q6s%nei3(target)%near(i), i=1, q6s%nei3(target)%nnear )
#endif
    !
    !これらの結合を近傍とする分子の、qlmsumからqlm値を差し引く
    !
    do ii=1, q6s%nei3(target)%nnear
       i = q6s%nei3(target)%near(ii)
       !
       !target-i結合の近傍の分子kについて
       !
       do kk=1, q6s%neib( ii, target )%size
          k = q6s%neib( ii, target )%set(kk)
          !
          !nbondとqlmsumはq6sの値を直接変更する(rejectされた場合に復旧が必要)
          !
          q6s%nbond(k) = q6s%nbond(k) - 2
          do m=-6,6
             q6s%qlmsum(m,k) = q6s%qlmsum(m,k) - q6s%qlm( m, ii, target ) * 2
#ifdef DEBUG20050309
             if(  m==-6 ) then
                write(6,*) -target, i,k,q6s%qlm( m, ii, target )
             endif
#endif
          enddo
       enddo
    enddo
    !
    !移動後のtargetの位置で、隣接ベクトルを準備する。
    !
#ifdef DEBUG20050309
    write(6,*) "d",deltax,deltay,deltaz
#endif
    call neighbors_displacement( q6s%nei6, target, deltax, deltay, deltaz, q6d%nei6modified, q6d%newnei6, q6d%newnei6idx )
    do i=1, q6d%nei6modified
       call neighbor1_compress2( q6d%newnei3(i), q6d%newnei6(i), 0, ORDER_BY_LABEL, 3d0 )
       call neighbor1_compress2( q6d%newnei5(i), q6d%newnei6(i), 0, ORDER_BY_LABEL, 5d0 )
    enddo
    !
    !移動後に、targetに水素結合している相手は
    !
#ifdef DEBUG20050309
    write(STDERR,*) "NEW NEAR", target, ":", ( q6d%newnei3(1)%near(i), i=1, q6d%newnei3(1)%nnear )
#endif

    !
    !移動後におこること
    !まずtargetが結合する相手が替わる。
    !変化した結合に隣接する分子のリストが替わる
    !隣接分子のqlmsum値が変化する。
    
    !targetの隣接分子集合を再構成する。
    !ほかの分子は、移動していないので、ほかの分子の隣接分子集合は不変
    call new( q6d%neis, nmol, q6d%newnei5(1)%nnear, q6d%newnei5(1)%near, .false. )
    call join( q6d%neis, target )
    !
    !新たに近傍に入ってきた分子の集合
    !
    call difference( q6d%neis, q6s%neis(target), q6d%plus )
    !
    !近傍から出ていった分子の集合
    !
    call difference( q6s%neis(target), q6d%neis, q6d%minus )

    !
    !まず、targetの直近分子との間にある結合が移動することによる変化の計算
    !
    !
    !1.targetと、直近隣接分子の隣接分子集合の交差を更新する。
    do jj=1, q6d%newnei3(1)%nnear
       j = q6d%newnei3(1)%near(jj)
       call intersection( q6d%neis, q6s%neis(j), q6d%neib( jj ) )
       !
       !この結合のqlm値をあらかじめ計算しておく。
       !
       pp = 1d0 / sqrt( q6d%newnei3(1)%dd(jj) )
       px = q6d%newnei3(1)%dx(jj) * pp
       py = q6d%newnei3(1)%dy(jj) * pp
       pz = q6d%newnei3(1)%dz(jj) * pp
       do m=-6,6
#ifdef USE_Q6CACHE
          call qlm_cache( 6, m, px, py, pz, q6d%qlm( m, jj ) )
#else
          q6d%qlm( m, jj )  = qlm( 6, m, px, py, pz )
#endif
#ifdef DEBUG20050309
          if ( m == -6 ) then
             write(6,*) target, j
             write(6,*) px,py,pz
             do k=1,q6s%nei3(target)%nnear
                if ( q6s%nei3(target)%near(k) == j ) then
                   write(6,*) q6s%nei3(target)%dx(k)/sqrt(q6s%nei3(target)%dd(k)),&
                        q6s%nei3(target)%dy(k)/sqrt(q6s%nei3(target)%dd(k)),&
                        q6s%nei3(target)%dz(k)/sqrt(q6s%nei3(target)%dd(k))
                endif
             enddo
          endif
#endif
       enddo
    enddo
    !
    !2.各分子のqlmsumに、新しいqlm値を加算する。
    !
    do jj=1, q6d%newnei3(1)%nnear
       !
       !jはtargetと水素結合している。
       !
       j = q6d%newnei3(1)%near(jj)
       !
       !target-j結合から5A以内の分子kについて
       !
       do kk=1, q6d%neib( jj )%size
          k = q6d%neib( jj )%set(kk)
          q6s%nbond(k) = q6s%nbond(k) + 2
          do m=-6,6
             q6s%qlmsum( m, k ) = q6s%qlmsum( m, k ) + q6d%qlm( m, jj ) * 2
             !if( 32 == k .and. m==-6 ) then
#ifdef DEBUG20050309
             if(  m==-6 ) then
                write(6,*) target, j,k,q6d%qlm( m, jj )
             endif
#endif
          enddo
       enddo
    enddo
    !
    !近傍に出現した分子は、もしかしたら結合も持ちこむかもしれない。
    !新近傍分子i⊂{Q}について
    do ii=1, q6d%plus%size
       i = q6d%plus%set(ii)
#ifdef DEBUG20050309
       write(6,*) "PLUS", i
#endif
       !
       !それと結合する分子jについて
       !
       do jj=1, q6s%nei3(i)%nnear
          j = q6s%nei3(i)%near(jj)
          !
          !もし、jもtargetの5A近傍であれば、i-j結合の近傍にtargetがあることになる。
          !
          if ( set_exist( q6d%neis, j ) ) then
             !
             !その場合は、targetのq値に加算する。
             !
             q6s%nbond(target) = q6s%nbond(target) + 2
             do m=-6,6
                q6s%qlmsum( m, target ) = q6s%qlmsum( m, target ) + q6s%qlm( m, jj, i ) * 2
             enddo
          endif
       enddo
    enddo
    !
    !近傍から消失した分子は、結合も確実に失う。
    !消失分子i⊂{P}について
    do ii=1, q6d%minus%size
       i = q6d%minus%set(ii)
#ifdef DEBUG20050309
       write(6,*) "MINUS", i
#endif
       !
       !それと結合する分子j⊂{P,R}について
       !
       do jj=1, q6s%nei3(i)%nnear
          j = q6s%nei3(i)%near(jj)
          !
          !もし、jももともとtargetの5A近傍であれば、i-j結合の近傍にtargetがあったことになる。
          !
          if ( exist( q6s%neis(target), j ) ) then
             !
             !その場合は、targetのq値から減算する。
             !
#ifdef DEBUG20050309
             write(6,*) "MINUS-", j
#endif
             q6s%nbond(target) = q6s%nbond(target) - 2
             do m=-6,6
                q6s%qlmsum( m, target ) = q6s%qlmsum( m, target ) - q6s%qlm( m, jj, i ) * 2
             enddo
          endif
       enddo
    enddo
    !
    !Average all q6
    !
    q6d%sumq6 = 0d0
    do i=1, nmol
       newq6 = 0d0
       do m=-6,6
          newq6 = newq6 + abs(q6s%qlmsum(m,i))**2
       enddo
       newq6 = sqrt( newq6 * 4d0 * pi / ( 2d0 * 6 + 1d0 ) ) / q6s%nbond(i)
#ifdef DEBUG20050309
       write(6,*) i, "nbond", q6s%nbond(i),"Q6", newq6
#endif
       q6d%q6(i) = newq6
       q6d%sumq6 = q6d%sumq6 + q6d%q6(i)
    enddo
    
    op = q6d%sumq6 / nmol
    do ii=1, q6d%nei6modified
       i = q6d%newnei6idx(ii)
#ifdef DEBUG20050309
       write(STDOUT,'(a3,2i3,f8.5)') "Q6 ",target,i,q6d%q6(i)
#endif
    enddo

  end subroutine sq6system_Differentiate
  
  subroutine sq6system_Integrate( q6s, q6d, nmol, target )
    type( sq6system ) :: q6s
    type( sq6diff )   :: q6d
    integer, intent(IN) :: nmol, target
    
    type(sSet) :: oldnei3,newnei3,lost,kept,found
    integer :: i,j,k,ii,jj,kk
    complex(kind=8) :: qlmtemp( -6:6, MAXBOND ), qlmtemp2( -6:6 )
#ifdef DEBUG5
    type(sNeighbors) :: test(MAXMOL)
#endif    
    !
    !q6を更新する
    !
    do i=1, nmol
       q6s%q6(i) = q6d%q6(i)
    enddo
    q6s%sumq6 = q6d%sumq6
    !
    !target自身のqlm配列を更新する。
    !
    do jj=1, q6d%newnei3(1)%nnear
       !q6d%qlm( m, jj )はtargetのjj番目の水素結合のqlm。
       !
       !jはtargetと水素結合している。
       !
#ifdef DEBUG20050309
       write(6,*) "UPDATED", target, q6d%newnei3(1)%near(jj)
#endif
       q6s%qlm( -6:6, jj, target ) = q6d%qlm( -6:6, jj )
    enddo
    !
    !jの隣りについて。
    !targetは動くので、結合が切れたりあらたに結合したりする。
    !過去に結合していたもので、結合が切れた場合はqlm表をつめなければいけないし、
    !新たに結合した場合はqlm表に追加しなければいけない。
    !どちらの場合も、newnei3内の順序がnei3内の順序とは何ら関係がないので、
    !並べ替えをおこなわなければいけない。
    !
    !
    !neibを安全に更新するのはかなり大変。neisから再生しよう。
    !
    !1.target自身のneisの更新
    !
    q6s%neis(target) = q6d%neis
    !
    !2.targetの移動により、targetが視野に入った/出ていった、辺縁部の分子のneisも変化する。
    !
    do ii=1, q6d%plus%size
       !
       !iは新たに近傍になった分子。
       !
       i = q6d%plus%set(ii)
       !targetをiの近傍に登録する。
       call join( q6s%neis(i), target )
    enddo
    do ii=1, q6d%minus%size
       !
       !iは近傍から転出した分子。
       !
       i = q6d%minus%set(ii)
       !targetをiの近傍から除く。
       call drop( q6s%neis(i), target )
    enddo
    !
    !neibは集合演算だけなので処理は速い。
    !ここで再計算してしまう。
    !
    !targetに結合していたが、移動により切れてしまった結合のneib
    !
    call new( oldnei3, MAXBOND, q6s%nei3(target)%nnear,q6s%nei3(target)%near,.false.)
    call new( newnei3, MAXBOND, q6d%newnei3(1)%nnear,  q6d%newnei3(1)%near,  .false.)
    call difference( oldnei3, newnei3, lost )
    call difference( newnei3, oldnei3, found )
    call intersection( oldnei3, newnei3, kept )
    !あとの計算は、nei3の更新のあとのほうが便利なのでうしろで行う。

    !
    !新たに近傍に入った/出た分子の、neibを更新する？
    !
    do ii=1, q6d%plus%size
       !
       !iは新たに近傍になった分子。
       !
       i = q6d%plus%set(ii)
       do jj=1, q6s%nei3(i)%nnear
          !
          !jはiと結合している分子
          !
          j = q6s%nei3(i)%near(jj)
          !
          !もしjもtargetの近傍であれば
          !
          if ( exist( q6d%neis, j ) ) then
             !
             !結合i-jのneibには、targetを含む。
             !
             call join( q6s%neib( jj, i ), target )
             kk = sNeighbor1_Lookup( q6s%nei3(j), i )
             if ( 0 < kk ) then
                call join( q6s%neib( kk, j ), target )
             endif
          endif
       enddo
    enddo
    do ii=1, q6d%minus%size
       !
       !iは近傍から転出した分子
       !
       i = q6d%minus%set(ii)
       do jj=1, q6s%nei3(i)%nnear
          !
          !jはiと結合している分子
          !
          j = q6s%nei3(i)%near(jj)
          !
          !jが、targetの移動前には近傍にあったなら
          !
          if ( exist( q6s%neis(target), j ) ) then
             call drop( q6s%neib( jj, i ), target )
             kk = sNeighbor1_Lookup( q6s%nei3(j), i )
             if ( 0 < kk ) then
                call drop( q6s%neib( kk, j ), target )
             endif
          endif
       enddo
    enddo

    !
    !隣接表を更新する。
    !
    do i=1, q6d%nei6modified
       j = q6d%newnei6idx(i)
       q6s%nei6(j) = q6d%newnei6(i)
    enddo


    
    !
    !G⊂{lost}はもともとtargetと結合したが移動により結合が切れた。その結合のqlmは除き、Gの残りの結合のqlmを再構築する必要がある。targetとの結合以外は動いていないので、値は旧い値をそのまま使えるはずだが、対応関係が失われているのが難儀。
    !
    !
    !targetの隣接分子のqlmの更新。nei3がラベル順にソートされているものとみなして処理を簡約にする。
    !
    !まず、もともと結合していたが、結合が切れてしまった、lostのメンバーについて
    !
#ifdef DEBUG20050309
    write(6,*) "oldnei3"
    call sset_write(oldnei3, stdout )
    write(6,*) "newnei3"
    call sset_write(newnei3, stdout )
    write(6,*) lost%size, "lost, ", kept%size, "kept, ", found%size, " found."
#endif
    do ii=1, lost%size
       i = lost%set(ii)
       jj = 1
       !
       !nnearは移動前の隣接分子数。
       !

       do kk=1, q6s%nei3(i)%nnear
          k = q6s%nei3(i)%near(kk)
#ifdef DEBUG20050309
          write(6,*) "LOST", i,k
#endif
          if ( k /= target ) then
             q6s%qlm( -6:6, jj, i ) = q6s%qlm( -6:6, kk, i )
             jj = jj + 1
          endif
       enddo
    enddo
    !
    !結合が切れも継がりもしない隣接分子について
    !
    do ii=1, kept%size
       i = kept%set(ii)
       !
       !target-i結合のqlmをあらかじめ調べておく。
       !
       jj = sNeighbor1_Lookup( q6d%newnei3(1), i )
       if ( 0 < jj ) then
          qlmtemp2( -6:6 ) = q6d%qlm( -6:6, jj )
       endif
       kk = sNeighbor1_Lookup( q6s%nei3(i), target )
       if ( 0 < kk ) then
          k = q6s%nei3(i)%near(kk)
#ifdef DEBUG20050309
          write(6,*) "KEPT", i,k
#endif
          !
          !target自身用に計算したqlm値を使用する。
          !
          q6s%qlm( -6:6, kk, i ) = qlmtemp2( -6:6 )
       endif
    enddo
    !
    !新たに結合ができた場合は？
    !qlmの旧い値をコピーしつつ、targetとの結合に関してはq3dの値を使う。
    !
    !
    !iは、targetの6A近傍の分子
    !
    do ii=1, q6d%nei6modified
       i = q6d%newnei6idx(ii)
       !
       !もしi⊂{found}なら
       !
       if ( exist( found, i ) ) then
          !qlmを再配置するので、一時配列に逃がす。
          !
          qlmtemp( -6:6, 1:q6d%newnei3(ii)%nnear-1 ) = q6s%qlm( -6:6, 1:q6d%newnei3(ii)%nnear-1, i )
          !
          !target-i結合のqlmをあらかじめ調べておく。
          !
          jj = sNeighbor1_Lookup( q6d%newnei3(1), i )
          if ( 0 < jj ) then
             qlmtemp2( -6:6 ) = q6d%qlm( -6:6, jj )
          endif
          !
          !移動後の配置でのiの隣接分子{k}(targetを含む)について
          !
          jj = 1
          do kk=1, q6d%newnei3(ii)%nnear
             k = q6d%newnei3(ii)%near(kk)
#ifdef DEBUG20050309
             write(6,*) "FOUND", i,k
#endif
             !
             !もともとはtargetは隣接していなかったので、qlm値は計算されていない
             !
             if ( k /= target ) then
#ifdef DEBUG20050309
                write(6,*) "(FOUND)=", kk,i,jj
#endif
                q6s%qlm( -6:6, kk, i ) = qlmtemp( -6:6, jj )
                jj = jj + 1
             else
                !
                !target自身用に計算したqlm値を使用する。
                !
#ifdef DEBUG20050309
                write(6,*) "(FOUND)+", kk,i
#endif
                q6s%qlm( -6:6, kk, i ) = qlmtemp2( -6:6 )
             endif
          enddo
       endif
    enddo
                
    !
    !nei3を更新する。
    !更新する必要のないものまで更新している。
    !
    do i=1, q6d%nei6modified
       j = q6d%newnei6idx(i)
       q6s%nei3(j) = q6d%newnei3(i)
    enddo
    do i=1, q6d%nei6modified
       j = q6d%newnei6idx(i)
       q6s%nei5(j) = q6d%newnei5(i)
    enddo

    !
    !移動後の座標で、targetに直結した分子のneib
    !nei3が更新されてからでないと面倒なのでここにもってきた。
    !
    do jj=1, q6s%nei3(target)%nnear
       j = q6s%nei3(target)%near(jj)
       call intersection( q6s%neis(target), q6s%neis(j), q6s%neib( jj, target) )
       do kk=1, q6s%nei3(j)%nnear
          k = q6s%nei3(j)%near(kk)
          call intersection( q6s%neis(j), q6s%neis(k), q6s%neib( kk, j) )
       enddo
    enddo
    !
    !もともとtargetと結合していたが移動により結合が切れた相手はneibを再生する必要がある。
    !
    do ii=1, lost%size
       i = lost%set(ii)
       do jj=1, q6s%nei3(i)%nnear
          j = q6s%nei3(i)%near(jj)
          call intersection( q6s%neis(i), q6s%neis(j), q6s%neib( jj, i) )
       end do
    end do
  end subroutine sq6system_Integrate

  subroutine sq6system_Reject( q6s, q6d, target, nmol )
    type( sq6system ) :: q6s
    type( sq6diff )   :: q6d
    integer, intent(IN) :: target, nmol

    integer :: i,j, k,m, ii,jj,kk

    !
    !qlmsumの復旧
    !
    do jj=1, q6d%newnei3(1)%nnear
       j = q6d%newnei3(1)%near(jj)
       do kk=1, q6d%neib( jj )%size
          k = q6d%neib( jj )%set(kk)
          q6s%nbond(k) = q6s%nbond(k) - 2
          do m=-6,6
             q6s%qlmsum( m, k ) = q6s%qlmsum( m, k ) - q6d%qlm( m, jj ) * 2
          enddo
       enddo
    enddo
    do ii=1, q6s%nei3(target)%nnear
       i = q6s%nei3(target)%near(ii)
       !
       !target-i結合の近傍の分子について
       !
       do kk=1, q6s%neib( ii, target )%size
          k = q6s%neib( ii, target )%set(kk)
          q6s%nbond(k) = q6s%nbond(k) + 2
          do m=-6,6
             q6s%qlmsum(m,k) = q6s%qlmsum(m,k) + q6s%qlm( m, ii, target ) * 2
          enddo
       enddo
    enddo
    do ii=1, q6d%plus%size
       i = q6d%plus%set(ii)
       do jj=1, q6s%nei3(i)%nnear
          j = q6s%nei3(i)%near(jj)
          if ( set_exist( q6d%neis, j ) ) then
             q6s%nbond(target) = q6s%nbond(target) - 2
             do m=-6,6
                q6s%qlmsum( m, target ) = q6s%qlmsum( m, target ) - q6s%qlm( m, jj, i ) * 2
             enddo
          endif
       enddo
    enddo
    do ii=1, q6d%minus%size
       i = q6d%minus%set(ii)
       do jj=1, q6s%nei3(i)%nnear
          j = q6s%nei3(i)%near(jj)
          if ( set_exist( q6s%neis(target), j ) ) then
             q6s%nbond(target) = q6s%nbond(target) + 2
             do m=-6,6
                q6s%qlmsum( m, target ) = q6s%qlmsum( m, target ) + q6s%qlm( m, jj, i ) * 2
             enddo
          endif
       enddo
    enddo
#ifdef DEBUG20050309
    call sq6system_write( q6s, nmol, STDOUT )
    write(STDOUT,*) "######################################"
#endif
  end subroutine sq6system_Reject

  subroutine sq6system_EndAStep( q6s )
    type( sq6system ) :: q6s
  end subroutine sq6system_EndAStep

  subroutine sQ6System_Write( q6s, nmol, file )
    use neighbors_module
    type( sq6system ) :: q6s
    integer :: file
    integer :: nmol

    integer :: i,j,k,ii,jj,kk
    write(file,*) "nei6========================="
    do i=1, nmol
       write(file,*) i
       call sneighbor1_write( q6s%nei6(i), file )
    enddo
    write(file,*) "nei5========================="
    do i=1, nmol
       write(file,*) i
       call sneighbor1_write( q6s%nei5(i), file )
    enddo
    write(file,*) "nei3========================="
    do i=1, nmol
       write(file,*) i
       call sneighbor1_write( q6s%nei3(i), file )
    enddo
    write(file,*) "q6()=========================", q6s%sumq6
    do i=1, nmol
       write(file,*) i, q6s%q6(i)
    enddo
    write(file,*) "qlmsum()========================="
    do i=1, nmol
       write(file,*) i, q6s%qlmsum( -6, i )
    enddo
    write(file,*) "qlm()========================="
    do i=1, nmol
       do j=1,q6s%nei3(i)%nnear
          write(file,*) i, j, q6s%qlm( -6, j, i )
       enddo
    enddo
    write(file,*) "neis()========================="
    do i=1, nmol
       write(file,*) i
       call sset_write( q6s%neis(i), file )
    enddo
    write(file,*) "neib()========================="
    do i=1, nmol
       write(file,*) i
       do j=1,q6s%nei3(i)%nnear
          write(file,*) i,j
          call sset_write( q6s%neib(j,i), file )
       enddo
    enddo
    write(file,*) "nbond()========================="
    do i=1, nmol
       write(file,*) i,q6s%nbond(i)
    enddo
    !
    !ある分子のq6に寄与している結合のリストがほしい。
    !
    do i=1, nmol
       do jj=1,q6s%nei3(i)%nnear
          j = q6s%nei3(i)%near(jj)
          do kk=1, q6s%neib(jj,i)%size
             k = q6s%neib(jj,i)%set(kk)
             if ( k == 5 ) then
                write(file,*) i,j,k
                !if ( ( i==4 .and. j==34 ) .or. ( i==34 .and. j==4 ) ) then
                !   write(file,*) "4,17",distance(q6s,4,17)
                !   write(file,*) "34,17",distance(q6s,34,17)
                !   write(file,*) "4,34",distance(q6s,4,34)
                !endif
             endif
          enddo
       enddo
    enddo
  end subroutine sQ6System_Write

  function distance( q6s, i, j )
    type( sq6system ) :: q6s
    integer :: i,j,k,kk
    real(kind=8) :: distance
    distance = 0d0
    do kk=1, q6s%nei6(i)%nnear
       k = q6s%nei6(i)%near(kk)
       if ( k == j ) then
          distance = sqrt(q6s%nei6(i)%dd(kk))
       endif
    enddo
  end function distance

  subroutine sQ6System_Compare( q1, q2, nmol, file )
    use neighbors_module
    type( sq6system ) :: q1,q2
    integer :: file, nmol

    integer :: i,j,k,ii,jj,kk
    write(file,*) "nei6++++++++++++++++++++++++="
    !
    !nei6は、そもそも緩衝帯として設けてあるので、少々不一致でも構わな
    !い。nei6からnei5やneisが正しく生成されればよい。
    !
    do i=1, nmol
       write(file,*) i
       call sneighbor1_compare( q1%nei6(i), q2%nei6(i), file )
    enddo
    write(file,*) "nei5++++++++++++++++++++++++="
    do i=1, nmol
       write(file,*) i
       call sneighbor1_compare( q1%nei5(i), q2%nei5(i), file )
    enddo
    write(file,*) "nei3++++++++++++++++++++++++="
    do i=1, nmol
       write(file,*) i
       call sneighbor1_compare( q1%nei3(i), q2%nei3(i), file )
    enddo
    write(file,*) "q6()++++++++++++++++++++++++="
    if ( 1d-12 < real_compare( q1%sumq6, q2%sumq6 ) ) then
       write(file,*) real_compare( q1%sumq6, q2%sumq6 ), "differs(sumq6)."
    endif
    do i=1, nmol
       if ( 1d-12 < real_compare( q1%q6(i), q2%q6(i) ) ) then
          write(file,*) i, real_compare( q1%q6(i), q2%q6(i) ), "differs(q6(i))."
       endif
    enddo
    write(file,*) "qlmsum()++++++++++++++++++++++++="
    do i=1, nmol
       if ( 1d-12 < complex_compare( q1%qlmsum( -6, i), q2%qlmsum( -6, i) ) ) then
          write(file,*) i, complex_compare( q1%qlmsum( -6, i), q2%qlmsum( -6, i) ), "differs(qlmsum( -6, i))."
       endif
    enddo
    write(file,*) "qlm()++++++++++++++++++++++++="
    do i=1, nmol
       k = min(q1%nei3(i)%nnear,q2%nei3(i)%nnear)
       do jj=1,k
          if ( 1d-12 < complex_compare( q1%qlm( -6, jj, i), q2%qlm( -6, jj, i) ) ) then
             write(file,*) i, jj, q1%nei3(i)%near(jj), q2%nei3(i)%near(jj), complex_compare( q1%qlm( -6, jj, i), q2%qlm( -6, jj, i) ), "differs(qlm( -6,jj,i)).",q1%qlm( -6, jj, i), q2%qlm( -6, jj, i)
          endif
       enddo
    enddo
    !      write(file,*) i, j, q6s%qlm( -6, j, i )
    !   enddo
    !enddo
    write(file,*) "neis()++++++++++++++++++++++++="
    do i=1, nmol
       write(file,*) i
       call sset_compare( q1%neis(i), q2%neis(i), file )
    enddo
    write(file,*) "neib()++++++++++++++++++++++++="
    do i=1, nmol
       write(file,*) i
       write(file,*) q1%nei3(i)%nnear, q2%nei3(i)%nnear, "(nei3(i)%nnear)."
       if ( q1%nei3(i)%nnear /= q2%nei3(i)%nnear ) then
          write(file,*) "differs."
       endif
       k = min(q1%nei3(i)%nnear,q2%nei3(i)%nnear)
       do j=1,k
          write(file,*) i,j,q1%nei3(i)%near(j),q2%nei3(i)%near(j)
          call sset_compare( q1%neib(j,i), q2%neib(j,i), file )
       enddo
    enddo
    write(file,*) "nbond()++++++++++++++++++++++++="
    do i=1, nmol
       if ( q1%nbond(i) /= q2%nbond(i) ) then
          write(file,*) q1%nbond(i), q2%nbond(i), " differs(nbond(i))."
       endif
    enddo
  end subroutine sQ6System_Compare

  function real_compare( x,y )
    real(kind=8) :: x,y,real_compare
    real_compare = abs((x-y)/y)
  end function real_compare

  function complex_compare( x,y )
    complex(kind=8) :: x,y
    real(kind=8) :: complex_compare
    complex_compare = abs((x-y)/y)
  end function complex_compare

end module q6local2
