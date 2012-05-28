! -*- f90 -*-
!
!純物質剛体系をquenchするプログラムのサンプル。
!
!他のモジュールと異なり、内部変数をもつので、多重呼びだしなどはできな
!い。(汎用でない)。剛体混合系以外に使うなら、同様のモジュールをカスタ
!ム製作する必要がある。
!
module rigidquench_module
  use common_module
  use rigid_module
  use mol_module
  use site_module
  use grid_module
  use box_module
  use interaction_module
  use interaction2_module
  use cutoff_module
  use standard_interaction_module

  implicit none
  integer, parameter:: RQ0=0,RQ_EULER=1, RQ_QUAT=2, RQ_EULER0=3
  
  type(sRigid),  pointer, private :: rigid(:)
  type(sMol),    pointer, private :: mol(:)
  type(sStdInt), pointer, private :: si(:)
  type(sSite),   pointer, private :: site
  type(sGrid),   pointer, private :: grid
  type(sBox),    pointer, private :: box
  type(sCutoff), pointer, private :: cutoff

  type(sInteraction)    , private :: inter(MAXCOMBI)
  type(sAddress)        , private :: addr(MAXCOMPO)
  type(sNeighborCell)   , private :: homo, hetero
  integer               , private :: ncompo, ncombi
  integer               , private :: combi(MAXCOMBI)
  integer               , private :: combj(MAXCOMBI)
  integer               , private :: mode
     
contains

  subroutine rigidquench_initialize(ncompo0, rigid0, mol0,&
       site0, grid0, cutoff0, box0, mode0, si0, ivscale)
    type(sRigid), target :: rigid0(1:)
    type(sMol),   target :: mol0(1:)
    type(sSite),  target :: site0
    type(sGrid),  target :: grid0
    type(sBox),   target :: box0
    type(sCutoff),target :: cutoff0
    integer,intent(in)   :: ncompo0
    type(sStdInt),target,optional :: si0(1:)
    real(kind=8),intent(IN), optional :: ivscale(MAXCOMPO, MAXCOMPO)

    real(kind=8)         :: rc
    integer              :: i,j,k
    integer,intent(in)   :: mode0
    !とりあえず、モジュール内で使用する変数を受け取って、ポインタを保
    !管しておく？
    rigid  => rigid0
    mol    => mol0
    site   => site0
    grid   => grid0
    box    => box0
    cutoff => cutoff0
    ncompo =  ncompo0
    mode   =  mode0
    if ( present( si0 ) ) then
       si  => si0
    else
       nullify( si )
    endif
    if(box%mode.ne.BOX_NONE)then
       do i=1, ncompo
          call Rigid_Relocate(rigid(i),mol(i),box)
       enddo
    endif
    if(box%mode.ne.BOX_NONE.and.cutoff%mode.eq.CUTOFF_NONE)then
       cutoff%mode=CUTOFF_POLY
       rc=dmin1(box%size%vec(1),box%size%vec(2),box%size%vec(3))*0.5d0
       call cutoff_initialize(cutoff,rc-2d0,rc)
    endif
    if(box%mode.eq.BOX_NONE)write(STDERR,*) "CLUSTER SYSTEM"
    if(cutoff%mode.eq.CUTOFF_NONE)then
       write(STDERR,*) "INTERACTION RANGE IS NOT LIMITED"
    endif
    !
    !相互作用対リストを初期化
    !
    !
    !同種分子間の相互作用(分子数をnとするとn(n+1)/2対)
    !
    do j=1,ncompo
       !
       !第2パラメータを省略した場合は自己相互作用とみなす。
       !
       call interaction_initialize(inter(j),mol(j)%nmol,-1)
       combi( j ) = j
       combj( j ) = j
    enddo
    !
    !異種分子間の相互作用(分子数をn,mとするとnm対)
    !
    ncombi = ncompo
    do i=1,ncompo
       do j=i+1,ncompo
          ncombi = ncombi + 1
          combi( ncombi ) = i
          combj( ncombi ) = j
          call interaction_initialize(inter(ncombi),&
               mol(i)%nmol, mol(j)%nmol)
       enddo
    enddo
    if ( present( ivscale ) ) then
       do k=1,ncombi
          i = combi( k )
          j = combj( k )
          inter(k)%scale = ivscale(i,j)
       enddo
    endif
    !
    !グリッド分割する場合
    !
    if( active(grid) )then
       !
       !各分子が所属するセルの表を格納するaddrを準備
       !
       do i=1,ncompo
          call Address_Initialize(addr(i),mol(i)%nmol)
       enddo
       call Grid_Update( grid, box, 0d0 )
    else
       !
       !全ての分子はほかの全ての分子と相互作用しうる。
       !
       do i=1,ncombi
          call interaction_alltoall(inter(i))
       enddo
    endif
  end subroutine rigidquench_initialize

  function rigidquench_func(p, noforce)
    use vector_module
    use tip4p_module
    use tip5p_module
    use st2_module
    use nvde_module
    use oplsmeoh_module
    use oplsmeoh_tip4p_module
    use matrix_module
    use interaction2_module
    implicit none
    real(kind=8), intent(in) :: p(*)
    logical,      intent(in) :: noforce
    real(kind=8) :: rigidquench_func
    real(kind=8) :: th,ph,ps
    real(kind=8) :: qa,qb,qc,qd,r
    real(kind=8) :: sina,sinb,sinc
    real(kind=8) :: cosa,cosb,cosc
    real(kind=8) :: epsum
    type(vector3):: com(MAXMOL, MAXCOMPO)
    integer      :: i, j, k, offset, mode
    real(kind=8), allocatable :: ep(:)
    type(sTensor) :: vr
    mode = rigidquench_mode()
    !
    !シリアライズされた配列から、剛体分子の位置と配向を抽出する。
    !
    offset = 0
    if ( mode .eq. RQ_EULER ) then
       do j=1,ncompo
          DO I=1,mol(j)%nmol
             rigid(j)%mol(i)%com%vec(1) = p(i+mol(j)%nmol*0+offset)
             rigid(j)%mol(i)%com%vec(2) = p(i+mol(j)%nmol*1+offset)
             rigid(j)%mol(i)%com%vec(3) = p(i+mol(j)%nmol*2+offset)
             TH                         = p(i+mol(j)%nmol*3+offset)
             PH                         = p(i+mol(j)%nmol*4+offset)
             PS                         = p(i+mol(j)%nmol*5+offset)
             SINA=SIN(TH)
             SINB=SIN(PH)
             SINC=SIN(PS)
             COSA=COS(TH)
             COSB=COS(PH)
             COSC=COS(PS)
             rigid(j)%mol(i)%t11        = COSC*COSB-COSA*SINB*SINC
             rigid(j)%mol(i)%t12        =-(SINC*COSB+COSA*SINB*COSC)
             rigid(j)%mol(i)%t13        = SINA*SINB
             rigid(j)%mol(i)%t21        = COSC*SINB+COSA*COSB*SINC
             rigid(j)%mol(i)%t22        =-SINC*SINB+COSA*COSB*COSC
             rigid(j)%mol(i)%t23        =-SINA*COSB
             rigid(j)%mol(i)%t31        = SINA*SINC
             rigid(j)%mol(i)%t32        = SINA*COSC
             rigid(j)%mol(i)%t33        = cosa
          enddo
          offset = offset + Mol_DoF( mol(j) )
       enddo
    else if ( mode .eq. RQ_QUAT ) then
       do j=1,ncompo
          DO I=1,mol(j)%nmol
             !write(*,*) "si10",i,j,offset
             rigid(j)%mol(i)%com%vec(1) = p(i+mol(j)%nmol*0+offset)
             rigid(j)%mol(i)%com%vec(2) = p(i+mol(j)%nmol*1+offset)
             rigid(j)%mol(i)%com%vec(3) = p(i+mol(j)%nmol*2+offset)
             qa                         = p(i+mol(j)%nmol*3+offset)
             qb                         = p(i+mol(j)%nmol*4+offset)
             qc                         = p(i+mol(j)%nmol*5+offset)
             qd                         = p(i+mol(j)%nmol*6+offset)
             !
             !本来なら、四元数の変化分を、四元数速度の定数倍にする計算は
             !quaternion_mulを使わなければいけないが、変化分が小さいもの
             !と考え、ここで規格化する。(どこかでp()にフィードバックしな
             !ければいけない。)
             !
             r = 1d0 / dsqrt(qa**2 + qb**2 + qc **2 + qd**2)
             qa = qa * r
             qb = qb * r
             qc = qc * r
             qd = qd * r
             !write(*,*) "si7",r,ncompo
             !
             rigid( j )%mol(i)%t11      = (qa*qa+qb*qb-(qc*qc+qd*qd))
             rigid( j )%mol(i)%t12      = -2.0d0*(qa*qd+qb*qc)
             rigid( j )%mol(i)%t13      = 2.0d0*(qb*qd-qa*qc)
             rigid( j )%mol(i)%t21      = 2.0d0*(qa*qd-qb*qc)
             rigid( j )%mol(i)%t22      = qa*qa+qc*qc-(qb*qb+qd*qd)
             rigid( j )%mol(i)%t23      = -2.0d0*(qa*qb+qc*qd)
             rigid( j )%mol(i)%t31      = 2.0d0*(qa*qc+qb*qd)
             rigid( j )%mol(i)%t32      = 2.0d0*(qa*qb-qc*qd)
             rigid( j )%mol(i)%t33      = qa*qa+qd*qd-(qb*qb+qc*qc)
          enddo
          offset = offset + mol(j)%nmol * 7
       enddo
    endif
    !
    !サイトに加わる力をリセットする。
    !
    call site_resetforce(site)
    !
    !剛体分子の座標から、サイトの座標を算出する。
    !rigid_setsiteposition()を使うと、quaternionから回転行列を再計算
    !してしまうので、代わりにrigid_setsiteposition0を使用している。
    !
    do j=1,ncompo
       call rigid_setsiteposition0(rigid(j),mol(j),site)
    enddo
    !
    !重心位置を一時配列にコピーしておく。
    !
    do j=1, ncompo
       com( 1:mol(j)%nmol, j ) = rigid( j )%mol(1:mol(j)%nmol)%com
    enddo
    if( active( grid ) )then
       !write(STDERR,*) (com(i)%vec(1),i=1,10),mol%nmol
       !
       !相互作用対リストを初期化する。
       !
       do i=1,ncombi
          call Interaction_Reset( inter( i ) )
       enddo
       !
       !分子をセルに割りあてる。
       !
       do i=1,ncompo
          call Address_Assign(addr(i),box,grid,mol(i)%nmol,com(1,i))
       enddo
       !
       !隣接するセルの表から、相互作用対リストを生成する
       !
       !
       !同種分子間相互作用対リストの作成
       !
       do i=1,ncompo
          call Grid_NeighborList(grid,grid%homo,inter(i),addr(i),addr(i))
       enddo
       !
       !異種分子間相互作用対リストの作成
       !
       do i=ncompo+1,ncombi
          call Grid_NeighborList(grid, grid%hetero, inter(i),&
               addr(combi(i)),addr(combj(i)))
       enddo
       !
       !カットオフ外の対を除いて、相互作用対リストを圧縮する。
       !
       do k=1,ncombi
          i=combi(k)
          j=combj(k)
          call Interaction_Compress2(inter(k),box%mode.ne.BOX_NONE,&
               box ,cutoff%out**2, com(1,i),com(1,j),.false.)
       enddo
    else
       do k=1,ncombi
          i = combi( k )
          j = combj( k )
          !
          !カットオフ外の対を除いて、相互作用対リストを圧縮する。
          !相互作用対リストを圧縮して、あとの力計算処理を高速化する
          !(ベクトル機では有効)
          !
          if ( cutoff%mode.ne.CUTOFF_NONE) then
             call Interaction_Compress(inter(k),com(1,i),com(1,j),&
                  .false. ,CutOff_OuterLimit(cutoff),box)
          else
             call Interaction_Compress(inter(k),com(1,i),com(1,j),&
                  .false. ,-1d0,box)
          endif
       enddo
    endif
    !
    !力を計算。ここではエネルギーしか使用しないが、dfuncで力も使う。
    !
    epsum = 0d0
    !vrsum = 0d0
    do k=1,ncombi
       i = combi( k )
       j = combj( k )
       !
       !Smooth cutoffによる減衰係数を事前計算
       !
       call CutOff_Dimmer2(cutoff, inter(k), com(1,i), com(1,j) )
       !
       !妙案が思いつかないので、その場で分子種を見ながらどの相互作用
       !関数を使うかを判断する。将来別の方法を採用する可能性あり。
       !
       allocate( ep(inter(k)%npair) )
       if ( associated( si ) ) then
          if ( noforce ) then
             call PairPotential(inter(k),mol(i),mol(j),site,ep,.true.,si(i),si(j))
          else
             call force(inter(k),mol(i),mol(j),site,ep,vr%mat33,.true.,si(i),si(j))
          endif
       else
          if ( noforce ) then
             call PairPotential(inter(k),mol(i),mol(j),site,ep,.true.)
          else
             call force(inter(k),mol(i),mol(j),site,ep,vr%mat33,.true.)
          endif
       endif
       do i=1,inter(k)%npair
          epsum = epsum + ep(i)
       enddo
       deallocate( ep )
       !vrsum = vrsum + vr(k)
    enddo
    !write(STDERR,*) ep(1:MAXCOMBI)
    rigidquench_func=epsum
    !write(STDERR,*) epsum
  end function rigidquench_func
  !
  !Quaternionを使用するバージョン(試作)
  !
  real(kind=8) function rigidquench_funcq(p)
    use vector_module
    use tip4p_module
    use tip5p_module
    use oplsmeoh_module
    use oplsmeoh_tip4p_module
    implicit none
    real(kind=8), intent(in) :: p(*)
    real(kind=8) :: qa,qb,qc,qd,r
    real(kind=8) :: epsum
    type(vector3):: com(MAXMOL, MAXCOMPO)
    integer      :: i, j, k, offset
    real(kind=8),allocatable :: ep(:)
    real(kind=8) :: vrmat(3,3)
    !
    !シリアライズされた配列から、剛体分子の位置と配向を抽出する。
    !
    offset = 0
    !write(*,*) "si11     ",ncompo,mol(1)%nmol
    do j=1,ncompo
       DO I=1,mol(j)%nmol
          !write(*,*) "si10",i,j,offset
          rigid(j)%mol(i)%com%vec(1) = p(i+mol(j)%nmol*0+offset)
          rigid(j)%mol(i)%com%vec(2) = p(i+mol(j)%nmol*1+offset)
          rigid(j)%mol(i)%com%vec(3) = p(i+mol(j)%nmol*2+offset)
          qa                         = p(i+mol(j)%nmol*3+offset)
          qb                         = p(i+mol(j)%nmol*4+offset)
          qc                         = p(i+mol(j)%nmol*5+offset)
          qd                         = p(i+mol(j)%nmol*6+offset)
          !
          !本来なら、四元数の変化分を、四元数速度の定数倍にする計算は
          !quaternion_mulを使わなければいけないが、変化分が小さいもの
          !と考え、ここで規格化する。(どこかでp()にフィードバックしな
          !ければいけない。)
          !
          r = 1d0 / dsqrt(qa**2 + qb**2 + qc **2 + qd**2)
          qa = qa * r
          qb = qb * r
          qc = qc * r
          qd = qd * r
          !write(*,*) "si7",r,ncompo
         !
          rigid( j )%mol(i)%t11      = (qa*qa+qb*qb-(qc*qc+qd*qd))
          rigid( j )%mol(i)%t12      = -2.0d0*(qa*qd+qb*qc)
          rigid( j )%mol(i)%t13      = 2.0d0*(qb*qd-qa*qc)
          rigid( j )%mol(i)%t21      = 2.0d0*(qa*qd-qb*qc)
          rigid( j )%mol(i)%t22      = qa*qa+qc*qc-(qb*qb+qd*qd)
          rigid( j )%mol(i)%t23      = -2.0d0*(qa*qb+qc*qd)
          rigid( j )%mol(i)%t31      = 2.0d0*(qa*qc+qb*qd)
          rigid( j )%mol(i)%t32      = 2.0d0*(qa*qb-qc*qd)
          rigid( j )%mol(i)%t33      = qa*qa+qd*qd-(qb*qb+qc*qc)
       enddo
       offset = offset + mol(j)%nmol * 7
    enddo
    !
    !サイトに加わる力をリセットする。
    !
    call site_resetforce(site)
    !
    !剛体分子の座標から、サイトの座標を算出する。
    !rigid_setsiteposition()を使うと、quaternionから回転行列を再計算
    !してしまうので、代わりにrigid_setsiteposition0を使用している。
    !
    do j=1,ncompo
       call rigid_setsiteposition0(rigid(j),mol(j),site)
    enddo
    !
    !重心位置を一時配列にコピーしておく。
    !
    do j=1, ncompo
       com( 1:mol(j)%nmol, j ) = rigid( j )%mol(1:mol(j)%nmol)%com
    enddo
    if( active( grid ) )then
       !write(STDERR,*) (com(i)%vec(1),i=1,10),mol%nmol
       !
       !相互作用対リストを初期化する。
       !
       do i=1,ncombi
          call Interaction_Reset( inter( i ) )
       enddo
       !
       !分子をセルに割りあてる。
       !
       do i=1,ncompo
          call Address_Assign(addr(i),box,grid,mol(i)%nmol,com(1,i))
       enddo
       !
       !隣接するセルの表から、相互作用対リストを生成する
       !
       !
       !同種分子間相互作用対リストの作成
       !
       do i=1,ncompo
          call Grid_NeighborList(grid,homo,inter(i),addr(i),addr(i))
       enddo
       !
       !異種分子間相互作用対リストの作成
       !
       do i=ncompo+1,ncombi
          call Grid_NeighborList(grid, hetero, inter(i),&
               addr(combi(i)),addr(combj(i)))
       enddo
       !
       !カットオフ外の対を除いて、相互作用対リストを圧縮する。
       !
       do k=1,ncombi
          i=combi(k)
          j=combj(k)
          call Interaction_Compress2(inter(k),box%mode.ne.BOX_NONE,&
               box ,cutoff%out**2, com(1,i),com(1,j),.false.)
       enddo
    else
       do k=1,ncombi
          i = combi( k )
          j = combj( k )
          !
          !カットオフ外の対を除いて、相互作用対リストを圧縮する。
          !相互作用対リストを圧縮して、あとの力計算処理を高速化する
          !(ベクトル機では有効)
          !
          if ( cutoff%mode.ne.CUTOFF_NONE) then
             call Interaction_Compress(inter(k),com(1,i),com(1,j),&
                  .false. ,CutOff_OuterLimit(cutoff),box)
          else
             call Interaction_Compress(inter(k),com(1,i),com(1,j),&
                  .false. ,-1d0,box)
          endif
       enddo
    endif
    !
    !力を計算。ここではエネルギーしか使用しないが、dfuncで力も使う。
    !
    epsum = 0d0
    do k=1,ncombi
       i = combi( k )
       j = combj( k )
       !
       !Smooth cutoffによる減衰係数を事前計算
       !
       call CutOff_Dimmer2(cutoff, inter(k), com(1,i), com(1,j) )
       !
       !妙案が思いつかないので、その場で分子種を見ながらどの相互作用
       !関数を使うかを判断する。将来別の方法を採用する可能性あり。
       !
       allocate( ep(inter(k)%npair) )
       if ( associated( si ) ) then
          call force(inter(k),mol(i),mol(j),site,ep,vrmat, .true., si(i),si(j))
       else
          if ( mol(i)%id(1:8) .eq. "TIP4P   " )then
             if ( mol(j)%id(1:8) .eq. "TIP4P   " )then
                call interaction_force_tip4p(&
                     inter(k),mol(i),mol(j),site,ep,vrmat)
             else if ( mol(j)%id(1:8) .eq. "OPLSMEOH" )then
                call interaction_force_tip4p_oplsmeoh(&
                     inter(k),mol(i),mol(j),site,ep,vrmat)
             else
                write(STDERR,*) "Unknown interaction", mol(i)%id(1:8)&
                     & ,":",mol(j)%id(1:8)
                stop
             endif
          else if ( mol(i)%id(1:8) .eq. "OPLSMEOH" )then
             if ( mol(j)%id(1:8) .eq. "TIP4P   " )then
                call interaction_force_tip4p_oplsmeoh(&
                     inter(k),mol(j),mol(i),site,ep,vrmat)
             else if ( mol(j)%id(1:8) .eq. "OPLSMEOH" )then
                call interaction_force_oplsmeoh(&
                     inter(k),mol(i),mol(j),site,ep,vrmat)
             else
                write(STDERR,*) "Unknown interaction", mol(i)%id(1:8)&
                     & ,":",mol(j)%id(1:8)
                stop
             endif
          else
             write(STDERR,*) "Unknown interaction", mol(i)%id(1:8)&
                  & ,":"&
                  & ,mol(j)%id(1:8)
             stop
          endif
       endif
       do i=1,inter(k)%npair
          epsum = epsum + ep(i)
       enddo
       deallocate( ep )
    enddo
    !write(STDERR,*) ep(1:MAXCOMBI)
    rigidquench_funcq = epsum
    !write(STDERR,*) epsum
  end function rigidquench_funcq

  subroutine do_nothing(p)
    real(kind=8), intent(in) :: p(*)
  end subroutine do_nothing

  subroutine q_normalize(p)
    implicit none
    real(kind=8), intent(inout) :: p(*)
    real(kind=8) :: qa,qb,qc,qd,r
    integer      :: i, j, offset
    !
    !シリアライズされた配列から、剛体分子の位置と配向を抽出する。
    !
    offset = 0
    do j=1,ncompo
       DO I=1,mol(j)%nmol
          qa                         = p(i+mol(j)%nmol*3+offset)
          qb                         = p(i+mol(j)%nmol*4+offset)
          qc                         = p(i+mol(j)%nmol*5+offset)
          qd                         = p(i+mol(j)%nmol*6+offset)
          r = 1d0 / dsqrt(qa**2 + qb**2 + qc **2 + qd**2)
          p(i+mol(j)%nmol*3+offset) = qa * r
          p(i+mol(j)%nmol*4+offset) = qb * r
          p(i+mol(j)%nmol*5+offset) = qc * r
          p(i+mol(j)%nmol*6+offset) = qd * r
       enddo
       offset = offset + mol( j )%nmol * 7
    enddo
  end subroutine q_normalize
  !
  !こちらは変位が力に直接比例する。(質量で割っていない)
  !
  subroutine rigidquench_dfunc(p,xi)
    real(kind=8) :: p(*),xi(*)
    real(kind=8) :: th,ph,ps
    real(kind=8) :: sina,sinb,sinc
    real(kind=8) :: cosa,cosb,cosc
    integer      :: i, j, offset
    real(kind=8) :: tx,ty,tz
    !みょうな処理順序だが、実際にはここでcomxなどに代入しなくても直前の
    !func()で設定されているはず。
    offset = 0
    do j=1,ncompo
       if ( mol(j)%isFixed ) then
          !
          !座標が固定されている分子には力を及ぼさない。
          !
          xi( i + mol( j )%nmol * 0 + offset ) = 0d0
          xi( i + mol( j )%nmol * 1 + offset ) = 0d0
          xi( i + mol( j )%nmol * 2 + offset ) = 0d0
          xi( i + mol( j )%nmol * 3 + offset ) = 0d0
          xi( i + mol( j )%nmol * 4 + offset ) = 0d0
          xi( i + mol( j )%nmol * 5 + offset ) = 0d0
       else
          !
          !各剛体分子に加わる並進/回転の力を算出
          !
          call rigid_collectforce(rigid(j),mol(j),site)
          DO I=1,mol(j)%nmol
             TH=p( i + mol( j )%nmol * 3 + offset )
             PH=p( i + mol( j )%nmol * 4 + offset )
             PS=p( i + mol( j )%nmol * 5 + offset )
             SINA=SIN(TH)
             SINB=SIN(PH)
             SINC=SIN(PS)
             COSA=COS(TH)
             COSB=COS(PH)
             COSC=COS(PS)
             xi( i + mol( j )%nmol * 0 + offset ) = &
                  -rigid(j)%mol(i)%forcex
             xi( i + mol( j )%nmol * 1 + offset ) = &
                  -rigid(j)%mol(i)%forcey
             xi( i + mol( j )%nmol * 2 + offset ) = &
                  -rigid(j)%mol(i)%forcez
             tx = rigid( j )%mol(i)%torquex
             ty = rigid( j )%mol(i)%torquey
             tz = rigid( j )%mol(i)%torquez
             xi( i + mol( j )%nmol * 3 + offset ) = (-tx*cosc + ty*sinc)
             xi( i + mol( j )%nmol * 4 + offset ) = &
                  (-tx*sina*sinc - ty*sina*cosc  - tz*cosa)
             xi( i + mol( j )%nmol * 5 + offset ) = (-tz)
          enddo
       endif
       offset = offset + Mol_DoF( mol( j ) )
    enddo
  end subroutine rigidquench_dfunc
  !
  !こちらは変位が力/質量に比例するようにしたい。(開発中)
  !
  subroutine rigidquench_dfunc2(p,xi)
    real(kind=8) :: p(*),xi(*)
    real(kind=8) :: th,  ph,  ps
    real(kind=8) :: dth, dph, dps
    real(kind=8) :: sina,sinb,sinc
    real(kind=8) :: cosa,cosb,cosc
    integer      :: i, j, offset
    real(kind=8) :: tx, ty, tz
#undef DEBUG
#ifdef DEBUG
    real(kind=8) :: pp, qq, rr, ss, tt, uu
    real(kind=8) :: dp, dq, dr, ds, dt, du, da,db,dc,dd
#endif
    !みょうな処理順序だが、実際にはここでcomxなどに代入しなくても直前の
    !func()で設定されているはず。
    offset = 0
    do j=1,ncompo
       if ( mol(j)%isFixed ) then
          !
          !座標が固定されている分子には力を及ぼさない。
          !
          xi( i + mol( j )%nmol * 0 + offset ) = 0d0
          xi( i + mol( j )%nmol * 1 + offset ) = 0d0
          xi( i + mol( j )%nmol * 2 + offset ) = 0d0
          xi( i + mol( j )%nmol * 3 + offset ) = 0d0
          xi( i + mol( j )%nmol * 4 + offset ) = 0d0
          xi( i + mol( j )%nmol * 5 + offset ) = 0d0
       else
          !
          !各剛体分子に加わる並進/回転の力を算出
          !
          call rigid_collectforce(rigid(j),mol(j),site)
          DO I=1,mol(j)%nmol
             TH=p( i + mol( j )%nmol * 3 + offset )
             PH=p( i + mol( j )%nmol * 4 + offset )
             PS=p( i + mol( j )%nmol * 5 + offset )
             SINA=SIN(TH)
             SINB=SIN(PH)
             SINC=SIN(PS)
             COSA=COS(TH)
             COSB=COS(PH)
             COSC=COS(PS)
             xi( i + mol( j )%nmol * 0 + offset ) = &
                  -rigid(j)%mol(i)%forcex * rigid(j)%massi
             xi( i + mol( j )%nmol * 1 + offset ) = &
                  -rigid(j)%mol(i)%forcey * rigid(j)%massi
             xi( i + mol( j )%nmol * 2 + offset ) = &
                  -rigid(j)%mol(i)%forcez * rigid(j)%massi
             !
             !剛体の回転の運動方程式はGoldsteinのP.268に記載されている。
             !静止状態でトルクを与えると、各軸回りの角加速度は、各軸まわ
             !りのトルクに比例する。
             !
             tx = rigid( j )%mol(i)%torquex * rigid(j)%ixxi
             ty = rigid( j )%mol(i)%torquey * rigid(j)%iyyi
             tz = rigid( j )%mol(i)%torquez * rigid(j)%izzi
             dth = (-tx*cosc + ty*sinc)
             dph = (-tx*sina*sinc - ty*sina*cosc - tz*cosa)
             dps = -tz
             xi( i + mol( j )%nmol * 3 + offset ) = dth
             xi( i + mol( j )%nmol * 4 + offset ) = dph
             xi( i + mol( j )%nmol * 5 + offset ) = dps
#ifdef DEBUG
             pp = dcos( th / 2d0 )
             qq = dsin( th / 2d0 )
             rr = dcos( ( ph + ps ) / 2d0 )
             ss = dcos( ( ps - ph ) / 2d0 )
             tt = dsin( ( ps - ph ) / 2d0 )
             uu = dsin( ( ph + ps ) / 2d0 )
             dp = -dth / 2d0 * qq
             dq =  dth / 2d0 * pp
             dr = ( dph + dps ) / 2d0 * (-1) * uu
             ds = ( dps - dph ) / 2d0 * (-1) * tt
             dt = ( dps - dph ) / 2d0 * ss
             du = ( dph + dps ) / 2d0 * rr
             da = dp * rr + pp * dr
             db = dq * ss + qq * ds
             dc = dq * tt + qq * dt
             dd = dp * uu + pp * du
             write(STDOUT,*) da,db,dc,dd, tx, ty, tz
             stop
#endif
          enddo
       endif
       offset = offset + Mol_DoF( mol( j ) )
    enddo
  end subroutine rigidquench_dfunc2
  !
  !Quaternionのままで共役勾配法を行う。linminなどの処理で、Quaternion
  !の定数倍計算を行う必要がでてくる。また、配列pの大きさに注意する。
  !
  subroutine rigidquench_dfuncq(p,xi)
    real(kind=8) :: p(*),xi(*)
    real(kind=8) :: qa,qb,qc,qd,da,db,dc,dd
    integer      :: i, j, offset
    real(kind=8) :: tx, ty, tz
    !みょうな処理順序だが、実際にはここでcomxなどに代入しなくても直前の
    !func()で設定されているはず。
    offset = 0
    do j=1,ncompo
       if ( mol(j)%isFixed ) then
          !
          !座標が固定されている分子には力を及ぼさない。
          !
          xi( i + mol( j )%nmol * 0 + offset ) = 0d0
          xi( i + mol( j )%nmol * 1 + offset ) = 0d0
          xi( i + mol( j )%nmol * 2 + offset ) = 0d0
          xi( i + mol( j )%nmol * 3 + offset ) = 0d0
          xi( i + mol( j )%nmol * 4 + offset ) = 0d0
          xi( i + mol( j )%nmol * 5 + offset ) = 0d0
          xi( i + mol( j )%nmol * 6 + offset ) = 0d0
       else
          !
          !各剛体分子に加わる並進/回転の力を算出
          !
          call rigid_collectforce(rigid(j),mol(j),site)
          DO I=1,mol(j)%nmol
             !
             !並進では、力が変位と直接関係し、座標は必要ない。
             !
             xi( i + mol( j )%nmol * 0 + offset ) = &
                  -rigid(j)%mol(i)%forcex * rigid( j )%massi
             xi( i + mol( j )%nmol * 1 + offset ) = &
                  -rigid(j)%mol(i)%forcey * rigid( j )%massi
             xi( i + mol( j )%nmol * 2 + offset ) = &
                  -rigid(j)%mol(i)%forcez * rigid( j )%massi
             !
             !Unserializeしてquaternionを入手する。
             !
             qa = p( i + mol( j )%nmol * 3 + offset )
             qb = p( i + mol( j )%nmol * 4 + offset )
             qc = p( i + mol( j )%nmol * 5 + offset )
             qd = p( i + mol( j )%nmol * 6 + offset )
             !
             !トルクを入手
             !トルクから、各軸回りの角加速度が求められる。
             !静止している物体であれば、角加速度は単純にトルクに比例する。
             !
             tx = rigid( j )%mol(i)%torquex * rigid( j )%ixxi
             ty = rigid( j )%mol(i)%torquey * rigid( j )%iyyi
             tz = rigid( j )%mol(i)%torquez * rigid( j )%izzi
             !
             !角加速度を、角速度とよみかえる→四元数の速度に変換
             !
             da = -( -tx*qb + ty*qc - tz*qd ) * 0.5d0
             db = -(  tx*qa - ty*qd - tz*qc ) * 0.5d0
             dc = -( -tx*qd - ty*qa + tz*qb ) * 0.5d0
             dd = -(  tx*qc + ty*qb + tz*qa ) * 0.5d0
             !
             !Serialize
             !
             xi( i + mol( j )%nmol * 3 + offset ) = da
             xi( i + mol( j )%nmol * 4 + offset ) = db
             xi( i + mol( j )%nmol * 5 + offset ) = dc
             xi( i + mol( j )%nmol * 6 + offset ) = dd
#ifdef DEBUG
             write(STDOUT,*) da,db,dc,dd, tx, ty, tz, qa, qb, qc, qd
             stop
#endif
          enddo
       endif
       offset = offset + mol( j )%nmol * 7
    enddo
  end subroutine rigidquench_dfuncq

  function rigidquench_mode()
    integer :: rigidquench_mode
    rigidquench_mode = mode
  end function rigidquench_mode

  function rigidquench_ncompo()
    integer :: rigidquench_ncompo
    rigidquench_ncompo = ncompo
  end function rigidquench_ncompo

end module rigidquench_module

!
!quench_moduleは、Numerical Recipeのプログラムをそのまま使用。commonな
!どf90にそぐわない点だけを変更した。
!
module quench_module
  use rigidquench_module
  implicit none
  integer, parameter, private ::  NMAX=MAXMOL*7
  real(kind=8),private        ::  pcom(NMAX),xicom(NMAX)
  integer,private             ::  ncom

  integer               , private :: ecount, fcount

contains

  SUBROUTINE linmin(p,xi,n,fret)
    INTEGER      :: n
    real(kind=8) :: fret,p(n),xi(n)
    real(kind=8), parameter :: TOL=1.d-4
    INTEGER      :: j
    real(kind=8) :: ax,bx,fa,fb,fx,xmin,xx
#ifdef VERBOSE
    write(STDERR,*) "linmin",n
#endif
    ncom=n
    do j=1,n
       pcom(j)=p(j)
       xicom(j)=xi(j)
    enddo
    ax=0.d0
    !     xx=1.d0
    !     fit for the case of water
    xx=0.001d0
    call mnbrak(ax,xx,bx,fa,fx,fb)
    fret=brent(ax,xx,bx,TOL,xmin)
    do j=1,n
       xi(j)=xmin*xi(j)
       p(j)=p(j)+xi(j)
    enddo
    !      write(7,*) xmin
  END SUBROUTINE linmin

  function f1dim(x)
    real(kind=8)         :: f1dim
    real(kind=8),intent(in) :: x
    real(kind=8)         :: xt(NMAX)
    INTEGER              :: j
#ifdef VERBOSE
    write(STDERR,*) "f1dim",ncom
#endif
    do j=1,ncom
       xt(j)=pcom(j)+x*xicom(j)
#ifdef VERBOSE
       if(j.le.10)write(STDERR,*) j,pcom(j),xicom(j)
#endif
    enddo
    !ここでは、forceを計算する必要がない。
    f1dim=func(xt, .true.)
  END function f1dim

  SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc)
    real(kind=8) :: ax,bx,cx,fa,fb,fc
    real(kind=8), parameter :: GOLD=1.618034d0, GLIMIT=100.d0, TINY=1.d-20
    real(kind=8) :: dum,fu,q,r,u,ulim
#ifdef VERBOSE
    write(STDERR,*) "mnbrak"
#endif
    fa=f1dim(ax)
    fb=f1dim(bx)
    if(fb.gt.fa)then
       dum=ax
       ax=bx
       bx=dum
       dum=fb
       fb=fa
       fa=dum
    endif
    cx=bx+GOLD*(bx-ax)
    fc=f1dim(cx)
1   if(fb.ge.fc)then
       r=(bx-ax)*(fb-fc)
       q=(bx-cx)*(fb-fa)
       u=bx-((bx-cx)*q-(bx-ax)*r)/(2.d0*sign(max(abs(q-r),TINY),q-r))
       ulim=bx+GLIMIT*(cx-bx)
       if((bx-u)*(u-cx).gt.0.d0)then
          fu=f1dim(u)
          if(fu.lt.fc)then
             ax=bx
             fa=fb
             bx=u
             fb=fu
             return
          else if(fu.gt.fb)then
             cx=u
             fc=fu
             return
          endif
          u=cx+GOLD*(cx-bx)
          fu=f1dim(u)
       else if((cx-u)*(u-ulim).gt.0.d0)then
          fu=f1dim(u)
          if(fu.lt.fc)then
             bx=cx
             cx=u
             u=cx+GOLD*(cx-bx)
             fb=fc
             fc=fu
             fu=f1dim(u)
          endif
       else if((u-ulim)*(ulim-cx).ge.0.d0)then
          u=ulim
          fu=f1dim(u)
       else
          u=cx+GOLD*(cx-bx)
          fu=f1dim(u)
       endif
       ax=bx
       bx=cx
       cx=u
       fa=fb
       fb=fc
       fc=fu
       goto 1
    endif
  END SUBROUTINE mnbrak

  FUNCTION brent(ax,bx,cx,tol,xmin)
    INTEGER ITMAX
    real(kind=8) :: brent,ax,bx,cx,tol,xmin,CGOLD,ZEPS
    PARAMETER (ITMAX=100,CGOLD=.3819660d0,ZEPS=1.0d-10)
    INTEGER iter
    real(kind=8) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
#ifdef VERBOSE
    write(STDERR,*) "brent"
#endif
    a=min(ax,cx)
    b=max(ax,cx)
    v=bx
    w=v
    x=v
    e=0.d0
    fx=f1dim(x)
    fv=fx
    fw=fx
    do iter=1,ITMAX
       xm=0.5d0*(a+b)
       tol1=tol*abs(x)+ZEPS
       tol2=2.d0*tol1
       if(abs(x-xm).le.(tol2-.5d0*(b-a))) goto 3
       if(abs(e).gt.tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.d0*(q-r)
          if(q.gt.0.d0) p=-p
          q=abs(q)
          etemp=e
          e=d
          if(abs(p).ge.abs(.5d0*q*etemp).or.&
               p.le.q*(a-x)             .or.&
               p.ge.q*(b-x))then
             goto 1
          endif
          d=p/q
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
       endif
1      if(x.ge.xm) then
          e=a-x
       else
          e=b-x
       endif
       d=CGOLD*e
2      if(abs(d).ge.tol1) then
          u=x+d
       else
          u=x+sign(tol1,d)
       endif
       fu=f1dim(u)
       if(fu.le.fx) then
          if(u.ge.x) then
             a=x
          else
             b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
       else
          if(u.lt.x) then
             a=u
          else
             b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
             v=w
             fv=fw
             w=u
             fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
             v=u
             fv=fu
          endif
       endif
    enddo
    write(STDERR,*) 'brent exceed maximum iterations'
    stop
3   xmin=x
    brent=fx
  end FUNCTION brent

  SUBROUTINE frprmn(p,n,ftol,iter,fret)
    INTEGER iter,n,ITMAX
    DOUBLE PRECISION fret,ftol,p(n)
    real(kind=8), parameter :: EPS=1.d-10
    INTEGER its,j,itn
    DOUBLE PRECISION dgg,fp,gam,gg,g(NMAX),h(NMAX),xi(NMAX)
    !     iterに繰り返し回数の上限を指示できるようにする。
#ifdef VERBOSE
    write(STDERR,*) "frprmn"
#endif
    ecount = 0
    fcount = 0
    itmax=iter
    fp=func(p, .false.)
    call dfunc(p,xi)
    do j=1,n
       g(j)=-xi(j)
       h(j)=g(j)
       xi(j)=h(j)
    enddo
    itn=1
    do its=1,ITMAX
       iter=its
       if(its.eq.itn)then
          write(STDERR,*) its,fp
          itn=itn*2
       endif
       call linmin(p,xi,n,fret)
       if(2.d0*abs(fret-fp).le.ftol*(abs(fret)+abs(fp)+EPS))then
          !          write(STDERR,*) "end",fret,fp
          return
       endif
       call normalize(p)
       fp=func(p, .false.)
       call dfunc(p,xi)
       gg=0.d0
       dgg=0.d0
       do j=1,n
          gg=gg+g(j)**2
          !     dgg=dgg+xi(j)**2
          dgg=dgg+(xi(j)+g(j))*xi(j)
       enddo
       if(gg.eq.0.d0)return
       gam=dgg/gg
       do j=1,n
          g(j)=-xi(j)
          h(j)=g(j)+gam*h(j)
          xi(j)=h(j)
       enddo
    enddo
    write(STDERR,*) 'frprmn maximum iterations exceeded'
  END SUBROUTINE frprmn

  function func(p, noforce)
    real(kind=8)             :: func
    real(kind=8), intent(in) :: p(*)
    logical,      intent(in) :: noforce
    ecount = ecount + 1
    func = rigidquench_func( p, noforce )
  end function func

  subroutine dfunc(p,xi)
    real(kind=8) :: p(*),xi(*)
    fcount = fcount + 1
    if ( rigidquench_mode() .eq. RQ_EULER0 ) then
       call rigidquench_dfunc( p, xi )
    else if ( rigidquench_mode() .eq. RQ_EULER ) then
       call rigidquench_dfunc2( p, xi )
    else if ( rigidquench_mode() .eq. RQ_QUAT ) then
       !
       !四元数を使って急冷する場合はこちらを使う。
       !(RigidModuleにも変更必要)
       !
       call rigidquench_dfuncq( p, xi )
    else
       stop
    endif
  end subroutine dfunc

  subroutine normalize(p)
    use rigidquench_module
    real(kind=8) :: p(*)
    if ( rigidquench_mode() .eq. RQ_QUAT ) then
       !
       !四元数を使って急冷する場合はこちらを使う。
       !(RigidModuleにも変更必要)
       !
       call q_normalize( p )
    endif
  end subroutine normalize

  subroutine getcount( e, f )
    integer, intent(out) :: e,f
    e = ecount
    f = fcount
  end subroutine getcount

end module quench_module
