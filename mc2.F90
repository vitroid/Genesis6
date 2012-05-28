! -*- f90 -*-
#undef DEBUG
#undef DEBUG2
#undef DEBUG3
#undef DEBUG4
!
!DEBUG20050322aは必須。計算精度がぜんぜんちがう。
!
#define DEBUG20050322a

#undef DEBUG20050309
!
!
!今のところ、ifort(8.1)+最適化(-O3)では正しい結果が得られないようだ。
!
!
!
#define OP
#ifdef OP
#define USE_FAKEMPI
#endif



program mc2
  use common_module
  use site_module
  use rigid_module
  use box_module
  use error_module
  use cutoff_module
  use standard_interaction_module
  use vector_module
  use book_module
  use grid_module
#ifdef OP
  use distrib_module
  use q6local2
  use neighbors_module
  use triangle
  use graph2
  use replica
#endif
  use set_module
  use mc2util
  use random_module

  use oplsmeoh_module
  use tip4p_module
  use tip5p_module
  use nvde_module
  use pairinfo_module
  use quat_module
  use interaction_module
  use mol_module
  use interaction2_module
  use monitor_module
  implicit none 
  external getnodeid
  integer :: getnodeid
  real(kind=8), parameter :: INCREMENT=1.05
  integer, parameter :: DEFAULT_SEED=5678
  type(sSystem) :: sys
  type(vector3) :: com( MAXMOL, MAXCOMPO ), comt(1)
  integer :: combi( MAXCOMBI ), combj( MAXCOMBI )

  type(sInteraction)  :: wwt
  type(sInteraction)  :: ww( 3 ) !あいかわらずおかしい。 ( MAXCOMBI )と書くとSegる。
  type(sMol) :: molt
  type(sRigid) :: rigidt
  type(sStdInt) :: sit
  type(sPairInfo)    :: pairinfo( 3 ),pairtest(MAXCOMBI)
  type(sAddress) :: addr(MAXCOMPO)
  integer          :: i,j,k,loop, outerloop, nouterloop, ninnerloop, innerloop
  real(kind=8)     :: epsum,oldsum, testsum
  !real(kind=8), allocatable :: random(:)
  !real(kind=8), allocatable :: delta(:), rtratio(:)
  real(kind=8) :: debop, sumdeltat, sumdeltar, dt,dr, avgrot, avgtrans
  type(sQ6System), pointer :: debq6
  !
  !組み合わせの数
  !
  integer      :: ncombi
  !
  !分子の総数
  !
  integer      :: nmol
  !
  !エネルギー差分
  !
  real(kind=8) :: deltae, deltaope, deltaesum
  !
  !微小回転前の四元子
  !
  real(kind=8) :: olda,oldb,oldc,oldd
  !
  !微小回転後の四元子
  !
  real(kind=8) :: qa,qb,qc,qd,qq
  !
  !各軸回りの微小回転角度
  !
  real(kind=8) :: rotx,roty,rotz
  !
  !内側ループ変数
  !
  integer      :: trial
  !
  !動かす分子の番号
  !
  integer      :: target
  !
  !微小並進量
  !
  real(kind=8) :: deltax, deltay, deltaz
  !
  !default radius to cut off interaction
  !
  real(kind=8) :: rc
  !
  !pair energy after trial
  !
  !real(kind=8),dimension(:),allocatable :: ep !,vr
  real(kind=8) :: ep(MAXMOL) !,vr
  !
  !true if accepted
  !
  logical      :: accept
  !
  !metropolis's ratio
  !
  real(kind=8) :: ratio
  !
  !number of trials in a loop ( usually set to twice of number of
  !molecules )
  !
  integer      :: ntrial
  integer          :: ii,jj,kk

  integer      :: compo

  logical :: verbose
  character(len=100) :: argv
#ifdef USE_FAKEMPI
  !
  !プロセス時間の均衡
  !
  integer :: deltatime, avgtime
  real    :: speed
  integer :: processtime(300), processloop(300)
#endif
#ifdef OP
  real(kind=8) :: op(MAXUMBRELLA), newop(MAXUMBRELLA)
  !type(pgraph), pointer :: subgraph
#endif
  !
  !replicaの状況
  !
  real(kind=8) :: rep_energy(300), rep_beta(300), rep_op(300*3)
  type(umbrella_parabolic) :: rep_umbrella(300*3)
  integer :: nexchange, n0, n00
  !
  !just for debug
  !
  real(kind=8) :: dum1, dum2
  !
  !replica exchange
  !
  !real(kind=8) :: acceptance0, acceptance1, dummy, a0,a1
  !
  !計算すべきOP領域1..Nを、どのノードが担っているか。読み込み時点では、領域順にノードにわりふられているものとする。(当面)平成17年2月11日(金)
  !
  integer :: node_in_charge(300)
  character(len=100)  ::  configfile
#ifdef USE_FAKEMPI
  integer :: ierr
#endif
  !write(6,*) box%mode
  !call new(nose)
  call new( sys )
  verbose = .false.
  !verbose = .true.

  do i=1,300
     node_in_charge(i) = i
  enddo
  !read common config from fort.10
  call getarg( 1, configfile )
  write(STDERR,*)"COMMON CONFIG FILE:", configfile
  open(10, file=trim(configfile))
  call ReadParams( sys, 10 )
  close (10)
  if ( 0 <= sys%myrank ) then
     open( 12, file=sys%infile )
     call ReadParams( sys, 12 )
     close( 12 )
  endif
#ifdef OP
  if ( .not. associated( sys%rxgraph ) ) then
     !
     !rxgraphが指定されていない場合は、従来通り線形に交換するのではなく、完全グラフを与える。
     !
     call new( sys%rxgraph )
     call graph_perfect( sys%rxgraph, sys%nprocs )
  end if
  call graph2pair( sys%rxgraph, sys%rxpair )
#endif

  do i=1,sys%mol(1)%nmol
     call qnormalize( &
          sys%rigid(1)%mol(i)%quat(1)%vec(1),&
          sys%rigid(1)%mol(i)%quat(1)%vec(2),&
          sys%rigid(1)%mol(i)%quat(1)%vec(3),&
          sys%rigid(1)%mol(i)%quat(1)%vec(4) )
  enddo
  if( sys%grid%mode .ne. GRID_NONE)then
     if(.not.active(sys%book))then
        call Book_Initialize(sys%book, BOOK_FIX, 0d0, 1)
     endif
  endif
  !
  !箱の大きさが指定されていてCutOffが指定されていない場合
  !
  if(sys%box%mode.ne.BOX_NONE.and.sys%cutoff%mode.eq.CUTOFF_NONE)then
     sys%cutoff%mode=CUTOFF_POLY
     rc=dmin1(sys%box%size%vec(1),sys%box%size%vec(2),sys%box%size%vec(3))*0.5d0
     call cutoff_initialize(sys%cutoff,rc-2d0,rc)
  endif
  !
  !いちおう変な設定をチェックしておく。
  !
  if( active(sys%book) .and. sys%grid%mode .eq. GRID_NONE)then
     write(STDERR,*) "ERROR: GRID IS NECESSARY FOR BOOKING."
     call die( 0, "mctest 2" )
  endif
  if( active(sys%book) .and. sys%cutoff%mode .eq. CUTOFF_NONE)then
     write(STDERR,*) "ERROR: CUTOFF IS NECESSARY FOR BOOKING."
     call die( 0, "mctest 3" )
  endif
  !
  !同種分子間の相互作用(分子数をnとするとn(n+1)/2対)
  !
  do j=1,sys%ncompo
     !
     !第2パラメータを省略した場合は自己相互作用とみなす。
     !
     !DECではoptionalの扱いがうまくいかないらしいので明示的に-1を与える。
     !
     call interaction_initialize(ww(j),sys%mol(j)%nmol,-1)
     combi( j ) = j
     combj( j ) = j
  enddo
  !
  !異種分子間の相互作用(分子数をn,mとするとnm対)
  !
  ncombi = sys%ncompo
  do i=1,sys%ncompo
     do j=i+1,sys%ncompo
        ncombi = ncombi + 1
        combi( ncombi ) = i
        combj( ncombi ) = j
        call interaction_initialize(ww(ncombi),sys%mol(i)%nmol,sys%mol(j)%nmol)
     enddo
  enddo
  write(*,*) "j",ncombi,sys%ncompo
  !
  !系が小さい場合を想定;全分子同士が相互作用する(グリッド分割しない)
  !
  if( active(sys%book) )then
     if( sys%book%margin .le. 0d0 )then
        sys%book%mode = BOOK_FIX
     endif
     if(sys%book%mode .eq. BOOK_AUTOINTERVAL)then
        call die( 0, "Main 27" )
        !
        !まだ考えていない。帳簿法はほとんど使わないから無視しても構わ
        !ないのだが・・・!!!
        !
        !com0(1:water%nmol)=rwater%mol(1:water%nmol)%com
        !call Book_ResetInterval(bk,co%Out+bk%margin,water%nmol,com0)
     endif
     !
     !初期エネルギー計算のために、とりあえず挿入(あとで除去しろ)
     !
     do j=1,ncombi
        call interaction_alltoall(ww(j))
     enddo
  else
     !
     !系が小さい場合を想定;全分子同士が相互作用する(グリッド分割しない)
     !
     do j=1,ncombi
        call interaction_alltoall(ww(j))
        !write(STDERR,*) "PAIR:", j, ww(j)%npair
     enddo
  endif
  if( sys%grid%mode.ne.GRID_NONE )then
     do i=1,sys%ncompo
        call Address_Initialize( addr(i),sys%mol(i)%nmol )
     enddo
  endif
  !
  !make interaction table
  !
  do j=1, sys%ncompo
     com( 1:sys%mol(j)%nmol, j ) = sys%rigid( j )%mol(1:sys%mol(j)%nmol)%com
  enddo
  !
  !set all site position
  !
  do j=1,sys%ncompo
     call rigid_setsiteposition(sys%rigid(j),sys%mol(j),sys%site)
  enddo
  pairinfo(1)%ep(:,:) = 0d0
  !pairinfo(1)%vr(:,:) = 0d0
  epsum = 0d0
  do k=1,ncombi
     i = combi( k )
     j = combj( k )
     !
     !相互作用対表の圧縮(相対座標が遠いものを落とす)
     !
     if ( sys%cutoff%mode.ne.CUTOFF_NONE) then
        call Interaction_Compress(ww(k),com(1,i),com(1,j),&
             .false. ,CutOff_OuterLimit(sys%cutoff),sys%box)
     else
        call Interaction_Compress(ww(k),com(1,i),com(1,j),&
             .false. ,-1d0,sys%box)
     endif
     !
     !Smooth cutoffによる減衰係数を事前計算
     !
     !call CutOff_Dimmer(co, ww(k), site, mol(i), mol(j) )
     call CutOff_Dimmer2( sys%cutoff, ww(k), com(1,i), com(1,j) )
     !
     !LJ+Coulombの標準相互作用の場合
     !
     call pairenergy(ww(k),sys%mol(i),sys%mol(j),sys%site,pairinfo(k), .false. ,sys%stdint(i),sys%stdint(j))
     epsum = epsum + pairinfo(k)%epsum
  enddo
  !epsum = epsum

  !
  !Count total number of molecules
  !
  nmol = 0
  do j=1,sys%ncompo
     nmol = nmol + sys%mol(j)%nmol
  enddo
  !
  !Trial専用に設けた分子セットの初期化を行う。
  !
  call new( molt, sys%mol(1)%nsite, sys%mol(1)%dof, RIGID_MODE, .not. FIXED )
  if ( sys%mol(1)%id(1:8) == "TIP4P   " ) then
     call rigid_tip4p_constructor( rigidt )
     call tip4p_setinteraction( sit )
  else if ( sys%mol(1)%id(1:8) == "NVDE____" ) then
     call rigid_nvde_constructor( rigidt )
     call nvde_setinteraction( sit )
  endif
  call mol_allocate(molt, sys%site, 1, RIGID_MODE )
  call interaction_initialize( wwt, sys%mol(1)%nmol, molt%nmol )
  call interaction_alltoall( wwt )
  !allocate( ep(wwt%npair0) )
  !allocate( vr(wwt%npair0) )
  if ( .not. associated( sys%delta ) ) then
     allocate( sys%delta(sys%mol(1)%nmol) )
     allocate( sys%rtratio(sys%mol(1)%nmol) )
     sys%delta(:)   = 1d-1
     sys%rtratio(:) = 1d0
  endif
  !delta(:) = 0d0
  if ( .not. associated( sys%rand ) ) then
     call random_initialize( sys%rand, DEFAULT_SEED )
  endif
  !
  !主loop
  !
  ntrial = 2*sys%mol(1)%nmol
  !allocate(random(ntrial*8+1))
  !call mtrngi( sys%seed )
  !
  !Fileのオープン
  !
  open( TRAJOUT, file=sys%outfile )
  !
  !sys%ninnerloopstepに一回程度、Replica Exchangeする。
  !ただし、あくまで平均値であり、実際はReplica Exchangeのタイミングが一致するように
  !ループ回数を動的に制御する。
  !
  loop = 0
  nouterloop = (sys%nloop-1) / sys%ninnerloop + 1
  ninnerloop = sys%nloop
  if ( sys%ninnerloop < sys%nloop ) then
     ninnerloop = sys%ninnerloop
  endif
  do outerloop=1, nouterloop
     !
     !reset timer
     !
#ifdef USE_FAKEMPI
     call fakempif_timer_delta( deltatime )
#endif
     sumdeltat = 0d0
     sumdeltar = 0d0
     do innerloop=1, ninnerloop
        !if ( 2 <= innerloop ) then
        !   pairtest(1) = pairinfo(1)
        !   oldsum = epsum
        !endif
        !write(STDERR,*) "LOOP", loop
        !
        ![0,1)乱数
        !
        !call mtrndv(random,ntrial*8+1)
        !
        !とりあえず、1MC stepに一回、グリッドを再設計する。
        !グリッドのサイズにはすこし余裕をみておく。
        !
        !
        !Serialize the center-of-molecules.
        !
        do compo=1, sys%ncompo
           nmol = sys%mol(compo)%nmol
           if ( sys%mol(compo)%isRigid ) then
              com( 1:nmol, compo ) = sys%rigid( compo )%mol(1:nmol)%com
           else
              do i=1,nmol
                 k = sys%mol( compo )%offset + i*sys%mol( compo )%nsite
                 com( i, compo )%vec(1) = sys%site%x(k)
                 com( i, compo )%vec(2) = sys%site%y(k)
                 com( i, compo )%vec(3) = sys%site%z(k)
              enddo
           endif
        enddo
        if ( sys%grid%mode .ne. GRID_NONE ) then
           !
           !grid and addr are local in this block
           !
           !
           !Decide the size of grid
           !
           call Grid_Update( sys%grid, sys%box, sys%book%outer )
           !
           !Assign molecules to the cell
           !
           do i=1,sys%ncompo
              call Address_Assign( addr(i), sys%box, sys%grid, sys%mol(i)%nmol, com(1,i) )
           enddo
           !
           !Flush pairlist
           !
           do k=1,ncombi
              call Interaction_Reset(ww(k))
              i=combi(k)
              j=combj(k)
              !
              !分子間相互作用対リストの作成
              !
              if ( ww(k)%isomol ) then
                 call Grid_NeighborList(sys%grid,sys%grid%homo,ww(k),addr(i),addr(j))
              else
                 call Grid_NeighborList(sys%grid,sys%grid%hetero,ww(k),addr(i),addr(j))
              endif
           enddo
           if( sys%book%mode.eq.BOOK_AUTOINTERVAL )then
              !
              !帳簿のマージン半径内に分子がどのくらい侵入してきたかを
              !チェックして、次に帳簿を更新すべきタイミングを決定する。
              !
              !call Book_AssumeInterval( sys%book, ww%nmol_i, com0, Book_BestMargin(sys%book,ww%nmol_i,com0) )
           endif
           sys%book%countdown=sys%book%interval
        endif
        testsum = 0d0
        do k=1,ncombi
           i=combi(k)
           j=combj(k)
           !call Interaction_Compress2(ww(k), active( box ),&
           !           box ,co%out**2, com(1,i),com(1,j),.false.)
           !   enddo
           if ( sys%cutoff%mode.ne.CUTOFF_NONE) then
              call Interaction_Compress(ww(k),com(1,i),com(1,j),&
                   .false. ,CutOff_OuterLimit(sys%cutoff),sys%box)
           else
              call Interaction_Compress(ww(k),com(1,i),com(1,j),&
                   .false. ,-1d0,sys%box)
           endif
#undef DEBUG20050322
#ifdef DEBUG20050322
           !
           !どうもトータルエネルギーがちゃんと計算できていない(保存されていない)気がするので、ここで完全に計算しなおしてみる。平成17年3月21日(月)ほぼ問題ない。
           !
           call CutOff_Dimmer2( sys%cutoff, ww(k), com(1,i), com(1,j) )
           !
           !LJ+Coulombの標準相互作用の場合
           !
           call pairenergy(ww(k),sys%mol(i),sys%mol(j),sys%site,pairinfo(k), .false. ,sys%stdint(i),sys%stdint(j))
           testsum = testsum + pairinfo(k)%epsum
           if ( k == ncombi ) then
              if ( 1d-12 < real_compare( testsum, epsum ) )then
                 write(STDERR,*) "EPSUM", real_compare( testsum, epsum ), testsum, epsum
              endif
           endif
#endif
        enddo

        
        !call box_save(box, STDERR, BOX_BOX3 )
#ifdef OP
        do i=1, sys%numbrella
           if ( sys%optype(i) == OP_TRIANGLE ) then
              call sTriangle_prepareall( sys%tri(i)%p, sys%mol(1)%nmol, ww(1), com(1,1), sys%refhist, op(i) )
              !write(STDERR,*) i,":TRIANGLE"
           else if ( sys%optype(i) == OP_Q6T ) then
              call sQ6System_prepareall( sys%q6s(i)%p, sys%mol(1)%nmol, ww(1), com(1,1), op(i) )
              ! call new( q6d(i)%p )
              !write(STDERR,*) i,":Q6"
           else if ( sys%optype(i) == OP_EP ) then
              op(i) = epsum
           endif
        enddo
        !write(STDERR,fmt1) trial, (op(i),i=1,sys%numbrella)
#ifdef DEBUG4
        if ( 2 <= innerloop ) then
           !dump all suspicious variables
           do i=1, sys%mol(1)%nmol
              do j=1, sys%mol(1)%nmol
                 if ( 1d-10 < real_compare( pairinfo(1)%ep(i,j), pairtest(1)%ep(i,j) ) ) then
                    write(STDERR,*) i,j,"PAIRINFO",real_compare( pairinfo(1)%ep(i,j), pairtest(1)%ep(i,j) )
                 end if
              enddo
           enddo
           if ( 1d-10 < real_compare( epsum, oldsum ) ) then
              write(STDERR,*) "EPSUM"
           endif
        endif
#endif
#endif /*OP*/



        
        do trial=0,ntrial-1
           !
           !
           !trial
           !
           !分子を選ぶ
           !
           target = random_getnext( sys%rand )*sys%mol(1)%nmol + 1
#ifdef DEBUG20050309
           write(6,*) "TRIAL", trial, target, (op(i),i=1,sys%numbrella)
#endif
           !if( trial==0 ) then
           !   write(STDERR,*) "TRIAL", trial, target
           !endif
           !write(STDERR,*) "TARGET:", target
           dt = sys%delta(target)
           dr = sys%delta(target) * sys%rtratio(target)
           sumdeltat = sumdeltat + dt
           sumdeltar = sumdeltar + dr
           !write(6,*) "target",target
           !
           !各軸回りの微小回転から四元数を求めているので、異方性は生まれないはず。
           !
           !rotx   = dr * ( random_getnext( sys%rand ) - 0.5d0 )
           !roty   = dr * ( random_getnext( sys%rand ) - 0.5d0 )
           !rotz   = dr * ( random_getnext( sys%rand ) - 0.5d0 )
           rotx   = random_normal( 0d0, dr**2, random_getnext( sys%rand ), random_getnext( sys%rand ) )
           roty   = random_normal( 0d0, dr**2, random_getnext( sys%rand ), random_getnext( sys%rand ) )
           rotz   = random_normal( 0d0, dr**2, random_getnext( sys%rand ), random_getnext( sys%rand ) )
           avgrot = ( abs(rotx) + abs(roty) + abs(rotz) ) / 3d0
           !こいつがちゃんと四元数をnormalizeしている。
           call qinfrotator( rotx,roty,rotz, qa,qb,qc,qd )
           olda   = sys%rigid(1)%mol( target )%quat(1)%vec(1)
           oldb   = sys%rigid(1)%mol( target )%quat(1)%vec(2)
           oldc   = sys%rigid(1)%mol( target )%quat(1)%vec(3)
           oldd   = sys%rigid(1)%mol( target )%quat(1)%vec(4)
           call qadd( qa,qb,qc,qd, olda,oldb,oldc,oldd )
           !
           !##rigidt
           !
           rigidt%mol( 1 )%quat( 1 )%vec(1) = qa
           rigidt%mol( 1 )%quat( 1 )%vec(2) = qb
           rigidt%mol( 1 )%quat( 1 )%vec(3) = qc
           rigidt%mol( 1 )%quat( 1 )%vec(4) = qd
           !qq = qa**2 + qb**2 + qc**2 + qd**2
           !if ( 1.001 < qq ) write(STDERR,*) "QQ",qq
           !deltax = dt * ( random_getnext( sys%rand ) - 0.5d0 )
           !deltay = dt * ( random_getnext( sys%rand ) - 0.5d0 )
           !deltaz = dt * ( random_getnext( sys%rand ) - 0.5d0 )
           !use box-muller method for generating normal random numbers
           deltax   = random_normal( 0d0, dt**2, random_getnext( sys%rand ), random_getnext( sys%rand ) )
           dum1 = random_getnext( sys%rand )
           dum2 = random_getnext( sys%rand )
           !if ( dum1 == 0d0 ) call die( 0, "Conguraturations!!!" )
           !deltay   = random_normal( 0d0, dt**2, random_getnext( sys%rand ), random_getnext( sys%rand ) )
           deltay   = random_normal( 0d0, dt**2, dum1, dum2 )
           deltaz   = random_normal( 0d0, dt**2, random_getnext( sys%rand ), random_getnext( sys%rand ) )
           avgtrans = ( abs(deltax) + abs(deltay) + abs(deltaz) ) / 3d0
           !if ( trial==0 ) write(STDERR,fmt1) sys%rand%next, sys%delta(target), sys%rtratio(target), deltax, deltay, deltaz
           !
           !##comt
           !
           comt(1)%vec(1) = deltax
           comt(1)%vec(2) = deltay
           comt(1)%vec(3) = deltaz
           comt(1)%vec(:) = comt(1)%vec(:) + sys%rigid(1)%mol( target )%com%vec(:)
           !
           !##rigidt
           !
           rigidt%mol( 1 )%com = comt(1)
           !
           !com()のほうも新しい座標にしておく
           !(Tetrahedralityの計算に必要)
           !
           !!!Reject時に復旧必要。
           !
           com( target, 1 ) = comt(1)
           call rigid_setsiteposition(rigidt,molt,sys%site)
           !
           !まず、rc+delta以内の隣接分子の表を作る。Compressが使えるはず。
           !本当はこの段階でtargetをpairlistから除外しておいた方がよい。
           !
           !
           !グリッド分割する場合は、グリッドの大きさを、相互作用距離+最大移
           !動距離よりも大きくとっておく。こうすれば、分子が移動した結果隣の
           !グリッドに移動してしまっても、移動前の隣接分子が突然見えなくなっ
           !てしまうことはない。
           !
           !グリッドを使わない場合も、半径(cutoff距離
           ! +移動距離)内の分子の一覧を作る。
           !
           !
           !##wwt
           !
           if ( sys%cutoff%mode.ne.CUTOFF_NONE) then
              call Interaction_Compress(wwt,com(1,1),comt,&
                   .false. ,CutOff_OuterLimit(sys%cutoff)+dt,sys%box)
           else
              call Interaction_Compress(wwt,com(1,1),comt,&
                   .false. ,-1d0,sys%box)
           endif
           !
           !Smooth cutoffによる減衰係数を事前計算
           !
           !
           !##wwt
           !
           !call CutOff_Dimmer(co, ww(k), site, mol(i), mol(j) )
           call CutOff_Dimmer2( sys%cutoff, wwt, com(1,1), comt )
           !
           !LJ+Coulombの標準相互作用の場合。
           !
           !
           !##ep
           !
           call pairpotential(wwt,sys%mol(1),molt,sys%site,ep,.false., sys%stdint(1), sit )

           deltae = 0d0
#ifdef DEBUG20050322a
           !
           !相互作用しない対のpairinfoが0であることを保証する。
           !
           do i=1,sys%mol(1)%nmol
              deltae = deltae - pairinfo(1)%ep(i,target)
           enddo
#endif
           do k=1,wwt%npair
              i = wwt%pair_i(k)
              if ( i /= target ) then
                 !!!相互作用していない対のpairinfoが0になっている保証はない。
#ifdef DEBUG20050322a
                 deltae = deltae + ep(k)
#else
                 deltae = deltae + ep(k) - pairinfo(1)%ep(i,target)
#endif
              endif
           enddo
           if ( verbose ) write(STDOUT,*) dt, dr, deltae*I2J*J2K

           !
           !Umbrella potentialのような外場はここで計算する。
           !
           !---------------------------------------------------------------------
#ifdef OP
           do i=1, sys%numbrella
              if ( sys%optype(i) == OP_TRIANGLE ) then
                 call sTriangle_Differentiate( sys%tri(i)%p, target, deltax, deltay, deltaz, sys%refhist, newop(i) )
              else if ( sys%optype(i) == OP_Q6T ) then
                 call sQ6System_Differentiate( sys%q6s(i)%p, target, deltax, deltay, deltaz, sys%mol(1)%nmol, sys%q6d(i)%p, newop(i) )
              else if ( sys%optype(i) == OP_EP ) then
                 newop(i) = epsum + deltae
              endif
           enddo
#endif
           !
           !---------------------------------------------------------------------
           !
           deltaesum = deltae
#ifdef OP
           deltaope = 0d0
           do i=1, sys%numbrella
              deltaope = deltaope + umbrella_potential( sys%umbrella(i), newop(i) ) - umbrella_potential( sys%umbrella(i), op(i) )
           enddo
           deltaesum = deltae + deltaope
#endif
           !
           !新しい分子配置をうけいれるかどうかを判定する。
           !
           accept = .false.
           if ( deltaesum < 0d0 ) then
              accept = .true.
           else
              ratio = dexp(-( deltaesum * I2J * J2K )*sys%beta)
              if ( random_getnext( sys%rand ) < ratio ) then
                 accept = .true.
                 !write(STDERR,*) "Ratio",ratio,deltae * I2J * J2K,"dE"
              endif
              !write(sys%myrank+70,*) ratio,deltae * I2J * J2K, accept
           endif
           !
           !うけいれた場合
           !
           if ( accept ) then
              if ( verbose ) write(6,*) "accept motion"
              !
              !座標を更新する。
              !
              sys%rigid(1)%mol( target ) = rigidt%mol(1)
              !
              !com(), Site Positionも1分子分だけ更新する。
              !
              !com()はtrial後の座標になってるはず。だから更新しなくていい。
              !com( target, 1 ) = rigidt%mol(1)%com
              i=sys%mol(1)%offset + (target-1)* sys%mol(1)%nsite
              j=molt%offset
              do k=1,sys%mol(1)%nsite
                 sys%site%x(i+k)=sys%site%x(j+k)
                 sys%site%y(i+k)=sys%site%y(j+k)
                 sys%site%z(i+k)=sys%site%z(j+k)
              enddo
              !
              !相互作用表を更新する。
              !
              epsum = epsum + deltae
#ifdef OP
#ifdef DEBUG3
              do i=1, sys%numbrella
                 if ( sys%optype(i) == OP_TRIANGLE ) then
                    call sTriangle_prepareall( sys%tri(i)%p, sys%mol(1)%nmol, ww(1), com(1,1), sys%refhist, op(i) )
                 else if ( sys%optype(i) == OP_Q6T ) then
                    call sQ6System_prepareall( sys%q6s(i)%p, sys%mol(1)%nmol, ww(1), com(1,1), op(i) )
                 else if ( sys%optype(i) == OP_EP ) then
                    op(i) = epsum
                 endif
              enddo
              !write(STDERR,fmt1) trial, epsum, op(1)
              if ( 1d-12 < dabs( op(1) - newop(1) ) ) stop
#else /*DEBUG3*/
              do i=1, sys%numbrella
                 if ( sys%optype(i) == OP_TRIANGLE ) then
                    call sTriangle_Integrate( sys%tri(i)%p )
                 else if ( sys%optype(i) == OP_Q6T ) then
                    call sQ6System_Integrate( sys%q6s(i)%p, sys%q6d(i)%p, sys%mol(1)%nmol, target )
                 else if ( sys%optype(i) == OP_EP ) then
                 endif
              enddo
#endif
              op(1:sys%numbrella) = newop(1:sys%numbrella)
              !write(STDERR,fmt1) trial, epsum, (op(i),i=1,sys%numbrella)
#endif /*OP*/
              pairinfo(1)%epsum = epsum
#ifdef DEBUG20050322a
              !
              !相互作用しない対のpairinfoが0であることを保証する。
              !
              do i=1,sys%mol(1)%nmol
                 pairinfo(1)%ep(i,target) = 0d0
                 pairinfo(1)%ep(target,i) = 0d0
              enddo
#endif
              do k=1,wwt%npair
                 i = wwt%pair_i(k)
                 if ( i /= target ) then
                    pairinfo(1)%ep(i,target) = ep(k)
                    pairinfo(1)%ep(target,i) = ep(k)
                 endif
              enddo
              !
              !移動量を更新する。
              !
              sys%delta( target ) = sys%delta( target ) * INCREMENT
              !
              !回転の大きさが適正になるように制御する。
              !
              if ( avgtrans * sys%rtratio(target) < avgrot ) then
                 sys%rtratio(target) = sys%rtratio(target) * INCREMENT
              else
                 sys%rtratio(target) = sys%rtratio(target) / INCREMENT
              endif
           else
              if ( verbose ) write(6,*) "reject motion"
              !
              !com()の内容をもとにもどす。
              !
              com( target, 1 )%vec(:) = sys%rigid(1)%mol( target )%com%vec(:)
              sys%delta( target ) = sys%delta( target ) / INCREMENT
              if ( sys%delta( target ) < 0.05 ) sys%delta( target ) = 0.05
              !
              !大きめに回転することでrejectされたなら、回転の比率をより減らす。
              !
              !if ( avgtrans * rtratio(target) < avgrot ) then
              !   rtratio(target) = rtratio(target) / INCREMENT
              !else
              !   rtratio(target) = rtratio(target) * INCREMENT
              !endif
#ifdef OP
              !
              !tri%newhistを捨てる
              !
              do i=1, sys%numbrella
                 if ( sys%optype(i) == OP_TRIANGLE ) then
                    call sTriangle_Reject( sys%tri(i)%p )
                 else if ( sys%optype(i) == OP_Q6T ) then
                    call sQ6System_Reject( sys%q6s(i)%p, sys%q6d(i)%p, target, sys%mol(1)%nmol )
                 else if ( sys%optype(i) == OP_EP ) then
                 endif
              enddo
#endif
           endif
#ifdef DEBUG20050309
           do i=1, sys%numbrella
              if ( sys%optype(i) == OP_Q6T ) then
                 write(6,*)"DEBUG5"
                 call new( debq6 )
                 call sQ6System_prepareall( debq6, sys%mol(1)%nmol, ww(1), com(1,1), debop )
                 call sQ6System_compare( debq6, sys%q6s(i)%p, sys%mol(1)%nmol, STDOUT )
                 call done( debq6 )
                 ! call new( q6d(i)%p )
                 write(6,*)"DEBUG5"
              else if ( sys%optype(i) == OP_EP ) then
              endif
           enddo
#endif
#ifdef DEBUG20050309
           do i=1, sys%numbrella
              if ( sys%optype(i) == OP_Q6T ) then
                 write(6,*)"DEBUG5"
                 call new( debq6 )
                 call sQ6System_prepareall( debq6, sys%mol(1)%nmol, ww(1), com(1,1), debop )
                 call sQ6System_compare( debq6, sys%q6s(i)%p, sys%mol(1)%nmol, STDOUT )
                 call done( debq6 )
                 ! call new( q6d(i)%p )
                 write(6,*)"DEBUG5"
              else if ( sys%optype(i) == OP_EP ) then
              endif
           enddo
#endif

        enddo
        !
        !calculate HB network
        !
        !
        !Save configuration
        !@SNPOで指定すべき。
        !
        if( 0 < sys%snpo )then
           if(mod(loop,sys%snpo).eq.0)then
              call writetag( TRAJOUT, "@NCMP" )
              write(TRAJOUT,*) sys%ncompo
              call Box_Save( sys%box, TRAJOUT, BOX_BOX3 )
              call Rigid_Save( sys%rigid(1), sys%mol(1), TRAJOUT, RIGID_NX4A )
           endif
        endif
        if ( 0 < sys%nlog ) then
           if ( mod(loop, sys%nlog) == 0 ) then
              call writetag( TRAJOUT, "@MC2O" )
              !
              !MC2O型式では、エネルギーはすべて内部単位系になっている。
              !本当は、エネルギーはK単位とし、最後にumbrella重みをいっしょに表示してほしい。
              !
              write(TRAJOUT,fmt1) loop, epsum, 1d0/sys%beta, ( op(i), sys%umbrella(i)%center, sys%umbrella(i)%coeffi, i=1, sys%numbrella )
           endif
        endif
        
#ifdef OP
        do i=1, sys%numbrella
           if ( sys%optype(i) == OP_TRIANGLE ) then
              call sTriangle_EndAStep( sys%tri(i)%p )
           else if ( sys%optype(i) == OP_EP ) then
           endif
        enddo
#endif
        loop = loop + 1
     enddo ! innerloop
     sumdeltat = sumdeltat / ( ninnerloop * ntrial )
     sumdeltar = sumdeltar / ( ninnerloop * ntrial )
#ifdef USE_FAKEMPI
     call fakempif_timer_delta( deltatime )
     !
     !計算時間の均衡化
     !
     speed = ninnerloop
     speed = speed / deltatime
     write( STDERR,* ) sys%myrank, "(", getnodeid(), ") took", deltatime, "ms for", ninnerloop, "loops", speed, sumdeltat, sumdeltar
     call fakempif_gather( 4, ninnerloop, processloop, ierr )
     call fakempif_gather( 4, deltatime , processtime, ierr )
#endif
#ifdef OP
     if ( sys%myrank == 0 ) then
        write( STDERR,* ) "Load balancing..."
        !
        !まず、全ノードの計算時間を求め、sys%ninnerloopステップ分計算するのに必要な平均時間を求める。
        !
        avgtime = 0
        do i=1, sys%nprocs
           avgtime = avgtime + sys%ninnerloop * processtime(i) / processloop(i)
        enddo
        avgtime = avgtime / sys%nprocs
        !
        !今度は、各ノードの計算速度をもとに、avgtime間で計算できるループ数を求める。
        !
        do i=1, sys%nprocs
           deltatime = processtime(i)
           if ( 100000 < deltatime ) deltatime = 100000 !safety valve for suspend
           processloop(i) = processloop(i) * avgtime / processtime(i)
        enddo
        write( STDERR,* ) "Done."
     endif
#endif
#ifdef USE_FAKEMPI
     call fakempif_scatter( 4, ninnerloop, processloop, ierr )
     call fakempif_timer_delta( deltatime )
     write( STDERR,* ) deltatime, "ms for scheduling task."
     !write(STDERR,*) "Next", sys%myrank, ninnerloop
     !
     !次にreplica exchange。交換するか否かの判断は親ノードが行い、結果をscatterする。
     !無駄な通信がある(親が保持していれば、gatherする必要がないもの)が、当面安全性を重視する。
     !
     call fakempif_gather( 8, epsum, rep_energy, ierr )
     call fakempif_gather( 8*sys%numbrella, op,    rep_op,     ierr )
     call fakempif_gather( 8, sys%beta,   rep_beta,     ierr )
     call fakempif_gather( umbrella_size*sys%numbrella, sys%umbrella, rep_umbrella,     ierr )
#endif
     !write(STDERR,*) sys%myrank,">", epsum, op, sys%beta, sys%center, sys%coeffi
#ifdef OP
     if ( sys%myrank == 0 ) then
        do i=1, sys%nprocs
           !
           !iはノードの番号ではなく、OP領域の番号。
           !その領域を担っているnodeはnode_in_charge(i)
           !
           n0 = node_in_charge(i)
           n00 = (n0-1) * sys%numbrella
           write(STDERR,'(i10,i3,1x,f13.4,1x,f5.0,1x,99(g10.3,1x,g10.3,1x,g10.3,1x,g10.3,1x))') &
                loop, &
                n0, &
                rep_energy(n0), &
                1d0/rep_beta(n0), &
                ( &
                  rep_op(n00+j), &
                  rep_umbrella(n00+j)%center, &
                  rep_umbrella(n00+j)%coeffi, &
                  umbrella_potential( rep_umbrella(n00+j), rep_op(n00+j) ), &
                  j=1, sys%numbrella &
                )
        enddo
        !
        !node 0でまとめて作業を行うのなら、Replica Exchangeを並列に実行する必要はない。
        if ( 1 < sys%nprocs ) then
           do nexchange = 1, sys%nexchange
              k  = random_getnext( sys%rand ) * sys%rxpair%nbond + 1
              j  = sys%rxpair%a(k)
              jj = sys%rxpair%b(k)
              !j,jjはノードの番号ではなく、OP領域の番号。
              !その領域を担っているnodeはnode_in_charge(j)
!              accept = replica_exchange_trial( j, jj, node_in_charge, sys%rand, sys%numbrella, rep_energy, rep_umbrella, rep_op, rep_beta )
              accept = replica_exchange_trial( j, jj, node_in_charge, sys%rand, rep_beta, rep_energy, sys%numbrella, rep_umbrella, rep_op )
              call writetag( TRAJOUT, "@MC2R" )
              n0 = node_in_charge(j)
              n00 = (n0-1) * sys%numbrella
              write(TRAJOUT,fmt1) nexchange, rep_energy(n0), 1d0/rep_beta(n0), ( rep_op(n00+i), rep_umbrella(n00+i)%center, rep_umbrella(n00+i)%coeffi, i=1, sys%numbrella )
              call writetag( TRAJOUT, "@MC2R" )
              n0 = node_in_charge(jj)
              n00 = (n0-1) * sys%numbrella
              write(TRAJOUT,fmt1) nexchange, rep_energy(n0), 1d0/rep_beta(n0), ( rep_op(n00+i), rep_umbrella(n00+i)%center, rep_umbrella(n00+i)%coeffi, i=1, sys%numbrella )
           enddo
        endif
     endif
#endif
#ifdef USE_FAKEMPI
     if ( 1 < sys%nprocs ) then
        call fakempif_scatter( 8, sys%beta,   rep_beta,     ierr )
        call fakempif_scatter( 8*2*sys%numbrella, sys%umbrella, rep_umbrella,     ierr )
     endif
     !write(STDERR,*) sys%myrank,">>", epsum, op, sys%beta, sys%center, sys%coeffi
     write( STDERR,* ) deltatime, "ms for replica exchanges."
#endif /*USE_FAKEMPI*/
  enddo
  !deallocate( random )
  !deallocate( ep )
  !deallocate( vr )
  !deallocate( delta )
  !deallocate( rtratio )
  !最終配置の出力、乱数の内部状態の出力
  call writetag( TRAJOUT, "@CONT" )  ! 継続計算のための情報の開始
  call writetag( TRAJOUT, "@DELT" )
  write(TRAJOUT,*) sys%mol(1)%nmol
  do i=1, sys%mol(1)%nmol
     write(TRAJOUT, fmt0 ) sys%delta(i), sys%rtratio(i)
  enddo
  call random_saveMTRN( sys%rand, TRAJOUT )
  do k=1, sys%numbrella
     call writetag( TRAJOUT, "@GEUM" )
     write(TRAJOUT,"(a8)") sys%id08(k)
     write(TRAJOUT,*) sys%umbrella(k)%center, sys%umbrella(k)%coeffi
  enddo
  call writetag( TRAJOUT, "@TEMP" )
  sys%temp = 1d0 / sys%beta
  write(TRAJOUT,*) sys%temp
  call Box_Save( sys%box, TRAJOUT, BOX_BOX3 )
  call Rigid_Save( sys%rigid(1), sys%mol(1), TRAJOUT, RIGID_NX4A )
  !
  !本当は乱数を保存すべき。
  !
#ifdef FAKEMPI
  call fakempi_finalize( ierr )
#endif
  
end program mc2

