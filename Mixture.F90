! -*- f90 -*-

!1成分(input512.snoda)、能勢あり、OK 平成14年12月10日(火)

!
!Rigidな分子の混合系汎用。もちろんNeatな系でも使える。
!
!ToDo
!Feedback temperature from replica exchange
!Final data in binary
!
!
!平成15年9月30日(火)
!  isRigidの周りの処理が、flex分子でもちゃんと動くかどうか。
!  flex moleculeのカットオフの種別
!     特定サイトを基準とするカットオフ
!     サイト群の重心を基準とするカットオフ
!        重心をカットオフの基準点とするなら、運動方程式の変更が必要では？


!it should work, but removed by matto 2012-5-28
!#define USE_REPLICA
!It fails to compile on gfortran.
!#define JARSYNSKI

program main
  use common_module

  use error_module
  use vector_module
  use rigid_average_module
  use box_module
  use cutoff_module
  use interaction_module
  use interaction2_module
  use mol_module
  use nose_module
  use rigid_module
  use flex_module
  use site_module
  use time_module
  use andersen_module
  use external_module
  use book_module
  use grid_module
  use monitor_module
  use standard_interaction_module
  use random_module
  use matrix_module
#ifdef RATTLE
  use rattle_module
#endif
#ifdef JARZYNSKI
  !inactivated 2012-5-28
  !use rigidmorph_module
#endif
  use stack_module
  use rigid_setup_module
  use flex_setup_module
  use triplet_module
  ! for scaling interaction
  use swsilicon_module
#ifdef USE_REPLICA
  use MDreplicaexchange
  use replica
#endif
  
  implicit none

  integer :: err
  integer :: combi( MAXCOMBI ), combj( MAXCOMBI )
  type(sTime) :: t
  type(sRigid) :: rigid( MAXCOMPO )
  type(sFlex) :: flex( MAXCOMPO )
  type(sMol) :: mol( MAXCOMPO )
  type(sInteraction) :: ww(MAXCOMBI)
  type(sBind) :: bind( MAXCOMBI )   ! 水同士をゆるく結びつけ、クラスタが壊れないようにする。
  type(sBind) :: jmon( MAXCOMBI )  ! 自由エネルギー計算用に、別のパラメータでも計算する。
#ifdef SOLVATIONTEST
  !
  !水素結合を固定する。水に特化した。
  !
  type(sJointHB) :: joint( MAXCOMBI )
  real(kind=8), allocatable :: ep1sol(:)
  real(kind=8)   :: ejointhb
#endif
  type(sAddress) :: addr(MAXCOMPO)
  type(sSite) :: site
  type(sNose) :: nose
  type(sAndersen) :: an
  type(sBox)  :: box
  type(sCutoff) :: co
  type(sRigidAverage) :: av
  type(sBook) :: bk
  type(sGrid) :: g
  type(sStdInt) :: si(MAXCOMPO)
  type(pRandom), pointer :: random
  !
  !SWSilicon 3-body interaction
  !
  type(sTriplet)              :: triplet
  real(kind=8)                :: vrmat3(3,3)
  real(kind=8), allocatable   :: ep3(:)
  integer                     :: tri
#ifdef RATTLE
  type(sRattle) :: rattle
#endif
  type(sTensor) :: virialp( MAXCOMBI )
  type(sTensor) :: virialk( MAXCOMPO )
  !圧力テンソルの運動量成分とポテンシャル成分。See JPCB 104, 1332(2000)
  type(sTensor) :: vpk, vpp
#ifdef DEBUG
  real(kind=8)  :: virialpsum, virialksum
#endif
  !
  !Volume x Pressure Tensor
  !
  type(sTensor) :: vp
  !
  !分子間相互作用に比例係数を掛ける
  !
  real(kind=8) :: ivscale( MAXCOMPO, MAXCOMPO ), value
  real(kind=8) :: ScaleSWSI
  !
  !分子を特定座標にバネでつなぐ。
  !
  type(sTie)  :: tie(MAXCOMPO)
  integer     :: ntie
  real(kind=8) :: exsum,etie
  logical      :: ftie  !dummy
  !
  !共通バネ定数
  !
  real(kind=8) :: tiea
  integer      :: tieb
  real(kind=8),allocatable :: tmpa(:),tmpx(:)
  !
  !Maximum number of molecules in a component!!
  !
  integer     :: MNOMIAC
  !
  !第0ステップでの温度の強制設定(座標だけあたえて、速度をランダムにし
  !たい場合)
  !
  real(kind=8) :: settemp
  integer      :: seed
  !
  !系全体の回転の抑止の頻度(クラスタの場合に必須)
  !
  integer      :: stoprotation
  type(vector3) :: moment,totalmoment,omega
  real(kind=8)  :: tensor(3,3),totaltensor(3,3)
  !
  !系全体の並進の抑止の頻度(クラスタの場合に必須)
  !
  integer      :: stopdrift
  !
  !重心を座標原点に
  !
  integer      :: shifttoorigin
  real(kind=8) :: totalmass,mass
  !
  !柔軟分子の質点の質量
  !
  real(kind=8) :: molmass(MAXSITE)
#ifdef MEMUSAGE
  type(sMonitor_integer)      :: mon(12)
#endif
  integer :: i,j,k,loop,nloop,nstep,nlog,snpo,snpv,compo,mdvw, mdvwsite
  integer :: ptensor
  integer :: pair
  !
  !分子を初期座標のまま固定。直後に読みこむ座標にのみ有効。継続データに含まれる。
  !
  logical :: isFixed
  !
  !@NCMPで明示的に指定される成分数
  !
  integer :: incompo
  real(kind=8) :: ekt(MAXCOMPO),ekr(MAXCOMPO),eksum,temperature,ektsum,ekrsum
  real(kind=8) :: andersen_ep,andersen_ek,nose_ep,nose_ek,pressure
  real(kind=8) :: ep(MAXCOMBI), epsum
  !real(kind=8) :: ep(MAXCOMBI), virial(MAXCOMBI), epsum, virialsum
  real(kind=8) :: rc,volumepressure
  real(kind=8), allocatable :: ep1(:)
  !real(kind=8) :: ep1(MAXPAIR)
  character(len=5) :: tag
  character(len=8) :: lastid
  logical :: fDt,fLoop,fVcnt
  logical :: fAverage
  logical :: ffast
  integer :: dof
  integer :: iVcnt
  real(kind=8)  :: vcnt
  !
  !成分の数(とりあえず上限は2)
  !
  integer      :: ncompo
  !
  !組み合わせの数
  !
  integer      :: ncombi
  !type(vector4),dimension(MAXMOL) :: quat
  !type(vector3),dimension(MAXMOL) :: force
  type(vector3) :: com( MAXMOL, MAXCOMPO )
  !
  !外場によって重心に加わる力
  !
  type(vector3) :: fcom(MAXMOL)
  !
  !分子対を結びつける。
  !
  logical :: press, stretch
  integer :: ii,jj
  real(kind=8) :: ebind, ejmon, emonsum
  !
  !体積を再設定する。
  !
  real(kind=8) :: vscale
  !
  !分子数。一時変数
  !
  integer :: nmol
#ifdef RATTLE
  !
  !Rattle内部でのサイト番号
  !
  integer :: rsite
#endif
  !
  integer             :: morphcompo = 0
#ifdef JARZYNSKI
  !Jarzynskiによる自由エネルギー計算のための、Morphing
  !内部でのとりあつかいはFIXCと同一にする。
  !とりあえず単一成分のみmorphできるものとする。
  !
  type( sRigidMorph ) :: rigidmorph
  real(kind=8)        :: mixi, mixf, mix

  real(kind=8)        :: dudl, dudlsum
  !for debug
#ifdef DEBUGTX
  real(kind=8)        :: tx(1000), ty(1000), tz(1000)
#endif
  !
  !Jarzynski. どうやって汎化すればいいかよくわからないので、当面morphingだけを考慮する。
  !
  logical             :: fJar
  integer             :: iJar
  real(kind=8)        :: epjar     ! work in a step
  real(kind=8)        :: jarwork( MAXCOMBI )   ! total work
#endif ! JARZYNSKI
  !
  !write文を整理する。
  !
  type(real8stack)    :: output

  !
  !全エネルギーを所定の値に収束させる。
  !
  real(kind=8) :: hset_value, hset_ratio, hset_ep, hset_ek, hset_param
  integer      :: hset_freq

#ifdef USE_REPLICA
  !
  !プロセス時間の均衡
  !
  integer :: deltatime, avgtime
  real    :: speed
  integer :: processtime(300), processloop(300), processaccm(300)
  integer :: ierr, intercept, ninnerloop
  !
  !replicaの状況
  !
  real(kind=8) :: rep_energy(300), rep_beta(300)
  type(sMDReplica) :: mdreplica
  !
  !計算すべきOP領域1..Nを、どのノードが担っているか。読み込み時点では、領域順にノードにわりふられているものとする。(当面)平成17年2月11日(金)
  !
  integer :: node_in_charge(300)
  character(len=100)  ::  configfilename
  integer :: n0, nexchange 
  !
  !true if accepted
  !
  logical      :: accept
  logical      :: finallap
  !
  !ep in k
  !
  real(kind=8) :: beta0
#endif
  integer :: LOGFILE, CONFIGFILE
  !real(kind=8), dimension(MAXMOL) :: comx0,comy0,comz0,qa0,qb0,qc0,qd0
  !real(kind=8), dimension(MAXMOL) :: comx1,comy1,comz1,qa1,qb1,qc1,qd1
  !real(kind=8), dimension(MAXMOL) :: forcex,forcey,forcez
  !microcanonicalの全エネルギー。
  !     initialization
  !     うまく算程化したい。とりあえず、初期化すべき構造体をすべてそのまま
  !     わたすことにする。
  call new(nose)
  call new(site, MAXsite_total)
  call new(an)
  call new(av)
  call new(box)
  call new(co)
  call new(bk)
  call new(g)
  call new( output, 100 ) !出力用スタック
#ifdef MEMUSAGE
  call new( mon(1), MON_MAXIMUM, MAXGRIDX,&
       MON_HALT, "MAXGRIDX                  ", 0 )
  call new( mon(2), MON_MAXIMUM, MAXGRIDY,&
       MON_HALT, "MAXGRIDY                  ", 0 )
  call new( mon(3), MON_MAXIMUM, MAXGRIDZ,&
       MON_HALT, "MAXGRIDZ                  ", 0 )
  call new( mon(4), MON_MAXIMUM, MAXINGRID,&
       MON_HALT, "MAXINGRID                 ", 0 )
  call new( mon(5), MON_MAXIMUM, MAXNEIBORCELL,&
       MON_HALT, "MAXNEIBORCELL             ", 0 )
  call new( mon(6), MON_MAXIMUM, MAXMOL,&
       MON_HALT, "MAXMOL                    ", 0 )
  call new( mon(7), MON_MAXIMUM, MAXsite_total,&
       MON_HALT, "MAXSITE                   ", 0 )
  call new( mon(8), MON_MAXIMUM, MAXPAIR,&
       MON_HALT, "MAXPAIR                   ", 0 )
  call new( mon(9), MON_MAXIMUM, MAXCELL,&
       MON_HALT, "MAXCELL                   ", 0 )
  call new( mon(10), MON_MAXIMUM, MAXPARTNER,&
       MON_HALT, "MAXPARTNER                ", 0 )
#endif
#ifdef USE_REPLICA
  call new(mdreplica)
#endif
  !
  !相互作用に比例係数を掛ける。初期値は1
  !
  ivscale(:,:) = 1d0
  ScaleSWSI    = 1d0
  dof=0d0
  tiea = 0d0
  hset_freq = 0
  !
  !3つ組の初期化
  !
  triplet%size = 0

  lastid=""
  fDt=.false.
  fLoop=.false.
  fVcnt=.false.
  snpo = 0
  snpv = 0
  ptensor = 0
  mdvw = 0
  fAverage=.false.
  ffast=.false.
  vscale=1.0d0
  nlog=0
  nstep=0
  isFixed = .false.
#ifdef JARZYNSKI
  mixi = 0
  mixf = 0
  morphcompo  = 0
  dudlsum = 0d0
  rigidmorph%mode = rigidmorph_none
  fJar = .false.
#endif !JARZYNSKI
  !
  !全体の回転の停止
  !
  stoprotation = 0
  !
  !
  !
  stopdrift = 0
  !
  !重心を座標原点に
  !
  shifttoorigin = 0
  !
  !温度の再設定(<0の場合は設定しない)
  !
  settemp = -1d0
  !
  !binaryから読みこんだtieの成分数
  !
  ntie = 0
  mnomiac=0
  !
  !念のため初期化
  !
  do i=1,MAXCOMPO
     call Tie_Constructor(tie(i))
  enddo

  ncompo=0
  incompo=-1
  !     中途半端にまぜるとややこしいので、
  !     基本的にすべてbinary I/Oとする。
  !     変更必要
  !
  !leave this for compatibility 2005-8-31
  !
  do
     read(TRAJIN,END=99, ERR=99) tag
     write(STDERR,*) tag
     if(tag.eq."@NLOG")then
        read(TRAJIN) nlog
        cycle
     endif
     call Box_BinaryLoader(box,TRAJIN,tag)
     call Cutoff_BinaryLoader(co,TRAJIN,tag)
     call Nose_BinaryLoader(nose,TRAJIN,tag)
     call Andersen_BinaryLoader(an,TRAJIN,tag)
     call Grid_BinaryLoader(g,trajin,tag)
     call Book_BinaryLoader(bk,trajin,tag)
     if(tag.eq."@VCNT")then
        read(TRAJIN) ivcnt,vcnt
        fVcnt=(ivcnt.gt.0)
        cycle
     endif
     if(tag.eq."@MDLP")then
        read(TRAJIN) nloop
        fLoop=.true.
        cycle
     endif
     if(tag.eq."@FIXC")then
        read(TRAJIN) isFixed
        cycle
     endif
     !
     !Morphingの初期座標をファイルから読みこむ場合、どの成分を使うかを指定。番号は1が最初の成分を表わす。
     !
     if(tag.eq."@DTPS")then
        call Time_ReadBinaryDTPS(t,TRAJIN)
        fDt=.true.
        cycle
     endif
     if(tag.eq."@DTSC")then
        call Time_ReadBinaryDTSC(t,TRAJIN)
        fDt=.true.
        cycle
     endif
     if(tag.eq."@TIE4")then
        ntie = ntie+1
        ftie = Tie_ReadBinaryTIE4(tie(ntie),TRAJIN)
        cycle
     endif
     if(tag.eq."@ASTP")then
        read(TRAJIN) nstep
        cycle
     endif
     !
     !i番目の成分とj番目の成分の相互作用に比例定数を掛ける。
     !
     if(tag.eq."@IVSC")then
        read(TRAJIN) i,j,value
        ivscale(i,j) = value
        ivscale(j,i) = value
        cycle
     endif
     !
     !Silicon専用; 三体力のみをscaleする。Volatile
     !
     if(tag.eq."@IVSI")then
        read(TRAJIN) ScaleSWSI
        cycle
     endif
     if(tag.eq."@ID08")then
        read(TRAJIN) lastid
        cycle
     endif
     if(tag.eq."@NCMP")then
        read(TRAJIN) incompo
        if(incompo.gt.MAXCOMPO)then
           call die( error_too_many_components, "Main 1" )
        endif
        cycle
     endif
     if(tag.eq.'@WTG2'.or.tag.eq.'@RMMC')then
        !
        !指定された成分数を読みこんだら入力をうちきる。
        !
        if(ncompo.eq.incompo)then
           cycle
        endif
        !
        !RigidMorphの場合は必ず固定座標扱いとする。
        !
        if ( tag .eq. '@RMMC' ) then
           isFixed = .true.
        endif
        if(lastid.eq."")then
           write(STDERR,*) "Error: Anonymous data"
           call die( 0, "Main 2" )
        endif
        ncompo = ncompo + 1
        if(ncompo.gt.MAXCOMPO)then
           write(STDERR,*) "Number of componets exceeded the&
                & limit."
           call die( 0, "Main 3" )
        endif
        !mol%id = id
        mol(ncompo)%id=lastid
        err = rigid_setup( rigid(ncompo), lastid, isFixed, mol(ncompo), si(ncompo) )
        if ( err .ne. error_none ) call die( err, "Main 4" )
        if ( mol(ncompo)%isRigid ) then
           if ( tag .eq. '@RMMC' ) then
#ifdef JARZYNSKI
              call RigidMorph_ReadBinaryNX4A( rigidmorph, rigid( ncompo ), mol( ncompo ), site, TRAJIN )
              morphcompo = ncompo
#else
#endif !JARZYNSKI
           else
              call Rigid_ReadBinaryWTG2(rigid(ncompo),mol(ncompo),site,TRAJIN)
           endif
        endif
        if ( mnomiac < mol(ncompo)%nmol )then
           mnomiac = mol(ncompo)%nmol
        endif
        if ( .not. isFixed ) dof = dof + Mol_DoF(mol(ncompo))
        lastid=""
        isFixed = .false.
        cycle
     endif
     if(tag.eq.'@APC5')then
        !
        !指定された成分数を読みこんだら入力をうちきる。
        !
        if(ncompo.eq.incompo)then
           cycle
        endif
        if(lastid.eq."")then
           write(STDERR,*) "Error: Anonymous data"
           call die( 0, "Main 5" )
        endif
        ncompo = ncompo + 1
        if(ncompo.gt.MAXCOMPO)then
           write(STDERR,*) "Number of componets exceeded the&
                & limit."
           call die( 0, "Main 6" )
        endif
        mol(ncompo)%id=lastid
        err = flex_setup( lastid, isFixed, mol(ncompo), si(ncompo), molmass )
        if ( err .eq. error_id_not_found ) then
           write( STDERR, * ) "Unknown molecule type: ",lastid
           call die( err, "Main 6.5" )
        endif
        if ( err .ne. error_none ) then
           call die( err, "Main 6.6" )
        endif
        call Flex_ReadBinaryAPC5(flex(ncompo),mol(ncompo),site,TRAJIN)
        call Flex_SetMass(flex(ncompo),mol(ncompo),molmass)
        if ( mnomiac < mol(ncompo)%nmol )then
           mnomiac = mol(ncompo)%nmol
        endif
        if ( .not. isFixed ) dof = dof + Mol_DoF(mol(ncompo))
        !dof = dof + Mol_DoF(mol(ncompo))
        lastid=""
        isFixed = .false.
        cycle
     endif
     if(tag.eq.'@NORT')then
        !only a dummy
        READ(TRAJIN)
        cycle
     endif
     if(tag.eq.'@VSTR')then
        !only a dummy
        READ(TRAJIN) av%interval,av%width
        fAverage=(av%interval.gt.0)
        cycle
     endif
  end do
99 continue
  
   !     Binaryを全部読みおえたあとで、変更する部分をstdinから読む。
   !     これがいちばんすっきりしている。
#ifdef USE_REPLICA
  CONFIGFILE=90
  call getarg( 1, configfilename )
  write(STDERR,*)"COMMON CONFIG FILE:", configfilename
  open( CONFIGFILE, FILE=configfilename )
#else
  CONFIGFILE=STDIN
#endif
  do
     read(CONFIGFILE,'(a5)',END=999) tag
     write(STDERR,*) tag
     if(tag.eq."@DTPS")then
        call Time_ReadDTPS(t,CONFIGFILE)
        fDt=.true.
        cycle
     endif
     if(tag.eq."@SNPO")then
        read(CONFIGFILE,*) snpo
        cycle
     endif
     if(tag.eq."@SNPV")then
        read(CONFIGFILE,*) snpv
        cycle
     endif
     !
     !圧力テンソルの出力
     !
     if(tag.eq."@PTSO")then
        read(CONFIGFILE,*) ptensor
        cycle
     endif
     if(tag.eq."@MDVW")then
        read(CONFIGFILE,*) mdvw
        cycle
     endif
     !
     !i番目の成分とj番目の成分の相互作用に比例定数を掛ける。
     !
     if(tag.eq."@IVSC")then
        read(CONFIGFILE,*) i,j,value
        ivscale(i,j) = value
        ivscale(j,i) = value
        cycle
     endif
     !
     !Silicon専用; 三体力のみをscaleする。Volatile
     !
     if(tag.eq."@IVSI")then
        read(CONFIGFILE,*) ScaleSWSI
        cycle
     endif
     if(tag.eq."@ID08")then
        read(CONFIGFILE,*) lastid
        cycle
     endif
     if(tag.eq."@NCMP")then
        read(CONFIGFILE,*) incompo
        if(incompo.gt.MAXCOMPO)then
           call die( error_too_many_components, "Main 8" )
        endif
        cycle
     endif
     if(tag.eq.'@NX4A' .or. tag.eq.'@AR3A' .or. tag.eq.'@FL3A' .or. tag.eq.'@RMMC' )then
        !
        !指定された成分数を読みこんだら入力をうちきる。
        !
        if( ncompo.eq.incompo )then
           cycle
        endif
        !
        !RigidMorphの場合は必ず固定座標扱いとする。
        !
        if ( tag .eq. '@RMMC' ) then
           isFixed = .true.
        endif
        if(lastid.eq."")then
           write(STDERR,*) "Error: Anonymous data"
           call die( 0, "Main 9" )
        endif
        ncompo = ncompo + 1
        if(ncompo.gt.MAXCOMPO)then
           write(STDERR,*) "Number of componets exceeded the&
                & limit."
           call die( 0, "Main 10" )
        endif
        !
        !Normal
        !
        mol(ncompo)%id=lastid
        err = rigid_setup( rigid(ncompo), lastid, isFixed, mol(ncompo), si(ncompo) )
        if ( err .eq. error_id_not_found ) then
           err = flex_setup( lastid, isFixed, mol(ncompo), si(ncompo), molmass )
           if ( err .eq. error_id_not_found ) then
              write( STDERR, * ) "Unknown molecule type: ",lastid
              call die( err, "Main 11.5" )
           endif
        endif
        if ( err .ne. error_none ) then
           call die( err, "Main 12" )
        endif
        if ( .not. mol(ncompo)%isRigid ) then
           mol(ncompo)%id=lastid
           if ( tag .eq. '@AR3A' ) then
              call Flex_ReadAR3A(flex(ncompo),mol(ncompo),site,.true.,CONFIGFILE)
           else if ( tag .eq. '@FL3A' ) then
              call Flex_ReadFL3A(flex(ncompo),mol(ncompo),site,.true.,CONFIGFILE)
           endif
           call Flex_SetMass(flex(ncompo),mol(ncompo),molmass)
#ifdef RATTLE
           if ( lastid.eq."SPC_E_F" ) then
              call Flex_Rattle_Register_SPCE( flex(ncompo), mol(ncompo), rattle )
           endif
#endif
        endif
        if ( mol(ncompo)%isRigid ) then
           if ( tag .eq. '@RMMC' ) then
#ifdef JARZYNSKI
              call RigidMorph_ReadNX4A( rigidmorph, rigid( ncompo ), mol( ncompo ), site, CONFIGFILE )
              morphcompo = ncompo
#else
#endif ! JARZYNSKI
           else
              call Rigid_ReadNX4A( rigid(ncompo),mol(ncompo),site,CONFIGFILE)
           endif
        endif
        if ( .not. isFixed ) dof = dof + Mol_DoF(mol(ncompo))
        !dof = dof + Mol_DoF(mol(ncompo))
        lastid=""
        isFixed = .false.
        cycle
     endif
     call Box_Loader(box,CONFIGFILE,tag)
     call Cutoff_Loader(co,CONFIGFILE,tag)
     call Nose_Loader(nose,CONFIGFILE,tag)
     call Andersen_Loader(an,CONFIGFILE,tag)
     call Grid_Loader(g,CONFIGFILE,tag)
     call Book_Loader(bk,CONFIGFILE,tag)
     if(tag.eq."@VCNT")then
        read(CONFIGFILE,*) ivcnt,vcnt
        fVcnt=(ivcnt.gt.0)
        cycle
     endif
     if(tag.eq."@MDLP")then
        read(CONFIGFILE,*) nloop
        fLoop=.true.
        cycle
     endif
     !
     !高速フォース処理。分子種ごとの関数を準備する必要がある。
     !FAST=2でmdgrapeを使うようにしようか。
     !
     if(tag.eq."@FAST")then
        read(CONFIGFILE,*) i
        ffast = ( i .ne. 0 )
        cycle
     endif
     !
     !trueなら次の成分の座標を固定
     !
     if(tag.eq."@FIXC")then
        read(CONFIGFILE,*) isFixed
        cycle
     endif
#ifdef JARZYNSKI
     if(tag.eq."@MRPH")then
        ! 最初と最後での、2構造の混合比を指定する。Volatile.
        read(CONFIGFILE,*) mixi, mixf
        cycle
     endif
     !
     !Jarzynskiのworkの計算を行うかどうか
     !
     if ( tag .eq. "@JZYN" ) then
        read(CONFIGFILE,*) ijar
        fjar = ( ijar.ne.0 )
        cycle
     endif
#endif !JARZYNSKI
     if(tag.eq."@NLOG")then
        read(CONFIGFILE,*) nlog
        cycle
     endif
     if(tag.eq.'@VSTR')then
        !only a dummy
        READ(CONFIGFILE,*) av%interval,av%width
        fAverage=(av%interval.gt.0)
        cycle
     endif
     if(tag.eq."@VSCA")then
        !体積のみスケーリング.体積比を与える。
        read(CONFIGFILE,*) vscale
        vscale=vscale**(1d0/3d0)
        cycle
     endif
     if(tag.eq."@TIEI")then
        !
        !全分子を重心に結びつける。
        !
        read(CONFIGFILE,*) tiea,tieb
        cycle
     endif
     if(tag.eq."@VSET")then
        read(CONFIGFILE,*) settemp,seed
        call new( random, seed )
        !@VSETは結果に出力しない。
        cycle
     endif
     !
     !定期的に、全エネルギーを計測し、それが所定の値になるように徐々に速度をスケールする。
     !
     if(tag.eq."@HSCA")then
        read(CONFIGFILE,*) hset_freq, hset_value, hset_ratio
        !@HSET is volatile.
        cycle
     endif
     if(tag.eq."@RSTO")then
        read(CONFIGFILE,*) stoprotation
        cycle
     endif
     if(tag.eq."@TSTO")then
        read(CONFIGFILE,*) stopdrift
        cycle
     endif
     if(tag.eq."@STOR")then
        read(CONFIGFILE,*) shifttoorigin
        cycle
     endif
     if(tag.eq."@JOIN")then
        !
        !分子対をバネで結ぶ。
        !
        if ( incompo .lt. 0 ) then
           write(STDERR,*) "(JOIN) Number of components must be specified in advance explicitly."
           call die( 0, "Main 14" )
        endif
        read(CONFIGFILE,*) i,j
        !
        !異種分子間のbindは必要があれば実装する
        !
        if ( i /= j ) call die( 0, "Main 15" )
        read(CONFIGFILE,*) press, stretch
        if ( .not. ( press .or. stretch ) ) then
           bind( i )%mode = BIND_NONE
        else
           bind( i )%mode    = BIND_ACTIVE
           bind( i )%press   = press
           bind( i )%stretch = stretch
           read(CONFIGFILE,*) bind( i )%a, bind( i )%b, bind( i )%balance
           bind( i )%a = bind( i )%a * J2I
           k = 0
           do
              k = k + 1
              read(CONFIGFILE,*) ii,jj
              if ( ii < 0 ) then
                 bind( i )%npair = k - 1
                 exit
              endif
              bind( i )%pair_i(k) = ii + 1
              bind( i )%pair_j(k) = jj + 1
           enddo
        endif
        cycle
     endif
     if(tag.eq."@JMON")then
        !
        !分子対をバネで結んだと仮定した場合のエネルギーを計算するが、全エネルギーには加算しない。
        !
        if ( incompo .lt. 0 ) then
           write(STDERR,*) "(JOIN) Number of components must be specified in advance explicitly."
           call die( 0, "Main 16" )
        endif
        read(CONFIGFILE,*) i,j
        !
        !異種分子間のbindは必要があれば実装する
        !
        if ( i /= j ) call die( 0, "Main 17" )
        read(CONFIGFILE,*) press, stretch
        if ( .not. ( press .or. stretch ) ) then
           jmon( i )%mode = BIND_NONE
        else
           jmon( i )%mode    = BIND_ACTIVE
           jmon( i )%press   = press
           jmon( i )%stretch = stretch
           read(CONFIGFILE,*) jmon( i )%a, jmon( i )%b, jmon( i )%balance
           jmon( i )%a = jmon( i )%a * J2I
           k = 0
           do
              k = k + 1
              read(CONFIGFILE,*) ii,jj
              if ( ii < 0 ) then
                 jmon( i )%npair = k - 1
                 exit
              endif
              jmon( i )%pair_i(k) = ii + 1
              jmon( i )%pair_j(k) = jj + 1
           enddo
        endif
        cycle
     endif
#ifdef SOLVATIONTEST
     if(tag.eq."@HBJO")then
        !
        !分子対をバネで結ぶ。
        !
        if ( incompo .lt. 0 ) then
           write(STDERR,*) "(JOIN) Number of components must be specified in advance explicitly."
           call die( 0, "Main 18" )
        endif
        read(CONFIGFILE,*) i,j
        !
        !異種分子間のbindは必要があれば実装する
        !
        if ( i /= j ) call die( 0, "Main 19" )

        read(CONFIGFILE,*) joint(i)%nmol, joint(j)%nmol
        joint( i )%mode = JOINTHB_ACTIVE
        read(CONFIGFILE,*) joint( i )%a, joint( i )%b, joint( i )%balance
        joint( i )%a = joint( i )%a * J2I
        allocate( joint(i)%members(joint(i)%nmol) )
        joint(i)%members(:) = 0
        allocate( joint(i)%bonds(joint(i)%nmol, joint(j)%nmol) )
        joint(i)%bonds(:,:) = 0
        allocate( joint(i)%partner( 2, joint(i)%nmol) )
        joint(i)%partner(:,:) = 0
        do
           read(CONFIGFILE,*) ii,jj
           if ( ii < 0 ) then
              exit
           endif
           !
           !iiからの2本目の出結合であれば、bondsの値を2とする。3本目以上が出現する可能性は考えなくてよい。
           !このやりかただと、i,jが異種分子のときに破綻するので注意。
           !
           if ( 10 <= joint( i )%members( ii+1 ) ) then
              joint( i )%bonds( ii+1, jj+1 ) = 2
              joint( i )%partner( 2, ii+1 ) = jj+1
           else
              joint( i )%bonds( ii+1, jj+1 ) = 1
              joint( i )%partner( 1, ii+1 ) = jj+1
           endif
           joint( i )%members( ii+1 )     = joint( i )%members( ii+1 ) + 10
           joint( i )%members( jj+1 )     = joint( i )%members( jj+1 ) + 1
        enddo
        cycle
     endif
#endif
#ifdef USE_REPLICA
     call MDReplica_Loader( mdreplica, CONFIGFILE, tag )
#endif
  enddo
999 continue
#ifdef USE_REPLICA
  !read common config from fort.10
  if ( 0 <= mdreplica%myrank ) then
     open( TRAJIN, file=mdreplica%intraj, form="UNFORMATTED", err=997 )
     do
        read(TRAJIN,END=997) tag
        write(STDERR,*) tag
        if(tag.eq."@NLOG")then
           read(TRAJIN) nlog
           cycle
        endif
        call Box_BinaryLoader(box,TRAJIN,tag)
        call Cutoff_BinaryLoader(co,TRAJIN,tag)
        call Nose_BinaryLoader(nose,TRAJIN,tag)
        call Andersen_BinaryLoader(an,TRAJIN,tag)
        call Grid_BinaryLoader(g,trajin,tag)
        call Book_BinaryLoader(bk,trajin,tag)
        if(tag.eq."@VCNT")then
           read(TRAJIN) ivcnt,vcnt
           fVcnt=(ivcnt.gt.0)
           cycle
        endif
        !
        !Removed for REPLICA
        !
        if(tag.eq."@MDLP")then
           !dummy
           read(TRAJIN) ierr
           cycle
        endif
        if(tag.eq."@FIXC")then
           read(TRAJIN) isFixed
           cycle
        endif
        !
        !Morphingの初期座標をファイルから読みこむ場合、どの成分を使うかを指定。番号は1が最初の成分を表わす。
        !
        if(tag.eq."@DTPS")then
           call Time_ReadBinaryDTPS(t,TRAJIN)
           fDt=.true.
           cycle
        endif
        if(tag.eq."@DTSC")then
           call Time_ReadBinaryDTSC(t,TRAJIN)
           fDt=.true.
           cycle
        endif
        if(tag.eq."@TIE4")then
           ntie = ntie+1
           ftie = Tie_ReadBinaryTIE4(tie(ntie),TRAJIN)
           cycle
        endif
        if(tag.eq."@ASTP")then
           read(TRAJIN) nstep
           cycle
        endif
        !
        !i番目の成分とj番目の成分の相互作用に比例定数を掛ける。
        !
        if(tag.eq."@IVSC")then
           read(TRAJIN) i,j,value
           ivscale(i,j) = value
           ivscale(j,i) = value
           cycle
        endif
        !
        !Silicon専用; 三体力のみをscaleする。Volatile
        !
        if(tag.eq."@IVSI")then
           read(TRAJIN) ScaleSWSI
           cycle
        endif
        if(tag.eq."@ID08")then
           read(TRAJIN) lastid
           cycle
        endif
        if(tag.eq."@NCMP")then
           read(TRAJIN) incompo
           if(incompo.gt.MAXCOMPO)then
              call die( error_too_many_components, "Main 1" )
           endif
           cycle
        endif
        if(tag.eq.'@WTG2'.or.tag.eq.'@RMMC')then
           !
           !指定された成分数を読みこんだら入力をうちきる。
           !
           if(ncompo.eq.incompo)then
              cycle
           endif
           !
           !RigidMorphの場合は必ず固定座標扱いとする。
           !
           if ( tag .eq. '@RMMC' ) then
              isFixed = .true.
           endif
           if(lastid.eq."")then
              write(STDERR,*) "Error: Anonymous data"
              call die( 0, "Main 2" )
           endif
           ncompo = ncompo + 1
           if(ncompo.gt.MAXCOMPO)then
              write(STDERR,*) "Number of componets exceeded the&
                   & limit."
              call die( 0, "Main 3" )
           endif
           !mol%id = id
           mol(ncompo)%id=lastid
           err = rigid_setup( rigid(ncompo), lastid, isFixed, mol(ncompo), si(ncompo) )
           if ( err .ne. error_none ) call die( err, "Main 4" )
           if ( mol(ncompo)%isRigid ) then
              if ( tag .eq. '@RMMC' ) then
#ifdef JARZYNSKI
                 call RigidMorph_ReadBinaryNX4A( rigidmorph, rigid( ncompo ), mol( ncompo ), site, TRAJIN )
                 morphcompo = ncompo
#else
#endif
              else
                 call Rigid_ReadBinaryWTG2(rigid(ncompo),mol(ncompo),site,TRAJIN)
              endif
           endif
           if ( mnomiac < mol(ncompo)%nmol )then
              mnomiac = mol(ncompo)%nmol
           endif
           if ( .not. isFixed ) dof = dof + Mol_DoF(mol(ncompo))
           lastid=""
           isFixed = .false.
           cycle
        endif
        if(tag.eq.'@APC5')then
           !
           !指定された成分数を読みこんだら入力をうちきる。
           !
           if(ncompo.eq.incompo)then
              cycle
           endif
           if(lastid.eq."")then
              write(STDERR,*) "Error: Anonymous data"
              call die( 0, "Main 5" )
           endif
           ncompo = ncompo + 1
           if(ncompo.gt.MAXCOMPO)then
              write(STDERR,*) "Number of componets exceeded the&
                   & limit."
              call die( 0, "Main 6" )
           endif
           mol(ncompo)%id=lastid
           err = flex_setup( lastid, isFixed, mol(ncompo), si(ncompo), molmass )
           if ( err .eq. error_id_not_found ) then
              write( STDERR, * ) "Unknown molecule type: ",lastid
              call die( err, "Main 6.5" )
           endif
           if ( err .ne. error_none ) then
              call die( err, "Main 6.6" )
           endif
           call Flex_ReadBinaryAPC5(flex(ncompo),mol(ncompo),site,TRAJIN)
           call Flex_SetMass(flex(ncompo),mol(ncompo),molmass)
           if ( mnomiac < mol(ncompo)%nmol )then
              mnomiac = mol(ncompo)%nmol
           endif
           if ( .not. isFixed ) dof = dof + Mol_DoF(mol(ncompo))
           !dof = dof + Mol_DoF(mol(ncompo))
           lastid=""
           isFixed = .false.
           cycle
        endif
        if(tag.eq.'@NORT')then
           !only a dummy
           READ(TRAJIN)
           cycle
        endif
        if(tag.eq.'@VSTR')then
           !only a dummy
           READ(TRAJIN) av%interval,av%width
           fAverage=(av%interval.gt.0)
           cycle
        endif
     end do
997  close( TRAJIN )

     open( CONFIGFILE, file=mdreplica%infile )
     do
        read( CONFIGFILE, *, END=998 ) tag
        call MDReplica_Loader( mdreplica, CONFIGFILE, tag )
        call Nose_Loader(nose,CONFIGFILE,tag)
        call Andersen_Loader(an,CONFIGFILE,tag)
        if(tag.eq."@VCNT")then
           read(CONFIGFILE,*) ivcnt,vcnt
           fVcnt=(ivcnt.gt.0)
           cycle
        endif
     enddo
998  close( CONFIGFILE )


  endif
  if ( .not. associated( mdreplica%rxgraph ) ) then
     !
     !rxgraphが指定されていない場合は、従来通り線形に交換するのではなく、完全グラフを与える。
     !
     call new( mdreplica%rxgraph )
     call graph_perfect( mdreplica%rxgraph, mdreplica%nprocs )
  end if
  call graph2pair( mdreplica%rxgraph, mdreplica%rxpair )
  !call pair_show( mdreplica%rxpair, STDERR )
#endif
  !
  !整合性のチェック
  !
  do compo=1,ncompo
     if ( tie(compo)%n > 0 .and. tie(compo)%n /= mol(compo)%nmol ) then
        write(STDERR,*) "FATAL(TIE):",compo,tie(compo)%n, mol(compo)%nmol
        call die( 0, "Main 20" )
     endif
  enddo
  !
  !tieaが設定されていたら、TIE4で読みこんだバネ定数を上書きしてよい。
  !
  if ( tiea > 0d0 ) then
     do compo=1,ncompo
        nmol = mol(compo)%nmol
        allocate( tmpa( nmol ) )
        allocate( tmpx( nmol ) )
        tmpa(:) = tiea
        tmpx(:) = 0d0
        call tie_initialize( tie( compo ), nmol, tmpx, tmpx, tmpx, tmpa, tieb )
        deallocate( tmpa )
        deallocate( tmpx )
     enddo
     ntie = ncompo
  endif
  !Hardcoding
  !体積変換平成14年4月25日(木)
  !
  !vscaleを読みこんだ場合は、体積をスケールするとともに、分子の重心位
  !置もスケールする。
  !
  !write(STDERR,*) "VSCALE: ",vscale
  call scale(box,vscale)
  do compo=1,ncompo
     nmol = mol(compo)%nmol
     if ( mol(compo)%isRigid ) then
        do i=1,nmol
           call scale(rigid(compo)%mol(i)%com,vscale)
        enddo
     else
        !
        !Scale position of each atom
        !(should be done in site_module)
        !
        do i=1, nmol * mol(compo)%nsite
           j = i + mol(compo)%offset
           site%x(j) = site%x(j) * vscale
           site%y(j) = site%y(j) * vscale
           site%z(j) = site%z(j) * vscale
        enddo
     endif
  enddo
  !write(6,*) t%dt
  !内部単位への変換(平成１２年４月３日(月)追加)
  !write(6,*) water%nmol,rwater%wx(1,1)
  do compo=1,ncompo
     if ( mol(compo)%isRigid ) then
        call rigid_toInternalUnit( rigid(compo), mol(compo), t )
     else
        call flex_toInternalUnit( flex(compo), mol(compo), t )
     endif
  enddo

  if ( nose%active .and. fVcnt ) then
     write(STDERR,*) "Do not specify two temperature-contorolling methods at a time."
     call die( 0, "Main 21" )
  endif
  !内部単位へ変換平成１２年４月６日(木)
  if(nose%active)nose%q=nose%q * t%dt
  if(an%mode.ne.noandersen)then
     if( .not. active( box ) )then
        write(STDERR,*)"Box is required to control pressure."
        call die( 0, "Main 22" )
     endif
     !dofで割ればたしかにシステムサイズに依存しなくはなるが、いいのか？
#ifdef TEST
     an%mass=an%mass/(t%dt**2)
     an%dump=an%dump
#else
     if ( an%compat_pcn23 ) then
        !
        !@PCN2/3互換モード。平成16年8月30日(月)廃止
        !
        write(STDERR,*) "PCN3 compat"
        an%mass=an%mass / (t%dt**2*dof )
     else
        !
        !分子あたりの質量が与えられるものとする。
        !
        
        an%mass=an%mass * dof / t%dt**2
     endif
#endif
     call Andersen_Initialize(an,box)
  endif
  if( 0 /= ptensor )then
     if( .not. active( box ) )then
        write(STDERR,*)"Box is required to calculate pressure tensor."
        call die( 0, "Main 23" )
     endif
  endif
  !
  !Grid分割を使用する時は帳簿も使う(ただしmarginは0とし、intervalを1に
  !するので、実質は使わないのと変わらない。
  !
  if(g%mode.ne.GRID_NONE)then
     if(.not.active(bk))then
        call Book_Initialize(bk, BOOK_FIX, 0d0, 1)
     endif
  endif
  if( active( box ) .and.co%mode.eq.CUTOFF_NONE)then
     co%mode=CUTOFF_POLY
     rc=dmin1(box%size%vec(1),box%size%vec(2),box%size%vec(3))*0.5d0
     write(STDERR,*) "Interaction truncated at ", rc, "A"
     call cutoff_initialize(co,rc-2d0,rc)
  endif
  !
  !いちおう変な設定をチェックしておく。
  !
  if(active(bk).and.g%mode.eq.GRID_NONE)then
     write(STDERR,*) "ERROR: GRID IS NECESSARY FOR BOOKING."
     call die( 0, "Main 24" )
  endif
  if(active(bk).and.co%mode.eq.CUTOFF_NONE)then
     write(STDERR,*) "ERROR: CUTOFF IS NECESSARY FOR BOOKING."
     call die( 0, "Main 25" )
  endif
#ifdef JARZYNSKI
  if ( fJar .and. rigidmorph%mode .ne. RigidMorph_Active ) then
     write(STDERR,*) "ERROR: RIGIDMORPH IS NECESSARY FOR FE ESTIMATION."
     call die( 0, "Main 26" )
  endif
  if ( fjar ) then
     jarwork(:) = 0d0
  endif
#endif !JARZYNSKI
  !
  !同種分子間の相互作用(分子数をnとするとn(n+1)/2対)
  !
  do compo=1,ncompo
     !
     !第2パラメータを省略した場合は自己相互作用とみなす。
     !
     !DECではoptionalの扱いがうまくいかないらしいので明示的に-1を与える。
     !
     call interaction_initialize(ww(compo),mol(compo)%nmol,-1)
     combi( compo ) = compo
     combj( compo ) = compo
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
        call interaction_initialize(ww(ncombi),mol(i)%nmol,mol(j)%nmol)
     enddo
  enddo
  do k=1,ncombi
     i = combi( k )
     j = combj( k )
     ww(k)%scale = ivscale(i,j)
     if ( mol(i)%id(1:8) .eq. "SWSILICO" .and. ww(k)%isomol ) then
        call SWSilicon_ScaleLambda( ScaleSWSI )
     endif
  enddo
  
  !write(*,*) "j",ncombi,ncompo
  if(active(bk))then
     if(bk%margin.le.0d0)then
        bk%mode = BOOK_FIX
     endif
     if(bk%mode .eq. BOOK_AUTOINTERVAL)then
        call die( 0, "Main 27" )
        !
        !まだ考えていない。帳簿法はほとんど使わないから無視しても構わ
        !ないのだが・・・!!!
        !
        !com0(1:water%nmol)=rwater%mol(1:water%nmol)%com
        !call Book_ResetInterval(bk,co%Out+bk%margin,water%nmol,com0)
     endif
  else
     !
     !系が小さい場合を想定;全分子同士が相互作用する(グリッド分割しない)
     !
     do j=1,ncombi
        call interaction_alltoall(ww(j))
     enddo
  endif
  if(g%mode.ne.GRID_NONE)then
     do i=1,ncompo
        call Address_Initialize(addr(i),mol(i)%nmol)
     enddo
  endif
  !call Rigid_initialize(rwater,water)
  if(settemp.ge.0d0)then
     !
     !settempが設定されていた場合は、強制的に温度を調節する。
     !
     do compo=1,ncompo
        if ( mol(compo)%isRigid ) then
           call Rigid_VReset2(rigid(compo),mol(compo),t,settemp, random)
        else
           call Flex_VReset2(flex(compo),mol(compo),t,settemp, random)
        endif
     enddo
  endif
  !
  !シミュレーションの条件をSTDERRに表示する。
  !
  write(STDERR,*) "Loop                                (steps):",nloop
  write(STDERR,*) "Interval                               (ps):",t%dt
  write(STDERR,*) "Statistics Logging           (steps, 0=off):",nlog
  write(STDERR,*) "Snapshot of coordinate       (steps, 0=off):",snpo
  write(STDERR,*) "Snapshot of coord and velocity (steps, 0=off):",snpv
  write(STDERR,*) "Output pressure tensor       (steps, 0=off):",ptensor
  write(STDERR,*) "Snapshot for MDView          (steps, 0=off):",mdvw
  write(STDERR,*) "Stop rotation at every n steps(steps, 0=off):",stoprotation
  write(STDERR,*) "Stop drift at every n steps  (steps, 0=off):",stopdrift
  write(STDERR,*) "Forced shift to origin       (steps, 0=off):",shifttoorigin
  write(STDERR,*) "Forced Temperature Renormalization         :",fVcnt
  if ( fVcnt ) then
  write(STDERR,*) "    Interval          (steps):",ivcnt
  write(STDERR,*) "    Temperature (K):",vcnt
  endif
  write(STDERR,*) "Bookkeeping method                         :",active(bk)
  if ( active(bk) )then
  write(STDERR,*) "    Interval of update              (steps):",bk%interval
  write(STDERR,*) "    Margin to avoid penetration         (A):",bk%margin
  endif
  if ( g%mode.eq.GRID_NONE ) then
  write(STDERR,*) "Cell division                              :",.false.
  endif
  if ( g%mode.eq.GRID_FIX ) then
  write(STDERR,*) "Cell division                              : Fixed"
  write(STDERR,*) "    Division                               :", g%ndivx, g%ndivy, g%ndivz
  endif
  if ( g%mode.eq.GRID_ADAPTIVE ) then
  write(STDERR,*) "Cell division                              : Adaptive"
  write(STDERR,*) "    Minimal cell size                   (A):", g%requested%vec(1), g%requested%vec(2), g%requested%vec(3)
  endif
  write(STDERR,*) "Nose's constant temperature method         :", nose%active
  if ( nose%active ) then
  write(STDERR,*) "    Temperature                         (K):", nose%temp
  write(STDERR,*) "    Mass                                   :", nose%q * 0021.8852721616497249235d0**2
  endif
  if ( an%mode.eq.noandersen ) then
  write(STDERR,*) "Andersen's constant pressure method        : none"
  endif
  if ( an%mode.eq.orthorhombic ) then
  write(STDERR,*) "Andersen's constant pressure method        : orthorhombic"
  endif
  if ( an%mode.eq.orz ) then
  write(STDERR,*) "Andersen's constant pressure method        : orz"
  endif
  if ( an%mode.ne.noandersen ) then
  write(STDERR,*) "    External pressure                 (atm):", an%pex
  write(STDERR,*) "    Mass                                   :", an%mass
  endif
  if ( co%mode.eq.CUTOFF_NONE ) then
  write(STDERR,*) "Interaction truncation                     : none"
  endif
  if ( co%mode.eq.CUTOFF_POLY ) then
  write(STDERR,*) "Interaction truncation                     : polynomial"
  write(STDERR,*) "    Interaction range                   (A):", co%out
  write(STDERR,*) "    Transition range                    (A):", co%out-co%in
  endif
  if ( .not. active( box ) ) then
  write(STDERR,*) "Simulation cell                            : none"
  endif
  if ( box%mode.eq.BOX_ORTHO ) then
  write(STDERR,*) "Simulation cell                            : orthorhombic"
  write(STDERR,*) "    Size                                (A):", box%size%vec(1),box%size%vec(2),box%size%vec(3)
  endif
  write(STDERR,*) "Fast algorithm                             :",ffast
  write(STDERR,*) "Volume scaling at initialization           :",vscale**3
  write(STDERR,*) "Random initial velocity                    :",0d0 <= settemp
  if ( 0d0 <= settemp ) then
  write(STDERR,*) "    Temperature                         (K):",settemp
  write(STDERR,*) "    Random Seed                            :",seed
  endif
  write(STDERR,*) "Tie molecules at specified positions       :",tiea.ne.0d0
  if ( tiea .ne. 0 ) then
  write(STDERR,*) "    Force constants                        :",tiea, tieb
  endif
  write(STDERR,*) "Total Energy Scaling                       :",hset_freq
  if ( hset_freq .ne. 0 ) then
  write(STDERR,*) "    Intended value                         :",hset_value
  write(STDERR,*) "    Feedback ratio                         :",hset_ratio
  endif

#ifdef JARZYNSKI
  write(STDERR,*) "Morphing (Rigid)   g                       :",rigidmorph%mode
  if ( rigidmorph%mode .eq. rigidmorph_active ) then
  write(STDERR,*) "    From, To                               :",mixi,mixf
  write(STDERR,*) "    Component                              :",morphcompo
  endif
#endif !JARZYNSKI
  write(STDERR,*) "Specified number of components             :",incompo
  do compo=1,ncompo
     write(STDERR,*) "Component ",compo
     if ( mol(compo)%isRigid ) then
        write(STDERR,*) "    Rigid molecule"
     else
        write(STDERR,*) "    Flexible/Atomic molecule"
     endif
     write(STDERR,*) "    ID                                     :", mol(compo)%id
     write(STDERR,*) "    Number of interaction sites            :", mol(compo)%nsite
     write(STDERR,*) "    Number of molecules                    :", mol(compo)%nmol
     write(STDERR,*) "    Degree of freedom                      :", mol(compo)%dof
     write(STDERR,*) "    Fixed to the initial coordinate        :", mol(compo)%isFixed
  enddo
  do k=1,ncombi
     i = combi( k )
     j = combj( k )
     write(STDERR,*) "Combination ", i,j
     if ( ww(k)%scale .ne. 1d0 ) then
        write(STDERR,*) "Interaction Scaling       : ",ww(k)%scale
     endif
     write(STDERR,*) "Bind molecular pairs                       :", bind( k )%mode
     if ( bind(k)%mode == BIND_ACTIVE ) then
        write(STDERR,*) "    Force constants                        :",bind(k)%a, bind(k)%b
        write(STDERR,*) "    Balancing distance                     :",bind(k)%balance
        write(STDERR,*) "    Number of pairs                        :",bind(k)%npair
        write(STDERR,*) "    Stretching force                       :",bind(k)%stretch
        write(STDERR,*) "    Pressing force                         :",bind(k)%press
     endif
     write(STDERR,*) "Simulate external field                       :", jmon( k )%mode
     if ( jmon(k)%mode == BIND_ACTIVE ) then
        write(STDERR,*) "    Force constants                        :",jmon(k)%a, jmon(k)%b
        write(STDERR,*) "    Balancing distance                     :",jmon(k)%balance
        write(STDERR,*) "    Number of pairs                        :",jmon(k)%npair
        write(STDERR,*) "    Stretching force                       :",jmon(k)%stretch
        write(STDERR,*) "    Pressing force                         :",jmon(k)%press
     endif
  enddo
  if ( ScaleSWSI .ne. 1d0 ) then
     write(STDERR,*) "Lambda of SWSilicon is scaled              :",scaleswsi
  endif
     
#ifdef USE_REPLICA
  do i=1,300
     node_in_charge(i) = i
  enddo
  !
  !reset timer
  !
  call fakempif_timer_delta( deltatime )
  !
  !dynamically updated interval between interceptions
  !
  ninnerloop = mdreplica%ninnerloop
  intercept = ninnerloop
  LOGFILE = 90
  open( LOGFILE, file=mdreplica%outfile )

  if ( nose%active ) then
     mdreplica%temp = nose%temp
  else if ( fVcnt ) then
     mdreplica%temp = vcnt
  else
     call die( 0, "Should be canonical." )
  endif

  call writetag( LOGFILE, "@MDRT" )
  write(LOGFILE,fmt0) mdreplica%temp
#else
  LOGFILE = STDOUT  
#endif
  !
  !主loop
  !
  do loop=1,nloop
#ifdef USE_REPLICA
     if ( loop == intercept ) then
        call fakempif_timer_delta( deltatime )
        write( STDERR,* ) deltatime, "ms for MD loops."
        !
        !計算時間の均衡化
        !
        speed = ninnerloop
        speed = speed / deltatime
        !write( STDERR,* ) mdreplica%myrank, "(", getnodeid(), ") took", deltatime, "ms for", ninnerloop, "loops", speed, sumdeltat, sumdeltar
        call fakempif_gather( 4, ninnerloop, processloop, ierr )
        call fakempif_gather( 4, deltatime , processtime, ierr )
        call fakempif_gather( 4, loop      , processaccm, ierr )

        if ( mdreplica%myrank == 0 ) then
           write( STDERR,* ) "Load balancing..."
           !
           !まず、全ノードの計算時間を求め、mdreplica%ninnerloopステップ分計算するのに必要な平均時間を求める。
           !
           avgtime = 0
           do i=1, mdreplica%nprocs
              avgtime = avgtime + mdreplica%ninnerloop * processtime(i) / processloop(i)
           enddo
           avgtime = avgtime / mdreplica%nprocs
           !
           !今度は、各ノードの計算速度をもとに、avgtime間で計算できるループ数を求める。
           !
           finallap = .false.
           do i=1, mdreplica%nprocs
              deltatime = processtime(i)
              if ( 100000 < deltatime ) deltatime = 100000 !safety valve for suspend
              processloop(i) = processloop(i) * avgtime / processtime(i)
              if ( nloop <= processaccm(i) + processloop(i) ) then
                 !
                 !final lap
                 !
                 finallap = .true.
              endif
           enddo
           write( STDERR,* ) "Done."
           if ( finallap ) then
              write( STDERR,* ) "This is the final lap."
              do i=1, mdreplica%nprocs
                 !processloop(i) = nloop
                 processloop(i) = 0
              enddo
           endif
        endif
        call fakempif_scatter( 4, ninnerloop, processloop, ierr )
        call fakempif_timer_delta( deltatime )
        write( STDERR,* ) deltatime, "ms for scheduling task."
        !write(STDERR,*) "Next", mdreplica%myrank, ninnerloop
        !
        !次にreplica exchange。交換するか否かの判断は親ノードが行い、結果をscatterする。
        !無駄な通信がある(親が保持していれば、gatherする必要がないもの)が、当面安全性を重視する。
        !
        !
        !convert ep from internal unit to kelvin
        !
        !epink = epsum * i2j * j2k
        !
        !Get temperature
        !
        if ( nose%active ) then
           mdreplica%temp = nose%temp
        else if ( fVcnt ) then
           mdreplica%temp = vcnt
        endif
        mdreplica%beta = 1d0 / mdreplica%temp
        beta0 = mdreplica%beta
        !
        !epsum --> rep_energy(myrank)
        !
        call fakempif_gather( 8, epsum, rep_energy, ierr )
        call fakempif_gather( 8, mdreplica%beta,   rep_beta,     ierr )
        if ( mdreplica%myrank == 0 ) then
           do i=1, mdreplica%nprocs
              !
              !iはノードの番号ではなく、OP領域の番号。
              !その領域を担っているnodeはnode_in_charge(i)
              !
              n0 = node_in_charge(i)
              write(STDERR,'(i10,i3,1x,f13.4,1x,f5.0,1x,99(g10.3,1x,g10.3,1x,g10.3,1x,g10.3,1x))') &
                   processaccm(n0), &
                   n0, &
                   rep_energy(n0), &
                   1d0/rep_beta(n0)
           enddo
           !
           !node 0でまとめて作業を行うのなら、Replica Exchangeを並列に実行する必要はない。
           if ( 1 < mdreplica%nprocs ) then
              do nexchange = 1, mdreplica%nexchange
                 k  = random_getnext( mdreplica%rand ) * mdreplica%rxpair%nbond + 1
                 j  = mdreplica%rxpair%a(k)
                 jj = mdreplica%rxpair%b(k)
                 !j,jjはノードの番号ではなく、OP領域の番号。
                 !その領域を担っているnodeはnode_in_charge(j)
                 accept = replica_exchange_trial( j, jj, node_in_charge, mdreplica%rand, rep_energy, rep_beta, 0 )
                 !
                 !When accepted, temp must be reflected to nose.
                 !
                 !call writetag( LOGFILE, "@MC2R" )
                 !n0 = node_in_charge(j)
                 !write(LOGFILE,fmt1) nexchange, rep_energy(n0), 1d0/rep_beta(n0) 
                 !call writetag( LOGFILE, "@MC2R" )
                 !n0 = node_in_charge(jj)
                 !write(LOGFILE,fmt1) nexchange, rep_energy(n0), 1d0/rep_beta(n0) 
              enddo
           endif
        endif
        if ( 1 < mdreplica%nprocs ) then
           call fakempif_scatter( 8, mdreplica%beta,   rep_beta,     ierr )
        endif
        !
        !Set temperature
        !
        if ( beta0 /= mdreplica%beta ) then
           mdreplica%temp = 1d0 / mdreplica%beta
           call writetag( LOGFILE, "@MDRT" )
           write(LOGFILE,fmt0) mdreplica%temp
           if ( nose%active ) then
              do compo=1,ncompo
                 if ( .not. mol(compo)%isFixed ) then
                    if ( mol(compo)%isRigid ) then
                       call Rigid_ScaleVelocity(rigid(compo),mol(compo),dsqrt( mdreplica%temp / nose%temp ))
                    else
                       call Flex_ScaleVelocity(flex(compo),mol(compo),dsqrt( mdreplica%temp / nose%temp ))
                    endif
                 endif
              enddo
              nose%temp = mdreplica%temp
              nose%zeta0 = 0d0
              nose%zeta1 = 0d0
              nose%zeta2 = 0d0
              nose%zeta3 = 0d0
              nose%zeta4 = 0d0
           else if ( fVcnt ) then
              vcnt = mdreplica%temp
           endif
        endif
        call fakempif_timer_delta( deltatime ) 
        write( STDERR,* ) deltatime, "ms for replica exchanges."
        if ( ninnerloop == 0 ) exit
        intercept = intercept + ninnerloop
     endif
#endif
     !
     !重心を座標原点にスライドさせる。
     !
     if ( 0 < shifttoorigin ) then
      if( mod(loop-1,shifttoorigin) == 0 ) then
        totalmoment%vec(:) = 0d0
        totalmass = 0d0
        do compo=1,ncompo
           if ( .not. mol(compo)%isFixed ) then
              if ( mol(compo)%isRigid ) then
                 call Rigid_GetOffset(rigid(compo),mol(compo),moment,mass)
              else
                 call Flex_GetOffset(flex(compo),mol(compo),site,moment,mass)
              endif
              totalmoment%vec(:) = totalmoment%vec(:) + moment%vec(:)
              totalmass          = totalmass + mass
           endif
        enddo
        totalmoment%vec(:) = -totalmoment%vec(:) / totalmass
        do compo=1,ncompo
           if ( .not. mol(compo)%isFixed ) then
              if ( mol(compo)%isRigid ) then
                 call Rigid_ShiftPosition(rigid(compo),mol(compo),totalmoment)
              else
                 !
                 !単にシフトさせるだけなので、flex()の情報は不要。
                 !
                 do i=1,mol(compo)%nmol*mol(compo)%nsite
                    j = i + mol(compo)%offset
                    site%x(j) = site%x(j) + totalmoment%vec(1)
                    site%y(j) = site%y(j) + totalmoment%vec(2)
                    site%z(j) = site%z(j) + totalmoment%vec(3)
                 enddo
              endif
           endif
        enddo
        write(STDERR,*) "OFFSET:",(totalmoment%vec(i),i=1,3)
     endif
     endif
     !
     !クラスタ全体の並進を止める
     !
     if ( 0 < stopdrift ) then
       if ( mod(loop-1,stopdrift) == 0 ) then
        totalmoment%vec(:) = 0d0
        totalmass          = 0d0
        do compo=1,ncompo
           if ( .not. mol(compo)%isFixed ) then
!!!
              if ( mol(compo)%isRigid ) then
                 call Rigid_GetMomentum( rigid(compo), mol(compo), moment, mass )
              else
                 call Flex_GetMomentum( flex(compo), mol(compo), moment, mass )
              endif
              totalmoment%vec(:) = totalmoment%vec(:) + moment%vec(:)
              totalmass = totalmass + mass
           endif
        enddo
        totalmoment%vec(:) = -totalmoment%vec(:) / totalmass
        do compo=1,ncompo
           if ( .not. mol(compo)%isFixed ) then
              if ( mol(compo)%isRigid ) then
                 call Rigid_AddMomenta(rigid(compo),mol(compo),totalmoment)
              else
                 call Flex_AddMomenta(flex(compo),mol(compo),totalmoment)
              endif
           endif
        enddo
        write(STDERR,*) "TRANS:",(totalmoment%vec(i),i=1,3)
     endif
     endif
     !
     !  クラスタ全体の回転を止める。(Bowl型の外場を与えれば、並進は止
     !  める必要がなくなる。
     !
     if ( 0 < stoprotation ) then
     if ( mod(loop-1,stoprotation) == 0 ) then
        !
        !系全体の回転の運動量を求める。
        !
        totalmoment%vec(:) = 0d0
        totaltensor(:,:) = 0d0
        do compo=1,ncompo
           if ( .not. mol(compo)%isFixed ) then
              if ( mol(compo)%isRigid ) then
                 call Rigid_GetAngularMomenta(rigid(compo),mol(compo),moment)
                 call Rigid_GetInertiaTensor(rigid(compo),mol(compo),tensor)
              else
                 call Flex_GetAngularMomenta(flex(compo),mol(compo),site,moment)
                 call Flex_GetInertiaTensor(flex(compo),mol(compo),site,tensor)
              endif
              totalmoment%vec(:) = totalmoment%vec(:) + moment%vec(:)
              totaltensor(:,:)   = totaltensor(:,:)   + tensor(:,:)
           endif
        enddo
        !
        !全体の角運動量はw = I**(-1) Lで求められる。
        !
        call InvSym33(totaltensor(1,1),totaltensor(2,2),totaltensor(3,3),&
             totaltensor(1,2),totaltensor(2,3),totaltensor(3,1))
        totaltensor(2,1) = totaltensor(1,2)
        totaltensor(3,2) = totaltensor(2,3)
        totaltensor(1,3) = totaltensor(3,1)
        do i=1,3
           omega%vec(i) = 0d0
           do j=1,3
              omega%vec(i) = omega%vec(i)-totaltensor(i,j)*totalmoment%vec(j)
           enddo
        enddo
        !
        !系全体の回転運動量を差しひく。混合物の場合、どう分配すればよ
        !いか。
        !
        !
        do compo=1,ncompo
           if ( .not. mol(compo)%isFixed ) then
              if ( mol(compo)%isRigid ) then
                 call Rigid_AddAngularVelocity(rigid(compo),mol(compo),omega)
              else
                 call Flex_AddAngularVelocity(flex(compo),mol(compo),site,omega)
              endif
           endif
        enddo
        write(STDERR,*) "ROT:",(omega%vec(i),i=1,3)
     endif
     endif
     !
     !予測子修正子法の予測子、あるいは分子位置の決定
     !
     do compo=1,ncompo
        if ( mol(compo)%isFixed ) then
#ifdef JARZYNSKI
           !
           !if rigidmorph is active,
           !
           if ( compo .eq. morphcompo ) then
              !
              !set new coordinate of the molecules
              !
              mix = dble(loop) / dble(nloop)
              mix = mixi + mix * ( mixf - mixi )
              call RigidMorph_Interpolate( rigidmorph, mix, rigid(compo) )
              !
              !mixを微小変化させた時の各相互作用点の変位を記録しておく。
              !自由エネルギー計算を速くできるかもしれない。平成16年8月19日(木)
              !
              call RigidMorph_SetTangentVectors0( rigidmorph, mix, rigid(compo), mol(compo), site )
#ifdef DEBUGTX
              !
              !alternative
              !
              do i=1,mol(compo)%nmol * mol(compo)%nsite
                 tx(i) = rigidmorph%tx(i)
                 ty(i) = rigidmorph%ty(i)
                 tz(i) = rigidmorph%tz(i)
              enddo
              call Rigid_SetSitePosition(rigid(compo),mol(compo),site)
              call RigidMorph_SetTangentVectors( rigidmorph, rigid(compo), mol(compo) )
              do i=1,mol(compo)%nmol * mol(compo)%nsite
                 write(STDERR,*) i,0,rigidmorph%tx(i),rigidmorph%ty(i),rigidmorph%tz(i)
                 write(STDERR,*) i,1,tx(i),ty(i),tz(i)
                 tx(i) = tx(i) - rigidmorph%tx(i)
                 ty(i) = ty(i) - rigidmorph%ty(i)
                 tz(i) = tz(i) - rigidmorph%tz(i)
                 if ( 0.001 .lt. dabs( tx(i)/rigidmorph%tx(i) ) .or. 0.001 .lt. dabs( ty(i)/rigidmorph%ty(i) ) .or. 0.001 .lt. dabs( tz(i)/rigidmorph%tz(i) ) )then
                    write(STDERR,*) "LARGE ERR", loop,i
                 endif
              enddo
#endif
           endif
#endif !JARZYNSKI
        else
           if ( mol(compo)%isRigid ) then
              call rigid_predict(rigid(compo),mol(compo))
              !write(STDERR,*) "check:", (rigid(compo)%mol(1)%com%vec(i),i=1,3)
              !write(STDERR,*) "check:", (rigid(compo)%mol(2)%com%vec(i),i=1,3)
              !write(STDERR,*) "check:", (rigid(compo)%mol(1)%intrax(i),i=1,3)
           else
              call flex_predict(flex(compo),mol(compo),site)
           endif
        endif
     enddo
     if(nose%active)call Nose_Predict(nose)
     if(an%mode.ne.noandersen) then
        call Andersen_Predict(an,box)
        do compo=1,ncompo
           if ( .not. mol(compo)%isFixed ) then
              !!!圧力一定と分子固定は両立するか
           endif
        enddo
     endif
     do compo=1,ncompo
        if ( mol(compo)%isRigid ) then
           !
           !剛体の位置と向きから、サイトの位置を決める。
           !
           call Rigid_SetSitePosition(rigid(compo),mol(compo),site)

        else
           !
           !座標に拘束を与える。
           !
#ifdef RATTLE
           call rattle_position( rattle, site )
           do i=1, mol(compo)%nmol * mol(compo)%nsite
              rsite = flex( compo )%xinfo(i)%rattlesite
              flex(compo)%xinfo(i)%dx(1) = flex(compo)%xinfo(i)%dx(1) + rattle%vx(rsite)
              flex(compo)%xinfo(i)%dy(1) = flex(compo)%xinfo(i)%dy(1) + rattle%vy(rsite)
              flex(compo)%xinfo(i)%dz(1) = flex(compo)%xinfo(i)%dz(1) + rattle%vz(rsite)
           enddo
#endif
           !
           !柔軟分子のサイト位置から、重心の位置を決める。
           !
           call flex_setcom( flex( compo ), mol( compo ), site )
        endif
     enddo
     !
     !重心座標だけを直列化
     !(重心を別配列で計算するのであれば、TIP4Pに第5のサイトを与える必要はないのでは？)
     !
     do compo=1, ncompo
        nmol = mol(compo)%nmol
        if ( mol(compo)%isRigid ) then
           com( 1:nmol, compo ) = rigid( compo )%mol(1:nmol)%com
        else
           do i=1,nmol
              k = mol( compo )%offset + i*mol( compo )%nsite
              com( i, compo )%vec(1) = site%x(k)
              com( i, compo )%vec(2) = site%y(k)
              com( i, compo )%vec(3) = site%z(k)
           enddo
        endif
     enddo
     if(active(bk))then
        bk%countdown=bk%countdown-1
        if(bk%countdown.le.0)then
           ! グリッドの分割のしかたを変える必要があるかどうか判断する。
           ! グリッド分割を更新
           call Grid_Update(g, box, bk%outer)
           do i=1,ncombi
              call Interaction_Reset(ww(i))
           enddo
           !
           !分子をセルに帰属させる。
           !
           do i=1,ncompo
              call Address_Assign(addr(i),box,g,mol(i)%nmol,com(1,i))
           enddo
           !
           !同種分子間相互作用対リストの作成
           !
           do i=1,ncompo
              call Grid_NeighborList(g,g%homo,ww(i),addr(i),addr(i))
           enddo
           !
           !異種分子間相互作用対リストの作成
           !
           do i=ncompo+1,ncombi
              call Grid_NeighborList(g,g%hetero,ww(i),&
                   addr(combi(i)),addr(combj(i)))
           enddo
           do k=1,ncombi
              i=combi(k)
              j=combj(k)
              call Interaction_Compress2(ww(k), active( box ),&
                   box ,co%out**2, com(1,i),com(1,j),.false.)
              !write(6,*) k, ww(k)%npair0, ww(k)%npair
           enddo
           if(bk%mode.eq.BOOK_AUTOINTERVAL)then
              !帳簿のマージン半径内に分子がどのくらい侵入してきたかを
              !チェックして、次に帳簿を更新すべきタイミングを決定する。
              call die( 0, "Main 28" )
              !call Book_AssumeInterval(bk,ww%nmol_i,com0&
              !& ,Book_BestMargin(bk,ww%nmol_i,com0))
           endif
           bk%countdown=bk%interval
        endif
     else
        do k=1,ncombi
           i = combi( k )
           j = combj( k )
           !
           !相互作用対リストを圧縮して、あとの力計算処理を高速化する
           !(ベクトル機では有効)
           !
           if ( co%mode.ne.CUTOFF_NONE) then
              call Interaction_Compress(ww(k),com(1,i),com(1,j),&
                   .false. ,CutOff_OuterLimit(co),box)
           else
              call Interaction_Compress(ww(k),com(1,i),com(1,j),&
                   .false. ,-1d0,box)
           endif
        enddo
     endif
     !
     !全分子のサイトに加わる力を初期化する。
     !
     call site_resetforce(site)

     !
     !ep: Potential energy
     !
#ifdef SOLVATIONTEST
     ejointhb = 0d0
#endif
     ep(1:ncombi) = 0d0
     do k=1,ncombi
        i = combi( k )
        j = combj( k )
        !
        !Smooth cutoffによる減衰係数を事前計算
        !
        !call CutOff_Dimmer(co, ww(k), site, mol(i), mol(j) )
        call CutOff_Dimmer2( co, ww(k), com(1,i), com(1,j) )
        !
        !妙案が思いつかないので、その場で分子種を見ながらどの相互作用
        !関数を使うかを判断する。将来別の方法を採用する可能性あり。
        !
        allocate( ep1(1:ww(k)%npair) )
        !
        !外場として作用する(固定された)分子同士の相互作用は計算しない。平成16年8月19日(木)これにより、外場サイトには、外場サイト以外から受ける力が蓄積される。これが、λ変化に対する外場の応答を計算するのに必要。
        !
        if ( .not. ( i .eq. morphcompo .and. j .eq. morphcompo ) ) then
           call force(ww(k),mol(i),mol(j),site,ep1,virialp(k)%mat33, ffast, si(i),si(j) )
        !else
        !   write( STDERR, * ) "Morph", i, j
        endif
        !
        !Siの3体相互作用。一応ここに展開して書く。
        !
        if ( mol(i)%id(1:8) .eq. "SWSILICO" .and. ww(k)%isomol ) then
           call ComposeTriplets( ww(k), SWSILICON_A_REAL, mol(i)%nmol, triplet )
           allocate( ep3( triplet%num) )
           call force_swsilicon_triplet(ww(k),triplet,mol(i),site,ep3,vrmat3)
           virialp(k)%mat33(:,:) = virialp(k)%mat33(:,:) + vrmat3(:,:)
           do tri=1, triplet%num
              ep(k) = ep(k) + ep3(tri)
           enddo
           deallocate( ep3 )
        endif
#ifdef SOLVATIONTEST
        if ( .not. ffast .and. joint(k)%mode .eq. JOINTHB_ACTIVE ) then
           allocate( ep1sol(1:ww(k)%npair) )
           call force_sol(ww(k),mol(i),mol(j),site,joint(k),ep1sol)
           do pair=1,ww(k)%npair
              ejointhb = ejointhb + ep1sol(pair)
           enddo
           deallocate( ep1sol )
        endif
#endif
        !write(STDERR,*) k,ww(k)%npair,ep1(1)
        if ( .not. ( i .eq. morphcompo .and. j .eq. morphcompo ) ) then
           do pair=1,ww(k)%npair
              ep(k) = ep(k) + ep1(pair)
           enddo
        endif
        deallocate( ep1 )
     enddo
     !
     !λ変化による外場ポテンシャルの変化量は、外場サイトに加わる力のベクトルと、λ変化による外場サイトの変位ベクトルの内積で計算される。See note 2004-08-20
     !
#ifdef JARZYNSKI
     if ( 0 < morphcompo ) then
        compo = morphcompo
        dudl = 0d0
        !write(STDERR,*) morphcompo, mol(compo)%nsite, mol(compo)%nmol
        do j=1,mol(compo)%nsite-1
           do i=1,mol(compo)%nmol
              k=(i-1)*mol(compo)%nsite+j
              dudl = dudl + rigidmorph%tx(k) * ( -site%fx( k + mol(compo)%offset ) )
              dudl = dudl + rigidmorph%ty(k) * ( -site%fy( k + mol(compo)%offset ) )
              dudl = dudl + rigidmorph%tz(k) * ( -site%fz( k + mol(compo)%offset ) )
           enddo
        enddo
        dudlsum = dudlsum + dudl
     endif
#endif !JARZYNSKI

     vpp%mat33(:,:) = 0d0
     epsum = 0d0
     do k=1,ncombi
        vpp%mat33(:,:) = vpp%mat33(:,:) + virialp(k)%mat33(:,:)
        epsum = epsum + ep(k)
     enddo
     !write(6,*) virialp(1)%mat33(1,1) + virialp(1)%mat33(2,2) + virialp(1)%mat33(3,3)  !NG
#ifdef DEBUG
     virialpsum = 0d0
     do i=1,3
        virialpsum = virialpsum + vpp%mat33(i,i)
     enddo
     !平成16年2月8日(日)virialsumの値が正確にvirialpsumの1/3になることを確認。
     write(LOGFILE,*) "VIRIAL",virialsum, virialpsum
#endif
     !
     !外場からの力の計算
     !
     exsum = 0d0
     do compo=1,ncompo
        if ( tie(compo)%mode == TIE_ACTIVE ) then
           write(STDERR,*) "Tie ",compo,"is active"
           call die( 0, "Main 29" )
           do i=1,mol(compo)%nmol
              fcom(i)%vec(:) = 0d0
           enddo
           call interaction_force_tie(tie(compo),box,mol(compo)%nmol,com(1,compo),fcom,etie)
           !
           !計算した力を重心に戻す。
           !
           do i=1,mol(compo)%nmol
              !
              !TIP4Pなら5番目のサイトが重心。
              !
              k = i*mol(compo)%nsite + mol(compo)%offset
              site%fx(k) = site%fx(k) + fcom(i)%vec(1)
              site%fy(k) = site%fy(k) + fcom(i)%vec(2)
              site%fz(k) = site%fz(k) + fcom(i)%vec(3)
           enddo
           exsum = exsum + etie
        endif
     enddo
     !
     !分子対形成力
     !
     emonsum = 0d0
     do j = 1, ncombi
        if ( bind( j )%mode == BIND_ACTIVE ) then
           compo = j
           do i=1,mol( compo )%nmol
              fcom(i)%vec(:) = 0d0
           enddo
           call interaction_force_bind( bind(j), box, com( 1, compo ), fcom, ebind )
           !
           !計算した力を重心に戻す。
           !
           do i=1,mol(compo)%nmol
              !
              !TIP4Pなら5番目のサイトが重心。
              !
              k = i*mol(compo)%nsite + mol(compo)%offset
              site%fx(k) = site%fx(k) + fcom(i)%vec(1)
              site%fy(k) = site%fy(k) + fcom(i)%vec(2)
              site%fz(k) = site%fz(k) + fcom(i)%vec(3)
           enddo
           exsum = exsum + ebind
        endif
        if ( jmon( j )%mode == BIND_ACTIVE ) then
           compo = j
           call interaction_force_bind( jmon(j), box, com( 1, compo ), fcom, ejmon )
           emonsum = emonsum + ejmon
        endif
     enddo
     !
     !Jarzynskiの方法で、workを計算する。
     !外場を微小変化させた場合のポテンシャルエネルギーの変化の積算
     !
#ifdef JARZYNSKI
     if ( fJar ) then
        do compo=1,ncompo
           if ( mol(compo)%isRigid .and. mol(compo)%isFixed .and. compo .eq. morphcompo ) then
              !
              !まず外場=クラスタの形をすこしだけ変化させる。
              !
              !
              !set new coordinate of the molecules
              !
              mix = dble(loop+1) / dble(nloop)
              mix = mixi + mix * ( mixf - mixi )
              call RigidMorph_Interpolate( rigidmorph, mix, rigid(compo) )
              !
              !剛体の位置と向きから、サイトの位置を決める。
              !
              call rigid_setsiteposition(rigid(compo),mol(compo),site)
              !
              !重心位置を更新する
              !
              nmol = mol(compo)%nmol
              com( 1:nmol, compo ) = rigid( compo )%mol(1:nmol)%com
           endif
        enddo
        do k=1,ncombi
           epjar = 0d0
           i = combi( k )
           j = combj( k )
           if ( i .eq. morphcompo .or. j .eq. morphcompo ) then
              call CutOff_Dimmer2( co, ww(k), com(1,i), com(1,j) )
              allocate( ep1(ww(k)%npair) )
              call pairpotential(ww(k),mol(i),mol(j),site,ep1,ffast,si(i),si(j) )
              do pair=1,ww(k)%npair
                 epjar = epjar + ep1(pair)
              enddo
              epjar = epjar - ep(k)
              deallocate( ep1 )
           endif
           jarwork( k )  = jarwork( k ) + epjar
           !
           !直接Uの変量を求めるこの方法は、正しく動いているようだが、コストが高いのでいずれ廃止する。
           !
        enddo
     endif
#endif !JARZYNSKI
           
     eksum = 0d0
     ektsum = 0d0
     ekrsum = 0d0
     do compo=1,ncompo
        if ( .not. mol(compo)%isFixed ) then
           if ( mol(compo)%isRigid ) then
              !
              !剛体分子単位で力を集計し、並進と回転の力に分類する。
              !
              call rigid_collectforce(rigid(compo),mol(compo),site)
              !
              !予測子修正子法の修正子
              !
              call rigid_correct2(rigid(compo),mol(compo),t,nose,an)
              !
              !運動エネルギーの算出
              !
              call rigid_ek(rigid(compo),mol(compo),t,ekt(compo),ekr(compo))
              ektsum = ektsum + ekt(compo)
              ekrsum = ekrsum + ekr(compo)
              !
              !圧力テンソルの運動量成分
              !
              call Rigid_PressureTensor( rigid(compo),mol(compo),t,virialk(compo)%mat33 )
           else
              !
              !CORSに加わった力は、カットオフの作用しているサイトに質量比で分配する。
              !
              call flex_distributeforce( flex( compo ), mol( compo ), site )
              !
              !予測子修正子法の修正子
              !
              call flex_correct2(flex(compo),mol(compo),site,t,nose,an)
              !
              !座標に拘束を与える。
              !
              !call rattle_position( rattle, site )
              !do i=1, mol(compo)%nmol * mol(compo)%nsite
              !   rsite = flex( compo )%xinfo(i)%rattlesite
              !   flex(compo)%xinfo(i)%dx(1) = flex(compo)%xinfo(i)%dx(1) + rattle%vx(rsite)
              !   flex(compo)%xinfo(i)%dy(1) = flex(compo)%xinfo(i)%dy(1) + rattle%vy(rsite)
              !   flex(compo)%xinfo(i)%dz(1) = flex(compo)%xinfo(i)%dz(1) + rattle%vz(rsite)
              !enddo
              !
              !速度に拘束を与える。
              !
              !do i=1, mol(compo)%nmol * mol(compo)%nsite
              !   rsite = flex( compo )%xinfo(i)%rattlesite
              !   rattle%vx(rsite) = flex(compo)%xinfo(i)%dx(1)
              !   rattle%vy(rsite) = flex(compo)%xinfo(i)%dy(1)
              !   rattle%vz(rsite) = flex(compo)%xinfo(i)%dz(1)
              !enddo
              !call rattle_velocity( rattle, site )
              !do i=1, mol(compo)%nmol * mol(compo)%nsite
              !   rsite = flex( compo )%xinfo(i)%rattlesite
              !   flex(compo)%xinfo(i)%dx(1) = rattle%vx(rsite)
              !   flex(compo)%xinfo(i)%dy(1) = rattle%vy(rsite)
              !   flex(compo)%xinfo(i)%dz(1) = rattle%vz(rsite)
              !enddo
              !
              !運動エネルギーの算出
              !
              call flex_ek(flex(compo),mol(compo),t,ekt(compo))
              ektsum = ektsum + ekt(compo)
              !ekrsum = ekrsum + ekr(compo)
              call Flex_PressureTensor( flex(compo),mol(compo),t,virialk(compo)%mat33 )
           endif
        endif
     enddo
     eksum = ektsum + ekrsum
     vpk%mat33(:,:)=0d0
     do compo=1,ncompo
        vpk%mat33(:,:) = vpk%mat33(:,:) + virialk(compo)%mat33(:,:)
     enddo
     vp%mat33(:,:) = vpk%mat33(:,:) + vpp%mat33(:,:)
     !write(6,*) vpk%mat33(1,1) + vpk%mat33(2,2) + vpk%mat33(3,3)  !OK
     !write(6,*) vpp%mat33(1,1) + vpp%mat33(2,2) + vpp%mat33(3,3)  !NG
     !
     !Volume x Pressure = vpの対角成分の和

     !
     volumepressure = vp%mat33(1,1)+vp%mat33(2,2)+vp%mat33(3,3)
#ifdef DEBUG
     !平成16年2月8日(日)前者の値が正確に後者の1/3になることを確認。
     write(LOGFILE,*) "VP", (2d0/3d0)*ektsum + virialsum, volumepressure
#endif
     !
     !運動エネルギーを温度に換算(固定されている分子は無視される。)
     !
     temperature = I2J*2d0*eksum/(rgas*dof)
     if( active( box ) )then
        !pressure = ((2d0/3d0)*ektsum+virialsum)*I2J*p2a/(Box_Volume(box)*NA*1d-30)
        !もしかしたら圧力が3倍になっているかも。
        pressure = (volumepressure/3d0)*I2J*p2a/(Box_Volume(box)*NA*ANG2M**3)
        !write(STDERR,*) volumepressure, volumepressure/3d0, pressure
        !stop
     else
        !
        !クラスタの圧力は計算できない。出力は体積圧力積で代用する。
        !何か示量変数で割りたいのだが、適当なものが思いあたらない。
        !
        pressure = (volumepressure/3d0)*I2J*p2a
     endif
     !
     !圧力テンソルP_{alpha beta}の瞬間値を出力。単位はPa。
     !See JPCB 104, 1332(2000)
     !
     if ( 0 < ptensor ) then
        if ( mod( nstep, ptensor ) .eq. 0 ) then
           !
           !
           !
           write(LOGFILE,'("@PTSR")')
           write(LOGFILE,'(9(e17.10,1x))') &
                (( vp%mat33(i,j)*I2J/(Box_Volume(box) * NA * ANG2M**3 ), i=1,3),j=1,3)
        endif
     endif
     !
     !拡張MD
     !
     if(nose%active) call Nose_Correct(nose,temperature)
     if(an%mode.ne.noandersen)then
        !pressure = 100.0
        call Andersen_Correct(an,box,pressure)
     endif
     !
     !温度を強制的に設定する。
     !
     if(fVcnt)then
        if(mod(nstep,ivcnt).eq.0)then
           do compo=1,ncompo
              if ( .not. mol(compo)%isFixed ) then
                 if ( mol(compo)%isRigid ) then
                    call Rigid_ScaleVelocity(rigid(compo),mol(compo),dsqrt(vcnt/temperature))
                 else
                    call Flex_ScaleVelocity(flex(compo),mol(compo),dsqrt(vcnt/temperature))
                 endif
              endif
           enddo
        endif
     endif
     !
     !全エネルギーを再設定する。
     !
     if( hset_freq .ne. 0 )then
        if( mod( nstep, hset_freq ) .eq. 0 )then
           hset_ep = 0d0
           do j=1,ncombi
              hset_ep = hset_ep + ep(j)
           enddo
           hset_ek = 0d0
           do j=1,ncompo
              hset_ek = hset_ek + ( ekt(j) + ekr(j) )
           enddo
           ! 全エネルギーに加えるべき量。
           hset_param = (hset_value * J2I * KX2X - ( hset_ep + hset_ek ) ) * hset_ratio
           ! 運動エネルギーのスケーリング
           hset_param = ( hset_ek + hset_param ) / hset_ek
           if ( hset_param .lt. 0d0 ) call die( 0, "Negative temperature." )
           do compo=1,ncompo
              if ( .not. mol(compo)%isFixed ) then
                 if ( mol(compo)%isRigid ) then
                    call Rigid_ScaleVelocity(rigid(compo),mol(compo),dsqrt( hset_param ))
                 else
                    call Flex_ScaleVelocity(flex(compo),mol(compo),dsqrt( hset_param ))
                 endif
              endif
           enddo
        endif
     endif


#ifdef NOSEDBG
     write(6,*) (nose%zetasum-nose%zeta0*0.5d0)*t%dt,dlog(nose%s0),nose%q/t%dt
#endif
     !     write(6,*) ep,ekt,ekr,ep+ekt+ekr,temperature
     if( 0 < snpo )then
        if(mod(nstep,snpo).eq.0)then
           write(LOGFILE,'("@ASTP")')
           write(LOGFILE,*)nstep
           call save(box,LOGFILE, BOX_BOX3)
           write(LOGFILE,'("@NCMP")')
           write(LOGFILE,*)ncompo
           do k=1,ncombi
              i = combi( k )
              j = combj( k )
              if ( ww(k)%scale .ne. 1d0 ) then
                 write(LOGFILE,'("@IVSC")')
                 write(LOGFILE,*) i,j,ww(k)%scale
              endif
           enddo
           do compo=1,ncompo
              if ( mol(compo)%isRigid ) then
                 call save( rigid(compo), mol(compo), LOGFILE, RIGID_NX4A)
              else
                 call save( flex(compo), mol(compo), site, LOGFILE, FLEX_AR3A)
              endif
           enddo
        endif
     endif
     if( 0 < snpv )then
        if(mod(nstep,snpv).eq.0)then
           write(LOGFILE,'("@ASTP")')
           write(LOGFILE,*)nstep
           call save(box,LOGFILE, BOX_BOX3)
           write(LOGFILE,'("@NCMP")')
           write(LOGFILE,*)ncompo
           do k=1,ncombi
              i = combi( k )
              j = combj( k )
              if ( ww(k)%scale .ne. 1d0 ) then
                 write(LOGFILE,'("@IVSC")')
                 write(LOGFILE,*) i,j,ww(k)%scale
              endif
           enddo
           do compo=1,ncompo
              if ( mol(compo)%isRigid ) then
                 call Rigid_Save( rigid(compo), mol(compo), LOGFILE, RIGID_WTG3, t)
              else
                 !call save( flex(compo), mol(compo), site, LOGFILE, FLEX_AR3A)
              endif
           enddo
        endif
     endif
     if( 0 < mdvw )then
        if(mod(nstep,mdvw).eq.0)then
           write(LOGFILE,'("@ASTP")')
           write(LOGFILE,*)nstep
           !MDViewは直方体セルをサポートしていないので、別個に情報を保
           !存しておく必要がある(と思う)平成15年10月27日(月)
           call save(box,LOGFILE, BOX_BOX3)
           write(LOGFILE,'("@MDVW")')
           call save(box,LOGFILE, BOX_MDVW)
           mdvwsite = 0
           do compo=1,ncompo
              mdvwsite = mdvwsite + Mol_CountSiteMDVW( mol(compo) )
           enddo
           write(LOGFILE,*) mdvwsite
           do compo=1,ncompo
              if ( mol(compo)%isRigid ) then
                 call save( rigid(compo), mol(compo), LOGFILE, RIGID_MDVW )
              else
                 call save( flex(compo), mol(compo), site, LOGFILE, FLEX_MDVW)
              endif
           enddo
        endif
     endif
     !不完全。運動する分子しか出力されない。
     if(fAverage)then
        !if(mod(nstep,av%interval).eq.0)then
        !   do i=1,water%nmol
        !      qa0(i)=rwater%mol(i)%qa(1)
        !      qb0(i)=rwater%mol(i)%qb(1)
        !      qc0(i)=rwater%mol(i)%qc(1)
        !      qd0(i)=rwater%mol(i)%qd(1)
        !   enddo
        !   call Average_Rigid(av,water%nmol,comx0,comy0,comz0,qa0,qb0,qc0,qd0)
        !   if(av%count.eq.0)then
        !      write(vstrout,21)
        !      write(vstrout,22)
        !      write(vstrout,23)
        !      write(vstrout,*)nstep
        !      call Box_Save(box,vstrout)
        !      call Average_RigidWriteNX4A(av,vstrout,water%nmol)
        !   endif
        !endif
     endif
     if(nlog.ne.0)then
        if(mod(nstep,nlog).eq.0)then
           nose_ek=(nose%zeta0**2*nose%q*t%Dt*dof*0.5d0)*RGAS
           !dz/dt=zeta0だから、zeta0を台形則で積分する。(シンプソンに
           !してもそんなに手間は変わらないけど)
           nose_ep=(dof*nose%temp*(nose%zetasum-0.5d0*nose%zeta0)*t%Dt)*RGAS
           ! N/m^2 * m^3 * NA 
           andersen_ep = an%pex   * an%v0         * a2p * NA * ANG2M**3
           ! A^6 * atm/A^3 * 0.5 * a2p * NA * 1d-30
           ! = m^3 Pa * NA
           andersen_ek = an%v1**2 * an%mass*0.5d0 * a2p * NA * ANG2M**3
           k = 0
           do j=1,ncompo
              k = k + mol(j)%nmol
           enddo
           !
           !混合物のエネルギーを、全分子数で割って表示するのはおかしい
           !ので、総和を表示する。
           !

           !
           !可変長データを一行に出力するための工夫。
           !
           call empty( output )
           call push( output, temperature )      !2
           call push( output, pressure )         !3
           call push( output, Box_Volume(box) )  !4
           do j=1,ncombi
              call push( output, ( ep(j) * I2J * X2KX ) )
           enddo
           do j=1,ncompo
              call push( output, ( ekt(j) + ekr(j) ) * I2J * X2KX )
           enddo
           call push( output, epsum*I2J/dble(k) * X2KX )  !7
           call push( output, eksum*I2J/dble(k) * X2KX )  !8
           call push( output, nose_ep * X2KX )              !9
           call push( output, nose_ek * X2KX )              !10
           call push( output, andersen_ep * X2KX )          !11
           call push( output, andersen_ek * X2KX )          !12
           call push( output, (epsum + eksum) * I2J + nose_ep + nose_ek + andersen_ep + andersen_ek )
           call push( output, exsum*I2J )
           call push( output, (epsum + eksum + exsum) * I2J + nose_ep + nose_ek +andersen_ep + andersen_ek )
           call push( output, emonsum*I2J )
#ifdef SOLVATIONTEST
           call push( output, (epsum + eksum + exsum + ejointhb) * I2J + nose_ep + nose_ek + andersen_ep + andersen_ek )
#endif
#ifdef JARZYNSKI
           if ( fJar ) then
              call push( output, epJar * I2J )
              do j=1, ncombi
                 call push( output, Jarwork(j) * I2J )
              enddo
           endif
           if ( morphcompo .ne. 0 ) then
              mix = ( mixf - mixi ) / dble( nloop )
              call push( output, dudl * I2J )          !これを平均すれば<du/dl>_lが得られる。
              call push( output, dudlsum * I2J * mix ) !mixi!=mixfの場合の検算用。
           endif
#endif !JARZYNSKI
           call writetag( LOGFILE, '@LOGD' )
           write(LOGFILE,'(i10,30(1x,E24.17))') nstep, ( output%value(j), j=1, output%depth )
        endif
     endif
#ifdef MEMUSAGE
     !
     !配列が溢れないか、逐一チェックする。
     !
     call update( mon(1), g%ndivx )
     call update( mon(2), g%ndivy )
     call update( mon(3), g%ndivz )
     do j=1,ncompo
        call update( mon(4), addr(j)%nmax )
     enddo
     call update( mon(5), g%homo%n )
     call update( mon(5), g%hetero%n )
     call update( mon(9), g%ndivx*g%ndivy*g%ndivz )
     do j=1,ncompo
        call update( mon(6), mol(j)%nmol )
     enddo
     call update( mon(7), site%nsite )
     do j=1,ncombi
        call update( mon(8), ww(j)%npair0 )
     enddo
#ifdef VPOPTIMIZE
     do j=1,ncombi
        call update( mon(10), ww(j)%maxpartner_i )
        call update( mon(10), ww(j)%maxpartner_j )
     enddo
#endif
     if(mod(nstep,MEMUSAGE).eq.0)then
        do i=1,10
           call showstatus( mon(i) )
        enddo
     endif
#endif
     nstep=nstep+1
  enddo
  !     結果の出力
#ifdef USE_REPLICA
  open( TRAJOUT, file=mdreplica%outtraj, action='WRITE', form='UNFORMATTED' )
#endif
  call save_binary(co,TRAJOUT)
  call save_binary(box,TRAJOUT)
  call save_binary(bk,TRAJOUT)
  call save_binary(g,TRAJOUT)
  if(fVcnt)then
     write(TRAJOUT) "@VCNT"
     write(TRAJOUT) ivcnt,vcnt
  endif
  write(TRAJOUT)"@NLOG"
  write(TRAJOUT) nlog
  call Time_WriteBinaryDTPS(t,TRAJOUT)
  write(TRAJOUT)"@ASTP"
  write(TRAJOUT)nstep
  if(nose%active)then
     !内部単位へ変換平成１２年４月６日(木)
     nose%q=nose%q / t%dt
     call Nose_WriteBinaryNOS2(nose,TRAJOUT)
  endif
  if(an%mode.ne.noandersen)then
     an%mass=an%mass*(t%dt**2) / dof
     call Andersen_WriteBinaryPCN5(an,TRAJOUT)
  endif
  write(TRAJOUT) "@TEST"
  write(TRAJOUT) "@MDLP"
  write(TRAJOUT) nloop
  !互換単位系にもどす
  !面倒をさけるため、できれば内部単位をぜんぶ変更した方がよい。
  do compo=1,ncompo
     if ( mol(compo)%isRigid ) then
        call rigid_toExternalUnit( rigid(compo), mol(compo), t )
     else
        call flex_toExternalUnit( flex(compo), mol(compo), t )
     endif
  enddo
  !
  !Number of components
  !
  write(TRAJOUT) "@NCMP"
  write(TRAJOUT) ncompo
  do compo=1,ncompo
     call save_binary(tie(compo), TRAJOUT)
  enddo
  !
  !morphの場合は、初座標と終座標と混合比を出力すべきでは？
  !平成16年10月8日(金)
  !
  do compo=1,ncompo
     if ( mol(compo)%isRigid ) then
#ifdef JARZYNSKI
        if ( rigidmorph%mode .eq. rigidmorph_active .and. morphcompo .eq. compo ) then
           call save_binary(rigidmorph, mol( compo ), TRAJOUT )
        else
           call save_binary(rigid(compo), mol(compo), TRAJOUT, RIGID_WTG2)
        endif
#else
        call save_binary(rigid(compo), mol(compo), TRAJOUT, RIGID_WTG2)
#endif !JARZYNSKI
     else
        call save_binary(flex(compo), mol(compo), site, TRAJOUT, FLEX_APC5)
     endif
  enddo
  do k=1,ncombi
     i = combi( k )
     j = combj( k )
     if ( ww(k)%scale .ne. 1d0 ) then
        write(TRAJOUT) '@IVSC'
        write(TRAJOUT) i,j,ww(k)%scale
     endif
  enddo
  !     Done
  !output for book and grid
  !
  !本当は乱数を保存すべき。
  !
#ifdef FAKEMPI
  call fakempi_finalize( ierr )
#endif
end program main
