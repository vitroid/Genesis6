#ifcはreal*8の定数値にrealの値を与えると、realの有効精度以下のぶぶんにひどい雑音が入る。
BINTARGET       = mc2 traj2 morph surface Mixture MixtureQuench # mctest NeatRigid2 MeOHWater MeOHWaterSmall Mixture MixtureSmall MixtureQuench NeatRigidQuench traj2
#BINTARGET       = traj2 morph Mixture MixtureQuench

#
#デバッグを行う場合
#-DMEMUSAGE=100 で、宣言した配列よりも大きいデータを扱おうとしていない
#かどうかを100ステップごとにチェックする。
#-DVERBOSEで頻繁にチェック文字列を表示する。
#
COMMONCPPFLAGS	= #-DMEMUSAGE=1 #-DVERBOSE #-DSOLVATIONTEST #-DVERBOSE##-DVPOPTIMIZE##-DEWALD#-Wp,-DVERBOSE

#include Makefile.osf1
#include Makefile.ifc
include Makefile.gfortran


#
#SAMPLES OF COMMONMODULE ARE PUT IN ./UNUSED/
#

include Makefile.rules

depend:: .depend
.depend: *.F90
	perl f90makedepend *.F90 > .depend

include .depend

