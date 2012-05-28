module replica
  implicit none
contains


  !  logical function replica_exchange_trial( j, jj, node_in_charge, rand, numbrella, rep_energy, rep_umbrella, rep_op, rep_beta )
  logical function replica_exchange_trial( j, jj, node_in_charge, rand, rep_energy, rep_beta, numbrella, rep_umbrella, rep_op )
    use physconst_module
    use random_module
    use umbrella
    type(pRandom),pointer :: rand
    integer, intent(IN) :: j,jj, numbrella
    integer, intent(INOUT) :: node_in_charge(*)
    real(kind=8), intent(IN) :: rep_energy(*)
    type(Umbrella_parabolic), intent(INOUT), optional :: rep_umbrella(*)
    real(kind=8), intent(IN), optional :: rep_op(*)
    real(kind=8), intent(INOUT) :: rep_beta(*)
    
    integer :: n0, n00,n1,n11,k
    real(kind=8) :: acceptance0, acceptance1, a0, a1, ratio, dummy
    logical :: accept
    type(Umbrella_parabolic) :: dumb
    !
    !現在の状態の実現確率
    !
    n0 = node_in_charge(j)
    n00 = (n0-1) * numbrella
    n1 = node_in_charge(jj)
    n11 = (n1-1) * numbrella
    
    a0 = rep_energy(n0)
    a1 = rep_energy(n1)
    do k=1, numbrella
       a0 = a0 + umbrella_potential( rep_umbrella(n00+k), rep_op(n00+k) )
       a1 = a1 + umbrella_potential( rep_umbrella(n11+k), rep_op(n11+k) )
    enddo
    acceptance0 = - ( a0 * rep_beta(n0) + a1 * rep_beta(n1) ) * I2J * J2K
    !
    !新しい状態の実現確率
    !
    a0 = rep_energy(n0)
    a1 = rep_energy(n1)
    do k=1, numbrella
       a0 = a0 + umbrella_potential( rep_umbrella(n11+k), rep_op(n00+k) )
       a1 = a1 + umbrella_potential( rep_umbrella(n00+k), rep_op(n11+k) )
    enddo
    acceptance1 = - ( a0 * rep_beta(n1) + a1 * rep_beta(n0) ) * I2J * J2K
    
    accept = .false.
    if ( acceptance0 < acceptance1 ) then
       accept = .true.
    else
       ratio = dexp( acceptance1 - acceptance0 )
       if ( random_getnext( rand ) < ratio ) then
          accept = .true.
          !write(STDERR,*) "Ratio",ratio,deltae * I2J * J2K,"dE"
       endif
    endif
    write(STDERR,*) n0,n1, acceptance1 - acceptance0, accept
    if ( accept ) then
       do k=1, numbrella
          dumb = rep_umbrella(n00+k)
          rep_umbrella(n00+k) = rep_umbrella(n11+k)
          rep_umbrella(n11+k) = dumb
       enddo
       dummy = rep_beta(n0)
       rep_beta(n0) = rep_beta(n1)
       rep_beta(n1) = dummy
       node_in_charge(j)   = n1
       node_in_charge(jj)  = n0
    endif
    replica_exchange_trial = accept
  end function replica_exchange_trial

end module replica
