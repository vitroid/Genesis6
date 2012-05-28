***********************************************************************
*          mtrnd (Mersenne Twister random numbers generator)
*
***********************************************************************
* (Reference) M.Matsumoto, & T.Nishimura, ACM Transactions on Modeling
* and Computer Simulation Vol.8, No.1, Pages 3-30 (1998)
*
***********************************************************************
* A C-program for MT19937, with initialization improved 2002/1/26.
* Coded by Takuji Nishimura and Makoto Matsumoto.
*
* Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
*
*   1. Redistributions of source code must retain the above copyright
*      notice, this list of conditions and the following disclaimer.
*
*   2. Redistributions in binary form must reproduce the above copyright
*      notice, this list of conditions and the following disclaimer in the
*      documentation and/or other materials provided with the distribution.
*
*   3. The names of its contributors may not be used to endorse or promote
*      products derived from this software without specific prior written
*      permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
* A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
* EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
* PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
* LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
************************************************************************
* Fortran translation and modification by BABA Akinori.  Apr. 25, 2002.
*
************************************************************************
      subroutine mtrnu
      parameter(N = 624, M = 397)
      parameter(NVEC = N-M-1, NLOOP = 50, NMAX = NVEC*NLOOP+N)
      parameter(MATA = -1727483681)
      parameter(MASKTB = -1658038656)
      parameter(MASKTC = -272236544)
      integer mt(NMAX)
      common /mtrnc/ mti, mt
      save   /mtrnc/
      integer mag01(0:1)
      MASKL = 2147483647
      MASKU = ieor(MASKL,-1)
      mag01(0) = 0
      mag01(1) = MATA
      do 1000 i = 1, N
        mt(i) = mt(NMAX-N+i)
 1000 continue
      do 2000 iloop = 1, NLOOP
        k = NVEC*(iloop-1)
        do 2100 i = 1, NVEC
          ki = k+i
          iy = ior(iand(mt(ki),MASKU),iand(mt(ki+1),MASKL))
          mt(ki+N) = ieor(ieor(mt(ki+M),ishft(iy,-1)),mag01(iand(iy,1)))
 2100   continue
 2000 continue
      do 3000 i = 1, NMAX-N
        iy = mt(i)
        iy = ieor(iy,ishft(iy,-11))
        iy = ieor(iy,iand(ishft(iy,7),MASKTB))
        iy = ieor(iy,iand(ishft(iy,15),MASKTC))
        iy = ieor(iy,ishft(iy,-18))
        mt(i) = iy
 3000 continue
      mti = 0
      return
      end
*
      subroutine mtrng0(iseed)
      integer iseed
      parameter(N = 624, M = 397)
      parameter(NVEC = N-M-1, NLOOP = 50, NMAX = NVEC*NLOOP+N)
      integer mt(NMAX)
      common /mtrnc/ mti, mt
      save   /mtrnc/
      do 1000 i = 1, NMAX-N
        mt(i) = 0
 1000 continue
      mt(NMAX-N+1) = iseed
      do 1100 i = 2, N
        iv = mt(NMAX-N+i-1)
        ix = ieor(iv,ishft(iv,-30))
        ih = ishft(ix,-16)
        il = iand(ix,65535)
        mt(NMAX-N+i) = ishft(iand(il*27655,65535),16)
     *    +ishft(iand(ih*35173,65535),16)+il*35173+i-1
 1100 continue
      return
      end
*
      subroutine mtrngi(iseed)
      integer iseed
      parameter(N = 624, M = 397)
      parameter(NVEC = N-M-1, NLOOP = 50, NMAX = NVEC*NLOOP+N)
      integer mt(NMAX)
      common /mtrnc/ mti, mt
      save   /mtrnc/
      call mtrng0(iseed)
      call mtrnu
      mti = N
      return
      end
*
      subroutine mtrngv(larray, narray)
      integer narray
      integer larray(narray)
      parameter(N = 624, M = 397)
      parameter(NVEC = N-M-1, NLOOP = 50, NMAX = NVEC*NLOOP+N)
      integer mt(NMAX)
      common /mtrnc/ mti, mt
      save   /mtrnc/
      call mtrng0(19650218)
      i = 2
      j = 1
      k = N
      if (k .lt. narray) k = narray
      do 1000 ik = 1, k
        iv = mt(NMAX-N+i-1)
        ix = ieor(iv,ishft(iv,-30))
        ih = ishft(ix,-16)
        il = iand(ix,65535)
        mt(NMAX-N+i) = ieor(ishft(iand(il*25,65535),16)+
     *    ishft(iand(ih*26125,65535),16)+il*26125,mt(NMAX-N+i))
     *    +larray(j)+j-1
        i = i+1
        j = j+1
        if (i .gt. N) then
          mt(NMAX-N+1) = mt(NMAX)
          i = 2
        endif
        if (j .gt. narray) j = 1
 1000 continue
      do 1100 ik = 1, N-1
        iv = mt(NMAX-N+i-1)
        ix = ieor(iv,ishft(iv,-30))
        ih = ishft(ix,-16)
        il = iand(ix,65535)
        mt(NMAX-N+i) = ieor(ishft(iand(il*23896,65535),16)+
     *    ishft(iand(ih*35685,65535),16)+il*35685,mt(NMAX-N+i))-i+1
        i = i+1
        if (i .gt. N) then
          mt(NMAX-N+1) = mt(NMAX)
          i = 2
        endif
 1100 continue
      mt(NMAX-N+1) = ieor(2147483647,-1)
      call mtrnu
      mti = N
      return
      end
*
      subroutine mtrni(irand)
      integer irand
      parameter(N = 624, M = 397)
      parameter(NVEC = N-M-1, NLOOP = 50, NMAX = NVEC*NLOOP+N)
      integer mt(NMAX)
      common /mtrnc/ mti, mt
      save   /mtrnc/
      if (mti .ge. NMAX-N) call mtrnu
      mti = mti+1
      irand = ishft(mt(mti),-1)
      return
      end
*
      subroutine mtrniv(lbuf, nbuf)
      integer nbuf
      integer lbuf(nbuf)
      parameter(N = 624, M = 397)
      parameter(NVEC = N-M-1, NLOOP = 50, NMAX = NVEC*NLOOP+N)
      integer mt(NMAX)
      common /mtrnc/ mti, mt
      save   /mtrnc/
      ibuf = 0
 1000 if (nbuf-ibuf .ge. NMAX-N-mti) then
        do 1100 i = 1, NMAX-N-mti
          lbuf(ibuf+i) = ishft(mt(mti+i),-1)
 1100   continue
        ibuf = ibuf + NMAX-N-mti
        call mtrnu
        goto 1000
      endif
      if (ibuf .lt. nbuf) then
        do 1200 i = 1, nbuf-ibuf
          lbuf(ibuf+i) = ishft(mt(mti+i),-1)
 1200   continue
        mti = mti + nbuf-ibuf
      endif
      return
      end
*
      subroutine mtrnj(irand)
      integer irand
      parameter(N = 624, M = 397)
      parameter(NVEC = N-M-1, NLOOP = 50, NMAX = NVEC*NLOOP+N)
      integer mt(NMAX)
      common /mtrnc/ mti, mt
      save   /mtrnc/
      if (mti .ge. NMAX-N) call mtrnu
      mti = mti+1
      irand = mt(mti)
      return
      end
*
      subroutine mtrnjv(lbuf, nbuf)
      integer nbuf
      integer lbuf(nbuf)
      parameter(N = 624, M = 397)
      parameter(NVEC = N-M-1, NLOOP = 50, NMAX = NVEC*NLOOP+N)
      integer mt(NMAX)
      common /mtrnc/ mti, mt
      save   /mtrnc/
      ibuf = 0
 1000 if (nbuf-ibuf .ge. NMAX-N-mti) then
        do 1100 i = 1, NMAX-N-mti
          lbuf(ibuf+i) = mt(mti+i)
 1100   continue
        ibuf = ibuf + NMAX-N-mti
        call mtrnu
        goto 1000
      endif
      if (ibuf .lt. nbuf) then
        do 1200 i = 1, nbuf-ibuf
          lbuf(ibuf+i) = mt(mti+i)
 1200   continue
        mti = mti + nbuf-ibuf
      endif
      return
      end
*
      subroutine mtrnd(drand)
      implicit double precision(a-h,o-z)
      double precision drand
      parameter(N = 624, M = 397)
      parameter(NVEC = N-M-1, NLOOP = 50, NMAX = NVEC*NLOOP+N)
      parameter(UNIT = 1.0d0/2147483648.0d0)
      integer mt(NMAX)
      common /mtrnc/ mti, mt
      save   /mtrnc/
      if (mti .ge. NMAX-N) call mtrnu
      mti = mti+1
      drand = dble(ishft(mt(mti),-1))*UNIT
      return
      end
*
      subroutine mtrnf(drand)
      implicit double precision(a-h,o-z)
      double precision drand
      parameter(N = 624, M = 397)
      parameter(NVEC = N-M-1, NLOOP = 50, NMAX = NVEC*NLOOP+N)
      parameter(UNIT2 = 67108864.0d0)
      parameter(UNIT3 = 1.0d0/9007199254740992.0d0)
      integer ia, ib
      integer mt(NMAX)
      common /mtrnc/ mti, mt
      save   /mtrnc/
      if (mti .ge. NMAX-N) call mtrnu
      mti = mti+1
      ia = ishft(mt(mti),-5)
      if (mti .ge. NMAX-N) call mtrnu
      mti = mti+1
      ib = ishft(mt(mti),-6)
      drand = (dble(ia)*UNIT2 + dble(ib))*UNIT3
      return
      end
*
      subroutine mtrndv(dbuf, nbuf)
      implicit double precision(a-h,o-z)
      integer nbuf
      double precision dbuf(nbuf)
      parameter(N = 624, M = 397)
      parameter(NVEC = N-M-1, NLOOP = 50, NMAX = NVEC*NLOOP+N)
      parameter(UNIT = 1.0d0/2147483648.0d0)
      integer mt(NMAX)
      common /mtrnc/ mti, mt
      save   /mtrnc/
      ibuf = 0
 1000 if (nbuf-ibuf .ge. NMAX-N-mti) then
        do 1100 i = 1, NMAX-N-mti
          dbuf(ibuf+i) = dble(ishft(mt(mti+i),-1))*UNIT
 1100   continue
        ibuf = ibuf + NMAX-N-mti
        call mtrnu
        goto 1000
      endif
      if (ibuf .lt. nbuf) then
        do 1200 i = 1, nbuf-ibuf
          dbuf(ibuf+i) = dble(ishft(mt(mti+i),-1))*UNIT
 1200   continue
        mti = mti + nbuf-ibuf
      endif
      return
      end
*
      subroutine mtrnfv(dbuf, nbuf)
      implicit double precision(a-h,o-z)
      integer nbuf
      double precision dbuf(nbuf)
      parameter(NWORK = 5000)
      parameter(UNIT2 = 67108864.0d0)
      parameter(UNIT3 = 1.0d0/9007199254740992.0d0)
      integer work(NWORK*2)
      ibuf = 0
 1000 if (nbuf-ibuf .ge. NWORK) then
        call mtrnjv(work, NWORK*2)
        do 1100 i = 1, NWORK
          dbuf(ibuf+i) = (dble(ishft(work(2*i-1),-5))*UNIT2
     *      +dble(ishft(work(2*i),-6)))*UNIT3
 1100   continue
        ibuf = ibuf + NWORK
        goto 1000
      endif
      if (ibuf .lt. nbuf) then
        call mtrnjv(work, (nbuf-ibuf)*2)
        do 1200 i = 1, nbuf-ibuf
          dbuf(ibuf+i) = (dble(ishft(work(2*i-1),-5))*UNIT2
     *      +dble(ishft(work(2*i),-6)))*UNIT3
 1200   continue
      endif
      return
      end
*
      subroutine mtrns(lbuf)
      parameter(N = 624, M = 397)
      parameter(MASKTC = -272236544)
      parameter(MASKB1 = 5760)
      parameter(MASKB2 = 802816)
      parameter(MASKB3 = 220200960)
      parameter(MASKB4 = -1879048192)
      parameter(MASKA1 = 2096128)
      parameter(MASKA2 = 1023)
      parameter(NVEC = N-M-1, NLOOP = 50, NMAX = NVEC*NLOOP+N)
      integer lbuf(N)
      integer mt(NMAX)
      common /mtrnc/ mti, mt
      save   /mtrnc/
      do 1000 i = 1, N
        lbuf(i) = mt(mti+i)
 1000 continue
      imax = N
      if (mti+N .gt. NMAX-N) imax = NMAX-N - mti
      do 2000 i = 1, imax
        iy = lbuf(i)
        iy = ieor(iy,ishft(iy,-18))
        iy = ieor(iy,iand(ishft(iy,15),MASKTC))
        iy = ieor(iy,iand(ishft(iy,7),MASKB1))
        iy = ieor(iy,iand(ishft(iy,7),MASKB2))
        iy = ieor(iy,iand(ishft(iy,7),MASKB3))
        iy = ieor(iy,iand(ishft(iy,7),MASKB4))
        iy = ieor(iy,iand(ishft(iy,-11),MASKA1))
        iy = ieor(iy,iand(ishft(iy,-11),MASKA2))
        lbuf(i) = iy
 2000 continue
      return
      end
*
      subroutine mtrnl(lbuf)
      parameter(N = 624, M = 397)
      parameter(NVEC = N-M-1, NLOOP = 50, NMAX = NVEC*NLOOP+N)
      integer lbuf(N)
      integer mt(NMAX)
      common /mtrnc/ mti, mt
      save   /mtrnc/
      do 1000 i = 1, N
        mt(NMAX-N+i) = lbuf(i)
 1000 continue
      call mtrnu
      return
      end
