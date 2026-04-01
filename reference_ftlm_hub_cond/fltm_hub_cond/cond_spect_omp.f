c======================================================================
c     Thermodynamics calculation
c     JK 12.4.2012 - each dkx,dky (or twisted boundary phase separately
c          and then final quantites averaged over - should reproduce
c          accurately NIE result. Different phases not just sumed
c          togehter in thermodynamic sum (Z).
c          From thermodyn_old.f (macbook))
c======================================================================
      program cond_spect
      implicit real*8 (a-h, o-z)
      implicit integer (i-n)
c-----------------------------------------------------------------------
c     INPUT FILE
c     number different idk's files to read (001,002,....,nfididk)
c     for sets (dkx,dky) and (dkxd,dkyd)
      parameter (nfididk=2) !number of diff. phases in file names
      parameter (nfidsmp=4) !number of diff. samples for phase (file names)
c-----------------------------------------------------------------------
c     name of the data file for results
      character ndf*11
c     name of the data file for results
      data ndf /'n08u120tp00'/
      character fid*3
      character fidd*3
      character fidsmp*3
c     for number of ndks in the name of resultfiles
      character cnndk*3
c-----------------------------------------------------------------------
c     wave vector
      real*8 kx,ky 
c-----------------------------------------------------------------------
c     twisted bounary condition average
c      (first result for each pase, than average)
c     maximal number of phases
      parameter(ndkmax=nfididk)
      real*8 dklist(4,ndkmax)
c-----------------------------------------------------------------------
c    chemical potential (mu) variable and values
      parameter (nmu=61) !number of diff. values
      parameter(dmumin=-4.d0, dmumax=10.d0) !interval of mu values
      real*8 mulist(nmu)
c     to calculate \tilde Z: min H-mu N_e or reference energy for each mu
      real*8 hmmumin(nmu,ndkmax)
c-----------------------------------------------------------------------
c    temperature (tt) variables and values
      parameter (ntt=50) !number of diff. temp. values
      parameter (ttmin=0.1d0, ttmax=10.1d0) !interval of temp. values
      real*8 ttlist(ntt)
c-----------------------------------------------------------------------
c    conductivity channel variables
      parameter (nchan=4001) !number of channels forcoduct (even so one chanel centered around om=0)
      parameter (ommin=-20.d0, ommax=20.0d0) !interval of temp. values
      parameter (eta=0.1d0) !eta broadening
      parameter (nchan2=4001) !number points in omega for eta-broadened conduct
      parameter (ommin2=-20.d0, ommax2=20.0d0) !interval of temp. values
c-----------------------------------------------------------------------
c     lanczos basis variables
      parameter (lstmax=200)
c     eigenenergies
      real*8 en1(lstmax)
c     overlap with random vector
      real*8 vka1(lstmax)
      character fmt*3
c-----------------------------------------------------------------------
      complex*16 crvka
c-----------------------------------------------------------------------
c     optical conductivity
c     second lanczos
      real*8 enb(lstmax)
      real*8 vkb1(lstmax)
      complex*16 crvkavib8(lstmax,lstmax)
      complex*16 crvkavib10(lstmax)
      complex*16 cvkb1e1(lstmax)
      complex*16 cdcond
c-----------------------------------------------------------------------
c     tau oper (kin en like oper) matrix elements
      real*8 rvkavib9(lstmax,lstmax)
      real*8 tau(nmu,ntt,ndkmax)
c-----------------------------------------------------------------------
c     curr oper (curr expectation value) for charge stiffnes
       complex*16 curr(nmu,ntt,ndkmax)
       complex*16 dcurr
c     curr oper^2 (curr expectation value) for charge stiffnes
       complex*16 curr2(nmu,ntt,ndkmax)
       complex*16 dcurr2,dcurr20,dcurr2a
c-----------------------------------------------------------------------
c     quantities, temp. and mu dependent
c     Z tilde
      real*8 zt(nmu,ntt,ndkmax)
c     nf average (number of fermions, electrons)
      real*8 nfa(nmu,ntt,ndkmax)
c-----------------------------------------------------------------------
c     optical conductivity 
      complex*16 cond(nmu,ntt,ndkmax)
c     channels
      complex*16 condchan(nmu,ntt,ndkmax,nchan)
      complex*16 dcondchan(nchan)
      real*8 dcondeta(nchan2)
      real*8 dcondeta1(nchan2)
      real*8 dcondeta2(nchan2)
      real*8 dcondeta3(nchan2)
      real*8 dcondeta4(nchan2)
      real*8 dcondeta5(nchan2)
      complex*16 dcond,ddcond
c-----------------------------------------------------------------------
c     list of nf - electron densities to calcualte T-dependence of
c             - chemical potential
c             - z sum
      parameter(nnfd=7)         !number of densitis
      real*8 dlnfd(nnfd) !list of density values
      data dlnfd / 1.d0,0.95d0,0.90d0,0.85d0,0.80d0,0.75d0,0.70d0/  !density values
      character str*4 !string for file name number
c     folder name for data files (on /mnt/lustre)
      character fldr*100
c-----------------------------------------------------------------------
c     TIMING
      real t0,ttime(2)
c-----------------------------------------------------------------------
      data  eps1, eps2/ 1.d-7, 1.d-5/
c=======================================================================
      pi=4.d0*datan(1.d0)      
c     chemical potential variables: values of mu in mulist
      write(0,*) "Chemical pot. nmu,dmumin,dmumax:",nmu,dmumin,dmumax
      write(0,*) "Temperature   ntt,ttmin,ttmax:",ntt,ttmin,ttmax
      ddmu=(dmumax-dmumin)/nmu
      do i=1,nmu
         mulist(i)=dmumin +i*ddmu
      enddo
c-----------------------------------------------------------------------     
c     temperature variables: values of mu in mulist
      dtt=(ttmax-ttmin)/ntt
      do i=1,ntt
         ttlist(i)=ttmin +(i-1)*dtt
      enddo
c-----------------------------------------------------------------------     
c     seting hmmumin to starting value (min H-mu*Ne)
      do i=1,nmu
       do idk=1,ndkmax  
         hmmumin(i,idk)=1.d10
       enddo  
      enddo
c-----------------------------------------------------------------------           
c     setting dklist to 0.d0
      do jk=1,4
         do idk=1, ndkmax
            dklist(jk,idk)=0.0d0
         enddo
      enddo
c     number of different phases determined from dkx,dky
      nndk=0
c-----------------------------------------------------------------------
c     READING FROM FILE
c-----------------------------------------------------------------------
      nf0=0
      ikold=0
c     switch if something wrong in data files
      ierr=0
      write(0,*) "starting first files reading "
c     loop over data files, idk, idkd
      do ifid=1,nfididk
      do ifidsmp=1,nfidsmp
      write(fid,250) ifid/100,(ifid-ifid/100*100)/10,(ifid-ifid/10*10)
c      write(fidd,250) ifidd/100,(ifidd-ifidd/100*100)/10,
c     * (ifidd-ifidd/10*10)
      write(fidsmp,250) ifidsmp/100,(ifidsmp-ifidsmp/100*100)/10,
     * (ifidsmp-ifidsmp/10*10)
 250   format(3i1)
      in=20
      write(0,*)'./hubTri2Dconddata/'//ndf//'.dat_idk'//fid//
     * '_ismp'//fidsmp
c     Local data folder
      write(fldr,*)'.'
       open(unit=in,file=trim(adjustl(fldr))//'/hubTri2Dconddata/'
     *  //ndf//'.dat_idk'//fid//'_ismp'//fidsmp,
     *  status='old',access='sequential',form='unformatted',iostat=kode)
c-----------------------------------------------------------------------
c     LOOP OVER DATA IN FILE
c-----------------------------------------------------------------------
      t0=dtime(ttime)
      write(0,*) "starting reading file"
c-----------------------------------------------------------------------
 100  continue
        read(in,end=900,err=800) nn0, nf, nu
c        write(0,*) 'nn0, nf , nu = ',nn0, nf, nu
        if(nf.ne.nf0) then
           write(0,*) "   reading data for nf =",nf,nf0
c          checking right order of nf data           
           if ((nf.ne.(nf0+1)).and.(nf0.ne.2*nn0)) then
             write(0,*) "   warning: from nf=",nf0," to nf=",nf
c             stop "ERROR in data file"
             ierr=1
           endif
           nf0=nf
        endif
        read(in) ik, kx, ky, dkx, dky, dkxd, dkyd, iismp     
c       checking right order of index ik
        if(ik.ne.ikold) then
c           write(0,*) "    reading data for ik =",ik
c          checking right order of ik in data           
           if (ik.ne.(ikold+1).and.nf.ne.1.and.ik.ne.1.and.
     *       (nu.ne.1.and.nf.ne.(nn0).and.idk.ne.1).and.
     *       nf.ne.2*nn0
     *       ) then
               write(0,*) "   warning: from idk=",ikold," to idk=",ik
c               stop "ERROR in data file"
               ierr=1
           endif
           ikold=ik
        endif
c      saving (dkx,dky) and (dkxd,dkyd) into dklist and determining idk
c        write(0,*) dkx,dky,dkxd,dkyd
        if(nndk.eq.0) then
           nndk=1
           idk=1
           dklist(1,idk)=dkx
           dklist(2,idk)=dky
           dklist(3,idk)=dkxd
           dklist(4,idk)=dkyd
        else
           idk=0
           do iidk=1, nndk
              if(dklist(1,iidk).eq.dkx.and.dklist(2,iidk).eq.dky.and.
     *           dklist(3,iidk).eq.dkxd.and.dklist(4,iidk).eq.dkyd) then
                 idk=iidk
              endif   
           enddo
           if(idk.eq.0) then
              nndk=nndk+1
              idk=nndk
              dklist(1,idk)=dkx
              dklist(2,idk)=dky
              dklist(3,idk)=dkxd
              dklist(4,idk)=dkyd
           endif   
        endif
        read(in) npw, nsmp, ndkl
         read(in,end=900,err=800) nex
c        loop over lanczos eigenstates
         do j=1,nex
          read(in,end=900,err=800) en, vn
c          write(0,*) en, vn

c        finding minimal H -mu*nf
          if (j.eq.1) then
             do jj=1,nmu
                if(hmmumin(jj,idk).gt.(en-mulist(jj)*dble(nf))) then
                   hmmumin(jj,idk)= en-mulist(jj)*dble(nf)
                endif
             enddo
          endif                                   
            
         enddo !j  
c        enddo !ismp

c        reading ba(-a) - norm of starting random vector
         read(in,end=900,err=800) bam1
c-----------------------------------------------------------------------
c        reading optical conductivity parameters
         read(in,end=900,err=800) nexb
         read(in,end=900,err=800) bbm1
         if (nexb.gt.lstmax) stop "ERROR, lstmax <nexb"
         read(in,end=900,err=800) (enb(ia),vkb1(ia),ia=1,nexb)
         read(in,end=900,err=800) ((crvka,kb=1,nexb),ja=1,nex)
         read(in,end=900,err=800) (crvka,ka=1,nex) !<lan1|j|lan1>
c-----------------------------------------------------------------------
c        reading tau oper (like kin en) matrix elements
         read(in,end=900,err=800) ((rij,ja=1,nex),ia=1,nex)
c-----------------------------------------------------------------------
      goto 100  !read another set of data - data loop

      
 800    write(0,*) "error in reading file"
 900    write(0,*) "end of reading file"
        close(in)

        enddo ! ifidsmp (over different sampes)
        enddo ! ifid (over files from different idk)
        write(0,*) "end of reading files"
c       stop if problem in data files
        write(0,*) "ierr for testing ",ierr
        if(ierr.eq.1) then
           write(0,*) "ERORR in data files - missing data?"
           stop "ERORR in data files - missing data?"
        endif
c-----------------------------------------------------------------------
c     write dklist
       if ((nndk).ne.(ndkl)) then
         write(0,*) "error:wrong num. of phases found",nndk,ndkl
       endif
       write(0,*) "  phases indices:" 
       do idk=1,nndk
          write(0,*) "  idk dkx dky dkxd dkyd", idk,
     *      dklist(1,idk),dklist(2,idk),dklist(3,idk),dklist(4,idk)
       enddo
c-----------------------------------------------------------------------
c     nndk to string
      write(cnndk,251)nndk/100,(nndk-nndk/100*100)/10,(nndk-nndk/10*10)
251   format(3i1)       
c-----------------------------------------------------------------------
        t0=dtime(ttime)
        write(0,*) "full first reading time in seconds:", ttime(1)
c-----------------------------------------------------------------------
c======================================================================
c     CALCULATIONG THERMODYNAMICS
c-----------------------------------------------------------------------
      do imu=1, nmu
         do itt=1, ntt
           do idk=1, ndkmax
            zt(imu,itt,idk)=0.d0
            nfa(imu,itt,idk)=0.d0
c           optical conductivity
            cond(imu,itt,idk)=dcmplx(0.d0,0.d0)
            do ichan=1,nchan
              condchan(imu,itt,idk,ichan)=0.d0
            enddo
            tau(imu,itt,idk)=0.d0
            curr(imu,itt,idk)=0.d0
            curr2(imu,itt,idk)=0.d0
           enddo
         enddo !itt
       enddo !imu
      
c     SECOND READING FROM FILE
c-----------------------------------------------------------------------     
c     loop over all files or treads ids
      do ifid=1, nfididk
      do ifidsmp=1, nfidsmp
      write(fid,250) ifid/100,(ifid-ifid/100*100)/10,(ifid-ifid/10*10)
      write(fidsmp,250) ifidsmp/100,(ifidsmp-ifidsmp/100*100)/10,
     * (ifidsmp-ifidsmp/10*10)
c 250   format(3i1)
      write(0,*) 'fid ',fid
      in=20
      write(0,*)'./hubTri2Dconddata/'//ndf//'.dat_idk'//fid//
     * '_ismp'//fidsmp 
c     Local data folder
      write(fldr,*)'.'
       open(unit=in,file=trim(adjustl(fldr))//'/hubTri2Dconddata/'
     *  //ndf//'.dat_idk'//fid//'_ismp'//fidsmp,
     *  status='old',access='sequential',form='unformatted',iostat=kode)
c     ifile - flag to be zero if no data in file - or file doesnt exsist
      ifile=0
c-----------------------------------------------------------------------
c     LOOP OVER DATA IN FILE
      write(0,*) "starting second file reading, ",ndf
      nf0=0
 1100 continue
        read(in,end=1900,err=1800) nn0, nf, nu
c        write(0,*) 'nn0, nf , nu = ',nn0, nf, nu
        ifile=1
        if(nf.ne.nf0) then
           write(0,*) "   reading data for nf =",nf
           nf0=nf
        endif
        read(in) ik, kx, ky, dkx, dky, dkxd, dkyd, iismp     
c      finding phase index idk
        idk=0
        do iidk=1, nndk
           if(dklist(1,iidk).eq.dkx.and.dklist(2,iidk).eq.dky.and.
     *      dklist(3,iidk).eq.dkxd.and.dklist(4,iidk).eq.dkyd) then
              idk=iidk
c              write(0,*) "idk ", idk
           endif   
        enddo
        if(idk.eq.0) then
           write(0,*) "error: phase index not found in sec. file read"
        endif
c        write(0,*) "   idk dkx dky", idk,dkx,dky
        read(in) npw, nsmp, ndkl
         read(in,end=1900,err=1800) nex
c         write(0,*) 'nex = ',nex
c        loop over lanczos eigenstates
         do j=1,nex
            read(in,end=1900,err=1800) en, vn
            en1(j)=en
            vka1(j)=vn
         enddo !j

c        reading ba(-1) - norm of starting random vector
         read(in,end=1900,err=1800) bam1

c           CALCULATION FOR EACH mu AND tt SEPARATELY
!$OMP PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(imu,dmu,emin,itt,tt,wgh,j)
            do imu=1, nmu
             dmu=mulist(imu)
             emin=hmmumin(imu,idk)
             do itt=1, ntt
               tt=ttlist(itt)
              do j=1,nex
               wgh=1.d0*npw*(vka1(j)**2)/nsmp
     *             *dexp(-(en1(j)-dmu*nf-emin)/tt)
               
c    spin symmetry
c    since result for exchanged nu and nd is the same and not explicitly calculated
               if(nu.ne.(nf-nu)) wgh=wgh*2.d0
               zt(imu,itt,idk)=zt(imu,itt,idk)+wgh
               nfa(imu,itt,idk)=nfa(imu,itt,idk)+nf*wgh              
               enddo !j  
               
             enddo !itt
            enddo !imu            
!$OMP END PARALLEL DO
         
c-----------------------------------------------------------------------
c        reading optical conductivity parameters
         read(in,end=1900,err=1800) nexb
         read(in,end=1900,err=1800) bbm1
         read(in,end=1900,err=1800) (enb(ia),vkb1(ia),ia=1,nexb)
c      current oper. matrix between second and first lanczos eigen. basis
        read(in,end=1900,err=1800)
     *    ((crvkavib8(kb,ja),kb=1,nexb),ja=1,nex)
c      current oper. matrix between first and first lanczos eigen. basis
        read(in,end=1900,err=1800)
     *    (crvkavib10(ka),ka=1,nex) !<lan1|j|lan1>
c-----------------------------------------------------------------------
c        reading tau oper (kin en like oper) matrix elements
         read(in,end=1900,err=1800)((rvkavib9(ia,ja),ja=1,nex),ia=1,nex)
c-----------------------------------------------------------------------
c           CALCULATION FOR EACH mu AND tt SEPARATELY
!$OMP PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(imu,dmu,emin,itt,tt,ia,fac1,ja,fac,ddocc,
!$OMP&  ddoccen,do1,den1,dborder1,dborder2,     
!$OMP&  lb,kb,fac4,fac3,ddcond,dom,ichan,cdcond,dtau,dcurr,dcurr2)     
         do imu=1, nmu
          dmu=mulist(imu)
          emin=hmmumin(imu,idk)
         do itt=1, ntt
            tt=ttlist(itt)
            do ia=1,nex
c    fac1: since result for exchanged nu and nd is the same and not explicitly calculated
              fac1=1.d0
              if(nu.ne.(nf-nu)) fac1=2.d0
              do ja=1,ia
                fac=2.d0 ! non-diagonal matrix elements twice 2 (ij+ji)
                if(ja.eq.ia) fac=1.d0
                
              enddo !ja
            enddo !ia

c          tau oper 
            do ia=1,nex
              fac1=1.d0
              if(nu.ne.(nf-nu)) fac1=2.d0
                 ja=ia  !vka1(i)*vka1(j)=delta_ij
                 fac=1.d0
c               tau oper
                dtau=1.d0*npw/nsmp
     *            *vka1(ia)*rvkavib9(ia,ja)*vka1(ja)
     *            *dexp(-(en1(ia)-dmu*nf-emin)/tt)
     *            *fac*fac1
                tau(imu,itt,idk)=tau(imu,itt,idk)+dtau
            enddo !ia
            
c          curr oper 
            do ia=1,nex
              fac1=1.d0
              if(nu.ne.(nf-nu)) fac1=2.d0
                 ja=ia !vka1(i)*vka1(j)=delta_ij  if N_random large
                 fac=1.d0
c               curr oper
                dcurr=1.d0*npw/nsmp
     *            *vka1(ia)*dimag(crvkavib10(ia))*vka1(ja)
     *            *dexp(-(en1(ia)-dmu*nf-emin)/tt)
     *            *fac*fac1
                dcurr2=dcurr*dimag(crvkavib10(ia))
                curr(imu,itt,idk)=curr(imu,itt,idk)+dcurr
                curr2(imu,itt,idk)=curr2(imu,itt,idk)+dcurr2
            enddo !ia
            
c---------
c           optical conductivity 
            do ia=1,nex
              fac1=1.d0
              if(nu.ne.(nf-nu)) fac1=2.d0
                fac=1.d0
                cdcond=dcmplx(0.d0,0.d0)
                do lb=1,nexb
                 if (dabs(enb(lb)-en1(ia))/tt.gt.eps1) then    
                  fac3=(dexp(-(en1(ia)-dmu*nf-emin)/tt)
     *             -dexp(-(enb(lb)-dmu*nf-emin)/tt))
     *            /(enb(lb)-en1(ia))
                  else
                  fac3=dexp(-(en1(ia)-dmu*nf-emin)/tt)/tt
c                  fac3=0.d0
                 endif
                 cdcond=+1.d0*npw/nsmp
     *            *vka1(ia)
     *            *dconjg(crvkavib8(lb,ia))
     *            *vkb1(lb)*bbm1
     *            *fac*fac1*fac3*pi
     *            /bam1
                 ddcond=cdcond
                 cond(imu,itt,idk)=cond(imu,itt,idk)+ddcond
c               optical conduct to channels
                 dom=(ommax-ommin)/nchan
                 ichan= int((enb(lb)-en1(ia)-ommin)/dom)
                 if(ichan.lt.1) ichan=1
                 if(ichan.gt.nchan) ichan=nchan
                 condchan(imu,itt,idk,ichan)=
     *             condchan(imu,itt,idk,ichan)+ddcond
              enddo             !lb
            enddo !ia

          enddo !itt
         enddo !imu            
!$OMP END PARALLEL DO
         
      goto 1100  !read anoter set of data - data loop

 1800   write(0,*) "error in second reading file"
 1900   write(0,*) "end of second reading file"
        close(in)
        
      enddo    !ifidsmp
      enddo    !ifid
c-----------------------------------------------------------------------
        t0=dtime(ttime)
        write(0,*) "full sec. reading time in seconds:", ttime(1)
c-----------------------------------------------------------------------
c     Writing results into file
c----------------------------------------------------------------------
      open(unit=40,file= './hubTri2Dcondres/'//ndf//'_ndk'//cnndk//
     * '_zt.dat',
     *  status='unknown',access='sequential',form='formatted')
      open(unit=41,file= './hubTri2Dcondres/'//ndf//'_ndk'//cnndk//
     * '_nf.dat',
     *  status='unknown',access='sequential',form='formatted')
      open(unit=55,file= './hubTri2Dcondres/'//ndf//'_ndk'//cnndk//
     * '_cond.dat',
     *  status='unknown',access='sequential',form='formatted')
      open(unit=56,file= './hubTri2Dcondres/'//ndf//'_ndk'//cnndk//
     * '_tau.dat',
     *  status='unknown',access='sequential',form='formatted')

       do imu=1, nmu
         dmu=mulist(imu)
         do itt=1, ntt
            tt=ttlist(itt)
            dzt=0.d0
            dnfa=0.d0
            dcond=dcmplx(0.d0,0.d0)
            dtau=0.d0
            dcurr=dcmplx(0.d0,0.d0)
            dcurr2=dcmplx(0.d0,0.d0)
            do idk=1,nndk
             dzt=dzt+zt(imu,itt,idk)/nndk
             dnfa=dnfa+nfa(imu,itt,idk)/zt(imu,itt,idk)/nndk
             dcond=dcond+ cond(imu,itt,idk)/zt(imu,itt,idk)/nndk
             dtau=dtau+
     *                tau(imu,itt,idk)/zt(imu,itt,idk)/nndk
            enddo
            write(40,*) dmu, tt, dzt
            write(41,*) dmu, tt, dnfa
            write(fmt,'(I3)') nn0
            write(55,*) dmu, tt, dcond
            write(56,*) dmu, tt, dtau
            ent1=0.d0
            dzt=zt(imu,itt,1)
            emin=hmmumin(imu,1)
         enddo !itt
         write(40,*) ' '
         write(41,*) ' '
         write(55,*) ' '
         write(56,*) ' '
       enddo !imu
      close(40)      
      close(41)      
      close(55)      
      close(56)      
c-----------------------------------------------------------------------
c     Writing results for particular nf into file 
c     chemical potential calculated from average nf
c-----------------------------------------------------------------------
c     Loop over all densities
c     file number nfile= starts with 100
      write(0,*) "starting to write to _ndf10x files"
      nfile=100
      do infd=1,nnfd
       write(str,1000) nfile
 1000   format (I3)
c       write(0,*) str
       write(0,*) "   writing to: ",'./hubTri2Dcondres/'//ndf//'_ndk'
     *  //cnndk//'_nfd'//str
       open(unit=80,file= './hubTri2Dcondres/'//ndf//'_ndk'//cnndk//
     *  '_nfd'//str,
     *  status='unknown',access='sequential',form='formatted')
c     storing results for conductivity
       open(unit=85,file= './hubTri2Dcondres/'//ndf//'_ndk'//cnndk//
     *  '_condchan_nfd'//str,
     *  status='unknown',access='sequential',form='formatted')
c     storing broadened conductivity to file condeta
       open(unit=86,file= './hubTri2Dcondres/'//ndf//'_ndk'//cnndk//
     *  '_condeta_nfd'//str,
     *  status='unknown',access='sequential',form='formatted')
c     storing tau oper expectaion value to file _tau_ndk
       open(unit=87,file= './hubTri2Dcondres/'//ndf//'_ndk'//cnndk//
     *  '_tau_nfd'//str,
     *  status='unknown',access='sequential',form='formatted')
       
       dnfd=dlnfd(infd)
       do itt=1, ntt
         tt=ttlist(itt)
         dnfdiff=100.d0
         imuflag=0
         iimulow=0
         iimuhgh=0
         iimu=0
         do imu=2, nmu-1
           dnfalow=0.d0
           dnfaup=0.d0
           do idk=1,nndk
            dnfalow=dnfalow+nfa(imu,itt,idk)/zt(imu,itt,idk)/nn0/nndk
            dnfaup=dnfaup+nfa(imu+1,itt,idk)/zt(imu+1,itt,idk)/nn0/nndk
           enddo
           if((dabs(dnfalow-dnfd).lt.dabs(dnfdiff)).or.
     *        (dabs(dnfalow-dnfd).lt.1.d-9)) then
              iimu=imu
              dnfdiff=dnfalow-dnfd
           endif
c          set imu for lowest value of flat (gaped) regime
           if((dabs(dnfalow-dnfd).lt.1.d-9).and.imuflag.eq.0) then
              iimulow=imu
              imuflag=1
           endif
c          set imu for highest value of flat (gaped) regime
           if((dabs(dnfalow-dnfd).lt.1.d-9).and.imuflag.eq.1) then
              iimuhgh=imu
           endif
c          set imu for dnf and quantites calc.
         enddo !imu
         if(iimulow.ne.0.and.iimuhgh.ne.0) iimu=(iimulow+iimuhgh)/2
         imu=iimu
           if (imu.ne.0.and.dabs(dnfdiff).lt.1.d-1) then
              dmu=mulist(imu)
              dzt=0.d0
              dnfa=dnfalow*nn0
              dnfa1=dnfalow*nn0              
              dena=0.d0
              demin=0.d0
              do ichan=1,nchan
                dcondchan(ichan)=dcmplx(0.d0,0.d0)
              enddo
              dtau=0.d0
              dcurr=dcmplx(0.d0,0.d0)
              dcurr2=dcmplx(0.d0,0.d0)
              dcurr20=dcmplx(0.d0,0.d0)
              dcurr2a=dcmplx(0.d0,0.d0)
              do idk=1,nndk
                dzt=dzt+zt(imu,itt,idk)/nndk
                demin=demin+hmmumin(imu,idk)/nndk
c               conductivity channels
                do ichan=1,nchan
                   dcondchan(ichan)=dcondchan(ichan)+
     *               condchan(imu,itt,idk,ichan)/zt(imu,itt,idk)/nndk
                enddo
                dtau=dtau+
     *                   tau(imu,itt,idk)/zt(imu,itt,idk)/nndk
                dcurr=dcurr+
     *                   curr(imu,itt,idk)/zt(imu,itt,idk)/nndk
c                D_c (charge stiffnes from 1st lanczos. TBC averged eq 6 in zotos1997 )
                dcurr2=dcurr2+0.5d0*(curr2(imu,itt,idk)/zt(imu,itt,idk)
     *           -(curr(imu,itt,idk)/zt(imu,itt,idk))**2)/tt/nndk
c               sum_i |<i|j|i>|^2  ... i -first lanczos basis
                dcurr20=dcurr20+
     *            curr2(imu,itt,idk)/zt(imu,itt,idk)/nndk
                
c             (sum_i <i|j|i>)^2  ... i -first lanczos basis (to correct om=0 chan)
                dcurr2a=dcurr2a+
     *            (curr(imu,itt,idk)/zt(imu,itt,idk))**2/nndk
              enddo
              ent1=0.d0
              fen1=0.d0
              dcv=0.d0
              dchis= 0.d0
              dchic= 0.d0
              write(80,1001) dnfd, tt, dmu, dzt, dena, ent1, dcv,
     *         dchis, dchic, fen1
 1001         format (3F10.4,7E14.5)

c            om=0 channel
              om=0.d0
              dom=(ommax-ommin)/nchan
              ichan0= int((om-ommin)/dom)

c            integral of optical conudctivity
              dsum=0.d0
              do ichan=1,nchan
                 dsum=dsum+dreal(dcondchan(ichan))
              enddo
              dsum=dsum/pi
c            integral of optical conudctivity (just om>0)
              dsumpoz=0.d0
              do ichan=ichan0+1,nchan
                 dsumpoz=dsumpoz+dreal(dcondchan(ichan))
              enddo
              dsumpoz=2*dsumpoz/pi
c            integral of optical conudctivity without om=0 channel
              dsum1=dsum-dreal(dcondchan(ichan0))/pi
c            charge stiffness
              dcstiff=0.5d0*(dtau-(dsum-dreal(dcondchan(ichan0))/pi))              
c           D_c from phase sum_i<J_i>^2-(sum_i<J_i>)^2 
            dcstiff1=dreal(dcurr2)   
            write(87,*)dnfd,tt,dmu,dtau,dsum,dsum1,dsumpoz,dcstiff
     *            ,dreal(dcondchan(ichan0)),dreal(dcurr20) 
     *            ,dreal(dcurr),dimag(dcurr),dreal(dcurr2) 
     *            ,dreal(dcurr2a)*4*pi/tt

 1002         format (3F10.4,6E14.5)
              write(fmt,'(I3)') nn0
              
c             write conductivity chan to file
              dom=(ommax-ommin)/nchan
              do ichan=1,nchan
                 om=ommin+ichan*dom +dom/2
c              --with D_c from tau
                 ddres=dreal(dcondchan(ichan)) 
                 if (ichan.eq.ichan0) ddres=2.d0*pi*dcstiff
c              --with subtracted om=0 part for finite current (dcurr2a)
                 ddres1=dreal(dcondchan(ichan))
                 if (ichan.eq.ichan0)
     *            ddres1=ddres1-dreal(dcurr2a)*4*pi/tt   
c              --with subtracted om=0 part for finite current (dcurr2a)
c                and subtracted D_C from phases (dcurr2) (just regular part)
                 ddres2=dreal(dcondchan(ichan))
                 if (ichan.eq.ichan0)
     *           ddres2=ddres2-dreal(dcurr2a)*4*pi/tt-dreal(dcurr2)*2*pi   
c              --with om=0 set to 0
                 ddres3=dreal(dcondchan(ichan))
                 if (ichan.eq.ichan0)  ddres3=0.d0
c              --with D_C from phase (dcurr2)
                 ddres4=dreal(dcondchan(ichan))
                 if (ichan.eq.ichan0)  ddres4=dreal(dcurr2)*2*pi
                 write(85,*) tt,ichan,om,dreal(dcondchan(ichan)),ddres
     *            ,ddres1,ddres2,ddres3,ddres4
              enddo   
              write(85,*) " "
              write(85,*) " "
c           --broadened conductiviy to dcondeta
              do ichan2=1,nchan2
                 dcondeta(ichan2)=0.d0
              enddo !ichan2
              dom=(ommax-ommin)/nchan
              dom2=(ommax2-ommin2)/nchan2
c             spread each channel over om-20eta,om+20eta frequencies
              do ichan=1,nchan
               om=ommin+ichan*dom+dom/2
               ichan2min= (om-20*eta -ommin2)/dom2
               ichan2max= (om+20*eta -ommin2)/dom2
               if (ichan2min.lt.1) ichan2min=1
               if (ichan2min.gt.nchan2) ichan2min=nchan2
               if (ichan2max.lt.1) ichan2max=1
               if (ichan2max.gt.nchan2) ichan2max=nchan2
               do ichan2=ichan2min,ichan2max
                 om2=ommin2+ichan2*dom2+dom2/2
                 dcondeta(ichan2)=dcondeta(ichan2)
     *            +dreal(dcondchan(ichan))*brf(om-om2,eta,1) ! broadening; 1- gausian, 2 -lorenzian
               enddo !ichan2
              enddo !ichan
c           --broadended conductivity with D_c from tau
              do ichan2=1,nchan2
                 dcondeta1(ichan2)=dcondeta(ichan2)
              enddo
c             fix just part from om=0 chanel
              ichan=ichan0
              om=ommin+ichan*dom+dom/2
              ichan2min= (om-20*eta -ommin2)/dom2
              ichan2max= (om+20*eta -ommin2)/dom2
              if (ichan2min.lt.1) ichan2min=1
              if (ichan2min.gt.nchan2) ichan2min=nchan2
              if (ichan2max.lt.1) ichan2max=1
              if (ichan2max.gt.nchan2) ichan2max=nchan2
              do ichan2=ichan2min,ichan2max
                om2=ommin2+ichan2*dom2+dom2/2
                dcondeta1(ichan2)=dcondeta1(ichan2)
     *           +(2.d0*pi*dcstiff-dreal(dcondchan(ichan)))
     *           *brf(om-om2,eta,1) ! broadening; 1- gausian, 2 -lorenzian
              enddo !ichan2
              do ichan2=1,nchan2
                 dcondeta2(ichan2)=dcondeta(ichan2)
              enddo
c             fix just part from om=0 chanel
              ichan=ichan0
              om=ommin+ichan*dom+dom/2
              ichan2min= (om-20*eta -ommin2)/dom2
              ichan2max= (om+20*eta -ommin2)/dom2
              if (ichan2min.lt.1) ichan2min=1
              if (ichan2min.gt.nchan2) ichan2min=nchan2
              if (ichan2max.lt.1) ichan2max=1
              if (ichan2max.gt.nchan2) ichan2max=nchan2
              do ichan2=ichan2min,ichan2max
                om2=ommin2+ichan2*dom2+dom2/2
                dcondeta2(ichan2)=dcondeta2(ichan2)
     *           +(-dreal(dcurr2a)*4*pi/tt)
     *           *brf(om-om2,eta,1) ! broadening; 1- gausian, 2 -lorenzian
              enddo !ichan2
c           --broadended conductivity with subtracted <J>^2(om=0) and D_c(phase) 
c             ddres2=dreal(dcondchan(ichan))-dreal(dcurr2a)*4*pi/tt-dreal(dcurr2)*2*pi   
              do ichan2=1,nchan2
                 dcondeta3(ichan2)=dcondeta(ichan2)
              enddo
c             fix just part from om=0 chanel
              ichan=ichan0
              om=ommin+ichan*dom+dom/2
              ichan2min= (om-20*eta -ommin2)/dom2
              ichan2max= (om+20*eta -ommin2)/dom2
              if (ichan2min.lt.1) ichan2min=1
              if (ichan2min.gt.nchan2) ichan2min=nchan2
              if (ichan2max.lt.1) ichan2max=1
              if (ichan2max.gt.nchan2) ichan2max=nchan2
              do ichan2=ichan2min,ichan2max
                om2=ommin2+ichan2*dom2+dom2/2
                dcondeta3(ichan2)=dcondeta3(ichan2)
     *           +(-dreal(dcurr2a)*4*pi/tt-dreal(dcurr2)*2*pi)
     *           *brf(om-om2,eta,1) ! broadening; 1- gausian, 2 -lorenzian
              enddo !ichan2
c           --broadended conductivity with om=0 chanel set to 0
              do ichan2=1,nchan2
                 dcondeta4(ichan2)=dcondeta(ichan2)
              enddo
c             fix just part from om=0 chanel
              ichan=ichan0
              om=ommin+ichan*dom+dom/2
              ichan2min= (om-20*eta -ommin2)/dom2
              ichan2max= (om+20*eta -ommin2)/dom2
              if (ichan2min.lt.1) ichan2min=1
              if (ichan2min.gt.nchan2) ichan2min=nchan2
              if (ichan2max.lt.1) ichan2max=1
              if (ichan2max.gt.nchan2) ichan2max=nchan2
              do ichan2=ichan2min,ichan2max
                om2=ommin2+ichan2*dom2+dom2/2
                dcondeta4(ichan2)=dcondeta4(ichan2)
     *           +(-dreal(dcondchan(ichan)))
     *           *brf(om-om2,eta,1) ! broadening; 1- gausian, 2 -lorenzian
              enddo !ichan2
c           --broadended conductivity D_c from phase (dcurr2)
              do ichan2=1,nchan2
                 dcondeta5(ichan2)=dcondeta(ichan2)
              enddo
c             fix just part from om=0 chanel
              ichan=ichan0
              om=ommin+ichan*dom+dom/2
              ichan2min= (om-20*eta -ommin2)/dom2
              ichan2max= (om+20*eta -ommin2)/dom2
              if (ichan2min.lt.1) ichan2min=1
              if (ichan2min.gt.nchan2) ichan2min=nchan2
              if (ichan2max.lt.1) ichan2max=1
              if (ichan2max.gt.nchan2) ichan2max=nchan2
              do ichan2=ichan2min,ichan2max
                om2=ommin2+ichan2*dom2+dom2/2
                dcondeta5(ichan2)=dcondeta5(ichan2)
     *           +(2.d0*pi*dreal(dcurr2)-dreal(dcondchan(ichan)))
     *               *brf(om-om2,eta,1) ! broadening; 1- gausian, 2 -lorenzian
              enddo !ichan2
c           --write broadened conductivity (dcondeta) to file
              dom2=(ommax2-ommin2)/nchan2
              write(86,*) "#tt =",tt,"eta=",eta
              write(86,*) "#tt,ichan2,om2,dcondeta(ichan2),D_c(om=0)"
     *        ,"-<J>^2(om=0),-<J>^2(om=0)-D_c(pahse), om=0 set to 0,"
     *        ," D_c phase"
              do ichan2=1,nchan2
                 om2=ommin2+ichan2*dom2+dom2/2
                 write(86,*) tt,ichan2,om2,dcondeta(ichan2),
     *           dcondeta1(ichan2),dcondeta2(ichan2),dcondeta3(ichan2),
     *           dcondeta4(ichan2),dcondeta5(ichan2)
              enddo   
              write(86,*) " "
              write(86,*) " "
              
           endif
      enddo                     !itt
       nfile=nfile+1
       close(80)
       close(85) !condchan
       close(86) !condeta
       close(87) !tau expect val.
      enddo ! indf - end loop over all densities   
c-----------------------------------------------------------------------
        t0=dtime(ttime)
        write(0,*) "time writing res in files in sec.:", ttime(1)
        t0=etime(ttime)
        write(0,*) "time full calculation in seconds :", ttime(1)
c-----------------------------------------------------------------------
      end
c======================================================================
c=======================================================================
c     broadening function : gaussian - i=1 : lorenzian - i=2
      real(8) function brf(om,sig,i)
      implicit real(8) (a-h, o-y), complex(8) (z)
      implicit integer (i-n)

      if(i.ne.1.and.i.ne.2) then
         stop "ERORR in brf function. Not right choice."
      endif
      
      pi=4.d0*datan(1.d0)
      brf=0.d0
      
      if (i.eq.2) then
         brf=sig/(om**2+sig**2)/pi
      endif
      
      if(i.eq.1) then
         brf=1.d0/dsqrt(2*pi)/sig*dexp(-om**2/(2*sig**2))
      endif
       
      return
      end
c-----------------------------------------------------------------------
   
