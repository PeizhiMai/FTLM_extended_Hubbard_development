c=======================================================================
c
c     Hubbard 2D, T>0, optical conductivity, FTLM
c     triangular lattice
c     J.K 29.5.2014 - implem. optical conductivity
c=======================================================================
      program hubTri2Dcond
      use IFPORT     ! to use rand with ifort compiler
      implicit real*8 (a-h, o-z)
      implicit integer (i-n)
c-----------------------------------------------------------------------
c   PARAMETERS for triangular lattice:
c      lx1,ly1 =postition of the right-down corner of the romboid *
c      lx2,ly2 =position of the left-up conrner of the romboid *
c            * - in units of lattice vectors a1=(a,0)
c              and a2=(a/2,a*sqrt(3)/2)
c              lx1*ly2 must be higer than ly1*lx2   
c      nnf   =number of fermions (N_up + N_down)
c      nca   =maximum total number of configurations for system
c             with fixed number of holes when S_z is varied
c      npa   =maximum number of parent configurations for system
c             with fixed number of holes when S_z is varied
c-----------------------------------------------------------------------
c   since calculatin goes over all number of ferminos, parameter nnf
c    should be large (at least nn0)
c-----------------------------------------------------------------------
c  4 sites systems
c      parameter (lx1=2,ly1=0, lx2=0, ly2=2)
c      parameter (nnf=4,  nca=36, npa=12, rhf=1.0, rnop=1.0)
c      parameter (nnfq=4,  ncaq=36, npaq=12, rhfq=1.0 )
c-----------------------------------------------------------------------
c  6 sites systems
c      parameter (lx1=3,ly1=0, lx2=-1, ly2=2)   !frustrated
c      parameter ( nnf=6,  nca=400,  npa=70, rhf=1.0,rnop=1.0 )
c      parameter ( nnfq=6, ncaq=400, npaq=70, rhfq=1.0 )
c-----------------------------------------------------------------------
c   8 sites systems
      parameter (lx1=4,ly1=0, lx2=-1, ly2=2)   
      parameter ( nnf=8,  nca=4900,  npa=628, rhf=1.0,rnop=1.0 )
      parameter ( nnfq=8, ncaq=4900, npaq=628, rhfq=1.0 )
c-----------------------------------------------------------------------
c  12 sites systems (imperfection as in kent2005)
c      parameter (lx1=4,ly1=0, lx2=0, ly2=3)    
c      parameter (nnf=12,  nca=853776, npa=71300 , rhf=1.0, rnop=1.0)
c      parameter (nnfq=12,  ncaq=853776, npaq=71300, rhfq=1.0 )
c-----------------------------------------------------------------------
c  16 sites systems 
c      parameter (lx1=4,ly1=0, lx2=0, ly2=4) ! symet. imp=1 
c      parameter (nnf=16,nca=165636900, npa=10360000, rhf=0.6, rnop=1.0)
c      parameter (nnfq=16,ncaq=165636900, npaq=10360000, rhfq=0.6 )
c-----------------------------------------------------------------------
c   CALCULATED CONSTANTS: 
c      nn0   =number of sites in the romb
c      nho   =maximum number of hops per configuration (array size)
c-----------------------------------------------------------------------
      parameter (nn0=lx1*ly2-ly1*lx2, nho=6*nnf, nhoq=6*nnfq )
      parameter (nham0=npa*nho*rhf+1, nham0q=npaq*nhoq*rhfq+1 )
      parameter (noper0=npaq*nn0*rnop)
c-----------------------------------------------------------------------
c   LANCZOS AND WAVE FUNCTION CONSTANTS:
c      lstmax=maximum number of Lanczos steps
c-----------------------------------------------------------------------
      parameter (lstmax=200)
c-----------------------------------------------------------------------
c   VARIABLES: Indexation
c-----------------------------------------------------------------------
      integer ix(40),iy(40),nii(-10:10,-10:10),nni(40,6)
      integer ibin(40,40),nshft(40,40)
      integer*2 irels(nca),irelm(nca)
      integer irelp(nca)
      integer iuconfn(40),idconf(40)
      integer ipar(npa)
      integer*2 ideg(npa)
      common /iii/ nshft,ibin,n0,nni
      common /ixiy/ ix,iy
c-----------------------------------------------------------------------
c   VARIABLES: Hamiltonian, operators
c-----------------------------------------------------------------------
c      kx, ky - wave vector of 1st Lanczos (or ground state)
c      q - wave vector of operator
      real*8 kx,ky
      real*8 kall(40,2)
      integer hamhopp(nham0)
      integer iham(npa+1)
      integer hamhop(nham0)
      integer ihamm
      integer ia
      complex*16 hamhopm(nham0)
      real*8 dnf(npa)
c     number of double occupied sites
      integer*2 hamdia(npa)
      complex*16 phase(nn0)
c     complex (due to twisted boundary conditions) hopping parameters
c     postive indexes  (1:40)    -  up  spin hopping (new TBC)
c     negative indexes (-1:-40)  - down spin hopping (new TBC)
      complex*16 ctc(-40:40)
      complex*16 ctc0
      complex*16 phsum
c     list of hopping parameters
      real*8 ctl(40)
c-----------------------------------------------------------------------
c     bond order operator (to store direction of hopping)
      integer*2 hamhopdir(nham0)
c-----------------------------------------------------------------------
c   CURRENT OPERATOR
c-----------------------------------------------------------------------
c     curdir: current operator hopping direction prefactor:
c        1 if hopping in positive direciton of current
c       -1 if hopping in negative direciton of current
c       0.5 if hopping for half lattice distance in positive direciton of current
c      -0.5 if hopping for half lattice distance in negative direciton of current
      real*8 curdir(-6:6)
c-----------------------------------------------------------------------
c     SECOND LANCZOS
c-----------------------------------------------------------------------
c     SECOND LANCZOS in the same sector
c     for optical conductivity
      real*8 b(0:lstmax-1)
      real*8 bb(-2:lstmax-1)
      complex*16 phib(0:npa,0:lstmax+1)
      real*8 db(lstmax),eb(lstmax)
      real*8 enb(lstmax)
      real*8 vkb(lstmax,lstmax)
      complex*16 crvkavib8(lstmax,lstmax) !current oper. (second lanc|J|first lanc)
      complex*16 crvkavib10(lstmax) !current oper. diag.(first lanc|J|first lanc)
      complex*16 cvkb1e1(lstmax) !explicit vkb(1) (symetrization of dynamic cal.)
c-----------------------------------------------------------------------
c   TAU (KINETIC ENERGY LIKE TERM) OPERATOR
c-----------------------------------------------------------------------
c     taudir: curdir**2 current operator hopping direction prefactor:
      real*8 taudir(-6:6)
      complex*16 crvkavib9(lstmax,lstmax) !tau operator
c-----------------------------------------------------------------------
c   VARIABLES: Wave functions, Lanczos
c-----------------------------------------------------------------------
      real*8 a(0:lstmax-1)
      real*8 da(lstmax),ea(lstmax)
      complex*16 phia(0:npa,0:lstmax+1)
      logical inromb
      real*8 ba(-2:lstmax-1)
      real*8 vka(lstmax,lstmax)
      real*8 en(lstmax)
c     vector used only in subroutine lanstep, but must be translated
c     there form main program, because of memory (16 site sistem)
      complex*16 psix(npa)
      complex*16 psi(npa)
      complex*16 crkji,crkji1
c-----------------------------------------------------------------------
c   VARIABLE: STATIC, THERMODYNAMIC
c-----------------------------------------------------------------------
c     name of the data file for results
      character ndf*15
c     number of idk or file id for data file name
      character fid*3 !for idk for up-spins
      character fidd*3 !for idk for down-spins
      character fidsmp*3 !for ismp in data filename
c     list of phases in units of reciprocal wave vectors (b1,b2)
      real*8 dklist(512,2)
c     folder name for data files 
      character fldr*100
c-----------------------------------------------------------------------
c   TIMING
c-----------------------------------------------------------------------
      real t0,tt0,ttime(2)
c-----------------------------------------------------------------------
c OPEN MP VARIABLES
c-----------------------------------------------------------------------
      integer id,nthrds
      INTEGER OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
c-----------------------------------------------------------------------      
c   DATA INITIALIZATION
c   ct - hopping parameter in (a,0) and (a/2,a*sqrt(3)/2) directions
c   ct1 - hopping in diagonal direction (-a/2,a*sqrt(3)/2)
c   cu - onsite repulsion paramater u - interaction
c   mq - index of wave vect of operator
c   eps1 - criterium for cancled parents and convergence of lanczos (full space)
c   eps2 - criterium for reality of lanczos parameters a and for
c          diffrences in b_k
c        - criterium for convergence of gs energy (lanczos)
c   nsmp - number of samples for FTLM 
c    ndkl - number of sampling over dphi or dkx and dky shifts (from list),
c          total number of samples 
c-----------------------------------------------------------------------
      data ct, ct1/ 1.0d0,  0.0d0/
      data cu /12.0d0/
      data eps1, eps2/ 1.d-7, 1.d-11/
      data ndf /'n08u120tp00.dat'/
      data nsmp /4/ !should be 1 for use of ismpmin,ismpmax
      data ndkl /2/ !number of phases from list (1,2,4,8,16,32,64,128)
      data idkmin, idkmax /1,2/ !TBC phases
      data ismpmin, ismpmax /1,4/ !samples over random vectors.
c-----------------------------------------------------------------------
c=======================================================================
c=======================================================================
c   BEGIN MAIN
c=======================================================================
      pi=4.d0*datan(1.d0)
      n0=nn0
      nca0=nca+2
c-----------------------------------------------------------------------
c     reseting random number series according to first idk,ismpmin
      call srand((idkmin*512+ismpmin))
c-----------------------------------------------------------------------
      write(0,*) "data file name ", ndf      
      write(0,*) 'nn0 ct ct1 nsmp ndkl ', nn0, ct, ct1, nsmp, ndkl
!$OMP PARALLEL PRIVATE(id,nthrds)
      nthrds = OMP_GET_NUM_THREADS()
      id = OMP_GET_THREAD_NUM()
      write (0,*) '---MP: Tread ID', id, nthrds
!$OMP END PARALLEL            
      do idk=idkmin,idkmax
       idkd=idk
       do ismp=ismpmin,ismpmax
c       clearing (deleting) old files
       iunit=20
       write(fid,100) idk/100,(idk-idk/100*100)/10,(idk-idk/10*10)
       write(fidsmp,100) ismp/100,(ismp-ismp/100*100)/10,
     *   (ismp-ismp/10*10)
 100   format(3i1)
c     local folder for raw data
      write(fldr,*)'.'
c     ifort bufferd option
       open(unit=iunit,file=trim(adjustl(fldr))//'/hubTri2Dconddata/'
     *    //ndf//'_idk'//fid//'_ismp'//fidsmp,
     *    status='unknown',access='sequential',form='unformatted',
     *    buffered='yes')
       write(0,*) trim(adjustl(fldr))//'/hubTri2Dconddata/'
     *    //ndf//'_idk'//fid//'_ismp'//fidsmp
       close(iunit,status='delete')
       enddo !ismp
      enddo !idk
c=======================================================================
c   SITE NUMERATION
c     site i has coordinates ix(i),iy(i);  site number nii(ix,iy);
c-----------------------------------------------------------------------
      isite=1
      do iix=min(0,lx1,lx2,lx1+lx2),max(0,lx1,lx2,lx1+lx2)
        do iiy=min(0,ly1,ly2,ly1+ly2),max(0,ly1,ly2,ly1+ly2)
c   point within the romb ?
           if (inromb(iix,iiy,lx1,ly1,lx2,ly2)) then
            nii(iix,iiy)=isite
            ix(isite)=iix
            iy(isite)=iiy
            isite=isite+1
      	  endif
        enddo
      enddo
c.......................................................................
      write(0,*) 'Site numbering:  i, x, y, i(x,y) '
      write(0,1000) (i,ix(i),iy(i),nii(ix(i),iy(i)),i=1,n0)
1000  format(4(i6))
      if ((isite-1).ne.n0)
     * write(0,*)'ERROR: wrong number of sites in romb', isite, n0
c-----------------------------------------------------------------------
c   NEAREST NEIGHBOURS:
c     nearest neighbours of site i in
c     (a,0),(a/2,a*sgrt(3)/2),(-a/2,a*sgrt(3)/2),(-a,0),
c      (-a/2, -a*sqrt(3)/2),(a/2,-a*sqrt(3)/2) direction
c     are nni(i,1), nni(i,2), nni(i,3), nni(i,4), nni(i,5)
c      and nni(i,6) respectively. 
c-----------------------------------------------------------------------
      do i=1,n0
c   loop through directions
        do j=1,6
c  moves for every direction
           if(j.eq.1) then
              mx=1
              my=0
             else if (j.eq.2) then
                mx=0
                my=1
             else if (j.eq.3) then
                mx=-1
                my=1
             else if (j.eq.4) then
                mx=-1
                my=0
             else if (j.eq.5) then
                mx=0
                my=-1
             else if (j.eq.6) then
                mx=1
                my=-1
             endif              
c   check 9 neighbouring squares for neighbour
          if (inromb(ix(i)+mx,iy(i)+my,lx1,ly1,lx2,ly2)) then
            iix=ix(i)+mx
            iiy=iy(i)+my           
           else if (inromb(ix(i)+mx+lx1,iy(i)+my+ly1,lx1,ly1,lx2,ly2))
     *      then
            iix=ix(i)+mx+lx1
            iiy=iy(i)+my+ly1       
           else if (inromb(ix(i)+mx-lx2,iy(i)+my-ly2,lx1,ly1,lx2,ly2))
     *      then
            iix=ix(i)+mx-lx2
            iiy=iy(i)+my-ly2        
           else if (inromb(ix(i)+mx-lx1,iy(i)+my-ly1,lx1,ly1,lx2,ly2))
     *      then
            iix=ix(i)+mx-lx1
            iiy=iy(i)+my-ly1        
           else if (inromb(ix(i)+mx+lx2,iy(i)+my+ly2,lx1,ly1,lx2,ly2))
     *      then
            iix=ix(i)+mx+lx2
            iiy=iy(i)+my+ly2        
           else if
     *      (inromb(ix(i)+mx+lx1+lx2,iy(i)+my+ly1+ly2,lx1,ly1,lx2,ly2))
     *      then
            iix=ix(i)+mx+lx1+lx2
            iiy=iy(i)+my+ly1+ly2        
           else if
     *      (inromb(ix(i)+mx-lx1-lx2,iy(i)+my-ly1-ly2,lx1,ly1,lx2,ly2))
     *      then
            iix=ix(i)+mx-lx1-lx2
            iiy=iy(i)+my-ly1-ly2        
           else if
     *      (inromb(ix(i)+mx+lx1-lx2,iy(i)+my+ly1-ly2,lx1,ly1,lx2,ly2))
     *      then
            iix=ix(i)+mx+lx1-lx2
            iiy=iy(i)+my+ly1-ly2
           else if
     *      (inromb(ix(i)+mx-lx1+lx2,iy(i)+my-ly1+ly2,lx1,ly1,lx2,ly2))
     *      then
            iix=ix(i)+mx-lx1+lx2
            iiy=iy(i)+my-ly1+ly2        
      	  endif
          nni(i,j)=nii(iix,iiy)
c          write(0,*)' i, j, mx, my , ', i, j, mx, my 
       enddo
      enddo
c-----------------------------------------------------------------------
c      write(0,*) 'Nearest neighbours: i, nn(1-6)'
c      write(0,1021) ((nni(i,j),j=1,6),i=1,27)
c1021  format(6(i6))
c-----------------------------------------------------------------------
c   SHIFT MAPPINGS:
c     shift s is defined as one which brings site 1 to site s;
c     for shift s, site j is brought into site nshft(s,j)
c-----------------------------------------------------------------------
      do i=1,n0
        do j=1,n0
          nshft(i,j)=j
        enddo
        isx=-ix(i)+ix(1)
        isy=-iy(i)+iy(1)
        if (isx.ge.0) then
          do j=1,isx
            do k=1,n0
              nshft(i,k)=nni(nshft(i,k),4)
            enddo
          enddo
         else
          do j=1,-isx
            do k=1,n0
              nshft(i,k)=nni(nshft(i,k),1)
            enddo
          enddo
        endif
        if (isy.ge.0) then
          do j=1,isy
            do k=1,n0
              nshft(i,k)=nni(nshft(i,k),5)
            enddo
          enddo
        else
          do j=1,-isy
            do k=1,n0
              nshft(i,k)=nni(nshft(i,k),2)
            enddo
          enddo
        endif
      enddo    
c      write(0,*)'nshift'
c      write(0,*)(nshft(2,i),i=1,n0)
c      do is=1,n0
c         do i=1,n0
c            write(0,*) "shift", is, i," -> ",nshft(is,i)
c     *        ," dx dy ",ix(i)-ix(nshft(is,i)),iy(i)-iy(nshft(is,i))   
c         enddo !i
c      enddo !is
c-----------------------------------------------------------------------
c   BINOMIAL COEFICIENTS:
c      sums of binomial symbols for indexation
c   ibin(i,j) is not binomial koeficient (i;j)
c             it holds that  ibin(i,j)= (j-2+i;j-2)
c             but i and j must be less than n0 and positive 
c-----------------------------------------------------------------------
      do i=1,n0
        ibin(i,1)=0
      enddo
      do j=1,n0
        ibin(1,j)=j-1
      enddo
      do i=2,n0
        do j=2,n0
          ibin(i,j)=ibin(i-1,j)+ibin(i,j-1)
        enddo
      enddo
c        write(0,*)'ibin',((ibin(i,j),i=1,n0),j=1,n0)      
c-----------------------------------------------------------------------
c     ALLOWED WAVE VECTORS CALCULATION
      call  wavevec(lx1,ly1,lx2,ly2,n0,kall)
      do i=1, n0
      write(0,*)i, kall(i,1),kall(i,2)
      enddo
c-----------------------------------------------------------------------
c     setting phase shift in the list by priority (dqa,dqb) [0:1)
      dklist(1,1)=0.d0
      dklist(1,2)=0.d0
      dklist(2,1)=0.5d0
      dklist(2,2)=0.5d0
      dklist(3,1)=0.5d0
      dklist(3,2)=0.d0
      dklist(4,1)=0.d0
      dklist(4,2)=0.5d0
      do ii=5,8
        dklist(ii,1)=dklist(ii-4,1)+0.25d0
        dklist(ii,2)=dklist(ii-4,2)+0.25d0
      enddo
      do ii=9,16
        dklist(ii,1)=dklist(ii-8,1)+0.25d0
        dklist(ii,2)=dklist(ii-8,2)
        if (dklist(ii,1).eq.1.d0) dklist(ii,1)=dklist(ii,1)-1.d0 
      enddo
      do ii=17,32
        dklist(ii,1)=dklist(ii-16,1)+0.125d0
        dklist(ii,2)=dklist(ii-16,2)+0.125d0
        if (dklist(ii,1).eq.1.d0) dklist(ii,1)=dklist(ii,1)-1.d0 
        if (dklist(ii,2).eq.1.d0) dklist(ii,2)=dklist(ii,2)-1.d0 
      enddo
      do ii=33,64
        dklist(ii,1)=dklist(ii-32,1)+0.125d0
        dklist(ii,2)=dklist(ii-32,2)
        if (dklist(ii,1).eq.1.d0) dklist(ii,1)=dklist(ii,1)-1.d0 
        if (dklist(ii,2).eq.1.d0) dklist(ii,2)=dklist(ii,2)-1.d0 
      enddo
      do ii=65,128
        dklist(ii,1)=dklist(ii-64,1)+1.d0/16
        dklist(ii,2)=dklist(ii-64,2)+1.d0/16
        if (dklist(ii,1).eq.1.d0) dklist(ii,1)=dklist(ii,1)-1.d0 
        if (dklist(ii,2).eq.1.d0) dklist(ii,2)=dklist(ii,2)-1.d0 
      enddo
      do ii=129,256
        dklist(ii,1)=dklist(ii-128,1)+1.d0/16
        dklist(ii,2)=dklist(ii-128,2)
        if (dklist(ii,1).eq.1.d0) dklist(ii,1)=dklist(ii,1)-1.d0 
        if (dklist(ii,2).eq.1.d0) dklist(ii,2)=dklist(ii,2)-1.d0 
      enddo
      do ii=257,512
        dklist(ii,1)=dklist(ii-256,1)+1.d0/32
        dklist(ii,2)=dklist(ii-256,2)+1.d0/32
        if (dklist(ii,1).eq.1.d0) dklist(ii,1)=dklist(ii,1)-1.d0 
        if (dklist(ii,2).eq.1.d0) dklist(ii,2)=dklist(ii,2)-1.d0 
      enddo
      write(0,*) "Phase shifts in units of reciprocal vectors b1,b2"
      do i=1,ndkl
        write(0,*)  i, (dklist(i,j),j=1,2)
      enddo
c-----------------------------------------------------------------------
c     list of hopping parameters for each direction of nearest neighbors
c-----------------------------------------------------------------------
      ctl(1)=-ct
      ctl(2)=-ct
      ctl(3)=-ct1
      ctl(4)=-ct
      ctl(5)=-ct
      ctl(6)=-ct1
c=======================================================================
c     LOOP OVER nf
c-----------------------------------------------------------------------
      nfmin=1
      nfmax=2*nn0
      do nf=nfmin,nfmax,1
      write(0,*)'==================================================='
      write(0,*)'=============  NF   ',nf,'   ========================='
      write(0,*)'==================================================='
c=======================================================================
      t0=dtime(ttime)
      write(0,2003) nn0,nf
 2003 format('  nn0 =',i3,'  nf = ',i3 )
      write(0,2010) ct,cu,ct1
 2010 format('  ct =',f7.2,'   cu =',f7.2,'   ct1 =',f7.2 )
c=======================================================================
      if(nf.le.n0) then
         numin=0
         numax=nf/2
      else
         numin=nf-n0
         numax=nf/2
      endif
c=======================================================================
c     LOOP OVER nu - number of up spins
c-----------------------------------------------------------------------
      do nu=numin,numax
         nd=nf-nu
         write(0,*)'nu,nd',nu,nd
c--------------------------------------------------------------------
         call parham(nu,nd,nca,npa,np,ncu,irels,irelm,irelp,
     *             ideg,ipar,nham0,nham,iham,hamhop,hamdia,
     *             hamhopdir)
c=======================================================================
c     LOOP OVER ik - wave vector
c-----------------------------------------------------------------------
         ikmin=1
         ikmax=nn0
         do ik=ikmin,ikmax
c          write(0,*)'Wave vector loop ik', ik
c-----------------------------------------------------------------------
c   select allowed wave vector for particular grid
c-----------------------------------------------------------------------
         kx=kall(ik,1)
         ky=kall(ik,2)
c         write(0,*)'nu, nd, ik, kx, ky ',nu, nd, ik, kx, ky
c-----------------------------------------------------------------------
c   WAVE VECTOR DEPENDENT PART:
c     the phases for each of the shifts k
c-----------------------------------------------------------------------
          do is=1,n0
           p=(ix(is)-ix(1))*kx+(iy(is)-iy(1))*(kx/2+ky*sqrt(3d0)/2)
           phase(is)=dcmplx(dcos(p),dsin(p))
          enddo
c          write(0,*)'***ik,kx,ky',ik,kx,ky
c          write(0,*)'phase',phase
c-----------------------------------------------------------------------
c   FIND CANCELED PARENT WAVE-VECTOR STATES:
c     check for destructive interference 
c-----------------------------------------------------------------------
          npw=np
          do i=1,np
c   Only parent states that translate into itself can be canceled.
           if (ideg(i).gt.1) then
             phsum=dcmplx(0.d0,0.d0)
c   Try all shifts
             do is=1,n0
c   Sum only over those shifts that duplicate the parent
               indx=ishift(ipar(i),is,iperm,nu,nd,ncu)
               if (indx.eq.ipar(i)) phsum=phsum+iperm*phase(is)
             enddo
c   Inverse of the normalization factor 
             dnf(i)=abs(phsum)/dsqrt(dble(ideg(i)))
             if (dnf(i).lt.eps1) npw=npw-1
           else
             dnf(i)=1.d0
           endif
          enddo !i
c     if all parents wave vector states cancel go to next ik
c     this is posible if nf=0, or nf=2*n0
          if (npw.eq.0) then
             goto 7110
          endif
c.......................................................................
c       basic reciprocal vectors
        bk1x=2*pi*ly2/n0
        bk1y=2*pi*(-2*lx2-ly2)/(sqrt(3d0)*n0)
        bk2x=-2*pi*ly1/n0
        bk2y=2*pi*(2*lx1+ly1)/(sqrt(3d0)*n0)
c=======================================================================
c     loop over idk for up and down spins
c=======================================================================          
       do idk=idkmin,idkmax !loop over TBC for up-spins
          idkd=idk  !for same phase for ud and down 
c.......................................................................
c        up-spin phase
         dqa=dklist(idk,1)
         dqb=dklist(idk,2)
         dkx=(dqa)*bk1x+(dqb)*bk2x 
         dky=(dqa)*bk1y+(dqb)*bk2y
c        down-spin phase
         dqa=dklist(idkd,1)
         dqb=dklist(idkd,2)
         dkxd=(dqa)*bk1x+(dqb)*bk2x 
         dkyd=(dqa)*bk1y+(dqb)*bk2y
c         do i=1, n0
c          write(0,*)'aaa', i, kall(i,1)+dkx,kall(i,2)+dky
c         enddo
c.......................................................................
c        up spins
         fixi=dkx*pi
         fiyi=dky*pi
c        down spins
         fixid=dkxd*pi
         fiyid=dkyd*pi
         write(0,*) 'dkx, dky = ', dkx, dky
c=======================================================================
c-----------------------------------------------------------------------
c   Hamiltonian k
c-----------------------------------------------------------------------
c     Hopping part, phase for up spin fermions (TBC)
         ctc(1)=dcmplx(dcos(fixi),dsin(fixi))
         ctc(2)=dcmplx(dcos(fixi/2+fiyi*sqrt(3d0)/2)
     *                      ,dsin(fixi/2+fiyi*sqrt(3d0)/2)) 
         ctc(3)=dcmplx(dcos(-fixi/2+fiyi*sqrt(3d0)/2)
     *                      ,dsin(-fixi/2+fiyi*sqrt(3d0)/2))
         ctc(4)=dcmplx(dcos(-fixi), dsin(-fixi))
         ctc(5)=dcmplx(dcos(-fixi/2-fiyi*sqrt(3d0)/2)
     *                     ,dsin(-fixi/2-fiyi*sqrt(3d0)/2))
         ctc(6)=dcmplx(dcos(fixi/2-fiyi*sqrt(3d0)/2)
     *                     ,dsin(fixi/2-fiyi*sqrt(3d0)/2))
c     Hopping part, phase for down spin ferminons (TBC)
         ctc(-1)=dcmplx(dcos(fixid),dsin(fixid))
         ctc(-2)=dcmplx(dcos(fixid/2+fiyid*sqrt(3d0)/2)
     *                      ,dsin(fixid/2+fiyid*sqrt(3d0)/2)) 
         ctc(-3)=dcmplx(dcos(-fixid/2+fiyid*sqrt(3d0)/2)
     *                      ,dsin(-fixid/2+fiyid*sqrt(3d0)/2))
         ctc(-4)=dcmplx(dcos(-fixid), dsin(-fixid))
         ctc(-5)=dcmplx(dcos(-fixid/2-fiyid*sqrt(3d0)/2)
     *                     ,dsin(-fixid/2-fiyid*sqrt(3d0)/2))
         ctc(-6)=dcmplx(dcos(fixid/2-fiyid*sqrt(3d0)/2)
     *                     ,dsin(fixid/2-fiyid*sqrt(3d0)/2))

!$OMP    PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ih,isig,idir,ia)         
         do i=1,np
           do ih=iham(i),iham(i+1)-1
c   index of the configuration
             if(hamhop(ih).ge.0) then
                isig=1
             else
                isig=-1
             endif
             idir=abs(hamhopdir(ih))
             if (idir.lt.1.or.idir.gt.6) then
             write(0,*) "idir error: idir, hamhop(ih) ",idir,hamhop(ih)
              stop ' error: stoping program due to wrong idir'
             endif                                
             ia=isig*hamhop(ih)
c     hopping matrix element between parents <hamhopp(i,l)| H_t |i>
             if (ia.le.0) then
                stop 'ERROR:ind. of conf.< 1 (ia=0 or <0)'
             endif
             if (dnf(i).gt.eps1.and.ia.ne.0) then
               hamhopp(ih)=irelp(ia)
c              hop matrix elemet without constants ct and phase due to
c              TBC (in ctc and ctl)
               hamhopm(ih)=isig*irelm(ia)*
     *           dnf(irelp(ia))/dnf(i)*phase(irels(ia))
             else
               hamhopp(ih)=0
               hamhopm(ih)=dcmplx(0.d0,0.d0)
             endif
           enddo ! ih
cc          write(0,*)"hamhopp(ih)",hamhopp(ih)
          enddo ! i
!$OMP     END  PARALLEL DO           
c=======================================================================
c     LOOP OVER RANDOM SAMPLES
c-----------------------------------------------------------------------
          write(0,*) "Starting loop over random samples"
          do ismp=ismpmin,ismpmax
           write(0,*)'Sample ', ismp,' ik ',ik
c-----------------------------------------------------------------------
c   Random state (non-normalized)
c-----------------------------------------------------------------------
           do i=1,np
            if (dnf(i).gt.eps1) then
             ran2=dble(rand(0))
             ran1=dble(rand(0))
             phia(i,1)=dcmplx(ran2-0.5d0,ran1-0.5d0)
            else
             phia(i,1)=dcmplx(0.d0,0.d0)
            endif
           enddo ! i
c-----------------------------------------------------------------------
c   Lanczos on random state (wave vector k)
c-----------------------------------------------------------------------
          en(1)=1.d4
          call lanczos(hamdia,hamhopp,hamhopm,hamhopdir,
     *                iham,nham0,nham,
     *                phia,a,ba,np,npa,lstmax,ldima,
     *                eps1,eps2,da,ea,en,vka,psix,cu,ctc,ctl)
          nex=min(ldima,lstmax)
c-----------------------------------------------------------------------
c   Operator on random state in this sym. sector 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c    OPTICAL CONDUCTIVITY
c-----------------------------------------------------------------------
c       second lanczos from J on random vector (J-current operator)
c       operator (current oper.) on random vector
        write(0,*) "   calculating optical cond. - sec. Lan"
c     projection of hoppings to current direction
c     defines the current direction
        curdir(1)=1.d0 
        curdir(2)=0.5d0
        curdir(3)=-0.5d0
        curdir(4)=-1.d0
        curdir(5)=-0.5d0
        curdir(6)=0.5d0
        curdir(-1)=1.d0 
        curdir(-2)=0.5d0
        curdir(-3)=-0.5d0
        curdir(-4)=-1.d0
        curdir(-5)=-0.5d0
        curdir(-6)=0.5d0
        curdir(0)=0.d0
c     square lattice test (in x+y direction for better averaging)
        if(ct1.eq.0.d0) then
        curdir(1)=dsqrt(2.d0)/2.d0 
        curdir(2)=dsqrt(2.d0)/2.d0 
        curdir(3)=0.0d0
        curdir(4)=-dsqrt(2.d0)/2.d0 
        curdir(5)=-dsqrt(2.d0)/2.d0 
        curdir(6)=0.0d0
        curdir(-1)=dsqrt(2.d0)/2.d0  
        curdir(-2)=dsqrt(2.d0)/2.d0 
        curdir(-3)=0.0d0
        curdir(-4)=-dsqrt(2.d0)/2.d0 
        curdir(-5)=-dsqrt(2.d0)/2.d0 
        curdir(-6)=0.0d0
        curdir(0)=0.d0
        endif
c     current operator on lanczos state phia(:,k)
c         |psi>=A |phia(:,k)>, A = current operator
         k=1
         call opercur(hamhopp,hamhopm,hamhopdir,
     *                   iham,nham0,nham,
     *                   phia,k,np,npa,lstmax,psi,
     *                   ctc,curdir)
c     second lanczos
         do i=1,np
         if (dnf(i).gt.eps1) then
          phib(i,1)=psi(i)
         else
          phib(i,1)=dcmplx(0.d0,0.d0)
         endif
        enddo ! i
        enb(1)=1.d4
        call lanczos(hamdia,hamhopp,hamhopm,hamhopdir,
     *                iham,nham0,nham,
     *                phib,b,bb,np,npa,lstmax,ldimb,
     *                eps1,eps2,db,eb,enb,vkb,psix,cu,ctc,ctl)
        nexb=min(ldimb,lstmax)
c        write(0,*) "Second lanczos nexb",nexb
c-----------------------------------------------------------------------
c     current oper. matrix elements between first and second lanczos basis 
        write(0,*) "   calculating current oper.-mat.elem"
        do ja=1,nex
          do kb=1,nexb
              crvkavib8(kb,ja)=dcmplx(0.d0,0.d0)
          enddo      
        enddo
c     <firs.lanc.|j|sec.lanc> just diagonal
        do ka=1,nex
            crvkavib10(ka)=dcmplx(0.d0,0.d0)
        enddo      
c     loop over all lanczos basis states phia(:,k)
        do k=1,ldima
c          operator (current) on lanczos state phia(:,k)
         call opercur(hamhopp,hamhopm,hamhopdir,
     *                   iham,nham0,nham,
     *                   phia,k,np,npa,lstmax,psi,
     *                   ctc,curdir)
         do l=1,ldimb
c           matrix element between <phib(:,l)|A|phia(:,k)>
            crkji=dcmplx(0.d0,0.d0)
!$OMP       PARALLEL DEFAULT(SHARED) PRIVATE(i,crkji1)
            crkji1=dcmplx(0.d0,0.d0)
!$OMP       DO 
            do i=1,np
              crkji1=crkji1+dconjg(phib(i,l))*psi(i)
            enddo
!$OMP       END DO
!$OMP       ATOMIC
            crkji=crkji+crkji1
!$OMP       END PARALLEL
            
            do kb=1,nexb
               do ja=1,nex
                  crvkavib8(kb,ja)=crvkavib8(kb,ja)+
     *              vkb(l,kb)/bb(l-2)*crkji*vka(k,ja)/ba(k-2)
              enddo ! ja
            enddo ! kb
          enddo ! l

c----- start <first lanczos|j|first lanczos> (just some test)
         do l=1,ldima
c           matrix element between <phia(:,l)|A|phia(:,k)>
            crkji=dcmplx(0.d0,0.d0)
c            do i=1,np
c              crkji=crkji+dconjg(phib(i,l))*psi(i)
c            enddo
            
!$OMP       PARALLEL DEFAULT(SHARED) PRIVATE(i,crkji1)
            crkji1=dcmplx(0.d0,0.d0)
!$OMP       DO 
            do i=1,np
              crkji1=crkji1+dconjg(phia(i,l))*psi(i)
            enddo
!$OMP       END DO
!$OMP       ATOMIC
            crkji=crkji+crkji1
!$OMP       END PARALLEL
            
            do ka=1,nex
                crvkavib10(ka)=crvkavib10(ka)+
     *            vka(l,ka)/ba(l-2)*crkji*vka(k,ka)/ba(k-2)
            enddo ! ka
          enddo ! l
c----- end <first lanczos|j|first lanczos> 
          
       enddo! k            
c-----------------------------------------------------------------------
c     TAU or KINETIC ENERGY LIKE OPERATOR
c-----------------------------------------------------------------------
c   Transformation from Lanczos to wave functions and matrix elements
c-----------------------------------------------------------------------
c     tau operator: direction of 1st nearest neighbour (idir=1)
c-----------------------------------------------------------------------
         write(0,*) "   calculating tau oper. - kinetic energy like"
c     each direction multiplied with taudir(idir) which is related to
c         curdir: taudir=curdir**2
        do idir=-6,6
           taudir(idir)=(curdir(idir))**2.d0
        enddo
c     matrix element beween eigenstates for tau operator
         do ia=1,nex
          do ja=1,nex
            crvkavib9(ia,ja)=dcmplx(0.d0,0.d0)
          enddo      
        enddo
        if(ldima.ne.nex) then
           write(0,*) "error nex ldim ne", ldima,nex
           stop
        endif
c     loop over all lanczos basis states phia(:,k)
        do k=1,ldima
c          operator (tau) on lanczos state phia(:,k)
         call opertau(hamhopp,hamhopm,hamhopdir,
     *                   iham,nham0,nham,
     *                   phia,k,np,npa,lstmax,psi,
     *                   ctc,taudir)
         do l=1,ldima
c           matrix element between <phia(:,l)|A|phia(:,k)>
            crkji=dcmplx(0.d0,0.d0)

!$OMP       PARALLEL DEFAULT(SHARED) PRIVATE(i,crkji1)
            crkji1=dcmplx(0.d0,0.d0)
!$OMP       DO 
            do i=1,np
              crkji1=crkji1+dconjg(phia(i,l))*psi(i)
            enddo
!$OMP       END DO
!$OMP       ATOMIC
            crkji=crkji+crkji1
!$OMP       END PARALLEL
            
            do ia=1,nex
c               crkji=crkji*vka(l,ia)/ba(l-2) !?
               do ja=1,nex
                  crvkavib9(ia,ja)=crvkavib9(ia,ja)+
     *              vka(l,ia)/ba(l-2)*crkji*vka(k,ja)/ba(k-2)
              enddo ! ja
            enddo ! ia
          enddo ! l
        enddo ! k
c-----------------------------------------------------------------------
c   WRITING SYM. SECTOR DATA TO FILE
c-----------------------------------------------------------------------
          iunit=20
          write(fid,8000) idk/100,(idk-idk/100*100)/10,(idk-idk/10*10)
          write(fidsmp,8000) ismp/100,(ismp-ismp/100*100)/10,
     *     (ismp-ismp/10*10)
 8000     format(3i1)
          write(0,*)'   writing to file ','./hubTri2Dconddata/'
     *    //ndf//'_idk'//fid//'_ismp'//fidsmp
c     local data folder
      write(fldr,*)'.'
      open(unit=iunit,file=trim(adjustl(fldr))//'/hubTri2Dconddata/'
     *    //ndf//'_idk'//fid//'_ismp'//fidsmp,
     *    status='unknown',access='append',form='unformatted',
     *    buffered='yes')

          kx=kall(ik,1)
          ky=kall(ik,2)
          write(iunit) nn0, nf, nu
          write(iunit) ik, kx, ky, dkx, dky, dkxd, dkyd, ismp     
          write(iunit) npw, nsmp, ndkl
c-----------------------------------------------------------------------
c   WRITING RESULTS TO FILE (SAVING)
c-----------------------------------------------------------------------
c         saving vka(1,i) and en(i)
          write(iunit) nex
          do i=1, nex
             write(iunit) en(i), vka(1,i)
          enddo
          write(iunit) ba(-1)
c-----------------------------------------------------------------------
c         storing optical conductivity parameters
          write(iunit) nexb
          write(iunit) bb(-1)
          write(iunit) (enb(i),vkb(1,i),i=1,nexb)          
          write(iunit) ((crvkavib8(kb,ja),kb=1,nexb),ja=1,nex)
          write(iunit) (crvkavib10(ka),ka=1,nex)
C-----------------------------------------------------------------------
C         Storing tau oper (like kin en) matrix elements (upper tridiagonal)
c          write(iunit) ((dreal(crvkavib9(ia,ja)),ja=1,ia),ia=1,nex)
          write(iunit) ((dreal(crvkavib9(ia,ja)),ja=1,nex),ia=1,nex)
c-----------------------------------------------------------------------
          close(iunit)
c-----------------------------------------------------------------------
       t0=dtime(ttime)
       write(0,*)'   gscal nu,nd,ik,idk,idkd,en0,en0b,nex,nexb,time',
     *     nu,nd,ik,idk,idkd,en(1),enb(1),nex,nexb,t0-tt0
       tt0=t0
c=======================================================================
c     END OF LOOP OVER SAMPLES ismp
c-----------------------------------------------------------------------
        enddo ! ismp
c=======================================================================
c     END OF LOOP OVER (dkx,dky) (idk) [phase]
c-----------------------------------------------------------------------
        enddo ! idk for up-spin phase (in list dklist)
 7110   continue          
c=======================================================================
c     END WAVE VECTOR LOOP
c-----------------------------------------------------------------------
       enddo ! ik
c=======================================================================
c     END OF nu LOOP 
c-----------------------------------------------------------------------
      enddo ! nu loop
c=======================================================================
c     END OF LOOP OVER nf
c-----------------------------------------------------------------------
       t0=dtime(ttime)
       write(0,*)'                     TIME: nf ',nf,' loop',ttime(1)
      enddo !nf
c=======================================================================
      End
c=======================================================================
c   END MAIN
c=======================================================================
c=======================================================================
c=======================================================================
c   Find if point lies within the romb
c-----------------------------------------------------------------------
      logical function inromb(iix,iiy,lx1,ly1,lx2,ly2)
      implicit real*8 (a-h, o-z)
      implicit integer (i-n)
      x=iix+(lx1+lx2)/1000d0
      y=iiy+(ly1+ly2)/1000d0
      d1=lx1*y-ly1*x
      d2=(x-lx2)*ly1-(y-ly2)*lx1
      d3=x*ly2-y*lx2
      d4=-(x-lx1)*ly2+(y-ly1)*lx2
      inromb=d1.gt.0.and.d2.gt.0.and.d3.gt.0.and.d4.gt.0
      return
      end
c=======================================================================
c     FINDS ALL ALLOWED WAVE VECTORS IN 1BC
c-----------------------------------------------------------------------
      subroutine wavevec(lx1,ly1,lx2,ly2,n0,kall)
      implicit real*8 (a-h, o-z)
      implicit integer (i-n)
      real*8 kall(40,2)
c-----------------------------------------------------------------------
      pi=4.d0*datan(1.d0)
c     edges of 1BC
      b1x=2*pi/2
      b1y=-2*pi/(sqrt(3d0)*2)
      b2x=0
      b2y=4*pi/(sqrt(3d0)*2)
      b3x=2*pi/2
      b3y=2*pi/(sqrt(3d0)*2)
c     basis vectors for k
      bk1x=2*pi*ly2/n0
      bk1y=2*pi*(-2*lx2-ly2)/(sqrt(3d0)*n0)
      bk2x=-2*pi*ly1/n0
      bk2y=2*pi*(2*lx1+ly1)/(sqrt(3d0)*n0)
c     max integers for 1BC
      iimax=int(1+n0/(sqrt(3d0*(lx2**2+lx2*ly2+ly2**2))))
      jjmax=int(1+n0/(sqrt(3d0*(lx1**2+lx1*ly1+ly1**2))))
      iik=0   
      do ii=-iimax,iimax
         do jj=-jjmax,jjmax
            bckx=(ii)*bk1x+(jj)*bk2x -0.001
            bcky=(ii)*bk1y+(jj)*bk2y -0.0003
            xx1=(bckx*b1x+bcky*b1y)/(b1x*b1x+b1y*b1y)
            xx2=(bckx*b2x+bcky*b2y)/(b2x*b2x+b2y*b2y)
            xx3=(bckx*b3x+bcky*b3y)/(b3x*b3x+b3y*b3y)
c     if in 1BC
            if (-1d0.lt.xx1.and.xx1.le.1d0
     *          .and.-1d0.lt.xx2.and.xx2.le.1d0
     *          .and.-1.lt.xx3.and.xx3.le.1) then
c                write(0,*) 'wave vect in 1BC', xx1,xx2,xx3,ii,jj
                bckx1=ii*bk1x+jj*bk2x
                bcky1=ii*bk1y+jj*bk2y
                iik=iik+1
c                write(0,*)'bck1=',bckx1,bcky1,'    ',iik
                kall(iik,1)=bckx1
                kall(iik,2)=bcky1
           endif
         enddo
      enddo
      if (iik.ne.n0) then
         Write(0,*) 'ERROR:wrong num. of  allowed  wave vectors',n0,iik
      endif
c      write(0,*)' Wave vectors, all in 1BC'
c      write(0,2513) (i,kall(i,1),kall(i,2),i=1,iik)
c 2513 format(I4,2F20.10)
      return
      end
c=======================================================================
c   Sort array and calculate sign of the permutation  
c-----------------------------------------------------------------------
      integer function isort(n,arr)
      integer arr(n),a,iper
      iper=1
      do 12 j=2,n
        a=arr(j)
        do 11 i=j-1,1,-1
          if(arr(i).le.a) go to 10
          iper=-iper
          arr(i+1)=arr(i)
 11    continue
        i=0
10      arr(i+1)=a
 12   continue
      isort=iper
      return
      end
c=======================================================================
c   Find index of the configuration (up + down spins)
c-----------------------------------------------------------------------
      integer function ind(iuconf,idconf,nu,nd,ncu)
      implicit integer (i-n)
      integer iuconf(40),idconf(40)
      integer icu(40)
      common /iii/ nshft(40,40),ibin(40,40),n0,nni(40,6)
c-----------------------------------------------------------------------
      ju=0
      do i=1,n0
        if (iuconf(i).gt.0) then
          ju=ju+1
          icu(ju)=i
        endif
      enddo
      indu=ind0(icu,nu)
      jd=0
      do i=1,n0
        if (idconf(i).gt.0) then
          jd=jd+1
          icu(jd)=i
        endif
      enddo
      indd=ind0(icu,nd)
      ind=(indd-1)*ncu+indu
      return
      end
c=======================================================================
c   Find index of configuration (one kind of particles only)
c-----------------------------------------------------------------------
      integer function ind0(iconf,np)
      implicit integer (i-n)
      integer iconf(40)
      common /iii/ nshft(40,40),ibin(40,40),n0,nni(40,6)
c-----------------------------------------------------------------------
      indx=1
      do j=1,np
        indx=indx+ibin(j,iconf(j)-j+1)
      enddo
      ind0=indx
      return
      end
c=======================================================================
c   Construct configuration with given index (up + down spins)
c-----------------------------------------------------------------------
      subroutine conf(indx,iuconf,idconf,nu,nd,ncu)
      implicit integer (i-n)
      integer iuconf(40),idconf(40)
      integer icu(40),icd(40)
      common /iii/ nshft(40,40),ibin(40,40),n0,nni(40,6)
c-----------------------------------------------------------------------
      ii=indx-1
      i=ii/ncu
      indd=i+1
      indu=indx-i*ncu
      call conf0(indu,n0,nu,icu)
      call conf0(indd,n0,nd,icd)
      do i=1,n0
        iuconf(i)=0
        idconf(i)=0
      enddo
      do i=1,nu
        iuconf(icu(i))=1
      enddo
      do i=1,nd
        idconf(icd(i))=1
      enddo
      return
      end
c=======================================================================
c   Construct configuration with given index (one kind of particles)
c-----------------------------------------------------------------------
      subroutine conf0(indx,nn,np,iconf)
      implicit integer (i-n)
      integer iconf(40)
      common /iii/ nshft(40,40),ibin(40,40),n0,nni(40,6)
c-----------------------------------------------------------------------
      ind=indx-1
      j=nn-np+1
      do i=np,1,-1
80      if (ind.ge.ibin(i,j)) go to 90
        j=j-1
        go to 80
90      continue
        ind=ind-ibin(i,j)
        iconf(i)=j+i-1
      enddo
      return
      end
c=======================================================================
c   Shift configuration and find new index and relative fermion sign
c-----------------------------------------------------------------------
      integer function ishift(indx,is,iper,nu,nd,ncu)
      integer iuconf0(40),idconf0(40),iuconf(40),idconf(40),ip(40)
      common /iii/ nshft(40,40),ibin(40,40),n0,nni(40,6)
c-----------------------------------------------------------------------
      call conf(indx,iuconf0,idconf0,nu,nd,ncu)
      j=0
      do i=1,n0
        iuconf(nshft(is,i))=iuconf0(i)
        if (iuconf0(i).ne.0) then
          j=j+1
          ip(j)=nshft(is,i)
        endif
      enddo
      iper=isort(j,ip)
      j=0
      do i=1,n0
        idconf(nshft(is,i))=idconf0(i)
        if (idconf0(i).ne.0) then
          j=j+1
          ip(j)=nshft(is,i)
        endif
      enddo
      iper=iper*isort(j,ip)
      ishift=ind(iuconf,idconf,nu,nd,ncu)
      return
      end
c=======================================================================
c   Lanczos diagonalization
c-----------------------------------------------------------------------
      subroutine lanczos(hamdia,hamhopp,hamhopm,hamhopdir,
     *                   iham,nham0,nham,
     *                   phi,a,b,np,npa,lstmax,ldim,
     *                   eps1,eps2,da,ea,en,vlan,psix,cu,ctc,ctl)
      implicit real*8 (a-h, o-z)
      implicit integer (i-n)
      integer hamhopp(nham0),iham(npa+1)
      complex*16 hamhopm(nham0)
      integer*2 hamhopdir(nham0)
      complex*16 phi(0:npa,0:lstmax+1)
      real*8 a(0:lstmax-1),b(-2:lstmax-1)
      integer*2 hamdia(npa)
      real*8 da(lstmax),ea(lstmax)
      real*8 vlan(lstmax,lstmax)
      real*8 en(lstmax)
      complex*16 ctc(-40:40)
      real*8 ctl(40)
      
c.......................................................................
c     vector used only in subroutine lanstep, but must be translated
c     there form main program, because of space (16 site sistem)
      complex*16 psix(npa)
c-----------------------------------------------------------------------
c8      write(0,*)'new lanczos'
      do i=1,np
        phi(i,0)=dcmplx(0.d0,0.d0)
      enddo
      b(-2)=1.d0
      bc=0.d0
      do i=1,np
        bc=bc+dreal(phi(i,1))**2+dimag(phi(i,1))**2
      enddo
      b(-1)=dsqrt(bc)
      if (b(-1).lt.eps1) then
        ldim=0
        return
      endif
      j0=0
      lsteps=min(np,lstmax)
      do k=0,lsteps-1
c     for saving lanczos basis wavefunctions
         call lanstep(hamdia,hamhopp,hamhopm,hamhopdir,
     *              iham,nham0,nham,
     *              phi,a,b,j0,j0+1,j0+2,k,np,npa,lstmax,eps2,psix,
     *              cu,ctc,ctl)
        j0=1+j0   !saving lanczos basis wavefunc.
c   dimension of Hilbert subspace
        ldim=k+1
c.......................................................................
c calculates Lanczos eigenvectors and energies to check convergence
c stops Lanczos if energy of ground state has converged 
        do i=1,ldim
          da(i)=a(i-1)
          ea(i)=b(i-2)
        enddo
        ifail=0
        call tql1(ldim,da,ea,ifail)
        if (ifail.ne.0) stop 'ERROR in tql1'        
        if(da(1).ne.0) then
          if (dabs(en(1)/da(1)-1.d0).lt.eps2) then
             goto 9000
          endif
        endif
        en(1)=da(1)
c.......................................................................
        if (b(k).lt.eps1) then
c8           write(0,*) 'Hilbert space exhausted. ldim',ldim
           go to 9000
        endif 
      enddo
c      write(0,*) 'Required number of iterations performed. '
9000  continue
c.......................................................................
c   calculate Lanczos eigenvectors
        do i=1,ldim
          do j=1,ldim
            vlan(i,j)=0.d0
          enddo
        enddo
        do i=1,ldim
          vlan(i,i)=1.d0
          da(i)=a(i-1)
          ea(i)=b(i-2)
        enddo
        ifail=0
        call tql2(lstmax,ldim,da,ea,vlan,ifail)
c   the lowest energies
        do i=1,ldim
          en(i)=da(i)
        enddo
        return
      end
c=======================================================================
c   Lanczos Step
c-----------------------------------------------------------------------
       subroutine lanstep(hamdia,hamhopp,hamhopm,hamhopdir,
     *                   iham,nham0,nham,
     *                   phi,a,b,j0,j1,j2,k,np,npa,lstmax,eps,psi,
     *                   cu,ctc,ctl)
      implicit real*8 (a-h, o-z)
      implicit integer (i-n)
      integer hamhopp(nham0),iham(npa+1)
      complex*16 hamhopm(nham0)
      integer*2 hamhopdir(nham0)
      complex*16 phi(0:npa,0:lstmax+1),psi(npa)
      real*8 a(0:lstmax-1),b(-2:lstmax-1)
      integer*2 hamdia(npa)
      complex*16 ak,bt0
      complex*16 ctc(-40:40)
      real*8 ctl(40)
      complex*16 ctc0
c     omp
      integer id1
      INTEGER OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
      complex*16 ctemp,ak1
      
c-----------------------------------------------------------------------
c   act by hamiltonian on vector
      phi(0,j1)=dcmplx(0.d0,0.d0)
!$OMP PARALLEL DO
!$OMP&   DEFAULT(SHARED) PRIVATE(i)     
      do i=1,np
c        psi(i)=hamdia(i)*phi(i,j1)
        psi(i)=dcmplx(cu*hamdia(i),0.d0)*phi(i,j1)
      enddo
!$OMP END PARALLEL DO
      
!$OMP PARALLEL DO
!$OMP&  DEFAULT(SHARED) PRIVATE(i,ih,idir,ctc0,ctemp)     
      do i=1,np
c       ctemp - to write to shared phi(i) less times
        ctemp=dcmplx(0.d0,0.d0)
        do ih=iham(i),iham(i+1)-1
c         phase for new TBC (up and down different)
          idir=hamhopdir(ih)
          ctc0=ctl(abs(idir))*ctc(idir)
          ctemp=ctemp+phi(hamhopp(ih),j1)*ctc0*hamhopm(ih)
c          psi(i)=psi(i)+phi(hamhopp(ih),j1)*ctc0*hamhopm(ih)
        enddo
          psi(i)=psi(i)+ctemp
      enddo
!$OMP END PARALLEL DO
      
c   calculate new vector and new coefficients
      ak=dcmplx(0.d0,0.d0)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,ak1)
      ak1=dcmplx(0.d0,0.d0)
!$OMP DO 
      do i=1,np
        ak1=ak1+dconjg(phi(i,j1))*psi(i)
c!$OMP ATOMIC
c        ak=ak+dconjg(phi(i,j1))*psi(i)
      enddo
!$OMP END DO
!$OMP ATOMIC
      ak=ak+ak1
!$OMP END PARALLEL

c   checks reality of matrix element a(k)
      if ( abs(ak).gt.eps ) then
        diff=dimag(ak)/abs(ak)
        if (diff.gt.eps) write(0,*) 'imag. part too large: ',
     *        diff,abs(ak)
      endif
      a(k)=dreal(ak)/(b(k-1)**2)
c   report too large diffrence between  old b(k) and new one
      if( k.gt.0 ) then
        bt0=dcmplx(0.d0,0.d0)
        do i=1,np
          bt0=bt0+dconjg(phi(i,j0))*psi(i)
        enddo
        bt0=bt0/b(k-2)/b(k-1)
        diff=abs(bt0-dcmplx(b(k-1),0.d0))
        if( diff.gt.eps ) write(0,*) 'diff. in b_k too large: ',
     *        diff,b(k-1)
      endif
      c1=1/b(k-1)
      c2=-a(k)/b(k-1)
      c3=-b(k-1)/b(k-2)
!$OMP PARALLEL DO
!$OMP&  DEFAULT(SHARED) PRIVATE(i)           
      do i=1,np
        phi(i,j2)=c1*psi(i)+c2*phi(i,j1)+c3*phi(i,j0)
      enddo
!$OMP END PARALLEL DO
      bt=0.d0
!$OMP PARALLEL   DEFAULT(SHARED) PRIVATE(i,bt1)
      bt1=0.d0
!$OMP DO
      do i=1,np
        bt1=bt1+dreal(phi(i,j2))**2+dimag(phi(i,j2))**2
c!$OMP ATOMIC
c        bt=bt+dreal(phi(i,j2))**2+dimag(phi(i,j2))**2
      enddo
!$OMP END DO
!$OMP ATOMIC
      bt=bt+bt1
!$OMP END PARALLEL
      b(k)=dsqrt(bt)
      return
      end
c=======================================================================
c  operator - tau operator - like kinetic enery or bond order
c-----------------------------------------------------------------------
       subroutine opertau(hamhopp,hamhopm,hamhopdir,
     *                   iham,nham0,nham,
     *                   phi,j1,np,npa,lstmax,psi,
     *                   ctc,taudir)
      implicit real*8 (a-h, o-z)
      implicit integer (i-n)
      integer hamhopp(nham0),iham(npa+1)
      complex*16 hamhopm(nham0)
      integer*2 hamhopdir(nham0)
      complex*16 phi(0:npa,0:lstmax+1), psi(npa)
      complex*16 ctc(-40:40)
      real*8 taudir(-6:6)
      complex*16 ctc0
c     omp
      complex*16 ctemp
c-----------------------------------------------------------------------
!$OMP PARALLEL DO
!$OMP&   DEFAULT(SHARED) PRIVATE(i,ih,idir,ctc0,ctemp)     
      do i=1,np
        ctemp=dcmplx(0.d0,0.d0)
        do ih=iham(i),iham(i+1)-1
          idir=hamhopdir(ih)
          ctc0=ctc(idir)*taudir(idir)
          ctemp=ctemp+phi(hamhopp(ih),j1)*ctc0*hamhopm(ih)
c          psi(i)=psi(i)+phi(hamhopp(ih),j1)*ctc0*hamhopm(ih)
        enddo
        psi(i)=ctemp
      enddo
!$OMP END PARALLEL DO

      return
      end
c=======================================================================
c  operator - current opperator
c-----------------------------------------------------------------------
       subroutine opercur(hamhopp,hamhopm,hamhopdir,
     *                   iham,nham0,nham,
     *                   phi,j1,np,npa,lstmax,psi,
     *                   ctc,curdir)
      implicit real*8 (a-h, o-z)
      implicit integer (i-n)
      integer hamhopp(nham0),iham(npa+1)
      complex*16 hamhopm(nham0)
      integer*2 hamhopdir(nham0)
      complex*16 phi(0:npa,0:lstmax+1), psi(npa)
      complex*16 ctc(-40:40)
      real*8 curdir(-6:6)
      complex*16 ctc0
c     omp
      complex*16 ctemp
c-----------------------------------------------------------------------
!$OMP PARALLEL DO
!$OMP&   DEFAULT(SHARED) PRIVATE(i,ih,idir,ctc0,ctemp)     
      do i=1,np
        ctemp=dcmplx(0.d0,0.d0)
        do ih=iham(i),iham(i+1)-1
          idir=hamhopdir(ih)
c         setting prefactor for directon
             ctc0=ctc(idir)*curdir(idir)
             ctemp=ctemp+phi(hamhopp(ih),j1)*ctc0*hamhopm(ih)
c             psi(i)=psi(i)+phi(hamhopp(ih),j1)*ctc0*hamhopm(ih)
        enddo
        psi(i)=ctemp
      enddo
!$OMP END PARALLEL DO

      return
      end
c=======================================================================
c   PARENTS + HAMILTONIAN 
c-----------------------------------------------------------------------
      subroutine parham(nu,nd,nca,npa,np,ncu,irels,irelm,irelp,
     *                 ideg,ipar,nham0,nham,iham,hamhop,hamdia,
     *                 hamhopdir)           
c-----------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
      implicit integer (i-n)
      integer iham(npa+1),irelp(nca)
      integer hamhop(nham0)
      integer ihamm
      integer*2 hamdia(npa)
      integer*2 irels(nca),irelm(nca)
      integer iuconf(40),idconf(40),ipar(npa)
      integer*2 ideg(npa)
c     direciton of hopping for bond order and TBC
      integer*2 hamhopdir(nham0)     
      common /iii/ nshft(40,40),ibin(40,40),n0,nni(40,6)
c-----------------------------------------------------------------------
c   NUMBER OF CONFIGURATIONS:
c-----------------------------------------------------------------------
      ncu=1
      do i=1,nu
        ncu=ncu+ibin(i,n0-nu+1)
      enddo
      ncd=1
      do i=1,nd
        ncd=ncd+ibin(i,n0-nd+1)
      enddo
      nc=ncu*ncd
c      write(0,*) '   Total configurations : ',nc
      if (nc.gt.nca) stop 'Error in number of configurations. nc>nca'
c-----------------------------------------------------------------------
c   INDEXING:
c     particles:
c       hole      =  0
c       fermion   =  1
c     configuration in which up spin is on site i:   iuconf(i)=1
c     configuration in which down spin is on site i: idconf(i)=1
c     index i of the configuration iuconf + idconf:
c       i=ind(iuconf,idconf,nu,nd,ncu)
c     configuration iuconf + idconf corresponding to index i:
c       call conf(i,iuconf,idconf,nu,nd,ncu)
c-----------------------------------------------------------------------
c   PARENT CONFIGURATIONS:
c     there are np parent configurations with indices ipar(1:np);
c     translational degeneracy of j-th parent is ideg(j);
c     for each configuration i, the following is defined:
c        irelp(i)=related parent
c        irels(i)=shift with respect to related parent
c        irelm(i)=translational fermionic sign 
c-----------------------------------------------------------------------
      np=0
      ndeg=0
      do i=1,nc
        irelp(i)=0
      enddo
c   All configurations
      do indx0=1,nc
        if (irelp(indx0).eq.0) then
          np=np+1
          ipar(np)=indx0
          irelp(indx0)=np
          irels(indx0)=1
          irelm(indx0)=1
          id=0
c   Try all shifts
          do is=1,n0
            indx=ishift(indx0,is,iperm,nu,nd,ncu)
c   If shifted configuration equals itself, degeneracy
            if (indx.eq.indx0) then
              id=id+1
            else
              irelp(indx)=np
              irels(indx)=is
              irelm(indx)=iperm
            endif
          enddo
          ideg(np)=id
          if (id.gt.1) ndeg=ndeg+1
        endif
      enddo
      if (np.gt.npa) stop 'ERROR in number of parents.'
c-----------------------------------------------------------------------
c   HAMILTONIAN:
c     loop through parent configurations
c-----------------------------------------------------------------------
      nca0=nca+2
      nham=0
      do indp=1,np
        iham(indp)=nham+1
        call conf(ipar(indp),iuconf,idconf,nu,nd,ncu)
c-----------------------------------------------------------------------
c   HOPPING PART:
c     hamhop(ih) stores index of the configuration that is 
c     produced from the original parent state (indp) by l-th hop:
c        l=1,5 left; l=2,6 down; l=3,7 up; l=4,8 right
c     the sign of hamhop is equal to the fermion sign due to
c     rearangement of operators after hopping
c-----------------------------------------------------------------------
c8888   idir - direction of hopping  
c     first hopping of spin-up fermions
c   loop through sites occupied by up-spin
        do i=1,n0
          if( iuconf(i).eq.1 ) then
c   probe all neighbours, consider only those occupied by a fermion
      	    do idir=1,6  ! 2D 
c      	    do idir=1,2  ! 1D
              j=nni(i,idir)
              if (iuconf(j).eq.0) then
c   find change of fermion sign due to rearangement of operators
                isig=1
                do jj=min(i,j)+1,max(i,j)-1
                  if (iuconf(jj).ne.0) isig=-isig
                enddo
c   hop !
                iuconf(i)=0
                iuconf(j)=1
c   index of hop-related configuration 
                ij=ind(iuconf,idconf,nu,nd,ncu)
c   store related configuration and relative fermion sign
                nham=nham+1               
                ihamm=ij
                hamhop(nham)=ihamm*isig
c     store that hopping of up-spin (positive hamhopdir)
c     store direcition of hopping for bond order, positive for up spin
                hamhopdir(nham)=idir
c   restore parent configurations
                iuconf(j)=0
                iuconf(i)=1
              endif
            enddo ! idir
          endif
        enddo ! i
c   hopping of spin-down fermions
c     loop through sites occupied by down-spin
        do i=1,n0
          if( idconf(i).eq.1 ) then
c   probe all neighbours, consider only those unoccupied
            do idir=1,6  ! 2D  + nnn
c            do idir=1,2   ! 1D
              j=nni(i,idir)
              if (idconf(j).eq.0) then
c   find change of fermion sign due to rearrangement of operators
                isig=1
                do jj=min(i,j)+1,max(i,j)-1
                  if (idconf(jj).ne.0) isig=-isig
                enddo
c   hop !
                idconf(i)=0
                idconf(j)=1
c   index of hop-related configuration 
                ij=ind(iuconf,idconf,nu,nd,ncu)
c   store related configuration and relative fermion sign
                nham=nham+1
                ihamm=ij
                hamhop(nham)=ihamm*isig
c     store that hopping of down-spin negative hamhopdir
c     store direcition of hopping for bond order, negative for down spin
                hamhopdir(nham)=-idir
c   restore parent configuration
                idconf(j)=0
                idconf(i)=1
              endif
            enddo ! idir
          endif
        enddo !i
c        write(0,*)'indp',indp
c        write(0,1113)(hamhop(ih),ih=iham(indp),nham)
c        write(0,1113)(irelp(iabs(hamhop(ih))),ih=iham(indp),nham)
c888        write(0,*)' **** iham', iham
 1113   format(16i6)        
c-----------------------------------------------------------------------
c   DIAGONAL PART:
c     hamdia(indp) stores number of doubly occupied sites
c-----------------------------------------------------------------------
        hamdia(indp)=0
        do i=1,n0
            hamdia(indp)=hamdia(indp)+iuconf(i)*idconf(i)
        enddo
c-----------------------------------------------------------------------
      enddo ! indp
      iham(np+1)=nham+1
      write(0,3456)nham/(npa*6.*(nu+nd)+1.d0)
 3456 format('    rhh=nham/npa/nho =',f7.4)
      if (nham.gt.nham0) stop ' nham > nham0 '
      return
      end
c=======================================================================
c   EISPACK 
c    eigenvalues only of symetric tridiagonal matrix
c    http://www.netlib.org/cgi-bin/netlibget.pl/eispack/tql1.f
c-----------------------------------------------------------------------
      subroutine tql1(n,d,e,ierr)
      integer*4 i,j,l,m,n,ii,l1,l2,mml,ierr
      real*8 d(n),e(n)
      real*8 c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
      ierr = 0
      if (n .eq. 1) go to 1001
      do 100 i = 2, n
      e(i-1) = e(i)
 100  continue
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
      do 290 l = 1, n
         j = 0
         h = abs(d(l)) + abs(e(l))
         if (tst1 .lt. h) tst1 = h
         do 110 m = l, n
            tst2 = tst1 + abs(e(m))
            if (tst2 .eq. tst1) go to 120
  110    continue
  120    if (m .eq. l) go to 210
  130    if (j .eq. 30) go to 1000
         j = j + 1
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
         do 140 i = l2, n
         d(i) = d(i) - h
 140     continue
  145    f = f + h
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
  200    continue
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + abs(e(l))
         if (tst2 .gt. tst1) go to 130
  210    p = d(l) + f
         if (l .eq. 1) go to 250
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
  230    continue
  250    i = 1
  270    d(i) = p
  290 continue
      go to 1001
 1000 ierr = l
 1001 return
      end
c=======================================================================
c   EISPACK 
c    eigenvalues and eigenvectors of symetric tridiagonal matrix
c    http://www.netlib.org/cgi-bin/netlibget.pl/eispack/tql2.f
c    http://www.netlib.org/eispack/tql2.f
c-----------------------------------------------------------------------
      subroutine tql2(nm,n,d,e,z,ierr)
      integer*4 i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      real*8 d(n),e(n),z(nm,n)
      real*8 c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
      ierr = 0
      if (n .eq. 1) go to 1001
      do 100 i = 2, n
      e(i-1) = e(i)
  100 continue
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
      do 240 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
  110    continue
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
         do 140 i = l2, n
         d(i) = d(i) - h
  140    continue
  145    f = f + h
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
  200    continue
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
  300 continue
      go to 1001
 1000 ierr = l
 1001 return
      end
c-----------------------------------------------------------------------
      real*8 function pythag(a,b)
      real*8 a,b
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
c=======================================================================
