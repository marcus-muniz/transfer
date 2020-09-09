c
      parameter (nnei=400, nmaxnd=700000, nmaxvl=1400000, nmaxgh=65000,
     &           nmaxsd=1000000)
      implicit real*8 (a-h,o-z)
c 
      common /idata/  niter,nstart,nmiter,nnodes,nvol,nvisc,iflow,
     &                iflow2,idof,nitj,imove,itime,idissip,ivmesh,igcl,
     &                imotion
      common /save/   nsave,iresp,ipress,nhist,nintsave,npressave
      common /rdata/  fsmach,alpha,rey,gamma,h,pg,akc,ymax,alpham,xea,
     &                tmax
      common /connec/ itable(nmaxvl,3),neighbor(nmaxvl,3)
      common /conne2/ nsides,iedges(nmaxsd,4),iedaux(nmaxsd,2)
      common /nodes/  xynode(nmaxnd,2)
      common /cells/  q(nmaxvl,4), vol(nmaxvl,2), qaux(nmaxvl,4)
      common /fluxes/ qav(4), evec(4), fvec(4),
     &                uedge, vedge,
     &                uvolg(nmaxgh), vvolg(nmaxgh)
      common /ueflux/ eplus(nmaxvl,4), eminus(nmaxvl,4)
      common /ufflux/ fplus(nmaxvl,4), fminus(nmaxvl,4)
      common /ugeflx/ egplus(nmaxgh,4), egminus(nmaxgh,4)
      common /ugfflx/ fgplus(nmaxgh,4), fgminus(nmaxgh,4)
      common /operat/ cpen(nmaxvl,4), dpen(nmaxvl,4)
      common /igdata/ nvolg,ighost(nmaxgh,4)
      common /rgdata/ rghost(nmaxgh,6)
      common /ichdat/ ivtstep,ibcchr,irhsmax,ischeme
      common /rchdat/ epscon,epsblw,rhsmax,rdl2n
      common /rhside/ rhs(nmaxvl,4)
      common /rbcond/ pexit
      common /artdis/ ak2,ak4
      common /vtstep/ deltat(nmaxvl),clength(nmaxvl)
      common /idsmth/ irsmth,ljacrs,irstag
      common /conacc/ epsirs
      common /dynm/   ht,alphai,alphat,x0(nmaxnd),y0(nmaxnd),nref
      common /dynmh/  xhist(nmaxnd,2),yhist(nmaxnd,2),xytemp(nmaxnd,2)
      common /vel/    vside(nmaxsd,2)
      common /genf1/  icount,in1,in2,in3,inode,itneigh,it
      common /genf2/  neigh(nmaxnd,nnei),iwn(nmaxnd),iobn(nmaxnd),
     &                iin(nmaxnd),jtot(nmaxnd)
      common /genf3/  nwn,nobn,nin
      common /hist/   cltemp,clsdtemp,cdtemp,cmtemp,cmsdtemp,cmextemp
      common /time/   finaltime
c
