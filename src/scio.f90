!
!     SCIO  by Xi (Rossi) LUO <xi.rossi.luo@gmail.com>
!     created on May 20, 2011
!     version 0.1.4 on May 28, 2011
!     version 0.1.5 on June 13, 2011
!     version 0.2  on August 26, 2011
!     version 0.3  on March 23, 2012
!     version 0.4 
      subroutine solveColA(n,ic,S,b,rhocol,thr,maxit, niter)
!     Add active set strategy
      implicit double precision(a-h,o-z)
      parameter(eps=1e-6)
      double precision S(n,n), b(n), rhocol(n)
      double precision, dimension (:), allocatable :: ss
      double precision, dimension (:), allocatable :: negb
      logical, dimension(:), allocatable :: ja
      logical ia
      ia = .false.
      allocate(ss(1:n), stat=ierr)
      allocate(ja(1:n), stat=ierr)
      allocate(negb(1:n), stat=ierr)
!     Loop through
      ss(1:n) = 0.0
      negb(1:n) = -1.0*b(1:n)
      ss = matmul(S,negb)
!      ia = .false.
      ja(1:n) = .true.
      
      thr50=50.0*thr
      

      do 101 niter=1,maxit
!     Through Loop
         dlx=0.0
         do 102 j=1,n
            if (ia .and. ja(j)) goto 102
            bj=b(j)
            v = ss(j) + S(j,j)*bj
!     if (niter .lt. 5) call dblepr("v",1,v,1)
!     if (abs(v).lt.thr) v=0.0
            if (j .eq. ic) v = v + 1.0
            b(j) = 0.0
            vth = abs(v) - rhocol(j)
            if (vth .gt. 0.0) b(j) = sign(vth, v)/S(j,j) 
!     if (niter .lt. 5) call intpr("j",1,j,1)
            if ( abs(bj-b(j)).lt.eps ) goto 102
            del = b(j)-bj
            dlx = max(dlx, abs(del))
            ss = ss - del*S(:,j)
 102     continue       
         if (ia) then 
            if (dlx .lt. thr) ia = .false.
         else 
            if (dlx .lt. thr) goto 103
            if (dlx .lt. thr50) then
               ia = .true.
               ja(1:n) = abs(b(1:n)).lt.eps
            endif
         endif
 101  continue
 103  continue
      if (allocated(ss)) deallocate(ss)
      if (allocated(ja)) deallocate(ja)
      return
      end              

      
      subroutine scio(n, S, w, rhomat,  thr, maxit, nniter, jerr, isym)
      implicit double precision(a-h,o-z)
      double precision S(n,n), w(n,n), rhomat(n,n)
      double precision, dimension (:), allocatable :: b
      allocate(b(1:n), stat=jerr)
      nniter = 0
      niter = 0
      do 104 j=1,n
!     cold start
!         b(1:n)=0.0
!     warm start
         b(1:n)=w(1:n,j)
!         call intpr("col",3,j,1)
!     call solveCol(n,j, S, b, rho, thr, maxit, niter)
         call solveColA(n,j, S, b, rhomat(:,j), thr, maxit, niter)
         nniter = max(nniter, niter)
!         w(1:n,j) = b(1:n)
         w(:,j) = b
!         call dblepr("b",1,b,n)
 104  continue
      if (isym>0) then
!     symmetrize
         do 201 ji=1,(n-1)
            do 202 jj =(ji+1),n
               if (abs(w(ji,jj))>abs(w(jj,ji))) then
                  w(ji,jj) = w(jj,ji)
               else
                  w(jj,ji) = w(ji,jj)
               end if
 202     continue
 201  continue
      end if
      if (allocated(b)) deallocate(b)
      return
      end

      subroutine sciopath(wlist, n, S, w, rholist, nrho, thr, maxit, nniter, jerrlist, idiag, isym)
      implicit double precision(a-h,o-z)
      integer jerrlist(nrho)
      double precision S(n,n), wlist(n,n,nrho), rholist(nrho)
      double precision w(n,n)
      double precision, dimension (:,:), allocatable :: rhomat
      rho = rholist(nrho)
      nniter =1 
      niter = 1
      jerr = 0
      allocate(rhomat(1:n, 1:n), stat=ierr)
      jerr = jerr + ierr
      if (ierr.gt.0) goto 301
      rhomat = 0.0
      rhomat = rho
      if (idiag .eq. 0) then
         do ijj = 1, n
            rhomat(ijj, ijj) = 0.0
         end do 
      end if
      call scio(n, S, w, rhomat, thr, maxit, nniter, jerr, isym)
      call assignw(wlist, nrho, nrho, w, n)
      jerrlist(nrho) = jerr
 !     call dblepr("w",1,w,n*n)
      do 105 k=(nrho-1), 1, -1
         rho = rholist(k)
         rhomat = rho
         if (idiag .eq. 0) then
            do ijj = 1, n
               rhomat(ijj, ijj) = 0.0
            end do
         end if
         call scio(n, S, w, rhomat, thr, maxit, niter, jerr, isym)
         ! scio does not need input on whether diagonal should be penalized
         ! built in rhomat diagonal nonzero or not
         call assignw(wlist, k, nrho, w, n)
!         call dblepr("w",1,w,n*n)
         jerrlist(k) = jerr
         nniter = max(nniter, niter)
 105  continue
!     call dblepr("wlist",5,wlist(:,:,10),n*n)
 301  if (allocated(rhomat)) deallocate(rhomat)
      return
      end
      
      subroutine assignw(wlist, jrho, nrho, w, n)
      implicit double precision(a-h,o-z)
      double precision wlist(n,n,nrho), w(n,n)
      do 106 i=1,n
         do 107 j=1,n
            wlist(i,j,jrho) = w(i,j)
 107     continue
 106  continue
      return
      end
