c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c                      REQUIRED UTILITY FUNCTIONS
!   47:        thomas
!   141:       splc_pr
!   214:       thomas_pr
!   375:       fast1
!   421:       pentdag
!   484:       gaussj
!   643:       delta_fun
!   693:       Driv_delta_fun
!   742:       deltaCKINSIDE
!   786:       DrivdeltaCKINSIDE
!   854:       delta_width
!   895:       bracketA
!   933:       invert
!   995:       invertelement
!   1022:      bicgstab2 
!   1463:      ecmat
!   1488:      cjgr
!   1553:      crossvec
!   1568:      crossvec_noscale
!   1577:      shapen_dn_ddn
!   1768:      GMRES_CSR
!   1826:      dhels
!   1849:      dheqr
!   1884:      dorth
!   1957:      dscalm
!   1967:      daxpy
!   1977:      daxpy_x
!   1987:      Matvec
!   2000:      dMatzero
!   2012:      dVeczero
!   2023:      amux 
!   2100:      bicgstab
!   2330:      matrix_product 
!   2361:      ClosestPtPointEle2Ele
!   2453:      ClosestPtPointTriangle
!   2539:      dist3D_Segment_to_Segment
!   2636:      lubksb
!   2684:      ClosestPointHighOrder
!   2785:      S_Alpha
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine thomas 
     +
     +     (n
!     matrix size
     +     ,a
!     diagonal 
     +     ,b
!     super-diagonal row
     +     ,c
!     sub-diagonal row      
     +     ,s
!     rhs
     +     ,x
!     solution 
     +  ,nmax
     +  )
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c thomas algorithm for tridiagonal systems
c coefficient matrix:
c
c  | a1 b1  0   0  ...  0   0    0    |
c  | c2 a2  b2  0  ...  0   0    0    |
c  | 0  c3  a3  b3 ...  0   0    0    |
c  | ..............................   |
c  | 0  0   0   0  ... cn-1 an-1 bn-1 |
c  | 0  0   0   0  ...  0   cn   an   |
c
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      implicit double precision (a-h,o-z)
      integer nmax
      dimension a(nmax),b(nmax),c(nmax),s(nmax),x(nmax)
      dimension d(nmax),y(nmax)

      parameter (tol=0.00000001d0)

c--------
c prepare
c--------
c      print*,'n', n;pause
c      print*,'a', a;pause
c      print*,'b', b;pause
c      print*,'c', c;pause
c      print*,'s', s;pause
c      print*,'x', x;pause

      na = n-1

c------------------------------
c reduction to upper bidiagonal
c------------------------------

      d(1) = b(1)/a(1)
      y(1) = s(1)/a(1)

      do i=1,na
       i1 = i+1
       den   = a(i1)-c(i1)*d(i)
       d(i1) = b(i1)/den
       y(i1) = (s(i1)-c(i1)*y(i))/den
      end do

c------------------
c back substitution
c------------------

      x(n) = y(n)

      do i=na,1,-1
        x(i)= y(i)-d(i)*x(i+1)
      end do

c-----------------------
c verification and alarm
c-----------------------

      res = s(1)-a(1)*x(1)-b(1)*x(2)

      if(abs(res).gt.tol) write (6,*) " thomas: alarm"

      do i=2,na
        res = s(i)-c(i)*x(i-1)-a(i)*x(i)-b(i)*x(i+1)
        if(abs(res).gt.tol) write (6,*) " thomas: alarm"
      end do

      res = s(n)-c(n)*x(n-1)-a(n)*x(n)

      if(abs(res).gt.tol) write (6,*) " thomas: alarm"

c-----
c done
c-----

 100  format (1x,f15.10)



      return
      end
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine splc_pr(n,xs,fp,a,b,c,nmax)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c  cubic-spline with periodic-end conditions
c
c  n  .... number of intervals
c  xs .... s coord. of prescribed data
c  fp .... y or x coord. of prescribed data
c  a .... polynomial coefficient related to 3rd derivative
c  b .... polynomial coefficient related to 2nd derivative
c  c .... polynomial coefficient related to 1st derivative
c  h .... interval between prescribed data
c  amat.. diagonal of tridiagonal matrix
c  bmat.. superdiagonal of tridiagonal matrix
c  cmat.. subdiagonal of tridiagonal matrix
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	implicit none
	integer nmax,n,na,n1
	real  xs(nmax),fp(nmax),slope1,slope2,a(nmax),b(nmax),c(nmax)
	real  h(nmax)
	real  amat(n),bmat(n),cmat(n)
	real  rhs(n),sln(n)
	integer i,i1,i2	
c-----------------------
c compute intervals h(i)
c-----------------------
      do 10, i=1,n
        h(i) = xs(i+1)-xs(i)
10	continue 
c--------------------------------
c tridiagonal matrixlike as:
c
c | amat1 bmat1   0   	 0  . . .  0          cmat1    |
c | cmat2 amat2   bmat2  0  . . .  0          0        |
c | 0     cmat3   amat3  bmat3  .  0          0        |
c | 0     0       ... 	           0          0        |
c | 0     0       ...              bmat(n-2)  0        |
c | 0     0     	...    cmat(n-1) amat(n-1)  bmat(n-1)|
c | bmatn 0     ...  	 0         cmatn      amatn    |
c--------------------------------
      do 11, i=1,n
        amat(i) = 2.0d0*(h(i)+h(i+1))
        bmat(i) = h(i+1)
        cmat(i) = h(i)
11	continue
c-----
c generate the right-hand side
c-----
      do 12, i=1,n-1
        rhs(i) = 3.0d0*((fp(i+2)-fp(i+1))/h(i+1)-(fp(i+1)-fp(i) )/h(i))
12	continue
      rhs(n) = 3.0d0*( (fp(2) -fp(1))/h(1)-(fp(n+1)-fp(n))/h(n) )
c-----
c solve the tridiagonal system for b_i
c-----
      call thomas_pr (n,amat,bmat,cmat,rhs,sln,nmax)
c-------------------------
c extract the coefficients
c-------------------------
      do 13, i=1,n
        b(i+1) = sln(i)
13	continue
      b(1) = b(n+1)
c---------------------------------
c compute the coefficients a and c
c---------------------------------
      do 15, i=1,n
        a(i) = (b(i+1)-b(i))/(3.0d0*h(i))
        c(i) = (fp(i+1)-fp(i))/h(i) - h(i)*(b(i+1)+2.0d0*b(i))/3.0d0
15	continue
      return
      end
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine thomas_pr (n,a,b,c,s,x,nmax)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c thomas algorithm for tridiagonal systems
c
c  t . x = s
c with a periodic condition
c
c coefficient matrix:
c  | a1 b1  0   0  ...  0      0       c1   |
c  | c2 a2  b2  0  ...  0      0       0    |
c  | 0  c3  a3  b3 ...  0      0       0    |
c  | ....................................   |
c  | 0  0   0   0  ... c(n-1) a(n-1) b(n-1) |
c  | bn 0   0   0  ...  0      c(n)   a(n)  |
c
c------------------------------------------
	implicit none
	integer n,nmax,na,nb
	real  a(nmax),b(nmax),c(nmax),s(nmax),x(nmax),tol
	real  d(n),y(n),den,res,r1,r0,x0(nmax)
	real  save1,savena
	integer i,i1
      parameter (tol=0.00000001d0)
c--------
c parameter
c--------
      na = n-1
      nb = n-2
      save1  = s(1)
      savena = s(na)
c------------------------------
c first assume that x(n) = 0
c and solve the first n-1 equations
c neglecting the last column
c and the last row
c------------------------------
      x(n) = 0.0d0
c-- regular thomas

c-----
c reduction to upper bidiagonal
c-----
      d(1) = b(1)/a(1)
      y(1) = s(1)/a(1)
      do 11,i=1,nb
       i1 = i+1
       den   = a(i1)-c(i1)*d(i)
       d(i1) = b(i1)/den
       y(i1) = (s(i1)-c(i1)*y(i))/den
11	continue
c----
c back substitution
c----
      x(na) = y(na)
      do 12, i=nb,1,-1
        x(i)= y(i)-d(i)*x(i+1)
12	continue
c-----
c compute the first residual:
c-----
      r0 = a(n)*x(n) + b(n)*x(1) + c(n)*x(na) - s(n)
c-----
c save the solution
c-----
      do 13, i=1,n
        x0(i) = x(i)
13	continue

c------------------------------
c second, assume that x(n) = 1
c and solve the first n-1 equations
c with a modified rhs
c------------------------------

      x(n) = 1.0d0

      s(1)  = s(1)  - c(1)  * x(n)
      s(na) = s(na) - b(na) * x(n)


c-- regular thomas

c-----
c reduction to upper bidiagonal
c-----


      d(1) = b(1)/a(1)
      y(1) = s(1)/a(1)

      do i=1,nb
       i1 = i+1
       den   = a(i1)-c(i1)*d(i)
       d(i1) = b(i1)/den
       y(i1) = (s(i1)-c(i1)*y(i))/den
      end do

c-----
c back substitution
c-----

      x(na) = y(na)

      do i=nb,1,-1
        x(i)= y(i)-d(i)*x(i+1)
      end do

c------
c compute the second residual:
c-----

      r1 = a(n)*x(n) + b(n)*x(1) + c(n)*x(na) - s(n)

c----------------------------
c rectify the right-hand side
c----------------------------

      s(1)  = save1
      s(na) = savena

c---------------------------------
c compute the correct value of x(n)
c---------------------------------

      x(n) = -r0/(r1-r0)

c----------------------------
c compose the solution vector
c----------------------------

      do i=1,na
       x(i) = (x(i)-x0(i)) * x(n) + x0(i)
      end do

c-----------------------
c verification and alarm
c-----------------------

      res = s(1) - a(1)*x(1) -b(1)*x(2)-c(1)*x(n)

      if(abs(res).gt.tol) write (6,*) " thomas_pr: alarm, 1",res

      do i=2,na
        res = s(i)-c(i)*x(i-1)-a(i)*x(i)-b(i)*x(i+1)
        if(abs(res).gt.tol) write (6,*) " thomas_pr: alarm ",i,res
      end do

      res = s(n)-c(n)*x(na)-a(n)*x(n)-b(n)*x(1)

      if(abs(res).gt.tol) write (6,*) " thomas_pr: alarm ",n,res

c-----
c done
c-----

 100  format (1x,f15.10)

      return
      end
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine fast1(x,ndim,ione)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c
c implement the fast fourier transform
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
 
      implicit real(a-h,o-z)
      complex  t,u,w,x(ndim),ci
      data pi1,ci/3.1415926535897932384d0,(0.d0,1.d0)/
      npow=int(log(ndim+1.d0)/log(2.d0))
      n1=2**npow
      n2=n1/2
      n3=n1-1
      j= 1
      do 300 i=1,n3
      if(i.ge.j) go to 100

      t=   x(j)
      x(j)=x(i)
      x(i)=t
100   k=n2
200   if(k.ge.j) go to 300
      j=j-k
      k=k/2
      go to 200
300   j=j+k
      do 500 l=1,npow
      le1=2**l
      le2=le1/2
      ang=pi1/le2
      u=1.d0
      w=cos(ang)+ci*ione*sin(ang)
      do 500 j=1,le2
      do 400 i=j,n1,le1
      ip=i+le2
      t=    x(ip)*u
      x(ip)=x(i)-t
400   x(i)= x(i)+t
500   u=u*w
      if(ione.eq.1)return
      scl= 1.d0/real(n1)
      do 600 i=1,n1
600   x(i)=x(i)*scl
      return
      end
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine pentdag(a,b,c,d,e,f,u,n) 
      implicit none
      save

c..solves for a vector u of length n in the pentadiagonal linear system 
c.. a_i u_(i-2) + b_i u_(i-1) + c_i u_i + d_i u_(i+1) + e_i u_(i+2) = f_i 
c..input are the a, b, c, d, e, and f and they are not modified 

c..in its clearest incarnation, this algorithm uses three storage arrays 
c..called p, q and r. here, the solution vector u is used for r, cutting 
c..the extra storage down to two arrays. 

c..declare the pass
      integer          n
      double precision a(n),b(n),c(n),d(n),e(n),f(n),u(n)

c..local variables
      integer          nmax,i 

      integer        i2

      parameter        (nmax=500) 
      double precision p(nmax),q(nmax),bet,den 


244    format(24(f20.9,'  '))


c..initialize elimination and backsubstitution arrays 
      if (c(1) .eq. 0.0)  stop 'eliminate u2 trivially' 
      bet  = 1.0d0/c(1) 
      p(1) = -d(1) * bet 
      q(1) = -e(1) * bet 
      u(1) = f(1)  * bet 

      bet = c(2) + b(2)*p(1) 
      if (bet .eq. 0.0) stop 'singular 1 in pentdag' 
      bet = -1.0d0/bet 
      p(2) = (d(2) + b(2)*q(1)) * bet 
      q(2) = e(2) * bet 
      u(2) = (b(2)*u(1) - f(2)) * bet 


c..reduce to upper triangular 
      do i=3,n 
       bet = b(i) + a(i) * p(i-2) 
       den = c(i) + a(i)*q(i-2) + bet*p(i-1) 
       if (den .eq. 0.0) stop 'singular 2 in pentdag' 
       den = -1.0d0/den 
       p(i) = (d(i) + bet*q(i-1)) * den 
       q(i) = e(i) * den 
       u(i) = (a(i)*u(i-2) + bet*u(i-1) - f(i)) * den 
      enddo

c..backsubstitution 
      u(n-1) = u(n-1) + p(n-1) * u(n) 
      do i=n-2,1,-1 
       u(i) = u(i) + p(i) * u(i+1) + q(i) * u(i+2) 
      enddo
      return 

      end 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine gaussj(a,n,np,b,m,mp,ierr)

c  purpose: solution of the system of linear equations ax = b by
c     gauss-jordan elimination, where a is a matrix of order n and b is
c     an n x m matrix.  on output a is replaced by its matrix inverse
c     and b is preplaced by the corresponding set of solution vectors.

c  source: w.h. press et al, "numerical recipes," 1989, p. 28.

c  modifications: 
c     1. double  precision.
c     2. error parameter ierr included.  0 = no error. 1 = singular 
c        matrix encountered; no inverse is returned.

c  prepared by j. applequist, 8/17/91.

      implicit real*8(a-h,o-z)

c        set largest anticipated value of n.

      parameter (nmax=500)
      dimension a(-3:np+3,-3:np+3),b(np,mp),
     &          ipiv(nmax),indxr(nmax),indxc(nmax)
      ierr=0
      do 11 j=1,n
      ipiv(j)=0
 11   continue
      do 22 i=1,n
      big=0.d0
      do 13 j=1,n
      if (ipiv(j).ne.1) then
      do 12 k=1,n
      if (ipiv(k).eq.0) then
      if (dabs(a(j,k)).ge.big) then
      big=dabs(a(j,k))
      irow=j
      icol=k
      endif
      else if (ipiv(k).gt.1) then
      ierr=1
      return
      endif
 12   continue
      endif
 13   continue
      ipiv(icol)=ipiv(icol)+1
      if (irow.ne.icol) then
      do 14 l=1,n
      dum=a(irow,l)
      a(irow,l)=a(icol,l)
      a(icol,l)=dum
 14   continue
      do 15 l=1,m
      dum=b(irow,l)
      b(irow,l)=b(icol,l)
      b(icol,l)=dum
 15   continue
      endif
      indxr(i)=irow
      indxc(i)=icol
      if (a(icol,icol).eq.0.d0) then
      ierr=1
      return
      endif
      pivinv=1.d0/a(icol,icol)
      a(icol,icol)=1.d0
      do 16 l=1,n
      a(icol,l)=a(icol,l)*pivinv
 16   continue
      do 17 l=1,m
      b(icol,l)=b(icol,l)*pivinv
 17   continue
      do 21 ll=1,n
      if (ll.ne.icol) then
      dum=a(ll,icol)
      a(ll,icol)=0.d0
      do 18 l=1,n
      a(ll,l)=a(ll,l)-a(icol,l)*dum
 18   continue
      do 19 l=1,m
      b(ll,l)=b(ll,l)-b(icol,l)*dum
 19   continue
      endif
 21   continue
 22   continue
      do 24 l=n,1,-1
      if (indxr(l).ne.indxc(l)) then
      do 23 k=1,n
      dum=a(k,indxr(l))
      a(k,indxr(l))=a(k,indxc(l))
      a(k,indxc(l))=dum
 23   continue
      endif
 24   continue
      return
      end
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      function s_function(x)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c function for fringe region s(x)
c inputs: x value 
c outputs: value of function
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      real x
      if (x .le. 0.0) then 
         s_function=0.0
      else if (x .lt. 1.0) then
         s_function=1.0/(1.0+exp(1.0/(x-1.0)+1.0/(x)))
      else 
         s_function=1.0
      end if
      return
      end
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      function ds_function(x)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c function for fringe region s(x)
c inputs: x value 
c outputs: value of function
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      real x
      if (x .le. 1.0d-2) then 
         ds_function=0.0
      else if (x .lt. 1.0-1.0d-2) then
         ds_function=
     &     (exp(1.0/(x - 1.0) + 1.0/x)
     &    *(1.0/(x - 1.0)**2 + 1.0/x**2))
     &    /(exp(1.0/(x - 1.0) + 1.0/x) + 1)**2
      else 
         ds_function=0.0
      end if
      return
      end
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      function areafun(x,y,n)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c inputs: array x with the x coords, and y with the y coords. 
c n: number of verts.
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      real x(*),y(*),area
      area=0.0 
      if (n-3) 1,1,2 
 1    write(*,3) n 
 3    format(19hunusual value of n=,i10,17h in call to area.) 
      return 
 2    if(n-1000000) 7,1,1 
 7    do 4 i=3,n,2 
 4    area=area+(x(i-2)-x(i))*y(i-1)-(y(i-2)-y(i))*x(i-1) 
      if (mod(n,2)) 6,5,6
 5    area=area+x(n-1)*y(n)-y(n-1)*x(n) 
 6    area=area*0.5 
      areafun=area
      return
      end
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine delta_fun(delta_type,r,f)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      implicit none
      integer delta_type 
      real r,f,pi1,fi_ib_6,K,r2,r3,rp 
      pi1 = 4. * atan(1.0)
      f=0.0
      select case (delta_type)
        case(1)
            r=abs(r)
           if (r .lt. 2.d0) f=0.25d0*(1.0d0+dcos(pi1*r/2.d0))
        case(2)
            r=abs(r)
           if (r .lt. 1.d0) then
              f=0.125*(3.-2.*r+sqrt(1.+4.*r-4.*r**2.))
           else if (r .lt. 2.d0) then
              f=0.125*(5.-2.*r-sqrt(-7.+12.*r-4.*r**2.))
           end if
        case(3)
            r=abs(r)
           if (r .le. 1.d0) then
              f=fi_ib_6(r)
           else if (r .le. 2.d0) then
              f=21./16.+7./12.*r-7./8.*r**2+1./6.*r**3-1.5*fi_ib_6(r-1.)
           else if (r .le. 3.d0) then
              f=9./8.-23./12.*r+3./4.*r**2-1/12.*r**3+0.5*fi_ib_6(r-2.)
           end if
        case(4)
            r=abs(r)
           if (r .lt. 1.d0) then
              f=1.-0.5*r-r**2+0.5*r**3
           else if (r .lt. 2.d0) then
              f=1.-11./6.*r+r**2-1./6.*r**3
           end if
        case(6)
           K=0.714075092976608   !59.0d0/60.0d0-sqrt(29.0d0)/20.0d0
           if((-3.0.ge. r) .and. (r .le. 3)) then
              Rp = r - ceiling(r) + 1.0  !!Rp between [0,1] 
              R2 = Rp * Rp; R3 = R2*Rp
              call deltaCKINSIDE(f,rp,r2,r3,K,r)
            end if  
        case default
           write(*,*) 'wrong delta type'
          stop
        end select

      return
      end 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine Driv_delta_fun(delta_type,r,f)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      implicit none
      integer delta_type 
      real r,f,pi1,dfi_ib_6 ,r2,r3,rp, K
      pi1 = 4. * atan(1.0)
      f=0.0
      select case (delta_type)
        case(1)
            r=abs(r)
           if (r .lt. 2.d0) f=-(pi1*sin((pi1*r)/2.0))/8.0
        case(2)
            r=abs(r)
           if (r .lt. 1.d0) then
              f=- (8.0*r-4.0)/(16.0*sqrt(-4.0*r**2+4.0*r+1.0)) - 0.25
           else if (r .lt. 2.d0) then
              f=  (8.0*r-12.0)/(16*sqrt(-4.0*r**2+12.0*r-7.0)) - 0.25
           end if
        case(3)
            r=abs(r)
           if (r .le. 1.d0) then
              f=dfi_ib_6(r)
           else if (r .le. 2.d0) then
              f=0.5*r**2 - 1.75*r + 0.58333333-1.5*dfi_ib_6(r-1.)
           else if (r .le. 3.d0) then
              f=- 0.25*r**2 + 1.5*r - 1.9166667+0.5*dfi_ib_6(r-2.)   
           end if
        case(4)
            r=abs(r)
           if (r .lt. 1.d0) then
              f=1.5*r**2 - 2.0*r - 0.5
           else if (r .lt. 2.d0) then
              f=- 0.5*r**2 + 2.0*r - 1.8333333
           end if
        case(6)
           K=0.714075092976608   !59.0d0/60.0d0-sqrt(29.0d0)/20.0d0
           if((-3.0.ge. r) .and. (r .le. 3)) then
              Rp = r - ceiling(r) + 1.0  !!Rp between [0,1] 
              R2 = Rp * Rp; R3 = R2*Rp
              call DrivdeltaCKINSIDE(f,rp,r2,r3,K,r)
            end if  
        case default
           write(*,*) 'wrong delta type'
          stop
        end select
      return
      end 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine deltaCKINSIDE(phi,r,R2,R3,K,rc)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      real r,K,r2,r3,rc
      real alpha,beta,gamma,discr,sgnK,phi

      alpha=28.0

      beta = 9./4. - 1.5 * (K + R2) + (22./3.-7.*K)*R - 7./3.*R3

      gamma = 0.25 * ( 0.5*(161./36. - 59./6.*K + 5.*K*K)*R2 
     &          + 1./3.*(-109./24. + 5.*K)*R2*R2 + 5./18.*R3*R3  )

      discr=alpha*gamma*(-4.0D0)+beta**2

      sgnK=-1.0
      if((1.5-K) .gt. 0) sgnK=1.0
          if ( (-3.0 < rc) .and. (rc <= -2.0) ) then
            phi = 1./(2*alpha) * ( -beta + sgnK * sqrt(discr) )
          else if ( (-2.0 < rc) .and. (rc <= -1.0) ) then
            phi = -3./(2*alpha) * ( -beta + sgnK * sqrt(discr) ) 
     &          - 1./16 + 1./8*( K+(rc+2)*(rc+2) )  
     &          + 1./12*(3*K-1)*(rc+2) + 1./12*(rc+2)*(rc+2)*(rc+2)
          else if ( (-1.0 < rc) .and. (rc <= 0.0) ) then
            phi = 2./(2*alpha) * ( -beta + sgnK * sqrt(discr) )  
     &          + 1./4 + 1./6*(4-3*K)*(rc+1)  
     &          - 1./6*(rc+1)*(rc+1)*(rc+1)
          else if ( (0.0 < rc) .and. (rc <= 1.0) ) then
            phi = 2./(2*alpha) * ( -beta + sgnK * sqrt(discr) )  
     &          +5./8 - 1./4 * ( K+rc*rc )
          else if ( (1.0 < rc) .and. (rc <= 2.0) ) then
            phi = -3./(2*alpha) * ( -beta + sgnK * sqrt(discr) )  
     &          +  1./4 - 1./6*(4-3*K)*(rc-1)  
     &          + 1./6*(rc-1)*(rc-1)*(rc-1)
          else if ( (2.0 < rc) .and. (rc <= 3.0) ) then
            phi = 1./(2*alpha) * ( -beta + sgnK * sqrt(discr) )  
     &          -   1./16 + 1./8*(K+(rc-2)*(rc-2))  
     &          - 1./12*(3*K-1)*(rc-2) - 1./12*(rc-2)*(rc-2)*(rc-2)
          else
            phi=0.0
          endif     
        return
        end subroutine
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine DrivdeltaCKINSIDE(phi,r,R2,R3,K,rc)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      real r,K,r2,r3,rc
      real alpha,beta,gamma,discr,sgnK,phi,dbeta,dgamma,pm3,dpm3

      alpha=28.0

      beta = 9./4. - 1.5 * (K + R2) + (22./3.-7.*K)*R - 7./3.*R3
      dbeta =                 -3.*R + (22./3-7*K)   - 7*R2  

      gamma = 0.25 * ( 0.5*(161./36. - 59./6.*K + 5.*K*K)*R2 
     &          + 1./3.*(-109./24. + 5.*K)*R2*R2 + 5./18.*R3*R3  )
      dgamma= 0.25 * (     (161./36 - 59./6*K + 5*K*K)*R  
     &          + 4./3*(-109./24 + 5*K)*R3    + 5./3 *R3*R2 )
      discr=alpha*gamma*(-4.0D0)+beta**2


      sgnK=-1.0
      if((1.5-K) .gt. 0) sgnK=1.0

      pm3  =  (-beta + sgn*sqrt(discr)) / (2*alpha)
      dpm3 = -(dbeta*pm3+dgamma) / (2*alpha*pm3+beta)

          if ( (-3.0 < rc) .and. (rc <= -2.0) ) then
            phi = dpm3
          else if ( (-2.0 < rc) .and. (rc <= -1.0) ) then
            phi = -3*dpm3 + 1./12*(3*K-1) + 1./4*R + 1./4*R2
          else if ( (-1.0 < rc) .and. (rc <= 0.0) ) then
            phi = 2*dpm3 +  1./6*(4-3*K)          - 1./2*R2
          else if ( (0.0 < rc) .and. (rc <= 1.0) ) then
            phi = 2*dpm3                 - 1./2*R
          else if ( (1.0 < rc) .and. (rc <= 2.0) ) then
            phi = -3*dpm3 - 1./6 *(4-3*K)          + 1./2*R2
          else if ( (2.0 < rc) .and. (rc <= 3.0) ) then
            phi =dpm3 - 1./12*(3*K-1) + 1./4*R - 1./4*R2
          else
            phi=0.0
          endif     
        return
        end subroutine
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      real function fi_ib_6(r)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      implicit none
      real r
      fi_ib_6=61./112.-11./42.*r-11./56.*r**2.+1./12.*r**3
     &         +sqrt(3.)/336.*sqrt(243.+1584.*r-748.*r**2-1560.*r**3
     &         +500.*r**4+336.*r**5-112.*r**6)
      return
      end
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      real function dfi_ib_6(r)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      implicit none
      real r
      dfi_ib_6=0.25*r**2
     &        -(0.0025774566
     &          *(672.0*r**5 - 1680.0*r**4 - 2000.0*r**3 
     &            + 4680.0*r**2 + 1496.0*r - 1584.0))
     &          /sqrt(- 112.0*r**6 + 336.0*r**5 + 500.0*r**4
     &                - 1560.0*r**3 - 748.0*r**2 + 1584.0*r + 243.0)
     &        - 0.39285714*r - 0.26190476
      return
      end
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine delta_width(del_x,del_y,del_z
     &                ,delta_typex,delta_typey,delta_typez)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c width of delta fuction ( in grid) in each direction
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      implicit none
      integer del_x,del_y,del_z,delta_typex,delta_typey,delta_typez
      select case (delta_typex)
        case(1)
           del_x=3
        case(2)
           del_x=3
        case(3)
           del_x=4
        case(4)
           del_x=3
      end select
      select case (delta_typey)
        case(1)
           del_y=3
        case(2)
           del_y=3
        case(3)
           del_y=4
        case(4)
           del_y=3
      end select
      select case (delta_typez)
        case(1)
           del_z=3
        case(2)
           del_z=3
        case(3)
           del_z=4
        case(4)
           del_z=3
      end select
      return
      end 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine bracketA(A,x) 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
         real x(3) ! input
         real A(3,3) ! input
         A=0.0
         A(1,2)=-x(3)
         A(1,3)= x(2)
         A(2,1)= x(3)
         A(2,3)=-x(1)
         A(3,1)=-x(2)
         A(3,2)= x(1)
      end subroutine
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	real function dist(xa,xe)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c    computes the distance between two points
c    with coordinates (xa,ya) and (xe,ye)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	implicit none
c cartesian dimension
	integer cdim 
	parameter(cdim=3)      
c coords of point 1
	real	 xa(cdim)
c coords of point 2      
	real	 xe(cdim)      
	integer   :: n
	real	:: sums
	sums= 0.0d0
	do 10,  n=1,cdim
         sums= sums + (xa(n)-xe(n))**2.d0
10	continue
	dist= sqrt(sums)
	return
	end 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine invert(matrix,ndim)
!
! this subroutine inverts a small square matrix onto itself.
!
       implicit none
       integer ndim,i,k
       real matrix(1:ndim,1:ndim)
       real det,j11,j12,j13,j21,j22,j23,j31,j32,j33,con
       if(ndim .eq. 2)then
         det=matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1)
         j11=matrix(1,1)
         matrix(1,1)=matrix(2,2)
         matrix(2,2)=j11
         matrix(1,2)=-matrix(1,2)
         matrix(2,1)=-matrix(2,1)
         matrix=matrix/det
       else if(ndim .eq. 3)then
         det=matrix(1,1)*(matrix(2,2)*matrix(3,3)
     &        -matrix(3,2)*matrix(2,3))
         det=det-matrix(1,2)
     &             *(matrix(2,1)*matrix(3,3)
     &        -matrix(3,1)*matrix(2,3))
         det=det+matrix(1,3)
     &             *(matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2))
         j11=matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3)
         j21=-matrix(2,1)*matrix(3,3)+matrix(3,1)*matrix(2,3)
         j31=matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2)
         j12=-matrix(1,2)*matrix(3,3)+matrix(3,2)*matrix(1,3)
         j22=matrix(1,1)*matrix(3,3)-matrix(3,1)*matrix(1,3)
         j32=-matrix(1,1)*matrix(3,2)+matrix(3,1)*matrix(1,2)
         j13=matrix(1,2)*matrix(2,3)-matrix(2,2)*matrix(1,3)
         j23=-matrix(1,1)*matrix(2,3)+matrix(2,1)*matrix(1,3)
         j33=matrix(1,1)*matrix(2,2)-matrix(2,1)*matrix(1,2)
         matrix(1,1)=j11
         matrix(1,2)=j12
         matrix(1,3)=j13
         matrix(2,1)=j21
         matrix(2,2)=j22
         matrix(2,3)=j23
         matrix(3,1)=j31
         matrix(3,2)=j32
         matrix(3,3)=j33
         matrix=matrix/det      
       else
         do k=1,ndim
           con=matrix(k,k)
           matrix(k,k)=1.0
           matrix(k,:)=matrix(k,:)/con
           do i=1,ndim
             if(i .ne. k)then
               con=matrix(i,k)
               matrix(i,k)=0.0
               matrix(i,:)=matrix(i,:)-matrix(k,:)*con
             end if

           end do
         end do
       end if
      return
      end 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine invertelement(matrix,ndim)
!
! this subroutine inverts a small square matrix onto itself.
!
       implicit none
       integer ndim,i,k
       real matrix(1:ndim,1:ndim)
       real det,j11,j12,j13,j21,j22,j23,j31,j32,j33,con
       do k=1,ndim
          do i=1,ndim
            if (matrix(i,k) .le. 1.0d-12) then
                matrix(i,k)=1.0d-12
            endif
            matrix(i,k)=1.0/matrix(i,k)
          enddo  
       enddo
      return
      end 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c
c      call bicgstab2 (.false.,2, ns_ibm_r(i),
c     &       gx_ibm_massive(i,1: ns_ibm_r(i)), fibm1(i,1: ns_ibm_r(i)),
c     &        .false., tol,
c     &      'rel',itmax, info,aamat(1: ns_ibm_r(i),1: ns_ibm_r(i)))
c
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine bicgstab2 (okprint,l, n, x, rhs,  nonzero, tol,
     &      typestop,mxmv, info,aamat)

      implicit none
c      external matvec
      integer l, n, mxmv,  info
      real  x(n), rhs(n), tol,aamat(n,n)
      logical   okprint,nonzero
      character*3   typestop
      real  work(n,5+2*l)
c     -----------------------------------------
      integer   lmax
      parameter(lmax=2)
      real rwork(lmax+1,3+2*(lmax+1))

      logical goon, rcmp, xpdt
      integer ii, i1, jj, kk, nmv
      real alpha,beta,hatgamma,kappa0, kappal,maxval1,
     &     mxnrmr,mxnrmx,omega,rho0,rho1,rnrm0,rnrm,rnrmmax,
     &     sigma,sum1,varrho
      integer z, zz, y0, yl, y

      integer rr, r, u, xp, bp

      real    zero, one, delta
      parameter(zero=0e0,one=1e0,delta=1.0e-2)


c 
c simple bicgstab(\ell) iterative method, \ell <= 2
c by m.a.botchev, jan.'98 
c report bugs to botchev@cwi.nl or botchev@excite.com
c
c copyright (c) 1998 by m.a.botchev
c permission to copy all or part of this work is granted,
c provided that the copies are not made or distributed
c for resale, and that the copyright notice and this
c notice are retained.
c
c this is the "vanilla" version of bicgstab(\ell) as described
c in phd thesis of d.r.fokkema, chapter 3.  it includes two enhancements 
c to bicgstab(\ell) proposed by g.sleijpen and h.van der vorst in
c 1) g.sleijpen and h.van der vorst "maintaining convergence 
c    properties of bicgstab methods in finite precision arithmetic",
c    numerical algorithms, 10, 1995, pp.203-223
c 2) g.sleijpen and h.van der vorst "reliable updated residuals in
c    hybrid bicg methods", computing, 56, 1996, pp.141-163
c
c {{ this code based on:
c subroutine bistbl v1.0 1995
c
c copyright (c) 1995 by d.r. fokkema.
c permission to copy all or part of this work is granted,
c provided that the copies are not made or distributed
c for resale, and that the copyright notice and this
c notice are retained.  }}
c
c okprint == (input) logical. if okprint=.true. residual norm
c            will be printed to *
c l       == (input) integer bicgstab's dimension <= 2
c            set l=2 for highly nonsymmetric problems
c n       == (input) integer size of the system to solve 
c x       == (input/output) double precision array dimension n
c            initial guess on input, solution on output
c rhs     == (input) double precision array dimension n
c            right-hand side (rhs) vector
c matvec  == (input) external name of matrix vector subroutine
c            to deliver y:=a*x by call matvec(n,x,y)
c nonzero == (input) logical tells
c            bicgstab(\ell) if initial guess in x is zero or not. 
c            if nonzero is .false., one matvec call is saved.
c tol     == (input/output) double precision tolerance for all possible
c            stopping criteria (see the 'typestop' parameter)
c            on output, if info=0 or 1, tol is actually achieved
c            residual reduction or whatever (see the 'typestop' parameter)
c typestop== (input) character*3 stopping criterion (||.|| denotes 
c            the 2-norm):
c            typestop='rel' -- relative stopping crit.: ||res|| < tol*||res0||
c            typestop='abs' -- absolute stopping crit.: ||res||<tol
c            typestop='max' -- maximum  stopping crit.: max(abs(res))<tol
c note(for typestop='rel' and 'abs'): to save computational work, the value of 
c            residual norm used to check the convergence inside the main iterative 
c            loop is computed from 
c            projections, i.e. it can be smaller than the true residual norm
c            (it may happen when e.g. the 'matrix-free' approach is used).
c            thus, it is possible that the true residual does not satisfy
c            the stopping criterion ('rel' or 'abs').
c            the true residual norm (or residual reduction) is reported on 
c            output in parameter tol -- this can be changed to save 1 matvec
c            (see comments at the end of the subroutine)
c mxmv   ==  (input/output) integer.  on input: maximum number of matrix 
c            vector multiplications allowed to be done.  on output: 
c            actual number of matrix vector multiplications done
c work   ==  (workspace) double precision array dimension (n,2*l+5))
c ldw    ==  (input) integer size of work, i.e. ldw >= n*(2*l+5)
c info   ==  (output) integer.  info = 0 in case of normal computations
c            and 
c            info = -m (<0) - means paramater number m has an illegal value
c            info = 1 - means no convergence achieved (stopping criterion
c            is not fulfilled)
c            info = 2 - means breakdown of the algorithm (try to enlarge
c            parameter l=\ell to get rid of this)
c ----------------------------------------------------------


      info = 0
      
      if (l.gt.lmax .or. l.lt.1) info = -2
      if (tol.le.zero) info = -8
      if (mxmv.lt.0) info = -10
      
      rr = 1
      r = rr+1
      u = r+(l+1)
      xp = u+(l+1)
      bp = xp+1
!      if (bp*n.gt.ldw) info = -12
      
      z = 1
      zz = z+(l+1)
      y0 = zz+(l+1)
      yl = y0+1
      y = yl+1
      
      if (info.ne.0) return

c      initialize first residual
      if (nonzero) then
!         call matvec (n, x, work(1,r) )
         work(1:n,r)=matmul(aamat,x(1:n))
         do ii=1,n
            work(ii,r) = rhs(ii) - work(ii,r)
         enddo
         nmv = 1
      else
         do ii=1,n
            work(ii,r) = rhs(ii)
         enddo
         nmv = 0
      endif

c     initialize iteration loop
      sum1 = zero
      do ii=1,n
         work(ii,rr) = work(ii,r)
         work(ii,bp) = work(ii,r)
         work(ii,xp) = x(ii)
         x(ii) = zero
         sum1=sum1+ work(ii,r)**2 
      enddo

      rnrm0 = sqrt( sum1 )
      rnrm = rnrm0
      if (typestop.eq.'max') then
         maxval1=zero
         do ii=1,n
            maxval1=max( maxval1, abs( work(ii,r) ) )
         enddo
         rnrmmax = maxval1
      endif

      mxnrmx = rnrm0

      mxnrmr = rnrm0
      rcmp = .false.
      xpdt = .false.

      alpha = zero
      omega = one
      sigma = one
      rho0 =  one

c     iterate
      if(typestop.eq.'rel')then
         goon = rnrm.ge.tol*rnrm0 .and. nmv.lt.mxmv
      else if(typestop.eq.'abs')then
         goon = rnrm.ge.tol       .and. nmv.lt.mxmv
      else if(typestop.eq.'max')then
         goon = rnrmmax.ge.tol    .and. nmv.lt.mxmv
      else
         info = -9
         return
      end if

      do while (goon)
c     =====================
c     --- the bicg part ---
c     =====================
      	 rho0 = -omega*rho0
      	 do kk=1,l
            sum1 = zero
            do ii=1,n
               sum1=sum1+ work(ii,rr)*work(ii,r+kk-1) 
            enddo
            rho1 = sum1

      	    if (rho0.eq.zero) then
      	       info = 2
      	       return
      	    endif
      	    beta = alpha*(rho1/rho0)
      	    rho0 = rho1
      	    do jj=0,kk-1
               do ii=1,n
                  work(ii,u+jj) = work(ii,r+jj) - beta*work(ii,u+jj)
               enddo
      	    enddo

!            call matvec(n, work(1,u+kk-1), work(1,u+kk))
            work(1:n,u+kk)=matmul(aamat,work(1:n,u+kk-1))
      	    nmv = nmv+1

            sum1 = zero
            do ii=1,n
               sum1=sum1+ work(ii,rr)*work(ii,u+kk) 
            enddo
            sigma = sum1
            
      	    if (sigma.eq.zero) then
      	       info = 2
      	       return
      	    endif

      	    alpha = rho1/sigma
            do ii=1,n
               x(ii) = alpha*work(ii,u) + x(ii)
            enddo
            
      	    do jj=0,kk-1
               do ii=1,n
                  work(ii,r+jj) = -alpha*work(ii,u+jj+1) + work(ii,r+jj)
               enddo
      	    enddo

!           call matvec (n, work(1,r+kk-1), work(1,r+kk))
            work(1:n,r+kk)=matmul(aamat,work(1:n,r+kk-1))
      	    nmv = nmv+1

            sum1 = zero
            do ii=1,n
               sum1=sum1+ work(ii,r)**2 
            enddo
            rnrm = sqrt( sum1 )
            if (typestop.eq.'max') then
               maxval1=zero
               do ii=1,n
                  maxval1=max( maxval1, abs( work(ii,r) ) )
               enddo
               rnrmmax = maxval1
            endif
      	    mxnrmx = max (mxnrmx, rnrm)
      	    mxnrmr = max (mxnrmr, rnrm)
      	 enddo

c     ==================================
c     --- the convex polynomial part ---
c     ================================== 

c        z = r`r
      	 do i1=1,l+1 
            do jj=i1-1,l
               sum1 = zero
               do ii=1,n
                  sum1=sum1+ work(ii,r+jj)*work(ii,r+i1-1) 
               enddo
               rwork(jj+1,z+i1-1) = sum1
               rwork(z+i1-1,jj+1) = rwork(jj+1,z+i1-1) 
            enddo
      	 enddo 

         do i1=zz,zz+l
            do ii=1,l+1
               rwork(ii,i1)   = rwork(ii,i1+(z-zz)) 
            enddo
         enddo
c        tilde r0 and tilde rl (small vectors)

       	 rwork(1,y0) = -one
         rwork(2,y0) = rwork(2,z) / rwork(2,zz+1)
      	 rwork(l+1,y0) = zero

      	 rwork(1,yl) = zero
         rwork(2,yl) = rwork(2,z+l) / rwork(2,zz+1)
      	 rwork(l+1,yl) = -one

c        convex combination
         do ii=1,l+1
            rwork(ii,y) = zero
         enddo

         do jj=1,l+1
            do ii=1,l+1
               rwork(ii,y) = rwork(ii,y) + rwork(jj,yl)*
     &              rwork(ii,z+jj-1)
            enddo
         enddo
         sum1 = zero
         do ii=1,l+1
            sum1=sum1+ rwork(ii,yl)*rwork(ii,y) 
         enddo
         kappal = sqrt( sum1 )

         do ii=1,l+1
            rwork(ii,y) = zero
         enddo
         do jj=1,l+1
            do ii=1,l+1
               rwork(ii,y) = rwork(ii,y) + rwork(jj,y0)*
     &              rwork(ii,z+jj-1)
            enddo
         enddo
         sum1 = zero
         do ii=1,l+1
            sum1=sum1+ rwork(ii,y0)*rwork(ii,y)  
         enddo
         kappa0 = sqrt( sum1 )

         sum1 = zero
         do ii=1,l+1
            sum1=sum1+ rwork(ii,yl)*rwork(ii,y) 
         enddo
         varrho = sum1  
         varrho = varrho / (kappa0*kappal)
         
      	 hatgamma = sign(1d0,varrho)*max(abs(varrho),7d-1) *
     &        (kappa0/kappal)

         do ii=1,l+1
            rwork(ii,y0) = -hatgamma*rwork(ii,yl) +      rwork(ii,y0)
         enddo

c        update
      	 omega = rwork(l+1,y0)
         do jj=1,l
            do ii=1,n
               work(ii,u) = work(ii,u) - rwork(1+jj,y0)*work(ii,u+jj)
               x(ii)      = x(ii)      + rwork(1+jj,y0)*work(ii,r+jj-1)
               work(ii,r) = work(ii,r) - rwork(1+jj,y0)*work(ii,r+jj)
            enddo
         enddo

         do ii=1,l+1
            rwork(ii,y) = zero
         enddo
         do jj=1,l+1
            do ii=1,l+1
               rwork(ii,y) = rwork(ii,y) + rwork(jj,y0)*
     &              rwork(ii,z+jj-1)
            enddo
         enddo

         sum1 = zero
         do ii=1,l+1
            sum1=sum1+ rwork(ii,y0)*rwork(ii,y) 
         enddo
         rnrm = sqrt( sum1 )

c     ================================
c     --- the reliable update part ---
c     ================================
      	 mxnrmx = max (mxnrmx, rnrm)
      	 mxnrmr = max (mxnrmr, rnrm)
      	 xpdt = (rnrm.lt.delta*rnrm0.and.rnrm0.lt.mxnrmx)
      	 rcmp = ((rnrm.lt.delta*mxnrmr.and.rnrm0.lt.mxnrmr) .or.xpdt)
      	 if (rcmp) then
!      	    call matvec (n, x, work(1,r) )
            work(1:n,r)=matmul(aamat,x(1:n))
            nmv = nmv + 1
            do ii=1,n
               work(ii,r) = work(ii,bp) - work(ii,r)
            enddo
      	    mxnrmr = rnrm
      	    if (xpdt) then
               do ii=1,n
                  work(ii,xp) = x(ii) + work(ii,xp)
                  x(ii) = zero
                  work(ii,bp) = work(ii,r)
               enddo
      	       mxnrmx = rnrm
      	    endif
      	 endif

         if(typestop.eq.'rel')then
            goon = rnrm.ge.tol*rnrm0 .and. nmv.lt.mxmv
            if(okprint) print *, nmv,' ',rnrm/rnrm0
         else if(typestop.eq.'abs')then
            goon = rnrm.ge.tol       .and. nmv.lt.mxmv
            if(okprint) print *, nmv,' ',rnrm
         else if(typestop.eq.'max')then
            maxval1=zero
            do ii=1,n
               maxval1=max(maxval1, abs( work(ii,r) ) )
            enddo
            rnrmmax = maxval1
            goon = rnrmmax.ge.tol    .and. nmv.lt.mxmv
            if(okprint) print *, nmv,' ',rnrmmax
         end if

      enddo

c     =========================
c     --- end of iterations ---
c     =========================

      do ii=1,n
         x(ii) = work(ii,xp) + x(ii)
      enddo

c --------------------- one matvec can be saved by commenting out this:
c                       (this is to compute the true residual)
!     call matvec (n, x, work(1,r) )
      work(1:n,r)=matmul(aamat,x(1:n))
     
      do ii=1,n
         work(ii,r) = rhs(ii) - work(ii,r)   
      enddo

      sum1 = zero
      do ii=1,n
         sum1=sum1+ work(ii,r)**2 
      enddo
      rnrm = sqrt( sum1 )
c --------------------- one matvec can be saved by commenting out this^

      if(typestop.eq.'rel')then
         if (rnrm.gt.tol*rnrm0) info = 1
         tol = rnrm/rnrm0
      else if(typestop.eq.'abs')then
         if (rnrm.gt.tol) info = 1
         tol = rnrm
      else if(typestop.eq.'max')then
         if (rnrmmax.gt.tol) info = 1
         tol = rnrmmax
      end if

      mxmv = nmv

      return
      end
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
       subroutine ecmat(ecm,fun,ndof,nodof)
!       call ecmat(ecm,fun,ndof,nodof)
!
! this subroutine forms the element consistent mass matrix.
!
       implicit none
       integer  ndof,nodof
       real fun(ndof)
       real ecm(ndof,ndof)
       integer  nod,i,j
       real nt(ndof,nodof),tn(nodof,ndof)
       nod=ndof
       nt=0.0
       tn=0.0
       do i=1,nod 
         do j=1,nodof
           nt((i-1)*nodof+j,j)=fun(i)
           tn(j,(i-1)*nodof+j)=fun(i)
         end do
       end do
       ecm=matmul(nt,tn)
      return
      end ! subroutine ecmatdel_s,betatemp_ibm,tsibmin
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine cjgr( x, a, b ,cgtol,cg_max_iter,n) 

       !***************************************************************!
       !*   conjugate gradient solver for a x = b                     *!
       !*   uses the routine a_times(x) below                         *!
       !*   no preconditioner (yet)
       !***************************************************************!

      implicit none
      integer n
      real b(n)
      real cgtol
      
      real a(n,n)
      real x(n)
      real r(n),d(n),q(n),s(n)
      real delta_new,delta_old,eps, alpha, beta
      integer  i,iter,cg_max_iter

      iter = 0

      r = b - matmul(a,x) 
      delta_new = sqrt( dot_product(r,r) ) ! for check
      if (delta_new < cgtol*cgtol ) then
         !      write(99,*) '  cg: no iter required '
      end if

      d = r !m_inverse_times(r)
      delta_new = dot_product(r,d) ! same as above for check
      eps = cgtol*cgtol
!         write(*,*) "  cg: initial residual", delta_new

      do while ((iter.lt.cg_max_iter).and.(delta_new.gt.eps))

       iter=iter+1
!       write(*,*) 'cg...iter, delta_new',iter,delta_new
       q = matmul(a, d )
       alpha = delta_new / dot_product( d, q )
       x = x + alpha * d

       if (mod(iter,75).eq.0) then
          r = b - matmul(a, x ) 
       else
          r = r - alpha * q
       end if

       s = r !m_inverse_times(r)
       delta_old = delta_new
       delta_new = dot_product(r,s)

       beta = delta_new/delta_old
       d = s + beta * d

      end do

      if (iter.eq.cg_max_iter) then
        write(111,*)  "......warning, cg used max iterations"
        write(*,*)    "......itmax = ",iter,", res = ",delta_new
      end if
      write(111,*) "......cg used ",iter," iterations"

      write(*,*) "......cg used ",iter,delta_new  !" iterations"
      end subroutine
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine crossvec(cross,a, b)
       Real*8  a(3), b(3),cross(3),valcross

        cross(1) = a(2) * b(3) - a(3) * b(2)
        cross(2) = a(3) * b(1) - a(1) * b(3)
        cross(3) = a(1) * b(2) - a(2) * b(1)
        valcross=sqrt(cross(1)**2+cross(2)**2+cross(3)**2)+1.0d-9
        cross(1)=cross(1)/valcross
        cross(2)=cross(2)/valcross
        cross(3)=cross(3)/valcross


      END subroutine
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine crossvec_noscale(cross,a, b)
       Real*8  a(3), b(3),cross(3)

        cross(1) = a(2) * b(3) - a(3) * b(2)
        cross(2) = a(3) * b(1) - a(1) * b(3)
        cross(3) = a(1) * b(2) - a(2) * b(1)
      END subroutine
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
       subroutine shapen_dn_ddn(v,w,n,dndv,ddndv2)
       real v,w,u
       real n(12),dndv(12,2),ddndv2(12,3)

      u=1.0-v-w
       n(1:12)= (/
     & u**4/12.0 + (v*u**3)/6.0,
     & u**4/12.0 + (w*u**3)/6.0, 
     & u**4/12.0 + (u**3*v)/2.0  
     &    + (w*u**3)/6.0 + u**2*v**2  
     &    + (w*u**2*v)/2.0 + (u*v**3)/2.0  
     &    + (w*u*v**2)/2.0 + v**4/12.0 + (w*v**3)/6.0, 
     & u**4/2.0 + 2.0*u**3*v + 2.0*u**3*w  
     &    + 2.0*u**2*v**2 + 5.0*u**2*v*w  
     &    + 2.0*u**2*w**2 + (2.0*u*v**3)/3.0  
     &    + 3.0*u*v**2*w + 3.0*u*v*w**2  
     &    + (2.0*u*w**3)/3.0 + v**4/12.0  
     &    + (v**3*w)/2.0 + v**2*w**2  
     &    + (v*w**3)/2.0 + w**4/12.0, 
     & u**4/12.0 + (u**3*w)/2.0  
     &    + (v*u**3)/6.0 + u**2*w**2  
     &    + (v*u**2*w)/2.0 + (u*w**3)/2.0  
     &    + (v*u*w**2)/2.0 + w**4/12.0  
     &    + (v*w**3)/6.0, 
     & v**4/12.0 + (u*v**3)/6.0, 
     & u**4/12.0 + (2.0*u**3*v)/3.0  
     &    + (u**3*w)/2.0 + 2.0*u**2*v**2  

     &    + 3.0*u**2*v*w + u**2*w**2  
     &    + 2.0*u*v**3 + 5.0*u*v**2*w  
     &    + 3.0*u*v*w**2 + (u*w**3)/2.0  
     &    + v**4/2.0 + 2.0*v**3*w  
     &    + 2.0*v**2*w**2 + (2.0*v*w**3)/3.0  
     &    + w**4/12.0, 
     & u**4/12.0 + (u**3*v)/2.0  
     &    + (2.0*u**3*w)/3.0 + u**2*v**2  
     &    + 3.0*u**2*v*w + 2.0*u**2*w**2  
     &    + (u*v**3)/2.0 + 3.0*u*v**2*w  
     &    + 5.0*u*v*w**2 + 2.0*u*w**3  
     &    + v**4/12.0 + (2.0*v**3*w)/3.0  
     &    + 2.0*v**2*w**2 + 2.0*v*w**3  
     &    + w**4/2.0, 
     & w**4/12.0 + (u*w**3)/6.0, 
     & v**4/12.0 + (w*v**3)/6.0, 
     & v**4/12.0 + (v**3*w)/2.0  
     &    + (u*v**3)/6.0 + v**2*w**2  
     &    + (u*v**2*w)/2.0 + (v*w**3)/2.0  
     &    + (u*v*w**2)/2.0 + w**4/12.0  
     &    + (u*w**3)/6.0, 
     & w**4/12.0 + (v*w**3)/6.0 
     & /) 

      dndv(1:12,1)= (/ 
     & - u**3/6.0 - (v*u**2)/2.0, 
     & - u**3/3.0 - (w*u**2)/2.0, 
     & u**3/6.0 + (u**2*v)/2.0   
     &    - (u*v**2)/2 - v**3/6.0, 
     & - 2.0*u**2*v - u**2*w   
     &    - 2.0*u*v**2 - 4.0*u*v*w   
     &    - u*w**2 - v**3/3.0   
     &    - (3.0*v**2*w)/2.0 - v*w**2   
     &    - w**3/6.0, 
     & - u**3/6.0 - u**2*w   
     &    - (v*u**2)/2.0 - (3.0*u*w**2)/2.0   
     &    - v*u*w - w**3/3.0   
     &    - (v*w**2)/2.0, 
     & v**3/6.0 + (u*v**2)/2.0, 
     & u**3/3.0 + 2.0*u**2*v   
     &    + (3.0*u**2*w)/2.0 + 2.0*u*v**2   
     &    + 4.0*u*v*w + u*w**2   
     &    + v**2*w + v*w**2   
     &    + w**3/6.0, 
     & u**3/6.0 + (u**2*v)/2.0   
     &    + u**2*w - (u*v**2)/2.0   
     &    + u*w**2 - v**3/6.0   
     &    - v**2*w - v*w**2, 
     & -w**3/6.0, 
     & v**3/3.0 + (w*v**2)/2.0, 
     & v**3/6.0 + v**2*w   
     &    + (u*v**2)/2.0 + (3*v*w**2)/2.0   
     &    + u*v*w + w**3/3.0   
     &    + (u*w**2)/2.0, 
     & w**3/6.0
     & /)

      dndv(1:12,2)= (/ 
     & - u**3/3.0 - (v*u**2)/2.0,
     & - u**3/6.0 - (w*u**2)/2.0,
     & - u**3/6.0 - u**2*v  
     &    - (w*u**2)/2.0 - (3.0*u*v**2)/2.0
     &    - w*u*v - v**3/3.0  
     &    - (w*v**2)/2.0,
     & - u**2*v - 2.0*u**2*w  
     &    - u*v**2 - 4.0*u*v*w  
     &    - 2.0*u*w**2 - v**3/6.0  
     &    - v**2*w - (3.0*v*w**2)/2.0  
     &    - w**3/3.0,
     & u**3/6.0 + (u**2*w)/2.0  
     &    - (u*w**2)/2.0 - w**3/6.0,
     & -v**3/6.0,
     & u**3/6.0 + u**2*v  
     &    + (u**2*w)/2.0 + u*v**2  
     &    - (u*w**2)/2.0 - v**2*w  
     &    - v*w**2 - w**3/6.0,
     & u**3/3.0 + (3.0*u**2*v)/2.0  
     &    + 2.0*u**2*w + u*v**2  
     &    + 4.0*u*v*w + 2.0*u*w**2  
     &    + v**3/6.0 + v**2*w  
     &    + v*w**2,
     & w**3/6.0 + (u*w**2)/2.0,
     & v**3/6.0,
     & v**3/3.0 + (3.0*v**2*w)/2.0  
     &    + (u*v**2)/2.0 + v*w**2  
     &    + u*v*w + w**3/6.0  
     &    + (u*w**2)/2.0,
     & w**3/3.0 + (v*w**2)/2.0
     & /)

      ddndv2(1:12,1)= (/ 
     & u*v,
     & u**2 + w*u,
     & -2.0*u*v,
     & - 2.0*u**2   
     &    - 2.0*w*u + v**2   
     &    + w*v,
     & u*v + u*w   
     &    + v*w + w**2,
     & u*v,
     & u**2 + w*u   
     &    - 2.0*v**2 - 2.0*w*v,
     & - 2.0*u*v - 2.0*u*w   
     &    - 2.0*v*w - 2.0*w**2,
     & 0.0,
     & v**2 + w*v,
     & u*v + u*w   
     &    + v*w + w**2,
     & 0.0
     & /)



      ddndv2(1:12,2)= (/ 
     & u**2 + v*u,
     & u*w,
     & u*v + u*w  
     &    + v*w + v**2,
     & - 2.0*u**2 - 2.0*v*u  
     &    + w**2 + v*w,
     & - 2.0*u*w,
     & 0.0,
     & - 2.0*u*v - 2.0*u*w  
     &    - 2.0*v*w - 2.0*v**2,
     & u**2 + v*u  
     &    - 2.0*w**2 - 2.0*v*w,
     & u*w,
     & 0.0,
     & u*v + u*w  
     &    + v*w + v**2,
     & w**2 + v*w
     & /)



      ddndv2(1:12,3)= (/ 
     & u**2/2.0 + v*u,
     & u**2/2.0 + w*u,
     & - u**2/2.0 - u*v  
     &    + v**2/2.0,
     & - u**2 + v**2/2.0  
     &    + 2.0*v*w + w**2/2.0,
     & - u**2/2.0 - u*w  
     &    + w**2/2.0,
     & -v**2/2.0,
     & u**2/2.0 - u*w  
     &    - v**2 - 2.0*v*w  
     &    - w**2/2.0,
     & u**2/2.0 - u*v  
     &    - v**2/2.0 - 2.0*v*w  
     &    - w**2,
     & -w**2/2.0,
     & v**2/2.0,
     & v**2/2.0 + 2.0*v*w  
     &    + u*v + w**2/2.0  
     &    + u*w,
     & w**2/2.0
     & /)

      return
      end
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE GMRES_CSR(Wmat,Jmat,Imat,bvec,x,tol,np,nup2)
      IMPLICIT NONE
      INTEGER nup,nhu,nhup1,nup2
      PARAMETER(nup=1000)
      PARAMETER(nhu=1000)
      PARAMETER(nhup1=1001)
      INTEGER i,LL,LLp1,np,maxite
      REAL*8 dnormm,x(np),v(nup,nhup1),bnorm,prod,rho,snormw,tol
      REAL*8 hes(nhup1,nhu),q(2*nhu),bvec(np),r02(nhup1),Wmat(nup2)
      INTEGER Jmat(nup2),Imat(np+1)
 
      REAL*8 dwrk(nup)
      integer noi
      maxite=nhu

      DO i=1,np
!        print*, i,np, bvec(i)
         v(i,1)=bvec(i)
      ENDDO
      bnorm=dnormm(bvec,dwrk,np,nup)
      CALL dscalm(v(1,1),1.D0/bnorm,np,nup)
      CALL dMatzero(hes,nhup1,nhu)
      prod=1.D0
      DO LL=1,nhu
         LLp1=LL+1
CCC   *** Matrix-vector multiplication A*q_n ***
         CALL amux (np, v(1,LL),v(1,LLp1), Wmat,Jmat,Imat,nup,nup2) 
!        CALL MatvecCSC(Wmat,Jmat,Imat,v(1,LL),v(1,LLp1),dwrk,np,nup)
CCC   *** Forming new vector q_(n+1) and new column in "hes" ***
         CALL dorth(v(1,LLp1),v,dwrk,hes,LL,snormw,nhu,nhup1,np,nup)
CCC   *** last element of "hes" ***
         hes(LLp1,LL)=snormw
CCC   *** scale q_(n+1) ***
         CALL dscalm(v(1,LLp1),1.D0/snormw,np,nup)
CCC   *** "hes" becomes "R" and "q" is "Q" (as Givens rotations) ***
         CALL dheqr(hes,nhup1,LL,q,nhu)
CCC   *** residual updated ***
         prod=prod*q(2*LL)
         rho=ABS(prod)
         IF ((rho.LE.tol).OR.(LL.EQ.maxite)) THEN
            CALL dVeczero(r02,nhup1)
CCC      *** vector y formed (called r02) ***
            r02(1)=bnorm
            CALL dhels(hes,nhup1,LL,Q,r02,nhu)
            CALL dVeczero(x,np)
CCC      *** solution x computed as linear combination of basis functions ***
            DO i=1,LL
               CALL daxpy_x(r02(i),v(1,i),x,np,nup)
            ENDDO
            write(*,*) LL,rho
            RETURN
         ENDIF
!         write(*,*) LL,rho
      ENDDO
      write(*,*) 'nhu too small'
      RETURN
      END
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE dhels(A,LDA,N,Q,B,nhu)
      IMPLICIT NONE
      INTEGER IQ,K,KB,KP1,nhu,LDA,N
      REAL*8 A(LDA,nhu),B(LDA),q(2*nhu),C,S,T,T1,T2
      DO K=1,N
         KP1=K+1
         IQ=2*(K-1)+1
         C=q(IQ)
         S=q(IQ+1)
         T1=B(K)
         T2=B(KP1)
         B(K)=C*T1-S*T2
         B(KP1)=S*T1+C*T2
      ENDDO
      DO KB=1,N
         K=N+1-KB
         B(K)=B(K)/A(K,K)
         T=-B(K)
         CALL daxpy(T,A(1,K),B,K-1,LDA)
      ENDDO
      RETURN
      END
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE dheqr(A,LDA,N,Q,nhu)
      IMPLICIT NONE
      INTEGER LDA,N,I,IQ,K,NM1,nhu
      REAL*8 A(LDA,nhu),q(2*nhu),C,S,T,T1,T2
      NM1=N-1
      DO K=1,NM1
         I=2*(K-1)+1
         T1=A(K,N)
         T2=A(K+1,N)
         C=q(I)
         S=q(I+1)
         A(K,N)=C*T1-S*T2
         A(K+1,N)=S*T1+C*T2
      ENDDO
      T1=A(N,N)
      T2=A(N+1,N)
      IF (T2.EQ.0.0D0) THEN
         C=1.0D0
         S=0.0D0
      ELSEIF( ABS(T2).GE.ABS(T1) ) THEN
         T=T1/T2
         S=-1.0D0/DSQRT(1.0D0+T*T)
         C=-S*T
      ELSE
         T=T2/T1
         C=1.0D0/DSQRT(1.0D0+T*T)
         S=-C*T
      ENDIF
      IQ=2*N-1
      q(IQ)=C
      q(IQ+1)=S
      A(N,N)=C*T1-S*T2
      RETURN
      END
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE dorth(vnew,v,dwrk,hes,LL,snormw,nhu,nhup1,n,nup)
      IMPLICIT NONE
      INTEGER n,LL,nup,i,nhu,nhup1
      REAL*8 vnew(nup),v(nup,nhup1),hes(nhup1,nhu),snormw
      REAL*8 arg,sumdsq,tem,vnrm,dnormm,ddot,dwrk(nup)
      vnrm=dnormm(vnew,dwrk,n,nup)
      DO i=1,LL
         hes(i,LL)=ddot(v(1,i),vnew,dwrk,n,nup)
         tem=-hes(i,LL)
         CALL daxpy(tem,v(1,i),vnew,n,nup)
      ENDDO
      snormw=dnormm(vnew,dwrk,n,nup)
      IF (vnrm+0.001D0*snormw.NE.vnrm) RETURN
      sumdsq=0.0D0
      DO i=1,LL
         tem=-ddot(v(1,i),vnew,dwrk,n,nup)
         IF (hes(i,LL)+0.001D0*tem.NE.hes(i,LL)) THEN 
            hes(i,LL)=hes(i,LL)-tem
            CALL daxpy(tem,v(1,i),vnew,n,nup)
            sumdsq=sumdsq+tem**2
         ENDIF
      ENDDO
      IF (sumdsq.NE.0.D0) THEN
         arg=MAX(0.0D0,snormw**2-sumdsq)
         snormw=DSQRT(arg)
      ENDIF
      RETURN
      END
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      REAL*8 FUNCTION dcompsum(x,np,nup)
      IMPLICIT NONE
      INTEGER i,np,nup
      REAL*8 x(nup),temp,err,y
      dcompsum=0.D0
C     *** standard summation of vector ***
      DO i=1,np
         dcompsum=dcompsum+x(i)
      ENDDO
C     *** extra accurate summation of vector ***
ccc      dcompsum=0.D0
ccc      err=0.D0
ccc      DO i=1,np
ccc         temp=dcompsum
ccc         y=x(i)+err
ccc         dcompsum=temp+y
ccc         err=(temp-dcompsum)+y
ccc      ENDDO
ccc      dcompsum=dcompsum+err
      RETURN
      END
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      REAL*8 FUNCTION dnormm(vec,dwrk,n,nup)
      IMPLICIT NONE
      INTEGER i,n,nup
      REAL*8 vec(nup),dwrk(nup),dcompsum
      DO i=1,n
         dwrk(i)=vec(i)*vec(i)
      ENDDO
      dnormm=DSQRT(dcompsum(dwrk,n,nup))
      RETURN
      END
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      REAL*8 FUNCTION ddot(vec1,vec2,dwrk,n,nup)
      IMPLICIT NONE
      INTEGER i,n,nup
      REAL*8 vec1(nup),vec2(nup),dwrk(nup),dcompsum
      DO i=1,n
         dwrk(i)=vec1(i)*vec2(i)
      ENDDO
      ddot=dcompsum(dwrk,n,nup)
      RETURN
      END
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE dscalm(vec,fac,n,nup)
      IMPLICIT NONE
      INTEGER i,n,nup
      REAL*8 vec(nup),fac
      DO i=1,n
         vec(i)=fac*vec(i)
      ENDDO
      RETURN
      END
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE daxpy(a,x,y,n,nup)
      IMPLICIT NONE
      INTEGER i,n,nup
      REAL*8 x(nup),y(nup),a
      DO i=1,n
         y(i)=y(i)+a*x(i)
      ENDDO
      RETURN
      END
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE daxpy_x(a,x,y,n,nup)
      IMPLICIT NONE
      INTEGER i,n,nup
      REAL*8 x(nup),y(n),a
      DO i=1,n
         y(i)=y(i)+a*x(i)
      ENDDO
      RETURN
      END
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE Matvec(A,x,b,dwrk,np,nup)
      IMPLICIT NONE
      INTEGER m,n,np,nup
      REAL*8 A(nup,nup),x(nup),b(nup),dcompsum,dwrk(nup)
      DO m=1,np
         DO n=1,np
            dwrk(n)=A(m,n)*x(n)
         ENDDO
         b(m)=dcompsum(dwrk,np,nup)
      ENDDO
      RETURN
      END
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE dMatzero(M,iup1,iup2)
      IMPLICIT NONE
      INTEGER i,j,iup1,iup2
      REAL*8 M(iup1,iup2)
      DO i=1,iup1
         DO j=1,iup2
            M(i,j)=0.D0
         ENDDO
      ENDDO
      RETURN
      END
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE dVeczero(v,iup)
      IMPLICIT NONE
      INTEGER i,iup
      REAL*8 v(iup)
      DO i=1,iup
         v(i)=0.D0
      ENDDO
      RETURN
      END
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine amux (n, x, y, a,ja,ia,nup,nup2)
      integer nup,nup2,n 
      real*8  x(nup), y(nup), a(nup2) 
      integer ja(nup2), ia(n+1)
c-----------------------------------------------------------------------
c         A times a vector
c----------------------------------------------------------------------- 
c multiplies a matrix by a vector using the dot product form
c Matrix A is stored in compressed sparse row storage.
c
c on entry:
c----------
c n     = row dimension of A
c x     = real array of length equal to the column dimension of
c         the A matrix.
c a, ja,
c    ia = input matrix in compressed sparse row format.
c
c on return:
c-----------
c y     = real array of length n, containing the product y=Ax
c
c-----------------------------------------------------------------------
c local variables
c
      real*8 t
      integer i, k
c-----------------------------------------------------------------------
      do 100 i = 1,n
c
c     compute the inner product of row i with vector x
c 
         t = 0.0d0
         do 99 k=ia(i), ia(i+1)-1 
            t = t + a(k)*x(ja(k))
 99      continue
c
c     store result in y(i) 
c
         y(i) = t
 100  continue
c
      return
c---------end-of-amux---------------------------------------------------
      end
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|

!**************************************************************************
! bicgstab solves the linear system A x = b via the biconjugate           *
! gradient stabilized method.                                             *
!**************************************************************************
! m1 - number of subdiagonal elements                     input           *
! m2 - number of superdiagonal elements                   input           *
! n - number of rows of quadratic matrix                  input           *
! a(1:n,1:m1+m2+1) - known coefficient matrix             input           *
! b(1:n) - known right hand side                          input           *
! g(1:n) - initial guess for the vector x                 input           *
! x(1:n) - calculated x satisfying A x = b                output          *
! itmax - maximum number of iterations                    input           *
! tol - tolerance for convergence criterion               input           *
! precond - used preconditioner                           input           *
! success - 0 if converge, -1 otherwise                   output          *
! itconverge - number of iterations at convergence        output          *
!**************************************************************************
! a contains the matrix A in band compact form. The matrix element        *
! A(i,j) is stored in a(i,m1+1+j-i).                                      *
!                                                                         *
! tol should be something between floating point accuracy and 1.          *
!                                                                         *
! precond=1 - identity matrix as preconditioner                           *
!                                                                         *
! precond=2 - Jacobi matrix as preconditioner                             *
!             'Jacobi matrix' refers to the diagonal part of A.           *
!                                                                         *
! precond=3 - diagonal incomplete LU decomposition as preconditioner      *
!**************************************************************************
       subroutine bicgstab 
     &        (m1,m2,n,a,b,g,x,itmax,tol,precond,success,itconverge)
       implicit none
       integer, intent(in) :: m1,m2,n,itmax,precond
       integer, intent(out) :: success,itconverge
       real(kind=8), intent(in) :: a(1:n,1:m1+m2+1),b(1:n),g(1:n),tol
       real(kind=8), intent(out) :: x(1:n)
       logical :: converge(1:2),inside(1:n,1:m1+m2+1)
     & ,inside_l(1:n,1:m1+m2+1),inside_u(1:n,1:m1+m2+1)
       integer :: i1,i2,i3,i4,i5,i6
       real(kind=8) :: r(1:n),aux0(1:n),rs(1:n)
     & ,rho(0:itmax),scalar,p(0:1,1:n),m(1:n)
     & ,c(0:1,1:n),v(0:1,1:n),ph(1:n),l(1:n,1:m1+m2+1)
     & ,u(1:n,1:m1+m2+1),d(1:n),sy(1-m1:n+m2),sx(1-m1:n+m2)
     & ,sz(1-m1:n+m2),sum_array(1:m1+m2+1),beta(0:itmax)
     & ,alpha(0:itmax),omega(0:itmax),s(1:n),sh(1:n)
     & ,t(1:n),aux1,norm_r,norm_a,norm_x,norm_b,norm_s
       character*30, parameter :: message='breakdown, division by zero'
       success=-1
       converge=.false.
       select case (precond)
       case (2)
       forall(i2=1:n)
       m(i2)=a(i2,m1+1)
       end forall
       case (3)
       inside=.false.
       inside_l=.false.
       inside_u=.false.
       do i1=1,n
       do i2=1,n
       if (m1+1+i2-i1.gt.0.and.m1+1+i2-i1.lt.m1+m2+2) then
       inside(i1,m1+1+i2-i1)=.true.
       endif
       enddo
       enddo
       forall (i1=1:n)
       forall (i2=1:m1)
       inside_l(i1,i2)=.true.
       endforall
       forall (i2=m1+2:m1+m2+1)
       inside_u(i1,i2)=.true.
       endforall
       endforall
       l=0d0
       u=0d0
       where(inside_l.eqv..true.)
       l=a
       endwhere
       where(inside_u.eqv..true.)
       u=a
       endwhere
       d=0d0
       do i1=1,n
       d(i1)=a(i1,m1+1)
       enddo
       do i1=1,n
       d(i1)=1d0/d(i1)
       do i2=m1+2,m1+m2+1
       if (inside(i1,i2).eqv..true.) then
       if (a(i1,i2).ne.0d0.and.a(i1+i2-m1-1,2*m1+2-i2).ne.0d0) then
       d(i1+i2-m1-1)=d(i1+i2-m1-1)
     &              -a(i1+i2-m1-1,2*m1+2-i2)*d(i1)*a(i1,i2)
       endif
       endif
       enddo
       enddo
       end select
       forall (i2=1:n)
       c(0,i2)=g(i2)
       endforall
       aux0=g
       call matrix_product(m1,m2,n,a,aux0)
       forall (i2=1:n)
       r(i2)=b(i2)-aux0(i2)
       endforall
       rs=r
       do i1=1,itmax
       rho(i1-1)=scalar(rs,r,n)
       if (i1.eq.1) then
       forall (i2=1:n)
       p(1,i2)=r(i2)
       end forall
       else
       if (rho(i1-2).eq.0d0.or.omega(i1-1).eq.0d0) then
!       print *, message, '1'
       success=-1; return;! stop
       endif
       beta(i1-1)=(rho(i1-1)/rho(i1-2))*(alpha(i1-1)/omega(i1-1))
       forall (i2=1:n)
       p(1,i2)=r(i2)+beta(i1-1)*(p(0,i2)-omega(i1-1)*v(0,i2))
       endforall
       endif
       select case (precond)
       case (1)
       forall (i2=1:n)
       ph(i2)=p(1,i2)
       endforall
       case (2)
       forall (i2=1:n)
       ph(i2)=p(1,i2)/m(i2)
       endforall
       case (3)
       sz=0d0
       sx=0d0
       forall(i2=1:n)
       sy(i2)=p(1,i2)
       endforall
       do i2=1,n
       sum_array=0d0
       forall (i3=1:m1)
       sum_array(i3)=l(i2,i3)*sz(i2+i3-m1-1)
       endforall
       sz(i2)=1d0/d(i2)*(sy(i2)-sum(sum_array))
       enddo
       do i2=n,1,-1
       sum_array=0d0
       forall (i3=m1+1:m1+m2+1)
       sum_array(i3)=u(i2,i3)*sx(i2+i3-m1-1)
       endforall
       sx(i2)=sz(i2)-1d0/d(i2)*sum(sum_array)
       enddo
       forall(i2=1:n)
       ph(i2)=sx(i2)
       endforall
       end select
       aux0=ph
       call matrix_product(m1,m2,n,a,aux0)
       forall (i2=1:n)
       v(1,i2)=aux0(i2)
       endforall
       forall (i2=1:n)
       aux0(i2)=v(1,i2)
       endforall
       aux1=scalar(rs,aux0,n)
       if (aux1.eq.0d0) then
!       print *, message, '2' !; stop
       success=-1; return;! stop

       endif
       alpha(i1)=rho(i1-1)/aux1
       forall (i2=1:n)
       s(i2)=r(i2)-alpha(i1)*v(1,i2)
       endforall
       norm_s=maxval(abs(s))
       converge(1)=(norm_s.lt.tol)
       if (converge(1)) then
       forall (i2=1:n)
       x(i2)=c(0,i2)+alpha(i1)*ph(i2)
       endforall
       itconverge=i1
       success=0
       exit
       endif
       select case (precond)
       case (1)
       forall (i2=1:n)
       sh(i2)=s(i2)
       endforall
       case (2)
       forall (i2=1:n)
       sh(i2)=s(i2)/m(i2)
       endforall
       case (3)
       sz=0d0
       sx=0d0
       sy=1d0
       forall(i2=1:n)
       sy(i2)=s(i2)
       endforall
       do i2=1,n
       sum_array=0d0
       forall (i3=1:m1)
       sum_array(i3)=l(i2,i3)*sz(i2+i3-m1-1)
       endforall
       sz(i2)=1d0/d(i2)*(sy(i2)-sum(sum_array))
       enddo
       do i2=n,1,-1
       sum_array=0d0
       forall (i3=m1+1:m1+m2+1)
       sum_array(i3)=u(i2,i3)*sx(i2+i3-m1-1)
       endforall
       sx(i2)=sz(i2)-1d0/d(i2)*sum(sum_array)
       enddo
       forall(i2=1:n)
       sh(i2)=sx(i2)
       endforall
       end select
       aux0=sh
       call matrix_product(m1,m2,n,a,aux0)
       forall (i2=1:n)

       t(i2)=aux0(i2)
       endforall
       aux1=scalar(t,t,n)
       if (aux1.eq.0d0) then
!       print *, message, '3'  !; stop
       success=-1; return;! stop
       endif
       omega(i1)=scalar(t,s,n)/aux1
       forall (i2=1:n)
       c(1,i2)=c(0,i2)+alpha(i1)*ph(i2)+omega(i1)*sh(i2)
       r(i2)=s(i2)-omega(i1)*t(i2)
       endforall
       forall (i2=1:n)
       c(0,i2)=c(1,i2)
       p(0,i2)=p(1,i2)
       v(0,i2)=v(1,i2)
       endforall
       forall (i2=1:n)
       aux0(i2)=c(1,i2)
       endforall
       norm_x=maxval(abs(aux0))
       call matrix_product(m1,m2,n,a,aux0)
       aux0=aux0-b!residual
       norm_r=maxval(abs(aux0))
       norm_a=maxval(sum(abs(a),dim=2))
       norm_b=maxval(abs(b))
       converge(2)=norm_r.lt.tol*(norm_a*norm_x+norm_b)
       if (converge(2)) then
       forall (i2=1:n)
       x(i2)=c(1,i2)
       endforall
       itconverge=i1
       success=0
       exit
       endif
       enddo
       end 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
       subroutine matrix_product (m1,m2,n,a,x)
       implicit none
       integer, intent(in) :: m1,m2,n
       real(kind=8), intent(in) :: a(1:n,1:m1+m2+1)
       real(kind=8), intent(inout) :: x(1:n)
       logical :: inside,inside_1,inside_2,inside_3,inside_4
       integer :: i1,i2,i3,i4,matrix_index
       real(kind=8) :: prod(1:n,1:m1+m2+1),xx(-m1-m2-1:n+m1+m2+1)
       xx=0d0
       forall (i2=1:n)
       xx(i2)=x(i2)
       endforall
       prod=0d0
       do i1=1,n
       do i2=1,m1+m2+1
       prod(i1,i2)=a(i1,i2)*xx(i1+i2-m1-1)
       enddo
       enddo
       x=sum(prod,dim=2)
       end 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
       function scalar(a,b,n)
       integer :: n,i
       real(kind=8) :: scalar,a(1:n),b(1:n),prod(1:n)
       forall (i=1:n)
       prod(i)=a(i)*b(i)
       endforall
       scalar=sum(prod)
       end 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine ClosestPtPointEle2Ele(out1,XP1,out2,XP2,X1,X2) 
        real*8 out1(3),XP1(3),out2(3),XP2(3),x1(3,3),x2(3,3)
        integer i,i_2,i2,i2_2,inode,inode2,idotmin,ni
        real*8 pointin(3),outtmp1(3),outin(3)
     &        ,point1(3),point2(3)
     &        ,pointchck1(20,1:3),pointchck2(20,1:3)
     &        ,outinchck1(20,1:3),outinchck2(20,1:3)
     &        ,ddotchck(20),dx(3),ddot

!first loop over triangular points

        ii=0
        do i=1,3
             ii=ii+1
             call ClosestPtPointTriangle(outin,pointin,x1(i,1:3)
     &                      ,x2(1,1:3),x2(2,1:3),x2(3,1:3))
             dx=pointin(1:3)-x1(i,1:3)
             ddot=sqrt(dot_product(dx,dx))
             pointchck2(ii,1:3)=pointin(1:3)
             pointchck1(ii,1:3)=x1(i,1:3)

             ddotchck(ii)=ddot
             outinchck1(ii,1:3)=0; outinchck1(ii,i)=1;
             outinchck2(ii,1:3)=outin(1:3);

        enddo

        do i=1,3
             ii=ii+1
             call ClosestPtPointTriangle(outin,point1,x2(i,1:3)
     &                      ,x1(1,1:3),x1(2,1:3),x1(3,1:3))
             dx=point1(1:3)-x2(i,1:3)
             ddot=sqrt(dot_product(dx,dx))
             pointchck1(ii,1:3)=point1(1:3)
             pointchck2(ii,1:3)=x2(i,1:3)

             ddotchck(ii)=ddot
             outinchck1(ii,1:3)=outin(1:3);
             outinchck2(ii,1:3)=0; outinchck2(ii,i)=1;
        enddo
        do i=1,3
           i_2=i+1
           if(i_2 .gt. 3) i_2=1 
           inode=i_2+1
           if(inode .gt. 3) inode=1 
           do i2=1,3
           ii=ii+1
           i2_2=i2+1
           if(i2_2 .gt. 3) i2_2=1   
           inode2=i2_2+1
           if(inode2 .gt. 3) inode2=1 
              call dist3D_Segment_to_Segment(outtmp1,point1,point2
     &                                    ,x1(i,1:3) ,x1(i_2,1:3)
     &                                    ,x2(i2,1:3),x2(i2_2,1:3))   

             dx=point1(1:3)-point2(1:3)
             ddot=sqrt(dot_product(dx,dx))
             pointchck2(ii,1:3)=point2(1:3)
             pointchck1(ii,1:3)=point1(1:3)
             ddotchck(ii)=ddot

             outinchck1(ii,1:3)=0;
             outinchck1(ii,i_2)= outtmp1(1)
             outinchck1(ii,i)= 1.0d0-outtmp1(1)
             outinchck1(ii,inode)= 0.0

             outinchck2(ii,1:3)=0;
             outinchck2(ii,i2_2)= outtmp1(2)
             outinchck2(ii,i2)= 1.0d0-outtmp1(2)
             outinchck2(ii,inode2)= 0.0
                    
           enddo
        enddo

        ni=ii
        ddotmin=1.0d10; idotmin=-1000
        do i=1,ni
           if(ddotmin .gt. ddotchck(i)) then
               ddotmin=ddotchck(i)
               
               idotmin=i
           endif
        enddo
	out1=outinchck1(idotmin,1:3)
	out2=outinchck2(idotmin,1:3)
	XP1 =pointchck1(idotmin,1:3)
	XP2 =pointchck2(idotmin,1:3)

       return 
      end             
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
       subroutine ClosestPtPointTriangle(out1,point1,p, a, b, c)
       real*8 p(3),a(3),b(3),c(3),ab(3),ac(3),ap(3),bp(3),cp(3)
     &       ,va,vb,vc,u,v,w,d1,d2,d3,d4,d5,d6
     &       ,out1(3),point1(3)

       do i=1,3
          ab(i)=b(i)-a(i)
          ac(i)=c(i)-a(i)
          ap(i)=p(i)-a(i)
       enddo
       d1 = dot_product(ab, ap)
       d2 = dot_product(ac, ap)
       if ((d1 .le. 0.0).and.(d2 .le. 0.0)) then
          point1(1:3)=a(1:3)
          out1(1)=1.0d0   
          out1(2)=0.0d0   
          out1(3)=0.0d0   
          return
       endif

       bp(1:3) = p(1:3) - b(1:3)
       d3 = dot_product(ab, bp)
       d4 = dot_product(ac, bp)
       if ((d3 .ge. 0.0).and.(d4 .le. d3)) then
          point1(1:3)=b(1:3)
          out1(1)=0.0d0   
          out1(2)=1.0d0   
          out1(3)=0.0d0   
          return
       endif

       vc = d1*d4 - d3*d2
       if (vc .le. 0.0 .and. d1 .ge. 0.0 .and. d3 <= 0.0) then
          v = d1 / (d1 - d3)
          point1(1:3)=a(1:3)+v*ab(1:3)
          out1(1)=1.0d0 - v  
          out1(2)=v   
          out1(3)=0.0d0   
          return
       endif

       cp(1:3) = p(1:3) - c(1:3)
       d5 = dot_product(ab, bp)
       d6 = dot_product(ac, bp)               
       if ((d6 .ge. 0.0).and.(d5 .le. d6)) then
          point1(1:3)=c(1:3)
          out1(1)=0.0d0   
          out1(2)=0.0d0   
          out1(3)=1.0d0   
          return
       endif

       vb = d5*d2 - d1*d6
       if (vb .le. 0.0 .and. d2 .ge. 0.0 .and. d6 <= 0.0) then
          w = d2 / (d2 - d6)
          point1(1:3)=a(1:3)+w*ac(1:3)
          out1(1)=1.0d0 - w  
          out1(2)=0.0d0   
          out1(3)=w   
          return
       endif

       va = d3*d6 - d5*d4
       if (va .le. 0.0 .and. (d4 - d3) .ge. 0.0  
     &                 .and. (d5 - d6) .ge. 0.0 ) then
          w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
          point1(1:3)=b(1:3)+w*(c(1:3)-b(1:3))
          out1(1)=0.0d0  
          out1(2)=1.0d0 - w   
          out1(3)=w   
          return
       endif

       denom = 1.0d0 / (va +vb+ vc)
       v = vb * denom 
       w = vc * denom 
       u=1.0-v-w
       point1(1:3)=a(1:3)+ ab(1:3)* v + ac(1:3) * w
!     ! u*a + v*b + w*c, u va * denom = l.Of - v - w
       out1(1)=u  
       out1(2)=v   
       out1(3)=w   
       return 
      end             

c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine dist3D_Segment_to_Segment(out1,point1,point2,
     &                                     S1P0,S1P1,S2P0,S2P1)


       real*8 S1P0(3), S1P1(3), S2P0(3),S2P1(3),u(3),v(3),w(3)
     &        ,SMALL_NUM
     &        ,a,b,c,d,e,DD,sc, sN, sD,tc, tN, tD
     &       ,out1(3),point1(3),point2(3),dP(3)


       SMALL_NUM=0.00000001
       u = S1P1 - S1P0 
       v = S2P1 - S2P0 
       w = S1P0 - S2P0 
       a = dot_product(u,u)          ! always >= 0
       b = dot_product(u,v) 
       c = dot_product(v,v)          ! always >= 0
       d = dot_product(u,w) 
       e = dot_product(v,w) 
       DD = a*c - b*b        ! always >= 0
       sD = DD       !? sc = sN / sD, default sD = D >= 0
       tD = DD       !? tc = tN / tD, default tD = D >= 0

      !! compute the line parameters of the two closest points
      if (DD .lt. SMALL_NUM) then ! the lines are almost parallel
        sN = 0.0          ! force using point P0 on segment S1
        sD = 1.0          ! to prevent possible division by 0.0 later
        tN = e 
        tD = c 
      else
!     get the closest points on the infinite lines
        sN = (b*e - c*d) 
        tN = (a*e - b*d) 
        if (sN < 0.0) then        ! sc < 0 => the s=0 edge is visible
            sN = 0.0 
            tN = e 
            tD = c 
        elseif (sN > sD) then  ! sc > 1  => the s=1 edge is visible
            sN = sD 
            tN = e + b 
            tD = c 
        endif
      endif

      if (tN < 0.0) then            ! tc < 0 => the t=0 edge is visible
        tN = 0.0 
        ! recompute sc for this edge
        if (-d < 0.0) then
            sN = 0.0 
        elseif (-d .gt. a) then
            sN = sD 
        else  
            sN = -d 
            sD = a 
        endif
    
      else if (tN > tD) then      ! tc > 1  => the t=1 edge is visible
        tN = tD 
        ! recompute sc for this edge
        if ((-d + b) < 0.0) then
            sN = 0 
        elseif ((-d + b) > a) then
            sN = sD 
        else  
            sN = (-d +  b);
            sD = a;
        endif
      endif

      ! finally do the division to get sc and tc
      if(abs(sN) < SMALL_NUM) then
        sc =  0.0 
      else
        sc =  sN / sD
      endif
      if(abs(tN)< SMALL_NUM) then
        tc =  0.0 
      else
        tc= tN / tD
      endif

      ! get the difference of the two closest points


      point1(1:3)=  w(1:3) + (sc * u(1:3))+  S2P0(1:3)
      point2(1:3)=  tc * v(1:3)+  S2P0(1:3)

      out1(1)=sc  
      out1(2)=tc   
  


      dP (1:3) = point1(1:3)-point2(1:3) 
      out1(3)=sqrt(dot_product(dP,dP))

      return 
      end
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE lubksb(a,b)
      IMPLICIT NONE
      REAL, DIMENSION(2,2) :: a
      REAL, DIMENSION(2) :: b
      REAL, DIMENSION(2) :: x
      REAL :: det


      det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      x(1)=1.0d0/det*(a(2,2)*b(1) -a(1,2)*b(2))
      x(2)=1.0d0/det*(-a(2,1)*b(1)+a(1,1)*b(2))
      b(1:2)=x(1:2)

       END SUBROUTINE lubksb
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      REAL*8 FUNCTION ContactForce(x,v,MaxP,h,a,b,eta)
      IMPLICIT NONE
!!   @b    p=MaxP
!!   @a    k/2
!!   @k    k=MaxP/(h/2+a-b)
!!   
!!   h=4*dx
!!   a=2*dx
!!   b=dx
!!   x1=[-h:.01:0]-a;
!!   x2=[-a:0.001:-b];
!!   y2=k*h/2+k*(x2+a);
!!   x3=[-b:0.0002:0-eps];
!!   y1=k/2/h*(x1+a+h).^2;
!!   y3=b*k*(h./(2*b)+(a/b-1)+log(-b./x3));

           REAL, INTENT(IN) :: MAxP,x,v
           REAL :: h,a,b,k,eta
           k=MAxP/(h/2.0+a-b)
           if((x.gt.h+a)) then
              ContactForce=0.0
           elseif((x .ge. a).and.(x.le.h+a)) then 
              ContactForce=k/2.0/h*(-x+a+h)**2     
           elseif ((x .ge. b).and.(x.le.a)) then
              ContactForce=k*h/2.0+k*(-x+a)
           else
              ContactForce=b*k*(h/(2.0*b)+(a/b-1.0)+log(b/x))
           endif
            ContactForce=ContactForce+eta*ContactForce*v
      END FUNCTION ContactForce
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE ClosestPointHighOrder(xsi,X0,vpos
     &                               ,cross,xP1,xP,XPVV,ntype,jtype
     &                               ,itype,l,iv,phi
     &                               ,jtypemax,maxordern              
     &                               ,ntrial,tolx,tolf)
      IMPLICIT NONE
      real  l(jtypemax,maxordern)
     *       , iv(jtypemax,maxordern,maxordern) 
     *       , phi(jtypemax,0:2,maxordern,12)

      INTEGER :: itype,ntype,jtypemax,maxordern
      INTEGER, INTENT(IN) :: ntrial,jtype
      REAL, INTENT(IN) :: tolx,tolf
      REAL, DIMENSION(3) :: xsi
      REAL, DIMENSION(3) ::xP1
      
      REAL, DIMENSION(maxordern,3) :: xP,xPVV
      REAL, DIMENSION(maxordern) :: s
      REAL, DIMENSION(2,maxordern) :: dsdv
      REAL, DIMENSION(3,maxordern) :: ddsdvdv
      REAL, DIMENSION(3)  :: X0,cross,vpos
      INTEGER :: i,noj,noi
      REAL :: d
      REAL, DIMENSION(2) :: fvec,p
      REAL, DIMENSION(2,2) :: fjac
      REAL, DIMENSION(3) :: xpos,dx
      REAL, DIMENSION(2,3) :: dxpos
      REAL, DIMENSION(3,3) :: d2xpos

      !print*,',jtypemax,maxordern',jtypemax,maxordern 

      do i=1,ntrial
              !print*,'xsi(2),xsi(3)',xsi(2),xsi(3)  
              call eval_n_dn_ddn(xsi(2),xsi(3), 
     *                             s(1:maxordern),
     *                             dsdv(1,1:maxordern), 
     *                             dsdv(2,1:maxordern), 
     *                             ddsdvdv(1,1:maxordern), 
     *                             ddsdvdv(2,1:maxordern), 
     *                             ddsdvdv(3,1:maxordern), 
     *                             itype, phi, l, iv, jtype
     *                            ,jtypemax,maxordern)

             do noj=1,3
                xpos(noj)=dot_product(s(1:ntype),XP(1:ntype,noj))
                do noi=1,2
                   dxpos(noi,noj)=
     &               dot_product(dsdv(noi,1:ntype),XP(1:ntype,noj))
                end do
                do noi=1,3
                   d2xpos(noi,noj)=
     &             dot_product(ddsdvdv(noi,1:ntype),XP(1:ntype,noj))
                end do
             end do
             dx(1:3)=xP1(1:3)-xpos(1:3)
             fvec(1)=dot_product(dx(1:3),dxpos(1,1:3))
             fvec(2)=dot_product(dx(1:3),dxpos(2,1:3))
             
             fjac(1,1)=-dot_product(dxpos(1,1:3),dxpos(1,1:3))+
     &                  dot_product(dx(1:3),d2xpos(1,1:3))

             fjac(1,2)=-dot_product(dxpos(2,1:3),dxpos(1,1:3))+
     &                  dot_product(dx(1:3),d2xpos(3,1:3))

             fjac(2,1)=fjac(1,2)
             fjac(2,2)=-dot_product(dxpos(2,1:3),dxpos(2,1:3))+
     &                  dot_product(dx(1:3),d2xpos(2,1:3))
             if (sum(abs(fvec)) <= tolf) RETURN
!     Check function convergence.
         p = -fvec    !     Right-hand side of linear equations.

         d=fjac(1,1)*fjac(2,2)-fjac(1,2)*fjac(2,1)

         call lubksb(fjac,p)
         xsi(2:3) = xsi(2:3)+p 
         xsi(1)=1.0d0-xsi(2)-xsi(3)   
         X0(1:3)=xpos(1:3)    

         call crossvec(cross,dxpos(1,1:3), dxpos(2,1:3))
  
         !Update solution.
         if (sum(abs(p)) <= tolx) RETURN ! Check root convergence.
             do noj=1,3
                vpos(noj)=dot_product(s(1:ntype),XPVV(1:ntype,noj))
             enddo
         if((xsi(2) .le. -0.1) .or. 
     &      (xsi(2) .ge. 1.1) .or.
     &      (xsi(3) .le. -0.1) .or.
     &      (xsi(3) .ge. 1.1)  )then
            xsi(1)=-999.0d0

            xsi(2)=-999.0d0
            xsi(3)=-999.0d0   
            X0(1)= -999.0d0     
            X0(2)= -999.0d0 
            X0(3)= -999.0d0    
            return
          endif
      end do
      END SUBROUTINE ClosestPointHighOrder 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      FUNCTION S_Alpha(i,ns,fix_bc1,fix_bc2,istartS,iraiseS)
      integer istartS,iendS,iraiseS,ns,i
      real*8 S_function
      logical fix_bc1,fix_bc2


      iendS=ns+1-istartS-iraiseS
      if(fix_bc1 .and. fix_bc2) then
        S_Alpha=1.0
      elseif(fix_bc1 .and. (.not.(fix_bc2))) then
        S_Alpha=1.0
     &       -S_function(dble(i-iendS)/dble(iraiseS))
      elseif((.not.(fix_bc1)) .and. fix_bc2) then
        S_Alpha=S_function(dble(i-istartS)/dble(iraiseS))
      else
        S_Alpha=S_function(dble(i-istartS)/dble(iraiseS))
     &         -S_function(dble(i-iendS)/dble(iraiseS))
      endif

      return
      end
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
