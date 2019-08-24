
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <REAL.H>
#include <CONSTANTS.H>
#include <SPECIALIZE_F.H>
#include <ArrayLim.H>

      subroutine FORT_FASTCOPY ( &
       dest,DIMS(dest), &
       imin, jmin, kmin,imax, jmax, kmax, &
       src, &
       DIMS(src), &
       imn,  jmn,kmn, &
       ncomp)

      implicit none
      INTEGER_T imin, jmin, kmin, imax, jmax, kmax
      INTEGER_T DIMDEC(dest)
      INTEGER_T imn,  jmn,  kmn
      INTEGER_T DIMDEC(src)
      INTEGER_T ncomp

      REAL_T  dest(DIMV(dest),ncomp)
      REAL_T  src(DIMV(src),ncomp)
      INTEGER_T i,j,k,l,ioff,joff,koff

      ioff=imn-imin
      joff=jmn-jmin
      koff=kmn-kmin

      do l = 1, ncomp
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  dest(i,j,k,l) = src(i+ioff,j+joff,k+koff,l)
               end do
            end do
         end do
      end do

      end

      subroutine FORT_FASTSETVAL ( &
       val,lo,hi, &
       dest,DIMS(dest), &
       ncomp)

      implicit none

      INTEGER_T ncomp
      INTEGER_T lo(3), hi(3)
      INTEGER_T DIMDEC(dest)
      REAL_T  val
      REAL_T  dest(DIMV(dest),ncomp)
      INTEGER_T i,j,k,l
      INTEGER_T imin,jmin,kmin,imax,jmax,kmax

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)

      do l = 1, ncomp
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  dest(i,j,k,l) = val
               end do
            end do
         end do
      end do

      end

      subroutine FORT_FASTZERONORM ( &
       src,DIMS(src), &
       lo,hi,ncomp,nrm)

      implicit none

      INTEGER_T ncomp
      INTEGER_T lo(3), hi(3)
      INTEGER_T DIMDEC(src)
      REAL_T  src(DIMV(src),ncomp)
      REAL_T nrm
      INTEGER_T i,j,k,l
      INTEGER_T imin,jmin,kmin,imax,jmax,kmax

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)

      nrm = 0.0d0

!$omp parallel private(i,j,k,l)
      do l = 1, ncomp
!$omp do reduction(max : nrm)
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  nrm = max(nrm, abs(src(i,j,k,l)))
               end do
            end do
         end do
!$omp end do
      end do
!$omp end parallel

      end

      subroutine FORT_FASTONENORM ( &
       src,DIMS(src), &
       lo,hi,ncomp,nrm)

      implicit none

      INTEGER_T ncomp
      INTEGER_T lo(3), hi(3)
      INTEGER_T DIMDEC(src)
      REAL_T  src(DIMV(src),ncomp)
      REAL_T nrm
      INTEGER_T i,j,k,l
      INTEGER_T imin,jmin,kmin,imax,jmax,kmax

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)

      nrm = 0.0d0

!$omp parallel private(i,j,k,l)
      do l = 1, ncomp
!$omp do reduction(+ : nrm)
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  nrm = nrm + abs(src(i,j,k,l))
               end do ! i
            end do ! j
         end do ! k
!$omp end do
      end do ! l
!$omp end parallel

      end

      subroutine FORT_FASTPLUS( &
       dst,DIMS(dst), &
       imin, jmin, kmin,imax, jmax,kmax, &
       src, &
       DIMS(src), &
       imn,  jmn, kmn, &
       ncomp)

      implicit none
      INTEGER_T imin, jmin, kmin, imax, jmax, kmax
      INTEGER_T DIMDEC(dst)
      INTEGER_T imn,  jmn,  kmn
      INTEGER_T DIMDEC(src)
      INTEGER_T ncomp

      REAL_T  dst(DIMV(dst),ncomp)
      REAL_T  src(DIMV(src),ncomp)
      INTEGER_T i,j,k,l,ioff,joff,koff

      ioff=imn-imin
      joff=jmn-jmin
      koff=kmn-kmin

!$omp parallel private(i,j,k,l)
      do l = 1, ncomp
!$omp do
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  dst(i,j,k,l) = dst(i,j,k,l) + src(i+ioff,j+joff,k+koff,l)
               end do
            end do
         end do
!$omp end do nowait
      end do
!$omp end parallel

      end

      subroutine FORT_FASTMULT ( &
       dst,DIMS(dst), &
       imin, jmin, kmin,imax, jmax,kmax, &
       src, &
       DIMS(src), &
       imn,  jmn,kmn, &
       ncomp)

      implicit none
      INTEGER_T imin, jmin, kmin, imax, jmax, kmax
      INTEGER_T DIMDEC(dst)
      INTEGER_T imn,  jmn,  kmn
      INTEGER_T DIMDEC(src)
      INTEGER_T ncomp

      REAL_T  dst(DIMV(dst),ncomp)
      REAL_T  src(DIMV(src),ncomp)
      INTEGER_T i,j,k,l,ioff,joff,koff

      ioff=imn-imin
      joff=jmn-jmin
      koff=kmn-kmin

!$omp parallel private(i,j,k,l)
      do l = 1, ncomp
!$omp do
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  dst(i,j,k,l) = dst(i,j,k,l) * src(i+ioff,j+joff,k+koff,l)
               end do
            end do
         end do
!$omp end do nowait
      end do
!$omp end parallel

      end
      subroutine FORT_FASTMINUS ( &
       dst, &
       DIMS(dst), &
       imin, jmin,kmin, imax, jmax, kmax, &
       src, &
       DIMS(src), &
       imn,  jmn, kmn, &
       ncomp)

      implicit none
      INTEGER_T imin, jmin, kmin, imax, jmax, kmax
      INTEGER_T DIMDEC(dst)
      INTEGER_T imn,  jmn,  kmn
      INTEGER_T DIMDEC(src)
      INTEGER_T ncomp

      REAL_T  dst(DIMV(dst),ncomp)
      REAL_T  src(DIMV(src),ncomp)
      INTEGER_T i,j,k,l,ioff,joff,koff

      ioff=imn-imin
      joff=jmn-jmin
      koff=kmn-kmin

!$omp parallel private(i,j,k,l)
      do l = 1, ncomp
!$omp do
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  dst(i,j,k,l) = dst(i,j,k,l) - src(i+ioff,j+joff,k+koff,l)
               end do
            end do
         end do
!$omp end do nowait
      end do
!$omp end parallel

      end
      subroutine FORT_FASTDIVIDE ( &
       dst, &
       DIMS(dst), &
       imin, jmin,kmin, imax, jmax,kmax,  &
       src, &
       DIMS(src), &
       imn,  jmn, kmn, &
       ncomp)

      implicit none
      INTEGER_T imin, jmin, kmin, imax, jmax, kmax
      INTEGER_T DIMDEC(dst)
      INTEGER_T imn,  jmn,  kmn
      INTEGER_T DIMDEC(src)
      INTEGER_T ncomp

      REAL_T  dst(DIMV(dst),ncomp)
      REAL_T  src(DIMV(src),ncomp)
      INTEGER_T i,j,k,l,ioff,joff,koff

      ioff=imn-imin
      joff=jmn-jmin
      koff=kmn-kmin

!$omp parallel private(i,j,k,l)
      do l = 1, ncomp
!$omp do
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  dst(i,j,k,l) = dst(i,j,k,l) / src(i+ioff,j+joff,k+koff,l)
               end do
            end do
         end do
!$omp end do nowait
      end do
!$omp end parallel

      end
