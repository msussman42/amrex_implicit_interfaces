
      SUBROUTINE getinfo(whalein,Nodes,whaleout)
      
C      whalefile is the input data file e.g. 'whalenormal.txt'
      
      Character (len=*) whalein, whaleout
      Integer Nodes
      
      whaleout=whalein
      
      print *, whaleout
      
      open(unit=4,file=whalein)
      
      read(4,*) Nodes
      
      close(4)
      
      END SUBROUTINE
      
      
      
      
 
   
      SUBROUTINE runonce(R,S,T,U,V,W,X,Y,Z,List,angle,timestep,spring,
     &              counter,counter_real,Nodes,Cells,
     &              X_init,Y_init,Z_init,whaleout)


C     creates files list.dat and connect.txt
C
C     ALL arguments are output variables used by get_new_geometry. 
C
C     R,S,T -> accel. , U,V,W -> vel. , X,Y,Z -> pos.
C
C     List -> data structure that stores node connectivity
C
C     angle and timestep -> arrays for tail angle data
C
C     spring -> array of spring constants 
C
C     dt*steps = 8


      IMPLICIT NONE

      CHARACTER*35 filename, whaleout
      INTEGER Nodes, Cells, Shape, cycles
      Integer, DIMENSION(Nodes,20) :: List
      INTEGER i, j, k, n1, n2, n3, garbage, p, counter, count, n
      INTEGER, PARAMETER :: STEPS=800
      REAL, PARAMETER :: dt=.01
      Real, Dimension(Nodes) ::  X, Y, Z, R, S, T, X_init, Y_init
      Real, Dimension(Nodes) ::  spring, Z_init
      Real, Dimension(Nodes) ::  U, V, W, V_init, temp_S, temp_T
      Real, Dimension(Nodes) ::  temp_U,temp_V,temp_W,temp_R
      Real, Dimension(22) :: angle, timestep
      Real  a, b, c, sum_x, sum_y, sum_z, mu, nu, value,temp
      Real  counter_real, tailpos, mag, mag_init,L1,L2,L4,L5
      Real  xtailpos

      tailpos=7.77
      xtailpos=0.5

      
      open(unit=2, file= whaleout)
      
      read(2,*) Nodes

      close(2)
      

      call nodelist(Nodes, 20, List)

      open(unit=20, file='list.dat')

      DO i=1,Nodes
         write(20,*) (List(i,k), k=1,20)
      END DO

      Close(20)

      open(unit=2, file= whaleout)      
      
      read(2,*) Nodes
      read(2,*) Cells
      read(2,*) Shape

      

C     Create X, Y, and Z vectors from data file

      DO i=1,Nodes
         read(2,*) garbage, a, b, c
         X(i)=a
         Y(i)=b
         Z(i)=c
         X_init(i)=a
         Y_init(i)=b
         Z_init(i)=c
      END DO   


      open(unit=16, file='connect.txt')

      DO i=1,3
         read(2,*) garbage
      END DO

      DO i=1,Cells-1
         read(2,*) n1
         read(2,*) n2
         read(2,*) n3
         write(16,*) n1, n2, n3

         read(2,*) garbage
         read(2,*) garbage
      END DO 

      read(2,*) n1
      read(2,*) n2
      read(2,*) n3
      write(16,*) n1, n2, n3

      close(unit=16)         
      close(2)
 
C     define angles for tail from angles file

      open(unit=18, file='angles3.txt')

      DO i=1,21
         read(18,*) a
         angle(i)=a*1.5
         timestep(i)=(i-1.0)*(STEPS/20.0)*4.0
      END DO

      close(18)


C     Setup initial velocity and acceleration fields 
C     
C     will start with initial velocity field defined as 
C     the accel in primitive model above


      call tailup(R,S,T,U,V,W,X,Y,Z,Nodes,tailpos)


C     Define spring constants


      call springs(Z, spring, Nodes, tailpos)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  set fluke in place
C
             temp=36.28*0.0174539252/2.0

  
        DO i=1,Nodes

        If (Z_init(i).gt.Z_init(7355)) Then


        If (value.le.0.0) Then


             L1=sqrt((Z_init(i)-Z_init(7355))**2+(Y_init(7355)-Y_init(i)
     &          )**2)


             L4=L1*cos(temp)

             
             L5=sqrt(L1**2-L4**2)


             L2=sqrt(2.0*L5**2)




             Y(i)=Y_init(i)-L2+(Y(7355)-Y_init(7355))
             Z(i)=Z_init(i)-L2+(Z(7355)-Z_init(7355))


        End If

        If (value.gt.0.0) Then

             L1=sqrt((Y_init(i)-Y_init(7355))**2+(Z_init(i)-Z_init(7355)
     &           )**2)

             L4=L1*cos(temp)
 
             L5=sqrt(L1**2-L4**2)

             L2=sqrt(2.0*L5**2)

             Y(i)=Y_init(i)+L2+(Y(7355)-Y_init(7355))
             Z(i)=Z_init(i)-L2+(Z(7355)-Z_init(7355))



        End If


        END IF

        END DO
        
        counter=1
        counter_real=1.0
   
      End Subroutine
      
      
   
      SUBROUTINE new_geometry(R,S,T,U,V,W,X,Y,Z,List,angle,timestep,
     &              spring,counter,counter_real,Nodes,Cells,
     &              X_init,Y_init,Z_init)


C     generates geometry for new time time step


      IMPLICIT NONE

      CHARACTER*35 filename
      INTEGER Nodes, Cells, Shape, cycles
      Integer, DIMENSION(Nodes,20) :: List
      INTEGER i, j, k, n1, n2, n3, garbage, p, counter, count, n
      Integer count1, count2, count3
      INTEGER, PARAMETER :: STEPS=800
      REAL, PARAMETER :: dt=.01
      Real, Dimension(Nodes) ::  X, Y, Z, R, S, T, X_init, Y_init
      Real, Dimension(Nodes) ::  spring, Z_init
      Real, Dimension(Nodes) ::  U, V, W, V_init, temp_S, temp_T
      Real, Dimension(Nodes) ::  temp_U,temp_V,temp_W,temp_R
      Real, Dimension(22) :: angle, timestep
      Real  a, b, c, sum_x, sum_y, sum_z, mu, nu, value,temp
      Real  counter_real, tailpos, mag, mag_init,L1,L2,L4,L5

      tailpos=7.77


C       Here, will force tail nodes to obey tail angle data
      
          call plininterp(timestep, angle, counter_real, 21, value)
 

   
             temp=abs(value)*0.0174539252/2.0

  
        DO i=1,Nodes

        If (Z_init(i).gt.Z_init(7355)) Then


        If (value.le.0.0) Then



             L1=sqrt((Z_init(i)-Z_init(7355))**2+(Y_init(7355)-Y_init(i)
     &          )**2)


             L4=L1*cos(temp)

             
             L5=sqrt(L1**2-L4**2)


             L2=sqrt(2.0*L5**2)


             Y(i)=Y_init(i)-L2+(Y(7355)-Y_init(7355))
             Z(i)=Z_init(i)-L2+(Z(7355)-Z_init(7355))


          End If

        

        If (value.gt.0.0) Then



             L1=sqrt((Y_init(i)-Y_init(7355))**2+(Z_init(i)-Z_init(7355)
     &           )**2)

             L4=L1*cos(temp)
 
             L5=sqrt(L1**2-L4**2)

             L2=sqrt(2.0*L5**2)

             Y(i)=Y_init(i)+L2+(Y(7355)-Y_init(7355))
             Z(i)=Z_init(i)-L2+(Z(7355)-Z_init(7355))


          End If



        END IF

        END DO



C     update accel
C
       mu=0.05
       nu=0.0

       
        DO i=1,Nodes

          If (Z(i).lt.tailpos) Then
           sum_x=0.0
           sum_y=0.0
           sum_z=0.0
           DO j=1,20
              p=List(i,j)
              If (p.ne.0) Then
                mag=sqrt((X(p)-X(i))**2+(Y(p)-Y(i))**2+
     &                (Z(p)-Z(i))**2)
                mag_init=sqrt((X_init(p)-X_init(i))**2+
     &            (Y_init(p)-Y_init(i))**2+(Z_init(p)-Z_init(i))**2)
                sum_x=sum_x+mu*(U(p)-U(i))+spring(i)
     &                 *(X(p)-X(i))*(1.0-mag_init/mag)
                sum_y=sum_y+mu*(V(p)-V(i))+spring(i)
     &                *(Y(p)-Y(i))*(1.0-mag_init/mag)
                sum_z=sum_z+mu*(W(p)-W(i))+spring(i)
     &                *(Z(p)-Z(i))*(1.0-mag_init/mag)
              End If
           END DO
           R(i)=sum_x-nu*U(i)
           S(i)=sum_y-nu*V(i)
           T(i)=sum_z-nu*W(i)
          End If
        
        END DO
           

      
C     Step 1:  Explicitly update velocity

      DO i=1,Nodes
         U(i)=U(i)+dt*R(i)
         V(i)=V(i)+dt*S(i)
         W(i)=W(i)+dt*T(i)
      END DO


C     Step 2:  Explicitly update position

      DO i=1,Nodes
         X(i)=X(i)+dt*U(i)
         Y(i)=Y(i)+dt*V(i)
         Z(i)=Z(i)+dt*W(i)
      END DO


C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    
     
      counter=counter+1
      counter_real=counter_real+1.0
      

      count1=1*STEPS
      count2=2*STEPS
      count3=3*STEPS

      If (count1.eq.counter.or.count3.eq.counter) Then

      DO j=1,Nodes
         U(j)=0.0
         V(j)=-V(j)
         W(j)=0.0
         R(j)=0.0
         S(j)=0.0
         T(j)=0.0
      END DO

      End If

      If (count2.eq.counter) Then

         call taildown(R,S,T,U,V,W,X,Y,Z,Nodes,tailpos)


      END IF

      
      END SUBROUTINE
      
      



    
      SUBROUTINE nodelist(Nodes, Ext, List)

C
C     Creates a data structure containg node connections
C

      IMPLICIT NONE


      INTEGER Nodes, Ext
      INTEGER garbage, kk
      INTEGER N, Cells, Shape, i, j, k, n1, n2, n3
      REAL x, y, z

      INTEGER, DIMENSION(Nodes,Ext) :: List
      INTEGER, DIMENSION(Nodes) :: Counter

      open(unit=12, file='whalenormal.txt')

      read(12,*) garbage
      read(12,*) Cells
      read(12,*) Shape

      DO i=1,Nodes
         read(12,*) N, x, y, z
         List(i,1)=N
         Counter(i)=1
      END DO

      
      DO i=1,3 
         read(12,*) garbage
      END DO

      DO i=1,Cells-1
         read(12,*) n1
         read(12,*) n2
         read(12,*) n3
         
         j=Counter(n1)
         List(n1,j)=n2
         j=j+1
         List(n1,j)=n3
         Counter(n1)=j+1

         j=Counter(n2)
         List(n2,j)=n1
         j=j+1
         List(n2,j)=n3
         Counter(n2)=j+1

         j=Counter(n3)
         List(n3,j)=n1
         j=j+1
         List(n3,j)=n2
         Counter(n3)=j+1

         read(12,*) garbage
         read(12,*) garbage
      END DO


         read(12,*) n1
         read(12,*) n2
         read(12,*) n3

         j=Counter(n1)
         List(n1,j)=n2
         j=j+1
         List(n1,j)=n3
         Counter(n1)=j+1

         j=Counter(n2)
         List(n2,j)=n1
         j=j+1
         List(n2,j)=n3
         Counter(n2)=j+1

         j=Counter(n3)
         List(n3,j)=n1
         j=j+1
         List(n3,j)=n2
         Counter(n3)=j+1

      close(12)


      DO i=1,Nodes
         j=Counter(i)
         DO k=1,j-1
            n1=List(i,k)
            DO kk=k+1,j-1
               n2=List(i,kk)
               If (n1.eq.n2) Then
                  List(i,kk)=0
               End If
            END DO
         END DO
      END DO

      open(unit=10, file='testfile2.dat')

      DO i=1,Nodes
         j=Counter(i)
         write(10,*) i, (List(i,k), k=1,j-1)
      END DO

      close(10)

      END SUBROUTINE
         
               
  


      SUBROUTINE plininterp(xvect, yvect, x, n, value)

C     Performs piecewise linear interpolation of data
C     given in vectors xvect and yvect of length n
C     at a point x.  p(x)=value

      IMPLICIT NONE

      Integer i, j, k, n, m
      Real, Dimension(n) :: xvect, yvect
      Real x, value

      
      DO k=1,n-1
      
         If (x.gt.xvect(k).and.x.lt.xvect(k+1)) Then
            value=yvect(k)+(yvect(k+1)-yvect(k))*(x-xvect(k))
     &            /(xvect(k+1)-xvect(k))
         End If
         
         
         If (x.eq.xvect(k)) Then
            value=yvect(k)
         End if
       

      END DO


      END SUBROUTINE





      SUBROUTINE tailup(R,S,T,U,V,W,X,Y,Z, Nodes, tailpos)

C     Initializes velocity for upward motion. 

      IMPLICIT NONE

      Integer Nodes, i
      Real tailpos
      Real, Dimension(Nodes) :: R,S,T,U,V,W,X,Y,Z

      
      DO i=1,Nodes
         U(i)=0.0
         V(i)=0.0
         V(i)=0.1+0.1*tanh(0.5*Z(i)-0.75)

              If (Z(i).lt.0.0) Then
                 V(i)=V(i)-.01*Z(i)
              End If

             If (Z(i).ge.0.0) Then
                V(i)=V(i)+.003*Z(i)**2
             End If
              
              If (Z(i).ge.tailpos) Then
                 V(i)=V(i)-.006*Z(i)+.006*tailpos
              End If   

         W(i)=0.0

      END DO

      DO i=1,Nodes
         R(i)=0.0
         S(i)=0.0
         T(i)=0.0
      END DO

      
      END SUBROUTINE
      




      SUBROUTINE taildown(R,S,T,U,V,W,X,Y,Z,Nodes,tailpos)

C     Initializes velocity for downward motion

      IMPLICIT NONE

      Integer Nodes, i
      Real tailpos
      Real, Dimension(Nodes) :: R,S,T,U,V,W,X,Y,Z

      
      DO i=1,Nodes
         U(i)=0.0
         V(i)=0.0
         V(i)=-(0.1+0.1*tanh(0.55*Z(i)-2.0))

               If (Z(i).lt.0.0) Then
                  V(i)=V(i)+0.015*Z(i)
               End If

               If (Z(i).gt.0.0) Then
                  V(i)=V(i)-.003*Z(i)**2
               End If

              If (Z(i).gt.tailpos) Then
                   V(i)=V(i)+.004*Z(i)-.004*tailpos
               End If

         W(i)=0.0
         R(i)=0.0
         S(i)=0.0
         T(i)=0.0

      END DO

      END SUBROUTINE





      SUBROUTINE springs(Z, spring, Nodes,tailpos)


C     Defines spring constants for each node

      IMPLICIT NONE

      Integer Nodes, i
      Real tailpos
      Real, Dimension(Nodes) :: spring, Z
      Real const1, const2,const3,springcons


      springcons=0.5

      DO i=1,Nodes
           spring(i)=0.1

          If (Z(i).gt.tailpos) Then
            spring(i)=spring(i)+springcons*(10.0-Z(i))**2
     &                -springcons*(10.0-tailpos)**2
          End If

      END DO
  
            


      END SUBROUTINE
