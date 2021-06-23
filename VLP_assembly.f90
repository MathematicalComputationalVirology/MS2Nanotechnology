      PROGRAM vlp_assembly

        IMPLICIT NONE

        DOUBLE PRECISION :: f,b,f1,b1,f2,b2,felongab,felongcc,a(10840)
        DOUBLE PRECISION :: atot,tout,tmax,time,r1,tau,ftime,vol,kd1,kd2,u1
        DOUBLE PRECISION :: KB,bet,temp,deltG1,deltG2,split1,split2,split3,split4,split5
        DOUBLE PRECISION :: random,dtime

        INTEGER :: i,j,i1,i2,k,l,ii,n,iseed,iout,iseed0
        INTEGER :: AB,CC,AB2,AB3,AB4,AB5,AB5CC1,AB5CC2,AB5CC3,AB5CC4,AB5CC5
        INTEGER :: MnoT4(10,4),MT3D5(15,13),MT3(30,11),MD5(31,25),MD3(30,19)
        INTEGER :: MD3new(30,16),MD3A(16,25),MD3B(17,30),MnoT3(11,6),MT4(45,50)
        
        INTEGER, PARAMETER :: mxg = 700000000

        INTEGER :: nrun
        
        iseed = 679340

        !open (unit=1,file='seed.list',status='unknown')

        DO nrun=1,1

        !read(1,*)iseed
        iseed0 = iseed
        
        time = 0.0d0
        ftime = 0.0d0
        tout = 1.0d-6
        tmax = 0.31d0

        iout = 1
        dtime = 1.0d-6
        
        !!!!!Initial Condition!!!!!
        
        AB = 0
        CC = 10988790
        AB2 = 0
        AB3 = 0
        AB4 = 0
        AB5 = 0
        AB5CC1 = 0
        AB5CC2 = 0
        AB5CC3 = 0
        AB5CC4 = 0
        AB5CC5 = 0
        MnoT3(:,:) = 0
        MnoT4(:,:) = 0
        MT3D5(:,:) = 0
        MT3(:,:) = 0
        MD5(:,:) = 0
        MD3(:,:) = 0
        MD3new(:,:) = 0
        MD3A(:,:) = 0
        MD3B(:,:) = 0
        MT4(:,:) = 0
        a(:) = 0.0d0
        
        !!!!!Parameters Values!!!!!
        
        vol = 7d-16 !in Liter
        
        KB = 1.9872d-3 !Boltzmann constant
        temp = 293.15 !Kelvin temperature
        bet = 1.0d0/KB/temp
        deltG1 = -2.7d0 !kcal/molar
        kd1 = EXP(bet*deltG1)
        deltG2 = -5.4d0 !kcal/molar
        kd2 = EXP(bet*deltG2)
        
        f = 40.0d0 !CC conversion rate to AB in per sec
        b = 0.02d0 !AB conversion rate to CC in per sec
        f1 = 1d+3 * 1.66d-24 / vol !binding rate of AB to AB in the nucleation step
        b1 = 1d+3 * kd1 !unbinding rate of AB from AB in the nucleation step
        f2 = 1d+6 * 1.66d-24 / vol !binding rate of AB to 4AB
        b2 = 1d+6 * kd2 !unbindin rate of AB from 5AB
        felongab = 1d+6 * 1.66d-24 / vol !elongation rate for AB
        felongcc = 1d+6 * 1.66d-24 / vol !elongation rate for CC
        
        split1 = 0.14d0 !0.28d0
        split2 = 0.022d0 !0.022d0
        split3 = 0.031d0 !0.031d0
        split4 = 0.07d0 !0.11d0
        split5 = 0.11d0 !0.11d0
        
        DO WHILE ( time < tmax )
        
          !IF (AB == 0 .and. CC == 0) THEN
          !   GOTO 1
          !ENDIF
          
          !!!!!Transition Rates!!!!!
          
          a(1) = f * DBLE(CC)
          a(2) = b * DBLE(AB)
          a(3) = f1 * DBLE(AB) * (DBLE(AB) - 1.0d0)/2
          a(4) = b1 * DBLE(AB2)
          a(5) = f1 * DBLE(AB) * DBLE(AB2)
          a(6) = 2.0d0 * b1 * DBLE(AB3)
          a(7) = f1 * DBLE(AB3) * DBLE(AB)
          a(8) = 2.0d0 * b1 * DBLE(AB4)
          a(9) = f2 * DBLE(AB4) * DBLE(AB)
          a(10) = 5.0d0 * b2 * DBLE(AB5)
          a(11) = f2 * DBLE(AB5) * DBLE(CC)
          a(12) = b2 * DBLE(AB5CC1)
          a(13) = f2 * DBLE(AB5CC1) * DBLE(CC)
          a(14) = 2.0d0 * b2 * DBLE(AB5CC2)
          a(15) = f2 * DBLE(AB5CC2) * DBLE(CC)
          a(16) = 3.0d0 * b2 * DBLE(AB5CC3)         
          a(17) = f2 * DBLE(AB5CC3) * DBLE(CC)
          a(18) = 4.0d0 * b2 * DBLE(AB5CC4)         
          a(19) = f2 * DBLE(AB5CC4) * DBLE(CC)
          a(20) = 5.0d0 * b2 * DBLE(AB5CC5)
          a(21) = felongab * DBLE(AB5CC5) * DBLE(AB)
          a(22) = split1 * felongcc * DBLE(AB5CC5) * DBLE(CC)
          
          n = 22
          !write(*,*) n
          
          DO i = 1,10
             DO j = 1,3
                a(n+(i-1)*3+j) = felongcc * DBLE(MnoT4(i,j)) * DBLE(CC)
             ENDDO
          ENDDO
          
          n = n+10*3
          !write(*,*) n
          
          DO i = 1,9
             DO j = 1,4
                a(n+(i-1)*4+j) = felongab * DBLE(MnoT4(i,j)) * DBLE(AB)
             ENDDO
          ENDDO
          
          n = n+9*4
          !write(*,*) n
          
          a(n+1) = felongab * DBLE(MnoT4(10,4)) * DBLE(AB)
          a(n+2) = split2 * felongcc * DBLE(MnoT4(10,4)) * DBLE(CC)
          
          n = n+2
          !write(*,*) n
          
          DO i = 1,15
             DO j = 1,12
                a(n+(i-1)*12+j) = felongcc * DBLE(MT3D5(i,j)) * DBLE(CC)
             ENDDO
          ENDDO
          
          n = n+15*12
          !write(*,*) n
          
          DO i = 1,14
             DO j = 1,13
                a(n+(i-1)*13+j) = felongab * DBLE(MT3D5(i,j)) * DBLE(AB)
             ENDDO
          ENDDO
          
          n = n+14*13
          !write(*,*) n
          
          a(n+1) = felongab * DBLE(MT3D5(15,13)) * DBLE(AB)
          a(n+2) = split4 * felongcc * DBLE(MT3D5(15,13)) * DBLE(CC)
          
          n = n+2
          !write(*,*) n
          
          DO i = 1,30
             DO j = 1,10
                a(n+(i-1)*10+j) = felongcc * DBLE(MT3(i,j)) * DBLE(CC)
             ENDDO
          ENDDO
          
          n = n+300
          !write(*,*) n
          
          DO i = 1,29
             DO j = 1,11
                a(n+(i-1)*11+j) = felongab * DBLE(MT3(i,j)) * DBLE(AB)
             ENDDO
          ENDDO
          
          n = n+29*11
          !write(*,*) n
          
          DO i = 1,31
             DO j = 1,24
                a(n+(i-1)*24+j) = felongcc * DBLE(MD5(i,j)) * DBLE(CC)
             ENDDO
          ENDDO
          
          n = n+31*24
          !write(*,*) n
          
          DO i = 1,30
             DO j = 1,25
                a(n+(i-1)*25+j) = felongab * DBLE(MD5(i,j)) * DBLE(AB)
             ENDDO
          ENDDO
          
          n = n+30*25
          !write(*,*) n
          
          DO i = 1,30
             DO j = 1,18
                a(n+(i-1)*18+j) = felongcc * DBLE(MD3(i,j)) * DBLE(CC)
             ENDDO
          ENDDO
          
          n = n+30*18
          !write(*,*) n
          
          DO i = 1,29
             DO j = 1,19
                a(n+(i-1)*19+j) = felongab * DBLE(MD3(i,j)) * DBLE(AB)
             ENDDO
          ENDDO
          
          n = n+29*19
          !write(*,*) n
          
          DO i = 1,11
             DO j = 1,5
                a(n+(i-1)*5+j) = felongcc * DBLE(MnoT3(i,j)) * DBLE(CC)
             ENDDO
          ENDDO
          
          n = n+11*5
          !write(*,*) n
          
          DO i = 1,10
             DO j = 1,6
                a(n+(i-1)*6+j) = felongab * DBLE(MnoT3(i,j)) * DBLE(AB)
             ENDDO
          ENDDO
          
          n = n+10*6
          !write(*,*) n
          
          a(n+1) = split3 * felongcc * DBLE(MnoT3(11,6)) * DBLE(CC)
          a(n+2) = felongab * DBLE(MnoT3(11,6)) * DBLE(AB)
          
          n = n+2
          !write(*,*) n
          
          DO i = 1,45
             DO j = 1,49
                a(n+(i-1)*49+j) = felongcc * DBLE(MT4(i,j)) * DBLE(CC)
             ENDDO
          ENDDO
          
          n = n+45*49
          !write(*,*) n
          
          DO i = 1,44
             DO j = 1,50
                a(n+(i-1)*50+j) = felongab * DBLE(MT4(i,j)) * DBLE(AB)
             ENDDO
          ENDDO
          
          n = n+44*50

          DO i = 1,30
             DO j = 1,15
                a(n+(i-1)*15+j) = felongcc * DBLE(MD3new(i,j)) * DBLE(CC)
             ENDDO
          ENDDO
          
          n = n+30*15
          !write(*,*) n
          
          DO i = 1,29
             DO j = 1,16
                a(n+(i-1)*16+j) = felongab * DBLE(MD3new(i,j)) * DBLE(AB)
             ENDDO
          ENDDO
          
          n = n+29*16
          
          a(n+1) = felongab * DBLE(MD3(30,19)) * DBLE(AB)
          a(n+2) = split5 * felongcc * DBLE(MD3(30,19)) * DBLE(CC)
          
          n = n+2
          !write(*,*) n
          
          DO i = 1,16
             DO j = 1,24
                a(n+(i-1)*24+j) = felongcc * DBLE(MD3A(i,j)) * DBLE(CC)
             ENDDO
          ENDDO
          
          n = n+16*24
          !write(*,*) n
          
          DO i = 1,15
             DO j = 1,25
                a(n+(i-1)*25+j) = felongab * DBLE(MD3A(i,j)) * DBLE(AB)
             ENDDO
          ENDDO
          
          n = n+15*25
          !write(*,*) n
          
          DO i = 1,17
             DO j = 1,29
                a(n+(i-1)*29+j) = felongcc * DBLE(MD3B(i,j)) * DBLE(CC)
             ENDDO
          ENDDO
          
          n = n+17*29
          !write(*,*) n
          
          DO i = 1,16
             DO j = 1,30
                a(n+(i-1)*30+j) = felongab * DBLE(MD3B(i,j)) * DBLE(AB)
             ENDDO
          ENDDO
          
          n = n+16*30
          
          
          
          a(n+1) = b2 * DBLE(MnoT4(1,1))
          a(n+2) = b2 * DBLE(MnoT3(1,1))
          
          a(n+3) = b2 * DBLE(MT3D5(1,1))
          a(n+4) = b2 * DBLE(MD3(1,1))
          
          a(n+5) = b2 * DBLE(MT3(1,1))
          a(n+6) = b2 * DBLE(MD5(1,1))
          
          a(n+7) = b2 * DBLE(MD3A(1,1))
          a(n+8) = b2 * DBLE(MD3B(1,1))
          
          a(n+9) = b2 * DBLE(MD3new(1,1))
          a(n+11) = b2 * DBLE(MT4(1,1))
          
          a(n+11) = felongab * DBLE(MD3new(30,16)) * DBLE(AB)
          a(n+12) = split5 * felongcc * DBLE(MD3new(30,16)) * DBLE(CC)

          !write(*,*) n
          
          !GOTO 1
          
          !=== Total Transition Rate ===!

          atot = 0.0d0

          DO i1=1,10840
             atot = atot + a(i1)
          ENDDO

          !=== Compute Time Increment ===!

          r1 = random(iseed)

          tau = DLOG(1.0d0/r1)
          tau = tau / atot

          time = time + tau

          IF ( time >= tout ) THEN
          
           OPEN(UNIT=3, FILE='VLP.dat', ACCESS='APPEND', STATUS='UNKNOWN')
           WRITE(UNIT=3,FMT='(E16.8,I10,I10,I10,I10,I10,I10)') tout,iseed0,MT3(30,11),MD5(31,25),MD3A(16,25),MD3B(17,30),MT4(45,50)
           CLOSE(UNIT=3)
           !WRITE(16,*) tout,iseed0,MT3(30,11),MD5(31,25),MD3A(16,25),MD3B(17,30),MT4(35,51)

!            iout = iout + 1
!            tout = tout + dtime
            tout = tout + 0.1d0

!            IF ( iout == 10 ) THEN
!              iout = 1
!              IF ( dtime < 1.0d-3 ) dtime = dtime * 10.0d0
!            ENDIF

          ENDIF

          !=== Pick Reaction ===!

          r1 = random(iseed)

          r1 = r1 * atot

          i2 = 0
          atot = 0.0d0

          DO WHILE ( atot < r1 )

            i2 = i2 + 1
            atot = atot + a(i2)

          ENDDO

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          IF ( i2 == 1 ) THEN

            !=== CC to AB ===!

            CC = CC - 1
            AB = AB + 1

          ENDIF

          IF ( i2 == 2 ) THEN

            !=== AB to CC ===!

            CC = CC + 1
            AB = AB - 1

          ENDIF

          IF ( i2 == 3 ) THEN

            !=== AB binds to AB ===!

            AB = AB - 2
            AB2 = AB2 + 1

          ENDIF

          IF ( i2 == 4 ) THEN

            !=== 2AB breaks apart ===!

            AB2 = AB2 - 1
            AB = AB + 2

          ENDIF

          IF ( i2 == 5 ) THEN
          
            !=== 2AB binds to AB ===!

            AB2 = AB2 - 1
            AB = AB - 1
            AB3 = AB3 + 1

          ENDIF

          IF ( i2 == 6 ) THEN

            !=== AB unbinds from 3AB ===!

            AB3 = AB3 - 1
            AB2 = AB2 + 1
            AB = AB + 1

          ENDIF
          
          IF ( i2 == 7 ) THEN
          
            !=== 3AB binds to AB ===!

            AB3 = AB3 - 1
            AB = AB - 1
            AB4 = AB4 + 1

          ENDIF
          
          IF ( i2 == 8 ) THEN
          
            !=== AB unbinds from 4AB ===!

            AB4 = AB4 - 1
            AB3 = AB3 + 1
            AB = AB + 1

          ENDIF
          
          IF ( i2 == 9 ) THEN
          
            !=== 4AB binds to AB ===!

            AB4 = AB4 - 1
            AB = AB - 1
            AB5 = AB5 + 1

          ENDIF
          
          IF ( i2 == 10 ) THEN
          
            !=== AB unbinds from 5AB ===!

            AB5 = AB5 - 1
            AB4 = AB4 + 1
            AB = AB + 1

          ENDIF
          
          IF ( i2 == 11 ) THEN
          
            !=== 5AB binds to CC ===!

            AB5 = AB5 - 1
            CC = CC - 1
            AB5CC1 = AB5CC1 + 1

          ENDIF
          
          IF ( i2 == 12 ) THEN
          
            !=== CC unbinds from 5AB:1CC ===!

            AB5 = AB5 + 1
            CC = CC + 1
            AB5CC1 = AB5CC1 - 1

          ENDIF
          
          IF ( i2 == 13 ) THEN
          
            !=== 5AB:1CC binds to CC ===!

            AB5CC1 = AB5CC1 - 1
            CC = CC - 1
            AB5CC2 = AB5CC2 + 1

          ENDIF
          
          IF ( i2 == 14 ) THEN
          
            !=== CC unbinds from 5AB:2CC ===!

            AB5CC1 = AB5CC1 + 1
            CC = CC + 1
            AB5CC2 = AB5CC2 - 1

          ENDIF
          
          IF ( i2 == 15 ) THEN
          
            !=== 5AB:2CC binds to CC ===!

            AB5CC2 = AB5CC2 - 1
            CC = CC - 1
            AB5CC3 = AB5CC3 + 1

          ENDIF
          
          IF ( i2 == 16 ) THEN
          
            !=== CC unbinds from 5AB:3CC ===!

            AB5CC2 = AB5CC2 + 1
            CC = CC + 1
            AB5CC3 = AB5CC3 - 1

          ENDIF
          
          IF ( i2 == 17 ) THEN
          
            !=== 5AB:3CC binds to CC ===!

            AB5CC3 = AB5CC3 - 1
            CC = CC - 1
            AB5CC4 = AB5CC4 + 1

          ENDIF
          
          IF ( i2 == 18 ) THEN
          
            !=== CC unbinds from 5AB:4CC ===!

            AB5CC3 = AB5CC3 + 1
            CC = CC + 1
            AB5CC4 = AB5CC4 - 1

          ENDIF
          
          IF ( i2 == 19 ) THEN
          
            !=== 5AB:4CC binds to CC ===!

            AB5CC4 = AB5CC4 - 1
            CC = CC - 1
            AB5CC5 = AB5CC5 + 1

          ENDIF
          
          IF ( i2 == 20 ) THEN
          
            !=== CC unbinds from 5AB:5CC ===!

            AB5CC4 = AB5CC4 + 1
            CC = CC + 1
            AB5CC5 = AB5CC5 - 1

          ENDIF
          
          IF ( i2 == 21 ) THEN
          
            !=== 5AB:5CC binds to AB ===!

            AB5CC5 = AB5CC5 - 1
            AB = AB - 1
            MnoT4(1,1) = MnoT4(1,1) + 1

          ENDIF
          
          IF ( i2 == 22 ) THEN
          
            !=== 5AB:5CC binds to CC ===!

            AB5CC5 = AB5CC5 - 1
            CC = CC - 1
            MnoT3(1,1) = MnoT3(1,1) + 1

          ENDIF
          
          IF ( i2 .GE. 23 .AND. i2 .LE. 52 ) THEN
           
            ii = i2 - 22
            k = INT(ii/3)
            l = MOD(ii,3)
            
            IF (l > 0) THEN
            
               !=== CC binds to MnoT4(k+1,l) ===!
               
               MnoT4(k+1,l) = MnoT4(k+1,l) - 1
               CC = CC - 1
               MnoT4(k+1,l+1) = MnoT4(k+1,l+1) + 1
               
            ENDIF
            
            IF (l == 0) THEN
            
               !=== CC binds to MnoT4(k,3) ===!
               
               MnoT4(k,3) = MnoT4(k,3) - 1
               CC = CC - 1
               MnoT4(k,4) = MnoT4(k,4) + 1
               
            ENDIF

          ENDIF
          
          IF ( i2 .GE. 53 .AND. i2 .LE. 88 ) THEN
           
            ii = i2 - 52
            k = INT(ii/4)
            l = MOD(ii,4)
            
            IF (l > 0) THEN
            
               !=== AB binds to MnoT4(k+1,l) ===!
               
               MnoT4(k+1,l) = MnoT4(k+1,l) - 1
               AB = AB - 1
               MnoT4(k+2,l) = MnoT4(k+2,l) + 1
               
            ENDIF
            
            IF (l == 0) THEN
            
               !=== AB binds to MnoT4(k,4) ===!
               
               MnoT4(k,4) = MnoT4(k,4) - 1
               AB = AB - 1
               MnoT4(k+1,4) = MnoT4(k+1,4) + 1
               
            ENDIF

          ENDIF
          
          IF (i2 == 89) THEN
             
             !=== AB binds to 25AB:5CC (MnoT4(20,6)) ===!
             
             MnoT4(10,4) = MnoT4(10,4) - 1
             AB = AB - 1
             MT3D5(1,1) = MT3D5(1,1) + 1
          
          ENDIF
          
          IF (i2 == 90) THEN
             
             !=== CC binds to 25AB:5CC (MnoT4(20,6)) ===!
             
             MnoT4(10,4) = MnoT4(10,4) - 1
             CC = CC - 1
             MD3(1,1) = MD3(1,1) + 1
          
          ENDIF                   
          
          IF ( i2 .GE. 91 .AND. i2 .LE. 270 ) THEN
           
            ii = i2 - 90
            k = INT(ii/12)
            l = MOD(ii,12)
            
            IF (l > 0) THEN
            
               !=== CC binds to MT3D5(k+1,l) ===!
               
               MT3D5(k+1,l) = MT3D5(k+1,l) - 1
               CC = CC - 1
               MT3D5(k+1,l+1) = MT3D5(k+1,l+1) + 1
               
            ENDIF
            
            IF (l == 0) THEN
            
               !=== CC binds to MT3D5(k,12) ===!
               
               MT3D5(k,12) = MT3D5(k,12) - 1
               CC = CC - 1
               MT3D5(k,13) = MT3D5(k,13) + 1
               
            ENDIF

          ENDIF
          
          IF ( i2 .GE. 271 .AND. i2 .LE. 452 ) THEN
           
            ii = i2 - 270
            k = INT(ii/13)
            l = MOD(ii,13)
            
            IF (l > 0) THEN
            
               !=== AB binds to MT3D5(k+1,l) ===!
               
               MT3D5(k+1,l) = MT3D5(k+1,l) - 1
               AB = AB - 1
               MT3D5(k+2,l) = MT3D5(k+2,l) + 1
               
            ENDIF
            
            IF (l == 0) THEN
            
               !=== AB binds to MT3D5(k,13) ===!
               
               MT3D5(k,13) = MT3D5(k,13) - 1
               AB = AB - 1
               MT3D5(k+1,13) = MT3D5(k+1,13) + 1
               
            ENDIF

          ENDIF        
          
          IF (i2 == 453) THEN
             
             !=== AB binds to 30AB:20CC (MT3D5(5,11)) ===!
             
             MT3D5(15,13) = MT3D5(15,13) - 1
             AB = AB - 1
             MT3(1,1) = MT3(1,1) + 1
          
          ENDIF
          
          IF (i2 == 454) THEN
             
             !=== CC binds to 30AB:20CC (MT3D5(5,11)) ===!
             
             MT3D5(15,13) = MT3D5(15,13) - 1
             CC = CC - 1
             MD5(1,1) = MD5(1,1) + 1
          
          ENDIF                             
          
          IF ( i2 .GE. 455 .AND. i2 .LE. 754 ) THEN
           
            ii = i2 - 454
            k = INT(ii/10)
            l = MOD(ii,10)
            
            IF (l > 0) THEN
            
               !=== CC binds to MT3(k+1,l) ===!
               
               MT3(k+1,l) = MT3(k+1,l) - 1
               CC = CC - 1
               MT3(k+1,l+1) = MT3(k+1,l+1) + 1
               
            ENDIF
            
            IF (l == 0) THEN
            
               !=== CC binds to MT3(k,10) ===!
               
               MT3(k,10) = MT3(k,10) - 1
               CC = CC - 1
               MT3(k,11) = MT3(k,11) + 1
               
            ENDIF

          ENDIF
          
          IF ( i2 .GE. 755 .AND. i2 .LE. 1073 ) THEN
           
            ii = i2 - 754
            k = INT(ii/11)
            l = MOD(ii,11)
            
            IF (l > 0) THEN
            
               !=== AB binds to MT3(k+1,l) ===!
               
               MT3(k+1,l) = MT3(k+1,l) - 1
               AB = AB - 1
               MT3(k+2,l) = MT3(k+2,l) + 1
               
            ENDIF
            
            IF (l == 0) THEN
            
               !=== AB binds to MT3(k,11) ===!
               
               MT3(k,11) = MT3(k,11) - 1
               AB = AB - 1
               MT3(k+1,11) = MT3(k+1,11) + 1
               
            ENDIF

          ENDIF

          IF ( i2 .GE. 1074 .AND. i2 .LE. 1817 ) THEN
           
            ii = i2 - 1073
            k = INT(ii/24)
            l = MOD(ii,24)
            
            IF (l > 0) THEN
            
               !=== CC binds to MD5(k+1,l) ===!
               
               MD5(k+1,l) = MD5(k+1,l) - 1
               CC = CC - 1
               MD5(k+1,l+1) = MD5(k+1,l+1) + 1
               
            ENDIF
            
            IF (l == 0) THEN
            
               !=== CC binds to MD5(k,24) ===!
               
               MD5(k,24) = MD5(k,24) - 1
               CC = CC - 1
               MD5(k,25) = MD5(k,25) + 1
               
            ENDIF

          ENDIF
          
          IF ( i2 .GE. 1818 .AND. i2 .LE. 2567 ) THEN
           
            ii = i2 - 1817
            k = INT(ii/25)
            l = MOD(ii,25)
            
            IF (l > 0) THEN
            
               !=== AB binds to MD5(k+1,l) ===!
               
               MD5(k+1,l) = MD5(k+1,l) - 1
               AB = AB - 1
               MD5(k+2,l) = MD5(k+2,l) + 1
               
            ENDIF
            
            IF (l == 0) THEN
            
               !=== AB binds to MD5(k,25) ===!
               
               MD5(k,25) = MD5(k,25) - 1
               AB = AB - 1
               MD5(k+1,25) = MD5(k+1,25) + 1
               
            ENDIF

          ENDIF                  
          
          IF ( i2 .GE. 2568 .AND. i2 .LE. 3107 ) THEN
           
            ii = i2 - 2567
            k = INT(ii/18)
            l = MOD(ii,18)
            
            IF (l > 0) THEN
            
               !=== CC binds to MD3(k+1,l) ===!
               
               MD3(k+1,l) = MD3(k+1,l) - 1
               CC = CC - 1
               MD3(k+1,l+1) = MD3(k+1,l+1) + 1
               
            ENDIF
            
            IF (l == 0) THEN
            
               !=== CC binds to MD3(k,16) ===!
               
               MD3(k,18) = MD3(k,18) - 1
               CC = CC - 1
               MD3(k,19) = MD3(k,19) + 1
               
            ENDIF

          ENDIF
          
          IF ( i2 .GE. 3108 .AND. i2 .LE. 3658 ) THEN
           
            ii = i2 - 3107
            k = INT(ii/19)
            l = MOD(ii,19)
            
            IF (l > 0) THEN
            
               !=== AB binds to MD3(k+1,l) ===!
               
               MD3(k+1,l) = MD3(k+1,l) - 1
               AB = AB - 1
               MD3(k+2,l) = MD3(k+2,l) + 1
               
            ENDIF
            
            IF (l == 0) THEN
            
               !=== AB binds to MD3(k,17) ===!
               
               MD3(k,19) = MD3(k,19) - 1
               AB = AB - 1
               MD3(k+1,19) = MD3(k+1,19) + 1
               
            ENDIF

          ENDIF         

          IF ( i2 .GE. 3659 .AND. i2 .LE. 3713 ) THEN
           
            ii = i2 - 3658
            k = INT(ii/5)
            l = MOD(ii,5)
            
            IF (l > 0) THEN
            
               !=== CC binds to MnoT3(k+1,l) ===!
               
               MnoT3(k+1,l) = MnoT3(k+1,l) - 1
               CC = CC - 1
               MnoT3(k+1,l+1) = MnoT3(k+1,l+1) + 1
               
            ENDIF
            
            IF (l == 0) THEN
            
               !=== CC binds to MnoT3(k,4) ===!
               
               MnoT3(k,5) = MnoT3(k,5) - 1
               CC = CC - 1
               MnoT3(k,6) = MnoT3(k,6) + 1
               
            ENDIF

          ENDIF
          
          IF ( i2 .GE. 3714 .AND. i2 .LE. 3773 ) THEN
           
            ii = i2 - 3713
            k = INT(ii/6)
            l = MOD(ii,6)
            
            IF (l > 0) THEN
            
               !=== AB binds to MnoT3(k+1,l) ===!
               
               MnoT3(k+1,l) = MnoT3(k+1,l) - 1
               AB = AB - 1
               MnoT3(k+2,l) = MnoT3(k+2,l) + 1
               
            ENDIF
            
            IF (l == 0) THEN
            
               !=== AB binds to MnoT3(k,5) ===!
               
               MnoT3(k,6) = MnoT3(k,6) - 1
               AB = AB - 1
               MnoT3(k+1,6) = MnoT3(k+1,6) + 1
               
            ENDIF

          ENDIF

          IF (i2 == 3774) THEN
             
             !=== CC binds to 25AB:10CC (MnoT3(21,5)) ===!
             
             MnoT3(11,6) = MnoT3(11,6) - 1
             CC = CC - 1
             MD3new(1,1) = MD3(1,1) + 1
          
          ENDIF
          
          IF (i2 == 3775) THEN
             
             !=== AB binds to 25AB:10CC (MnoT3(21,5)) ===!
             
             MnoT3(11,6) = MnoT3(11,6) - 1
             AB = AB - 1
             MT4(1,1) = MT4(1,1) + 1
          
          ENDIF

          IF ( i2 .GE. 3776 .AND. i2 .LE. 5980 ) THEN
           
            ii = i2 - 3775
            k = INT(ii/49)
            l = MOD(ii,49)
            
            IF (l > 0) THEN
            
               !=== CC binds to MT4(k+1,l) ===!
               
               MT4(k+1,l) = MT4(k+1,l) - 1
               CC = CC - 1
               MT4(k+1,l+1) = MT4(k+1,l+1) + 1
               
            ENDIF
            
            IF (l == 0) THEN
            
               !=== CC binds to MT4(k,50) ===!
               
               MT4(k,49) = MT4(k,49) - 1
               CC = CC - 1
               MT4(k,50) = MT4(k,50) + 1
               
            ENDIF

          ENDIF
          
          IF ( i2 .GE. 5981 .AND. i2 .LE. 8180 ) THEN
           
            ii = i2 - 5980
            k = INT(ii/50)
            l = MOD(ii,50)
            
            IF (l > 0) THEN
            
               !=== AB binds to MT4(k+1,l) ===!
               
               MT4(k+1,l) = MT4(k+1,l) - 1
               AB = AB - 1
               MT4(k+2,l) = MT4(k+2,l) + 1
               
            ENDIF
            
            IF (l == 0) THEN
            
               !=== AB binds to MT4(k,50) ===!
               
               MT4(k,50) = MT4(k,50) - 1
               AB = AB - 1
               MT4(k+1,50) = MT4(k+1,50) + 1
               
            ENDIF

          ENDIF

          IF ( i2 .GE. 8181 .AND. i2 .LE. 8630 ) THEN
           
            ii = i2 - 8180
            k = INT(ii/15)
            l = MOD(ii,15)
            
            IF (l > 0) THEN
            
               !=== CC binds to MT4(k+1,l) ===!
               
               MD3new(k+1,l) = MD3new(k+1,l) - 1
               CC = CC - 1
               MD3new(k+1,l+1) = MD3new(k+1,l+1) + 1
               
            ENDIF
            
            IF (l == 0) THEN
            
               !=== CC binds to MT4(k,50) ===!
               
               MD3new(k,15) = MD3new(k,15) - 1
               CC = CC - 1
               MD3new(k,16) = MD3new(k,16) + 1
               
            ENDIF

          ENDIF
          
          IF ( i2 .GE. 8631 .AND. i2 .LE. 9094 ) THEN
           
            ii = i2 - 8630
            k = INT(ii/16)
            l = MOD(ii,16)
            
            IF (l > 0) THEN
            
               !=== AB binds to MT4(k+1,l) ===!
               
               MD3new(k+1,l) = MD3new(k+1,l) - 1
               AB = AB - 1
               MD3new(k+2,l) = MD3new(k+2,l) + 1
               
            ENDIF
            
            IF (l == 0) THEN
            
               !=== AB binds to MT4(k,50) ===!
               
               MD3new(k,16) = MD3new(k,16) - 1
               AB = AB - 1
               MD3new(k+1,16) = MD3new(k+1,16) + 1
               
            ENDIF

          ENDIF

          IF (i2 == 9095) THEN
             
             !=== AB binds to MD3(30,19) ===!
             
             MD3(30,19) = MD3(30,19) - 1
             AB = AB - 1
             MD3A(1,1) = MD3A(1,1) + 1
          
          ENDIF
          
          IF (i2 == 9096) THEN
             
             !=== CC binds to MD3(30,19) ===!
             
             MD3(30,19) = MD3(30,19) - 1
             CC = CC - 1
             MD3B(1,1) = MD3B(1,1) + 1
          
          ENDIF 

          IF ( i2 .GE. 9097 .AND. i2 .LE. 9480 ) THEN
           
            ii = i2 - 9096
            k = INT(ii/24)
            l = MOD(ii,24)
            
            IF (l > 0) THEN
            
               !=== CC binds to MD3A(k+1,l) ===!
               
               MD3A(k+1,l) = MD3A(k+1,l) - 1
               CC = CC - 1
               MD3A(k+1,l+1) = MD3A(k+1,l+1) + 1
               
            ENDIF
            
            IF (l == 0) THEN
            
               !=== CC binds to MD3A(k,24) ===!
               
               MD3A(k,24) = MD3A(k,24) - 1
               CC = CC - 1
               MD3A(k,25) = MD3A(k,25) + 1
               
            ENDIF

          ENDIF
          
          IF ( i2 .GE. 9481 .AND. i2 .LE. 9855 ) THEN
           
            ii = i2 - 9480
            k = INT(ii/25)
            l = MOD(ii,25)
            
            IF (l > 0) THEN
            
               !=== AB binds to MD3A(k+1,l) ===!
               
               MD3A(k+1,l) = MD3A(k+1,l) - 1
               AB = AB - 1
               MD3A(k+2,l) = MD3A(k+2,l) + 1
               
            ENDIF
            
            IF (l == 0) THEN
            
               !=== AB binds to MD3A(k,25) ===!
               
               MD3A(k,25) = MD3A(k,25) - 1
               AB = AB - 1
               MD3A(k+1,25) = MD3A(k+1,25) + 1
               
            ENDIF

          ENDIF

          IF ( i2 .GE. 9856 .AND. i2 .LE. 10348 ) THEN
           
            ii = i2 - 9855
            k = INT(ii/29)
            l = MOD(ii,29)
            
            IF (l > 0) THEN
            
               !=== CC binds to MD3B(k+1,l) ===!
               
               MD3B(k+1,l) = MD3B(k+1,l) - 1
               CC = CC - 1
               MD3B(k+1,l+1) = MD3B(k+1,l+1) + 1
               
            ENDIF
            
            IF (l == 0) THEN
            
               !=== CC binds to MD3B(k,29) ===!
               
               MD3B(k,29) = MD3B(k,29) - 1
               CC = CC - 1
               MD3B(k,30) = MD3B(k,30) + 1
               
            ENDIF

          ENDIF
          
          IF ( i2 .GE. 10349 .AND. i2 .LE. 10828 ) THEN
           
            ii = i2 - 10348
            k = INT(ii/30)
            l = MOD(ii,30)
            
            IF (l > 0) THEN
            
               !=== AB binds to MD3B(k+1,l) ===!
               
               MD3B(k+1,l) = MD3B(k+1,l) - 1
               AB = AB - 1
               MD3B(k+2,l) = MD3B(k+2,l) + 1
               
            ENDIF
            
            IF (l == 0) THEN
            
               !=== AB binds to MD3B(k,30) ===!
               
               MD3B(k,30) = MD3B(k,30) - 1
               AB = AB - 1
               MD3B(k+1,30) = MD3B(k+1,30) + 1
               
            ENDIF

          ENDIF

          IF (i2 == 10829) THEN
             
             !===  ===!
             
             MnoT4(1,1) = MnoT4(1,1) - 1
             AB = AB + 1
             AB5CC5 = AB5CC5 + 1
          
          ENDIF
          
          IF (i2 == 10830) THEN
             
             !===  ===!
             
             MnoT3(1,1) = MnoT3(1,1) - 1
             CC = CC + 1
             AB5CC5 = AB5CC5 + 1
          
          ENDIF
          
          IF (i2 == 10831) THEN
             
             !===  ===!
             
             MT3D5(1,1) = MT3D5(1,1) - 1
             AB = AB + 1
             MnoT4(10,4) = MnoT4(10,4) + 1
          
          ENDIF
          
          IF (i2 == 10832) THEN
             
             !===  ===!
             
             MD3(1,1) = MD3(1,1) - 1
             CC = CC + 1
             MnoT4(10,4) = MnoT4(10,4) + 1
          
          ENDIF         
          
          IF (i2 == 10833) THEN
             
             !===  ===!
             
             MT3(1,1) = MT3(1,1) - 1
             AB = AB + 1
             MT3D5(15,13) = MT3D5(15,13) + 1
          
          ENDIF
          
          IF (i2 == 10834) THEN
             
             !===  ===!
             
             MD5(1,1) = MD5(1,1) - 1
             CC = CC + 1
             MT3D5(15,13) = MT3D5(15,13) + 1
          
          ENDIF
          
          IF (i2 == 10835) THEN
             
             !===  ===!
             
             MD3A(1,1) = MD3A(1,1) - 1
             AB = AB + 1
             MD3(30,19) = MD3(30,19) + 1
          
          ENDIF
          
          IF (i2 == 10836) THEN
             
             !===  ===!
             
             MD3B(1,1) = MD3B(1,1) - 1
             CC = CC + 1
             MD3(30,19) = MD3(30,19) + 1
          
          ENDIF
          
          IF (i2 == 10837) THEN
             
             !===  ===!
             
             MD3new(1,1) = MD3new(1,1) - 1
             CC = CC + 1
             MnoT3(11,6) = MnoT3(11,6) + 1
          
          ENDIF
          
          IF (i2 == 10838) THEN
             
             !===  ===!
             
             MT4(1,1) = MT4(1,1) - 1
             AB = AB + 1
             MnoT3(11,6) = MnoT3(11,6) + 1
          
          ENDIF  
          
          IF (i2 == 10839) THEN
             
             !=== AB binds to MD3new(30,16) ===!
             
             MD3new(30,16) = MD3new(30,16) - 1
             AB = AB - 1
             MD3A(1,1) = MD3A(1,1) + 1
          
          ENDIF
          
          IF (i2 == 10840) THEN
             
             !=== CC binds to MD3new(30,16) ===!
             
             MD3new(30,16) = MD3new(30,16) - 1
             CC = CC - 1
             MD3B(1,1) = MD3B(1,1) + 1
          
          ENDIF        

        ENDDO

1       CONTINUE

        write(*,*)time,iseed0,MT3(30,11),MD5(31,25),MD3A(16,25),MD3B(17,30),MT4(45,50)

        ENDDO
      
      END PROGRAM
      
      !!!!! Creating a Random Number !!!!!
      
      DOUBLE PRECISION FUNCTION RANDOM (ISEED)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        INTEGER, INTENT(INOUT) :: iseed

        !=== VARIABLES ===!

        INTEGER :: hi,lo,test

        INTEGER, PARAMETER :: ae = 16807
        INTEGER, PARAMETER :: m = 2147483647
        INTEGER, PARAMETER :: q = 127773
        INTEGER, PARAMETER :: re = 2836


        hi = INT(iseed/q)
        lo = MODULO(iseed,q)

        test = ae * lo - re * hi

        IF ( test > 0 ) THEN

          iseed = test

        ELSE

          iseed = test + m

        ENDIF

        RANDOM = DBLE(iseed) / DBLE(m)

        RETURN

      END FUNCTION RANDOM
