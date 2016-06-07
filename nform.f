      SUBROUTINE NFORM(RG,QQG,GEP,GEN,GMP,GMN)
C
C  ----------------------------------------------------------
C
C   CALCULATE NUCLEON FORM FACTORS
C
C   IG = 1 - DIPOLE WITH FORM FACTOR SCALING AND GEN=0.0
C      = 2 - IJL FIVE PARAMETER MODEL FIT
C        3 - GEP AND GMP FROM IJL, GMN=DIPOLE, GEN=GALSTER
C            WE CALL THIS THE "BEST FIT" NUCLEON FORM FACTORS
C        4 - BEST FIT EXCEPT GEN = 0.0
C        5 - BLATNIK + ZOVKO VDM FIT
C        6 - JANNSENS 1966 STANDARD FIT
C        7 - DIPOLE + F1N = 0.0
C        8 - GES = 0.0,  GMS = 1.0
C        9 - GES = 1.0,  GMS = 0.0
C       10 - HOHLER1 - PROTON AND NEUTRON FROM FIT 8.2
C       11 - HOHLER2 - PROTON FROM FIT 5.3, NEUTRON FROM 8.2
C       12 - GARI + KRUMPELMANN, Z PHYS. A322,689(1985)
C       13 - KORNER + KURODA, PHYS. REV. D16,2165(1977)
C
C   QQG = INPUT Q SQUARED (GEV**2)
C
C  ------------------------------------------------------------
C
      implicit none

      integer*4 ig
!inputs
      real*8 rg,qqg,gep,gen,gmp,gmn
!parameters
      real*8 gam,br,bw,bf,af
      real*8 rmn2,rmw2,rmf2,rmr2,gamr,pi,rmpi,rmpi2
      real*8 tro,trop,tropp,tfi,tom,tomp
      real*8 rmus,rmuv,bs,bv
      real*8 rmrho2,crho,rkrho,rks,rmomg2,comg,rkomg,rkv
      real*8 rlam12,rlam22,rlamq2
      real*8 vrh1,vrh2,vrh3,vom1,vom2,vom3
      real*8 rmup,rmun
!internal variables
      real*8 qq,tau,gt,t1,t2,alph,top,bot,rho
      real*8 f1s,f2s,f1v,f2v,f1e,f2e,f1m,f2m
      real*8 gd,rs,rv,f1,f2,f3
      real*8 ges,gms,gev,gmv,f1rho,f2rho,f1p,f2p
      real*8 qqp,c1,c2,c3,c4,f2vk,f2sk,f1n,f2n

ccc      COMMON/IO/INPUT,IOUT
C
C
C
C   IJL PARAMETERS FOR 5 PARAMETER DIPOLE FIT (IN GEV UNITS)
C   PHYS LETT. 43B, 191(1973)
      DATA GAM  , BR    , BW   , BF   , AF
     *   /0.25  , 0.672 , 1.102, 0.112, -0.052 /
      DATA RMN2 , RMW2  , RMF2  , RMR2  ,  GAMR  ,  PI
     *    / 0.8817,.6146, 1.0384, 0.5852 ,  .112 , 3.14159  /
      DATA RMPI ,  RMPI2
     *   / .139  , .019321  /
C
C   PARAMETERS FOR BLATNIK AND ZOVKO VDM FIT
C   ACTA PHYSICA AUSTRIACA 39, 62(1974)
C     VECTOR MESON MASSES SQUARED (GEV UNITS)
      DATA   TRO,  TROP,  TROPP,  TFI   ,  TOM  ,  TOMP
     *   / 0.585, 1.30 ,  2.10 ,  1.039 , 0.614 , 1.40  /
C     FITTED PARAMETERS
      DATA  RMUS  , RMUV  ,  BS   ,   BV
     *   / -0.060 , 1.853 , -0.91 , -1.10  /
C
C  PARAMETERS FOR GARI AND KRUMPELMANN, Z.PHYS.A322,689(1985)
      DATA  RMRHO2, CRHO, RKRHO,   RKS, RMOMG2,  COMG, RKOMG, RKV
     *   /  0.6022, 0.377, 6.62, -0.12, 0.6147, 0.411, 0.163, 3.706/
      DATA  RLAM12, RLAM22, RLAMQ2
     *   /   0.632, 5.153,  0.0841/
C
C PARAMETERS FOR KORNER AND KURODA
C    VECTOR MESON MASSES SQUARED USING REGGE PARAMETER ALPHA=1.
      DATA VRH1,   VRH2,   VRH3,   VOM1,   VOM2,   VOM3
     *   / 0.593,  1.593,  2.593,  0.614,  1.614,  2.614/
C
      DATA RMUP   , RMUN
     *   / 2.792782, -1.913148 /
C
ccc      NAMELIST/JUNK/QQG1,QQG2,QQ,QQP,F1,F2,C1,C2,C3,C4,F1V,F2VK,F1S,
ccc     *              F2SK,F1P,F1N,F2P,F2N,IG
C   CONVERT TO FM**-2
      QQ  = QQG/.197328**2
      TAU = QQG/(4.*RMN2)
      IG  = NINT(RG)
      GO TO (110,120,120,120,150,160,170,180,190,200,200,220,230),IG
C
C   DIPOLE
  110 GEP = 1./(1. + QQG/0.71)**2
      GEN = 0.0
      GMP = RMUP*GEP
      GMN = RMUN*GEP
      GO TO 900
C
C   IJL 5 PARAMTER JOB
C
  120 GT = 0.5/(1. + GAM*QQG)**2
      T1 = SQRT(QQG + 4.*RMPI2)
      T2 = SQRT(QQG)
      ALPH = 2.*T1*LOG((T1 + T2)/(2.*RMPI) ) /(T2*PI)
      TOP = RMR2 + 8.*GAMR*RMPI/PI
      BOT = RMR2 + QQG + (4.*RMPI2 + QQG)*GAMR*ALPH/RMPI
      RHO = TOP/BOT
      F1S = GT*((1. - BW - BF) + BW/(1. + QQG/RMW2)
     *                         + BF/(1. + QQG/RMF2)  )
C
      F1V = GT*((1. - BR) + BR*RHO  )
C
      F2S = GT*((-0.12 - AF)/(1. + QQG/RMW2)
     *           + AF/(1. + QQG/RMF2)  )
C
      F2V = GT*(3.706*RHO  )
C
      GEP = F1V + F1S - TAU*(F2V + F2S)
      GEN = F1S - F1V - TAU*(F2S - F2V)
      GMP = F1V + F1S + F2V + F2S
      GMN = F1S - F1V + F2S - F2V
      IF(IG.EQ.2) GO TO 900
      GD = 1./(1.+QQG/.71)**2
      GMN = RMUN*GD
      GEN = -RMUN*TAU*GD/(1. + 5.6*TAU)
      IF(IG.EQ.3) GO TO 900
      GEN = 0.0
      GO TO 900
C
C
C   BLATNIK AND ZOVKO
  150 RS = 1./((1.+QQG/TOM)*(1.+QQG/TFI)*(1.+QQG/TOMP))
      RV = 1./((1.+QQG/TRO)*(1.+QQG/TROP)*(1.+QQG/TROPP))
      F1E = (0.5-TAU*(RMUS+2.*RMN2*BS))*RS
      F2E = (0.5-TAU*(RMUV+2.*RMN2*BV))*RV
      F1M = (0.5+RMUS-0.5*BS*QQG)*RS
      F2M = (0.5+RMUV-0.5*BV*QQG)*RV
      GEP = F1E + F2E
      GMP = F1M + F2M
      GEN = F1E - F2E
      GMN = F1M - F2M
      GO TO 900
C
C   JANNSSENS
  160 F1 = 1. + QQ/15.7
      F2 = 1. + QQ/26.7
      F3 = 1. + QQ/8.19
      GES = 0.5*(2.5/F1 - 1.6/F2 + 0.10)
      GMS = 0.44*(3.33/F1 - 2.77/F2 + 0.44)
      GEV = 0.5*(1.16/F3 - 0.16)
      GMV = 2.353*(1.11/F3 - 0.11)
      GEP = GES + GEV
      GMP = GMS + GMV
      GEN = GES - GEV
      GMN = GMS - GMV
      GO TO 900
C
C   DIPOLE + F1N = 0.0
  170 GEP = 1./(1. + QQG/0.71)**2
      GEN = -RMUN*TAU*GEP
      GMP = RMUP*GEP
      GMN = RMUN*GEP
      GO TO 900
 180  GEP = 0.0
      GEN = 0.0
      GMP = 1.0
      GMN = 0.0
      GO TO 900
C
C
  190 GEP = 1.0
      GEN = 0.0
      GMP = 0.0
      GMN = 0.0
      GO TO 900
C
C  HOHLER1 AND HOHLER2
  200 F1RHO = 0.5*(0.955 + 0.090/(1. + QQG/0.355)**2)/
     *            (1. + QQG/0.536)
C
      F2RHO = 0.5*(5.335 + 0.962/(1. + QQG/0.268))/
     *            (1. + QQG/0.603)
C
      F1S = 0.71/(0.6129 + QQG) - 0.64/(1.0404 + QQG)
     *      - 0.13/(3.240 + QQG)
C
      F1V = F1RHO + 0.05/(1.464 + QQG) - 0.52/(6.0025 + QQG)
     *      + 0.28/(8.7025 + QQG)
C
      F2S = -0.11/(0.6129 + QQG) + 0.13/(1.0404 + QQG)
     *      - 0.02/(3.240 + QQG)
C
      F2V = F2RHO - 1.99/(1.464 + QQG) + 0.20/(6.0025 + QQG)
     *      + 0.19/(8.7025 + QQG)
C
      GEP = F1V + F1S - TAU*(F2V + F2S)
      GEN = F1S - F1V - TAU*(F2S - F2V)
      GMP = F1V + F1S + F2V + F2S
      GMN = F1S - F1V + F2S - F2V
      IF(IG.EQ.10) GO TO 900
C
C   HOHLER2 - USE PROTON FIT 5.3
      F1P = F1RHO + 0.67/(0.6129 + QQG) - 0.39/(0.9216 + QQG)
     *      -0.54/(2.7556 + QQG)
C
      F2P = F2RHO + 0.04/(0.6129 + QQG) - 1.88/(1.2996 + QQG)
     *      + 0.24/(10.1761 + QQG)
      GEP = F1P - TAU*F2P
      GMP = F1P + F2P
C
      GO TO 900
C GARI AND KRUMPELMANN
 220  QQP = QQG*LOG(((RLAM22+QQG)/RLAMQ2))/LOG(RLAM22/RLAMQ2)
      C1 = RLAM12/(RLAM12 + QQP)
      C2 = RLAM22/(RLAM22 + QQP)
      F1 = C1*C2
      F2 = F1*C2
      C3 = RMRHO2/(RMRHO2 + QQG)
      C4 = RMOMG2/(RMOMG2 + QQG)
      F1V = (C3*CRHO + (1-CRHO))*F1
      F2VK = (C3*CRHO*RKRHO + (RKV-CRHO*RKRHO))*F2
      F1S = (C4*COMG + (1-COMG))*F1
      F2SK = (C4*COMG*RKOMG + (RKS-COMG*RKOMG))*F2
      F1P = 0.5*(F1S + F1V)
      F1N = 0.5*(F1S - F1V)
      F2P = 0.5*(F2SK + F2VK)
      F2N = 0.5*(F2SK - F2VK)
      GEP = F1P - TAU*F2P
      GMP = F1P + F2P
      GEN = F1N - TAU*F2N
      GMN = F1N + F2N
      GO TO 900
C
C KORNER AND KURODA
 230  F1S = (1/(1+QQG/VOM1))*(1/(1+QQG/VOM2))
      F1V = (1/(1+QQG/VRH1))*(1/(1+QQG/VRH2))
      F2S = F1S*(1/(1+QQG/VOM3))
      F2V = F1V*(1/(1+QQG/VRH3))
      F1P = 0.5*F1S + 0.5*F1V
      F1N = 0.5*F1S - 0.5*F1V
      F2P = (RMUP-1)*(-0.0335*F2S + 1.0335*F2V)
      F2N =    -RMUN*(-0.0335*F2S - 1.0335*F2V)
      GEP = F1P - TAU*F2P
      GMP = F1P + F2P
      GEN = F1N - TAU*F2N
      GMN = F1N + F2N

 900  RETURN
C
      END
