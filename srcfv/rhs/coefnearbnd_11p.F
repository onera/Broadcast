  !  coef are obtained with sympy for FD scheme
  !face 9/2:
  !{a8: -3/56, a1: -1/63, a5: 6/5, a7: 4/21, a2: 3/28, a9: 1/105, a10: -1/1260, a6: -1/2, a3: -4/7, a0: 1/840, a4: -11/30}
  a0  =   1.d0/840.d0
  a1  = - 1.d0/63.d0
  a2  =   3.d0/28.d0
  a3  = - 4.d0/7.d0
  a4  = -11.d0/30.d0
  a5  =   6.d0/5.d0
  a6  = - 0.5d0
  a7  =   4.d0/21.d0
  a8  = - 3.d0/56.d0
  a9  =   1.d0/105.d0
  a10 = - 1.d0/1260.d0
  ! Phi(9/2) = Phi(11/2) - Sum( ai*fi)
  c9d0  = c5 - a0
  c9d1  = c4 - a1
  c9d2  = c3 - a2
  c9d3  = c2 - a3
  c9d4  = c1 - a4
  c9d5  = c1 - a5
  c9d6  = c2 - a6
  c9d7  = c3 - a7
  c9d8  = c4 - a8
  c9d9  = c5 - a9
  c9d10 = -a10
  ! face 7/2
  ! {a8: 3/40, a1: 1/24, a5: -21/20, a7: -1/4, a2: -3/8, a9: -1/72, a10: 1/840, a6: 7/12, a3: -319/420, a0: -1/360, a4: 7/4}
  ! FD37
  a0  = -  1.d0/360.d0
  a1  =    1.d0/24.d0
  a2  = -  3.d0/8.d0
  a3  = -319.d0/420.d0
  a4  =    1.75d0
  a5  = - 21.d0/20.d0
  a6  =    7.d0/12.d0
  a7  = -  0.25d0
  a8  =    0.075d0
  a9  = -  1.d0/72.d0
  a10 =    1.d0/840.d0
  ! Phi(7/2) = Phi(9/2) - Sum( ai*fi)
  c7d0  = c9d0  - a0
  c7d1  = c9d1  - a1
  c7d2  = c9d2  - a2
  c7d3  = c9d3  - a3
  c7d4  = c9d4  - a4
  c7d5  = c9d5  - a5
  c7d6  = c9d6  - a6
  c7d7  = c9d7  - a7
  c7d8  = c9d8  - a8
  c7d9  = c9d9  - a9
  c7d10 = c9d10 - a10
  ! face 5/2
  ! {a8: -1/6, a1: -2/9, a5: 28/15, a7: 8/15, a2: -341/280, a9: 2/63, a10: -1/360, a6: -7/6, a3: 8/3, a0: 1/90, a4: -7/3}
  ! FD28
  a0=    1.d0/90.d0
  a1= -  2.d0/9.d0
  a2= -341.d0/280.d0
  a3=    8.d0/3.d0
  a4= -  7.d0/3.d0
  a5=   28.d0/15.d0
  a6= -  7.d0/6.d0
  a7=    8.d0/15.d0
  a8= -  1.d0/6.d0
  a9=    2.d0/63.d0
  a10= - 1.d0/360.d0
  ! Phi(5/2) = Phi(7/2) - Sum( ai*fi) given in Berland
  c5d0  = c7d0  - a0
  c5d1  = c7d1  - a1
  c5d2  = c7d2  - a2
  c5d3  = c7d3  - a3
  c5d4  = c7d4  - a4
  c5d5  = c7d5  - a5
  c5d6  = c7d6  - a6
  c5d7  = c7d7  - a7
  c5d8  = c7d8  - a8
  c5d9  = c7d9  - a9
  c5d10 = c7d10 - a10
  ! face 3/2
  ! {a8: 9/14, a1: -4609/2520, a5: -63/10, a7: -2, a2: 9/2, a9: -1/8, a10: 1/90, a6: 21/5, a3: -6, a0: -1/10, a4: 7}
  ! FD19
  a0= -   0.1d0
  a1= -4609.d0/2520.d0
  a2=     4.5d0
  a3= -   6.d0
  a4=     7.d0
  a5= -  6.3d0
  a6=    21.d0*0.2d0
  a7= -   2.d0
  a8=     9.d0/14.d0
  a9= -   0.125d0
  a10=    1.d0/90.d0
  ! Phi(3/2) = Phi(5/2) - Sum( ai*fi) given in Berland
  c3d0  = c5d0  - a0
  c3d1  = c5d1  - a1
  c3d2  = c5d2  - a2
  c3d3  = c5d3  - a3
  c3d4  = c5d4  - a4
  c3d5  = c5d5  - a5
  c3d6  = c5d6  - a6
  c3d7  = c5d7  - a7
  c3d8  = c5d8  - a8
  c3d9  = c5d9  - a9
  c3d10 = c5d10 - a10
