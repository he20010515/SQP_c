!// ForQuill v1.01 Beta www.fcode.cn
!
!      ALGORITHM 733, COLLECTED ALGORITHMS FROM ACM.
!      TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 20, NO. 3, SEPTEMBER, 1994, PP. 262-281.
!      https://doi.org/10.1145/192115.192124
!
!
!      https://web.archive.org/web/20170106155705/http://permalink.gmane.org/gmane.comp.python.scientific.devel/6725
!      ------
!      From: Deborah Cotton <cotton@hq.acm.org>
!      Date: Fri, 14 Sep 2007 12:35:55 -0500
!      Subject: RE: Algorithm License requested
!      To: Alan Isaac
!
!      Prof. Issac,
!
!      In that case, then because the author consents to [the ACM] releasing
!      the code currently archived at http://www.netlib.org/toms/733 under the
!      BSD license, the ACM hereby releases this code under the BSD license.
!
!      Regards,
!
!      Deborah Cotton, Copyright & Permissions
!      ACM Publications
!      2 Penn Plaza, Suite 701**
!      New York, NY 10121-0701
!      permissions@acm.org
!      212.869.7440 ext. 652
!      Fax. 212.869.0481
!      ------
!

!***********************************************************************
!                              optimizer                               *
!***********************************************************************

Subroutine slsqp(m, meq, la, n, x, xl, xu, f, c, g, a, acc, iter, mode, w, l_w, jw, l_jw, alpha, f0, gs, h1, h2, h3, h4, t, t0, tol, iexact, incons, ireset, itermx, line, n1, n2, n3)

    !   SLSQP       S EQUENTIAL  L EAST  SQ UARES  P ROGRAMMING
    !            TO SOLVE GENERAL NONLINEAR OPTIMIZATION PROBLEMS
    
    !***********************************************************************
    !*                                                                     *
    !*                                                                     *
    !*            A NONLINEAR PROGRAMMING METHOD WITH                      *
    !*            QUADRATIC  PROGRAMMING  SUBPROBLEMS                      *
    !*                                                                     *
    !*                                                                     *
    !*  THIS SUBROUTINE SOLVES THE GENERAL NONLINEAR PROGRAMMING PROBLEM   *
    !*                                                                     *
    !*            MINIMIZE    F(X)                                         *
    !*                                                                     *
    !*            SUBJECT TO  C (X) .EQ. 0  ,  J = 1,...,MEQ               *
    !*                         J                                           *
    !*                                                                     *
    !*                        C (X) .GE. 0  ,  J = MEQ+1,...,M             *
    !*                         J                                           *
    !*                                                                     *
    !*                        XL .LE. X .LE. XU , I = 1,...,N.             *
    !*                          I      I       I                           *
    !*                                                                     *
    !*  THE ALGORITHM IMPLEMENTS THE METHOD OF HAN AND POWELL              *
    !*  WITH BFGS-UPDATE OF THE B-MATRIX AND L1-TEST FUNCTION              *
    !*  WITHIN THE STEPLENGTH ALGORITHM.                                   *
    !*                                                                     *
    !*    PARAMETER DESCRIPTION:                                           *
    !*    ( * MEANS THIS PARAMETER WILL BE CHANGED DURING CALCULATION )    *
    !*                                                                     *
    !*    M              IS THE TOTAL NUMBER OF CONSTRAINTS, M .GE. 0      *
    !*    MEQ            IS THE NUMBER OF EQUALITY CONSTRAINTS, MEQ .GE. 0 *
    !*    LA             SEE A, LA .GE. MAX(M,1)                           *
    !*    N              IS THE NUMBER OF VARIBLES, N .GE. 1               *
    !*  * X()            X() STORES THE CURRENT ITERATE OF THE N VECTOR X  *
    !*                   ON ENTRY X() MUST BE INITIALIZED. ON EXIT X()     *
    !*                   STORES THE SOLUTION VECTOR X IF MODE = 0.         *
    !*    XL()           XL() STORES AN N VECTOR OF LOWER BOUNDS XL TO X.  *
    !*                   ELEMENTS MAY BE NAN TO INDICATE NO LOWER BOUND.   *
    !*    XU()           XU() STORES AN N VECTOR OF UPPER BOUNDS XU TO X.  *
    !*                   ELEMENTS MAY BE NAN TO INDICATE NO UPPER BOUND.   *
    !*    F              IS THE VALUE OF THE OBJECTIVE FUNCTION.           *
    !*    C()            C() STORES THE M VECTOR C OF CONSTRAINTS,         *
    !*                   EQUALITY CONSTRAINTS (IF ANY) FIRST.              *
    !*                   DIMENSION OF C MUST BE GREATER OR EQUAL LA,       *
    !*                   which must be GREATER OR EQUAL MAX(1,M).          *
    !*    G()            G() STORES THE N VECTOR G OF PARTIALS OF THE      *
    !*                   OBJECTIVE FUNCTION; DIMENSION OF G MUST BE        *
    !*                   GREATER OR EQUAL N+1.                             *
    !*    A(),LA,M,N     THE LA BY N + 1 ARRAY A() STORES                  *
    !*                   THE M BY N MATRIX A OF CONSTRAINT NORMALS.        *
    !*                   A() HAS FIRST DIMENSIONING PARAMETER LA,          *
    !*                   WHICH MUST BE GREATER OR EQUAL MAX(1,M).          *
    !*    F,C,G,A        MUST ALL BE SET BY THE USER BEFORE EACH CALL.     *
    !*  * ACC            ABS(ACC) CONTROLS THE FINAL ACCURACY.             *
    !*                   IF ACC .LT. ZERO AN EXACT LINESEARCH IS PERFORMED,*
    !*                   OTHERWISE AN ARMIJO-TYPE LINESEARCH IS USED.      *
    !*  * ITER           PRESCRIBES THE MAXIMUM NUMBER OF ITERATIONS.      *
    !*                   ON EXIT ITER INDICATES THE NUMBER OF ITERATIONS.  *
    !*  * MODE           MODE CONTROLS CALCULATION:                        *
    !*                   REVERSE COMMUNICATION IS USED IN THE SENSE THAT   *
    !*                   THE PROGRAM IS INITIALIZED BY MODE = 0; THEN IT IS*
    !*                   TO BE CALLED REPEATEDLY BY THE USER UNTIL A RETURN*
    !*                   WITH MODE .NE. IABS(1) TAKES PLACE.               *
    !*                   IF MODE = -1 GRADIENTS HAVE TO BE CALCULATED,     *
    !*                   WHILE WITH MODE = 1 FUNCTIONS HAVE TO BE CALCULATED
    !*                   MODE MUST NOT BE CHANGED BETWEEN SUBSEQUENT CALLS *
    !*                   OF SQP.                                           *
    !*                   EVALUATION MODES:                                 *
    !*        MODE = -1: GRADIENT EVALUATION, (G&A)                        *
    !*                0: ON ENTRY: INITIALIZATION, (F,G,C&A)               *
    !*                   ON EXIT : REQUIRED ACCURACY FOR SOLUTION OBTAINED *
    !*                1: FUNCTION EVALUATION, (F&C)                        *
    !*                                                                     *
    !*                   FAILURE MODES:                                    *
    !*                2: NUMBER OF EQUALITY CONSTRAINTS LARGER THAN N      *
    !*                3: MORE THAN 3*N ITERATIONS IN LSQ SUBPROBLEM        *
    !*                4: INEQUALITY CONSTRAINTS INCOMPATIBLE               *
    !*                5: SINGULAR MATRIX E IN LSQ SUBPROBLEM               *
    !*                6: SINGULAR MATRIX C IN LSQ SUBPROBLEM               *
    !*                7: RANK-DEFICIENT EQUALITY CONSTRAINT SUBPROBLEM HFTI*
    !*                8: POSITIVE DIRECTIONAL DERIVATIVE FOR LINESEARCH    *
    !*                9: MORE THAN ITER ITERATIONS IN SQP                  *
    !*             >=10: WORKING SPACE W OR JW TOO SMALL,                  *
    !*                   W SHOULD BE ENLARGED TO L_W=MODE/1000             *
    !*                   JW SHOULD BE ENLARGED TO L_JW=MODE-1000*L_W       *
    !*  * W(), L_W       W() IS A ONE DIMENSIONAL WORKING SPACE,           *
    !*                   THE LENGTH L_W OF WHICH SHOULD BE AT LEAST        *
    !*                   (3*N1+M)*(N1+1)                        for LSQ    *
    !*                  +(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ         for LSI    *
    !*                  +(N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1       for LSEI   *
    !*                  + N1*N/2 + 2*M + 3*N + 3*N1 + 1         for SLSQPB *
    !*                   with MINEQ = M - MEQ + 2*N1  &  N1 = N+1          *
    !*        NOTICE:    FOR PROPER DIMENSIONING OF W IT IS RECOMMENDED TO *
    !*                   COPY THE FOLLOWING STATEMENTS INTO THE HEAD OF    *
    !*                   THE CALLING PROGRAM (AND REMOVE THE COMMENT C)    *
    !#######################################################################
    !     INTEGER LEN_W, LEN_JW, M, N, N1, MEQ, MINEQ
    !     PARAMETER (M=... , MEQ=... , N=...  )
    !     PARAMETER (N1= N+1, MINEQ= M-MEQ+N1+N1)
    !     PARAMETER (LEN_W=
    !    $           (3*N1+M)*(N1+1)
    !    $          +(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ
    !    $          +(N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1
    !    $          +(N+1)*N/2 + 2*M + 3*N + 3*N1 + 1,
    !    $           LEN_JW=MINEQ)
    !     DOUBLE PRECISION W(LEN_W)
    !     INTEGER          JW(LEN_JW)
    !#######################################################################
    !*                   THE FIRST M+N+N*N1/2 ELEMENTS OF W MUST NOT BE    *
    !*                   CHANGED BETWEEN SUBSEQUENT CALLS OF SLSQP.        *
    !*                   ON RETURN W(1) ... W(M) CONTAIN THE MULTIPLIERS   *
    !*                   ASSOCIATED WITH THE GENERAL CONSTRAINTS, WHILE    *
    !*                   W(M+1) ... W(M+N(N+1)/2) STORE THE CHOLESKY FACTOR*
    !*                   L*D*L(T) OF THE APPROXIMATE HESSIAN OF THE        *
    !*                   LAGRANGIAN COLUMNWISE DENSE AS LOWER TRIANGULAR   *
    !*                   UNIT MATRIX L WITH D IN ITS 'DIAGONAL' and        *
    !*                   W(M+N(N+1)/2+N+2 ... W(M+N(N+1)/2+N+2+M+2N)       *
    !*                   CONTAIN THE MULTIPLIERS ASSOCIATED WITH ALL       *
    !*                   ALL CONSTRAINTS OF THE QUADRATIC PROGRAM FINDING  *
    !*                   THE SEARCH DIRECTION TO THE SOLUTION X*           *
    !*  * JW(), L_JW     JW() IS A ONE DIMENSIONAL INTEGER WORKING SPACE   *
    !*                   THE LENGTH L_JW OF WHICH SHOULD BE AT LEAST       *
    !*                   MINEQ                                             *
    !*                   with MINEQ = M - MEQ + 2*N1  &  N1 = N+1          *
    !*                                                                     *
    !*  THE USER HAS TO PROVIDE THE FOLLOWING SUBROUTINES:                 *
    !*     LDL(N,A,Z,SIG,W) :   UPDATE OF THE LDL'-FACTORIZATION.          *
    !*     LINMIN(A,B,F,TOL) :  LINESEARCH ALGORITHM IF EXACT = 1          *
    !*     LSQ(M,MEQ,LA,N,NC,C,D,A,B,XL,XU,X,LAMBDA,W,....) :              *
    !*                                                                     *
    !*        SOLUTION OF THE QUADRATIC PROGRAM                            *
    !*                QPSOL IS RECOMMENDED:                                *
    !*     PE GILL, W MURRAY, MA SAUNDERS, MH WRIGHT:                      *
    !*     USER'S GUIDE FOR SOL/QPSOL:                                     *
    !*     A FORTRAN PACKAGE FOR QUADRATIC PROGRAMMING,                    *
    !*     TECHNICAL REPORT SOL 83-7, JULY 1983                            *
    !*     DEPARTMENT OF OPERATIONS RESEARCH, STANFORD UNIVERSITY          *
    !*     STANFORD, CA 94305                                              *
    !*     QPSOL IS THE MOST ROBUST AND EFFICIENT QP-SOLVER                *
    !*     AS IT ALLOWS WARM STARTS WITH PROPER WORKING SETS               *
    !*                                                                     *
    !*     IF IT IS NOT AVAILABLE USE LSEI, A CONSTRAINT LINEAR LEAST      *
    !*     SQUARES SOLVER IMPLEMENTED USING THE SOFTWARE HFTI, LDP, NNLS   *
    !*     FROM C.L. LAWSON, R.J.HANSON: SOLVING LEAST SQUARES PROBLEMS,   *
    !*     PRENTICE HALL, ENGLEWOOD CLIFFS, 1974.                          *
    !*     LSEI COMES WITH THIS PACKAGE, together with all necessary SR's. *
    !*                                                                     *
    !*     TOGETHER WITH A COUPLE OF SUBROUTINES FROM BLAS LEVEL 1         *
    !*                                                                     *
    !*     SQP IS HEAD SUBROUTINE FOR BODY SUBROUTINE SQPBDY               *
    !*     IN WHICH THE ALGORITHM HAS BEEN IMPLEMENTED.                    *
    !*                                                                     *
    !*  IMPLEMENTED BY: DIETER KRAFT, DFVLR OBERPFAFFENHOFEN               *
    !*  as described in Dieter Kraft: A Software Package for               *
    !*                                Sequential Quadratic Programming     *
    !*                                DFVLR-FB 88-28, 1988                 *
    !*  which should be referenced if the user publishes results of SLSQP  *
    !*                                                                     *
    !*  DATE:           APRIL - OCTOBER, 1981.                             *
    !*  STATUS:         DECEMBER, 31-ST, 1984.                             *
    !*  STATUS:         MARCH   , 21-ST, 1987, REVISED TO FORTRAN 77       *
    !*  STATUS:         MARCH   , 20-th, 1989, REVISED TO MS-FORTRAN       *
    !*  STATUS:         APRIL   , 14-th, 1989, HESSE   in-line coded       *
    !*  STATUS:         FEBRUARY, 28-th, 1991, FORTRAN/2 Version 1.04      *
    !*                                         accepts Statement Functions *
    !*  STATUS:         MARCH   ,  1-st, 1991, tested with SALFORD         *
    !*                                         FTN77/386 COMPILER VERS 2.40*
    !*                                         in protected mode           *
    !*                                                                     *
    !***********************************************************************
    !*                                                                     *
    !*  Copyright 1991: Dieter Kraft, FHM                                  *
    !*                                                                     *
    !***********************************************************************
    
      Integer il, im, ir, is, iter, iu, iv, iw, ix, l_w, l_jw, jw(l_jw), la, m, meq, mineq, mode, n
    
      Double Precision acc, a(la, n+1), c(la), f, g(n+1), x(n), xl(n), xu(n), w(l_w)
    
      Integer iexact, incons, ireset, itermx, line, n1, n2, n3
    
      Double Precision alpha, f0, gs, h1, h2, h3, h4, t, t0, tol
    
    !     dim(W) =         N1*(N1+1) + MEQ*(N1+1) + MINEQ*(N1+1)  for LSQ
    !                    +(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ          for LSI
    !                    +(N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1        for LSEI
    !                    + N1*N/2 + 2*M + 3*N +3*N1 + 1           for SLSQPB
    !                      with MINEQ = M - MEQ + 2*N1  &  N1 = N+1
    
    !   CHECK LENGTH OF WORKING ARRAYS
    
      n1 = n + 1
      mineq = m - meq + n1 + n1
      il = (3*n1+m)*(n1+1) + (n1-meq+1)*(mineq+2) + 2*mineq + (n1+mineq)*(n1-meq) + 2*meq + n1*n/2 + 2*m + 3*n + 4*n1 + 1
      im = max(mineq, n1-meq)
      If (l_w<il .Or. l_jw<im) Then
        mode = 1000*max(10, il)
        mode = mode + max(10, im)
        Return
      End If
    
    !   PREPARE DATA FOR CALLING SQPBDY  -  INITIAL ADDRESSES IN W
    
      im = 1
      il = im + max(1, m)
      il = im + la
      ix = il + n1*n/2 + 1
      ir = ix + n
      is = ir + n + n + max(1, m)
      is = ir + n + n + la
      iu = is + n1
      iv = iu + n1
      iw = iv + n1
    
      Call slsqpb(m, meq, la, n, x, xl, xu, f, c, g, a, acc, iter, mode, w(ir), w(il), w(ix), w(im), w(is), w(iu), w(iv), w(iw), jw, alpha, f0, gs, h1, h2, h3, h4, t, t0, tol, iexact, incons, ireset, itermx, line, n1, n2, n3)
    
    End Subroutine slsqp
    
    Subroutine slsqpb(m, meq, la, n, x, xl, xu, f, c, g, a, acc, iter, mode, r, l, x0, mu, s, u, v, w, iw, alpha, f0, gs, h1, h2, h3, h4, t, t0, tol, iexact, incons, ireset, itermx, line, n1, n2, n3)
    
    !   NONLINEAR PROGRAMMING BY SOLVING SEQUENTIALLY QUADRATIC PROGRAMS
    
    !        -  L1 - LINE SEARCH,  POSITIVE DEFINITE  BFGS UPDATE  -
    
    !                      BODY SUBROUTINE FOR SLSQP
    
      Integer iw(*), i, iexact, incons, ireset, iter, itermx, k, j, la, line, m, meq, mode, n, n1, n2, n3
      Logical badlin
    
      Double Precision a(la, n+1), c(la), g(n+1), l((n+1)*(n+2)/2), mu(la), r(m+n+n+2), s(n+1), u(n+1), v(n+1), w(*), x(n), xl(n), xu(n), x0(n), ddot_sl, dnrm2_, linmin, acc, alfmin, alpha, f, f0, gs, h1, h2, h3, h4, hun, one, t, t0, ten, tol, two, zero
    
    !     dim(W) =         N1*(N1+1) + MEQ*(N1+1) + MINEQ*(N1+1)  for LSQ
    !                     +(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ
    !                     +(N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1       for LSEI
    !                      with MINEQ = M - MEQ + 2*N1  &  N1 = N+1
    
      Data zero/0.0D0/, one/1.0D0/, alfmin/1.0D-1/, hun/1.0D+2/, ten/1.0D+1/, two/2.0D0/
    
    !     The badlin flag keeps track whether the SQP problem on the current
    !     iteration was inconsistent or not.
      badlin = .False.
    
      If (mode) 260, 100, 220
    
      100 itermx = iter
      If (acc>=zero) Then
        iexact = 0
      Else
        iexact = 1
      End If
      acc = abs(acc)
      tol = ten*acc
      iter = 0
      ireset = 0
      n1 = n + 1
      n2 = n1*n/2
      n3 = n2 + 1
      s(1) = zero
      mu(1) = zero
      Call dcopy_(n, s(1), 0, s, 1)
      Call dcopy_(m, mu(1), 0, mu, 1)
    
    !   RESET BFGS MATRIX
    
      110 ireset = ireset + 1
      If (ireset>5) Goto 255
      l(1) = zero
      Call dcopy_(n2, l(1), 0, l, 1)
      j = 1
      Do i = 1, n
        l(j) = one
        j = j + n1 - i
      End Do
    
    !   MAIN ITERATION : SEARCH DIRECTION, STEPLENGTH, LDL'-UPDATE
    
      130 iter = iter + 1
      mode = 9
      If (iter>itermx) Goto 330
    
    !   SEARCH DIRECTION AS SOLUTION OF QP - SUBPROBLEM
    
      Call dcopy_(n, xl, 1, u, 1)
      Call dcopy_(n, xu, 1, v, 1)
      Call daxpy_sl(n, -one, x, 1, u, 1)
      Call daxpy_sl(n, -one, x, 1, v, 1)
      h4 = one
      Call lsq(m, meq, n, n3, la, l, g, a, c, u, v, s, r, w, iw, mode)
    
    !   AUGMENTED PROBLEM FOR INCONSISTENT LINEARIZATION
    !
    !   If it turns out that the original SQP problem is inconsistent,
    !   disallow termination with convergence on this iteration,
    !   even if the augmented problem was solved.
    
      badlin = .False.
      If (mode==6) Then
        If (n==meq) Then
          mode = 4
        End If
      End If
      If (mode==4) Then
        badlin = .True.
        Do j = 1, m
          If (j<=meq) Then
            a(j, n1) = -c(j)
          Else
            a(j, n1) = max(-c(j), zero)
          End If
        End Do
        s(1) = zero
        Call dcopy_(n, s(1), 0, s, 1)
        h3 = zero
        g(n1) = zero
        l(n3) = hun
        s(n1) = one
        u(n1) = zero
        v(n1) = one
        incons = 0
        150 Call lsq(m, meq, n1, n3, la, l, g, a, c, u, v, s, r, w, iw, mode)
        h4 = one - s(n1)
        If (mode==4) Then
          l(n3) = ten*l(n3)
          incons = incons + 1
          If (incons>5) Goto 330
          Goto 150
        Else If (mode/=1) Then
          Goto 330
        End If
      Else If (mode/=1) Then
        Goto 330
      End If
    
    !   UPDATE MULTIPLIERS FOR L1-TEST
    
      Do i = 1, n
        v(i) = g(i) - ddot_sl(m, a(1,i), 1, r, 1)
      End Do
      f0 = f
      Call dcopy_(n, x, 1, x0, 1)
      gs = ddot_sl(n, g, 1, s, 1)
      h1 = abs(gs)
      h2 = zero
      Do j = 1, m
        If (j<=meq) Then
          h3 = c(j)
        Else
          h3 = zero
        End If
        h2 = h2 + max(-c(j), h3)
        h3 = abs(r(j))
        mu(j) = max(h3, (mu(j)+h3)/two)
        h1 = h1 + h3*abs(c(j))
      End Do
    
    !   CHECK CONVERGENCE
    
      mode = 0
      If (h1<acc .And. h2<acc .And. .Not. badlin .And. f==f) Goto 330
      h1 = zero
      Do j = 1, m
        If (j<=meq) Then
          h3 = c(j)
        Else
          h3 = zero
        End If
        h1 = h1 + mu(j)*max(-c(j), h3)
      End Do
      t0 = f + h1
      h3 = gs - h1*h4
      mode = 8
      If (h3>=zero) Goto 110
    
    !   LINE SEARCH WITH AN L1-TESTFUNCTION
    
      line = 0
      alpha = one
      If (iexact==1) Goto 210
    
    !   INEXACT LINESEARCH
    
      190 line = line + 1
      h3 = alpha*h3
      Call dscal_sl(n, alpha, s, 1)
      Call dcopy_(n, x0, 1, x, 1)
      Call daxpy_sl(n, one, s, 1, x, 1)
      mode = 1
      Goto 330
      200 If (h1<=h3/ten .Or. line>10) Goto 240
      alpha = max(h3/(two*(h3-h1)), alfmin)
      Goto 190
    
    !   EXACT LINESEARCH
    
      210 If (line/=3) Then
        alpha = linmin(line, alfmin, one, t, tol)
        Call dcopy_(n, x0, 1, x, 1)
        Call daxpy_sl(n, alpha, s, 1, x, 1)
        mode = 1
        Goto 330
      End If
      Call dscal_sl(n, alpha, s, 1)
      Goto 240
    
    !   CALL FUNCTIONS AT CURRENT X
    
      220 t = f
      Do j = 1, m
        If (j<=meq) Then
          h1 = c(j)
        Else
          h1 = zero
        End If
        t = t + mu(j)*max(-c(j), h1)
      End Do
      h1 = t - t0
      Goto (200, 210) iexact + 1
    
    !   CHECK CONVERGENCE
    
      240 h3 = zero
      Do j = 1, m
        If (j<=meq) Then
          h1 = c(j)
        Else
          h1 = zero
        End If
        h3 = h3 + max(-c(j), h1)
      End Do
      If ((abs(f-f0)<acc .Or. dnrm2_(n,s,1)<acc) .And. h3<acc .And. .Not. badlin .And. f==f) Then
        mode = 0
      Else
        mode = -1
      End If
      Goto 330
    
    !   CHECK relaxed CONVERGENCE in case of positive directional derivative
    
      255 Continue
      h3 = zero
      Do j = 1, m
        If (j<=meq) Then
          h1 = c(j)
        Else
          h1 = zero
        End If
        h3 = h3 + max(-c(j), h1)
      End Do
      If ((abs(f-f0)<tol .Or. dnrm2_(n,s,1)<tol) .And. h3<tol .And. .Not. badlin .And. f==f) Then
        mode = 0
      Else
        mode = 8
      End If
      Goto 330
    
    !   CALL JACOBIAN AT CURRENT X
    
    !   UPDATE CHOLESKY-FACTORS OF HESSIAN MATRIX BY MODIFIED BFGS FORMULA
    
      260 Do i = 1, n
        u(i) = g(i) - ddot_sl(m, a(1,i), 1, r, 1) - v(i)
      End Do
    
    !   L'*S
    
      k = 0
      Do i = 1, n
        h1 = zero
        k = k + 1
        Do j = i + 1, n
          k = k + 1
          h1 = h1 + l(k)*s(j)
        End Do
        v(i) = s(i) + h1
      End Do
    
    !   D*L'*S
    
      k = 1
      Do i = 1, n
        v(i) = l(k)*v(i)
        k = k + n1 - i
      End Do
    
    !   L*D*L'*S
    
      Do i = n, 1, -1
        h1 = zero
        k = i
        Do j = 1, i - 1
          h1 = h1 + l(k)*v(j)
          k = k + n - j
        End Do
        v(i) = v(i) + h1
      End Do
    
      h1 = ddot_sl(n, s, 1, u, 1)
      h2 = ddot_sl(n, s, 1, v, 1)
      h3 = 0.2D0*h2
      If (h1<h3) Then
        h4 = (h2-h3)/(h2-h1)
        h1 = h3
        Call dscal_sl(n, h4, u, 1)
        Call daxpy_sl(n, one-h4, v, 1, u, 1)
      End If
      If (h1==0 .Or. h2==0) Then
    !         Singular update: reset hessian.
        Goto 110
      End If
      Call ldl(n, l, u, +one/h1, v)
      Call ldl(n, l, v, -one/h2, u)
    
    !   END OF MAIN ITERATION
    
      Goto 130
    
    !   END OF SLSQPB
    
    330 End Subroutine slsqpb
    
    
    Subroutine lsq(m, meq, n, nl, la, l, g, a, b, xl, xu, x, y, w, jw, mode)
    
    !   MINIMIZE with respect to X
    
    !             ||E*X - F||
    !                                      1/2  T
    !   WITH UPPER TRIANGULAR MATRIX E = +D   *L ,
    
    !                                      -1/2  -1
    !                     AND VECTOR F = -D    *L  *G,
    
    !  WHERE THE UNIT LOWER TRIDIANGULAR MATRIX L IS STORED COLUMNWISE
    !  DENSE IN THE N*(N+1)/2 ARRAY L WITH VECTOR D STORED IN ITS
    ! 'DIAGONAL' THUS SUBSTITUTING THE ONE-ELEMENTS OF L
    
    !   SUBJECT TO
    
    !             A(J)*X - B(J) = 0 ,         J=1,...,MEQ,
    !             A(J)*X - B(J) >=0,          J=MEQ+1,...,M,
    !             XL(I) <= X(I) <= XU(I),     I=1,...,N,
    !     ON ENTRY, THE USER HAS TO PROVIDE THE ARRAYS L, G, A, B, XL, XU.
    !     WITH DIMENSIONS: L(N*(N+1)/2), G(N), A(LA,N), B(M), XL(N), XU(N)
    !     THE WORKING ARRAY W MUST HAVE AT LEAST THE FOLLOWING DIMENSION:
    !     DIM(W) =        (3*N+M)*(N+1)                        for LSQ
    !                    +(N-MEQ+1)*(MINEQ+2) + 2*MINEQ        for LSI
    !                    +(N+MINEQ)*(N-MEQ) + 2*MEQ + N        for LSEI
    !                      with MINEQ = M - MEQ + 2*N
    !     ON RETURN, NO ARRAY WILL BE CHANGED BY THE SUBROUTINE.
    !     X     STORES THE N-DIMENSIONAL SOLUTION VECTOR
    !     Y     STORES THE VECTOR OF LAGRANGE MULTIPLIERS OF DIMENSION
    !           M+N+N (CONSTRAINTS+LOWER+UPPER BOUNDS)
    !     MODE  IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANINGS:
    !          MODE=1: SUCCESSFUL COMPUTATION
    !               2: ERROR RETURN BECAUSE OF WRONG DIMENSIONS (N<1)
    !               3: ITERATION COUNT EXCEEDED BY NNLS
    !               4: INEQUALITY CONSTRAINTS INCOMPATIBLE
    !               5: MATRIX E IS NOT OF FULL RANK
    !               6: MATRIX C IS NOT OF FULL RANK
    !               7: RANK DEFECT IN HFTI
    
    !     coded            Dieter Kraft, april 1987
    !     revised                        march 1989
    
      Double Precision l, g, a, b, w, xl, xu, x, y, diag, zero, one, ddot_sl, xnorm
    
      Integer jw(*), i, ic, id, ie, if, ig, ih, il, im, ip, iu, iw, i1, i2, i3, i4, la, m, meq, mineq, mode, m1, n, nl, n1, n2, n3, nancnt, j
    
      Dimension a(la, n), b(la), g(n), l(nl), w(*), x(n), xl(n), xu(n), y(m+n+n)
    
      Data zero/0.0D0/, one/1.0D0/
    
      n1 = n + 1
      mineq = m - meq
      m1 = mineq + n + n
    
    !  determine whether to solve problem
    !  with inconsistent linerarization (n2=1)
    !  or not (n2=0)
    
      n2 = n1*n/2 + 1
      If (n2==nl) Then
        n2 = 0
      Else
        n2 = 1
      End If
      n3 = n - n2
    
    !  RECOVER MATRIX E AND VECTOR F FROM L AND G
    
      i2 = 1
      i3 = 1
      i4 = 1
      ie = 1
      if = n*n + 1
      Do i = 1, n3
        i1 = n1 - i
        diag = sqrt(l(i2))
        w(i3) = zero
        Call dcopy_(i1, w(i3), 0, w(i3), 1)
        Call dcopy_(i1-n2, l(i2), 1, w(i3), n)
        Call dscal_sl(i1-n2, diag, w(i3), n)
        w(i3) = diag
        w(if-1+i) = (g(i)-ddot_sl(i-1,w(i4),1,w(if),1))/diag
        i2 = i2 + i1 - n2
        i3 = i3 + n1
        i4 = i4 + n
      End Do
      If (n2==1) Then
        w(i3) = l(nl)
        w(i4) = zero
        Call dcopy_(n3, w(i4), 0, w(i4), 1)
        w(if-1+n) = zero
      End If
      Call dscal_sl(n, -one, w(if), 1)
    
      ic = if + n
      id = ic + meq*n
    
      If (meq>0) Then
    
    !  RECOVER MATRIX C FROM UPPER PART OF A
    
        Do i = 1, meq
          Call dcopy_(n, a(i,1), la, w(ic-1+i), meq)
        End Do
    
    !  RECOVER VECTOR D FROM UPPER PART OF B
    
        Call dcopy_(meq, b(1), 1, w(id), 1)
        Call dscal_sl(meq, -one, w(id), 1)
    
      End If
    
      ig = id + meq
    
    !  RECOVER MATRIX G FROM LOWER PART OF A
    !  The matrix G(mineq+2*n,m1) is stored at w(ig)
    !  Not all rows will be filled if some of the upper/lower
    !  bounds are unbounded.
    
      If (mineq>0) Then
    
        Do i = 1, mineq
          Call dcopy_(n, a(meq+i,1), la, w(ig-1+i), m1)
        End Do
    
      End If
    
      ih = ig + m1*n
      iw = ih + mineq + 2*n
    
      If (mineq>0) Then
    
    !  RECOVER H FROM LOWER PART OF B
    !  The vector H(mineq+2*n) is stored at w(ih)
    
        Call dcopy_(mineq, b(meq+1), 1, w(ih), 1)
        Call dscal_sl(mineq, -one, w(ih), 1)
    
      End If
    
    !  AUGMENT MATRIX G BY +I AND -I, AND,
    !  AUGMENT VECTOR H BY XL AND XU
    !  NaN value indicates no bound
    
      ip = ig + mineq
      il = ih + mineq
      nancnt = 0
    
      Do i = 1, n
        If (xl(i)==xl(i)) Then
          w(il) = xl(i)
          Do j = 1, n
            w(ip+m1*(j-1)) = 0
          End Do
          w(ip+m1*(i-1)) = 1
          ip = ip + 1
          il = il + 1
        Else
          nancnt = nancnt + 1
        End If
      End Do
    
      Do i = 1, n
        If (xu(i)==xu(i)) Then
          w(il) = -xu(i)
          Do j = 1, n
            w(ip+m1*(j-1)) = 0
          End Do
          w(ip+m1*(i-1)) = -1
          ip = ip + 1
          il = il + 1
        Else
          nancnt = nancnt + 1
        End If
      End Do
    
      Call lsei(w(ic), w(id), w(ie), w(if), w(ig), w(ih), max(1,meq), meq, n, n, m1, m1-nancnt, n, x, xnorm, w(iw), jw, mode)
    
      If (mode==1) Then
    
    !   restore Lagrange multipliers (only for user-defined variables)
    
        Call dcopy_(m, w(iw), 1, y(1), 1)
    
    !   set rest of the multipliers to nan (they are not used)
    
        If (n3>0) Then
          y(m+1) = 0
          y(m+1) = 0/y(m+1)
          Do i = m + 2, m + n3 + n3
            y(i) = y(m+1)
          End Do
        End If
    
      End If
      Call bound(n, x, xl, xu)
    
    !   END OF SUBROUTINE LSQ
    
    End Subroutine lsq
    
    
    Subroutine lsei(c, d, e, f, g, h, lc, mc, le, me, lg, mg, n, x, xnrm, w, jw, mode)
    
    !     FOR MODE=1, THE SUBROUTINE RETURNS THE SOLUTION X OF
    !     EQUALITY & INEQUALITY CONSTRAINED LEAST SQUARES PROBLEM LSEI :
    
    !                MIN ||E*X - F||
    !                 X
    
    !                S.T.  C*X  = D,
    !                      G*X >= H.
    
    !     USING QR DECOMPOSITION & ORTHOGONAL BASIS OF NULLSPACE OF C
    !     CHAPTER 23.6 OF LAWSON & HANSON: SOLVING LEAST SQUARES PROBLEMS.
    
    !     THE FOLLOWING DIMENSIONS OF THE ARRAYS DEFINING THE PROBLEM
    !     ARE NECESSARY
    !     DIM(E) :   FORMAL (LE,N),    ACTUAL (ME,N)
    !     DIM(F) :   FORMAL (LE  ),    ACTUAL (ME  )
    !     DIM(C) :   FORMAL (LC,N),    ACTUAL (MC,N)
    !     DIM(D) :   FORMAL (LC  ),    ACTUAL (MC  )
    !     DIM(G) :   FORMAL (LG,N),    ACTUAL (MG,N)
    !     DIM(H) :   FORMAL (LG  ),    ACTUAL (MG  )
    !     DIM(X) :   FORMAL (N   ),    ACTUAL (N   )
    !     DIM(W) :   2*MC+ME+(ME+MG)*(N-MC)  for LSEI
    !              +(N-MC+1)*(MG+2)+2*MG     for LSI
    !     DIM(JW):   MAX(MG,L)
    !     ON ENTRY, THE USER HAS TO PROVIDE THE ARRAYS C, D, E, F, G, AND H.
    !     ON RETURN, ALL ARRAYS WILL BE CHANGED BY THE SUBROUTINE.
    !     X     STORES THE SOLUTION VECTOR
    !     XNORM STORES THE RESIDUUM OF THE SOLUTION IN EUCLIDIAN NORM
    !     W     STORES THE VECTOR OF LAGRANGE MULTIPLIERS IN ITS FIRST
    !           MC+MG ELEMENTS
    !     MODE  IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANINGS:
    !          MODE=1: SUCCESSFUL COMPUTATION
    !               2: ERROR RETURN BECAUSE OF WRONG DIMENSIONS (N<1)
    !               3: ITERATION COUNT EXCEEDED BY NNLS
    !               4: INEQUALITY CONSTRAINTS INCOMPATIBLE
    !               5: MATRIX E IS NOT OF FULL RANK
    !               6: MATRIX C IS NOT OF FULL RANK
    !               7: RANK DEFECT IN HFTI
    
    !     18.5.1981, DIETER KRAFT, DFVLR OBERPFAFFENHOFEN
    !     20.3.1987, DIETER KRAFT, DFVLR OBERPFAFFENHOFEN
    
      Integer jw(*), i, ie, if, ig, iw, j, k, krank, l, lc, le, lg, mc, mc1, me, mg, mode, n
      Double Precision c(lc, n), e(le, n), g(lg, n), d(lc), f(le), h(lg), x(n), w(*), t, ddot_sl, xnrm, dnrm2_, epmach, zero
      Data epmach/2.22D-16/, zero/0.0D+00/
    
      mode = 2
      If (mc>n) Goto 75
      l = n - mc
      mc1 = mc + 1
      iw = (l+1)*(mg+2) + 2*mg + mc
      ie = iw + mc + 1
      if = ie + me*l
      ig = if + me
    
    !  TRIANGULARIZE C AND APPLY FACTORS TO E AND G
    
      Do i = 1, mc
        j = min(i+1, lc)
        Call h12(1, i, i+1, n, c(i,1), lc, w(iw+i), c(j,1), lc, 1, mc-i)
        Call h12(2, i, i+1, n, c(i,1), lc, w(iw+i), e, le, 1, me)
        Call h12(2, i, i+1, n, c(i,1), lc, w(iw+i), g, lg, 1, mg)
      End Do
    
    !  SOLVE C*X=D AND MODIFY F
    
      mode = 6
      Do i = 1, mc
        If (abs(c(i,i))<epmach) Goto 75
        x(i) = (d(i)-ddot_sl(i-1,c(i,1),lc,x,1))/c(i, i)
      End Do
      mode = 1
      w(mc1) = zero
      Call dcopy_(mg-mc, w(mc1), 0, w(mc1), 1)
    
      If (mc==n) Goto 50
    
      Do i = 1, me
        w(if-1+i) = f(i) - ddot_sl(mc, e(i,1), le, x, 1)
      End Do
    
    !  STORE TRANSFORMED E & G
    
      Do i = 1, me
        Call dcopy_(l, e(i,mc1), le, w(ie-1+i), me)
      End Do
      Do i = 1, mg
        Call dcopy_(l, g(i,mc1), lg, w(ig-1+i), mg)
      End Do
    
      If (mg>0) Goto 40
    
    !  SOLVE LS WITHOUT INEQUALITY CONSTRAINTS
    
      mode = 7
      k = max(le, n)
      t = sqrt(epmach)
      Call hfti(w(ie), me, me, l, w(if), k, 1, t, krank, xnrm, w, w(l+1), jw)
      Call dcopy_(l, w(if), 1, x(mc1), 1)
      If (krank/=l) Goto 75
      mode = 1
      Goto 50
    !  MODIFY H AND SOLVE INEQUALITY CONSTRAINED LS PROBLEM
    
      40 Do i = 1, mg
        h(i) = h(i) - ddot_sl(mc, g(i,1), lg, x, 1)
      End Do
      Call lsi(w(ie), w(if), w(ig), h, me, me, mg, mg, l, x(mc1), xnrm, w(mc1), jw, mode)
      If (mc==0) Goto 75
      t = dnrm2_(mc, x, 1)
      xnrm = sqrt(xnrm*xnrm+t*t)
      If (mode/=1) Goto 75
    
    !  SOLUTION OF ORIGINAL PROBLEM AND LAGRANGE MULTIPLIERS
    
      50 Do i = 1, me
        f(i) = ddot_sl(n, e(i,1), le, x, 1) - f(i)
      End Do
      Do i = 1, mc
        d(i) = ddot_sl(me, e(1,i), 1, f, 1) - ddot_sl(mg, g(1,i), 1, w(mc1), 1)
      End Do
    
      Do i = mc, 1, -1
        Call h12(2, i, i+1, n, c(i,1), lc, w(iw+i), x, 1, 1, 1)
      End Do
    
      Do i = mc, 1, -1
        j = min(i+1, lc)
        w(i) = (d(i)-ddot_sl(mc-i,c(j,i),1,w(j),1))/c(i, i)
      End Do
    
    !  END OF SUBROUTINE LSEI
    
    75 End Subroutine lsei
    
    
    Subroutine lsi(e, f, g, h, le, me, lg, mg, n, x, xnorm, w, jw, mode)
    
    !     FOR MODE=1, THE SUBROUTINE RETURNS THE SOLUTION X OF
    !     INEQUALITY CONSTRAINED LINEAR LEAST SQUARES PROBLEM:
    
    !                    MIN ||E*X-F||
    !                     X
    
    !                    S.T.  G*X >= H
    
    !     THE ALGORITHM IS BASED ON QR DECOMPOSITION AS DESCRIBED IN
    !     CHAPTER 23.5 OF LAWSON & HANSON: SOLVING LEAST SQUARES PROBLEMS
    
    !     THE FOLLOWING DIMENSIONS OF THE ARRAYS DEFINING THE PROBLEM
    !     ARE NECESSARY
    !     DIM(E) :   FORMAL (LE,N),    ACTUAL (ME,N)
    !     DIM(F) :   FORMAL (LE  ),    ACTUAL (ME  )
    !     DIM(G) :   FORMAL (LG,N),    ACTUAL (MG,N)
    !     DIM(H) :   FORMAL (LG  ),    ACTUAL (MG  )
    !     DIM(X) :   N
    !     DIM(W) :   (N+1)*(MG+2) + 2*MG
    !     DIM(JW):   LG
    !     ON ENTRY, THE USER HAS TO PROVIDE THE ARRAYS E, F, G, AND H.
    !     ON RETURN, ALL ARRAYS WILL BE CHANGED BY THE SUBROUTINE.
    !     X     STORES THE SOLUTION VECTOR
    !     XNORM STORES THE RESIDUUM OF THE SOLUTION IN EUCLIDIAN NORM
    !     W     STORES THE VECTOR OF LAGRANGE MULTIPLIERS IN ITS FIRST
    !           MG ELEMENTS
    !     MODE  IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANINGS:
    !          MODE=1: SUCCESSFUL COMPUTATION
    !               2: ERROR RETURN BECAUSE OF WRONG DIMENSIONS (N<1)
    !               3: ITERATION COUNT EXCEEDED BY NNLS
    !               4: INEQUALITY CONSTRAINTS INCOMPATIBLE
    !               5: MATRIX E IS NOT OF FULL RANK
    
    !     03.01.1980, DIETER KRAFT: CODED
    !     20.03.1987, DIETER KRAFT: REVISED TO FORTRAN 77
    
      Integer i, j, le, lg, me, mg, mode, n, jw(lg)
      Double Precision e(le, n), f(le), g(lg, n), h(lg), x(n), w(*), ddot_sl, xnorm, dnrm2_, epmach, t, one
      Data epmach/2.22D-16/, one/1.0D+00/
    
    !  QR-FACTORS OF E AND APPLICATION TO F
    
      Do i = 1, n
        j = min(i+1, n)
        Call h12(1, i, i+1, me, e(1,i), 1, t, e(1,j), 1, le, n-i)
        Call h12(2, i, i+1, me, e(1,i), 1, t, f, 1, 1, 1)
      End Do
    
    !  TRANSFORM G AND H TO GET LEAST DISTANCE PROBLEM
    
      mode = 5
      Do i = 1, mg
        Do j = 1, n
          If (.Not. (abs(e(j,j))>=epmach)) Goto 50
          g(i, j) = (g(i,j)-ddot_sl(j-1,g(i,1),lg,e(1,j),1))/e(j, j)
        End Do
        h(i) = h(i) - ddot_sl(n, g(i,1), lg, f, 1)
      End Do
    
    !  SOLVE LEAST DISTANCE PROBLEM
    
      Call ldp(g, lg, mg, n, h, x, xnorm, w, jw, mode)
      If (mode/=1) Goto 50
    
    !  SOLUTION OF ORIGINAL PROBLEM
    
      Call daxpy_sl(n, one, f, 1, x, 1)
      Do i = n, 1, -1
        j = min(i+1, n)
        x(i) = (x(i)-ddot_sl(n-i,e(i,j),le,x(j),1))/e(i, i)
      End Do
      j = min(n+1, me)
      t = dnrm2_(me-n, f(j), 1)
      xnorm = sqrt(xnorm*xnorm+t*t)
    
    !  END OF SUBROUTINE LSI
    
    50 End Subroutine lsi
    
    Subroutine ldp(g, mg, m, n, h, x, xnorm, w, index, mode)
    
    !                     T
    !     MINIMIZE   1/2 X X    SUBJECT TO   G * X >= H.
    
    !       C.L. LAWSON, R.J. HANSON: 'SOLVING LEAST SQUARES PROBLEMS'
    !       PRENTICE HALL, ENGLEWOOD CLIFFS, NEW JERSEY, 1974.
    
    !     PARAMETER DESCRIPTION:
    
    !     G(),MG,M,N   ON ENTRY G() STORES THE M BY N MATRIX OF
    !                  LINEAR INEQUALITY CONSTRAINTS. G() HAS FIRST
    !                  DIMENSIONING PARAMETER MG
    !     H()          ON ENTRY H() STORES THE M VECTOR H REPRESENTING
    !                  THE RIGHT SIDE OF THE INEQUALITY SYSTEM
    
    !     REMARK: G(),H() WILL NOT BE CHANGED DURING CALCULATIONS BY LDP
    
    !     X()          ON ENTRY X() NEED NOT BE INITIALIZED.
    !                  ON EXIT X() STORES THE SOLUTION VECTOR X IF MODE=1.
    !     XNORM        ON EXIT XNORM STORES THE EUCLIDIAN NORM OF THE
    !                  SOLUTION VECTOR IF COMPUTATION IS SUCCESSFUL
    !     W()          W IS A ONE DIMENSIONAL WORKING SPACE, THE LENGTH
    !                  OF WHICH SHOULD BE AT LEAST (M+2)*(N+1) + 2*M
    !                  ON EXIT W() STORES THE LAGRANGE MULTIPLIERS
    !                  ASSOCIATED WITH THE CONSTRAINTS
    !                  AT THE SOLUTION OF PROBLEM LDP
    !     INDEX()      INDEX() IS A ONE DIMENSIONAL INTEGER WORKING SPACE
    !                  OF LENGTH AT LEAST M
    !     MODE         MODE IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING
    !                  MEANINGS:
    !          MODE=1: SUCCESSFUL COMPUTATION
    !               2: ERROR RETURN BECAUSE OF WRONG DIMENSIONS (N.LE.0)
    !               3: ITERATION COUNT EXCEEDED BY NNLS
    !               4: INEQUALITY CONSTRAINTS INCOMPATIBLE
    
      Double Precision g, h, x, xnorm, w, u, v, zero, one, fac, rnorm, dnrm2_, ddot_sl, diff
      Integer index, i, if, iw, iwdual, iy, iz, j, m, mg, mode, n, n1
      Dimension g(mg, n), h(m), x(n), w(*), index(m)
      diff(u, v) = u - v
      Data zero, one/0.0D0, 1.0D0/
    
      mode = 2
      If (n<=0) Goto 50
    
    !  STATE DUAL PROBLEM
    
      mode = 1
      x(1) = zero
      Call dcopy_(n, x(1), 0, x, 1)
      xnorm = zero
      If (m==0) Goto 50
      iw = 0
      Do j = 1, m
        Do i = 1, n
          iw = iw + 1
          w(iw) = g(j, i)
        End Do
        iw = iw + 1
        w(iw) = h(j)
      End Do
      if = iw + 1
      Do i = 1, n
        iw = iw + 1
        w(iw) = zero
      End Do
      w(iw+1) = one
      n1 = n + 1
      iz = iw + 2
      iy = iz + n1
      iwdual = iy + m
    
    !  SOLVE DUAL PROBLEM
    
      Call nnls(w, n1, n1, m, w(if), w(iy), rnorm, w(iwdual), w(iz), index, mode)
    
      If (mode/=1) Goto 50
      mode = 4
      If (rnorm<=zero) Goto 50
    
    !  COMPUTE SOLUTION OF PRIMAL PROBLEM
    
      fac = one - ddot_sl(m, h, 1, w(iy), 1)
      If (.Not. (diff(one+fac,one)>zero)) Goto 50
      mode = 1
      fac = one/fac
      Do j = 1, n
        x(j) = fac*ddot_sl(m, g(1,j), 1, w(iy), 1)
      End Do
      xnorm = dnrm2_(n, x, 1)
    
    !  COMPUTE LAGRANGE MULTIPLIERS FOR PRIMAL PROBLEM
    
      w(1) = zero
      Call dcopy_(m, w(1), 0, w, 1)
      Call daxpy_sl(m, fac, w(iy), 1, w, 1)
    
    !  END OF SUBROUTINE LDP
    
    50 End Subroutine ldp
    
    
    Subroutine nnls(a, mda, m, n, b, x, rnorm, w, z, index, mode)
    
    !     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY:
    !     'SOLVING LEAST SQUARES PROBLEMS'. PRENTICE-HALL.1974
    
    !      **********   NONNEGATIVE LEAST SQUARES   **********
    
    !     GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B, COMPUTE AN
    !     N-VECTOR, X, WHICH SOLVES THE LEAST SQUARES PROBLEM
    
    !                  A*X = B  SUBJECT TO  X >= 0
    
    !     A(),MDA,M,N
    !            MDA IS THE FIRST DIMENSIONING PARAMETER FOR THE ARRAY,A().
    !            ON ENTRY A()  CONTAINS THE M BY N MATRIX,A.
    !            ON EXIT A() CONTAINS THE PRODUCT Q*A,
    !            WHERE Q IS AN M BY M ORTHOGONAL MATRIX GENERATED
    !            IMPLICITLY BY THIS SUBROUTINE.
    !            EITHER M>=N OR M<N IS PERMISSIBLE.
    !            THERE IS NO RESTRICTION ON THE RANK OF A.
    !     B()    ON ENTRY B() CONTAINS THE M-VECTOR, B.
    !            ON EXIT B() CONTAINS Q*B.
    !     X()    ON ENTRY X() NEED NOT BE INITIALIZED.
    !            ON EXIT X() WILL CONTAIN THE SOLUTION VECTOR.
    !     RNORM  ON EXIT RNORM CONTAINS THE EUCLIDEAN NORM OF THE
    !            RESIDUAL VECTOR.
    !     W()    AN N-ARRAY OF WORKING SPACE.
    !            ON EXIT W() WILL CONTAIN THE DUAL SOLUTION VECTOR.
    !            W WILL SATISFY W(I)=0 FOR ALL I IN SET P
    !            AND W(I)<=0 FOR ALL I IN SET Z
    !     Z()    AN M-ARRAY OF WORKING SPACE.
    !     INDEX()AN INTEGER WORKING ARRAY OF LENGTH AT LEAST N.
    !            ON EXIT THE CONTENTS OF THIS ARRAY DEFINE THE SETS
    !            P AND Z AS FOLLOWS:
    !            INDEX(1)    THRU INDEX(NSETP) = SET P.
    !            INDEX(IZ1)  THRU INDEX (IZ2)  = SET Z.
    !            IZ1=NSETP + 1 = NPP1, IZ2=N.
    !     MODE   THIS IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANING:
    !            1    THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY.
    !            2    THE DIMENSIONS OF THE PROBLEM ARE WRONG,
    !                 EITHER M <= 0 OR N <= 0.
    !            3    ITERATION COUNT EXCEEDED, MORE THAN 3*N ITERATIONS.
    
      Integer i, ii, ip, iter, itmax, iz, izmax, iz1, iz2, j, jj, jz, k, l, m, mda, mode, n, npp1, nsetp, index(n)
    
      Double Precision a(mda, n), b(m), x(n), w(n), z(m), asave, diff, factor, ddot_sl, zero, one, wmax, alpha, c, s, t, u, v, up, rnorm, unorm, dnrm2_
    
      diff(u, v) = u - v
    
      Data zero, one, factor/0.0D0, 1.0D0, 1.0D-2/
    
    !     revised          Dieter Kraft, March 1983
    
      mode = 2
      If (m<=0 .Or. n<=0) Goto 290
      mode = 1
      iter = 0
      itmax = 3*n
    
    ! STEP ONE (INITIALIZE)
    
      Do i = 1, n
        index(i) = i
      End Do
      iz1 = 1
      iz2 = n
      nsetp = 0
      npp1 = 1
      x(1) = zero
      Call dcopy_(n, x(1), 0, x, 1)
    
    ! STEP TWO (COMPUTE DUAL VARIABLES)
    ! .....ENTRY LOOP A
    
      110 If (iz1>iz2 .Or. nsetp>=m) Goto 280
      Do iz = iz1, iz2
        j = index(iz)
        w(j) = ddot_sl(m-nsetp, a(npp1,j), 1, b(npp1), 1)
      End Do
    
    ! STEP THREE (TEST DUAL VARIABLES)
    
      130 wmax = zero
      Do iz = iz1, iz2
        j = index(iz)
        If (w(j)<=wmax) Goto 140
        wmax = w(j)
        izmax = iz
      140 End Do
    
    ! .....EXIT LOOP A
    
      If (wmax<=zero) Goto 280
      iz = izmax
      j = index(iz)
    
    ! STEP FOUR (TEST INDEX J FOR LINEAR DEPENDENCY)
    
      asave = a(npp1, j)
      Call h12(1, npp1, npp1+1, m, a(1,j), 1, up, z, 1, 1, 0)
      unorm = dnrm2_(nsetp, a(1,j), 1)
      t = factor*abs(a(npp1,j))
      If (diff(unorm+t,unorm)<=zero) Goto 150
      Call dcopy_(m, b, 1, z, 1)
      Call h12(2, npp1, npp1+1, m, a(1,j), 1, up, z, 1, 1, 1)
      If (z(npp1)/a(npp1,j)>zero) Goto 160
      150 a(npp1, j) = asave
      w(j) = zero
      Goto 130
    ! STEP FIVE (ADD COLUMN)
    
      160 Call dcopy_(m, z, 1, b, 1)
      index(iz) = index(iz1)
      index(iz1) = j
      iz1 = iz1 + 1
      nsetp = npp1
      npp1 = npp1 + 1
      Do jz = iz1, iz2
        jj = index(jz)
        Call h12(2, nsetp, npp1, m, a(1,j), 1, up, a(1,jj), 1, mda, 1)
      End Do
      k = min(npp1, mda)
      w(j) = zero
      Call dcopy_(m-nsetp, w(j), 0, a(k,j), 1)
    
    ! STEP SIX (SOLVE LEAST SQUARES SUB-PROBLEM)
    ! .....ENTRY LOOP B
    
      180 Do ip = nsetp, 1, -1
        If (ip==nsetp) Goto 190
        Call daxpy_sl(ip, -z(ip+1), a(1,jj), 1, z, 1)
        190 jj = index(ip)
        z(ip) = z(ip)/a(ip, jj)
      End Do
      iter = iter + 1
      If (iter<=itmax) Goto 220
      210 mode = 3
      Goto 280
    ! STEP SEVEN TO TEN (STEP LENGTH ALGORITHM)
    
      220 alpha = one
      jj = 0
      Do ip = 1, nsetp
        If (z(ip)>zero) Goto 230
        l = index(ip)
        t = -x(l)/(z(ip)-x(l))
        If (alpha<t) Goto 230
        alpha = t
        jj = ip
      230 End Do
      Do ip = 1, nsetp
        l = index(ip)
        x(l) = (one-alpha)*x(l) + alpha*z(ip)
      End Do
    
    ! .....EXIT LOOP B
    
      If (jj==0) Goto 110
    
    ! STEP ELEVEN (DELETE COLUMN)
    
      i = index(jj)
      250 x(i) = zero
      jj = jj + 1
      Do j = jj, nsetp
        ii = index(j)
        index(j-1) = ii
        Call dsrotg(a(j-1,ii), a(j,ii), c, s)
        t = a(j-1, ii)
        Call dsrot(n, a(j-1,1), mda, a(j,1), mda, c, s)
        a(j-1, ii) = t
        a(j, ii) = zero
        Call dsrot(1, b(j-1), 1, b(j), 1, c, s)
      End Do
      npp1 = nsetp
      nsetp = nsetp - 1
      iz1 = iz1 - 1
      index(iz1) = i
      If (nsetp<=0) Goto 210
      Do jj = 1, nsetp
        i = index(jj)
        If (x(i)<=zero) Goto 250
      End Do
      Call dcopy_(m, b, 1, z, 1)
      Goto 180
    ! STEP TWELVE (SOLUTION)
    
      280 k = min(npp1, m)
      rnorm = dnrm2_(m-nsetp, b(k), 1)
      If (npp1>m) Then
        w(1) = zero
        Call dcopy_(n, w(1), 0, w, 1)
      End If
    
    ! END OF SUBROUTINE NNLS
    
    290 End Subroutine nnls
    
    Subroutine hfti(a, mda, m, n, b, mdb, nb, tau, krank, rnorm, h, g, ip)
    
    !     RANK-DEFICIENT LEAST SQUARES ALGORITHM AS DESCRIBED IN:
    !     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY, 1973 JUN 12
    !     TO APPEAR IN 'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974
    
    !     A(*,*),MDA,M,N   THE ARRAY A INITIALLY CONTAINS THE M x N MATRIX A
    !                      OF THE LEAST SQUARES PROBLEM AX = B.
    !                      THE FIRST DIMENSIONING PARAMETER MDA MUST SATISFY
    !                      MDA >= M. EITHER M >= N OR M < N IS PERMITTED.
    !                      THERE IS NO RESTRICTION ON THE RANK OF A.
    !                      THE MATRIX A WILL BE MODIFIED BY THE SUBROUTINE.
    !     B(*,*),MDB,NB    IF NB = 0 THE SUBROUTINE WILL MAKE NO REFERENCE
    !                      TO THE ARRAY B. IF NB > 0 THE ARRAY B() MUST
    !                      INITIALLY CONTAIN THE M x NB MATRIX B  OF THE
    !                      THE LEAST SQUARES PROBLEM AX = B AND ON RETURN
    !                      THE ARRAY B() WILL CONTAIN THE N x NB SOLUTION X.
    !                      IF NB>1 THE ARRAY B() MUST BE DOUBLE SUBSCRIPTED
    !                      WITH FIRST DIMENSIONING PARAMETER MDB>=MAX(M,N),
    !                      IF NB=1 THE ARRAY B() MAY BE EITHER SINGLE OR
    !                      DOUBLE SUBSCRIPTED.
    !     TAU              ABSOLUTE TOLERANCE PARAMETER FOR PSEUDORANK
    !                      DETERMINATION, PROVIDED BY THE USER.
    !     KRANK            PSEUDORANK OF A, SET BY THE SUBROUTINE.
    !     RNORM            ON EXIT, RNORM(J) WILL CONTAIN THE EUCLIDIAN
    !                      NORM OF THE RESIDUAL VECTOR FOR THE PROBLEM
    !                      DEFINED BY THE J-TH COLUMN VECTOR OF THE ARRAY B.
    !     H(), G()         ARRAYS OF WORKING SPACE OF LENGTH >= N.
    !     IP()             INTEGER ARRAY OF WORKING SPACE OF LENGTH >= N
    !                      RECORDING PERMUTATION INDICES OF COLUMN VECTORS
    
      Integer i, j, jb, k, kp1, krank, l, ldiag, lmax, m, mda, mdb, n, nb, ip(n)
      Double Precision a(mda, n), b(mdb, nb), h(n), g(n), rnorm(nb), factor, tau, zero, hmax, diff, tmp, ddot_sl, dnrm2_, u, v
      diff(u, v) = u - v
      Data zero/0.0D0/, factor/1.0D-3/
    
      k = 0
      ldiag = min(m, n)
      If (ldiag<=0) Goto 270
    
    !   COMPUTE LMAX
    
      Do j = 1, ldiag
        If (j==1) Goto 20
        lmax = j
        Do l = j, n
          h(l) = h(l) - a(j-1, l)**2
          If (h(l)>h(lmax)) lmax = l
        End Do
        If (diff(hmax+factor*h(lmax),hmax)>zero) Goto 50
        20 lmax = j
        Do l = j, n
          h(l) = zero
          Do i = j, m
            h(l) = h(l) + a(i, l)**2
          End Do
          If (h(l)>h(lmax)) lmax = l
        End Do
        hmax = h(lmax)
    
    !   COLUMN INTERCHANGES IF NEEDED
    
        50 ip(j) = lmax
        If (ip(j)==j) Goto 70
        Do i = 1, m
          tmp = a(i, j)
          a(i, j) = a(i, lmax)
          a(i, lmax) = tmp
        End Do
        h(lmax) = h(j)
    
    !   J-TH TRANSFORMATION AND APPLICATION TO A AND B
    
        70 i = min(j+1, n)
        Call h12(1, j, j+1, m, a(1,j), 1, h(j), a(1,i), 1, mda, n-j)
        Call h12(2, j, j+1, m, a(1,j), 1, h(j), b, 1, mdb, nb)
      End Do
    
    !   DETERMINE PSEUDORANK
    
      Do j = 1, ldiag
        If (abs(a(j,j))<=tau) Goto 100
      End Do
      k = ldiag
      Goto 110
      100 k = j - 1
      110 kp1 = k + 1
    
    !   NORM OF RESIDUALS
    
      Do jb = 1, nb
        rnorm(jb) = dnrm2_(m-k, b(kp1,jb), 1)
      End Do
      If (k>0) Goto 160
      Do jb = 1, nb
        Do i = 1, n
          b(i, jb) = zero
        End Do
      End Do
      Goto 270
      160 If (k==n) Goto 180
    
    !   HOUSEHOLDER DECOMPOSITION OF FIRST K ROWS
    
      Do i = k, 1, -1
        Call h12(1, i, kp1, n, a(i,1), mda, g(i), a, mda, 1, i-1)
      End Do
      180 Do jb = 1, nb
    
    !   SOLVE K*K TRIANGULAR SYSTEM
    
        Do i = k, 1, -1
          j = min(i+1, n)
          b(i, jb) = (b(i,jb)-ddot_sl(k-i,a(i,j),mda,b(j,jb),1))/a(i, i)
        End Do
    
    !   COMPLETE SOLUTION VECTOR
    
        If (k==n) Goto 240
        Do j = kp1, n
          b(j, jb) = zero
        End Do
        Do i = 1, k
          Call h12(2, i, kp1, n, a(i,1), mda, g(i), b(1,jb), 1, mdb, 1)
        End Do
    
    !   REORDER SOLUTION ACCORDING TO PREVIOUS COLUMN INTERCHANGES
    
        240 Do j = ldiag, 1, -1
          If (ip(j)==j) Goto 250
          l = ip(j)
          tmp = b(l, jb)
          b(l, jb) = b(j, jb)
          b(j, jb) = tmp
        250 End Do
      End Do
      270 krank = k
    End Subroutine hfti
    
    Subroutine h12(mode, lpivot, l1, m, u, iue, up, c, ice, icv, ncv)
    
    !     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY, 1973 JUN 12
    !     TO APPEAR IN 'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974
    
    !     CONSTRUCTION AND/OR APPLICATION OF A SINGLE
    !     HOUSEHOLDER TRANSFORMATION  Q = I + U*(U**T)/B
    
    !     MODE    = 1 OR 2   TO SELECT ALGORITHM  H1  OR  H2 .
    !     LPIVOT IS THE INDEX OF THE PIVOT ELEMENT.
    !     L1,M   IF L1 <= M   THE TRANSFORMATION WILL BE CONSTRUCTED TO
    !            ZERO ELEMENTS INDEXED FROM L1 THROUGH M.
    !            IF L1 > M THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
    !     U(),IUE,UP
    !            ON ENTRY TO H1 U() STORES THE PIVOT VECTOR.
    !            IUE IS THE STORAGE INCREMENT BETWEEN ELEMENTS.
    !            ON EXIT FROM H1 U() AND UP STORE QUANTITIES DEFINING
    !            THE VECTOR U OF THE HOUSEHOLDER TRANSFORMATION.
    !            ON ENTRY TO H2 U() AND UP
    !            SHOULD STORE QUANTITIES PREVIOUSLY COMPUTED BY H1.
    !            THESE WILL NOT BE MODIFIED BY H2.
    !     C()    ON ENTRY TO H1 OR H2 C() STORES A MATRIX WHICH WILL BE
    !            REGARDED AS A SET OF VECTORS TO WHICH THE HOUSEHOLDER
    !            TRANSFORMATION IS TO BE APPLIED.
    !            ON EXIT C() STORES THE SET OF TRANSFORMED VECTORS.
    !     ICE    STORAGE INCREMENT BETWEEN ELEMENTS OF VECTORS IN C().
    !     ICV    STORAGE INCREMENT BETWEEN VECTORS IN C().
    !     NCV    NUMBER OF VECTORS IN C() TO BE TRANSFORMED.
    !            IF NCV <= 0 NO OPERATIONS WILL BE DONE ON C().
    
      Integer incr, ice, icv, iue, lpivot, l1, mode, ncv
      Integer i, i2, i3, i4, j, m
      Double Precision u, up, c, cl, clinv, b, sm, one, zero
      Dimension u(iue, *), c(*)
      Data one/1.0D+00/, zero/0.0D+00/
    
      If (0>=lpivot .Or. lpivot>=l1 .Or. l1>m) Goto 80
      cl = abs(u(1,lpivot))
      If (mode==2) Goto 30
    
    !     ****** CONSTRUCT THE TRANSFORMATION ******
    
      Do j = l1, m
        sm = abs(u(1,j))
        cl = max(sm, cl)
      End Do
      If (cl<=zero) Goto 80
      clinv = one/cl
      sm = (u(1,lpivot)*clinv)**2
      Do j = l1, m
        sm = sm + (u(1,j)*clinv)**2
      End Do
      cl = cl*sqrt(sm)
      If (u(1,lpivot)>zero) cl = -cl
      up = u(1, lpivot) - cl
      u(1, lpivot) = cl
      Goto 40
    !     ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C ******
    
      30 If (cl<=zero) Goto 80
      40 If (ncv<=0) Goto 80
      b = up*u(1, lpivot)
      If (b>=zero) Goto 80
      b = one/b
      i2 = 1 - icv + ice*(lpivot-1)
      incr = ice*(l1-lpivot)
      Do j = 1, ncv
        i2 = i2 + icv
        i3 = i2 + incr
        i4 = i3
        sm = c(i2)*up
        Do i = l1, m
          sm = sm + c(i3)*u(1, i)
          i3 = i3 + ice
        End Do
        If (sm==zero) Goto 70
        sm = sm*b
        c(i2) = c(i2) + sm*up
        Do i = l1, m
          c(i4) = c(i4) + sm*u(1, i)
          i4 = i4 + ice
        End Do
      70 End Do
    80 End Subroutine h12
    
    Subroutine ldl(n, a, z, sigma, w)
    !   LDL     LDL' - RANK-ONE - UPDATE
    
    !   PURPOSE:
    !           UPDATES THE LDL' FACTORS OF MATRIX A BY RANK-ONE MATRIX
    !           SIGMA*Z*Z'
    
    !   INPUT ARGUMENTS: (* MEANS PARAMETERS ARE CHANGED DURING EXECUTION)
    !     N     : ORDER OF THE COEFFICIENT MATRIX A
    !   * A     : POSITIVE DEFINITE MATRIX OF DIMENSION N;
    !             ONLY THE LOWER TRIANGLE IS USED AND IS STORED COLUMN BY
    !             COLUMN AS ONE DIMENSIONAL ARRAY OF DIMENSION N*(N+1)/2.
    !   * Z     : VECTOR OF DIMENSION N OF UPDATING ELEMENTS
    !     SIGMA : SCALAR FACTOR BY WHICH THE MODIFYING DYADE Z*Z' IS
    !             MULTIPLIED
    
    !   OUTPUT ARGUMENTS:
    !     A     : UPDATED LDL' FACTORS
    
    !   WORKING ARRAY:
    !     W     : VECTOR OP DIMENSION N (USED ONLY IF SIGMA .LT. ZERO)
    
    !   METHOD:
    !     THAT OF FLETCHER AND POWELL AS DESCRIBED IN :
    !     FLETCHER,R.,(1974) ON THE MODIFICATION OF LDL' FACTORIZATION.
    !     POWELL,M.J.D.      MATH.COMPUTATION 28, 1067-1078.
    
    !   IMPLEMENTED BY:
    !     KRAFT,D., DFVLR - INSTITUT FUER DYNAMIK DER FLUGSYSTEME
    !               D-8031  OBERPFAFFENHOFEN
    
    !   STATUS: 15. JANUARY 1980
    
    !   SUBROUTINES REQUIRED: NONE
    
      Integer i, ij, j, n
      Double Precision a(*), t, v, w(*), z(*), u, tp, one, beta, four, zero, alpha, delta, gamma, sigma, epmach
      Data zero, one, four, epmach/0.0D0, 1.0D0, 4.0D0, 2.22D-16/
    
      If (sigma==zero) Goto 280
      ij = 1
      t = one/sigma
      If (sigma>zero) Goto 220
    ! PREPARE NEGATIVE UPDATE
      Do i = 1, n
        w(i) = z(i)
      End Do
      Do i = 1, n
        v = w(i)
        t = t + v*v/a(ij)
        Do j = i + 1, n
          ij = ij + 1
          w(j) = w(j) - v*a(ij)
        End Do
        ij = ij + 1
      End Do
      If (t>=zero) t = epmach/sigma
      Do i = 1, n
        j = n + 1 - i
        ij = ij - i
        u = w(j)
        w(j) = t
        t = t - u*u/a(ij)
      End Do
      220 Continue
    ! HERE UPDATING BEGINS
      Do i = 1, n
        v = z(i)
        delta = v/a(ij)
        If (sigma<zero) tp = w(i)
        If (sigma>zero) tp = t + delta*v
        alpha = tp/t
        a(ij) = alpha*a(ij)
        If (i==n) Goto 280
        beta = delta/tp
        If (alpha>four) Goto 240
        Do j = i + 1, n
          ij = ij + 1
          z(j) = z(j) - v*a(ij)
          a(ij) = a(ij) + beta*z(j)
        End Do
        Goto 260
        240 gamma = t/tp
        Do j = i + 1, n
          ij = ij + 1
          u = a(ij)
          a(ij) = gamma*u + beta*z(j)
          z(j) = z(j) - v*u
        End Do
        260 ij = ij + 1
        t = tp
      End Do
      280 Return
    ! END OF LDL
    End Subroutine ldl
    
    Double Precision Function linmin(mode, ax, bx, f, tol)
    !   LINMIN  LINESEARCH WITHOUT DERIVATIVES
    
    !   PURPOSE:
    
    !  TO FIND THE ARGUMENT LINMIN WHERE THE FUNCTION F TAKES IT'S MINIMUM
    !  ON THE INTERVAL AX, BX.
    !  COMBINATION OF GOLDEN SECTION AND SUCCESSIVE QUADRATIC INTERPOLATION.
    
    !   INPUT ARGUMENTS: (* MEANS PARAMETERS ARE CHANGED DURING EXECUTION)
    
    ! *MODE   SEE OUTPUT ARGUMENTS
    !  AX     LEFT ENDPOINT OF INITIAL INTERVAL
    !  BX     RIGHT ENDPOINT OF INITIAL INTERVAL
    !  F      FUNCTION VALUE AT LINMIN WHICH IS TO BE BROUGHT IN BY
    !         REVERSE COMMUNICATION CONTROLLED BY MODE
    !  TOL    DESIRED LENGTH OF INTERVAL OF UNCERTAINTY OF FINAL RESULT
    
    !   OUTPUT ARGUMENTS:
    
    !  LINMIN ABSCISSA APPROXIMATING THE POINT WHERE F ATTAINS A MINIMUM
    !  MODE   CONTROLS REVERSE COMMUNICATION
    !         MUST BE SET TO 0 INITIALLY, RETURNS WITH INTERMEDIATE
    !         VALUES 1 AND 2 WHICH MUST NOT BE CHANGED BY THE USER,
    !         ENDS WITH CONVERGENCE WITH VALUE 3.
    
    !   WORKING ARRAY:
    
    !  NONE
    
    !   METHOD:
    
    !  THIS FUNCTION SUBPROGRAM IS A SLIGHTLY MODIFIED VERSION OF THE
    !  ALGOL 60 PROCEDURE LOCALMIN GIVEN IN
    !  R.P. BRENT: ALGORITHMS FOR MINIMIZATION WITHOUT DERIVATIVES,
    !              PRENTICE-HALL (1973).
    
    !   IMPLEMENTED BY:
    
    !     KRAFT, D., DFVLR - INSTITUT FUER DYNAMIK DER FLUGSYSTEME
    !                D-8031  OBERPFAFFENHOFEN
    
    !   STATUS: 31. AUGUST  1984
    
    !   SUBROUTINES REQUIRED: NONE
    
      Integer mode
      Double Precision f, tol, a, b, c, d, e, p, q, r, u, v, w, x, m, fu, fv, fw, fx, eps, tol1, tol2, zero, ax, bx
      Data c/0.381966011D0/, eps/1.5D-8/, zero/0.0D0/
    
    !  EPS = SQUARE - ROOT OF MACHINE PRECISION
    !  C = GOLDEN SECTION RATIO = (3-SQRT(5))/2
    
      Goto (10, 55), mode
    
    !  INITIALIZATION
    
      a = ax
      b = bx
      e = zero
      v = a + c*(b-a)
      w = v
      x = w
      linmin = x
      mode = 1
      Goto 100
    
    !  MAIN LOOP STARTS HERE
    
      10 fx = f
      fv = fx
      fw = fv
      20 m = 0.5D0*(a+b)
      tol1 = eps*abs(x) + tol
      tol2 = tol1 + tol1
    
    !  TEST CONVERGENCE
    
      If (abs(x-m)<=tol2-0.5D0*(b-a)) Goto 90
      r = zero
      q = r
      p = q
      If (abs(e)<=tol1) Goto 30
    
    !  FIT PARABOLA
    
      r = (x-w)*(fx-fv)
      q = (x-v)*(fx-fw)
      p = (x-v)*q - (x-w)*r
      q = q - r
      q = q + q
      If (q>zero) p = -p
      If (q<zero) q = -q
      r = e
      e = d
    
    !  IS PARABOLA ACCEPTABLE
    
      30 If (abs(p)>=0.5D0*abs(q*r) .Or. p<=q*(a-x) .Or. p>=q*(b-x)) Goto 40
    
    !  PARABOLIC INTERPOLATION STEP
    
      d = p/q
    
    !  F MUST NOT BE EVALUATED TOO CLOSE TO A OR B
    
      If (u-a<tol2) d = sign(tol1, m-x)
      If (b-u<tol2) d = sign(tol1, m-x)
      Goto 50
    
    !  GOLDEN SECTION STEP
    
      40 If (x>=m) e = a - x
      If (x<m) e = b - x
      d = c*e
    
    !  F MUST NOT BE EVALUATED TOO CLOSE TO X
    
      50 If (abs(d)<tol1) d = sign(tol1, d)
      u = x + d
      linmin = u
      mode = 2
      Goto 100
      55 fu = f
    
    !  UPDATE A, B, V, W, AND X
    
      If (fu>fx) Goto 60
      If (u>=x) a = x
      If (u<x) b = x
      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu
      Goto 85
      60 If (u<x) a = u
      If (u>=x) b = u
      If (fu<=fw .Or. w==x) Goto 70
      If (fu<=fv .Or. v==x .Or. v==w) Goto 80
      Goto 85
      70 v = w
      fv = fw
      w = u
      fw = fu
      Goto 85
      80 v = u
      fv = fu
      85 Goto 20
    
    !  END OF MAIN LOOP
    
      90 linmin = x
      mode = 3
      100 Return
    
    !  END OF LINMIN
    
    End Function linmin
    
    !## Following a selection from BLAS Level 1
    
    Subroutine daxpy_sl(n, da, dx, incx, dy, incy)
    
    !     CONSTANT TIMES A VECTOR PLUS A VECTOR.
    !     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
    !     JACK DONGARRA, LINPACK, 3/11/78.
    
      Double Precision dx(*), dy(*), da
      Integer i, incx, incy, ix, iy, m, mp1, n
    
      If (n<=0) Return
      If (da==0.0D0) Return
      If (incx==1 .And. incy==1) Goto 20
    
    !        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
    !        NOT EQUAL TO 1
    
      ix = 1
      iy = 1
      If (incx<0) ix = (-n+1)*incx + 1
      If (incy<0) iy = (-n+1)*incy + 1
      Do i = 1, n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
      End Do
      Return
    
    !        CODE FOR BOTH INCREMENTS EQUAL TO 1
    
    !        CLEAN-UP LOOP
    
      20 m = mod(n, 4)
      If (m==0) Goto 40
      Do i = 1, m
        dy(i) = dy(i) + da*dx(i)
      End Do
      If (n<4) Return
      40 mp1 = m + 1
      Do i = mp1, n, 4
        dy(i) = dy(i) + da*dx(i)
        dy(i+1) = dy(i+1) + da*dx(i+1)
        dy(i+2) = dy(i+2) + da*dx(i+2)
        dy(i+3) = dy(i+3) + da*dx(i+3)
      End Do
      Return
    End Subroutine daxpy_sl
    
    Subroutine dcopy_(n, dx, incx, dy, incy)
    
    !     COPIES A VECTOR, X, TO A VECTOR, Y.
    !     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
    !     JACK DONGARRA, LINPACK, 3/11/78.
    
      Double Precision dx(*), dy(*)
      Integer i, incx, incy, ix, iy, m, mp1, n
    
      If (n<=0) Return
      If (incx==1 .And. incy==1) Goto 20
    
    !        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
    !        NOT EQUAL TO 1
    
      ix = 1
      iy = 1
      If (incx<0) ix = (-n+1)*incx + 1
      If (incy<0) iy = (-n+1)*incy + 1
      Do i = 1, n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
      End Do
      Return
    
    !        CODE FOR BOTH INCREMENTS EQUAL TO 1
    
    !        CLEAN-UP LOOP
    
      20 m = mod(n, 7)
      If (m==0) Goto 40
      Do i = 1, m
        dy(i) = dx(i)
      End Do
      If (n<7) Return
      40 mp1 = m + 1
      Do i = mp1, n, 7
        dy(i) = dx(i)
        dy(i+1) = dx(i+1)
        dy(i+2) = dx(i+2)
        dy(i+3) = dx(i+3)
        dy(i+4) = dx(i+4)
        dy(i+5) = dx(i+5)
        dy(i+6) = dx(i+6)
      End Do
      Return
    End Subroutine dcopy_
    
    Double Precision Function ddot_sl(n, dx, incx, dy, incy)
    
    !     FORMS THE DOT PRODUCT OF TWO VECTORS.
    !     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
    !     JACK DONGARRA, LINPACK, 3/11/78.
    
      Double Precision dx(*), dy(*), dtemp
      Integer i, incx, incy, ix, iy, m, mp1, n
    
      ddot_sl = 0.0D0
      dtemp = 0.0D0
      If (n<=0) Return
      If (incx==1 .And. incy==1) Goto 20
    
    !        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
    !          NOT EQUAL TO 1
    
      ix = 1
      iy = 1
      If (incx<0) ix = (-n+1)*incx + 1
      If (incy<0) iy = (-n+1)*incy + 1
      Do i = 1, n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
      End Do
      ddot_sl = dtemp
      Return
    
    !        CODE FOR BOTH INCREMENTS EQUAL TO 1
    
    !        CLEAN-UP LOOP
    
      20 m = mod(n, 5)
      If (m==0) Goto 40
      Do i = 1, m
        dtemp = dtemp + dx(i)*dy(i)
      End Do
      If (n<5) Goto 60
      40 mp1 = m + 1
      Do i = mp1, n, 5
        dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) + dx(i+2)*dy(i+2) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
      End Do
      60 ddot_sl = dtemp
      Return
    End Function ddot_sl
    
    Double Precision Function dnrm1(n, x, i, j)
      Integer n, i, j, k
      Double Precision snormx, sum, x(n), zero, one, scale, temp
      Data zero/0.0D0/, one/1.0D0/
    
    !      DNRM1 - COMPUTES THE I-NORM OF A VECTOR
    !              BETWEEN THE ITH AND THE JTH ELEMENTS
    
    !      INPUT -
    !      N       LENGTH OF VECTOR
    !      X       VECTOR OF LENGTH N
    !      I       INITIAL ELEMENT OF VECTOR TO BE USED
    !      J       FINAL ELEMENT TO USE
    
    !      OUTPUT -
    !      DNRM1   NORM
    
      snormx = zero
      Do k = i, j
        snormx = max(snormx, abs(x(k)))
      End Do
      dnrm1 = snormx
      If (snormx==zero) Return
      scale = snormx
      If (snormx>=one) scale = sqrt(snormx)
      sum = zero
      Do k = i, j
        temp = zero
        If (abs(x(k))+scale/=scale) temp = x(k)/snormx
        If (one+temp/=one) sum = sum + temp*temp
      End Do
      sum = sqrt(sum)
      dnrm1 = snormx*sum
      Return
    End Function dnrm1
    
    Double Precision Function dnrm2_(n, dx, incx)
      Integer n, i, j, nn, next, incx
      Double Precision dx(*), cutlo, cuthi, hitest, sum, xmax, zero, one
      Data zero, one/0.0D0, 1.0D0/
    
    !     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE
    !     INCREMENT INCX .
    !     IF    N .LE. 0 RETURN WITH RESULT = 0.
    !     IF N .GE. 1 THEN INCX MUST BE .GE. 1
    
    !           C.L.LAWSON, 1978 JAN 08
    
    !     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
    !     HOPEFULLY APPLICABLE TO ALL MACHINES.
    !         CUTLO = MAXIMUM OF  SQRT(U/EPS)   OVER ALL KNOWN MACHINES.
    !         CUTHI = MINIMUM OF  SQRT(V)       OVER ALL KNOWN MACHINES.
    !     WHERE
    !         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
    !         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
    !         V   = LARGEST  NO.            (OVERFLOW  LIMIT)
    
    !     BRIEF OUTLINE OF ALGORITHM..
    
    !     PHASE 1    SCANS ZERO COMPONENTS.
    !     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
    !     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
    !     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
    !     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
    
    !     VALUES FOR CUTLO AND CUTHI..
    !     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
    !     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
    !     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
    !                   UNIVAC AND DEC AT 2**(-103)
    !                   THUS CUTLO = 2**(-51) = 4.44089E-16
    !     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
    !                   THUS CUTHI = 2**(63.5) = 1.30438E19
    !     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
    !                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
    !     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
    !     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
    !     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
      Data cutlo, cuthi/8.232D-11, 1.304D19/
    
      If (n>0) Goto 10
      dnrm2_ = zero
      Goto 300
    
      10 Assign 30 To next
      sum = zero
      nn = n*incx
    !                       BEGIN MAIN LOOP
      i = 1
      20 Goto next, (30, 50, 70, 110)
      30 If (abs(dx(i))>cutlo) Goto 85
      Assign 50 To next
      xmax = zero
    
    !                        PHASE 1.  SUM IS ZERO
    
      50 If (dx(i)==zero) Goto 200
      If (abs(dx(i))>cutlo) Goto 85
    
    !                        PREPARE FOR PHASE 2.
    
      Assign 70 To next
      Goto 105
    
    !                        PREPARE FOR PHASE 4.
    
      100 i = j
      Assign 110 To next
      sum = (sum/dx(i))/dx(i)
      105 xmax = abs(dx(i))
      Goto 115
    
    !                   PHASE 2.  SUM IS SMALL.
    !                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
    
      70 If (abs(dx(i))>cutlo) Goto 75
    
    !                   COMMON CODE FOR PHASES 2 AND 4.
    !                   IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
    
      110 If (abs(dx(i))<=xmax) Goto 115
      sum = one + sum*(xmax/dx(i))**2
      xmax = abs(dx(i))
      Goto 200
    
      115 sum = sum + (dx(i)/xmax)**2
      Goto 200
    
    !                  PREPARE FOR PHASE 3.
    
      75 sum = (sum*xmax)*xmax
    
    !     FOR REAL OR D.P. SET HITEST = CUTHI/N
    !     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
    
      85 hitest = cuthi/float(n)
    
    !                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
    
      Do j = i, nn, incx
        If (abs(dx(j))>=hitest) Goto 100
        sum = sum + dx(j)**2
      End Do
      dnrm2_ = sqrt(sum)
      Goto 300
    
      200 Continue
      i = i + incx
      If (i<=nn) Goto 20
    
    !              END OF MAIN LOOP.
    
    !              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
    
      dnrm2_ = xmax*sqrt(sum)
      300 Continue
      Return
    End Function dnrm2_
    
    Subroutine dsrot(n, dx, incx, dy, incy, c, s)
    
    !     APPLIES A PLANE ROTATION.
    !     JACK DONGARRA, LINPACK, 3/11/78.
    
      Double Precision dx(*), dy(*), dtemp, c, s
      Integer i, incx, incy, ix, iy, n
    
      If (n<=0) Return
      If (incx==1 .And. incy==1) Goto 20
    
    !       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
    !         TO 1
    
      ix = 1
      iy = 1
      If (incx<0) ix = (-n+1)*incx + 1
      If (incy<0) iy = (-n+1)*incy + 1
      Do i = 1, n
        dtemp = c*dx(ix) + s*dy(iy)
        dy(iy) = c*dy(iy) - s*dx(ix)
        dx(ix) = dtemp
        ix = ix + incx
        iy = iy + incy
      End Do
      Return
    
    !       CODE FOR BOTH INCREMENTS EQUAL TO 1
    
      20 Do i = 1, n
        dtemp = c*dx(i) + s*dy(i)
        dy(i) = c*dy(i) - s*dx(i)
        dx(i) = dtemp
      End Do
      Return
    End Subroutine dsrot
    
    Subroutine dsrotg(da, db, c, s)
    
    !     CONSTRUCT GIVENS PLANE ROTATION.
    !     JACK DONGARRA, LINPACK, 3/11/78.
    !                    MODIFIED 9/27/86.
    
      Double Precision da, db, c, s, roe, scale, r, z, one, zero
      Data one, zero/1.0D+00, 0.0D+00/
    
      roe = db
      If (abs(da)>abs(db)) roe = da
      scale = abs(da) + abs(db)
      If (scale/=zero) Goto 10
      c = one
      s = zero
      r = zero
      Goto 20
      10 r = scale*sqrt((da/scale)**2+(db/scale)**2)
      r = sign(one, roe)*r
      c = da/r
      s = db/r
      20 z = s
      If (abs(c)>zero .And. abs(c)<=s) z = one/c
      da = r
      db = z
      Return
    End Subroutine dsrotg
    
    Subroutine dscal_sl(n, da, dx, incx)
    
    !     SCALES A VECTOR BY A CONSTANT.
    !     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
    !     JACK DONGARRA, LINPACK, 3/11/78.
    
      Double Precision da, dx(*)
      Integer i, incx, m, mp1, n, nincx
    
      If (n<=0) Return
      If (incx==1) Goto 20
    
    
    !        CODE FOR INCREMENT NOT EQUAL TO 1
    
      nincx = n*incx
      Do i = 1, nincx, incx
        dx(i) = da*dx(i)
      End Do
      Return
    
    !        CODE FOR INCREMENT EQUAL TO 1
    
    !        CLEAN-UP LOOP
    
      20 m = mod(n, 5)
      If (m==0) Goto 40
      Do i = 1, m
        dx(i) = da*dx(i)
      End Do
      If (n<5) Return
      40 mp1 = m + 1
      Do i = mp1, n, 5
        dx(i) = da*dx(i)
        dx(i+1) = da*dx(i+1)
        dx(i+2) = da*dx(i+2)
        dx(i+3) = da*dx(i+3)
        dx(i+4) = da*dx(i+4)
      End Do
      Return
    End Subroutine dscal_sl
    
    Subroutine bound(n, x, xl, xu)
      Integer n, i
      Double Precision x(n), xl(n), xu(n)
      Do i = 1, n
    !        Note that xl(i) and xu(i) may be NaN to indicate no bound
        If (xl(i)==xl(i) .And. x(i)<xl(i)) Then
          x(i) = xl(i)
        Else If (xu(i)==xu(i) .And. x(i)>xu(i)) Then
          x(i) = xu(i)
        End If
      End Do
    End Subroutine bound
    