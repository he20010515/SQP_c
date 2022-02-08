# 基于c语言实现的SQP算法仓库
```fortran
C   SLSQP       S EQUENTIAL  L EAST  SQ UARES  P ROGRAMMING
C            TO SOLVE GENERAL NONLINEAR OPTIMIZATION PROBLEMS

C*
C*                                                                     
C*                                                                     
C* A NONLINEAR PROGRAMMING METHOD WITH QUADRATIC PROGRAMMING SUBPROBLEMS 
C*                  基于二次优化子问题的非线性优化算法
C*                                                                     
C*                                                                     
C*  THIS SUBROUTINE SOLVES THE GENERAL NONLINEAR PROGRAMMING PROBLEM   
C*                                                                     
C*            MINIMIZE    F(X)                                         
C*                                                                     
C*            SUBJECT TO  C (X) .==. 0  ,  J = 1,...,MEQ               
C*                         J                                           
C*                                                                     
C*                        C (X) .>=. 0  ,  J = MEQ+1,...,M             
C*                         J                                           
C*                                                                     
C*                        XL .<=. X .<=. XU , I = 1,...,N.             
C*                          I      I       I                           
C*                                                                     
C*  THE ALGORITHM IMPLEMENTS THE METHOD OF HAN AND POWELL              
C*  WITH BFGS-UPDATE OF THE B-MATRIX AND L1-TEST FUNCTION              
C*  WITHIN THE STEPLENGTH ALGORITHM.                                   
C*  该算法在步长算法中利用b矩阵的bfgs-update和l1-test函数实现了han和Powell方法。
C*    PARAMETER DESCRIPTION:                                           
C*    ( * MEANS THIS PARAMETER WILL BE CHANGED DURING CALCULATION )    
C*                                                                     
C*    M              IS THE TOTAL NUMBER OF CONSTRAINTS, M .>=. 0     
                     M为约束的总数量
C*    MEQ            IS THE NUMBER OF EQUALITY CONSTRAINTS, MEQ .>=. 0 
                     等式约束的数量
C*    LA             SEE A, LA .>=. MAX(M,1)                           
C*    N              IS THE NUMBER OF VARIBLES, N .>=. 1  
                     变量的数量             
C*  * X()            X() STORES THE CURRENT ITERATE OF THE N VECTOR X  
C*                   ON ENTRY X() MUST BE INITIALIZED. ON EXIT X()     
C*                   STORES THE SOLUTION VECTOR X IF MODE = 0.
                     X()在X()必须初始化的条目上存储了n个向量X的当前迭代。 退出时，如果mode = 0，则x()存储解决方案向量x。           
C*    XL()           XL() STORES AN N VECTOR OF LOWER BOUNDS XL TO X.  
C*                   ELEMENTS MAY BE NAN TO INDICATE NO LOWER BOUND.   
C*    XU()           XU() STORES AN N VECTOR OF UPPER BOUNDS XU TO X.  
C*                   ELEMENTS MAY BE NAN TO INDICATE NO UPPER BOUND. 
                     上下边界
C*    F              IS THE VALUE OF THE OBJECTIVE FUNCTION. 
                     目标函数的值          
C*    C()            C() STORES THE M VECTOR C OF CONSTRAINTS,         
C*                   EQUALITY CONSTRAINTS (IF ANY) FIRST.              
C*                   DIMENSION OF C MUST BE GREATER OR EQUAL LA,       
C*                   which must be GREATER OR EQUAL MAX(1,M).   
                     C()首先存储约束的向量C，等式约束(如果有的话)。 C的维数必须大于或等于LA，而LA必须大于或等于MA(1,M)。         
C*    G()            G() STORES THE N VECTOR G OF PARTIALS OF THE      
C*                   OBJECTIVE FUNCTION; DIMENSION OF G MUST BE        
C*                   GREATER OR EQUAL N+1.
                     G()存储目标函数偏导数的n个向量G; g的维数必须大于或等于n +1。                               
C*    A(),LA,M,N     THE LA BY N + 1 ARRAY A() STORES                  
C*                   THE M BY N MATRIX A OF CONSTRAINT NORMALS.        
C*                   A() HAS FIRST DIMENSIONING PARAMETER LA,          
C*                   WHICH MUST BE GREATER OR EQUAL MAX(1,M).          
C*    F,C,G,A        MUST ALL BE SET BY THE USER BEFORE EACH CALL.     
C*  * ACC            ABS(ACC) CONTROLS THE FINAL ACCURACY.             
C*                   IF ACC .LT. ZERO AN EXACT LINESEARCH IS PERFORMED,
C*                   OTHERWISE AN ARMIJO-TYPE LINESEARCH IS USED.      
C*  * ITER           PRESCRIBES THE MAXIMUM NUMBER OF ITERATIONS.      
C*                   ON EXIT ITER INDICATES THE NUMBER OF ITERATIONS.  
C*  * MODE           MODE CONTROLS CALCULATION:                        
C*                   REVERSE COMMUNICATION IS USED IN THE SENSE THAT   
C*                   THE PROGRAM IS INITIALIZED BY MODE = 0; THEN IT IS
C*                   TO BE CALLED REPEATEDLY BY THE USER UNTIL A RETURN
C*                   WITH MODE .NE. IABS(1) TAKES PLACE.               
C*                   IF MODE = -1 GRADIENTS HAVE TO BE CALCULATED,     
C*                   WHILE WITH MODE = 1 FUNCTIONS HAVE TO BE CALCULATE
C*                   MODE MUST NOT BE CHANGED BETWEEN SUBSEQUENT CALLS 
C*                   OF SQP.                                           
C*                   EVALUATION MODES:                                 
C*        MODE = -1: GRADIENT EVALUATION, (G&A)                        
C*                0: ON ENTRY: INITIALIZATION, (F,G,C&A)               
C*                   ON EXIT : REQUIRED ACCURACY FOR SOLUTION OBTAINED 
C*                1: FUNCTION EVALUATION, (F&C)                        
C*                                                                     
C*                   FAILURE MODES:                                    
C*                2: NUMBER OF EQUALITY CONSTRAINTS LARGER THAN N      
C*                3: MORE THAN 3*N ITERATIONS IN LSQ SUBPROBLEM        
C*                4: INEQUALITY CONSTRAINTS INCOMPATIBLE               
C*                5: SINGULAR MATRIX E IN LSQ SUBPROBLEM               
C*                6: SINGULAR MATRIX C IN LSQ SUBPROBLEM               
C*                7: RANK-DEFICIENT EQUALITY CONSTRAINT SUBPROBLEM HFTI
C*                8: POSITIVE DIRECTIONAL DERIVATIVE FOR LINESEARCH    
C*                9: MORE THAN ITER ITERATIONS IN SQP                  
C*             >=10: WORKING SPACE W OR JW TOO SMALL,                  
C*                   W SHOULD BE ENLARGED TO L_W=MODE/1000             
C*                   JW SHOULD BE ENLARGED TO L_JW=MODE-1000*L_W       
C*  * W(), L_W       W() IS A ONE DIMENSIONAL WORKING SPACE,           
C*                   THE LENGTH L_W OF WHICH SHOULD BE AT LEAST        
C*                   (3*N1+M)*(N1+1)                        for LSQ    
C*                  +(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ         for LSI    
C*                  +(N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1       for LSEI   
C*                  + N1*N/2 + 2*M + 3*N + 3*N1 + 1         for SLSQPB 
C*                   with MINEQ = M - MEQ + 2*N1  &  N1 = N+1          
C*        NOTICE:    FOR PROPER DIMENSIONING OF W IT IS RECOMMENDED TO 
C*                   COPY THE FOLLOWING STATEMENTS INTO THE HEAD OF    
C*                   THE CALLING PROGRAM (AND REMOVE THE COMMENT C)    
c######################################################################
C     INTEGER LEN_W, LEN_JW, M, N, N1, MEQ, MINEQ
C     PARAMETER (M=... , MEQ=... , N=...  )
C     PARAMETER (N1= N+1, MINEQ= M-MEQ+N1+N1)
C     PARAMETER (LEN_W=
c    $           (3*N1+M)*(N1+1)
c    $          +(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ
c    $          +(N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1
c    $          +(N+1)*N/2 + 2*M + 3*N + 3*N1 + 1,
c    $           LEN_JW=MINEQ)
C     DOUBLE PRECISION W(LEN_W)
C     INTEGER          JW(LEN_JW)
c######################################################################
C*                   THE FIRST M+N+N*N1/2 ELEMENTS OF W MUST NOT BE    
C*                   CHANGED BETWEEN SUBSEQUENT CALLS OF SLSQP.        
C*                   ON RETURN W(1) ... W(M) CONTAIN THE MULTIPLIERS   
C*                   ASSOCIATED WITH THE GENERAL CONSTRAINTS, WHILE    
C*                   W(M+1) ... W(M+N(N+1)/2) STORE THE CHOLESKY FACTOR
C*                   L*D*L(T) OF THE APPROXIMATE HESSIAN OF THE        
C*                   LAGRANGIAN COLUMNWISE DENSE AS LOWER TRIANGULAR   
C*                   UNIT MATRIX L WITH D IN ITS 'DIAGONAL' and        
C*                   W(M+N(N+1)/2+N+2 ... W(M+N(N+1)/2+N+2+M+2N)       
C*                   CONTAIN THE MULTIPLIERS ASSOCIATED WITH ALL       
C*                   ALL CONSTRAINTS OF THE QUADRATIC PROGRAM FINDING  
C*                   THE SEARCH DIRECTION TO THE SOLUTION X*           
C*  * JW(), L_JW     JW() IS A ONE DIMENSIONAL INTEGER WORKING SPACE   
C*                   THE LENGTH L_JW OF WHICH SHOULD BE AT LEAST       
C*                   MINEQ                                             
C*                   with MINEQ = M - MEQ + 2*N1  &  N1 = N+1          
C*                                                                     
C*  THE USER HAS TO PROVIDE THE FOLLOWING SUBROUTINES:                 
C*     LDL(N,A,Z,SIG,W) :   UPDATE OF THE LDL'-FACTORIZATION.          
C*     LINMIN(A,B,F,TOL) :  LINESEARCH ALGORITHM IF EXACT = 1          
C*     LSQ(M,MEQ,LA,N,NC,C,D,A,B,XL,XU,X,LAMBDA,W,....) :              
C*                                                                     
C*        SOLUTION OF THE QUADRATIC PROGRAM                            
C*                QPSOL IS RECOMMENDED:                                
C*     PE GILL, W MURRAY, MA SAUNDERS, MH WRIGHT:                      
C*     USER'S GUIDE FOR SOL/QPSOL:                                     
C*     A FORTRAN PACKAGE FOR QUADRATIC PROGRAMMING,                    
C*     TECHNICAL REPORT SOL 83-7, JULY 1983                            
C*     DEPARTMENT OF OPERATIONS RESEARCH, STANFORD UNIVERSITY          
C*     STANFORD, CA 94305                                              
C*     QPSOL IS THE MOST ROBUST AND EFFICIENT QP-SOLVER                
C*     AS IT ALLOWS WARM STARTS WITH PROPER WORKING SETS               
C*                                                                     
C*     IF IT IS NOT AVAILABLE USE LSEI, A CONSTRAINT LINEAR LEAST      
C*     SQUARES SOLVER IMPLEMENTED USING THE SOFTWARE HFTI, LDP, NNLS   
C*     FROM C.L. LAWSON, R.J.HANSON: SOLVING LEAST SQUARES PROBLEMS,   
C*     PRENTICE HALL, ENGLEWOOD CLIFFS, 1974.                          
C*     LSEI COMES WITH THIS PACKAGE, together with all necessary SR's. 
C*                                                                     
C*     TOGETHER WITH A COUPLE OF SUBROUTINES FROM BLAS LEVEL 1         
C*                                                                     
C*     SQP IS HEAD SUBROUTINE FOR BODY SUBROUTINE SQPBDY               
C*     IN WHICH THE ALGORITHM HAS BEEN IMPLEMENTED.                    
C*                                                                     
C*  IMPLEMENTED BY: DIETER KRAFT, DFVLR OBERPFAFFENHOFEN               
C*  as described in Dieter Kraft: A Software Package for               
C*                                Sequential Quadratic Programming     
C*                                DFVLR-FB 88-28, 1988                 
C*  which should be referenced if the user publishes results of SLSQP  
C*                                                                     
C*  DATE:           APRIL - OCTOBER, 1981.                             
C*  STATUS:         DECEMBER, 31-ST, 1984.                             
C*  STATUS:         MARCH   , 21-ST, 1987, REVISED TO FORTRAN 77       
C*  STATUS:         MARCH   , 20-th, 1989, REVISED TO MS-FORTRAN       
C*  STATUS:         APRIL   , 14-th, 1989, HESSE   in-line coded       
C*  STATUS:         FEBRUARY, 28-th, 1991, FORTRAN/2 Version 1.04      
C*                                         accepts Statement Functions 
C*  STATUS:         MARCH   ,  1-st, 1991, tested with SALFORD         
C*                                         FTN77/386 COMPILER VERS 2.40
C*                                         in protected mode           
C*                                                                     
C**********************************************************************
C*                                                                     
C*  Copyright 1991: Dieter Kraft, FHM                                  
C*                                                                     
C**********************************************************************
```
