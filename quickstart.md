# 快速使用指南

## 编译

```bash
git clone https://github.com/he20010515/SQP_c.git # 克隆项目代码
cd SQP_c        # 进入到项目目录
mkdir build     # 建立build目录
cd build        # 进入build目录
cmake ..        # 配置项目
cmkae --build . # 编译项目
```
项目配置时参考输出
```bash
[proc] Executing command: "D:\Program Files (x86)\VisualStudio\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.EXE" --no-warn-unused-cli -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE -DCMAKE_BUILD_TYPE:STRING=Debug -DCMAKE_C_COMPILER:FILEPATH=C:\mingw64\bin\x86_64-w64-mingw32-gcc.exe -DCMAKE_CXX_COMPILER:FILEPATH=C:\mingw64\bin\x86_64-w64-mingw32-g++.exe -Sd:/Workspace/HIT/SQP_c -Bd:/Workspace/HIT/SQP_c/build -G "Unix Makefiles"
[cmake] Not searching for unused variables given on the command line.
[cmake] -- The C compiler identification is GNU 8.1.0
[cmake] -- The CXX compiler identification is GNU 8.1.0
[cmake] -- Detecting C compiler ABI info
[cmake] -- Detecting C compiler ABI info - done
[cmake] -- Check for working C compiler: C:/mingw64/bin/x86_64-w64-mingw32-gcc.exe - skipped
[cmake] -- Detecting C compile features
[cmake] -- Detecting C compile features - done
[cmake] -- Detecting CXX compiler ABI info
[cmake] -- Detecting CXX compiler ABI info - done
[cmake] -- Check for working CXX compiler: C:/mingw64/bin/x86_64-w64-mingw32-g++.exe - skipped
[cmake] -- Detecting CXX compile features
[cmake] -- Detecting CXX compile features - done
[cmake] -- Found OpenMP_C: -fopenmp (found version "4.5") 
[cmake] -- Found OpenMP_CXX: -fopenmp (found version "4.5") 
[cmake] -- Found OpenMP: TRUE (found version "4.5")  
[cmake] -- found openmp
[cmake] -- Found OpenMP_C: -fopenmp (found version "4.5") 
[cmake] -- Found OpenMP_CXX: -fopenmp (found version "4.5") 
[cmake] -- Configuring done
[cmake] -- Generating done
[cmake] -- Build files have been written to: D:/Workspace/HIT/SQP_c/build
```
项目的编译时参考输出

```bash
[main] Building folder: SQP_c 
[build] Starting build
[proc] Executing command: "D:\Program Files (x86)\VisualStudio\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.EXE" --build d:/Workspace/HIT/SQP_c/build --config Debug --target all -j 14 --
[build] Consolidate compiler generated dependencies of target sqp
[build] [  1%] Building C object CMakeFiles/sqp.dir/src/function/function.c.obj
[build] [  2%] Building C object CMakeFiles/sqp.dir/src/linarg/sparse_,matrix.c.obj
[build] [  3%] Building C object CMakeFiles/sqp.dir/src/util/random.c.obj
[build] [  5%] Building C object CMakeFiles/sqp.dir/src/util/util.c.obj
[build] [  6%] Building C object CMakeFiles/sqp.dir/src/log/elog.c.obj
[build] D:/Workspace/HIT/SQP_c/src/util/util.c:21: warning: "LOG_TAG" redefined
[build]  #define LOG_TAG "util"
[build]  
[build] In file included from D:/Workspace/HIT/SQP_c/src/util/util.c:13:
[build] D:/Workspace/HIT/SQP_c/include/elog.h:216: note: this is the location of the previous definition
[build]      #define LOG_TAG          "NO_TAG"
[build]  
[build] [  7%] Building C object CMakeFiles/sqp.dir/src/log/elog_async.c.obj
[build] [  9%] Building C object CMakeFiles/sqp.dir/src/log/elog_buf.c.obj
[build] [ 10%] Building C object CMakeFiles/sqp.dir/src/log/elog_port.c.obj
[build] [ 11%] Building C object CMakeFiles/sqp.dir/src/log/elog_utils.c.obj
[build] [ 12%] Building C object CMakeFiles/sqp.dir/src/optimize/lp.c.obj
[build] [ 14%] Building C object CMakeFiles/sqp.dir/src/optimize/qp.c.obj
[build] D:/Workspace/HIT/SQP_c/src/optimize/lp.c:9: warning: "LOG_TAG" redefined
[build]  #define LOG_TAG "lp"
[build]  
[build] In file included from D:/Workspace/HIT/SQP_c/src/optimize/lp.c:7:
[build] D:/Workspace/HIT/SQP_c/include/elog.h:216: note: this is the location of the previous definition
[build]      #define LOG_TAG          "NO_TAG"
[build]  
[build] [ 15%] Building C object CMakeFiles/sqp.dir/src/optimize/simplex.c.obj
[build] D:/Workspace/HIT/SQP_c/src/optimize/qp.c:11: warning: "LOG_TAG" redefined
[build]  #define LOG_TAG "qp"
[build]  
[build] In file included from D:/Workspace/HIT/SQP_c/src/optimize/qp.c:9:
[build] D:/Workspace/HIT/SQP_c/include/elog.h:216: note: this is the location of the previous definition
[build]      #define LOG_TAG          "NO_TAG"
[build]  
[build] [ 16%] Building C object CMakeFiles/sqp.dir/src/optimize/sqp.c.obj
[build] D:/Workspace/HIT/SQP_c/src/optimize/simplex.c:7: warning: "LOG_TAG" redefined
[build]  #define LOG_TAG "simplex"
[build]  
[build] In file included from D:/Workspace/HIT/SQP_c/src/optimize/simplex.c:4:
[build] D:/Workspace/HIT/SQP_c/include/elog.h:216: note: this is the location of the previous definition
[build]      #define LOG_TAG          "NO_TAG"
[build]  
[build] D:/Workspace/HIT/SQP_c/src/optimize/sqp.c:9: warning: "LOG_TAG" redefined
[build]  #define LOG_TAG "SQP"
[build]  
[build] In file included from D:/Workspace/HIT/SQP_c/src/optimize/sqp.c:6:
[build] D:/Workspace/HIT/SQP_c/include/elog.h:216: note: this is the location of the previous definition
[build]      #define LOG_TAG          "NO_TAG"
[build]  
[build] [ 18%] Linking C static library libsqp.a
[build] [ 24%] Built target sqp
[build] [ 25%] Building C object test/CMakeFiles/test_function_function_with_openmp.dir/test_function_function_with_openmp.c.obj
[build] [ 27%] Building C object test/CMakeFiles/test_optimize_sqp2.dir/test_optimize_sqp2.c.obj
[build] [ 28%] Building C object test/CMakeFiles/test_function_function.dir/test_function_function.c.obj
[build] [ 29%] Building C object test/CMakeFiles/test_optimize_sqp.dir/test_optimize_sqp.c.obj
[build] [ 31%] Building C object test/CMakeFiles/test_linarg_linearequations4.dir/test_linarg_linearequations4.c.obj
[build] [ 32%] Building C object test/CMakeFiles/test_linarg_linearequations.dir/test_linarg_linearequations.c.obj
[build] [ 33%] Building C object test/CMakeFiles/test_optimize_sqp4.dir/test_optimize_sqp4.c.obj
[build] [ 35%] Building C object test/CMakeFiles/test_linarg_linearequations2.dir/test_linarg_linearequations2.c.obj
[build] [ 36%] Building C object test/CMakeFiles/test_linarg_linearequations3.dir/test_linarg_linearequations3.c.obj
[build] [ 37%] Building C object test/CMakeFiles/test_linarg_matrix2.dir/test_linarg_matrix2.c.obj
[build] [ 38%] Building C object test/CMakeFiles/test_linarg_matrix.dir/test_linarg_matrix.c.obj
[build] [ 40%] Building C object test/CMakeFiles/test_optimize_simplex6.dir/test_optimize_simplex6.c.obj
[build] [ 41%] Building C object test/CMakeFiles/test_openmp.dir/test_openmp.c.obj
[build] [ 42%] Building C object test/CMakeFiles/test_log_elog.dir/test_log_elog.c.obj
[build] [ 44%] Linking C executable test_function_function_with_openmp.exe
[build] [ 45%] Linking C executable test_optimize_sqp.exe
[build] [ 46%] Linking C executable test_function_function.exe
[build] [ 48%] Linking C executable test_linarg_linearequations.exe
[build] [ 49%] Linking C executable test_linarg_linearequations4.exe
[build] [ 50%] Linking C executable test_linarg_linearequations2.exe
[build] [ 51%] Linking C executable test_optimize_sqp2.exe
[build] [ 53%] Linking C executable test_optimize_simplex6.exe
[build] [ 54%] Linking C executable test_linarg_matrix.exe
[build] [ 55%] Linking C executable test_linarg_matrix2.exe
[build] [ 57%] Linking C executable test_openmp.exe
[build] [ 58%] Linking C executable test_optimize_sqp4.exe
[build] [ 59%] Linking C executable test_linarg_linearequations3.exe
[build] [ 61%] Linking C executable test_log_elog.exe
[build] [ 61%] Built target test_function_function_with_openmp
[build] [ 61%] Built target test_linarg_linearequations4
[build] [ 61%] Built target test_optimize_sqp
[build] [ 61%] Built target test_optimize_sqp2
[build] [ 61%] Built target test_function_function
[build] [ 61%] Built target test_optimize_sqp4
[build] [ 61%] Built target test_linarg_linearequations
[build] [ 61%] Built target test_linarg_linearequations2
[build] [ 61%] Built target test_linarg_linearequations3
[build] [ 61%] Built target test_linarg_matrix2
[build] [ 61%] Built target test_linarg_matrix
[build] [ 61%] Built target test_optimize_simplex6
[build] [ 61%] Built target test_log_elog
[build] [ 61%] Built target test_openmp
[build] [ 62%] Building C object test/CMakeFiles/test_optimize_lp.dir/test_optimize_lp.c.obj
[build] [ 63%] Building C object test/CMakeFiles/test_optimize_qp_active_set3.dir/test_optimize_qp_active_set3.c.obj
[build] [ 64%] Building C object test/CMakeFiles/test_optimize_qp_active_set.dir/test_optimize_qp_active_set.c.obj
[build] [ 66%] Building C object test/CMakeFiles/test_optimize_lp3.dir/test_optimize_lp3.c.obj
[build] [ 67%] Building C object test/CMakeFiles/test_optimize_qp_active_set2.dir/test_optimize_qp_active_set2.c.obj
[build] [ 68%] Building C object test/CMakeFiles/test_optimize_simplex.dir/test_optimize_simplex.c.obj
[build] [ 70%] Building C object test/CMakeFiles/test_optimize_qp.dir/test_optimize_qp.c.obj
[build] [ 71%] Building C object test/CMakeFiles/test_optimize_lp2.dir/test_optimize_lp2.c.obj
[build] [ 72%] Building C object test/CMakeFiles/test_optimize_simplex3.dir/test_optimize_simplex3.c.obj
[build] [ 74%] Building C object test/CMakeFiles/test_optimize_simplex2.dir/test_optimize_simplex2.c.obj
[build] [ 75%] Building C object test/CMakeFiles/test_optimize_simplex4.dir/test_optimize_simplex4.c.obj
[build] [ 76%] Building C object test/CMakeFiles/test_optimize_simplex5.dir/test_optimize_simplex5.c.obj
[build] [ 77%] Building C object test/CMakeFiles/test_optimize_sqp3.dir/test_optimize_sqp3.c.obj
[build] [ 79%] Building C object test/CMakeFiles/test_util_index_set.dir/test_util_index_set.c.obj
[build] [ 80%] Linking C executable test_optimize_lp.exe
[build] [ 84%] Linking C executable test_optimize_qp_active_set2.exe
[build] [ 83%] Linking C executable test_optimize_simplex.exe
[build] [ 84%] Linking C executable test_optimize_qp_active_set.exe
[build] [ 85%] Linking C executable test_optimize_qp.exe
[build] [ 87%] Linking C executable test_optimize_qp_active_set3.exe
[build] [ 88%] Linking C executable test_optimize_lp2.exe
[build] [ 90%] Linking C executable test_optimize_lp3.exe
[build] [ 89%] Linking C executable test_optimize_simplex5.exe
[build] [ 92%] Linking C executable test_optimize_simplex3.exe
[build] [ 93%] Linking C executable test_optimize_simplex2.exe
[build] [ 94%] Linking C executable test_optimize_simplex4.exe
[build] [ 96%] Linking C executable test_optimize_sqp3.exe
[build] [ 97%] Linking C executable test_util_index_set.exe
[build] [ 97%] Built target test_optimize_lp
[build] [ 97%] Built target test_optimize_qp_active_set
[build] [ 97%] Built target test_optimize_lp2
[build] [ 97%] Built target test_optimize_lp3
[build] [ 97%] Built target test_optimize_qp_active_set2
[build] [ 97%] Built target test_optimize_qp_active_set3
[build] [ 97%] Built target test_optimize_qp
[build] [ 97%] Built target test_optimize_simplex
[build] [ 97%] Built target test_optimize_simplex2
[build] [ 97%] Built target test_optimize_simplex3
[build] [ 97%] Built target test_optimize_simplex4
[build] [ 97%] Built target test_optimize_simplex5
[build] [ 97%] Built target test_optimize_sqp3
[build] [ 97%] Built target test_util_index_set
[build] [ 98%] Building C object test/CMakeFiles/test_util_random.dir/test_util_random.c.obj
[build] [100%] Linking C executable test_util_random.exe
[build] [100%] Built target test_util_random
[build] Build finished with exit code 0
```


## 运行测试



### 运行所有测试

```bash
cd build
ctest
```
参考输出
```bash
Test project D:/Workspace/HIT/SQP_c/build
      Start  1: test_function_function
 1/29 Test  #1: test_function_function ...............   Passed    0.03 sec
      Start  2: test_function_function_with_openmp
 2/29 Test  #2: test_function_function_with_openmp ...   Passed    0.03 sec
      Start  3: test_linarg_linearequations
 3/29 Test  #3: test_linarg_linearequations ..........   Passed    0.02 sec
      Start  4: test_linarg_linearequations2
 4/29 Test  #4: test_linarg_linearequations2 .........   Passed    0.02 sec
      Start  5: test_linarg_linearequations3
 5/29 Test  #5: test_linarg_linearequations3 .........   Passed    0.02 sec
      Start  6: test_linarg_linearequations4
 6/29 Test  #6: test_linarg_linearequations4 .........   Passed    0.02 sec
      Start  7: test_linarg_matrix
 7/29 Test  #7: test_linarg_matrix ...................   Passed    0.02 sec
      Start  8: test_linarg_matrix2
 8/29 Test  #8: test_linarg_matrix2 ..................   Passed    0.02 sec
      Start  9: test_log_elog
 9/29 Test  #9: test_log_elog ........................   Passed    0.01 sec
      Start 10: test_openmp
10/29 Test #10: test_openmp ..........................   Passed    0.02 sec
      Start 11: test_optimize_lp
11/29 Test #11: test_optimize_lp .....................   Passed    0.02 sec
      Start 12: test_optimize_lp2
12/29 Test #12: test_optimize_lp2 ....................   Passed    0.02 sec
      Start 13: test_optimize_lp3
13/29 Test #13: test_optimize_lp3 ....................   Passed    0.02 sec
      Start 14: test_optimize_qp
14/29 Test #14: test_optimize_qp .....................   Passed    0.03 sec
      Start 15: test_optimize_qp_active_set
15/29 Test #15: test_optimize_qp_active_set ..........   Passed    0.06 sec
      Start 16: test_optimize_qp_active_set2
16/29 Test #16: test_optimize_qp_active_set2 .........   Passed    0.02 sec
      Start 17: test_optimize_qp_active_set3
17/29 Test #17: test_optimize_qp_active_set3 .........   Passed    0.02 sec
      Start 18: test_optimize_simplex
18/29 Test #18: test_optimize_simplex ................   Passed    0.07 sec
      Start 19: test_optimize_simplex2
19/29 Test #19: test_optimize_simplex2 ...............   Passed    0.02 sec
      Start 20: test_optimize_simplex3
20/29 Test #20: test_optimize_simplex3 ...............   Passed    0.35 sec
      Start 21: test_optimize_simplex4
21/29 Test #21: test_optimize_simplex4 ...............   Passed    0.02 sec
      Start 22: test_optimize_simplex5
22/29 Test #22: test_optimize_simplex5 ...............   Passed    0.02 sec
      Start 23: test_optimize_simplex6
23/29 Test #23: test_optimize_simplex6 ...............   Passed    0.02 sec
      Start 24: test_optimize_sqp
24/29 Test #24: test_optimize_sqp ....................***Failed    0.22 sec
      Start 25: test_optimize_sqp2
25/29 Test #25: test_optimize_sqp2 ...................   Passed    0.06 sec
      Start 26: test_optimize_sqp3
26/29 Test #26: test_optimize_sqp3 ...................***Exception: SegFault  0.04 sec
      Start 27: test_optimize_sqp4
27/29 Test #27: test_optimize_sqp4 ...................***Exception: SegFault  0.10 sec
      Start 28: test_util_index_set
28/29 Test #28: test_util_index_set ..................   Passed    0.02 sec
      Start 29: test_util_random
29/29 Test #29: test_util_random .....................   Passed    0.02 sec

90% tests passed, 3 tests failed out of 29

Total Test time (real) =   1.46 sec

The following tests FAILED:
         24 - test_optimize_sqp (Failed)
         26 - test_optimize_sqp3 (SEGFAULT)
         27 - test_optimize_sqp4 (SEGFAULT)
Errors while running CTest
Output from these tests are in: D:/Workspace/HIT/SQP_c/build/Testing/Temporary/LastTest.log
Use "--rerun-failed --output-on-failure" to re-run the failed cases verbosely.
```


### 运行某个测试

例如测试test_optimize_simplex6

```bash
cd build
cd test
.\test_optimize_simplex6.exe
```

参考输出

```bash
I/elog            (D:/Workspace/HIT/SQP_c/src/log/elog.c:244 elog_start)EasyLogger V2.2.99 is initialize success.
I/simplex         (D:/Workspace/HIT/SQP_c/src/optimize/simplex.c:555 _linprog_simplex)simplex optimize complete
I/simplex         (D:/Workspace/HIT/SQP_c/src/optimize/simplex.c:556 _linprog_simplex)Optimization terminated successfully.
[2, 2, 0, 0, 0, 0, 4, 2, 2, ]
```
