# 基于C语言的SQP算法库

## 目录

```bash
SQP_c
├─include // 项目的包含文件
├─other   // 开发过程中的参考资料/参考代码等
├─src     // 源码树
|  ├─log      //有关日志的代码
│  ├─function //有关函数的代码 [自动梯度][自动Hession矩阵计算]
│  ├─linarg   //线性代数代码   [矩阵][线性方程组]
│  ├─optimize //有关优化的代码 [二次优化][SQP]
│  └─util     //工具          [error]
└─test        //测试文件
```

项目正在开发中,整个项目由Cmake以及Ctest组织和构建.

## quick start

```bash
git clone https://github.com/he20010515/SQP_c.git
cd ./SQP_c
mkdir build
cd build
cmkae ..
cmake --build ,
ctest
```



## TODOlist:

1. 基础数据结构
   1. 矩阵
   2. 函数
2. 线性方程组求解 done // 可能需要优化某些SOR方法不收敛的情况
   1. 直接法 高斯消元法
   2. 迭代法 SOR
   3. 迭代法 
3. 数值微分             
    1. 中心梯度        done
    2. 中心Hession矩阵 done
4. 二次优化
    1. 二次线性约束优化问题 done
    2. 二次不等式约束优化问题 ←  done
    3. 二次混合约束优化问题   ←  done
5. SQP   ← HERE
   1. 分解为子问题 done
   2. 求解子问题 done
   3. BFGS update doing
   4. 验证结果,写测试
