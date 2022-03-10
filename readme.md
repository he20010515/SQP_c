# 基于C语言的SQP算法库

## 目录

```bash
SQP_c
├─include // 项目的包含文件
├─other   // 开发过程中的参考资料/参考代码等
├─src     // 源码树
│  ├─function //有关函数的代码 [自动梯度][自动Hession矩阵计算]
│  ├─linarg   //线性代数代码   [矩阵][线性方程组]
│  ├─optimize //有关优化的代码 [二次优化][SQP]
│  └─util     //工具          [log][error]
└─test        //测试文件
```

项目正在开发中,整个项目由Cmake以及Ctest组织和构建.


TODOlist:

1. 线性方程组求解 done // 可能需要优化某些SOR方法不收敛的情况
2. 数值微分             
    1. 中心梯度        done
    2. 中心Hession矩阵 done
3. 二次优化
    1. 二次线性约束优化问题 done 线性方程组解法需要优化
    2. 二次不等式约束优化问题 ← HERE
    3. 二次混合约束优化问题 
4. SQP
