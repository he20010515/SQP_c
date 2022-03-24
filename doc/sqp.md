# Chapter 18 Sequential Quadratic Programming

本章考虑的所有方法都是活动集方法； 本章的一个更具描述性的标题可能是“非线性规划的活动集方法”。在第 14 章中，我们研究了非线性规划的内点方法，这是一种处理不等式约束问题的竞争方法。

有两种类型的活动集 SQP 方法。 在 IQP 方法中，在每次迭代中求解一般**不等式约束**的二次程序，具有计算步长和生成最优活动集估计的双重目标。  EQP 方法将这些计算解耦。 他们首先计算最优活动集的估计，然后求解等式约束的二次程序以找到步骤。 在本章中，我们研究 IQP 和 EQP 方法。

我们对 SQP 方法的开发分两个阶段进行。 首先，我们考虑激发 SQP 方法的局部方法，并允许我们在简单的设置中引入步进计算技术。 其次，我们考虑从远程起点实现收敛的实用线搜索和信任区域方法。 在整章中，我们都在考虑解决大问题的算法需求。

## 18.1 局部SQP方法

我们以一个等式约束问题开始:

$$
\min_x f(x) \\
s.t.\ \ c(x)=0
$$

其中$f:R^n \rightarrow R$ , $c:R^n\rightarrow R^m$是光滑的,SQP 方法背后的想法是通过二次规划子问题在当前迭代 xk 处对上述等式优化问题建模，然后使用该子问题的最小化器来定义新的迭代 xk+1。 挑战在于设计二次子问题，以便为非线性优化问题提供一个很好的步骤。 也许我们现在提出的 SQP 方法的最简单推导，将它们视为牛顿方法对 (18.1) 的 KKT 最优条件的应用。

上述问题的拉格朗日函数为$L(x,\lambda) = f(x)-\lambda ^Tc(x)$,我们用$A(x)$表示约束函数$c$的雅克比矩阵

$$
A(x)^{T}=\left[\nabla c_{1}(x), \nabla c_{2}(x), \ldots, \nabla c_{m}(x)\right]
$$

上述问题的一阶KKT条件可以表示为一个$m+n$阶的系统:

$$
F(x, \lambda)=\left[\begin{array}{c}
\nabla f(x)-A(x)^{T} \lambda \\
c(x)
\end{array}\right]=0
$$

使得$A(x^*)$满秩的任何解$(x^*,\lambda^*)$都满足上述系统.我们可以使用牛顿法求解上述非线性方程.

$$
\left[\begin{array}{cc}
\nabla_{x x}^{2} \mathcal{L}_{k} & -A_{k}^{T} \\
A_{k} & 0
\end{array}\right]\left[\begin{array}{c}
p_{k} \\
p_{\lambda}
\end{array}\right]=\left[\begin{array}{c}
-\nabla f_{k}+A_{k}^{T} \lambda_{k} \\
-c_{k}
\end{array}\right]
$$

**SQP框架**

从另一种角度看上述的迭代序列: 假设在迭代 (xk, λk) 时，我们使用二次规划对问题 (18.1) 建模
$$
\begin{aligned}
\min\ \  _{p} & f_{k}+\nabla f_{k}^{T} p+\frac{1}{2} p^{T} \nabla_{x x}^{2} \mathcal{L}_{k} p \\
\text { subject to }\ \ \ \  & A_{k} p+c_{k}=0 .
\end{aligned}
$$

如果假设 18.1 成立，则该问题有一个唯一解 (pk, lk)，满足

$$
\begin{aligned}
\nabla_{x x}^{2} \mathcal{L}_{k} p_{k}+\nabla f_{k}-A_{k}^{T} I_{k} &=0 \\
A_{k} p_{k}+c_{k} &=0
\end{aligned}
$$

向量 pk 和 lk 可以用牛顿方程 (18.6) 的解来识别。 如果我们从 (18.6) 中第一个方程的两边减去 AT k λk，我们得到