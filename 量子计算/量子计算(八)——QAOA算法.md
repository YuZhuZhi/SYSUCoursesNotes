# 量子计算(八)——QAOA算法

量子近似优化算法（Quantum Approximate Optimization Algorithm，QAOA）是一种混合量子–经典架构的变分算法，即用经典优化器迭代训练一个含参量子线路，主要用于求解组合优化问题。所谓组合优化问题，典型的如旅行商问题、背包问题等，是在有限且离散的解空间中寻找最优解的优化问题，这一解空间会随着问题规模指数增长，因此一般是NP-hard的，可以给出数学模型为：
$$\begin{equation}
\begin{aligned}
    \min\quad &\mathcal{C}(x) \\
    \text{s.t.}\quad &\mathcal{U}(x)\geq 0 \\
    &x\in\mathbb{D}
\end{aligned}
\end{equation}$$
其中$\mathcal{C}(x)$是代价函数，$\mathcal{U}(x)$是约束条件，$\mathbb{D}$表示离散空间。在量子算法中，还会进一步要求$x$是二值的。

量子绝热算法（Quantum Adiabatic Algorithm，QAA）指出，通过设计哈密顿量的演化过程，量子系统可以从易制备的初始基态绝热演化至问题哈密顿量的基态，从而获得问题的一个近似解。而QAOA借鉴于这一物理图像，将连续的绝热演化过程，经由Trotter分解处理，离散化为多层的参数化量子门序列，从而在有限的电路深度下近似求解组合优化问题。

QAOA一般要求问题的数学模型要进一步约化到二次无约束二值优化（Quadratic Unconstrained Binary Optimization，QUBO）问题，即要求代价函数$\mathcal{C}(x)$是最高次项不超过2的多项式，变量$x$取值为0或1，并且没有任何约束条件。实际上，QAOA亦可直接应用于多项式无约束二值优化（Polynomial Unconstrained Binary Optimization，PUBO）问题。相较于仅包含二次项的QUBO问题，PUBO问题推广了交互项的阶数，能够编码包含高阶关联的复杂目标函数。本节将从QUBO和PUBO问题的数学模型入手，详细说明QAOA算法的构造步骤与运行原理。

## 一、伊辛模型，QUBO问题与PUBO问题

伊辛模型是统计物理中最具代表性的模型之一，用于研究自旋系统中的相变行为。伊辛模型定义在一个离散格点系综上，系统的总能量由**伊辛哈密顿量**给出：
$$
\begin{equation}
    \mathcal{H} = -J \sum_{\langle i,j \rangle} s_i s_j - h \sum_{i} s_i
\end{equation}
$$
其中，$\langle i,j \rangle$ 表示所有相邻格点对的集合；$J$ 为耦合常数，决定自旋之间的相互作用强度和类型；每个格点 $i$ 上的自旋变量 $s_i$ 只能取 $+1$（向上）或 $-1$（向下）两种状态。**自旋变量**$s$的取值与计算机中常用的0-1布尔变量$b$不相同，但可以使用$s=2b-1$相互转换。注意到式(2)中的次数为2，如果该次数能取更大的值，那么称模型为高阶伊辛模型。
在组合优化和量子计算中，许多 NP-hard 问题都可以转化为寻找等效伊辛模型的基态。而QAOA算法的思想正是将问题编码为哈密顿量，使用变分量子线路模拟绝热演变，当系统能量最低时就得到了问题的近似最优解。

QAOA算法由Edward Farhi在2014年提出时，就已经用在了有界度数图的最大割问题中。该问题是典型的QUBO问题，即其代价函数具有数学形式：
$$
\begin{equation}
    \mathcal{C}(\bm{x}) \coloneqq \bm{x}^\top \bm{Q} \bm{x}, \quad \bm{x} \in \{0,1\}^n,
\end{equation}
$$
其中问题中有$n$个变量，$\bm{x}$是一个具有$n$个元素的0-1向量，$\bm{Q}$是$n\times n$的实数矩阵。这显然是一个二次型，但由于对任意二进制变量$b$都有$b^2=b$，因此$\bm{Q}$的对角线元素在QUBO问题中可以特指一次项系数。有时候式(3)会被写成求和的形式：
$$
\begin{equation}
    \mathcal{C}(x) = \sum_{\substack{i,j=1\\i\leq j}}^{n} Q'_{ij} x_i x_j + \sum_{i=1}^{n} c_i x_i,
\end{equation}
$$
从而更加方便地转换为伊辛模型。需要注意：
+ 式(4)中的$Q'_{ij}$不一定等于式(3)中二次型矩阵$\bm{Q}$的第$i$行第$j$列元素$Q_{ij}$，而应视$\bm{Q}$的形式决定。例如，若$\bm{Q}$是对称矩阵，那么$Q'_{ij}=2Q_{ij}$。
+ $Q'_{ij}$允许$i=j$，这似乎会造成与一次项求和部分的冗余。但允许$i=j$会使之后转换伊辛模型的计算更加简单。事实上，可以令$Q'_{ii}=0,c_i=Q_{ii}$以消除冗余，此时求和号下$i\leq j$改为$i<j$不影响结果。

为了将代价函数转换为伊辛模型，只需将布尔变量 $x_i$ 转换为自旋变量 $s_i$。通过代入变换关系$x_i = \frac{s_i + 1}{2}$，即可得到伊辛哈密顿量：
$$
\begin{align}
    \mathcal{H}(s) &= \sum_{\substack{i,j=1\\i\leq j}}^{n} \frac{Q'_{ij}}{4}(s_i + 1)(s_j + 1) + \sum_{i=1}^{n} \frac{c_i}{2}(s_i + 1) \\
    &= \sum_{\substack{i,j=1\\i< j}}^{n} \frac{Q'_{ij}}{4} s_i s_j + \sum_{i=1}^{n} \frac{1}{2}\left( c_i + \sum_{j=1}^{n}Q'_{ij} \right) s_i + \bcancel{\sum_{i=1}^{n} \frac{1}{2}\left( c_i + \sum_{j=1}^{n}\frac{Q'_{ij}}{2} \right)}.
\end{align}
$$
由于式(6)的最后一个求和项为常数，不影响优化问题的解，因此常常略去。令$J_{ij} = -\frac{Q'_{ij}}{4}, h_i = -\frac{1}{2}\left( c_i + \sum_{j=1}^{n}Q'_{ij} \right)$，就成功转化为伊辛哈密顿量：
$$
\begin{equation}
    \mathcal{H}(s) = -\sum_{i,j=1}^{n} J_{ij} s_i s_j - \sum_{i=1}^{n} h_i s_i .
\end{equation}
$$

当实际问题具有约束的时候，一般会通过添加惩罚项的方式将问题转换为QUBO问题。例如式(1)所示的原始问题，将会被转化为：
$$
\begin{equation}
    \min \quad \mathcal{C}(x) + \lambda\mathcal{P}(x),
\end{equation}
$$
其中$\mathcal{P}(x)$是约束条件$\mathcal{U}(x)$对应的惩罚函数；$\lambda>0$是任取的惩罚系数，一般不应取太小。下表给出了一些典型约束条件应转换到的惩罚函数。

| 序号 | 约束条件 | 惩罚函数 |
|:---:|:---|:---|
| 1 | $x_1 + x_2 \leq 1$ | $x_1 x_2$ |
| 2 | $x_1 + x_2 \geq 1$ | $1 - x_1 - x_2 + x_1 x_2$ |
| 3 | $x_1 + x_2 = 1$ | $1 - x_1 - x_2 + 2x_1 x_2$ |
| 4 | $x_1 = x_2$ | $x_1 + x_2 - 2x_1 x_2$ |
| 5 | $x_1 \leq x_2$ | $x_1 - x_1 x_2$ |
| 6 | $x_1 + x_2 + x_3 \leq 1$ | $x_1x_2 + x_1x_3 + x_2x_3$ |

以QUBO问题为基础，若不限定代价函数中最高次项的次数，则问题扩展为PUBO问题。此时代价函数可以写为：
$$
\begin{equation}
\begin{aligned}
    \mathcal{C}(x) &= \sum_{i_1}T_{i_1}x_{i_1} + \sum_{i_1 < i_2}T_{i_1i_2}x_{i_1}x_{i_2} + \sum_{i_1 < i_2 < i_3}T_{i_1i_2i_3}x_{i_1}x_{i_2}x_{i_3} + \dots \\
    &= \sum_{\delta=(\delta_1,\dots,\delta_n) \in \{0,1\}^n} T_{\delta}\prod_{i=1}^{n}x_i^{\delta_i},
\end{aligned}
\end{equation}
$$
其中$n$是变量的数量。通过变量代换、添加约束并转换为惩罚项的方式，PUBO问题总能被归约为QUBO问题。相关算法在文献中有详尽的总结。利用变换式$s=2x-1$能将PUBO问题转换为高阶伊辛模型。

## 二、将问题编码为厄米矩阵

在上节中已提到，使用自旋变量 $s$ 与布尔变量 $x$的变量代换$s=2x-1$即可将QUBO问题或PUBO问题转换为伊辛模型（包括高阶伊辛模型）。要应用QAOA算法，下一步就是将伊辛模型编码为厄米矩阵。现假设已经得到的伊辛哈密顿量为：
$$
\begin{equation}
    \mathcal{H}(s) = \sum_{\delta=(\delta_1,\dots,\delta_n) \in \{0,1\}^n} h_{\delta}\prod_{i=1}^{n}s_i^{\delta_i} ,
\end{equation}
$$
那么只需做一点改动，即可得到对应的厄米矩阵形式哈密顿量，称为**问题哈密顿量**：
$$
\begin{equation}
    \bm{\mathcal{H}}_{\mathcal{C}} = \sum_{\delta=(\delta_1,\dots,\delta_n) \in \{0,1\}^n} h_{\delta}\bigotimes_{i=1}^{n}\bm{Z}_i^{\delta_i},
\end{equation}
$$
其中$\bm{Z}$即泡利Z矩阵，且$\bm{Z}^0=\bm{Z}^2=\bm{I}$、$\bm{Z}^1=\bm{Z}$，$\bm{Z}_i$等价于将在第$i$量子比特上作用Z门。之所以这样替换，是因为可以自然地将物理可观测量——自旋沿 $z$ 轴的分量，即泡利 Z 矩阵与自旋变量 $s$ 对应起来。针对问题哈密顿量$\bm{\mathcal{H}}_{\mathcal{C}}$，有以下定理：

> **定理**  
> 问题哈密顿量$\bm{\mathcal{H}}_{\mathcal{C}}$的特征向量是从$\ket{0}$到$\ket{2^n-1}$的计算基态，且每个特征向量$\ket{x}$对应的特征值是代价函数在相应解$x$下的函数值$\mathcal{C}(x)$（在不忽略常数项的情况下）。

> **证明：**  
> 由于：
> $$
> \bm{I} = \begin{bmatrix}
>     1 & 0 \\ 0 & 1
> \end{bmatrix}, \quad \bm{Z} = \begin{bmatrix}
>     1 & 0 \\ 0 & -1
> \end{bmatrix}
> $$
> 都是对角矩阵，$\bigotimes_{i=1}^{n}\bm{Z}_i^{\delta_i}$自然也是对角矩阵，相应地求累加后$\bm{\mathcal{H}}_{\mathcal{C}}$还是对角的，故$\bm{\mathcal{H}}_{\mathcal{C}}$的对角元素就是特征值。同时，由于$\bm{I},\bm{Z}$在计算基下可选$\{ \ket{0},\ket{1} \}$作为共同特征向量，故任意矩阵$\bm{U}\in\{ \bm{I},\bm{Z} \}^{\otimes n}$都可选$\{ \ket{0},\ket{1} \}^{\otimes n}$即所有计算基态作为特征向量，进一步地$\bm{\mathcal{H}}_{\mathcal{C}}$作为$\{ \bm{I},\bm{Z} \}^{\otimes n}$中矩阵的实系数线性组合也能以所有计算基态作为特征向量。
> 另一方面，将特征向量$\ket{x}$写为二进制形式$\ket{x_1\dots x_n}$并令$s_i=2x_i-1$。一个显然的结果是对于任意$\delta\in\{0,1\}^n$有：
> $$
> \left( \bigotimes_{i=1}^{n}\bm{Z}_i^{\delta_i} \right) \ket{x} = \left( \prod_{i=1}^{n}s_i^{\delta_i} \right) \ket{x},
> $$
> 因为每有一个$x_i=1$且$\delta_i=1$就会使$\ket{x}$乘一次$-1$。那么：
> $$
> \bm{\mathcal{H}}_{\mathcal{C}} \ket{x} = \left( \sum_{\delta} h_{\delta}\bigotimes_{i=1}^{n}\bm{Z}_i^{\delta_i} \right) \ket{x}
> = \left( \sum_{\delta} h_{\delta}\prod_{i=1}^{n}s_i^{\delta_i} \right) \ket{x} = \mathcal{H}(s) \ket{x}.
> $$
> 在不忽略常数项的情况下，$\mathcal{H}(s) = \mathcal{C}(x)$，所以$\bm{\mathcal{H}}_{\mathcal{C}} \ket{x} = \mathcal{C}(x) \ket{x}$，故$\ket{x}$对应的特征值就是相应的代价$\mathcal{C}(x)$。
> **证毕.**

显然一种求解此组合优化问题的方法就是从$\bm{\mathcal{H}}_{\mathcal{C}}$中找到最小的特征值，再以此求出对应特征向量。但必须指出，$\bm{\mathcal{H}}_{\mathcal{C}}$是$2^n$阶的，因此该方法并没有摆脱指数增长带来的困难。然而另一方面，对于一个哈密顿量，其最小特征值对应的特征向量正是量子系统的基态。因此，量子绝热算法基于量子绝热定理，通过设计一条演化路径，将一个易于制备基态的哈密顿量连续地演化为问题哈密顿量。若演化足够缓慢且满足绝热条件，最终系统将处于问题哈密顿量的基态，此时对该态进行测量，即可得到优化问题的解。而QAOA借鉴量子绝热算法，以参数化量子电路来分段近似模拟绝热演化。

> **定理（量子绝热定理）**  
> 对于一个缓慢变化的、无简并的哈密顿量，如果系统初始处于其基态或某个瞬时能级的本征态，则系统在随时间演化的过程中将始终保持在相应的瞬时能级态上，而不会跃迁到其他能级。 

## 三、分段近似模拟绝热演化

接下来首先给出QAOA的整体结构与电路图，之后再分模块介绍。QAOA的算法过程可以用下式完整表示：
$$
\begin{equation}
    \ket{\bm{\gamma}, \bm{\beta}} = \underbrace{U(\bm{\mathcal{H}}_{B}, \beta_p)U(\bm{\mathcal{H}}_{\mathcal{C}}, \gamma_p)}_{\text{第$p$层}} \cdots \underbrace{U(\bm{\mathcal{H}}_{B}, \beta_2)U(\bm{\mathcal{H}}_{\mathcal{C}}, \gamma_2)}_{\text{第2层}} \underbrace{U(\bm{\mathcal{H}}_{B}, \beta_1)U(\bm{\mathcal{H}}_{\mathcal{C}}, \gamma_1)}_{\text{第1层}} \ket{s},
\end{equation}
$$
其对应量子电路如图\ref{fig:chap9-QAOAGlobal}所示。其中，$\bm{\mathcal{H}}_{\mathcal{C}}$即前文所述的问题哈密顿量，$\bm{\mathcal{H}}_{B}$是**混合哈密顿量**，定义为：
$$
\begin{equation}
    \bm{\mathcal{H}}_{B} = \sum_{i=1}^{n} \bm{X}_i,
\end{equation}
$$
$\bm{X}_i$等价于将在第$i$量子比特上作用X门。$\bm{\gamma} = (\gamma_1,\dots,\gamma_p),\, \bm{\beta} = (\beta_1,\dots,\beta_p)$分别是包含$p$个角度参数的参数向量，称参数$p$为**深度**。由图\ref{fig:chap9-QAOAGlobal}可知$\ket{s}$正是初态$\ket{0}^{\otimes n}$被作用$H^{\otimes n}$之后得到的均匀叠加态$\ket{+}^{\otimes n}$。而：
$$
\begin{equation}
    U(\bm{\mathcal{H}}, \theta) = \mathrm{e}^{-\mathrm{i}\theta \bm{\mathcal{H}}}
\end{equation}
$$
是一个酉算子，接受一个矩阵和角度作为输入，输出一个酉矩阵。简便起见，将任意一层记为酉算子$L_i(\gamma_i, \beta_i) = U(\bm{\mathcal{H}}_{B}, \beta_i)U(\bm{\mathcal{H}}_{\mathcal{C}}, \gamma_i)$，所有层记为酉算子$L(\bm{\gamma}, \bm{\beta}) = \prod_{i=p}^{1} L_i(\gamma_i, \beta_i)$。

\begin{figure}[htbp]
\centering
\begin{quantikz}[row sep = 0.4cm, column sep = 0.4cm]
    \lstick{$q_1 = \ket{0}$} & \gate{H}\slice{$\ket{s}$} & \gate[4]{U(\bm{\mathcal{H}}_{\mathcal{C}}, \gamma_1)} & \gate[4]{U(\bm{\mathcal{H}}_{B}, \beta_1)} & \ \ldots\  & \gate[4]{U(\bm{\mathcal{H}}_{\mathcal{C}}, \gamma_p)} & \gate[4]{U(\bm{\mathcal{H}}_{B}, \beta_p)}\slice{$\ket{\gamma, \beta}$} & \meter{} \\
    % \lstick{\vdots} \wave &&&&&&&  \\
    \lstick{$q_2 = \ket{0}$} & \gate{H} & & & \ \ldots\  & & & \meter{} \\
    \lstick{\vdots} \wave &&&&&&& \\
    \lstick{$q_n = \ket{0}$} & \gate{H} & & & \ \ldots\  & & & \meter{}
\end{quantikz}
\caption{QAOA量子电路示意}
\label{fig:chap9-QAOAGlobal}
\end{figure}

若在式(14)中代入$\bm{\mathcal{H}}_{B}$与$\beta$，那么可做如下计算：
$$
\begin{align}
    U(\bm{\mathcal{H}}_{B}, \beta) &= \mathrm{e}^{-\mathrm{i}\beta \sum_{i=1}^{n} \bm{X}_i} \\
    &= \prod_{i=1}^{n} \mathrm{e}^{-\mathrm{i}\beta \bm{X}_i} = \bigotimes_{i=1}^{n} \mathrm{e}^{-\mathrm{i}\beta \bm{X}} \\
    &= \bigotimes_{i=1}^{n} (\cos(\beta)\bm{I} - \mathrm{i}\sin(\beta)\bm{X}) \\
    &= \bigotimes_{i=1}^{n} RX_i(2\beta),
\end{align}
$$
这说明混合酉算子$U(\bm{\mathcal{H}}_{B}, \beta)$的作用就是在每一个量子比特上都作用一个$RX(2\beta)$门，如图\ref{fig:chap9-MixerUnitaryOperatorEffect}所示。$RX(\theta)$门是：
$$
\begin{equation}
    RX(\theta) = \begin{bmatrix}
        \cos\frac{\theta}{2} & -\mathrm{i}\sin\frac{\theta}{2} \\ -\mathrm{i}\sin\frac{\theta}{2} & \cos\frac{\theta}{2}
    \end{bmatrix}
\end{equation}
$$
因此，在量子电路图\ref{fig:chap9-QAOAGlobal}中，唯一会随着问题而变的结构是问题酉算子$U(\bm{\mathcal{H}}_{\mathcal{C}}, \gamma_*)$，故只要研究清楚问题哈密顿量$\bm{\mathcal{H}}_{\mathcal{C}}$会如何构建量子电路，QAOA的电路构建就手到擒来了。

\begin{figure}[htbp]
\centering
\begin{quantikz}[row sep = 0.7cm, align equals at=2.5]
    & \gate[4, disable auto height]{\raisebox{1.5em}{$U(\bm{\mathcal{H}}_{B}, \beta)$}} & \\
    & & \\
    \wave && \\
    & & 
\end{quantikz} \quad $=$ \quad \begin{quantikz}[row sep = 0.2cm, align equals at=2.5]
    & \gate{RX(2\beta)} & \\
    & \gate{RX(2\beta)} & \\
    & \midstick{$\vdots$} & \ghost{RX(2\beta)} \\
    & \gate{RX(2\beta)} & 
\end{quantikz}
\caption{混合酉算子$U(\bm{\mathcal{H}}_{B}, \beta)$的作用}
\label{fig:chap9-MixerUnitaryOperatorEffect}
\end{figure}

观察式~\eqref{eq:chap9-QuestionHamiltonian}，注意到$\{ \bm{I}, \bm{Z} \}^{\otimes n}$中的矩阵是两两对易的，因此在式(14)中代入问题哈密顿量$\bm{\mathcal{H}}_{\mathcal{C}}$与其对应角度参数$\gamma$即得：
$$
\begin{align}
    U(\bm{\mathcal{H}}_{\mathcal{C}}, \gamma) &= \mathrm{e}^{-\mathrm{i}\gamma \sum_{\delta} h_{\delta}\bigotimes_{i=1}^{n}\bm{Z}_i^{\delta_i}} \\
    &= \prod_{\delta} \mathrm{e}^{-\mathrm{i} (h_\delta\gamma) \bigotimes_{i=1}^{n}\bm{Z}_i^{\delta_i}}.
\end{align}
$$
显然可以认为系数是角度参数的一部分，之后再回过头来处理。上式主要透露出的信息是：只需要考虑$\mathrm{e}^{-\mathrm{i} \gamma \bigotimes_{i=1}^{n}\bm{Z}_i^{\delta_i}}$最终能转换成什么电路，就解决了QAOA整体电路构建的问题。现给出以下定理：

> **定理 ($Z_{i_1}$与$Z_{i_1}\otimes Z_{i_2}$对应的酉算子及其量子电路)：**  
> 若算子$\bm{\mathcal{H}} = \bm{Z}_{i_1}$，那么其对应酉算子是$U(\bm{\mathcal{H}}, \gamma) = \mathrm{e}^{-\mathrm{i} \gamma \bm{Z}_{i_1}} = RZ_{i_1}(2\gamma)$，即在第$i_1$量子比特上作用$RZ(2\gamma)$门，如下左图所示。$RZ(\theta)$门是：
> $$
> RZ(\theta) = \begin{bmatrix}
>     \mathrm{e}^{-\mathrm{i}\frac{\theta}{2}} & 0 \\ 0 & \mathrm{e}^{\mathrm{i}\frac{\theta}{2}}
> \end{bmatrix}
> $$
> 若算子$\bm{\mathcal{H}} = \bm{Z}_{i_1}\otimes Z_{i_2}$，那么其对应酉算子是$\mathrm{e}^{-\mathrm{i} \gamma \bm{Z}_{i_1}\otimes Z_{i_2}} = \text{CNOT}_{i_2}^{i_1}\cdot RZ_{i_2}(2\gamma)\cdot \text{CNOT}_{i_2}^{i_1}$，其量子电路如下右图所示，$\text{CNOT}_{i_2}^{i_1}$表示以第$i_1$量子比特为控制位、第$i_2$量子比特为受控位。
> **电路示意：**
> 左图（单比特 $RZ$）：  
> ```
> q_i1: ── RZ(2γ) ──
> ```
> 右图（双比特 $ZZ$ 演化）：  
> ```
> q_i1: ── ● ──────────── ● ──
>          │               │
> q_i2: ── ⊕ ── RZ(2γ) ── ⊕ ──
> ```
> （其中 ● 表示控制位，⊕ 表示受控非门目标位）
> **图：** $Z_{i_1}$与$Z_{i_1}\otimes Z_{i_2}$对应的量子电路

> **证明：**  
> 首先当$\bm{\mathcal{H}} = \bm{Z}$时做如下计算：
> $$
> \begin{aligned}
> U(\bm{\mathcal{H}}, \gamma) &= \mathrm{e}^{-\mathrm{i} \gamma \bm{Z}} \\
> &= \cos(\gamma)\bm{I} - \mathrm{i}\sin(\gamma)\bm{Z} \\
> &= \operatorname{diag}\{\mathrm{e}^{-\mathrm{i}\gamma}, \mathrm{e}^{\mathrm{i}\gamma}\} \\
> &= RZ(2\gamma) 
> \end{aligned}
> $$
> 在考虑将要作用到的量子比特后，易证定理的第一部分。而当$\bm{\mathcal{H}} = \bm{Z}\otimes Z$时，计算如下：
> $$
> \begin{aligned}
> U(\bm{\mathcal{H}}, \gamma) &= \mathrm{e}^{-\mathrm{i} \gamma \bm{Z}\otimes \bm{Z}} \\
> &= \operatorname{diag}\{ -\mathrm{e}^{\mathrm{i}\gamma}, \mathrm{e}^{\mathrm{i}\gamma}, \mathrm{e}^{\mathrm{i}\gamma}, \mathrm{e}^{-\mathrm{i}\gamma} \} \\
> &= \text{CNOT}\cdot RZ(2\gamma)\cdot \text{CNOT}
> \end{aligned}
> $$
> 在考虑将要作用到的量子比特后定理第二部分得证。 ∎

必须指出，图\ref{fig:chap9-SingleCoupleFoldZGate->QuantumCircuit}中的$q_{i_1}$与$q_{i_2}$并不一定是相邻的，$i_1,i_2$的取值只取决于式\eqref{eq:chap9-QuestionHamiltonian}中的$\delta$。由此定理得到以下推论：

> **推论 ($\bigotimes_{k=1}^{n} Z_{i_k}$对应的酉算子及其量子电路)：**  
> 若算子$\bm{\mathcal{H}} = \bigotimes_{k=1}^{n} Z_{i_k}$，那么其对应酉算子$U(\bm{\mathcal{H}}, \gamma)$的量子电路如下式及下图所示。
> $$
> \left( \prod_{k=1}^{n-1}\text{CNOT}_{i_n}^{i_k} \right)\cdot RZ_{i_n}(2\gamma)\cdot \left( \prod_{k=1}^{n-1}\text{CNOT}_{i_n}^{i_{n-k}} \right)
> $$
> **电路示意（$n$比特情况）：**
> ```
> q_i1:  ── ● ──────────────────────────── ● ──
>           │                               │
> q_i2:  ── ┼── ● ────────────────────── ● ──┼──
>           │    │                         │  │
> q_i3:  ── ┼── ┼── ● ──────────────── ● ──┼──┼──
>           │    │    │                 │   │  │
>          ...  ...  ...  ...  ...  ... ... ... ...
>           │    │    │                 │   │  │
> q_i{n-1}: ── ┼── ┼── ┼── ● ──────── ● ──┼──┼──┼──
>           │    │    │    │         │   │  │  │
> q_i{n}:  ── ⊕── ⊕── ⊕── ⊕── RZ(2γ) ── ⊕── ⊕── ⊕── ⊕──
> ```
> 其中：
> - 最底部线路为$q_{i_n}$（受控目标比特）
> - 上方各线路为控制比特$q_{i_1}, q_{i_2}, \dots, q_{i_{n-1}}$
> - ● 表示控制位，⊕ 表示CNOT的目标位
> - 电路结构：先将所有上方比特作为控制位，依次对$q_{i_n}$作用CNOT门；然后在$q_{i_n}$上作用$RZ(2\gamma)$门；最后再用相同的CNOT序列（顺序相反）进行还原。
> **图：** $\bigotimes_{k=1}^{n} Z_{i_k}$对应的量子电路

> **证明：**  
> 证明留给读者作为习题。提示：使用数学归纳法。 ∎

由定理\ref{theo:chap9-SingleCoupleFoldZGate->QuantumCircuit}，任意QUBO问题都能轻易地转化为QAOA量子电路了；进一步地，由推论\ref{coro:chap9-MultiFoldZGate->QuantumCircuit}，任意PUBO问题也能转化为量子电路。但千万别忘了：QAOA的终态$\ket{\bm{\gamma}, \bm{\beta}}$是由$\bm{\gamma}, \bm{\beta}$中共$2p$个角度参数决定的，这些参数直接影响了算法对量子绝热过程的近似模拟的好坏。因此QAOA的下一步，是寻找到一组足够好的参数、以尽可能地模拟真实的量子绝热过程，或者说输出足够好的结果。

