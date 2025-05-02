# Qubitization

## 一、基础

### 标准形式编码

现在将量子系统看做一个信号系统，不再仔细思考系统中的每一个量子位是如何产生作用的。我们只关心系统整体是如何运作的，或者系统内部的子系统如何产生耦合。或者说，现在系统的子系统会被当做一个量子位对待，这就是量子位化的由来。

假设系统的初态是$| \psi \rangle_s \in \mathcal{H}_s$——或者说，将之认为是系统的输入信号。那么显而易见的是，假如系统对信号的操作算符是$\hat{H}$，那么系统的输出是$\hat{H}| \psi \rangle_s$。然而在$\hat{H}$不是酉矩阵的情况下，我们就需要考虑将其嵌入到一个更大的希尔伯特空间中。

假设嵌入到空间$\mathcal{H}_a \otimes \mathcal{H}_s$，嵌入后得到的酉矩阵操作算符是$\hat{U}$，称$\mathcal{H}_a$为辅助空间而$\mathcal{H}_s$为系统空间。假设辅助态$| G \rangle_a \in \mathcal{H}_a$，显然可得：

$$\begin{equation}
    \hat{U} | G \rangle_a | \psi \rangle_s = | G \rangle_a \hat{H}| \psi \rangle_s + \sqrt{1-\Vert\hat{H}| \psi \rangle_s\Vert^2}| G_\psi^\perp \rangle_{as}
\end{equation}$$

其中，$| G^\perp \rangle_{as}$是$| G \rangle_a \otimes \hat{\mathcal{I}}_s$的正交补空间，而$| G_\psi^\perp \rangle_{as}$则是在上述运算后得到的属于此空间的一个分量，也就是存在：

$$\begin{equation}
    (\langle G |_a \otimes \hat{\mathcal{I}}_s)| G_\psi^\perp \rangle_{as}=0
\end{equation}$$

因此，将$(1)$式左边乘上$\langle G |_a \otimes \hat{\mathcal{I}}_s$，就可以消除右边加号后的第二项，得到：

$$\begin{equation}
    \begin{aligned}
        (\langle G |_a \otimes \hat{\mathcal{I}}_s)\hat{U} (| G \rangle_a | \psi \rangle_s) &= (\langle G |_a \otimes \hat{\mathcal{I}}_s)(| G \rangle_a \hat{H}| \psi \rangle_s) \\
        &=\langle G | G \rangle_a \cdot \hat{H}| \psi \rangle_s = \hat{H}| \psi \rangle_s
    \end{aligned}
\end{equation}$$

将$(3)$式两边同时消掉$| \psi \rangle_s$，即因为$| G \rangle_a | \psi \rangle_s=(| G \rangle_a \otimes \hat{I}_s)| \psi \rangle_s$，从而我们得到从$\hat{U}$中提取出$\hat{H}$的方式：

$$\begin{equation}
    \hat{H} = (\langle G |_a \otimes \hat{\mathcal{I}}_s)\hat{U} (| G \rangle_a \otimes \hat{I}_s)
\end{equation}$$

>$\hat{U}$的形式是：$$\hat{U} = \begin{bmatrix}
    \hat{H} & \cdot \\
    \cdot & \cdot \\
\end{bmatrix}$$

### 作用

那么，$(1)$式有什么作用？这里就是局部测量的一种应用。测量辅助空间，如果得到的结果是$| G \rangle_a$，那么可以断定成功地将$\hat{H}$作用到$| \psi \rangle_s$上。这一成功的概率是$\Vert\hat{H}| \psi \rangle_s\Vert^2$，整个系统的输出即$\hat{H}| \psi \rangle_s$被投影到$\dfrac{| G \rangle_a \hat{H}| \psi \rangle_s}{\Vert\hat{H}| \psi \rangle_s\Vert}$。

定义$\mathcal{H}_G=| G \rangle \otimes \mathcal{H}_s$，在明显可以辨别的情况下将$| G \rangle_a \otimes \hat{I}_s$简写为$| G \rangle$。再定义$| G_\psi \rangle = | G \rangle | \psi \rangle \in \mathcal{H}_G$(补充上一节的定义)。

>必须指出：虽然我们对子空间$\mathcal{H}_{G^\perp}$并不感兴趣，尽管$\hat{U}$是作用在$\mathcal{H}_G$空间上，但$\hat{U}$中未定义的部分会使$\hat{U}| G_\psi \rangle$中的一部分泄漏到$\mathcal{H}_{G^\perp}$中。这也意味着$\hat{U}| G_\psi^\perp \rangle$的结果是未定义的。

### 当$\hat{H}$是埃尔米特矩阵

埃尔米特矩阵是指$H=H^\dagger$，也就是复对称矩阵。类似于实对称矩阵，埃尔米特矩阵的特征值也都是实数。假设$\hat{H}$的特征值为$\lambda$，特征向量为$| \lambda \rangle$，第一显然的是$\hat{H}| \lambda \rangle = \lambda | \lambda \rangle$，于是改写$(1)$式，令其中的$| \psi \rangle$是$\hat{H}$的特征向量$| \lambda \rangle$，那么就可以得到：

$$\begin{equation}
    \hat{U} | G \rangle | \lambda \rangle = \hat{U}| G_\lambda \rangle = \lambda| G_\lambda \rangle + \sqrt{1-|\lambda|^2}| G_\lambda^\perp \rangle
\end{equation}$$

>根据上式，再定义子空间$\mathcal{H}_\lambda = span\{| G_\lambda \rangle, | G_\lambda^\perp \rangle\}$。

### 问题的全貌

已知原矩阵的的标准形式编码$\hat{H} = (\langle G |_a \otimes \hat{\mathcal{I}}_s)\hat{U} (| G \rangle_a \otimes \hat{I}_s)$。如果存在标量函数$f$，在此函数作用下，$f[\hat{H}] = \displaystyle\sum_{\lambda} f(\lambda)| \lambda \rangle\langle \lambda | = (\langle G |_a \otimes \hat{\mathcal{I}}_s)\hat{U}' (| G \rangle_a \otimes \hat{I}_s)$。那么，如何寻找最优的量子线路，将$\hat{U}$转变为$\hat{U}'$？这就是量子信号处理器的任务。量子信号处理器(Quantum Signal Proccessor)将会解决：

* 输入：1. 埃尔米特矩阵$\hat{H}$，但要求其谱范数$\Vert \hat{H} \Vert \leq 1$；
        2. 函数$f: [-1, 1] \rightarrow \mathbb{D}$，其中$\mathbb{D}$是指复数单位圆，包括圆上；
* 输出：标准形式编码的$f[\hat{H}]$，或者说，一个酉矩阵$\hat{Q}$使得$f[\hat{H}] = \langle0|\hat{Q}|0\rangle$；
* 需求：1. 求得$\hat{H} = \langle G |\hat{U}| G \rangle$；
        2. 常数级的辅助量子位；
        3. 任意调用的单量子门与双量子门。

### $\hat{H}$的高次项

由于$f$是标量函数，那么我们自然会思考其中的一些高次项要如何计算。例如，$\hat{H}^2$。从量子信号处理器的定义中，我们已知$\hat{H} = \langle G |\hat{U}| G \rangle$，正常的想法就会认为：

$$\hat{H}^2 = \langle G |\hat{U}| G \rangle\langle G |\hat{U}| G \rangle = \langle G |\hat{U}^2| G \rangle$$

也就是利用$\hat{U}$的高次项来计算$\hat{H}$的高次项。但实际上，这并不正确。如果对$(1)$的两端再作一次$\hat{U}$，那么所得结果中会包含$\hat{U}| G_\psi^\perp \rangle$。虽然不知道它的结果是什么，但一般来说，这就导致一部分振幅泄漏到了$\mathcal{H}_{\lambda}$外。

>个人理解：换句话说，由于不能保证$(\langle G |_a \otimes \hat{\mathcal{I}}_s)\hat{U}| G_\psi^\perp \rangle_{as}=0$，所以类似于$(1)\sim(4)$的推导不能保证成立，也就不能保证$\hat{H}^2 = \langle G |\hat{U}^2| G \rangle$。

然而$\hat{H}$的高次项依然还是可以计算的。做法是，将$\hat{U}$替换为$\hat{W}$迭代子矩阵。虽然这么做依然使编码为$\hat{H}=\langle G |\hat{W}| G \rangle$，但不同点在于，$\hat{W}$有更严格的结构，使其能够迭代，使得迭代之后振幅依然保持在$\mathcal{H}_{\lambda}$内。为了做到这一点，需要保证$\hat{W}| G_\lambda \rangle$所在的子空间恰好与$| G_\lambda^\perp \rangle$相同，也就是$\mathcal{H}_\lambda = span\{ | G_\lambda \rangle, \hat{W}| G_\lambda \rangle \} = span\{ | G_\lambda \rangle, | G_\lambda^\perp \rangle \}$。这也变相通过施密特正交化定义了$| G_\lambda^\perp \rangle$，以及一系列作用在子空间上的泡利算符：

$$\begin{equation}
    | G_\lambda^\perp \rangle = \dfrac{(\hat{W}-\lambda)| G_\lambda \rangle}{\sqrt{1-|\lambda|^2}}
\end{equation}$$

$$\begin{equation}
    \begin{cases}
        \hat{X}_\lambda | G_\lambda \rangle = | G_\lambda^\perp \rangle \\
        \hat{Y}_\lambda | G_\lambda \rangle = i| G_\lambda^\perp \rangle \\
        \hat{Z}_\lambda | G_\lambda \rangle = | G_\lambda \rangle
    \end{cases}
\end{equation}$$

>$$\begin{cases}
    \hat{X}_\lambda = |G_\lambda \rangle\langle G_\lambda^\perp| + | G_\lambda^\perp \rangle\langle G_\lambda| \\
    \hat{Y}_\lambda = i|G_\lambda \rangle\langle G_\lambda^\perp| - i| G_\lambda^\perp \rangle\langle G_\lambda| \\
    \hat{Z}_\lambda = |G_\lambda \rangle\langle G_\lambda| - | G_\lambda^\perp \rangle\langle G_\lambda^\perp|
\end{cases}$$

从而得到：

$$\begin{equation}
    \begin{aligned}
        \hat{W} &= \begin{bmatrix}
        \lambda & -\sqrt{1-|\lambda|^2} \\
        \sqrt{1-|\lambda|^2} & \lambda \\
        \end{bmatrix}_\lambda \\
        &= \lambda |G_\lambda \rangle\langle G_\lambda| - \sqrt{1-|\lambda|^2} |G_\lambda \rangle\langle G_\lambda^\perp| + \sqrt{1-|\lambda|^2} |G_\lambda^\perp \rangle\langle G_\lambda | + \lambda |G_\lambda^\perp \rangle\langle G_\lambda^\perp |
    \end{aligned}
\end{equation}$$

那么，对于输入$| G \rangle \displaystyle\sum_{\lambda} a_\lambda | \lambda \rangle \in \mathcal{H}_G$有：

$$\begin{equation}
    \hat{W} = \bigoplus_{\lambda} \begin{bmatrix}
        \lambda & -\sqrt{1-|\lambda|^2} \\
        \sqrt{1-|\lambda|^2} & \lambda \\
        \end{bmatrix}_\lambda = \bigoplus_{\lambda} e^{-i\hat{Y}_\lambda \theta_\lambda}
\end{equation}$$

其中，$\theta_\lambda = \arccos(\lambda)$。

## 二、量子位化(Qubitization)

经过以上过程，使用$\hat{W}$即可有效地实现$\hat{H}$的高次项。但是，构造$\hat{W}$并不简单，因为对于每个本征态，都需要计算$\sqrt{1-|\lambda|^2}$(可以使用相位估计解决)。为了系统性地从标准形式编码转换到$\hat{W}$，使用所谓量子位化(Qubitization)来实现。

假设$| G \rangle$已知(虽然这个条件不是必要的)，同时有一个酉矩阵$\hat{S}'$只作用在辅助空间上，使得$\hat{W}=\hat{S}'\hat{U}$。为了得到$\hat{S}'$，再构造$\hat{S}$使得：

$$\begin{equation}
   \hat{W} = ((2|G \rangle\langle G | - \hat{I})_a \otimes \hat{I}_s)\hat{S}\hat{U}
\end{equation}$$

$$\begin{equation}
   \begin{cases}
       | G_\lambda \rangle = | G \rangle | \lambda \rangle \\
       | G_\lambda^\perp \rangle = \dfrac{\lambda | G_\lambda \rangle - \hat{S}\hat{U}| G_\lambda \rangle}{\sqrt{1-\lambda^2}}
   \end{cases}
\end{equation}$$

### 条件与难点

允许量子位化的条件是$\langle G |_a \hat{S}\hat{U} | G \rangle_a = \hat{H}$的同时，$\langle G |_a \hat{S}\hat{U}\hat{S}\hat{U} | G \rangle_a = \hat{\mathcal{I}}$。第二个条件暗示了$\hat{S}\hat{U}$是一个类似于Grover算法的反射算符，也要求$\hat{S}\hat{U}\hat{S} = \hat{U}^\dagger$。以上这些是相当严苛的条件，对于一般的$\hat{U}$，$\hat{S}$很难存在。为了解决这个问题，又重新构造一个量子线路$\hat{U}'$，这个线路中包含$\hat{U}$且总是有简单的解$\hat{S}$。构造受控$\hat{U}$门：

$$\begin{cases}
    \hat{V}_1 = | 0 \rangle\langle 0 | \otimes \hat{I} + | 1 \rangle\langle 1 | \otimes \hat{U}^\dagger \\
    \hat{V}_2 = | 0 \rangle\langle 0 | \otimes \hat{U} + | 1 \rangle\langle 1 | \otimes \hat{I} \\
\end{cases}$$

然后相乘得到$\hat{U}'$：

$$\hat{U}' = \hat{V}_2 \hat{V}_1 = | 0 \rangle\langle 0 | \otimes \hat{U} + | 1 \rangle\langle 1 | \otimes \hat{U}^\dagger $$

也就是说，以上构造了一个根据额外引入的一个量子位，选择作用$\hat{U}$抑或$\hat{U}^\dagger$。现在再令这个引入的量子位为$\dfrac{1}{\sqrt{2}}(| 0 \rangle + | 1 \rangle)$、并并入到辅助空间中，此时辅助态为$| G' \rangle = \dfrac{1}{\sqrt{2}}(| 0 \rangle + | 1 \rangle)| G \rangle$。又令$\hat{S} = (| 0 \rangle\langle 1 | + | 1 \rangle\langle 0 |) \otimes \hat{\mathcal{I}}_{as}$，此时重新计算量子位化的条件的左边，将$\hat{U}$替换为$\hat{U}'$：

$$\begin{equation}
    \langle G' | \hat{S}\hat{U}' | G' \rangle = \dfrac{1}{2}(\hat{H} + \hat{H}^\dagger) = \hat{H}, \qquad
    \langle G' | \hat{S}\hat{U}'\hat{S}\hat{U}' | G' \rangle = \langle G' | \hat{U}'^\dagger\hat{U}' | G' \rangle = \hat{\mathcal{I}}
\end{equation}$$

使用上述方法，注意到对$| G \rangle$没有任何要求，且$(12)$是普适的，因此对于任何$\hat{G},\hat{U}$都能成功实现量子位化。

### 添加参数

在$(9)$的基础上，添加参数$\phi$，这时我们称$\hat{W}_\phi$为相位迭代子：

$$\begin{equation}
    \hat{W}_\phi = \bigoplus_{\lambda} \begin{bmatrix}
        \lambda & -ie^{-i\phi}\sqrt{1-|\lambda|^2} \\
        -ie^{-i\phi}\sqrt{1-|\lambda|^2} & \lambda \\
        \end{bmatrix}_\lambda = \bigoplus_{\lambda} e^{-i\hat{\phi}^\lambda \theta_\lambda}
\end{equation}$$

这里$\hat{\phi}^\lambda = cos\phi\ \hat{X}_\lambda + sin\phi\ \hat{Y}_\lambda$。如果仔细计算，会发现这一带相位的迭代子可以从不带相位的迭代子得到：

$$\begin{equation}
    \hat{W}_\phi = \hat{Z}_{\phi+\pi/2} \hat{W} \hat{Z}_{-\phi-\pi/2}
\end{equation}$$

而$\hat{Z}_\phi = (1+e^{-i\phi})| G \rangle\langle G | - \hat{\mathcal{I}}$，这个形式已是非常熟悉的了，就是关于$| G \rangle$做对称并附加相位，详细可见之前所讲的振幅放大。当然，也可以将$\hat{Z}_\phi$展开为矩阵形式：

$$\begin{equation}
    \hat{Z}_\phi = \bigoplus_{\lambda} \begin{bmatrix}
        e^{-i\phi} & 0 \\
        0 & -1
    \end{bmatrix}_\lambda
\end{equation}$$

### 奇偶性相同时


