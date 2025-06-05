第13章 协同进化中问题结构的影响  
13.1 协同进化问题结构及其对搜索的影响简介

协同进化计算已在众多不同设置下的各类问题中得到应用，例如战略决策 [1, 11, 12, 25]、机器学习问题（如分类） [26, 65]、优化 [26, 34, 44, 56, 66] 以及人工生命模拟 [23, 60]。这些方法已经被用于不同的问题领域（见第4章）。鉴于协同进化算法（CEA）所涉及的问题范围极其广泛，因此对在类似问题设置下应用CEA的相关文献进行仔细综述几乎成为必要。不过，实践者和算法设计者依然会合理地提出这样的问题：是否存在一些通用的原则，能够指导CEA的指定与定制以解决手头的问题？或者，换个角度看，有哪些相关且关键的问题结构会影响协同进化搜索的性能？本章旨在介绍并展示我们对一些相关但一般性问题结构所进行的正式研究，以及这些结构如何具体与协同进化搜索过程相关联，从而可以提取出一些用于指导CEA性能的通用原则。这需要首先明确所研究的协同进化问题的范围，以便更好地解释我们后续展示的结果。

如我们在第2.3节中已介绍，CEA主要应用于两种通用但不同的设置：竞争性和合作性。简而言之，竞争性CEA适用于测试型问题，此类问题需要通过测试用例来评估生成候选解的质量；而合作性CEA则适用于组合型问题，此时需要通过组装各个组件来构建候选解 [53]。需要注意的是，在组合型问题中，候选解的评估对于更具挑战性和复杂性的问题可能也涉及测试用例的使用，或者可用目标函数。本质上，首先要区分清楚哪些问题必须通过测试用例进行评估，哪些问题则可以通过查询某个函数来为候选解提供具体的目标质量度量，这一点至关重要。
特别地，协作型协同进化算法（cooperative CEAs）的主要原理是通过“分而治之”（Divide-and-Conquer）策略，以并行搜索的方式解决问题。如我们在第3.3.2节中所讨论，这要求优化问题可以被分解为若干子问题，每个子问题由相互隔离的种群分别独立求解。协作搜索方法的优势在于，这些种群之间的交互能够同步改进各自的子组件，从而协同推动全局搜索的进行。在连续优化问题[45]中，变量之间的相互作用是一个至关重要的问题结构。例如，不可分离问题会表现出不同程度的变量交互。文献[54]的研究表明，这种变量交互会对协作型协同进化算法的性能产生影响。特别地，文中采用了详细的案例研究，分析了在二维连续优化问题中，不可分离性如何影响协同进化算法的设计决策，比如精英策略（elitism）的使用和种群规模的设定。

在这种问题设置下，协作型协同进化算法自然需要采用双种群协同搜索（每个种群分别负责一个特定变量），而两个变量之间存在相互作用。这种相互作用可以通过最佳响应曲线（best response curves）来刻画，即优化函数关于某一变量的偏导数。假设以最小化为目标，对 $x_2$ 的最佳响应曲线是类型为 ${\mathbb{R}}_{x_2} \rightarrow {\mathbb{R}}_{x_1}$ 的函数，定义为
$$
g_{\mathrm{best}}(x_2) = \mathop{\mathrm{argmin}}_{x_1 \in {\mathbb{R}}_{x_1}} f(x_1, x_2) 
$$
对 $x_1$ 的最佳响应曲线同理可以定义。需要注意的是，该研究所选的优化函数使得最佳响应曲线始终是一一对应的关系。更重要的是，这两条最佳响应曲线的几何结构及其相互关系具有重要的解释意义。基于这一点，可以构建一个理想化的协同进化搜索过程。具体来看，首先将上述这两条分别类型为 ${\mathbb{R}}_{x_2} \rightarrow {\mathbb{R}}_{x_1}$ 和 ${\mathbb{R}}_{x_1} \rightarrow {\mathbb{R}}_{x_2}$ 的最佳响应曲线嵌入到 ${\mathbb{R}}^2$ 平面中。在理想情况下，一次进化搜索循环的终极结果就是只采样来自相应最佳响应曲线上的点。在一个协同进化周期内，协作型协同进化算法的搜索轨迹（即采用最优合作体策略best-of-generation）可以想象为：首先从某条最佳响应曲线采样得到一个点，再沿着曲线之间的方向（根据两条曲线的几何关系，选择水平或垂直）移动至另一条曲线上，再重复该过程，从而回到第一条曲线上。整个轨迹就是 ${\mathbb{R}}^2$ 空间中这些“跳跃”点的有序连接序列。

对于更简单的偏轴二次函数（Off Axis Quadratic function）来说，搜索轨迹类似于一个以全局最优为吸引子的负反馈系统[24]。然而，在如Rosenbrock等复杂的不可分离问题中，变量之间的强交互和非线性可能对这种协作型协同进化搜索的性能产生负面影响。此时，最佳响应曲线是非线性的，算法生成的搜索轨迹甚至可能导致协作搜索过程中，完整候选解 $\boldsymbol{x} = (x_1, x_2)$ 实际上远离全局最优解。
在本章中，我们主要关注于竞争性协同进化问题的结构分析，以及这些结构对协同进化搜索的影响，内容基于我们在文献[13, 14]中的研究。我们需要的数学工具主要针对大规模竞争性协同进化问题，这类问题中，解-测试用例之间的交互被建模为双人策略博弈。我们已经在第5.3.1节中介绍了协同进化有向图（coevolutionary digraph），该有向图能够完整地捕捉此类问题的相关结构。本章剩余部分组织如下：13.2节首先将抽象协同进化表述为在协同进化有向图上的随机游走（有限状态马尔可夫链），并据此提出关键的理论结果，为阐明协同进化问题对于搜索过程具有挑战性的原因提供了重要的定性见解。我们对协同进化马尔可夫链的理论同样被用于开发定量分析工具，以刻画给定问题（循环）复杂度下协同进化搜索的速度。第13.3节进一步提出了一项理论研究，揭示了PageRank与抽象协同进化之间的深刻联系。PageRank曾用于刻画网页网络中网页节点的重要性，我们建立了PageRank权威（PageRank authorities）的概念可用于指示协同进化有向图中顶点的重要性（性能）。至关重要的是，PageRank权威在带有重启机制的有向图上，可以自然地、二重地解释为协同进化搜索中对顶点的访问概率。第13.4节以简要评论结束全章，强调为更好理解协同进化系统而采用的多样化理论工具及其间的联系。

13.2 随机游走在有向图上作为协同进化搜索过程的模型

在本节中，我们[14]给出了一种理论化的种群一体（population-one）竞争性协同进化系统的表述，这一表述将使我们能够深入研究协同进化问题中结构——即在相关有向图中以顶点对关系完全捕捉的结构——如何影响搜索过程。第13.2.1节将引入抽象协同进化，将其视为在协同进化有向图上的随机游走，用自然的采样过程刻画协同进化搜索过程。这些马尔可夫链使我们得以研究协同进化有向图中的特定结构如何影响协同进化搜索过程的动态特性（第13.2.2节）。这些分析为竞争性问题中令协同进化搜索具有挑战性的具体结构提供了关键且易于解释的见解。随后在第13.2.3节，我们将展示如何应用我们所建立的协同进化马尔可夫链理论，来开发定量分析工具，从而用于判定协同进化在特定问题（循环结构）下找到主导解的速度。

13.2.1 
协同进化搜索过程作为有向图上的随机游走
我们在第5.3.1节中将有向图 $D = (V, A)$ 上的游走定义为顶点和弧交替出现的序列 $v_1 a_1 v_2 a_2 v_3 \ldots a_{k-1} v_k$，其中 $v_1, \ldots, v_k \in V$，$a_1, \ldots, a_{k-1} \in A$。需要注意的是，游走必须满足：对于所有 $i$，存在弧 $v_i a_i v_{i + 1}$。对于至少一个 $v_i$，其出邻居可能有多个，这些出邻居由集合 $N_{D}^{+}(v_i) = \{u \in V \setminus \{v_i\} : v_i \rightarrow u \in A\}$ 给出。这意味着，可能存在多个策略 $u$ 可以击败 $v_k$。因此，在有向图中，两个顶点 $v_1, v_k$ 之间可能存在多条不同的游走路线。这种不同顶点对之间游走集合的扩展同样适用于图中任意两顶点之间。随机游走的概念可以非正式地视为在有向图上的一种采样过程。更具体而言，该过程会根据某种潜在分布，从有向图中两顶点间定义的游走集合中进行采样。随机游走可以以过程化的方式实现。算法13.1描述了有标签有向图 $D$ 上的标准随机游走（此处“有标签”指的是对每个顶点都已知其出邻居的先验信息）。该算法从一个随机顶点出发，在每一步中，以均匀概率随机移动到当前顶点的一个出邻居。这一过程化内容在算法13.1中有具体描述。在实际应用中，人们通常会在无标签有向图如 $(1+1)$ CEA（算法13.2）上实现随机游走。在 $(1+1)$ CEA 中，每一步的随机移动通过变异（从顶点集 $V$ 中生成一个新的随机后代）和选择（当后代顶点支配父代顶点时移动至后代）两个阶段完成。因此，$(1+1)$ CEA 可能在某一步仍停留在当前顶点。只有在到达 $D$ 的支配顶点时，标准随机游走才会出现这种停留的情形。还要注意，实际操作中的 $(1+1)$ CEA 可以采用任意终止准则（如算法13.2中的有限步数）。

算法13.1 有向图 $D$ 上的随机游走 [14]
输入：共演有向图 $D = (V, A)$，$u$ 当前代表博弈策略的顶点
1: 过程 RWALKD($D$, $u$)
2: $v := \text{uniSamp}(D.V)$ // $v \in V$，以均匀分布随机选取
3: 重复
4: $N_{D}^{+}(v) := \text{outNeigh}(v)$ // $N_{D}^{+}(v) = \{u \in V \setminus \{v\} : v_i \rightarrow u \in A\}$
5: 如果 $N_{D}^{+}(v) \neq \emptyset$ // $N_{D}^{+}(v) \subset V$ 是击败 $v$ 的策略集合
6: $u := \text{uniSamp}(N_{D}^{+}(v))$ // $u \in N_{D}^{+}(v)$，以概率 $1/|N_{D}^{+}(v)|$ 随机选取
7: 否则 $u$ 保持不变
8: $v := u$
9: 直到 永远
10: 过程结束

算法13.2 $(1+1)$ CEA 在有向图 $D$ 上 [14]
DCoevolutionary Digraph $D = (V, A)$

当前顶点表示一个博弈策略
$b_{gen}$ 终止准则，固定的迭代代数

1: 过程 CEA($D, u, b_{gen}$)
2: $t := 0$ \hspace{1em} // 初始化时间戳
3: $v_t := \text{uniSamp}(D.V)$ \hspace{1em} // 从顶点集$V$中均匀随机选择$v \in V$
4: $u := v_t$ \hspace{1em} // 用单个顶点初始化父代$u$
5: 当 $t \leq b_{gen}$ 时，执行
6: \hspace{1em} $v_t := \text{mutate}(u)$  \hspace{1em} // 通过从$V \setminus \{u\}$均匀采样生成后代
7: \hspace{1em} 如果 $v_t \leftarrow u$ 则
8: \hspace{2em} $u := v_t$
9: \hspace{1em} 结束如果 \hspace{1em} // 选择操作在$\{u, v_t\}$中选出获胜策略
10: \hspace{1em} $t := t + 1$
11: 结束当
12: 过程结束

无论如何实际实现有向图上的随机游走，我们在本节关注的核心是生成这些游走的底层随机过程。更具体地说，我们希望将这些种群规模为1的共演搜索过程建模为离散时间、有限状态的马尔可夫链 [47, 50]。我们将抽象共演[14]作为有向共演图$\varPhi$上的随机游走，即$\{\varPhi_t : t \in \mathbb{N}_{\geq 0}\}$组成的随机变量序列，其中每个$\varPhi_t$在时间$t$取值于一个可数集合$X$，状态空间$X$对应于共演有向图$D = (V, A)$的有限顶点集$V(D)$。

然而，要将一般抽象共演过程构建为离散时间随机过程在技术上非常复杂，因为需要在概率空间上进行操作。一般来说，首先要对样本路径行为和有向图上随时间演化的步进动力学给予结构上的定义，以此刻画控制马尔可夫链$\varPhi$演化的概率律。即便在常见的情形下，序列$\{\varPhi_0, \varPhi_1, \varPhi_2,\ldots\}$在$X$中取离散值，如果有向图的结构不严格限制其样本路径的行为（例如无环图和正规锦标赛），想要在大多数有向图上对随机游走的分布给出完整的、具备操作性的描述仍然非常困难。

因此，我们注意到，大量受限的共演搜索过程满足马尔可夫性，即对未来状态的转移仅依赖于当前状态。例如，在$(1+1)$ CEA的情形下，从特定顶点$x_t = u \in V(D)$到其出邻居$x_{t+1} = v \in N^+_D(u)$的转移概率并不依赖于时间（世代）步$t$，即无论是在初始时$x_0 = u$还是在$t$代后的$x_t = u$，该转移概率均相同。我们利用这一马尔可夫性，形式化并定义共演过程中有向图上的马尔可夫链如下：

设$D(V_{S(n)}, A) \in \mathcal{D}(V_{S(n)})$为一个拥有有限顶点集$V_{S(n)}$（大小为$n$）的共演有向图。共演马尔可夫链（CMC）[14]是在$D(V_{S(n)}, A)$上的随机游走，其初始分布为$\boldsymbol{\mu}$，马尔可夫转移矩阵为$\boldsymbol{P}$，状态空间为$X = V_{S(n)}$，满足：

1. $\boldsymbol{\mu} = (\mu_x : x \in X)$，其中 $0 \leq \mu_x \leq 1$ 且 $\sum_{x \in X} \mu_x = 1$。
2. $\boldsymbol{P} = (\mathbb{P}(x,z) : x, z \in X)$，其中每一行都是概率分布，满足 $0 \leq \mathbb{P}(x, z) \leq 1$ 且 $\sum_{y \in X} \mathbb{P}(x, y) = 1$。
由于$V_{S(n)}$包含有限个元素，并且按定义是可数集，根据可数状态空间马尔可夫链的存在定理[18]，可以保证构建出定义在状态空间$X = V_{S(n)}$上的如下齐次马尔可夫链$\varPhi = \{\varPhi_t : t \in {\mathbb{N}}_0\}$，其初始分布为$\boldsymbol \mu$，转移矩阵为$\boldsymbol P$。用$\boldsymbol \mu$和$\boldsymbol P$可以计算描述齐次马尔可夫链（CMC）的分布。我们以初始分布为$x_0 \in X$的定点分布（点质量分布）启动马尔可夫链。随后，对于每个$x_0 \in X$，递归地定义该马尔可夫链的$\mathbb{P}_{x_0}$如下：
\[
\mathbb{P}_{x_0}(\varPhi_0 = x_0) = 1, \\
\mathbb{P}_{x_0}(\varPhi_0 = x_0, \varPhi_1 = x_1) = \mathbb{P}(x_0, x_1), \\
\mathbb{P}_{x_0}(\varPhi_0 = x_0, \varPhi_1 = x_1, \varPhi_2 = x_2) = \mathbb{P}(x_0, x_1)\mathbb{P}(x_1, x_2)，
\]
以此类推，利用$\boldsymbol P$给出该链的一步转移概率$\mathbb{P}(x_t, x_{t+1})$[47]。此外，还可以获得其他有用信息，例如，链起始于$x$，经过$n$步后到达$y$的概率为
\[
\mathbb{P}_{\mu}(\varPhi_n = y \ | \ \varPhi_0 = x) = \mathbb{P}^n(x,y).
\]
首先设$\boldsymbol P^0 = \mathbf{I}$（单位矩阵），即初始条件$\mathbb{P}^0(x,x) = 1$对所有$x \in X$均成立。则$n$步转移矩阵的递归定义为
\[
\boldsymbol P^n = (\mathbb{P}^n(x, y) : x, y \in X)
\]
对于$n \geq 0$，可通过迭代实现$\boldsymbol P^{n+1} = \boldsymbol P^n \boldsymbol P$，其递推公式为[14]：
\[
\mathbb{P}^{n+1}(x, z) = \sum_{y \in X} \mathbb{P}^n(x, y)\mathbb{P}(y, z).
\]
接下来，我们给出一些具体的CMC构造实例。首先考虑两个定义在传递锦标赛上的CMC例子。初始分布可以任意，因此只需构造适当的马尔可夫转移矩阵。第一个例子为带标签的传递共进化锦标赛上的随机游走。设$T(V_{S(m)}, A)$为一个传递共进化锦标赛，其顶点$(v_1, v_2, \ldots, v_m)$按其得分序列$(s_1, s_2, \ldots, s_m)$进行索引（注意：共进化有向图中一个顶点$v$的得分为其入弧数量，等于$v$所支配的顶点个数）。为简化记号，令$\boldsymbol P_{m \times m} = (p_{ij} : i, j \in \{1, 2, \ldots, m\})$为该CMC在$T(V_{S(m)}, A)$上的马尔可夫转移矩阵。$\boldsymbol P_{m \times m}$递归定义如下。基例（$n = 0$）为$T(V_{S(m)}, A)$阶数$m = n + 1 = 1$。此时$T(V_{S(m)}, A)$只有一个顶点，同时也是唯一的获胜顶点，即$v_m = v_1$，其自环转移概率为$p_{mm} = 1$，因此$\boldsymbol P_{1 \times 1} = [p_{11}] = [1]$。归纳步骤（$n+1$）增加一个顶点，使得锦标赛阶数$m = n+1$。新加入的顶点被索引为$v_1$。存在$n = m - 1$个顶点从$v_2$到$v_m$支配$v_1$（见图13.1a）。因此，从$v_1$出发到其他顶点的转移概率相等。

图13.1：(a) 带标签的传递锦标赛$T(V_{S(n+1)}, A)$上的随机游走的CMC； (b) 传递锦标赛$T(V_{S(n+1)}, A)$上的(1+1)共进化算法的CMC [14]。
对于每一个 $\{v_2, v_3, \ldots, v_m\}$，有如下概率：
$$
p_{11} = 0, \\
p_{1j} = \frac{1}{m-1}, \quad j = 2,3,\ldots,m, \\
p_{i1} = 0, \quad i = 2,3,\ldots,m.
$$
我们接着可以通过 $m = 1, 2, \ldots, n+1$ 的归纳方式，构造$n+1$阶标记传递锦标赛$T(V_{S(n+1)}, A)$上随机游走的马尔可夫转移矩阵，其概率为：
$$
p_{ii} = 0, \quad i = m \\
p_{ij} = \frac{1}{n+1-m}, \quad i = m, \ j = i+1, i+2, \ldots, n+1, \\
p_{ij} = 0, \quad i = j+1, j+2, \ldots, n+1, \ j = m,
$$
并且有 $p_{(n+1)(n+1)} = 1$。

这样便可以得到如下的马尔可夫转移矩阵$P$ $(n+1) \times (n+1)$：
$$
P = 
\begin{pmatrix}
0 & \frac{1}{n} & \frac{1}{n} & \cdots & \frac{1}{n} \\
\frac{1}{n-1} & 0 & \frac{1}{n-1} & \cdots & \frac{1}{n-1} & 0 \\
\vdots & \vdots & \ddots & & \vdots \\
\frac{1}{3} & \frac{1}{3} & \frac{1}{3} & 0 & 0 & \cdots & 0 \\
0 & \frac{1}{2} & \frac{1}{2} & 0 & 0 & \cdots & 0 \\
0 & 0 & 1 & 0 & 0 & \cdots & 0 \\
0 & 0 & 1
\end{pmatrix},
\tag{13.2}
$$
其中最后一行表示锦标赛中唯一支配顶点的转移概率，用水平线突出显示 [14]。

下一个例子是关于无标签传递协同进化锦标赛上的$(1+1)$ CEA。前述马尔可夫转移矩阵具有良好的归纳结构（右下子矩阵诱导了在相应标记子锦标赛上的随机游走）。而对于$(1+1)$ CEA，我们需要直接计算其在不同阶数的传递协同进化锦标赛上的马尔可夫转移矩阵。

幸运的是，转移概率的计算较为直接。考虑一个$n+1$阶的传递协同进化锦标赛$T(V_{S(n+1)},A)$。其转移概率如下给出 [14]：
上三角部分：$p_{ij} = 1/n, \quad i < j + 1, \quad i = 2, 3, \ldots, n+1$

对角线部分：$p_{ii} = m/n, \quad m = i - 1, \quad i = 2, 3, \ldots, n+1$

下三角部分：$p_{ij} = 0, \quad i > j - 1, \quad i = 2, 3, \ldots, n+1$

$P$ 为 $(n+1) \times (n+1)$的矩阵，其形式如下：

$$
P_{(n+1)\times(n+1)} = 
\begin{pmatrix}
0 & \frac{1}{n} & \frac{1}{n} & \cdots & \frac{1}{n} \\
\frac{1}{n} & 0 & \frac{1}{n} & \cdots & \frac{1}{n} \\
\frac{1}{n} & \frac{1}{n} & 0 & \cdots & \frac{2}{n} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
\frac{1}{n} & \frac{1}{n} & \cdots & 0 & \frac{n-2}{n} \\
0 & \frac{n-1}{n} & \frac{1}{n} & 0 & \cdots & 1 \\
0 & 0 & 1
\end{pmatrix}
\tag{13.3}
$$

其对应的马尔可夫链在图 13.1b 中进行了说明。我们现在考虑关于规则巡回赛（标签与非标签）的两个 CMC 示例。设 $T(V_{S(n)}, A)$ 为一个奇数阶 $n \geq (3 + 2m)$ 的规则共演进巡回赛，其中 $m = 0, 1, 2, \ldots$。它的邻接矩阵 $\mathbf{M}_{\mathrm{adj}}$ 可以按如下方式构造 [7]：

$$
\mathbf{M}_{\mathrm{adj}} = \mathbf{L}_n + \mathbf{L}_n^2 + \mathbf{L}_n^3 + \ldots + \mathbf{L}_n^{(n-1)/2}
$$

其中，$\mathbf{L}_n$ 是由排列向量 $(2, 3, 4, \ldots, n, 1)$ 得到的置换矩阵。排列向量 $(r_1, r_2, r_3, \ldots, r_m)$ 指定在 $m \times m$ 单位矩阵上的行置换操作。对第一行来说，第 1 列的 1 会与第 $r_1$ 列的 0 交换，依此类推，直到第 $m$ 行，其中第 $m$ 列的 1 会与第 $r_m$ 列的 0 交换。例如，排列向量 $(2, 3, 1)$ 所对应的置换矩阵为：

$$
\begin{bmatrix}
0 & 1 & 0 \\
0 & 0 & 1 \\
1 & 0 & 0 \\
\end{bmatrix}
$$

在邻接矩阵中，$\mathrm{adj}(i, j) = 1$ 表示对应的有向图中存在从 $i$ 到 $j$ 的弧（$i \rightarrow j$）。马尔可夫转移矩阵需要我们从该有向图的支配矩阵 $\mathbf{M}_{\mathrm{dom}}$ 计算，其中 $\mathrm{dom}(i, j) = 1$ 表示 $i \leftarrow j$。注意，$\mathbf{M}_{\mathrm{dom}}$ 是 $\mathbf{M}_{\mathrm{adj}}$ 的转置，即：

$$
\mathbf{M}_{\mathrm{dom}} = \mathbf{M}_{\mathrm{adj}}^T
$$

因此，对于标记规则巡回赛 $T(V_{S(n)}, A)$ 的随机游走，其马尔可夫转移矩阵可表示为 [14]：

$$
\mathbf{P} = \frac{2}{(n-1)}\mathbf{M}_{\mathrm{dom}}
$$

$\tag{13.4}$
而对于无标签的正则图$T(V_{S(n)}, A)$上的$(1 + 1)$CEA，其转移矩阵为[14]：$\mathbf{P} = \frac{1}{2}\mathbf{I} + \frac{1}{(n-1)}\mathbf{M}_{\mathrm{dom}}$。

这些例子突出了以下情形：我们能够形式化地进行结构构建，并在下一个小节进行分析，而无需像实际问题中常见的那样依赖采样和大规模计算。

### 13.2.2 协同进化马尔可夫链的刻画

马尔可夫链中通信类（communicating classes）的概念对于其性质刻画非常重要。设$x, y$为马尔可夫链的两个不同状态，$x, y \in X = V_{S(n)}$。我们有$x \leadsto y$当且仅当存在某个$n \geq 0$，使得从$x$出发，马尔可夫链经过$n$步可以以正概率$p_{xy}^{(n)} > 0$到达$y$。关系$\leadsto$按照定义是自反的（即$x \leadsto x$），并且根据$n$步转移概率的性质（由公式(13.1)给出）是传递的，因为$p_{xz}^{(m+n)} \geq p_{xy}^{(m)}p_{yz}^{(n)} > 0$，其中$m \geq 0$且$n \geq 0$。于是有$x \leftrightsquigarrow y$（即$x$与$y$“通信”）当且仅当$x \leadsto y$且$y \leadsto x$。关系$\leftrightsquigarrow$是对称且传递的，因此在$X$上构成一个等价关系，从而能够将$X$划分为不相交的子集，称为通信等价类（communicating classes）[18]。设$C(x)$为包含$x \in X$的通信等价类。若状态空间$X$仅包含一个通信类，即$C(x) = X$，则称该马尔可夫链是不可约的（irreducible）。

在文献[14]中，我们形式化地证明了有向图理论中的强连通性的概念、博弈论中的可解性概念以及马尔可夫链中的通信性之间的联系。这一点在引理13.1和13.2中总结，它们共同为三者的统一奠定基础，使我们能够对竞争性协同进化环境做出定性的分析，涵盖了从问题结构（以及通过两两关系确定的解的结构）到协同进化搜索过程的各个方面。在对竞争性协同进化的定性刻画中，这种二分性从问题结构转移到了搜索过程。具体而言，如果CMC是不可约的，则它作用于一个强连通（不可约）的协同进化有向图，该图对应的是不可解的协同进化问题；反之，吸收型CMC具有唯一的吸收类，对应协同进化有向图中包含可解协同进化问题解的主导子集。

**引理13.1 ([14])** 如果CMC作用于不可约的协同进化有向图，则该CMC是不可约的。

**引理13.2 ([14])** 对于可约的CMC，只有唯一的吸收类。吸收类由协同进化有向图中对应主导解的状态组成。
除了前文对抽象协同进化所作的一般性定性描述之外，当将其应用于具体问题时，可能还需要更为细致的定量表征。对此类协同马尔可夫链（CMC）行为的一项相关分析，是研究其在$X$空间中演化时，于某些随机时刻的分布。下文将对这类量进行精确定义与描述。特别地，其中一类量指的是CMC $\boldsymbol{\varPhi}$访问$X$中某些状态的随机时刻。设$\mathcal{A} \subset \mathcal{B}(X)$，其中$\mathcal{B}(X)$为$X$的所有子集的集合。集合$\mathcal{A}$的首次到达时间（first hitting time）是随机变量$\tau_{\mathcal{A}}: X \rightarrow \mathbb{N}_0 \cup \{\infty\}$，定义为
$$
\tau_{\mathcal{A}} = \mathrm{inf}\{t \geq 0 : \varPhi_t \in \mathcal{A}\}.
$$
对于状态$x \in X$的首次到达时间简写为$\tau_{x} = \mathrm{inf}\{t \geq 0 : \varPhi_t = x\}$。按照惯例[50]，对空集$\emptyset$的首次到达时间定义为$\infty$。此外，可以关注起始状态不属于$\mathcal{A}$的情形，从而定义对$\mathcal{A}$的首次访问时间为
$$
\tau_{\mathcal{A}}^{+} = \mathrm{inf}\{t \geq 1 : \varPhi_t \in \mathcal{A}\}
$$
[47]。类似于进化算法分析中的情形[31]，研究可吸收CMC $\boldsymbol{\varPhi}$ 的平均首次到达时间具有重要意义。可吸收CMC具有两点特性，而这在进化算法模型中并不总是成立。首先，其状态空间总能被划分为两类——一类为封闭的吸收类（即遍历类），另一类仅由跃迁态（transient states）构成。其次，协同进化有向图内的完全连通性意味着$\tau_{\mathcal{A}}^{+} = 1$，尽管$\boldsymbol{\varPhi}$一步到达$\mathcal{A}$通常是极小概率事件。尽管如此，形式上我们可以定义，从$x \in X$出发达到非空且可达的$\mathcal{A}$时的期望首次到达时间为
$$
\mathbb{E}_{x}(\tau_{\mathcal{A}}) = \sum_{t < \infty} t \mathbb{P}_x(\tau_{\mathcal{A}} = t),
$$
其中$\mathbb{P}_x(\tau_{\mathcal{A}} = t)$为起始于$x \in X$时$\boldsymbol{\varPhi}$在时刻$t$第一次到达$\mathcal{A}$的概率[50]。在假设$\mathcal{A}$非空且可达的前提下，定义上$\mathbb{P}_x(\tau_{\mathcal{A}} = \infty) = 0$。我们采用期望首次到达时间$\mathbb{E}_x(\tau_{\mathcal{A}})$作为协同进化搜索$\mathcal{A}$中主导解速度的量化指标[14]。在计算$\mathbb{E}_x(\tau_{\mathcal{A}})$之前，首先需对可约CMC的状态重新编号，并对其做适当简化，得到其马尔可夫转移矩阵的规范形式$\boldsymbol{P}^*$ [46]
$$
\boldsymbol{P}^* = \begin{bmatrix}
\mathbf{I} & \mathbf{0} \\
\mathbf{R} & \mathbf{Q}
\end{bmatrix}.
$$
其中，$\mathbf{R}$表示从跃迁态到吸收类的转移概率，$\mathbf{Q}$表示跃迁态之间的转移概率，$\mathbf{0}$为全零行向量。随后，文献[14]给出如下理论结果，使得可直接计算任意跃迁态$i \in Q$出发时的期望首次到达时间$\mathbb{E}_i(\tau_{\mathcal{A}})$。
引理 13.3 ([14]) 考虑一个在其标准形下具有马尔可夫转移矩阵 $P^*$ 的可约连续马尔可夫链（CMC）。令 $i \in \mathcal{Q}$ 为瞬态状态。定义列向量 $\mathbf{h} = [h_i]_{i \in \mathcal{Q}}$，其中 $h_i = \mathbb{E}_i(\tau_{\mathcal{A}})$ 表示从瞬态状态 $i$ 出发到达吸收类 $\mathcal{A}$ 的期望命中时间。则有：

$${\mathbf{h}} = ({\mathbf{I}} - {\mathbf{Q}})^{-1}\mathbf{1}, \quad h = (I − Q)^{-1}1, \tag{13.8}$$

其中 $\mathbf{1}$ 为元素全为1的列向量。该结论的技术性证明利用了每一个可约连续马尔可夫链始终是具有一个吸收类的吸收马尔可夫链（此时可应用 [36] 中的定理3.2）。对于吸收马尔可夫链，矩阵 $({\mathbf{I}} - {\mathbf{Q}})^{-1}$ 存在，并被称为马尔可夫转移矩阵 ${\boldsymbol P}$ 的基本矩阵。通过首先计算其基本矩阵，可以直接计算任意可约连续马尔可夫链的期望命中时间（见 [29] 中的定理11.4）：

$$({\mathbf{I}} - {\mathbf{Q}})^{-1} = {\mathbf{I}} + {\mathbf{Q}} + {\mathbf{Q}}^{2} + {\mathbf{Q}}^{3} + \cdots. $$
$$(I - Q)^{-1} = I + Q + Q^2 + Q^3 + \cdots $$

由于 $\mathbf{Q}^n \rightarrow \mathbf{0}$ 当 $n \rightarrow \infty$ 时，该级数收敛（见 [29] 中的定理11.3）。现在，我们可以为可传递共进化竞赛上的特定CMC分析性地给出期望命中时间的计算表达式。注意，公式 (13.2) 和 (13.3) 所给出的马尔可夫转移矩阵均为上三角矩阵。我们对相应的 ${\boldsymbol P}$ 重新标记，使其变为下三角矩阵。对于这两类特定的CMC，${\mathbf{Q}}$ 也呈现下三角矩阵形式。关键在于，每一个 $n \times n$ 的矩阵 $({\mathbf{I}} - {\mathbf{Q}})$ 都有 $\mathrm{rank}({\mathbf{I}} - {\mathbf{Q}}) = n$，因此其逆矩阵存在 [46]。对于阶数为 $n+1$ 的有标签可传递共进化竞赛上的随机游走，容易得到以下形式 [14]：

$${\mathbf{I}} - {\mathbf{Q}} = \left( \begin{array}{cccccc}
1 & 0 & 0 & \dots & 0 & 0 \\
-\frac{1}{2} & 1 & 0 & \dots & 0 & 0 \\
-\frac{1}{3} & -\frac{1}{3} & \ddots & \ddots & \vdots & \vdots \\
\vdots & \vdots & \ddots & \ddots & \ddots & \vdots \\
-\frac{1}{n-1} & -\frac{1}{n-1} & \ddots & -\frac{1}{n-1} & 1 & 0 \\
-\frac{1}{n} & -\frac{1}{n} & \dots & -\frac{1}{n} & -\frac{1}{n} & 1 \\
\end{array} \right).$$

$$(I - Q) = \begin{pmatrix}
1 & 0 & 0 & \cdots & 0 & 0 \\
-\frac{1}{2} & 1 & 0 & \cdots & 0 & 0 \\
-\frac{1}{3} & -\frac{1}{3} & \ddots & \ddots & \vdots & \vdots \\
\vdots & \vdots & \ddots & \ddots & \ddots & \vdots \\
-\frac{1}{n-1} & -\frac{1}{n-1} & \cdots & -\frac{1}{n-1} & 1 & 0 \\
-\frac{1}{n} & -\frac{1}{n} & \cdots & -\frac{1}{n} & -\frac{1}{n} & 1 \\
\end{pmatrix}.$$

（13.9）

基本矩阵 ${\boldsymbol P}$ 具有可归纳的结构，我们可以加以利用。具体而言，我们应用高斯-约旦消元法，先获得 $m \times m$ 的基本矩阵，再用该结果部分填充 $(m+1) \times (m+1)$ 基本矩阵的 $m \times m$ 子矩阵。从 $m=2$ 开始，迭代应用上述过程，我们得到如下形式 [14]：
$({\mathbf{I}} - {\mathbf{Q}})^{-1} = 
\left( 
\begin{array}{cccccc} 
1 & 0 & 0 & \dots & 0 & 0 \\ 
\frac{1}{2} & 1 & 0 & \dots & 0 & 0 \\ 
\frac{1}{2} & \frac{1}{3} & \ddots & \ddots & \vdots & \vdots \\ 
\vdots & \vdots & \ddots & \ddots & \ddots & \vdots \\ 
\frac{1}{2} & \frac{1}{3} & \ddots & \frac{1}{n-1} & 1 & 0 \\ 
\frac{1}{2} & \frac{1}{3} & \dots & \frac{1}{n-1} & \frac{1}{n} & 1 \\ 
\end{array} 
\right).$

递归结构可以如下展示。需要注意的是，随着基本矩阵的规模增加，将在底部添加一行、在左侧添加一列。可以观察到，$({\mathbf{I}} - {\mathbf{Q}})^{-1}$的第一列和最后一行的元素之和始终为$1/2$。对于$2 \times 2$矩阵，显然通过Gauss-Jordan消元运算得到$1/2$。当矩阵规模增加时：

\[
\begin{array}{lrcl}
3 \times 3: & \frac{1}{3} \times \frac{1}{2} + \frac{1}{3} \times 1 & = & \frac{1}{3} \times (\frac{3}{2}) \\
& & = & \frac{1}{2} \\
4 \times 4: & \frac{1}{4} \times \frac{1}{2} + \frac{1}{4} \times \frac{1}{2} + \frac{1}{4} \times 1 & = & \frac{1}{4} \times (\frac{4}{2}) \\
& & = & \frac{1}{2} \\
5 \times 5: & \frac{1}{5} \times \frac{1}{2} + \frac{1}{5} \times \frac{1}{2} + \frac{1}{5} \times \frac{1}{2} + \frac{1}{5} \times 1 & = & \frac{1}{5} \times (\frac{5}{2}) \\
& & = & \frac{1}{2}. \\
\end{array}
\]

对于$(n+1)\times(n+1)$的矩阵，可以推出

\[
\begin{array}{cl}
& \overbrace{\frac{1}{n+1} \times \frac{1}{2} + \frac{1}{n+1} \times \frac{1}{2} + \cdots + \frac{1}{n+1} \times \frac{1}{2}}^{n-1} + \frac{1}{n+1} \times 1 \\
= & \frac{1}{n+1} \times (\frac{n+1}{2}) \\
= & \frac{1}{2}. \\
\end{array}
\]

一个简单但繁琐的计算可以验证$({\mathbf{I}} - {\mathbf{Q}})({\mathbf{I}} - {\mathbf{Q}})^{-1} = {\mathbf{I}}$。我们应用引理13.3计算期望到达时间，得到[14]

\[
{\mathbf{h}} = ({\mathbf{I}} - {\mathbf{Q}})^{-1}\mathbf{1} = \left(1, \frac{1}{2} + 1, \frac{1}{2} + \frac{1}{3} + 1, \dots, \sum_{i=1}^n \frac{1}{i} + 1\right)^T.
\]
（13.11）

对于传递性锦标赛，$P^*$中的重新索引意味着${\mathbf{h}} = [h_i]_{i\in\mathcal{Q}} = [h_i]_{i=1}^n$，其中索引$i=1,\ldots,n$对应的是分数序列$(s_j : j = 1,\ldots,n)$逆序的顶点，分数是从小到大排列。因此，$[h_1] = [1]$就是顶点$s_n$（在根据分数序列原始索引下表示，分数为$n-1$）到主导顶点$s_{n+1}$（分数为$n$，在$n+1$阶锦标赛中）的期望到达时间$\mathbb{E}_{s_n}(\tau_{s_{n+1}})$。而从最弱顶点到主导顶点的期望到达时间为
\[
\mathbb{E}_{s_1}(\tau_{s_{n+1}}) = [h_n] = \left[ \sum_{i=1}^n \frac{1}{i} + 1 \right].
\]

尽管我们能够得出$(1+1)$CEA在$(n+1)$阶递归协同进化锦标赛上的解析表达式，并可通过直接计算得到${\mathbf{h}} = [h_i]_{i\in\mathcal{Q}} = [n,n,n,\ldots,n]^T$，但对吸收型$(1+1)$CEA的进一步理论研究揭示了一个令人惊讶的普遍性结论。相关内容在下述定理13.1中给出。其技术性证明利用了所有吸收型$(1+1)$CEA相关概率转移矩阵的独特不变量性质，即$Q$的全部对角元素值均为1。
2n−1n\frac{1}{2}\frac{n-1}{n}。

**定理13.1**（[14]）  
考虑在可约协同进化锦标赛 $T(V_{S(n+n_d)},A) \in \mathcal{T}(V_{S(n+n_d)})$ 上运行的 (1+1) 协同进化算法（CEA）。设 $V_{S(n+n_d)}$ 由两个不相交的子集 $V_{S(n)}^1$ 和 $V_{S(n_d)}^2$ 构成，其中 $V_{S(n)}^1 \cap V_{S(n_d)}^2 = \emptyset$。设 $V_{S(n_d)}^2 \Leftarrow V_{S(n)}^1$。令 CMC 的马尔可夫转移矩阵为其标准形式 ${\boldsymbol P}^{*}$。令 $i \in \mathcal{Q}$ 表示瞬态状态，且 $|\mathcal{Q}| = n$。定义列向量 ${\mathbf{h}} = [h_i]_{i \in \mathcal{Q}}$，其中 $h_i = {\mathbb{E}}_{i}(\tau_{\mathcal{A}})$ 表示从瞬态状态 $i$ 到吸收类 $\mathcal{A} = V_{S(n_d)}^2$ 的期望到达时间（hitting time）。则有  
${\mathbf{h}} = [h_i]_{i \in \mathcal{Q}} = [n, n, n, \ldots, n]^T$。

下面，我们将通过一个例子展示我们在[14]中为定理13.1给出的简化证明思路，以呈现相关的技术论证。考虑在一个未标记的协同进化锦标赛 $T(V_{S(n+1)},A)$ 上的 (1+1) CEA，该锦标赛有一个主导顶点和一个由单一受主导强分量组成的正则锦标赛，其中奇数阶 $n \geq (3 + 2k)$，$k = 0, 1, 2, \ldots$。

此时，我们有如下的马尔可夫转移矩阵：
$${\boldsymbol P}^{*} = \left(\begin{array}{cccccccc}
1 & 0 & 0 & 0 & 0 & \cdots & 0 & 0 \\
\frac{1}{n} & \frac{1}{2}\frac{n-1}{n} & 0 & \frac{1}{n} & 0 & \cdots & 0 & \frac{1}{n} \\
\frac{1}{n} & \frac{1}{n} & \frac{1}{2}\frac{n-1}{n} & 0 & \frac{1}{n} & \cdots & \frac{1}{n} & 0 \\
\frac{1}{n} & 0 & \frac{1}{n} & \frac{1}{2}\frac{n-1}{n} & 0 & \cdots & 0 & \frac{1}{n} \\
\vdots & \vdots & \vdots & \vdots & \ddots & \ddots & \vdots & \vdots \\
\vdots & \vdots & \vdots & \vdots & \vdots & \ddots & \ddots & \vdots \\
\frac{1}{n} & \frac{1}{n} & 0 & \frac{1}{n} & 0 & \cdots & \frac{1}{2}\frac{n-1}{n} & 0 \\
\frac{1}{n} & 0 & \frac{1}{n} & 0 & \frac{1}{n} & \cdots & \frac{1}{n} & \frac{1}{2}\frac{n-1}{n}
\end{array}\right)
$$

我们可以将 ${\mathbf{Q}}$ 用如下两个等式替换其中的元素 $a = \frac{1}{2}\frac{n-1}{n}$ 和 $b = \frac{1}{n}$，得到
$${\mathbf{Q}} = \left( \begin{array}{ccccccc}
a & 0 & b & 0 & \cdots & 0 & b \\
b & a & 0 & b & \cdots & b & 0 \\
0 & b & a & 0 & \cdots & 0 & b \\
\vdots & \vdots & \vdots & \ddots & \ddots & \vdots & \vdots \\
\vdots & \vdots & \vdots & \vdots & \ddots & \ddots & \vdots \\
b & 0 & b & 0 & \cdots & a & 0 \\
0 & b & 0 & b & \cdots & b & a
\end{array} \right)$$
我们观察到，$\mathbf{Q}$的第$(i+1)$行是第$i$行经过一次右移后的结果。此外，$\mathbf{Q}$的每一行的行和都相同，其值为 $a + \left(\frac{n-1}{2}\right)b = \frac{1}{2}\frac{n-1}{n} + \left(\frac{n-1}{2}\right)\frac{1}{n} = \frac{n-1}{n}$。这是由正则锦标赛（regular tournament）的定义所决定的，其中每个顶点$i$的入度和出度相同，即 $d_D^+(i) = d_D^-(i)$。进一步地，$\mathbf{Q}^{1}_{n \times n}\mathbf{1} = [\frac{n-1}{n}, \ldots, \frac{n-1}{n}]^T$。显然，也可以通过计算 $\mathbf{1} - \mathbf{R}\mathbf{1}$ 得到相同的结果，因为 $\mathbf{R} = [\frac{1}{n}, \ldots, \frac{1}{n}]^T$，且 $(\mathbf{R}\mathbf{1} + \mathbf{Q}\mathbf{1}) = \mathbf{1}$。

我们需要证明 $(\mathbf{I} - \mathbf{Q})^{-1} = \mathbf{I} + \mathbf{Q} + \mathbf{Q}^2 + \mathbf{Q}^3 + \cdots$ 的行和是几何级数。我们有 $\mathbf{Q}^m_{n \times n}\mathbf{1} = [Q_i^m: i = 1,\ldots,n]^T$, 其中 $Q_i^m = \left(\frac{n-1}{n}\right)^m$，可对$m$归纳证明。因此，$\mathbf{h} = (\mathbf{I} - \mathbf{Q})^{-1}\mathbf{1}$ 的每个元素都服从如下几何级数：$1 + \frac{n-1}{n} + (\frac{n-1}{n})^2 + (\frac{n-1}{n})^3 + \cdots = \frac{1}{1 - \frac{n-1}{n}} = n$，于是 $\mathbf{h} = [h_i]_{i \in \mathcal{Q}} = [n, n, n, \ldots, n]^T$。

我们可以验证上述计算，注意 $\mathbf{I} - \mathbf{Q} = \left( \begin{array}{ccccccc} 1 - a & 0 & -b & 0 & \cdots & 0 & -b \\ -b & 1 - a & 0 & -b & \cdots & -b & 0 \\ 0 & -b & 1 - a & 0 & \cdots & 0 & -b \\ \vdots & \vdots & \vdots & \ddots & \ddots & \vdots & \vdots \\ \vdots & \vdots & \vdots & \vdots & \ddots & \ddots & \vdots \\ -b & 0 & -b & 0 & \cdots & 1 - a & 0 \\ 0 & -b & 0 & -b & \cdots & -b & 1 - a \\ \end{array} \right)$。

$(\mathbf{I} - \mathbf{Q})\mathbf{h}$ 的第$i$个元素精确为1，即：
$$(1 - a)n + \left(\frac{n-1}{2}\right)(-b)n = \left(1 - \frac{1}{2}\frac{n-1}{n}\right)n + \left(\frac{n-1}{2}\right)(-\frac{1}{n})n = \frac{2n-n+1}{2} - \frac{n-1}{2} = \frac{2}{2} = 1$$

更重要的是，定理13.1和式(13.11)的理论结果为我们揭示了一个重要见解：在一个由其关联有向图表示的协同进化问题结构中，该结构将影响在此有向图上进行随机游走所建模的搜索过程。特别地，有向图中由每个顶点的出邻域 $N_D^+(v_i),~i=1,\ldots, n+1$ 所给定的标记信息可以被解释为一种领域知识。将这样的信息纳入协同进化搜索，将导致构造出具有快速搜索过程的特定CMC（即标准随机游走）。例如，在可递锦标赛（transitive tournaments）上的随机游走具有 $O(\ln n)$ 的时间复杂度来获得支配解 $v_{n+1}$ [20]。（注：从 $v_1$ 到 $v_{n+1}$ 的期望抵达时间呈调和级数 $H_n$，见式(13.11)）。定理13.1将为在可递锦标赛上的$(1+1)$协同进化算法（CEA）给出$O(n)$的时间复杂度。

然而，在CMC分析中只考虑从 $v_1$ 到 $v_{n+1}$ 的期望抵达时间并不能充分理解标准随机游走与 $(1+1)$ CEA 之间在速度上的差异。如果我们考虑从每个顶点 $v_1, \ldots, v_n$ 的期望抵达时间，那么各自的CMC找到 $v_{n+1}$ 的实际速度差异将会非常明显。特别地，对$(1+1)$ CEA而言，无论从 $v_1$ 还是 $v_n$ 开始，搜索过程的速度都没有变化。由于缺少形如 $N_D^+(v_i),~i=1,\ldots,n+1$ 的领域知识，（$1+1$）CEA 的行为就好像它仅仅在顶点集 $V_{S(n+1)}$ 上进行均匀采样。如果我们充分利用定理13.1的解释，那么在被支配集 $V_{S(n+1)}\setminus \{v_{n+1}\}$ 中采用不同结构的情况下，$(1+1)$ CEA 获得支配解的速度依然保持不变。换句话说，$(1+1)$ CEA 在任何可约的协同进化锦标赛上都执行的是一种盲目搜索。
13.2.3 协同进化马尔可夫链定量度量的计算研究

本节将以基于期望到达时间的定量度量为基础，对协同进化马尔可夫链（CMC）的分析进行计算研究。在前文中，我们能够推导出解析表达式，从而可以计算具有简单结构（如可传递锦标赛）的可约协同进化锦标赛特定CMC的到达时间。然而，实际中的可约协同进化锦标赛通常比传递性锦标赛更为复杂。因此，本计算研究的另一个方面是探索如何以可控方式生成随机协同进化锦标赛。这里的潜在应用是利用这些随机锦标赛作为基准，用于研究在更复杂且更现实的CMC中估计感兴趣量的其他方法。需要注意的是，简单地对底层图的边进行随机定向，通常会生成不可约的随机锦标赛 [48]。例如，阶数为5的随机锦标赛生成不可约结构的概率超过$50\%$，而对于阶数为16的随机锦标赛，该概率则超过$99.9\%$，此后生成的随机锦标赛几乎必然是不可约的。另一种方法是首先生成所有合法的得分序列 [32]，再据此重构相应的锦标赛 [48]。在这里，我们提出了一种通过单一控制参数，在不可约锦标赛与具有明显传递结构的锦标赛这两个复杂性极端之间生成随机协同进化锦标赛层次的方法。该方法受到复杂网络研究中通过优先连接的网络生长方法的启发 [58]。算法13.3中描述了一种通用的过程式实现。通过下述过程，可以以随机或任意方式（如无环排序）用一个初始小规模锦标赛作为种子来生成锦标赛：

算法13.3 基于生长的随机锦标赛生成算法 $T(V_{S(n)}, A)$ [14]
输入：共演锦标赛 $T(V, A)$，$n_{\mathrm{init}}$ 初始锦标赛的规模，$n_{\mathrm{final}}$ 最终锦标赛的规模  
输出：$T(V, A)$ 共演锦标赛 $T(V_S(n_{\mathrm{final}}), A)$ 

1: 过程 $\mathrm{GROW}(T, n_{\mathrm{init}}, n_{\mathrm{final}})$  
2: $n := 0$  
3: $T := \mathrm{initTournament}(n_{\mathrm{init}})$  // 初始化共演锦标赛 $T(V_S(n_{\mathrm{init}}), A)$  
4: $n := n_{\mathrm{init}}$  
5: 当 $n \leq n_{\mathrm{final}}$ 时循环  
6: $T_t := \mathrm{reindex}(T)$  // 对 $T_t(V_t, A_t)$，令 $V_t := V$，$A_t := A$，然后根据其得分序列 $(s_1, s_2, s_3, ..., s_n)$ 重新编号  
7: $p := \mathrm{attachPr}(T_t)$  // 计算 $V_S(n)$ 的连接概率，$p = \{P_i^-\}_{1}^{n}$  
8: $T.V := \mathrm{addVertex}(T.V_t, v_{n+1})$  // $V := V_t \cup v_{n+1}$  
9: $i := 0$  
10: 对所有 $i \leq n$ 进行循环  
11: $T.A := \mathrm{addArc}(T.A_t, v_{n+1}, v_i, p[i])$  // $A := A_t \cup (v_{n+1} \rightarrow v_i)$，概率为 $P_i^-$，否则 $A := A_t \cup (v_i \rightarrow v_{n+1})$  
12: $i := i + 1$  
13: 结束循环  
14: $n := n + 1$  
15: 结束循环  
16: 返回 $T$  
17: 结束过程  

$\mathrm{initTournament}(n_{\mathrm{init}})$。针对每一次迭代 $n_{\mathrm{init}} \leq n \leq n_{\mathrm{final}}$，在已有的规模为 $n$ 的锦标赛中添加一个新的节点（顶点）。然后，来自新节点的 $n + 1$ 条边将定向指向现有锦标赛中的节点。这个定向是随机且独立进行的，但可以通过受控的方式对得分较高的现有顶点施加倾向性。为此，需要计算连接概率
$$
\mathbb{P}_i^{-} = \frac{1}{1 + e^{-\beta(x_i - \overline{x})}}
$$
其中 $\{x_i\}_{1}^{n}$ 表示锦标赛的相关统计量。参数 $\beta > 0$ 控制伯努利分布中 S型连接函数的陡峭程度 [14]。较小的 $\beta$ 值会降低函数的陡峭程度，当 $\beta \to 0$ 时有 $\mathbb{P}_i^{-} \to 1/2$（这与上文讨论的随机定向情况相同）。该方法的动机来源于统计物理学中描述的现象。现有锦标赛的每个节点都具有与其得分（入度）相关的能量。新加入的节点 $n + 1$ 更可能将其边指向现有锦标赛中的高能量节点。参数 $\beta$ 起到逆温度的作用。在高温区，节点的能量趋于均等，在现有锦标赛中任意节点的边指向概率也是一样的。相反，当系统降温（增大 $\beta$）时，新节点 $n + 1$ 指向高能量节点的边定向趋势会变得更为明显。
以下的计算研究需要对共进化锦标赛中循环结构的复杂性进行度量。我们正式构造了如下指标。令一个阶数为 $n$ 的共进化锦标赛 $T(V_{S(n)}, A)$，其顶点依照得分序列 $(s_1, s_2, s_3, \ldots, s_n)$ 被标号为 $v_1, v_2, v_3, \ldots, v_n$。令两个整数 $i, j = 1, 2, 3, \ldots, n$。$T(V_{S(n)}, A)$ 关联有两个数列 $S = \left(s_1, \sum_{i=1}^2 s_i, \sum_{i=1}^3 s_i, \ldots, \sum_{i=1}^n s_i\right)$ 和 $L = \left(0, \binom{2}{2}, \binom{3}{2}, \ldots, \binom{n}{2}\right)$。我们用 $S_i$ 表示 $S$ 的第 $i$ 个元素（同理，$L$ 用 $L_i$）。然后，我们定义 Landau 指数 $\nu$ 如下 [14]：

\[
\nu = \sum_{i=1}^n \mathbb{I}\{S_i = L_i\}
\]
\[
\mathbb{I}\{S_i = L_i\} = \begin{cases}
1 & \text{若 } \sum_{j=1}^i s_j = \binom{i}{2} \\
0 & \text{若 } \sum_{j=1}^i s_j \neq \binom{i}{2}
\end{cases}
\]

该指数 $\nu$ 的取值范围为 $1$ 到 $n$。本质上，它统计了共进化锦标赛中强连通分量（strong components）的数量，并为共进化锦标赛中循环结构的复杂度提供了有用的宏观尺度指标。$\nu=1$ 表示该锦标赛为不可约锦标赛（irreducible tournament），而 $\nu=n$ 表示为传递性锦标赛（transitive tournament）。当 $\nu$ 靠近 $n$ 时，表示锦标赛中存在许多较短的循环，可以分别归类为不同的子锦标赛。当 $\nu$ 接近 $1$ 时，锦标赛中包含长度在 $n/2$ 到 $n-1$ 之间的较长循环。预先计算好 $L$ 后，该指标的计算复杂度与锦标赛规模线性相关。我们在文献 [14] 中已进行了初步实验，以确定相关参数的合适取值。
我们在初始化图增长时使用传递性锦标赛。对于统计量 $\{x_i\}_1^n$，我们使用节点排名 $(i)_1^n$ 且 $\overline{x} = n/2$。我们有 $\beta \in [0,4]$。最终阶数 $n_{\mathrm{final}} = 1000$ 的锦标赛由初始阶数 $n_{\mathrm{init}} \in \{10, 100\}$ 的锦标赛生成。每组实验重复 $100$ 次，因此随机且独立地生成 $100$ 个锦标赛。结果汇总于图 13.2，其中，图 13.2b 显示，当过程 GROW 以更大的初始锦标赛启动时，生成的锦标赛拥有更多的强连通分量。
接下来，我们在 1000 节点的锦标赛中，再引入一个额外的主导节点（第 1001 个节点），该节点从原有的 1000 个顶点接收所有指向自身的有向边。此结构将锦标赛变为可约（reducible）。随后我们计算从分数最低的节点 $x_1$ 到主导（吸收）节点 $x_{1001}$ 的期望到达时间（expected hitting time）。这样可以研究 Landau 指数 $\nu$ 所反映的传递性结构增加时，对于共进化搜索寻找到主导节点速度的影响。对于共进化，我们采用标准的随机游走（random walk）。我们观察到的特征如下：

0
20
40
60
80
100
图13.2的柱状图展示了在生成具有特定强连通分量数量（即Landau指数）的随机锦标赛实验中，不同实验运行次数的分布情况。其中（a）和（b）分别对应初始传递锦标赛的阶数为10和100的实验集[14]。横轴表示强连通分量的数量，纵轴表示实验运行次数，$B$取不同值：$B=0.1$、$B=0.5$、$B=1.0$、$B=2.0$。

图13.3展示了箱线图，其中中线表示中位数，两端分别对应第25百分位数和第75百分位数，离群点以“+”表示。参数$\beta$（逆温度）分别在(a) $[0,2]$和(b) $[1,4]$的区间内取值[14]，用于表示不同的高温和低温区间。

在图13.3(a)所示的高温区间（$\beta\approx[0,0.001]$），所有1000节点的基础锦标赛均为不可约类型，因此仅存在一个强连通分量。再加上主导的第1001个节点，总共存在两个强连通分量。此时击中时间（hitting time）的波动极小（范围约为492–497），这反映了被主导的1000节点不可约基础锦标赛内在结构的有限变化。作为对比，普通999节点基础锦标赛的期望击中时间恰好为500。

在低温区间（$\beta\approx[1,4]$）下，显著的传递结构主导基础锦标赛，击中时间的范围为$[7.486,10.659]$。作为对比，如果基础锦标赛为节点全游（pancyclic）结构（即具有最长传递结构），则其期望击中时间略小于10；对于传递结构的1000节点基础锦标赛，期望击中时间约为7.4855。

与大多数复杂系统类似，在这两个温度区间之间，还存在一个中等温度区间，该区间表现出最有趣的行为。此时可以观察到Landau指数和期望击中时间的全部可能取值范围。实际上，这一情况对应于传递结构与多样的击中时间共存，反映了强连通分量中循环结构程度的不同。图13.3展示了（逆）温度$\beta$与主导的第1001个节点击中时间之间的关系。如预期所示，极高温度会导致基础锦标赛不可约化，形成复杂的循环结构并带来较长的击中时间，而极低温区间则导致基础锦标赛呈现显著而平凡的传递结构，从而得到最短的击中时间。此外，可以预期（从生成大量随机锦标赛的角度来看），随着温度的降低，基础锦标赛中的传递结构逐渐增强，期望击中时间将随之缩短。最复杂和有趣的行为区间位于这两极之间，即$\beta\approx(0.005,0.5)$。
我们用$(1+1)$协同进化算法（CEA）重复了实验，并获得了验证第13.1定理所表述结果的数据。期望命中时间的计算依赖于与当前研究的协同进化马尔可夫链（CMC）相关的概率转移矩阵。协同进化搜索过程的不同构建方式可能导致期望命中时间分析出现不同的结果。期望命中时间具有实用的解释意义，即可作为协同进化搜索支配解速度的度量。然而，要合理利用这一定量指标，需事先获得这些有向图基础结构的定性知识。虽然可以在一个不可约有向图中形式化地以起始顶点和终止顶点来定义命中时间，但在这里我们将命中时间提出作为衡量协同进化搜索支配解速度的度量。只有当CMC作用在可约的锦标赛图（tournaments）上时，这一度量才具有有意义的解释。

13.3 竞争性协同进化中解性能和访问概率的对偶性

在本节中，我们呈现了一项关于CMC与PageRank之间深层联系的研究 [13]，以应对抽象协同进化的定量分析需要关于有向图基础结构的先验定性知识的问题。我们发展了一种有理论基础的方法，用于在给定的协同进化有向图中度量和排序解（顶点）的性能（重要性）。在PageRank形式主义中，若$u$支配$v$，即$v \rightarrow u \in A$，那么$v$会将其部分权威传递给$u$，其中$D = (V, A)$为有向图。因此，PageRank权威值表示顶点的重要性。我们将在第13.3.1节中，于竞争性协同进化的设定下形式化PageRank权威值。随后在第13.3.2节，我们将形式上证明，通过适当归一化后的PageRank权威值，与协同进化随机游走（带重启）在有向图上的长期访问概率具有自然的对偶解释。我们还可以精确量化当重启概率变化时这些权威值的变化。最后，在第13.3.3节，我们将进行计算实验，以展示如何利用PageRank权威值对具有不同基础结构的协同进化有向图进行特征刻画。
13.3.1 不变量测度与抽象协同进化的PageRank权威性

在13.2节中，我们已经引入并研究了从任一瞬态状态 $x \in \mathcal{Q}$ 出发至吸收类 $\tau_{\mathcal{A}}$ 的期望命中时间 $\mathbb{E}_x(\tau_{\mathcal{A}})$，作为分析可约马尔可夫链（CMC）的定量度量。瞬态状态对应于被支配的顶点，而吸收类对应于协同进化有向图中的主导子集。那么，是否存在适用于可约和不可约有向图CMC的相关定量表征？我们首先引入马尔可夫链上的不变量测度的概念[50]。

设行向量 $\boldsymbol{\pi} = (\pi_x : x \in X)$ 为 $X$ 上的概率测度，$\boldsymbol{P}$ 为马尔可夫转移矩阵。如果 $\boldsymbol{\pi}$ 满足 $\boldsymbol{\pi}\boldsymbol{P} = \boldsymbol{\pi}$，则称 $\boldsymbol{\pi}$ 为不变量。定义在 $X$ 上且满足该不变量性质的 $\boldsymbol{\pi}$ 便代表了相关马尔可夫链的长期极限分布。

在抽象协同进化的背景下[14]，此类分布 $\boldsymbol{\pi}$ 可以等价表述如下：

(i) 平稳性：若 $(\varPhi_n)_{n \geq 0}$ 是满足 $\boldsymbol{\mu} = \boldsymbol{\pi}$ 和 $\boldsymbol{P}$ 的CMC，则对一切 $x \in X$，$(\varPhi_{m+n})_{n \geq 0}$ 同样是CMC，且有 $\mathbb{P}(\varPhi_m = x) = (\boldsymbol{\pi} \boldsymbol{P}^m)_x = \pi_x$（见[50]中的定理1.7.1）。

(ii) 平衡性：假设协同进化有向图完全连通，则对所有 $x, y \in X$，有 $\mathbb{P}^n(x, y) \rightarrow \pi_y$，当 $n \rightarrow \infty$ 时（见[50]中的定理1.7.2）。

已知，对于定义在可数状态空间$X$上的不可约马尔可夫链$\boldsymbol{\varPhi}$，其平稳分布可由期望回归时间 $\mathbb{E}_x(\varrho_x)$ 计算得到[42]。其基本直觉在于，将 $\boldsymbol{\varPhi}$ 的轨迹按到某一任意状态的来回访问划分为同分布片段，每段在状态 $x \in X$ 上的平均时间占比，会近似等于 $\boldsymbol{\varPhi}$ 长期在 $x$ 上的占比。此外，无论初始分布 $\boldsymbol{\mu}$ 如何，链最终会收敛到其平稳分布[42, 50]。

在 $D_C$ 的全局连通结构的定性与 $\boldsymbol{\varPhi}$ 的定量表征之间存在着紧密联系[14]。在可约的 $D_C$ 上的随机游走会导致吸收性的 $\boldsymbol{\varPhi}$，而不可约的 $\boldsymbol{\varPhi}$ 则作用于强连通的 $D_C$。对于 $\boldsymbol{\varPhi}$，可以获得满足不变量性质 $\boldsymbol{\pi} \boldsymbol{P} = \boldsymbol{\pi}$ 的长期极限分布。对于吸收性的 $\boldsymbol{\varPhi}$，概率质量集中在吸收类，其它所有状态的概率均为零[47]。

针对CMC，我们有如下结果。
引理13.4（[14]）设有不可约的CMC，其马尔可夫转移矩阵为${\boldsymbol P}$。则在$X$上存在唯一的不变概率分布$\boldsymbol \pi = (\pi_x : x \in X)$，其中$0 \leq \pi_x \leq 1$且$\sum_{x \in X}\pi_x = 1$，满足$\boldsymbol \pi{\boldsymbol P} = \boldsymbol \pi$。此外，$\pi_x = 1/\mathbb{E}_x(\varrho_x)$。

引理13.5（[14]）设有可约的CMC，其马尔可夫转移矩阵为$\boldsymbol P$。对于所有非本质状态$y \in X$，有$\pi_y = 0$。此时，存在唯一的不变概率分布$\boldsymbol \pi$，该分布集中于$X$中本质状态的吸收类$C$。尽管CMC的平稳分布可以被形式化，但相关的定量分析与后续结论的解释，需要事先了解有向图的基础结构。例如，若CMC作用于可约有向图，则全部概率质量将集中于其吸收类上，继而导致对于支配子集顶点的结构没有任何了解。在此，我们对不变测度$\boldsymbol \pi$进行了细化，使其能够体现有关协同进化有向图结构的某些信息。这意味着要在抽象协同进化与PageRank（最初用于衡量网页重要性[6, 40]）之间建立直接联系，而PageRank正是大型网络分析方法之一[27]。在PageRank形式中，策略在成对交互时会将其一部分权威传递给支配它们的策略，这些权威值可以反映策略的性能。权威值与有向图的结构密切相关。我们利用有向图理论的符号对PageRank权威值[5]重表述，以明确其与有向图结构的关联。

设$u, v \in V_S$，则$u$的入邻域定义为$N^{-}_D(u) = \{v \in V_S\backslash\{u\} : vu \in A_R\}$，而$v$的出度（即发出的弧的个数）定义为$d^+_D(v) = |\{(u,v) \in A_R : v \in V_S\backslash\{u\}\}|$。$u$的PageRank权威值为
\[
\varphi_u = \alpha + (1 - \alpha) \sum_{v \in N^{-}_D(u)} \frac{\varphi_v}{d^+_D(v)},
\]
其中$\alpha \in (0,1)$。参数$\alpha$的含义可根据具体上下文调整，这里用于以凸组合的方式，将每个顶点（节点）固有的自由权威与其基于连通结构获得的权威相结合，从而计算PageRank权威值[5]。

设协同进化有向图$D_C(V_{S(n)})$的顶点标记为$v_1, v_2, v_3, \ldots, v_n$，则$D_C(V_{S(n)})$的PageRank权威值向量为
\[
\boldsymbol \varphi = \alpha \mathbf{e} + (1 - \alpha) \boldsymbol \varphi \mathbf{M}, \tag{13.15}
\]
其中$\mathbf{e} = [1, 1, 1, \ldots, 1]$，矩阵$\mathbf{M} = (m_{ij} : i, j \in \{1,2,3,\ldots,n\})$，且
\[
m_{ij} = 
\begin{cases}
\frac{1}{d^+_D(v_i)} & \text{若 } v_j \in N^{-}_D(v_i) \\
0 & \text{否则}
\end{cases}
\]
$m_{ij} = \frac{1}{d_D^+(v_i)}$ 当 $i \rightarrow j$，$m_{ij} = 0$ 当 $i \nrightarrow j$。$\boldsymbol\varphi = [\varphi_1, \ldots, \varphi_n]$（这里我们用简化记号$\varphi_i$代替每个$\varphi_{v_i}$，表示$\boldsymbol\varphi$的第$i$个分量）[13]。有向图结构的影响体现在$\mathbf{M}$中非零的元素上，这对应于与$D_C(V_{S(n)})$相关的邻接矩阵。

#### 13.3.2 PageRank与协同进化马尔科夫链的联系

实际的PageRank计算需要将其转化为特征系统或等价的线性系统，前提是$\boldsymbol\psi$经过适当归一化（$\boldsymbol\psi\mathbf{e}^{\mathsf{T}} = 1$）[27]。PageRank可以被视为一个马尔可夫链，其转移概率矩阵$\boldsymbol{P}_{\mathrm{PR}}$是原始的概率转移矩阵（非负、不可约，且其谱半径$\rho(\boldsymbol{P}_{\mathrm{PR}})$是其唯一的模为$r$的特征值[46]）。常用的PageRank修正规则为$\boldsymbol{P}_{\mathrm{PR}} = \alpha \mathbf{e}^{\mathsf{T}} \mathbf{s} + (1 - \alpha) \boldsymbol{P}$，其中$\mathbf{s}$为一般概率向量（如均匀分布$\frac{1}{n}\mathbf{e}$）[41]，这样任意CMC都变成不可约（例如，原本为吸收型CMC也变成不可约），并且形成与$D_C(V_{S(n)})$上的PageRank类似的马尔可夫链。关键的是，这一设定可以解释为以概率$\alpha$在CMC中引入重启。正是这种简单的改进，使得PageRank权威度可以适用于所有协同进化有向图（即既可约也不可约）。

接下来，$\boldsymbol{P}$的设置指示使用哪种协同进化搜索过程，例如$\boldsymbol{P} = \mathbf{M}$意味着CMC为标准的随机游走。需要注意的是，协同进化有向图的PageRank权威度$\boldsymbol\psi$是方程$\boldsymbol\psi\boldsymbol{P}_{\mathrm{PR}} = \boldsymbol\psi$的解。这为$\boldsymbol\psi$提供了两个自然的解释：（1）它指示协同进化如何度量个体$v\in V_{S(n)}$的重要性（权威度）并用于排序；（2）它给出了由特定CMC建模的协同进化搜索在有向图上的长期访问概率分布。

下文将介绍我们在[13]中建立的、连接CMC与PageRank的多个理论结果，这些结果反映在与CMC和PageRank（带重启概率$\alpha$的CMC）相关的平稳分布$\boldsymbol\pi$和$\boldsymbol\psi_\alpha$的关系上[37]。我们首先给出了与CMC相关的PageRank向量存在性的保证，摘要如下引理13.6所示。技术证明涉及对$\big(\beta \mathbf{I} + (\mathbf{I} - \boldsymbol{P})\big)^{-1}$的直接计算，这利用了“惰性随机游走”。这种方法将概率转移矩阵修改为$\boldsymbol{Z} = (\mathbf{I} + \boldsymbol{P})/2$，但仍保持原来的平稳分布$\boldsymbol\pi$，因为$\boldsymbol\pi = \boldsymbol\pi\boldsymbol{Z}$。更重要的是，惰性随机游走器的这种修改确保$\boldsymbol{Z}$是无周期的[16]。这使我们能够宣称，对于一个在强连通有向图$D_C(V_{S(n)})$上运行、且无周期的CMC，其平稳分布$\boldsymbol\pi$就是极限不变分布，这可以由Perron-Frobenius定理推出。最后，$\mathbf{W} = \beta\mathbf{I} + (\mathbf{I} - \boldsymbol{P})$是严格占优的对角矩阵，因此可逆[35]。

**引理13.6（[13]）**  
设$D_C(V_{S(n)}) \in \mathcal{D}_C(V_{S(n)})$为协同进化有向图。令$\boldsymbol{P}$为在$D_C(V_{S(n)})$上运行的CMC对应的概率转移矩阵，$\boldsymbol{Z} = (\mathbf{I} + \boldsymbol{P})/2$为其惰性版本。行向量$\mathbf{s}$是$V_{S(n)}$顶点集上的概率分布。给定$\alpha \in (0,1)$，以及$\beta = \frac{2\alpha}{1 - \alpha}$时，个性化PageRank向量是如下线性系统的唯一解：
${\boldsymbol \psi}_{\alpha}({\mathbf{s}}) = \alpha{\mathbf{s}} + (1 - \alpha){\boldsymbol \psi}_{\alpha}({\mathbf{s}}){\boldsymbol Z}$

并且可以被计算为
$${\boldsymbol \psi}_{\alpha}({\mathbf{s}}) = \beta{\mathbf{s}}\left(\beta{\mathbf{I}} + ({\mathbf{I}} - {\boldsymbol P})\right)^{-1}.$$

接下来，我们采用线性系统的表述来揭示关于系数矩阵${\mathbf{W}}$的若干性质，并在引理13.7中进行总结。设${\mathbf{R}}^{n \times n}$为$n \times n$的实方阵集合。若${\mathbf{A}} \in {\mathbf{R}}^{n \times n}$是一个M-矩阵，其形式为${\mathbf{A}} = c{\mathbf{I}} - {\mathbf{B}}$，其中${\mathbf{B}} \geq {\mathbf{0}} = (b_{ij} \geq 0 : 1 \leq i, j \leq n)$且$c \geq \rho({\mathbf{B}})$，$\rho(\cdot)$表示谱半径[4]。$||{\mathbf{A}}||_{\infty}$为${\mathbf{A}}$的$\infty$-范数，而$\kappa_{\infty}({\mathbf{A}}) = ||{\mathbf{W}}||_{\infty} \,||{\mathbf{W}}^{-1}||_{\infty}$是针对${\mathbf{A}}$的$\infty$-范数条件数[46]。技术性证明的主要思想利用了${\mathbf{W}} = c{\mathbf{I}} - {\mathbf{B}}$是一个M-矩阵，这里$c = 1 + \beta$且${\mathbf{B}} = {\boldsymbol P}$为一个随机矩阵。因此，${\mathbf{W}}$是逆正的，即${\mathbf{W}}^{-1}$存在且${\mathbf{W}}^{-1} \geq {\mathbf{0}}$（见定理2.3，[4]）。

引理13.7（[13]） 设${\boldsymbol \psi}_{\alpha}({\mathbf{s}})(\beta{\mathbf{I}} + ({\mathbf{I}} - {\boldsymbol P})) = \beta{\mathbf{s}}$为PageRank问题的线性系统表述。系数矩阵${\mathbf{W}} = \beta{\mathbf{I}} + ({\mathbf{I}} - {\boldsymbol P})$具有如下性质：

1. ${\mathbf{W}}$为M-矩阵。
2. ${\mathbf{W}}$非奇异。
3. ${\mathbf{W}}$的行和为$\beta$。
4. $||{\mathbf{W}}||_{\infty} = 2 + \beta$。
5. ${\mathbf{W}}^{-1} \geq {\mathbf{0}}$。
6. ${\mathbf{W}}^{-1}$的行和为$\frac{1}{\beta}$。
7. $||{\mathbf{W}}^{-1}||_{\infty} = \frac{1}{\beta}$。
8. $\kappa_{\infty}({\mathbf{W}}) = \frac{2 + \beta}{\beta} = \frac{1}{\alpha}$。

引理13.7中的$\kappa_{\infty}({\mathbf{W}})$指示了PageRank解对${\mathbf{W}}$扰动的敏感性。当采用基于高斯消元的直接计算法时，它量化了${\mathbf{W}}$相对于机器精度的病态程度。对于机器精度$m_p$，该量取为$m_p\kappa_{\infty}({\mathbf{W}})$（详见第3章，[28]）。例如，如果$m_p\kappa_{\infty}({\mathbf{W}}) \leq 1$，对于$\alpha = 0.1$，则表明只会丢失一位有效数字的精度。在实际应用中确定$\kappa_{\infty}({\mathbf{W}})$非常重要，因为PageRank的计算通常采用迭代方法（有关这些方法的误差特性，见定理2.2 [27]）。

现在，我们可以继续研究，具体地说，对CMC引入重启并改变$\alpha$，将如何影响CMC的长期极限分布。我们首先给出定性结果，具体陈述在定理13.2中。其技术性证明采用直接计算，感兴趣的读者可参考[13]附录部分。该证明应用了与随机游走的有向图Laplacian相关的结果，这些结果已被推广到有向图[43]，我们在下文中正式定义。

设${\boldsymbol P}$为作用于不可约$D_C(V_{S(n)})$上的CMC关联的转移矩阵，其平稳向量为${\boldsymbol \pi} = (\pi_1, \pi_2, \pi_3, \ldots, \pi_n)$，满足${\boldsymbol \pi}{\boldsymbol P} = {\boldsymbol \pi}$。设对角矩阵为${\boldsymbol \varPi} = \mathrm{diag}(\pi_i)$。归一化有向图Laplacian定义为
$$\widetilde{L} = {\boldsymbol \varPi}^{1/2}({\mathbf{I}} - {\boldsymbol P}){\boldsymbol \varPi}^{-1/2}$$
$2 \tilde{\boldsymbol{\mathcal{L}}} = {\boldsymbol\varPi}^{\frac{1}{2}}({\mathbf{I}} - {\boldsymbol P}){\boldsymbol \varPi}^{-\frac{1}{2}}$，其各元素定义为

$$
\tilde{\mathcal{L}}_{ij} = \left\{ 
\begin{array}{ll} 
1 - p_{ii} & \text{如果 } i = j \\ 
-\pi_i^{\frac{1}{2}} p_{ij} \pi_j^{\frac{1}{2}} & \text{如果 } (i,j) \in A \\ 
0 & \text{否则} 
\end{array} 
\right.
$$

有向图的格林函数为 $\tilde{\boldsymbol{\mathcal{Z}}} = \tilde{\boldsymbol{\mathcal{L}}}^+$，即 $\tilde{\boldsymbol{\mathcal{L}}}$ 的 Moore-Penrose 广义逆，其满足 $\tilde{\boldsymbol{\mathcal{Z}}} \tilde{\boldsymbol{\mathcal{L}}} = \tilde{\boldsymbol{\mathcal{L}}} \tilde{\boldsymbol{\mathcal{Z}}} = \mathbf{I} - \tilde{\boldsymbol{\mathcal{J}}}$，其中 $\tilde{\boldsymbol{\mathcal{J}}} = ({\boldsymbol \pi}^{\frac{1}{2}})^T {\boldsymbol \pi}^{\frac{1}{2}}$。

定理13.2（[13]）：设${\boldsymbol \pi}$为在不可约DC($V_{S(n)}$)上的CMC的平稳分布向量。个性化PageRank向量表达如下：

$${\boldsymbol \psi}_{\alpha}({\mathbf{s}}) = {\mathbf{s}} \left( {\mathbf{I}} - {\boldsymbol \varPi}^{-\frac{1}{2}}\tilde{\boldsymbol{\mathcal{L}}} ({\beta} {\mathbf{I}} + \tilde{\boldsymbol{\mathcal{L}}})^{-1} {\boldsymbol \varPi}^{\frac{1}{2}} \right).$$

同时有：
$$
\lim\limits_{\beta \rightarrow 0} \tilde{\boldsymbol{\mathcal{L}}}({\beta}{\mathbf{I}} + \tilde{\boldsymbol{\mathcal{L}}})^{-1} = \tilde{\boldsymbol{\mathcal{L}}}(\tilde{\boldsymbol{\mathcal{L}}})^+ = \tilde{\boldsymbol{\mathcal{L}}}\tilde{\boldsymbol{\mathcal{Z}}}
$$

其中，$\beta = \frac{2\alpha}{1 - \alpha}$，而$\tilde{\boldsymbol{\mathcal{Z}}} = \tilde{\boldsymbol{\mathcal{L}}}^+$是$\tilde{\boldsymbol{\mathcal{L}}}$的Moore-Penrose 广义逆。

当$\alpha$取较小值（即$\beta$较小）时，PageRank向量近似等于平稳分布向量，即${\boldsymbol \psi}_{\alpha}({\mathbf{s}}) \approx {\boldsymbol \pi}$，对于任意$\mathbf{s}$，总有等号成立（即${\boldsymbol \psi}_0({\mathbf{s}}) = {\boldsymbol \pi}$）。

有向图拉普拉斯算子明确建立了DC结构及DC上随机游走相关量之间的联系。定理13.2展示了在CMC中引入重启后的平稳分布${\boldsymbol \pi}$与PageRank 向量${\boldsymbol \psi}_{\alpha}({\mathbf{s}})$之间的直接关系，该关系通过有向图拉普拉斯算子$\tilde{\boldsymbol{\mathcal{L}}}$及其广义逆——即格林函数$\tilde{\boldsymbol{\mathcal{Z}}}$来表达。后者还包含了重启概率$\alpha$的信息。非正式而言，定理13.2表明， 当$\alpha > 0$且很小时，随机游走者长时间停留在各顶点的概率将根据摄动过的归一化有向图拉普拉斯算子$(\beta \mathbf{I} + \tilde{\boldsymbol{\mathcal{L}}})^{-1}$被重新分配。

接下来的引理13.8定量给出了由于以概率$\alpha$重启引起的${\boldsymbol \pi}$与${\boldsymbol \psi}_{\alpha}({\mathbf{s}})$之间的具体差异，并给出了差异的闭式表达。其直接计算的证明利用了有限马尔可夫链扰动理论（具体为[59]中的定理1和2）。然而，随后的推论13.1和13.2更具实际意义。第一个结果声明，方程(13.19)中表述的差异随$\alpha$单调递增，即$\|{\boldsymbol \pi} - {\boldsymbol \psi}_{\alpha}({\mathbf{s}})\|$随$\alpha$的增加而增加。第二个结果则给出了该差异富有直观意义的上界。

最大差异发生于重启将${\boldsymbol \pi}$重新分配为均匀分布$s = \mathbf{1}$的情形。
$\mathbf{s} = \frac{1}{n}\mathbf{e}$（注意，我们采用了中心性度量的方法 [27]，并在整体有向图分析中取 $\mathbf{s}$ 为均匀分布）。

引理 13.8（[13]） 设一个在不可约 $\mathrm{DC}(V_{S(n)}) \in \mathcal{D}_C(V_{S(n)})$ 上运行的 CMC 与概率转移矩阵 $\mathbf{P}$、平稳分布 $\boldsymbol{\pi}$ 以及基本矩阵 $Z$ 对应。类似地，$\mathrm{DC}(V_{S(n)})$ 上的个性化 PageRank 是具有概率转移矩阵 $\mathbf{P}_{\mathrm{PR}}$ 和平稳分布 $\boldsymbol{\psi}_{\alpha}(\mathbf{s})$ 的 CMC。则有
\[
\begin{array}{rcl}
\boldsymbol{\pi} - \boldsymbol{\psi}_{\alpha}(\mathbf{s}) & = & -\boldsymbol{\psi}_{\alpha}(\mathbf{s})\big(\mathbf{I} - \mathbf{P}\big)Z \\
 & = & \beta(\boldsymbol{\pi} - \mathbf{s})\big(\beta \mathbf{I} + (\mathbf{I} - \mathbf{P})\big)^{-1}.
\end{array}
\]
即
\[
\boldsymbol{\pi} - \boldsymbol{\psi}_{\alpha}(\mathbf{s}) = -\boldsymbol{\psi}_{\alpha}(\mathbf{s})(\mathbf{I} - \mathbf{P})Z = \beta(\boldsymbol{\pi} - \mathbf{s})(\beta \mathbf{I} + (\mathbf{I} - \mathbf{P}))^{-1}.
\]

推论 13.1（[13]） 对于在不可约 $\mathrm{DC}(V_{S(n)}) \in \mathcal{D}_C(V_{S(n)})$ 上运行的 CMC 所对应的 $\boldsymbol{\pi}$ 和 $\boldsymbol{\psi}_{\alpha}(\mathbf{s})$，对于 $0 < \alpha_1 \leq \alpha_2 < 1$，有如下不等式：
\[
\left\|\boldsymbol{\pi} - \boldsymbol{\psi}_{\alpha_1}(\mathbf{s})\right\| \leq \left\|\boldsymbol{\pi} - \boldsymbol{\psi}_{\alpha_2}(\mathbf{s})\right\|
\]
其中 $\alpha_1, \alpha_2 \in (0, 1)$。

推论 13.2（[13]） 对于每一个不可约的 $\mathrm{DC}(V_{S(n)}) \in \mathcal{D}_C(V_{S(n)})$，其对应的 CMC，在 $\alpha \in (0,1)$ 下具有平稳分布 $\boldsymbol{\psi}_{\alpha}(\mathbf{s})$ 和 $\boldsymbol{\pi}$，满足如下扰动界：
\[
\left\|\boldsymbol{\pi} - \boldsymbol{\psi}_{\alpha}(\mathbf{s})\right\|_{\infty} \leq \left\|\boldsymbol{\pi} - \mathbf{s}\right\|_{\infty}
\]
当 $\alpha = 1$ 时取等号。

尽管刻画 $\left\|\boldsymbol{\pi} - \boldsymbol{\psi}_{\alpha}(\mathbf{s})\right\|_1$ 随 $\alpha$ 的变化更为有用，但获得紧的界限较为困难（如 $2/\alpha$ [5, 41, 49]）。我们使用耦合方法 [49] 结合有向图理论的论证，将界限改进为
\[
\left\|\boldsymbol{\pi} - \boldsymbol{\psi}_{\alpha}(\mathbf{s})\right\|_1 \leq \frac{2}{1 - \alpha}\left(1 - \left(\pi_{v_n} - \frac{1}{n}\right)\right)
\]
对于秩为 $n$ 且具有特定环结构的不可约 $\mathrm{DC}$ 成立（详细内容见 [13] 的附录），其中 $v_n$ 具有最高分数。

设 $\mathrm{DC}$ 是具有最少 3-回路的顶点全回路图 [10]。其顶点可划分为三个互不相交的子集 $\{v_1\}$, $V_{S(n-2)}$, 和 $\{v_n\}$。注意到子锦标赛 $T(V_{S(n-2)})$ 是传递的。则 $v_1$ 和 $v_n$ 的长期极限访问概率相同，且为 [13]
\[
\pi_{v_n} = \pi_{v_1} = \frac{1}{2}\left(1 - \frac{\mathbb{E}(\eta_{V_{S(n-2)}})}{2 + \mathbb{E}(\eta_{V_{S(n-2)}})}\right)
\]
其中
\[
\mathbb{E}(\eta_{V_{S(n-2)}}) = \frac{1}{n - 2} \sum_{m=1}^{n-2} \mathbb{E}_{x_m}(\tau_{v_n})
\]
$n-2 \mathbb{E}_{m=1} \mathbb{E}_{x_m}(\tau_{v_n}),\quad (13.23)$

其中，$\mathbb{E}(\eta_{V_{S(n-2)}})$ 表示随机游走者在 $V_{S(n-2)}$ 中平均花费的时间比例（即期望），此游走者重复经历 $v_1 \leadsto V_{S(n-2)} \leadsto v_n \leadsto v_1$ 这样的循环。因此，$\mathbb{E}(\eta_{V_{S(n-2)}})$ 其实就是从 $V_{S(n-2)}$ 中每个顶点出发到 $v_n$ 的期望命中时长（hitting time）的平均值，即：

\[
\begin{align*}
\mathbb{E}_{v_{n-1}}(\tau_{v_n}) &= \mathbb{E}_{x_1}(\tau_{v_n}) = H_1 = 1 \\
\mathbb{E}_{v_{n-2}}(\tau_{v_n}) &= \mathbb{E}_{x_2}(\tau_{v_n}) = H_2 = 1 + \frac{1}{2} \\
&\ \vdots \\
\mathbb{E}_{v_2}(\tau_{v_n}) &= \mathbb{E}_{x_{n-2}}(\tau_{v_n}) = H_{n-2} = 1 + \frac{1}{2} + \frac{1}{3} + \cdots + \frac{1}{n-2}
\end{align*}
\]

其中 $H_n = \sum_{k=1}^n \frac{1}{k}$，即第 $n$ 个调和数（Harmonic number）。

对于与有标签（同构）锦标赛上的随机游走相关的概率向量，还可以获得 $\|{\boldsymbol \pi} - {\boldsymbol \psi}_{\alpha}({\mathbf{s}})\|_1$ 的其它一般性上界 [48, 61]，这些上界将其单向的和定向对偶性 [3, 9] 考虑在内，我们的表述如下 [13]：

\[
\|{\boldsymbol \pi}_1 - {\boldsymbol \pi}_2\|_1 \leq 2\Big(1 - \frac{1}{n}\Big)
\]

即 $\|{\boldsymbol \pi}_1 - {\boldsymbol \pi}_2\|_1 \leq 2(1 - \frac{1}{n})$。

#### 13.3.3 PageRank 竞争共进化中的计算研究

在本节，我们介绍在 [13] 中针对在具有特定结构 [3, 14, 48] 的共进化锦标赛 $T_C(V_{S(n)})$ 上运行的 CMC（循环马尔可夫链）所进行的系统计算实验。这些锦标赛的选取基于已知的内部环路结构，这些结构具有代表性并有助于我们更好地利用先前发展出的定量指标，评估不同结构对 CMC 访问顶点概率的影响。同时，我们还计算了 PageRank 排序 [30] 的差异，后续将作详细说明。我们采用标准的幂法（Power Method）[41] 来计算 PageRank 向量 ${\boldsymbol \psi}_{\alpha}({\mathbf{s}})$，其中 CMC 使用带有均匀跳转向量（teleportation vector）${\mathbf{s}}$ 的“懒惰”版本。

在本次研究中，我们采用以下几类共进化锦标赛 $T_C(V_{S(n)})$：（a）不可约型（Irreducible）：选用具有最少 3-环的全环型（极大传递子锦标赛）$T_C(V_{S(n)})$ [10]。它们由有序分数序列索引的传递型锦标赛生成，并只需将弧 $v_1 \rightarrow v_n$ 反向为 $v_n \rightarrow v_1$。我们还使用另一类顶点全环型锦标赛（vertex pancyclic 或简称 pancyclic），即进一步将 $v_2 \rightarrow v_{n-1}$ 反向为 $v_{n-1} \rightarrow v_2$。

（b）可约型（Reducible）：通过两种方法以控制方式生成具有不同环路结构的可约型 $T_C(V_{S(n)})$。第一种方式利用已知的有向图结构，如可约型 $T_C(V_{S(n)}), n \geq 2$ 具有强连通分解（详见第 5.3.1 节）。除了生成传递型锦标赛外，还产生包含奇数组分的可约型 $T_C(V_{S(n)})$，其中每个奇数下标组分为单一顶点，偶数下标组分为已知结构的强连通部分，例如顶点全环型（极大传递子锦标赛）以及正则组件 [3]。第二种方法对应于算法 13.3 的相同实现，即利用第 13.2.3 节已用到的基于择优连接（preferential attachment）的网络增长方法生成可约锦标赛。

（c）随机型（Random）：通过给每条边 $\{v_i, v_j\} \in E$ 随机定向生成随机型 $T_C(V_{S(n)})$。
在其底层完全图 $G_C(V_{S(n)}, E)$ 上的每条边的方向均以相等的概率随机选取。尽管对于本研究所考虑的大小，随机的 $T_C(V_{S(n)})$ 很可能是不可约的 [48]，但我们仍然采用基于得分序列的不可约性测试，对随机生成的 $T_C(V_{S(n)})$ 进行筛选，并仅选择那些不可约的情形 [3]。首先，我们在结构和规模受控变化的条件下，研究重启概率 $\alpha$ 对最大范数 $\|\boldsymbol{\pi} - \boldsymbol{\psi}_{\alpha}(\mathbf{s})\|_{\infty}$ 和一范数 $\|\boldsymbol{\pi} - \boldsymbol{\psi}_{\alpha}(\mathbf{s})\|_{1}$ 的影响。图13.4 展示了 $\|\boldsymbol{\pi} - \boldsymbol{\psi}_{\alpha}(\mathbf{s})\|_{\infty}$ 的实验结果。

对于可约的 $T_C(V_{S(n)})$，它们都存在一个吸收态 $v_n$，使得 $\boldsymbol{\pi} = (0, 0, 0, \dots, 1)$。因此，定理13.2意味着 $\boldsymbol{\psi}_{\alpha}(\mathbf{s})$ 是 $\boldsymbol{\pi}$ 的重新分布，本质上体现了原本集中于 $v_n$ 的概率质量的泄漏，即 $\|\boldsymbol{\pi} - \boldsymbol{\psi}_{\alpha}(\mathbf{s})\|_{\infty} = |\pi_{v_n} - \psi_{v_n}|$。由于概率质量泄漏至 $V_{S(n)} \setminus \{v_n\}$ 的 $n-1$ 个顶点，在没有泄漏的不可约情形相比，可约情形下 $\|\boldsymbol{\pi} - \boldsymbol{\psi}_{\alpha}(\mathbf{s})\|_{\infty}$ 会更大。随着3-环的数量增加，导致可约 $T_C(V_{S(n)})$ 中的传递分量减少，$\|\boldsymbol{\pi} - \boldsymbol{\psi}_{\alpha}(\mathbf{s})\|_{\infty}$ 变大并趋近于上界 $(1 - \frac{1}{n})$（图13.4中的粗实线）。这是通过将支配的强连通分量从具有最少3-环的全环型结构转变为在 $B$ 中3-环数量最多的正则结构（$C$中）来实现的 [48][13]。

对不可约情形，平稳分布和PageRank向量的各元素都是非负的（全部大于0），因此 $\|\boldsymbol{\pi} - \boldsymbol{\psi}_{\alpha}(\mathbf{s})\|_{\infty}$ 远小于 $(1 - \frac{1}{n})$。随着3-环数量的增加，即当结构更趋于具有最大数目3-环且 $\boldsymbol{\pi}$ 均匀分布的正则循环型竞赛图时 [14]，$\|\boldsymbol{\pi} - \boldsymbol{\psi}_{\alpha}(\mathbf{s})\|_{\infty}$ 将会减小。这可通过比较D（3-环最少的全环型 $T_C(V_{S(n)})$）、E以及尤其是F（随机 $T_C(V_{S(n)})$ 随顶点数增加逐渐趋近于正则分布，见F中的插图）中的实验结果观察到 [13]。

计算结果表明，$\|\boldsymbol{\pi} - \boldsymbol{\psi}_{\alpha}(\mathbf{s})\|_{\infty}$ 随 $\alpha$ 单调增大（见图13.4中较高虚线，对应更大$\alpha$时可约和不可约情形的结果）。$\|\boldsymbol{\pi} - \boldsymbol{\psi}_{\alpha}(\mathbf{s})\|_{\infty}$ 趋近于上界。
有向图阶数（对数刻度）  
$0 \quad 0.5 \quad 1$  
最大范数 A $10^1 \quad 10^2 \quad 10^3$  
有向图阶数（对数刻度）  
$0 \quad 0.5 \quad 1$  
最大范数 B $10^1 \quad 10^2 \quad 10^3$  
有向图阶数（对数刻度）  
$0 \quad 0.5 \quad 1$  
最大范数 C $10^1 \quad 10^2 \quad 10^3$  
有向图阶数（对数刻度）  
$0 \quad 0.5 \quad 1$  
最大范数 D $10^1 \quad 10^2 \quad 10^3$  
有向图阶数（对数刻度）  
$0 \quad 0.5 \quad 1$  
最大范数 E $10^1 \quad 10^2 \quad 10^3$  
有向图阶数（对数刻度）  
$0 \quad 0.5 \quad 1$  
最大范数 F $10^2 \quad 10^0$  

图 13.4  CMC 算法作用于不同锦标赛类型时的 $\|\boldsymbol{\pi} - \boldsymbol{\psi}_{\alpha}(\mathbf{s})\|_{\infty}$ ：（a）传递性，（b）包含阶数为 $9$ 的泛环（极大传递子锦标赛）分量，（c）包含阶数为 $9$ 的正则分量，（d）泛环（极大传递子锦标赛），（e）泛环，以及（f）阶数为奇数 $n \in [5,999]$ 的随机锦标赛。除 F 图中的插图为对数-对数坐标外，其余所有图均为半对数坐标。虚线表示 $\|\boldsymbol{\pi} - \boldsymbol{\psi}_{\alpha}(\mathbf{s})\|_{\infty}$，其中 $\alpha \in [0.05, 0.95]$，步长为 $0.05$。上界：(粗线) $(1 - \frac{1}{n})$，及 D-E 图中点划线 $(\frac{a_1}{a_2 + H_{m + 1}} - \frac{1}{m + 2})$ [13]，均对 $\|\boldsymbol{\pi} - \mathbf{s}\|_\infty \leq \|\boldsymbol{\pi}^{\mathrm{tran}} - \mathbf{s}\|_\infty = (1 - \frac{1}{n})$ 给出约束。
对于可约的A-C情形，有 $||{\boldsymbol \pi} - {\mathbf{s}}||_{\infty} \leq ||{\boldsymbol \pi}_{\mathrm{tran}} - {\mathbf{s}}||_{\infty} = (1 - \frac{1}{n})$，而对于不可约的D-E情形则有 $||{\boldsymbol \pi} - {\mathbf{s}}||_{\infty}$。尽管该性质在推论13.1中仅针对不可约的CMC得到了证明，我们的计算实验表明，从渗漏和设定均匀分布$\mathbf{s}$的角度来看，对可约CMC有 $||{\boldsymbol \pi} - {\boldsymbol \psi}_{\alpha}({\mathbf{s}})||_{\infty}$ 的单调性[13]。针对在 $T_C(V_{S(n)})$ 上运行的CMC，如果$D$具有由全环（极大传递子锦标赛）组成的特定结构，则可以获得较为精确的界。我们注意到 $||{\boldsymbol \pi}-{\boldsymbol \psi}_{\alpha}({\mathbf{s}})||_{\infty}= |\pi_{v_1}-\psi_{v_1}| \leq ||{\boldsymbol \pi} - {\mathbf{s}}||_{\infty}= |\pi_{v_1} - \frac{1}{n}|$。这与我们对该结构如何影响PageRank计算的理解是一致的。重新分布主要影响顶点$v_1$，而该顶点与$v_n$共同分担了大部分概率质量。如文献[5]所述，$v_1$的PageRank权威性由唯一的链接$v_n \rightarrow v_1$赋予，尽管$v_n$的权威性最高。我们可以利用有向图理论的方法直接计算该界。我们有如下等式：$H_n = \sum_{k = 1}^{n}\frac{1}{k}$（第$n$个调和数），$\sum_{i = 1}^{n}H_i = (n + 1)(H_{n + 1} - 1)$ [19]，以及$\mathbb{E}(\eta_{V_{S(n-2)}}) = \frac{1}{n-2}\sum_{i=1}^{n-2}H_i$。令 $m = n-2, n \geq 3$，将上述等式应用到式(13.22)并注意$\pi_{v_1} = \pi_{v_n}$，可得[13]：

\[
\begin{align*}
\pi_{v_1} &= \frac{1}{2}\left(1 - \frac{\frac{1}{m}\sum_{i=1}^m H_i}{2 + \frac{1}{m}\sum_{i=1}^m H_i}\right) \\
&= \frac{1}{2}\left(1 - \frac{\frac{1}{m}(m + 1)(H_{m+1} - 1)}{2 + \frac{1}{m}(m + 1)(H_{m+1} - 1)}\right)\\
&= \frac{a_1}{a_2 + H_{m+1}},
\end{align*}
\]

其中，$a_1 = \frac{m}{m+1}$，$a_2 = \frac{m-1}{m+1}$。于是，

\[
||{\boldsymbol \pi} - {\boldsymbol \psi}_{\alpha}({\mathbf{s}})||_{\infty} \leq \left(\frac{a_1}{a_2 + H_{m+1}} - \frac{1}{m+2}\right).
\]

对于这些包含大量但有限个顶点的锦标赛图（$n=m+2$），上述不等式右端由$H_{m+1}$的倒数主导（因为$H_{m+1}$呈对数增长），且始终大于零[19]。等号仅在$m=1$时成立，此时右端为零，对应于$n=3$顶点的同构全环锦标赛，该结构是正则的，并且 $||{\boldsymbol \pi} - {\boldsymbol \psi}_{\alpha}({\mathbf{s}})|| = 0$。图13.5给出了同一实验集下 $||{\boldsymbol \pi} - {\boldsymbol \psi}_{\alpha}({\mathbf{s}})||_1$ 的结果。上界 $||{\boldsymbol \pi} - {\boldsymbol \psi}_{\alpha}({\mathbf{s}})||_1 \leq 2(1 - \frac{1}{n})$ 以粗线标记，而 $||{\boldsymbol \pi} - {\boldsymbol \psi}_{\alpha}({\mathbf{s}})||_1 \leq \frac{2}{1 - \alpha}\left(1 - \left(\pi_{v_n} - \frac{1}{n}\right)\right)$ 以点划线标记在图13.5中。对于仅含一个吸收态$v_n$的可约A-C情况，可以数值验证 $||{\boldsymbol \pi} - {\boldsymbol \psi}_{\alpha}({\mathbf{s}})||_1 = 2||{\boldsymbol \pi} - {\boldsymbol \psi}_{\alpha}({\mathbf{s}})||_{\infty}$。任意由最大范数度量的渗漏都会重新分布到所有其他状态中；而1-范数则在最大范数下对$v_n$的渗漏之外，再累计其它渗漏。至关重要的是……
有向图阶数（对数刻度） 0 0.5 1 1.5 2  
$L_1$范数 A $10^1$ $10^2$ $10^3$  
有向图阶数（对数刻度） 0 0.5 1 1.5 2  
$L_1$范数 B $10^1$ $10^2$ $10^3$  
有向图阶数（对数刻度） 0 0.5 1 1.5 2  
$L_1$范数 C $10^1$ $10^2$ $10^3$  
有向图阶数（对数刻度） 0 0.5 1 1.5 2  
$L_1$范数 D $10^1$ $10^2$ $10^3$  
有向图阶数（对数刻度） 0 0.5 1 1.5 2  
$L_1$范数 E $10^1$ $10^2$ $10^3$  
有向图阶数（对数刻度） 0 0.5 1 1.5 2  
$L_1$范数 F $10^2$ $10^{-2}$ $10^0$  
图13.5   
$\|\boldsymbol{\pi} - \boldsymbol{\psi}_{\alpha}(\mathbf{s})\|_1$  
CMC（可组合马尔可夫链）在以下类型锦标赛下的表现：  
(a) 传递性，  
(b) 包含阶数为9的全环（极大传递子锦标赛）组分，  
(c) 包含阶数为9的正则组分，  
(d) 全环（极大传递子锦标赛），  
(e) 全环，  
(f) 阶数$n \in [5,999]$的随机锦标赛。  
除(f)的插图为对数-对数坐标外，其余均为半对数坐标。  
虚线表示$\|\boldsymbol{\pi} - \boldsymbol{\psi}_{\alpha}(\mathbf{s})\|_1$，其中$\alpha \in [0.05,0.95]$，步长为$0.05$。  
上界：$2(1-\frac{1}{n})$（粗线）以及$\frac{2}{1-\alpha}\left(1-\left(\pi_{v_n}-\frac{1}{n}\right)\right)$在$\alpha=0.05$时的取值（(d)-(f)中的点划线）[13]。

图13.6  
$\|\boldsymbol{\pi} - \boldsymbol{\psi}_{\alpha}(\mathbf{s})\|_{\infty}$  
CMC作用于基于网络生成方法、参数$\gamma$分别为(a) 2.0，(b) 0.1，(c) 0.05，(d) 0.01时生成的可约锦标赛的性能。  
蓝色（底部）、黑色（中间）和红色（顶部）的箱线图分别表示最终生成的锦标赛阶数为100、500和1000时的结果[13]。
“泄漏关系”同样扩展到了有向图理论推导的上界之间，因为 $||{\boldsymbol \pi} - {\boldsymbol \psi}_{\alpha}({\mathbf{s}})||_1 \leq 2(1 - \frac{1}{n})$，且 $||{\boldsymbol \pi} - {\boldsymbol \psi}_{\alpha}({\mathbf{s}})||_{\infty} \leq (1 - \frac{1}{n})$ [13]。接下来我们仅考虑通过基于网络增长的方法生成的可约 $T_C(V_{S(n)})$，并仅测量 $||{\boldsymbol \pi} - {\boldsymbol \psi}_{\alpha}({\mathbf{s}})||_{\infty}$。图13.6展示了系统温度 $1/{\gamma}$ 从A变化到D的实验结果。在A中温度较低时，由于在每次方法迭代时引入新节点所生成的边的随机定向，偏向于已有得分更高的节点。这导致生成的 $T_C(V_{S(n)})$ 中具有更加显著的传递结构，这一点可通过更大的Landau指数 $\nu$ 体现出来。最终生成的 $T_C(V_{S(n)})$ 在 $n=100$、$500$、$1000$ 时，$\nu$ 的范围分别为 $[8,49]$、$[8,39]$ 以及 $[7,47]$。在较高温度下，由B到D生成的 $T_C(V_{S(n)})$ 仅包含两个强分量（一个主导节点和一个由 $n-1$ 个节点组成的强支配分量）。更高的 $||{\boldsymbol \pi} - {\boldsymbol \psi}_{\alpha}({\mathbf{s}})||_{\infty}$ 值则是因为较大规模的强支配分量导致更多的“泄漏”，从而将更多的概率质量从吸收的主导节点“拉出”[13]。PageRank向量 ${\boldsymbol \psi}_{\alpha}({\mathbf{s}})$ 给出了顶点 $v_i \in V_{S(n)}$ 的访问概率 $\psi_{v_i}$，可用于衡量顶点的重要性（性能），反映出底层结构的代表性。
10^{-3} 10^{-2} 10^{-1} 10^{0} PageRank 概率（对数坐标）  
A 62 64 66 68 70 0.005 0.0055 0.006 0 20 40 60 80 100 顶点  
10^{-3} 10^{-2} 10^{-1} 10^{0} PageRank 概率（对数坐标）  
B 62 64 66 68 70 0.005 0.0055 0.006 0.0065  

图 13.7 CMC 带重启的 PageRank 向量在 101 阶可约锦标赛各顶点上的访问概率：(a) 具有 9 阶泛环（极大传递子锦标赛）分量的情况，(b) 具有 9 阶正则分量的情况。图中蓝色（下）、黑色（中）、红色（上）叉号分别代表 $\alpha$ 取 $0.05$、$0.15$ 和 $0.25$ 时的访问概率。内嵌图显示当 $\alpha = 0.25$ 时在其中一个强分量内的访问概率。顶点按照它们的得分序列 [13] 以及它们在共进化有向图 $D_C(V_{S(n)})$ 中的对偶关系结构排序。我们采用与得分序列相同的方法对这些顶点的重要性进行 PageRank 排序，以保证一致性。设 $\lambda : \{\psi_{v_i}\}_i^n \rightarrow {\mathbb{Z}}_{1 \leq z \leq n}$ 为 $D_C(V_{S(n)})$ 中顶点 $v_i \in V_{S(n)}$ 的 PageRank，其中访问概率为 $\{\psi_{v_i}\}_i^n$，${\mathbb{Z}}_{1 \leq z \leq n} = \{1,2,3,\ldots,n\}$。顶点 PageRank 的排序为升序，从最不重要的 $\lambda(v_1) = 1$ 到最重要的 $\lambda(v_n) = n$。若出现并列，则采用平均排名法（分数排名）[39]。在我们的实验中，生成了奇数阶的锦标赛与强分量，使得可以用中位数（整数）为并列结果赋值。结果总结如图 13.7 所示。

图 13.7 显示，对于传递型 $T_C(V_{S(n)})$，即使在较高 $\alpha$ 设置下，顶点的 PageRanks $(\lambda)_i^n$ 仍与其得分序列 $\big(d_T^-(v_i)\big)_i^n$ 保持一致。这并不令人意外，因为吸收节点的概率质量被均匀再分配。更有意义的是涉及不可约 $T_C(V_{S(n)})$ 并存在明显传递结构的情况。对于泛环（极大传递子锦标赛）$T_C(V_{S(n)})$，访问概率排序为 $\psi_{v_2} < \psi_{v_3} < \psi_{v_4} < \cdots < \psi_{v_{n-1}} < \psi_{v_1} < \psi_{v_n}$，在较低 $\alpha$ 设置下成立。与平稳概率向量情形不同，后者中最高的两个排名会并列（$\pi_{v_1} = \pi_{v_n}$），PageRank 能够区分 $v_1$ 和 $v_n$ [13]。若此结构嵌入到可约 $T_C(V_{S(n)})$（如实际复杂结构问题中常见）中，对于规模较小的 $T_C(V_{S(n)})$，PageRank 可用于揭示该强分量结构，因为不存在规模性问题。例如，当强分量为泛环（极大传递子锦标赛）时（见图 13.7A 内嵌图），节点 $v_{62}$ 和 $v_{63}$ 在全局上具有相同的出度（$d_T^-(v_{62}) = d_T^-(v_{63})$），且在该强分量内局部均仅有一个向外发出的链接，然而 $\psi_{v_{62}} > \psi_{v_{63}}$。这些信息都指示它们属于主导了顶点 $\{v_1, \ldots, v_{61}\}$ 的强分量。如果强分量是正则的，则分量内所有顶点的 PageRank 权重相同，因而 PageRanks 完全并列（见图 13.7b 的内嵌图）[13]。

已知选择合适的 $\alpha$ 设置可以提升幂法计算 PageRank 的收敛速度（见[41]定理 5.1），但代价是会放大 PageRank 权重波动，从而改变实际排名顺序 [30]。我们通过计算基准值 $\alpha = 0.05$ 与较大 $\alpha$ 时所得 PageRanks 间的排名差异，考察 $\alpha$ 对 PageRanks 的影响。排名差异采用秩聚合方式计算。我们使用两种不同的距离度量。设 $T_C(V_{S(n)})$ 为奇数阶 $n$，第一种度量是 Spearman 足规距离 [22]：
其中，斯皮尔曼距离定义为：

$$d_{\mathrm{spear}}(\lambda_1,\lambda_2) = \frac{4}{n^2 - 1} \sum_{i = 1}^n |\lambda_1(v_i) - \lambda_2(v_i)|$$

这里的归一化常数由绝对差之和的最大值给出：

$$\max d(\lambda_1,\lambda_2) = \frac{n^2 - 1}{4}$$

该归一化方法与斯皮尔曼等级相关系数 [39] 使用相同，旨在计算完全相反排序下获得的距离 [9]。最大距离对应于传递赛和规则锦标赛之间的情况。第二个衡量指标是平均等级差：

$$d_{\mathrm{avg}}(\lambda_1,\lambda_2) = \frac{1}{n} \sum_{i = 1}^n |\lambda_1(v_i) - \lambda_2(v_i)|$$

$d_{\mathrm{avg}}(\lambda_1,\lambda_2)$ 提供了与 $d_{\mathrm{spear}}(\lambda_1,\lambda_2)$ 的有用对比，适用于某些大型锦标赛实验中，当 $n^2$ 因子可能掩盖了小的排名差异影响的情况 [13]。如图 13.8 所总结的结果表明，较高的 $\alpha$ 设置会增强对 PageRank 的扰动，这可以通过更大的 $d(\lambda_1, \lambda_2)$ 值衡量（例如，图 13.8 中的粗线）。然而，我们观察到，随着所生成的 $T_C(V_{S(n)})$ 规模的增大，$d(\lambda_1, \lambda_2)$ 的值反而降低。这些结果说明 PageRank 的变化实际上局限于某些特定的顶点。例如，对于泛环（极大可传递子锦标赛）$T_C(V_{S(n)})$ 的 A-B 情况，通常是 $v_1$。除了 $T_C(V_{S(n)})$ 的特定结构外，我们的结果进一步表明，显著的可传递结构（A-B）及那些强制可约性的结构（C-D）的存在，会使 PageRank 的扰动局限于少数顶点，从而随着 $T_C(V_{S(n)})$ 规模增大，导致 $d(\lambda_1, \lambda_2)$ 的值下降。此外，我们普遍观察到，对于同阶 $n$，可约 $T_C(V_{S(n)})$ 的 $d(\lambda_1, \lambda_2)$ 值通常要低于不可约 $T_C(V_{S(n)})$ [13]。这并不出乎意料，因为重启操作会将概率质量从吸收节点均匀重新分配到瞬态节点，这也是 PageRank 随机游走中普遍的行为。在实际应用中，可以预期分析需要在可解共进化问题上进行，这些问题可由可约的共进化有向图表示。算法的 PageRank 计算性能由此受到显著影响。
图 13.8 不可约合作演化锦标赛的 PageRank 差异，包含 (a)-(b) 泛环（极大传递子锦标赛）结构和可约的 (c)-(d) 泛环（极大传递子锦标赛）分量，阶数为 9。图 (a) 和 (c) 展示了 $d_{\mathrm{spear}}(\lambda_1, \lambda_2)$，图 (b) 和 (d) 展示了 $d_{\mathrm{avg}}(\lambda_1, \lambda_2)$。虚线表示 $\alpha = 0.05$ 时获得的基线 PageRank 与设置更高 $\alpha \in [0.10, 0.90]$ 时获得的 PageRank 之间的差异。粗实线对应 $\alpha = 0.95$ 的情形[13]，此时算法的速度（例如迭代法的收敛速度）变得非常关键。PageRank 直接应用于竞争性合作演化中的一个重要问题在于计算成本。与 PageRank 通常应用于稀疏连接网络[27]并且已有提升计算速度方法[41]的情形不同，合作演化的有向图在定义上是完全连接的。尽管如此，我们的计算研究表明，PageRank 提供了有意义的定量分析，有助于揭示合作演化有向图中的问题结构，因为 PageRank 随机游动者在机制上体现了这些定性信息。该研究旨在奠定基础，使相关的定量分析成为可能，尽管对于大规模、现实世界的合作演化系统仍然存在真正的挑战。相比在实际中进行全局分析（计算量极大），可以将分析重点放在对应于搜索空间中局部邻域的特定子有向图的局部分析。这可以通过对瞬移向量 $\mathbf{s}$ 的适当设置来实现。此外，也有空间发展相关方法，实现在实际中对完整有向图的强连通子分量进行相对较短的 PageRank 计算，并将这些运行结果进行聚合，从而在这些子分量连通性的意义下构建相关且更宏观的整体图景。
13.4 讨论与结论

本章旨在解决竞争环境中的一个基本问题，即问题结构如何影响通过协同进化迭代方式发现的解的性能。在一些温和的假设下，单个种群协同进化过程可以被表示为在有向图上的马尔可夫链模型，这些有向图刻画了问题中相互竞争解之间的偏好关系。这种将抽象协同进化形式化为有向图上的随机游走的方法，使我们能够对竞争性协同进化系统进行深入研究。除了确保协同进化问题结构的二分性能够传递到协同进化搜索过程的总体定性结果外，即（不可）约 CMC 操作于（不可）约的协同进化有向图之外，这些对竞争性协同进化的定性刻画为构建用于协同进化搜索过程严格定量分析的工具提供了有益的信息。

特别地，我们提出了将“期望到达吸收类的时间”（expected hitting times to the absorbing class）的概念，作为衡量单种群协同进化搜索过程发现主导解子集速度的方法。我们还获得了一些重要见解，使我们能够精确地说明特定领域知识如何影响协同进化搜索的速度。例如，标准随机游走是一种特定的 CMC，它通过顶点的出邻域纳入了领域知识，而 $(1+1)$ CEA 并未利用任何领域知识，表现为盲目搜索。在这两种 CMC 的最坏情形比较中（即从最弱顶点出发），标准随机游走具有亚线性时间复杂度 $O(\ln n)$，但如果搜索起点更接近主导顶点，则其速度可显著提升。

我们为抽象协同进化提出的另一项定量度量是不可约 CMC 的平稳分布，其对应协同进化搜索访问顶点的概率。由于这种不变测度需要先验的定性知识（如协同进化有向图的不可约性等），我们建立了协同进化过程马尔可夫链模型与作用于这些有向图上的 PageRank 之间的直接且正式的联系。我们保证了对任何单种群协同进化系统都存在 PageRank。我们证明，PageRank 向量实际上是通过引入重启机制后 CMC 在协同进化有向图上的平稳向量的再分配。由具体重启概率 $\alpha$ 设定而导致的 PageRank 权威分数的变化可以被精确量化。更为关键的是，引入重启机制后，即使原本是可约的 CMC 也能够被有效转化为不可约 CMC。因此，PageRank 权威分数构成了对原始平稳分布度量的细化。这使我们能够为所有类型的协同进化有向图（包括可约和不可约）建立一种有原则的定量刻画方法。

PageRank 随机游走者通过权威度对协同进化有向图进行表征，这些权威度衡量各个解（顶点）的性能，并通过（递增的）权威度顺序对应其 PageRank。当 PageRank 被用来对协同进化有向图进行全局分析（即设置均匀的随机跳转向量）时，无论是在理论上还是计算上都已证明，PageRank 能保持可约有向图中顶点基于分值序列的原序。最后，经适当归一化后的 PageRank 权威度在协同进化搜索带有重启的情形下，自然地解释为搜索过程访问有向图中顶点的概率。
虽然我们研究的主要动机是解决在对协同进化系统进行定量分析时，面对问题结构几乎毫无先验定性知识而产生的基本挑战，但我们的技术贡献具有普适性，可应用于其他采用在半完全有向图（一个边可以是双向的）上随机游走模型的问题领域。在对抗博弈的语境下，这意味着交互关系包含了平局的情形。我们所考虑的情形是受限但广泛的一类单群体协同进化系统，其中涵盖了涉及自对弈的其他学习算法。然而，PageRank方法也已被应用于演化系统，其中底层问题的结构被建模为特定的可约有向图 [33]。此外，我们的技术方法还扩展了有关一般化为有向图的随机游走 [17, 43] 以及相应电网络模型 [2, 51] 的丰富研究领域。相比之下，本研究关注于在不可逆CMC（连续马尔可夫链）上运行的、既包括不可约也包括可约的协同进化半完全有向图，其中完全定向图（竞赛图）是其一个子集。我们的理论研究经常利用特定有向图的性质（如强连通性）来证明不可约有向图上的理论结果。重要的是，我们的方法受到有向图理论、尤其是竞赛图理论的启发，从而能够建立一种更为自然的形式体系，在交互可以建模为两人对策博弈的情形下，研究竞争性协同进化。我们以对于竞争性协同进化研究现状、特别是理论研究状况的一些说明作为本章的结语。

我们首先指出，协同进化系统极难被系统研究。其结果之一是，这些系统的不同方面往往（但非总是）被孤立地在多种常见理论框架下分别研究。即使在同一框架下，由于协同进化系统的设计、问题结构和动态性的相互影响，人们也会采用不同的方法进行研究。例如，研究 [24, 62] 基于进化博弈论的协同进化动力学，使用连续动力系统理论框架，从人口层面分析了纯策略混合体在无限体量主体群体中的表现。此时，协同进化搜索关注的焦点往往是那些策略处于某种均衡状态（如纳什均衡）的特定混合体。然而，不同于演化动力学 [63]，协同进化动力学可以表现为从单一吸引子（对应纳什均衡）到周期性甚至复杂混沌动力学的多种状态 [24]。在后者情形下，必须开展更为深入的研究，并引入动力系统理论中的一些复杂工具（如结构稳定性），以判定协同进化动力学系统是否具有影子性（Shadowing）[62]。自然和计算生成的轨迹，往往由于有限种群效应和选择及变异中的随机性，表现为伪轨道。确立协同进化中的影子性至关重要，因为如果伪轨道能够在时空上被真实轨道紧密跟踪（即被“影子化”），那么就可以认为这些伪轨道构成对协同进化过程相关行为的有效观测。尽管我们无法证明这些协同进化动力系统一般都具备影子性，但通过进一步的理论构造（例如在约化动力学设定下，混沌轨道存在于某一不变子集内），是有可能确立某种形式的影子性的。
另一个例子是在帕累托协同进化框架下的各种研究 [21]。在这些工作中，常见的做法是使用一组测试用例评估每个协同进化个体解的性能，并像多目标优化问题（MOP）那样，从帕累托支配的角度定义解的质量。然而，与进化多目标优化问题场景不同，在帕累托协同进化中，每一个测试用例都作为需要被优化的目标。如果进一步考虑单种群的帕累托协同进化问题，候选解本身也可以作为测试用例，这使问题更加复杂。尽管如此，已有若干基础性研究为帕累托协同进化框架的发展做出了贡献。文献 [64] 提出了一种问题生成方法，可以在测试目标中以受控方式引入特定属性（如非传递性）。在 [21] 中，发展了能够检测和利用个体层次或非支配解子集之间以帕累托覆盖顺序表现出来的梯度的协同进化算法（CEA）。文献 [8] 基于序理理论提出了理论基础，研究解集合的偏序集上几何结构的涌现，把此理论用于分析 [64] 中所引入问题家族，从而为更具洞察力的帕累托协同进化搜索（比较）提供理论支持。这些框架之间是否存在明显的联系？例如，帕累托协同进化与本章所用的协同进化问题的有向图表示及协同进化过程中的随机游走框架之间的关系如何？一个紧密的联系在于，[8] 中用于揭示成对比较空间结构的形式化方法，以及我们在 [14] 首次提出的对协同进化有向图中底层结构的定性刻画，这些对于有效帕累托协同进化搜索的设计与分析都是至关重要的。对于任何对抗性协同进化问题，关于问题的定性知识对于保证问题可解（即存在主导解子集）是至关重要的。其底层可约协同进化有向图允许把解集合分解为具有无环排序的强连通分量。当在帕累托协同进化中有意使用解集合本身进行评估时，解集合（搜索空间）上所诱导的帕累托覆盖顺序就对应于该可约有向图的强分解结构。这表明可在这些帕累托层（强分量）之间，进而在最终的唯一主导解子集之间，挖掘并利用结构性信息以支持协同进化搜索。
为了进一步说明这两种框架的互补特性，可以构建一个关于伪嵌入共进化解的程序的归纳性论证：

1. 取一个 $n \geq 4$ 个顶点的传递性锦标赛（transitive tournament），顶点按其得分序列 $\{v_i\}_1^n$ 标记，排列顺序为 $s_1 < s_2 < s_3 < \ldots < s_n$。

2. 对 $v_1v_n$ 进行一次有向弧反转，变为 $v_nv_1$，从而获得一个泛圈（pancyclic，极大传递子锦标赛）锦标赛。

3. 该泛圈锦标赛被伪嵌入到 $\mathbb{N}^2$ 上，其线性实现为 $L = \{()_X, ()_Y\} = \{([v_1], [v_2], [v_3], \ldots, [v_n]), ([v_n], [v_2], [v_3], \ldots, [v_1])\}$。我们用记号 $L = \{()_X, ()_Y\}$ 表示该伪嵌入可以在 $\mathbb{N}^2$ 空间中可视化，即 $x$ 轴由 $()_X$ 张成，$y$ 轴由 $()_Y$ 张成。

该泛圈结构的一个特定性质是，包含 $n-2$ 个顶点的子锦标赛 $\{v_2, v_3, v_4, \ldots, v_{n-1}\}$ 始终具有传递性，这使得可以设立归纳论证，并利用上述程序对 $n \geq 4$ 个顶点的锦标赛进行伪嵌入。归纳步骤可以被可视化为，将由两端点指定的嵌入点集从 $\{(v_1, v_n), (v_n, v_1)\}$ 移动到 $\{(v_1, v_{n+1}), (v_{n+1}, v_1)\}$——在后者中，添加了一个额外点，对应于形成传递子锦标赛的顶点 $\{v_2, v_3, v_4, \ldots, v_{n}\}$。

未来研究中，将其他强锦标赛的伪嵌入推广到超越 $\mathbb{N}^2$ 空间的情形，将是一个有趣的研究方向。除此之外，还有其他框架被用于对共进化系统的形式化研究。在第 12 章中，我们将竞争性共进化建模为共进化学习系统，并大量借鉴了机器学习方法论中的泛化概念与分析方法 [15, 52, 57]。竞争性共进化还可以在共优化（cooptimization）的背景下建模，通过在全域问题实例上对 CEA 的整体性能进行最坏情形优化分析进行比较 [52, 55]。文献 [38] 则从理论层面对共进化相互作用的本质进行了基础且深入的研究。其形式化了一个基于偏好函数 $f: \mathcal{A} \rightarrow 2^S$ 的共进化设定。其核心思想是，集合 $\mathcal{A}$ 表示共进化用以判定 $S$ 中哪些解被优选（并因此成为 $S$ 的幂集 $2^S$ 元素）的信息状态。该设定要求集合 $\mathcal{A}$ 具备代数结构，例如偏序集 $(\mathcal{A}, \subseteq)$，其中 $\mathcal{A} = 2^{S \times T}$，即 $S$ 和 $T$ 分别为解集与测试案例集时的幂集。

文献 [38] 的一个重要理论结果是，只要偏好函数 $f$ 涉及的是有限集合，其必然收敛。非正式地说，这一性质表明：对于 $S$ 中某些解在信息状态 $\mathcal{A}$ 累积到极限时已建立的偏好，其实在更早的累积阶段就已成立。感兴趣的读者可参考 [38] 获取更多技术细节。

在本章结尾，我们希望强调，面对广泛而深刻的研究框架，对于共进化系统的完全理解仍面临巨大挑战。这也消除了“共进化相关研究缺乏严谨性”的看法，事实上，情况恰恰相反。
参考文献

1. Axelrod, R.: 重复囚徒困境中策略的演化. 载于 L.D. Davis (编), 《遗传算法与模拟退火》，第3章, 第32–41页. Morgan Kaufmann, 纽约 (1987)
2. Balázs, M., Folly, A.: 非可逆马尔科夫链的电网络. 《美国数学月刊》 123(7), 657–682 (2016)
3. Bang-Jensen, J., Gutin, G.Z.: 有向图理论、算法与应用. Springer, 伦敦 (2009)
4. Berman, A., Plemmons, R.J.: 数学科学中的非负矩阵. SIAM, 费城 (1994)
5. Bianchini, M., Gori, M., Scarselli, F.: PageRank的内部机制. 《ACM互联网技术汇刊》 5(1), 92–128 (2005)
6. Brin, S., Page, L.: 大规模超文本Web搜索引擎的结构. 见第七届国际万维网大会（WWW7）论文集, 第107–117页. 澳大利亚布里斯班 (1998)
7. Brualdi, R.A.: 有向图的谱. 《线性代数与其应用》 432, 2181–2213 (2010)
8. Bucci, A.: 协同进化算法中涌现的几何组织与信息维数. 博士论文, 布兰迪斯大学, 马萨诸塞州 (2007)
9. Burns, T., Meeker, L.D.: 评价、决策与社会交互的数学模型. 见 J. Cochrane, M. Zeleny (编), 《多准则决策制定》，第141–163页. 南卡罗来纳大学出版社, 哥伦比亚 (1973)
10. Burzio, M., Pelant, J.: 关于最少3-圈的完全强连通有向图. 《离散数学》 155, 27–30 (1996)
11. Chellapilla, K., Fogel, D.B.: 演化、神经网络、博弈与智能. 《IEEE会刊》 87(9), 1471–1496 (1999)
12. Chong, S.Y., Tan, M.K., White, J.D.: 观察神经网络学习黑白棋博弈的演化过程. 《IEEE进化计算汇刊》 9(3), 240–251 (2005)
13. Chong, S.Y., Tiˇno, P., He, J.: 协同进化系统与PageRank. 《人工智能》 277, 103164 (2019)
14. Chong, S.Y., Tiˇno, P., He, J., Yao, X.: 协同进化系统分析的新框架——有向图表示与随机游走. 《进化计算》 27(2), 195–228 (2019)
15. Chong, S.Y., Tiˇno, P., Yao, X.: 协同进化学习中泛化性能的度量. 《IEEE进化计算汇刊》 12(4), 479–505 (2008)
16. Chung, F.: 有向图的拉普拉斯算子与Cheeger不等式. 《组合学年鉴》 9(1), 1–19 (2005)
17. Chung, F., Zhao, W.: PageRank与图上的随机游走. 见《组合学与计算机科学盛典》, Bolyai数学研究丛书，卷20，第43–62页. Springer, 柏林 (2010)
18. Chung, K.L.: 具有平稳转移概率的马尔可夫链. Springer, 柏林 (1960)
19. Conway, J.H., Guy, R.K.: 数字之书. Springer, 纽约 (1996)
20. Cormen, T.H., Leiserson, C.E., Rivest, R.L., Stein, C.: 算法导论（第2版）. MIT出版社, 剑桥 (2001)
21. de Jong, E.D., Pollack, J.B.: 来自协同进化的理想评价. 《进化计算》 12(2), 159–192 (2004)
22. Dwork, C., Kuman, R., Naor, M., Sivakumar, D.: Web的排序聚合方法. 见第10届国际万维网大会(WWW’01)论文集，第613–622页. 香港 (2001)
23. Falgueras-Cano, J., Falgueras-Cano, J., Moya, A.: 基于进化元胞自动机的数字生物协同进化研究. 《生物学》 10(11), 1147 (2021)
24. Ficici, S.G., Melnik, O., Pollack, J.B.: 协同进化选择方法的博弈论与动力系统分析. 《IEEE进化计算汇刊》 9(6), 580–602 (2005)
25. Fogel, D.B., Hays, T.J., Hahn, S.L., Quon, J.: 自学习进化国际象棋程序. 《IEEE会刊》 92(12), 1947–1954 (2004)
26. García-Pedrajas, N., Hervás-Martínez, C., Ortiz-Boyer, D.: 基于协作协同进化的人工神经网络集成在模式分类中的应用. 《IEEE进化计算汇刊》 9(3), 271–302 (2005)
27. Gleich, D.F.: PageRank在Web之外的应用. 《SIAM评论》 57(3), 321–363 (2015)
28. Golub, G.H., Van Loan, C.F.: 矩阵计算（第4版）. 约翰·霍普金斯大学出版社, 马里兰 (2013)
29. Grinstead, C.M., Snell, J.L.: 概率论导论. 美国数学学会, 普罗维登斯 (1997)
30. Haveliwala, T.H.: 主题敏感的PageRank：一种上下文相关的网页搜索排名算法. 《IEEE知识与数据工程汇刊》 15(4), 784–796 (2003)
31. He, J., Yao, X.: 面向分析的进化算法计算时间分析框架. 《人工智能》 145, 59–97 (2003)
32. Hemasinha, R.: 一种生成锦标赛得分序列的算法. 《数学与计算机建模》 37, 377–382 (2003)
33. Herrmann, S., Ochoa, G., Rothlauf, F.: 基于PageRank中心性进行性能预测：局部最优网络模型的影响. 《启发式方法杂志》24(3), 243–264 (2018)
34. Hillis, W.D.: 共进化寄生体作为优化过程提升模拟进化. 《物理学D》 42, 228–234 (1990)
35. Horn, R.A., Johnson, C.R.: 矩阵分析（第2版）. 剑桥大学出版社, 剑桥 (2012)
36. Iosifescu, M.: 有限马尔科夫过程及其应用. Wiley, 纽约 (1980)
37. Jeh, G., Widom, J.: 个性化Web搜索的扩展. 见第12届万维网大会(WWW’03)论文集，第271–279页. 匈牙利布达佩斯 (2003)
38. Jung, A., Rowe, J.E.: 偏好函数的收敛性. 《理论计算机科学》488, 66–77 (2013)
39. Kendall, M.G.: 排名相关方法（第4版）. Charles Griffin, 伦敦 (1975)
40. Kleinberg, J.M.: 超链接环境中的权威信息源. 《ACM学报》 46(5), 604–632 (1999)
41. Langville, A.N., Meyer, C.D.: PageRank的内部机制. 《ACM互联网技术汇刊》 1(3), 335–380 (2004)
42. Levin, D.A., Peres, Y., Wilmer, E.L.: 马尔科夫链与混合时间. 美国数学学会, 普罗维登斯 (2009)
43. Li, Y., Zhang, Z.: 有向图上的随机游走、广义有向图拉普拉斯算子与不对称度. 见《Web-Graph的算法与模型》(WAW 2010)，《计算机科学讲义丛书》，卷6516，第74–85页. Springer, 柏林 (2010)
44. Lu, X., Menzel, S., Tang, K., Yao, X.: 基于协同协同进化的设计优化：并行工程视角. 《IEEE进化计算汇刊》 22(2), 173–188 (2018)
45. Ma, X., Li, X., Zhang, Q., Tang, K., Liang, Z., Xie, W., Zhu, Z.: 协同协同进化算法综述. 《IEEE进化计算汇刊》 23(3), 421–441 (2019)
46. Meyer, C.D.: 矩阵分析与应用线性代数. SIAM, 费城 (2000)
47. Meyn, S., Tweedie, R.L.: 马尔科夫链与随机稳定性（第2版）. 剑桥大学出版社, 剑桥 (2009)
48. Moon, J.W.: 关于锦标赛的若干主题. Holt, Rinehart and Winston (1968)
