第5章 共演搜索中的问题与分析  
**5.1 共演的一般挑战**  

在第3章中，我们已经介绍并讨论了竞争性和协作性共演的两种框架，它们分别以算法$3.1$和算法$3.2$的形式呈现。在我们讨论竞争性和协作性共演算法（简称CEAs）时，尤其在关于其基本过程组件（如变异算子和选择算子）的问题上，经常提到它们与进化算法（EAs）的对应过程组件具有相似的设计考虑，例如算法$1.1$中所展示的内容。这意味着，在设计具有高效组件的CEAs中所遇到的大多数挑战，与设计性能优异的EAs的情况是相同的。特别是，在为解决问题选择合适的设计选项和具体规范时，人们需要考虑解决方案表示形式、变异算子和选择算子中的那些能够促进（共）演化搜索有效且高效运行的组件。

第1章已经详细介绍了这些EAs基本组件相关的总体细节和所面临的挑战，这些内容同时适用于CEAs。在第$1.2.2$节中，我们引入了EAs中使用的各种解决方案表示形式，并简要讨论了由于它们特定的构造所需的设计考虑。第3章中对应用于可加分离优化问题的协作性CEAs的讨论，多数情况下解决方案表示形式的实现比较直接，特别是在设定为连续优化问题时，考虑已知的基准问题[41, 63, 66]。然而，第4章介绍了一些更具挑战性的解决问题环境，其中包含在CEAs中更复杂的解决方案表示形式的应用。例如，研究[60]中使用协作性CEA与基因编程（GP）树表示形式来设计用于云计算资源分配规则的方法。

在竞争性CEAs的应用中，特别是用于涉及游戏的问题，许多研究探讨了使用复杂的人工神经网络（ANNs）来提供游戏策略的可能性[13, 16, 37, 40]。
任何解决方案表示形式的使用都与特定的变异算子紧密相关。这是因为变异算子的主要职责是根据当前种群中的父代候选解生成新的子代候选解。事实上，我们在1.2.3中已经介绍和讨论了解决方案表示形式与变异算子的组合如何共同指定搜索空间的邻域结构。因此，在演化算法（EAs）中，解决方案表示形式和相关变异算子的设计考虑之间存在着紧密联系。同样，可以发现，在协同进化算法（CEAs）中，也存在类似的解决方案表示形式和变异算子的研究与开发。

例如，当使用具有多个层的复杂人工神经网络（ANNs）时，就需要采用复杂的变异算子，并引入参数控制机制（如自适应调整self-adaptation [13, 16]），以便能够有效地协同进化搜索两人棋盘游戏的策略。然而，如果能够发现并利用领域知识所提供的有用的不变量和其他重要结构，从而实现更有效、简洁的解决方案表示，就可以使用较为简单的变异算子，同时依然能够进行有效的协同进化搜索[46, 48]。

鉴于CEAs和EAs中的选择算子的作用机制都基于与候选解关联的适应度值，以决定候选解在下一代中被复制的频率，因此选择算子的设计对（协同）进化过程的搜索性能会产生影响。如同之前，对于EAs设计的选择算子（参见1.2.4），也可以在CEAs中实现。然而，CEAs与EAs框架之间的主要概念区别在于生成候选解的适应度评估。因此，本章将着重于协同进化过程中固有的主要动态行为，特别是由于候选解相对的动态适应度评估而导致的这些行为。本章的目标是揭示协同进化过程中为何以及如何出现这些负面影响搜索性能的病理行为。

本章剩余内容分为两个主要部分。第5.2节首先会介绍并非正式地描述这些协同进化的病理现象。为了更清楚地展示这些病理现象如何影响协同进化搜索完整解的过程，我们将在讨论中使用竞争协同进化的设置。第5.3节将介绍近年来关于构建新框架的研究，这些框架能够捕获协同进化问题与过程中的特定结构。这进一步推动了用于正式分析协同进化的工具的发展，并因此更精确地描述了一些协同进化病理现象。最后，本章将总结如何解决这些病理现象，并减轻其对协同进化搜索过程的负面影响。

### 5.2 
协同进化搜索中的病理现象
在本节中，我们将介绍协演化病理学（coevolutionary pathologies），其统称为协演化过程中观察到的一系列特征行为或动态，这些行为或动态会影响协演化系统的搜索性能 [6, 8, 11, 14, 22, 30, 55]。前两个小节将介绍并讨论那些在文献中普遍已经被很好理解的协演化病理，并针对这些问题的研究开发了理论和实证方法。尽管需要研究简化且更加抽象的协演化系统，采用理论方法的一个关键优势在于它能够更精确地描述这些病理。第5.2.1节描述了第一组病理，包括过度专业化（overspecialization）、聚焦（focusing）和非传递性（intransitivity），这些问题可以通过研究协演化所应用解决的目标问题中的结构来理解。第5.2.2节描述了第二组病理，包括循环（cycling）和遗忘（forgetting），这些问题可以通过研究潜在协演化过程的结构来理解。协演化过程中循环动态与协演化问题解空间关系的非传递性之间存在一定的联系。第5.2.3节将介绍并讨论脱离（disengagement）和平庸稳定状态（mediocre stable states）的病理，这些问题的相关研究文献主要通过实证研究来探索。

### 5.2.1 过度专业化、聚焦和非传递性

在协演化系统中，过度专业化指的是种群中的候选解集中于问题的某些特定属性，而不是处理所有相关属性，从而解决整个问题。这种特定的协演化病理已经被大量研究，并且通常被很好地理解。有两个视角可以用来更好地理解协演化中的过度专业化问题。在竞争性协演化设置下，一个视角是从候选解之间的成对交互的水平来考察问题，例如将其作为两人竞争游戏。在其他方面，协演化试图解决的问题——协演化问题——实际上可以被看作是两人竞争游戏。作为一种解概念，可以考虑通过成对方式，某种策略对策略集合中的所有其他策略的支配性（dominance）。如果一个两人竞争游戏存在一种策略，它在策略集合中唯一地击败其他所有策略，则该游戏被认为具有主导策略（dominant strategy） [51]。然而，拥有非常大的策略集合的游戏可能在策略集合之上具有复杂的成对或二进制关系。特别地，即使存在一个唯一的主导策略，也不会否定在剩余被支配策略的子集中，必须存在一个策略能够支配该子集中的所有其他策略。

实际上，现在已经明确，这些游戏可以接受形成一个有序序列的二进制关系结构（策略子集）。特别地，越在高序位置的组件（策略子集）能够支配在较低序位置的组件，从而形成一个传递的支配组件链。然而，每个组件内并不存在单一的主导策略 [19, 50]。在这种情况下，任何策略仅专注于击败某些策略，而非击败组件中的所有策略 [24]。这种情况的出现源自于该组件中策略二进制关系的非传递性。给定一个策略集合$\boldsymbol{S}$，一个策略$s_i \in \boldsymbol{S}$支配（击败）另一个策略$s_j \in \boldsymbol{S}$可表示为$s_i \leftarrow s_j$。当对于所有$s_i, s_j, s_k \in \boldsymbol{S}$，只要$s_i \leftarrow s_j$且$s_j \leftarrow s_k$，便有$s_i \leftarrow s_k$时，则策略集合$\boldsymbol{S}$上的二进制关系$A$被认为是传递的。

非传递性发生在某些情况下存在$s_i, s_j, s_k \in \boldsymbol{S}$，使得$s_i \leftarrow s_j$, $s_j \leftarrow s_k$且$s_i \rightarrow s_k$。换而言之，策略集合$\{s_i, s_j, s_k\}$组成了一个循环，而不是一个传递的三元组，就像“石头、剪刀、布”游戏的简单例子一样。我们将在后续的第5.3节中更详细地重新审视这一问题，并提出一种基于有向图理论（Directed Graph Theory）的理论工具来分析协演化中的此类问题结构。
另一种视角是在交互或游戏规则中考虑问题。特别是，研究文献[65]提出了一种数字游戏，这种游戏为研究协演化系统提供了两个有用的特性。首先，这些游戏可以通过对象的成对比较（关系）来说明非传递性的本质，这种关系只能通过比较来决定哪些对象被认为是优越的。其次，通过以某种特定方式对集合中的对象进行组织，可以得到几何结构[8]，这种几何结构可以被利用来专门突出协演化系统中的聚焦问题[9]。

通常来说，数字游戏涉及自然数集合$N = {\mathbb{N}}$的$n$元组。一种特定的数字游戏通过确定两个$n$元组之间比较结果的方法以及$n \geq 2$的选择来具体化。例如，当考虑$n = 1$并禁止自比较（即不允许比较相同的数字时），那么在集合$N = {\mathbb{N}}$上的二元关系$<<$是传递的，从而该数字游戏表示的问题变得很简单。大多数研究[9, 65]考虑了涉及2元组（例如，一对数字）的双人数字游戏，因为这些游戏已经能够揭示协演化的病理现象，例如代理的过度专业化，即仅专注于解决（击败）对手游戏策略中的某些方面。

在对称双人数字游戏中，两位玩家可以选择从集合$\boldsymbol{S} \subset {\mathbb{N}}^2$中实施策略，这由2元组$(x, y) \in {\mathbb{N}}^2$指明。设第一位玩家实施策略$(x_1, y_1)$，而第二位玩家实施策略$(x_2, y_2)$。对于被称为非传递游戏（intransitive game）的双人数字游戏，其对于第一位玩家的收益函数定义如下：

\[
p_{\mathrm{IG}}\big((x_1,y_1),(x_2,y_2)\big) = \left\{ 
\begin{array}{ll} 
1 & \quad \text{如果 } \delta_x > \delta_y \land y_1 > y_2 \\ 
1 & \quad \text{如果 } \delta_y > \delta_x \land x_1 > x_2 \\ 
0 & \quad \text{否则}, 
\end{array} 
\right.
\]

其中，$\delta_x = |x_1 - x_2|$, $\delta_y = |y_1 - y_2|$[9]。

为了简化，可将游戏限制为仅赢输结果，即拥有较大收益值的玩家获胜，并确保当两位玩家具有相同收益值时不会出现平局。可以通过从集合$\boldsymbol{S}$中筛选一对策略来达到这一目的，使得：(i) $x_1 \neq x_2$ 且 $y_1 \neq y_2$（即无自比较）；(ii) $\delta_x \neq \delta_y$。

按照文献[65]中提供的示例，三个策略$s_1 = (1, 6), s_2 = (4, 5), s_3 = (2, 4)$形成了一个循环。该循环基于以下直接计算结果得到：

1. $s_1 \leftarrow s_2$，其中$p_{\mathrm{IG}}(s_1, s_2) = 1$，因为$\delta_x > \delta_y \land y_1 > y_2$，以及$p_{\mathrm{IG}}(s_2, s_1) = 0$，因为$\delta_x > \delta_y \land y_1 \ngtr y_2$。
2. $s_2 \leftarrow s_3$，其中$p_{\mathrm{IG}}(s_2, s_3) = 1$，因为$\delta_x > \delta_y \land y_1 > y_2$，以及$p_{\mathrm{IG}}(s_3, s_2) = 0$，因为$\delta_x > \delta_y \land y_1 \ngtr y_2$。
3. $s_1 \rightarrow s_3$，其中$p_{\mathrm{IG}}(s_1, s_3) = 0$，因为$\delta_y > \delta_x \land x_1 \ngtr x_2$，以及$p_{\mathrm{IG}}(s_3, s_1) = 1$，因为$\delta_y > \delta_x \land x_1 > x_2$。
在数字游戏的策略集合$\boldsymbol{S}$的结构中，基于策略间的成对关系所形成的几何组织可以被适当地利用，以支持某种形式的可视化分析[8, 9]。一个偏序集（poset）$(\boldsymbol{S}, \preceq)$ 用来比较$\boldsymbol{S}$中的元素，并根据二元关系$\preceq$对它们进行排序。（注意：通过禁止自比较，例如对于二人对战游戏指定的关系中禁止自体对战，可以得到一个简化的偏序集$(\boldsymbol{S}, \prec)$。）我们考虑一个非传递游戏，其中$\boldsymbol{S} = \{s_1, s_2, s_3\}$。如果所有游戏都在策略三元组$s_1, s_2, s_3$间进行，并且支配关系$\leftarrow$满足传递性，那么偏序集$(\boldsymbol{S}, \leftarrow)$是完全有序或线性有序的，其关系可以单调地表示为$s_1 \leftarrow s_2 \leftarrow s_3$，因为根据传递性有$s_1 \leftarrow s_3$。非正式地，通过可视化方式，$\boldsymbol{S}$的元素可以作为一维线上的不同点，其排序由点的位置指定（例如，可以按照递增顺序表示为$s_3 \rightarrow s_2 \rightarrow s_1$）。对于之前的例子$\boldsymbol{S}$，其中$s_1 = (1,6), s_2 = (4,5), s_3 = (2,4)$，由于非传递性，策略三元组形成了一个不可比较的循环。因此，这个三元组无法嵌入到一维线中进行线性排序。然而，这个$\boldsymbol{S}$可以嵌入到$\mathbb{N}^2$中，而$\mathbb{N}^2$本身可以进一步嵌入到二维空间$\mathbb{R}^2$中（因此可以在二维空间中可视化）[9]。对于一个策略$s_i = (x, y)$，每个维度上的值可以分别表示为$s_i.x$和$s_i.y$。

在当前例子中，$(\boldsymbol{S}, \leftarrow)$基于$s.x$维度有线性排序$s_1 \rightarrow s_3 \rightarrow s_2$，同时基于$s.y$维度有另一个线性排序$s_3 \rightarrow s_2 \rightarrow s_1$。 $(\boldsymbol{S}, \leftarrow)$嵌入在具有两个线性排序的二维空间中的这个事实表明，在$\boldsymbol{S}$上的成对关系中存在循环。假设现在在策略集合中添加一个新的策略$s_4 = (s_1.x + 10, s_1.y + 10) = (11, 16)$。对于新的策略集合$\boldsymbol{S} \cup \{s_4\}$，可以像以前一样可视化子集$\{s_1, s_2, s_3\}$，而后将策略$s_4$定位于二维空间的右上方。事实上，通过直接计算，可以确定策略$s_4$支配了子集$\{s_1, s_2, s_3\}$中的所有策略。然而，也可以利用$\{s_1, s_2, s_3\}$形成循环并嵌入于$\mathbb{R}^2$的几何结构，以及游戏规则的知识来断定策略$s_4$由于在每个维度上具有更高的值而支配了其他策略，即$s_4.x > \max \{s_1.x, s_2.x, s_3.x\}$和$s_4.y > \max \{s_1.y, s_2.y, s_3.y\}$。特别地，策略$s_4$不是仅在某一个维度上寻求更大的值，而是在两个维度上同时采用更大的值来支配剩余策略。这一点通过策略$s_4$相对于子集$\{s_1, s_2, s_3\}$的位置以及其位于远离$\mathbb{R}^2$原点对角线方向上的几何位置得到体现。
### 5.2.2 循环与遗忘

在我们目前的讨论中，我们使用**循环**（cycle）这个术语来指代二元关系中存在的非传递性现象，而不是那些满足传递性性质的关系。具体而言，循环出现在策略组 $s_i, s_j, s_k \in \boldsymbol{S}$ 内的非传递性关系中。这一现象与协同进化系统（co-evolutionary systems）中的种群动力学研究相关。相关研究重点在于通过反复应用变异算子（variation operator）和选择算子（selection operator），种群的配置状态如何发生变化。这里主要讨论的病理性问题就是**循环种群动力学**（cycling population dynamics），简称**循环**（cycling）[11, 14, 22, 30, 55]。

在第 1 章的 1.2 节中，我们描述了一个离散时间的动态系统 $\mathcal{F}: \boldsymbol{\mathcal{X}} \rightarrow \boldsymbol{\mathcal{X}}$，其中包括构成协同进化过程的变异算子 $\mathcal{V}$ 和选择算子 $\mathcal{S}$，这两个算子共同定义了协同进化过程 $\mathcal{F} = \mathcal{S} \circ \mathcal{V}$，其作用空间是种群配置空间 $\boldsymbol{\mathcal{X}}$。这一数学抽象和建模方法同样被应用到了协同进化研究中[35, 62]。通过在配置空间 $\boldsymbol{\mathcal{X}}$ 中生成轨迹（trajectories）作为迭代形式 $\boldsymbol{X}_1, \boldsymbol{X}_2, \boldsymbol{X}_3, \ldots$，这种协同进化动态系统能够分析复杂的协同进化种群动力学。

具体而言，周期循环的动态机制（regimes of periodic cycles）被用于更精确地描述协同进化中的循环种群动力学。在这里，周期为 $n \geq 1$ 的周期循环是指种群配置或状态的轨迹或序列 $\boldsymbol{X}_1, \boldsymbol{X}_2, \boldsymbol{X}_3, \ldots, \boldsymbol{X}_n$ 会重复，例如 $\boldsymbol{X}_n = \boldsymbol{X}_1$。在后续的 **5.3 节**中，我们将给出这种协同进化动态系统的精确理论表述，以便进行严谨的种群动力学分析。

在本节的剩余部分中，我们首先介绍并讨论与协同进化过程中出现的病理性搜索相关的动态机制。需要明确区分以下两种**循环结构**：一种是种群动力学中，由协同进化过程 $\mathcal{F}$ 在种群配置空间 $\boldsymbol{\mathcal{X}}$ 上产生的循环结构；另一种是关于策略组 $\boldsymbol{S}$ 中双向关系 $\leftarrow$ 所产生的循环结构。前者描述与协同进化过程相关的结构，而后者描述与协同进化问题相关的结构。

尽管如此，这两种循环结构之间存在联系。协同进化过程中的选择机制通过种群中的代理实施策略，生成一系列交互结果，而这些交互结果又受到策略组 $\boldsymbol{S}$ 内双向关系的底层结构支配 [19]。例如，对于“石头－剪刀－布”游戏，存在非传递性时，其抽象协同进化动态系统的相图（phase portrait，或称相位图）会指示出循环动力学。例如，研究 [30] 生成了这种相图，并以轨迹的曲线图形式表示。这些轨迹位于嵌入在 ${{\mathbb{R}}^{2 + 1}}$ 空间中的二维标准单形（Simplex） ${{\mathbb{S}}^{2}}$ 中。（非正式地说，减少一个维度的原因是，由于策略比例总和为 1，已知两种策略的比例就可以自动确定剩余策略的比例。）

**相空间**或种群配置空间 $\boldsymbol{\mathcal{X}}$ 是由 ${{\mathbb{S}}^{2}}$ 表示的，包含三种策略“石头－剪刀－布”的具体混合或比例分布。
周期性循环动力学的产生是由于相空间中存在固定点，这些固定点代表特定的种群配置处于不稳定的平衡状态并充当排斥点。这种动力学表明，相邻区域中的点（例如，由于突变算子作用产生的种群配置的微小扰动）的轨迹会远离这些排斥固定点。具有三种纯策略供代理选择实施的协进化动力系统可能具有固定点，这些固定点表示单型和多型平衡，它们分别代表只实现一种策略的种群配置和至少实现两种策略的种群配置。

对于“石头、剪刀、布”游戏，其协进化动力系统存在一个不稳定的多型平衡，该平衡中石头、剪刀和布策略的比例相等，以及三个单型平衡，分别对应所有代理都选择石头、剪刀或布策略。在图形上，简单形$\mathbb{S}^{2}$是一个等边三角形，每个角包含一个鞍点，而其中心位置为不稳定的多型平衡点。该中心点也对应纳什平衡，因为在一个代理种群中，如果所有代理均以石头、剪刀和布策略的相等比例进行游戏，在完全混合条件下（即每个代理与所有代理进行游戏，包括与自己对局），他们将获得相同的游戏收益。这三个角上的鞍点则因“石头、剪刀、布”策略之间的非传递关系而存在。鞍点的动力学既有吸引成分也有排斥成分。如果将动力学限制在简单形$\mathbb{S}^{2}$的各边上，则任何非角点的种群配置对应的代理只实现三种可选策略中的两种，其配置会被吸引到某一角点并从相对角点排斥离开。例如，任何非角点会被吸引到代表所有代理选择剪刀策略的角点，同时从代表所有代理选择石头策略的角点被排斥。因此，对于从简单形$\mathbb{S}^{2}$中非平衡点开始的轨迹，它们将在该相空间的区域中循环运动。通常，多个具有排斥动力学的不稳定固定点的存在会导致各种循环动力学，包括周期性循环动力学。因此，即便是存在仅两种纯策略的协进化动力系统，如“鹰-鸽”游戏系统，也可能产生循环动力学，因为游戏结构导致多型和单型平衡均不稳定。此外，通过对协进化中种群配置的单值测量（$S^n \rightarrow \mathbb{R}_{\geq 0}$）以及结合游戏的领域知识，例如涉及更复杂游戏（如迭代囚徒困境 (IPD)）的协进化也可以观察到循环动力学的存在。

在协进化动力学中，“遗忘病理”（forgetting pathology）指的是种群配置随着时间的变化，其中某些早期出现的行为特征消失，随后又被重新发现。这种现象有三种可能的情况 [30]：第一种情况是该行为特征在协进化过程中遭到选择性淘汰，因为实现该特征的策略会使代理获得较低的适应度值；后两种情况中，这种特征则被选择性保留下来。然而，由于变异算子的作用，选择压力可能较低，此时可能出现漂变（drift）；或者在更极端的情况下，由于变异算子机制本身存在对该特征的偏见，导致其无法在产生的后代中表现出来。在讨论行为特征时，可以考虑将纯策略直接视为完整体现这些特征的更简单设定。例如“石头、剪刀、布”游戏，协进化动力系统可以通过引入变异算子带来的特定扰动来展示表现“遗忘病理”的循环动力学。这种情况是指扰动轨迹在简单形$\mathbb{S}^{2}$的边上循环运动，使得周期性地有一种策略没有被协进化的种群所实施。然而，行为特征也可能涉及某些策略或策略子集在对抗入侵对手策略时采用的游戏机制。例如在数字游戏中，一个特征可以对应一个特定维度；在另一个例子“迭代囚徒困境”游戏中，行为特征可能指策略在游戏中采用某些特定选择（如合作、背叛或某些更复杂的游戏机制，如中间水平的合作以及合作与背叛的组合）的倾向。因此，对于实施中间水平合作的IPD游戏策略的协进化过程，循环动力学可能与行为特征的遗忘问题相一致。
### 5.2.3 脱离与平庸稳定状态  

在这一子节中，我们介绍并讨论与协同演化过程相关的另外两种主要病态现象。这两种病态是文献中大多数研究通过实证调查所探讨的对象。此外，大部分研究针对这两种病态现象所进行的调查通常发生在更复杂的双种群竞争协同演化环境中，这与之前所讨论的两组病态现象不同，后者通常发生在较简单的单种群竞争协同演化环境中。尽管如此，通过受控的实证研究，依然能够揭示这一最后一组协同演化病态现象并分离出其所依赖的条件。对于这些病态现象的某些解释来自其他研究领域，例如演化生物学。  

“脱离”（Disengagement）病态现象指的是在双种群竞争协同演化中，由于存在不对称交互（例如捕食者-猎物协同演化），导致一个种群的大部分表现普遍优于另一个种群的情形 $[11]$。当出现这一情况时，相对适应度评估（Relative Fitness Evaluation）便失去意义，因为对立种群无法有效区分当前种群的表现。一种典型的案例是捕食者-猎物协同演化，这通常是指机器人在追逐-逃避交互中的适应性。这种模拟对应着自然界中类似现象的演化机器人研究 $[21]$。对此，“红皇后效应”（Red Queen Effect）被提出，用于解释两个协同演化种群在演化时间线中持续的“军备竞赛”：两个种群轮流通过进化响应去侵蚀并抵消对立种群所发现的赋予适应性优势的性状，而这种适应性优势性状会促使当前种群产生反制该性状的响应性适应 $[64]$。只要在适应阶段一个种群具有足够的多样性，并在另一种群中实现反制适应，这种动态交互就可以在耦合的双种群之间持续存在。  
然而，当前种群成员可能会发现某一具有内在优势的特定适应性特征，从而相对于对立种群获得显著优势，并通过繁殖在整个种群中迅速扩散，这种扩散可能在一个短暂的演化周期内完成。这将导致脱离现象的发生，因为当前种群的大多数成员在表现上优于对立种群的成员。当这种情况出现时，相对适应度（relative fitness）将不再能够对种群成员进行区分以用于选择，它们可能开始发生基因漂移（例如通过变异操作符）。此外，同样的动态可能也会在对立种群中发生。在这种情况下，对立种群的成员多数表现较差并获得相等糟糕的适应度值，而这些适应度值无法用于区分成员。因此，对立种群也开始发生漂移。当这种情况发生时，两种协演化种群就被认为处于完全脱离（full disengagement）的状态。[11]中还识别出了两种其他程度的脱离现象。

**非对称脱离（Asymmetric disengagement）**发生在一个种群达到全局最优解，而对立种群开始基因漂移的情况下。以主机-寄生虫竞争协演化用于排序网络的例子[44]为例，[11]的一项研究对这种非对称脱离的实际情景及其对搜索性能的损害进行了阐释。需要注意的是，主机-寄生虫竞争协演化除了主要目标寻找最优排序网络外，还可能具备次要目标——寻找最佳测试序列。假设问题环境本身也表现出非对称特性，例如最佳排序网络能够排列难度最大的测试序列。如果排序网络种群演化到发现全局最优解，则非对称脱离将开始出现。全局最优解将不断扩散，而在测试序列种群中，基因漂移可能发生并停止寻找更困难的测试案例。在这种情况下，协演化系统可以实现主要目标，但无法实现次要目标。此外，如果协演化系统继续运行，持续的基因漂移可能会向排序网络种群引入有害突变，在选择操作符存在非精英主义机制的情况下，这可能导致它们丧失甚至原本的最优解。

**共谋脱离（Collusive disengagement）**发生在两个协演化种群之间，当其中一个种群进入某种特定的稳定配置或成员平衡状态，即便它们的表现是次优的。然而，这种种群整体被认为具有共谋现象，因为其成员具有类似的适应度值。[11]的一项研究提供了一个关于两种群协演化系统中的硬币投掷游戏的例子，该系统由生成器（generators）和预测器（predictors）种群组成。此处的最优解是生成器能够以 $50:50$ 比率随机产生正面与反面，而预测器则以相等的概率预测这两个可能结果。然而，也可能出现这样的情况：预测器种群中总是预测正面和总是预测反面的成员各占一半比例。虽然这些预测器作为个体来说是“专家”，拥有几乎相同的平均适应度或表现，但作为一个整体种群，它们形成了一种策略混合，能够在对抗多样化生成器种群时表现出最佳效果。后者种群将由于非对称脱离开始基因漂移。该种群配置的有效性意味着演化过程难以探索实际的最优个体预测器。种群中会存在强大的选择压力维持这样的种群配置或状态，在这种情况下，这种状态也被称为中等稳定状态的病态问题（pathological issue of mediocre stable states）[31]。

实际上，生成器种群也可能进入这种状态，即总是抛正面和总是抛反面的成员各占一半比例。如此，生成器和预测器的两种种群都被认为有共谋现象。
5.3 共演搜索分析  
本节将介绍两个用于分析竞争性共演系统的正式框架。这些理论框架能够准确描述共演系统用于执行搜索时的显著特征。它们揭示了多个与共演地病（pathologies）相关的重要结构，这些病理可能会在应用共演系统解决问题时出现。更重要的是，这两个框架提供了构建工具的手段，从而能够对共演系统进行严格分析。第5.3.1节将首先介绍基于**有向图理论**的竞争性共演问题结构分析的正式框架。接下来，第5.3.2节将介绍基于**动力系统理论**的竞争性共演过程分析的正式框架。

### 5.3.1 共演问题中的循环结构  
最近，我们 [19] 引入了一个框架，该框架利用**有向图（digraphs）**来捕捉任何共演问题中的整体结构，这些问题表现为具有确定性结果的两人策略游戏，这是竞争性共演环境中的典型问题 [1, 13, 15, 25, 36, 40]。该共演问题被表示为共演**有向图**$D_C = (V_S, A_R)$。  
顶点集（vertex set）表示解集（solution set，即策略集）$V_S = S$，其中每个独特的顶点$v \in V_S$对应一个独特的解（即纯策略）。  
弧集（arc set）捕获了解集$S$上的二元关系，该关系定义为$A_R = R$，其中每条弧（有向边）可以是单向的$u \leftarrow v$ 或 $u \rightarrow v$，或双向的$u \leftrightarrows v$，它表示解对$u, v \in V_S$之间的偏好关系（pairwise preference relationship）。  

通过这种方式，各种竞争性共演问题作为两人策略游戏，可以通过这些有向图来表示，从而完全捕获它们在对应策略集上的成对关系（优势关系）结构。
策略对之间的博弈通过有向图中的边 $\{u, v\}$ 将两个顶点 $u, v \in V_S$ 连接起来来捕获。注意，底层图 $D_C$ 是完全图，因此任何两个顶点都存在边连接，即任何策略对都存在博弈行为。博弈的胜负结果通过边的方向来表示，例如弧 $u \leftarrow v$ 表示 $u$ 支配（击败）$v$。以这种方式，协同进化锦标赛通过具有胜负结果的双人博弈模型化协同进化问题，其中给定特定策略集下所有策略之间相互竞争（自我博弈除外）。如果结果包括平局（即 $u \leftrightarrows v$），则这些带有胜负平局结果的协同进化问题被包含在协同进化半完整有向图的范畴中。接下来，在本小节的讨论中，我们将重点聚焦于协同进化锦标赛的范畴，这一范畴足以呈现有向图理论中的循环结构，这对于协同进化问题结构的完整刻画至关重要。我们在后续讨论中提供了有向图的一些概念和属性。

支配的概念可以扩展到 $V_S$ 的子集，其中顶点可以看作是由一个元素组成的 $V_S$ 的子集。对于两个不相交的子集 $V_S^1$ 和 $V_S^2$（满足 $V_S^1 \cap V_S^2 = \emptyset$），如果满足 $V_S^2 \Leftarrow V_S^1$，则意味着 $V_S^2$ 中的所有顶点 $v \in V_S^2$ 都支配 $V_S^1$ 中的所有顶点 $u \in V_S^1$。

在有向图 $D_C = (V_S, A_R)$ 中的**行走**定义为顶点和弧交替出现的序列 $v_1a_1v_2a_2v_3\ldots a_{k-1}v_k$，其中 $v_1, \ldots, v_k \in V_S$，且 $a_1, \ldots, a_{k-1} \in A_R$。通常该行走可以简单写为顶点序列 $v_1v_2v_3\ldots v_k$。行走的长度表示其中弧的数量，在上述定义中即为 $k$。一个行走是**闭合的**，如果满足 $v_k = v_1$。接着，一个 $(v_1, v_k)$-路径可以定义为行走，其中所有顶点均不重复。于是，一个 $k$-循环被定义为长度 $k \geq 3$ 的闭合 $(v_1, v_k)$-行走，该行走构成了 $(v_1, v_{k-1})$-路径。

如果对于 $V_S$ 中的每一对不同顶点 $u, v \in V_S$ 都分别存在 $(u, v)$-路径和 $(v, u)$-路径，则该有向图被称为**强连通**（简称强连通）。一个**可约有向图**满足其顶点集 $V_S$ 可以划分为两个互不相交且非空的子集 $V_S^1$ 和 $V_S^2$，满足 $V_S^1 \cap V_S^2 = \emptyset$ 且 $V_S^1 \cup V_S^2 = V_S$，并且进一步的性质是两个子集的方向是单向的，即满足 $V_S^1 \Rightarrow V_S^2 \lor V_S^1 \Leftarrow V_S^2$。否则，该有向图被称为不可约的。实际上，所有不可约的有向图根据定义都是强连通的。注意，可约有向图 $D$ 的顶点划分会引导其生成子有向图，例如 $D^1$ 和 $D^2$。我们用 $V(D^{(i)})$ 表示子有向图 $D^{(i)}$ 的顶点集，用 $A(D^{(i)})$ 表示其弧集。

现在，我们可以使用有向图理论的语言准确地描述游戏的一些特性，这些特性用于刻画与协同进化病理相关的循环结构，这些结构在过去通常是非正式解释的。例如，**石头、剪刀、布**游戏中的结构可以通过协同进化有向图 $D_C = (V_S, A_R)$ 捕获，其中 $V_S = \{v_1,v_2,v_3\} = \{\mathrm{Rock},\mathrm{Paper},\mathrm{Scissors}\}$ 且 $A_R = \{v_1 \rightarrow v_2, v_2 \rightarrow v_3, v_3 \rightarrow v_1\}$。对于只有三种策略的游戏，该三元组即构成策略集。这个三元组中的支配关系的不传递性正好被 $(v_1,v_3)$-路径上的 3-循环 $v_1v_2v_3v_1$ 精确捕获。由于闭合的循环和路径都包含 $V_S$ 中的所有策略，它们被称为**哈密顿结构**。
在另一个关于一维数字游戏的示例中，假设 $VS = \{v_1,v_2,v_3\} = \{1,2,3\}$ 且 $AR = \{v_1 \rightarrow v_2, v_2 \rightarrow v_3, v_1 \rightarrow v_3\}$，该三元组满足传递性性质。这里，协演化有向图 $D_C = (V_S, A_R)$ 是传递的。令 $v_1, v_2, v_3, \ldots, v_n$ 是 $V_S$ 的一个排序，使得对于所有弧 $v_i \rightarrow v_j \in A_R$ 都有 $i < j$。这样的有向图被称为具有无环排序的图，并且无环图没有循环。实际上，已知一个锦标赛图只有在无环的情况下才是传递的 [5]，因为可以将所有支配关系单调地表达为 $v_1 \rightarrow v_2 \rightarrow v_3 \rightarrow \cdots \rightarrow v_n$。

石头、剪刀、布游戏以及只有三个策略的一维数字游戏现在可以被视为通过有向图理论语言描述的简单例子。然而，通过循环结构和无环排序来表征锦标赛的概念可以推广到更大的锦标赛图。特别地，用于得到这些表征的理论结果利用了这样的事实：非哈密顿低阶循环可以被打包到 $V_S$ 的子集作为强连通分量，并且在这些强连通分量上的支配关系是单调的。实际上，使用有向图表示法的动机不仅在于它自然地捕获了协演化问题底层的所有相关结构，还在于可以使用有向图理论及其中特别是锦标赛理论 [50] 中的各种工具和强有力的结果来对其进行严格分析。

具体地，我们 [19] 已将不可约性的有向图理论概念与“可解性”的博弈理论概念联系起来，以对这些协演化锦标赛进行定性表征。实际上，一个主要结果是协演化锦标赛要么是可解的（即存在一个支配的策略子集），要么不是。一个可解的协演化锦标赛 $T_C = (V_S, A_R)$ 是可约的，即策略集合可以通过顶点划分出一个支配子集，例如 $V_S^1 \Rightarrow V_S^2$。协演化锦标赛中的循环结构决定了其是否可解（可约）。一个主要的技术结果是证明一个可约的协演化锦标赛 $T$ 可以被强分解为子图 $T^{(1)}, T^{(2)}, T^{(3)}, \ldots, T^{(l)}$，其中 $\cup_{i = 1}^{l} V\big(T^{(i)}\big) = V(T)$ 且当 $i \neq j$ 时 $V\big(T^{(i)}\big) \cap V\big(T^{(j)}\big) = \emptyset$。此外，每个子图都是一个强连通分量，并且它们形成一个支配关系链，可以单调地表达为 $V\big(T^{(1)}\big) \Rightarrow V\big(T^{(2)}\big) \Rightarrow V\big(T^{(3)}\big) \Rightarrow \cdots \Rightarrow V\big(T^{(l)}\big)$。

另一个重要的技术结果是不可约的协演化锦标赛 $T$ 是强连通的，并且顶点是泛循环的，即对于每个 $v \in V(T)$ 和每个 $k \in \{3, 4, 5, \ldots, n\}$，在锦标赛 $T$ 中包含 $v$ 的 $k$-循环皆存在。因此，我们能够利用有向图理论分析来揭示协演化问题中出现的所有循环结构。我们的分析能够揭示这些循环结构比以往文献研究显示的更复杂。

例如，一个这样的发现是低阶循环包含在各种程度的高阶循环中，即便在简单的一种群体协演化设置中，这些设置被建模为马尔可夫链，可能会影响搜索性能 [19]。此外，这些结果与其他理论研究一致，例如 [49] 中具体探讨的情境，即在协演化进化算法（CEA）中若构造出一个具有关键属性的合适目标度量，则能产生类似于进化算法（EA）的搜索行为。例如，一个具有赢-输结果（内部主观性能度量）的协演化问题，其解决方案构成一个传递链，其排名等同于由解决方案得分序列所给出的外部客观性能度量，其中解决方案的得分为其胜场数量。
我们将以几个共演竞赛的实例描述来结束本小节，以突出关于这些循环结构可以获得的关键特性。通过对一个传递竞赛（transitive tournament）中的单一弧的反向操作，可以很容易获得一个顶点全周期共演竞赛（vertex pancyclic coevolutionary tournament）。对于按照无环排序（acyclic ordering）标记顶点的传递竞赛，其顶点为$v_1, v_2, v_3, \ldots, v_n$，其中一端$v_1$表示表现最差的策略（无任何胜利），另一端$v_n$表示表现最佳的策略（获胜次数为最大值$n - 1$）。在顶点$v_1, v_2, v_n$之间具有以下优势关系：$v_1 \rightarrow v_2, v_2 \rightarrow v_n, v_1 \rightarrow v_n$。只需将弧$v_1 \rightarrow v_n$反向为$v_1 \leftarrow v_n$即可生成一个顶点全周期竞赛$T$。实际上，在竞赛$T$中，存在一条哈密顿循环（Hamiltonian cycle）：$v_1v_2v_3{\ldots}v_nv_1$，基于哈密顿路径$v_1v_2v_3{\ldots}v_n$。这种顶点全周期竞赛具有最少数量的三元循环（3-cycles）[10]。此外，一个规则竞赛（regular tournament），其中每个顶点均具有相同数量的入弧（表示该顶点获胜的游戏数量）和出弧（表示该顶点失败的游戏数量），则包含最多数量的三元循环。只有顶点数为奇数的竞赛可以是规则竞赛，并且它是最小规则竞赛（仅由三个顶点组成）的推广，例如“石头、剪刀、布”游戏中的非正式示例。一般而言，与包含最少数量三元循环的竞赛子组件相比，共演过程需要更多时间才能从包含最多数量三元循环的规则子组件中逃逸[19]。

### 5.3.2 共演过程中的种群动态

我们首先将介绍竞争性共演的具体数学表述，将其定义为离散时间动力系统，然后在这种数学框架内讨论生成于共演过程中的种群动态所涉及的复杂循环结构。研究[35]首次正式定义了共演的动力系统，尽管还有其他相关研究，例如较早的研究[38]以及后续的研究[39]，通过计算机生成的轨迹经验性调查了共演的种群动态。我们已经在第5.2.2节中引入了种群动态的概念，描述其为轨迹集合，每一条轨迹是特定的正向序列$\{X_n : n \in \mathbb{N}_{\geq 1}\} = \{ F_0, X_1, (F_1), X_1, (F_2), X_1, \ldots \}$。
\[
\{{\boldsymbol X}_{n} : n \in {\mathbb{N}}_{\geq 1}\} = \big\{\mathcal{F}^0\big({\boldsymbol X}_{1}\big),\mathcal{F}^1\big({\boldsymbol X}_{1}\big),\mathcal{F}^2\big({\boldsymbol X}_{1}\big),\ldots\big\}
\]
表示了通过协同演化过程 $\mathcal{F}: \boldsymbol{\mathcal{X}} \rightarrow \boldsymbol{\mathcal{X}}$ 从初始种群状态 ${\boldsymbol X}_{1} \in \boldsymbol{\mathcal{X}}$ 生成的一系列种群配置。采用动态系统作为协同演化的抽象数学模型有两方面的动机：(1)在某些设置下，协同演化过程的机制可以得到合适的形式化。(2)无论这些轨迹表现为简单还是复杂，其行为或模式都可以进行定性表征并进行精确描述。然而，在高维度度量空间 $\boldsymbol{\mathcal{X}} = ({\mathbb{R}}^D,d)$ 中，由动态系统 $\mathcal{F}$ 引发的种群动态十分丰富，该空间通常用于描述在来自大量纯策略集 $(D-1)$ 中选择的代理之间的竞争性协同演化交互。实际上，正如我们之前所强调的，即使仅存在两个纯策略的情形，协同演化也能生成复杂的种群动态。然而，我们仍然可以为这种系统提供易于理解的定性和定量表征。因此，大多数研究[35, 38, 39, 61, 62]集中于涉及两玩家、两纯策略游戏（如鹰鸽博弈）的协同演化场景。此外，动态系统分析也已在合作协同演化设置中得到应用[53, 54]，尽管该设置更为复杂，因此主要依赖计算机生成的轨迹作为分析的数据来源。

接下来，我们的讨论将重点放在涉及两玩家、两纯策略游戏的协同演化动态系统上，其模型将在下文中正式描述。这类系统属于一种一维离散时间动态系统家族，该系统由度量空间上的函数 $f: \boldsymbol{\mathcal{X}} \rightarrow \boldsymbol{\mathcal{X}}$ 指定，距离度量 $d$ 定义在实线的紧区间 $([0,1],d)$ 上。随后，我们将采用简化记号 $f: [0,1] \rightarrow [0,1]$。映射 $f$ 在区间 $[0,1]$ 上的作用生成一个轨迹（路径），该轨迹通过差分方程 
\[
x_{n + 1} = f(x_n) \in [0,1], \ n \in {\mathbb{N}}_{\geq 1}
\]
进行描述，初始条件为 $x_1 \in [0,1]$。点 $x_1$ 的真实轨迹是序列 
\[
\{x_n : n \in {\mathbb{N}}_{\geq 1}\} = \{x_1,f(x_1),f^2(x_1),\ldots\},
\]
其中函数的迭代 $f^m, \ m \in {\mathbb{N}}$ 定义为 $f^0(x_1) = x_1$ 和 $f^n(x_1) = f\big(f^{n-1}(x_1)\big)$。

在这样的设置中，[35] 的主要贡献是通过进化博弈理论(EGT)来识别竞争性协同演化和动态系统之间的主要建模接口。协同演化过程的模型假定：(1)存在无限数量的代理种群（每个代理实施两种游戏纯策略之一）。(2)每个代理与种群中的所有其他代理互动（即完全混合），并从游戏结果中积累收益。(3)代理根据累积的收益进行无性繁殖（仅生产克隆）。以这种方式，不同的游戏规格和选择机制会导出不同构建的复制器映射 $f: [0,1] \rightarrow [0,1]$ [45, 57]。

我们考虑的交互为两玩家、两纯策略的对称游戏，其游戏为同时进行的一次性交互，结果通过规范形式收益矩阵来指定，例如
\[
{\mathbf{A}} = (a_{ij} : 1 \leq i,j \leq 2)
\]
，其中每个项 $a_{ij}$ 为第一玩家指派一个非负收益。令 $p_1,p_2 \in [0,1]$ 分别表示无限种群中两策略的比例。由于策略是可遗传的，并且在任意世代（或迭代轮次）$t$ 中，所有代理为了获得决定其平均繁殖成功率的收益而竞争，复制器映射 $f$ 决定了种群状态（即 $p_1,p_2$）随时间 $t$ 的变化。为了构造 $f$，我们收集状态中的策略比例形成一个列向量 
\[
\mathbf{p}(t) = \big(p_1(t),p_2(t)\big){^{\mkern-1.5mu\mathsf{T}}},
\]
以及与之关联的策略累计收益向量 
\[
\mathbf{w}(t) = \big(w_1(t),w_2(t)\big){^{\mkern-1.5mu\mathsf{T}}} = {\mathbf{A}} \mathbf{p}(t),
\]
在每次迭代 $t$ 中满足上述关系。令 $p_1 = p$ 和 $p_2 = 1 - p_1 = 1 - p$，则：

$$w_1 = p_1a_{11} + p_2a_{12} = pa_{11} + (1 - p)a_{12}$$  
$$w_2 = p_1a_{21} + p_2a_{22} = pa_{21} + (1 - p)a_{22}$$  

复制因子为  
$$f(p_i) = \frac{p_i w_i}{\sum_{i=1}^2 p_i w_i}$$  
使用适应度比例选择策略。对于$p_1$，有：  
$$f(p) = \frac{p \cdot w_1}{p \cdot w_1 + (1 - p) \cdot w_2}$$  

请注意，共演动态系统是一维的，因为由$f$定义的动态在嵌入到$\mathbb{R}^2$中的标准一维单纯形$\mathbb{S}$上演化，且满足$\sum_{i=1}^2 p_i(t) = 1$，对于所有$t = 0, 1, 2, \ldots$。 由  
$$\mathbf{p}(t+1) = f\big(\mathbf{p}(t)\big$$  
描述的动态在使用以下形式代替时仍保持不变：  
$$p(t+1) = f\big(p(t)\big$$  
即观察时间中$p = p_1$的变化，因为$\mathbb{R}^2$的坐标轴中与$p_1$相关的线段与紧区间$[0,1]$的投影是一致的。  

通过构造与研究不同的复制因子映射，可以建模不同的选择机制[35, 62]。从非正式的角度看，共演动态系统可以作为一种模型来研究选择作用下的长期种群动态。变异的作用（如突变）可以在特定有限扰动的轨迹中进行考虑，同时结合其他对种群状态产生影响的因素，如有限种群效应[61,62]。更重要的是，现在可以将这些共演动态系统表述为特定的复制因子映射$f: [0, 1] \rightarrow [0, 1]$，这使得能够使用动态系统理论广泛的表达语言来精确研究和描述由这些映射生成的种群动态。  

设$\mathcal{O}_f = \{x_n: n \in \mathbb{N}\}$为由$f$生成的轨道。复制因子映射$f$所生成轨道的各种动态行为或机制将在此介绍。一种简单机制涉及到$f$的固定点$f(q) = q \in [0,1]$，我们用$\boldsymbol{\mathcal{O}}_f = \{q\}$来表示$f$的固定点轨道以示强调。此概念可以扩展到$f^k$的固定点，以描述周期为$k$的周期轨道，定义为$k$个不同点组成的序列：  
$$\boldsymbol{\mathcal{O}}_{f^k} = \{q_1, q_2, q_3, \ldots, q_k\} = \big\{f^0(q_1), f^1(q_1), f^2(q_1), \ldots, f^{k-1}(q_1)\big\}.$$  

请注意，这假设一种极小性条件，使得对于每个$q \in \boldsymbol{\mathcal{O}}_{f^k}$均满足$f^j(q) \neq q, 1 \leq j \leq k - 1$。尤其是，如果$f^k(q_i) = q_i \in (0,1)$对于所有$q_i \in \boldsymbol{\mathcal{O}}_{f^k}$成立，则$\boldsymbol{\mathcal{O}}_{f^k}$称为周期。如果系统从周期点$q_i \in \boldsymbol{\mathcal{O}}_{f^k}$的周期为$k$开始，则轨道将始终固定在$\boldsymbol{\mathcal{O}}_{f^k}$中（以$\boldsymbol{\mathcal{O}}_{f^k} = \{q_1, q_2, q_3, \ldots, q_k\}$形式无限循环，且$k \geq 2$）。  
动态系统理论精确描述了另一个概念，即点$p$的轨迹（orbit）累积到固定点（fixed point）[7, 28]。对于真实世界系统的测量或计算模拟中，经常观察到轨迹的长期行为——某些轨迹进入并停留在固定点的现象。特别地，周期为$k$的周期轨迹$\boldsymbol{\mathcal{O}}_{f^k}$被称为吸引的（attracting），如果点$p$在周期点附近的轨迹累积到该轨迹上，即$f^j(p) \rightarrow \boldsymbol{\mathcal{O}}_{f^k}$，当$j \rightarrow \infty$时[28]。因此，周期$k \geq 2$的吸引性周期轨迹也被称为周期为$k$的极限环（limit cycle），这一现象曾在竞争性共同演化（competitive coevolution）的研究中被报道[35, 38]。这种极限环的一个典型特征是，通过测量或计算生成的种群状态的轨迹，最初在状态空间$[0, 1]$中游走，但随后迅速进入周期性，因为它们进入了某一周期点。

除了这些被研究文献充分记录的吸引固定点和周期轨迹的状态，共同演化种群动态中还观察到更复杂的行为，这些被归因于混沌动态（chaotic dynamics）[35, 38]。然而，当前共同演化动态的文献中对混沌动态的描述仍不够精确。部分原因在于，大多数关于共同演化动态的研究是基于计算生成轨迹并且在短时间尺度上进行的。在文献[35]中，通过$10^4$次迭代的种群动态被归因于混沌，其实这些动态是计算生成的伪轨迹（pseudo-orbits），受到有限计算精度的误差影响。然而众所周知，混沌动态对初始条件非常敏感，相邻点的轨迹可以指数级分离（例如，初始间距为$10^{-14}$，典型的误差倍增率意味着大约经过50次迭代后会有一个接近$2^{50} \times 10^{-14}$的单位误差）[43]。但是，在动态系统的研究中，也有文献表明，通过混沌映射计算生成的伪轨迹可以可靠地跟踪超过$10^5$次迭代[20]。我们[61, 62]已经开始研究基于连续复制子映射（replicator maps）$f$的共同演化混沌动态，这种映射可能具有双曲性（hyperbolicity）。这是关键，因为双曲映射在结构上是稳定的，其轨迹结构丰富且稳健。尽管对初始条件有敏感依赖，但双曲映射具有“跟踪性”（shadowing property），这使得任何由于每次迭代的有限扰动而生成的伪轨迹都可以被真实轨迹紧密跟踪。在共同演化问题中使用图论（digraph theory）研究循环结构的同时，引入动态系统理论中的更多方法，能够严格分析共同演化系统中可能出现的复杂和混沌动态的精密结构，并且在特定情况下（例如跟踪动态）可以利用计算机辅助的工具。事实上，混沌动态中的这两个看似对立的特征实际上为研究具有非常复杂行为的动态系统提供了可能。

除了混沌轨迹以看似不规则和不可预测的方式穿越种群状态空间（例如对于共同演化动态系统，单位区间$[0, 1]$），它们还会接近其起始点。这种行为源于混沌映射在拓扑上具有可传递性（topologically transitive），即对于任意成对的不为空的开子集$U, V \subset X$，总存在一个$n > 0$使得$f^n(U) \cap V \neq \emptyset$[29]。一种非正式的方法来理解混沌映射$f: [0, 1] \rightarrow [0, 1]$中轨迹的丰富性是分析轨迹的各种状态和它们的密度。例如，移除$[0, 1]$中所有周期性的、所有周期点的轨迹，这些轨迹在有限次$f$迭代后最终成为周期性轨迹，以及那些轨迹渐近趋于周期性的点[47]。这样仍然有属于集合$U$的点，它们的轨迹不是渐近周期性的，但在回归意义上满足$f^n(U) \cap U \neq \emptyset$，其中$n > 0$。
**5.4 总结与进一步阅读建议**

在本章中，我们总结了共进化（co-evolution）这一研究领域中广泛存在的各种病理问题，这些问题已被证实会对共进化搜索性能产生负面影响。我们首先回顾了这些共进化的病理，并将它们分为三个类别。第一类涉及过度专门化（overspecialization）、聚焦（focusing）以及非传递性（intransitivity）的病理，这些病理与共进化问题中的循环结构相关，并且已经取得了深刻的理论理解。这些结构源于解决方案集合所依据的二元关系，通过对解决方案在竞争性共进化环境中的表现进行成对比较得以确定。

尽管我们讨论了如何利用有向图理论（Digraph Theory）来全面刻画共进化问题中的结构 [19]，我们也引用了其他使用序关系理论（Order Theory）来进行类似描述的研究 [8, 9]。这并不令人意外，因为连接强连通分量的非循环图（acyclic graphs）对应于有限集合的偏序关系 [2–4]。竞争性共进化中的循环结构可以描述为共进化有向图的谱（spectrum），其中一端是具有最大数量$3$-循环（$3$-cycles）的规则有向图（regular digraphs），另一端是不含任何循环的传递有向图（transitive digraphs）。在这两个极端之间，还有一大范围的不可约有向图（irreducible digraphs），它们可根据其中的$3$-循环数量进行区分；以及另一大范围的可约有向图（reducible digraphs），所有这些可约有向图都可以通过强连通分量的分解所形成的传递性支配关系链来表现。

当把竞争性共进化设置视为一个两玩家的策略博弈（strategic game），那么解决方案集合上的二元关系可以进一步表示为策略集合上的支配关系（dominance relations）。这种转换与博弈论（Game Theory）[51]相联系，并且许多过去的研究在竞争性共进化中都采用了这一方法。除了对共进化问题（有向图）是否可解（可约）或不可解（不可约）的当前定性理解之外，对这些可解的共进化问题在搜索最优解（即支配解）的速度方面也有定量刻画 [19]。
剩下的两类病态现象主要是关于协同进化过程动态结构的研究。具体而言，所考虑的基本对象是种群配置（状态），它不仅包括种群中存在的具体解，还包括这些解在种群中的数量。显而易见，此时所研究的基本对象变得更为复杂，因此该领域的理解部分依赖于受控的经验研究，这些研究采用了生物学领域的概念。例如，我们归类的最后一组病态现象包括“脱离”与“平庸的稳定状态”，这些病态现象出现在更为复杂的双种群竞争性协同进化中，其中具有非对称交互。然而，通过使用动力系统理论，已经获得了关于包括“循环”和“遗忘”病态现象的较深层次理论理解。在这里，博弈论，特别是演化博弈理论（EGT），作为建模界面起着重要作用，能够构建协同进化的动力系统，例如复制子映射（replicator maps）[35]。大多数关于“循环”和“遗忘”的观察可以表述为这些协同进化动力系统生成的周期轨迹。然而，也存在其他显著更复杂和混沌的动态行为，对其形式研究已经开始[61, 62]。需要注意的是，目前为止的这些协同进化形式研究主要涉及确定性动态。还有一些研究考虑了随机动态，这些研究直接纳入有限种群效应和进化过程中选择与变异的其他随机性，例如通过使用马尔科夫过程来研究协同进化中的种群动态[33]，以及研究不同问题结构中CEA的搜索性能[18, 19]。我们通过文献综述简要总结了几项关键研究，并集中讨论了旨在解决协同进化问题中某些属性（例如不可传递性）导致的各种病态现象的方法。这些方法直接作用于选择过程。需要注意的是，还有其他机制和启发式方法，例如利用外部档案影响选择的方式[26, 34, 56]。在此，我们重点关注那些直接改变或重新定义解的适应度值的的方法。这些方法通常被统称为“多样性维护技术”[17]，因为此前的研究声称这些方法对于协同进化搜索性能的影响来源于它们能够在种群中引入并保持多样性。我们首先介绍一个关于适应度共享（fitness sharing）的方法族。适应度共享的概念并非新创，它早期在遗传算法（GA）中被引入并用于更好地解决多峰优化问题[42]。对于竞争性CEA，适应度共享被用作一种选择过程中的机制，通过改变分配给竞争代理的适应度值，从而引入并维持一个多样化的测试实例种群。这些测试实例将共同为种群中的解提出不同的挑战需要解决。因此，适应度共享的使用可能会防止竞争性协同进化中过度专业化和集中化的出现。在接下来的内容中，我们将介绍并简要讨论两种不同的适应度共享在协同进化中的实现方式。
第一个方法涉及竞争适应度共享（competitive fitness sharing），用于捕食者-猎物的竞争共进化环境中，这是在文献[56]中提出的。虽然这一实现可以轻松应用于单一竞争共进化环境，但以下描述中的符号表现为文献[56]中提出的双种群（two-population）竞争共进化环境。每个演化的解代理$s \in S$与一个对手或测试代理的随机子集$T'_s \subset T$进行交互（如果考虑锦标赛环境，而非涉及整个种群$T$的循环赛环境）。与通常的交互模型（即获胜-失败的游戏结果为获胜者分配相同的胜利值）不同，竞争适应度共享允许根据对手分配不同的值。这是通过为对手$t \in T'_s$的胜利来分配一个特定权重$w_t = 1/n_t$实现的，其中$n_t$表示整个种群中击败对手$t$的解决方案代理数量。假设该种群的大小为$|S|$。如果只考虑获胜-失败的游戏结果，并且平局（包括自我对弈）对双方而言均算作失败，那么权重范围为$w_t \in [1/(|S|-1), 1]$。对于解代理$s$的共享适应度计算公式为：

$$f_{\mathrm{sh}}(s) = \sum_{t \in T'_s} w_t \cdot M(s, t)$$  

其中，$M(s, t)$是常规的赋值函数，若$s$击败$t$，则$M(s, t) = 1$，否则为$M(s, t) = 0$。本质上，这种方法鼓励选择独特的解决方案，即能够解决一些测试案例，而这些案例其他竞争解决方案难以解决。值$1/n_t$与$n_t$呈反比例关系，即当$n_t$较小时，$1/n_t$值较大。

隐式适应度共享（Implicit fitness sharing）[24]采用一种算法程序，该程序能够近似计算对于解代理$s \in S$的共享适应度值[59]，因此具有无需设计特定于问题的共享函数的优势[23]。我们将在单种群竞争共进化环境中介绍这种方法。设$g(s, t)$为策略$s$在与$t$的双人游戏中获得的游戏收益（类似地，$g(t, s)$是策略$t$获得的游戏收益）。对于种群中的每个代理$s \in S$，从种群中抽取随机的测试对手子集$T'_{c}, \ c = 1,\ldots,C$，其中$T'_c \subset S$。注意每个子集$T'_c$都排除$s$以避免自我对弈。对于每个随机子集$T'_c$，程序将游戏收益分配给在双人游戏中对$s$获得最大游戏收益的对手，这被定义为：

$$\max_{t \in T'_c} g(t, s)$$  

在出现平局时，对手会均分适应度（游戏收益除以出现平局的对手数量）。注意，这是用于非零和游戏的交互模型，解决方案概念或优化的目标是最大化累积收益。对于零和游戏，可以采用最大的获胜边际，即：

$$\max_{t \in T'_c} \delta_{\mathrm{win}}(t, s)$$  

其中，$\delta_{\mathrm{win}}(t, s) = g(t, s) - g(s, t)$，如果$t$击败$s$，则计算获胜边际；否则为$0$。

另一种直接改变适应度值以影响代理选择的方法是减少毒性（Reducing virulence）[11,12]。这种基于自然启发的方法基于以下观察：减少毒性可能解决生物共进化系统中的脱离问题。高毒性的寄生物会削弱宿主，从而降低宿主传播到新宿主的繁殖机会。相比之下，中度毒性的菌株能够更快繁殖，因此在演化上占优势。对于人工竞争共进化系统来说，采用降低毒性的方法可以避免因种群中高性能代理对测试案例的快速增殖而导致的脱离问题，进而促进更长期和持续的共进化搜索。

假设代理$s$的适应度值$f(s)$与其对解决测试案例的表现成正比，例如：

$$f(s) = \sum_{t \in T'} M(s, t)$$  

减少毒性的实现通过重新调整$f(s)$，定义为：

$$f(s, v) = 2f(s)/v - (f(s)/v)^2$$
$$f({\boldsymbol s}, v) = \frac{2f({\boldsymbol s})}{v} - \left(\frac{f({\boldsymbol s})}{v}\right)^2 \ [12]$$，其中参数 $v \in [0.5, 1.0]$ 控制代表倾向于最大化 $f({\boldsymbol s})$ 的病毒性水平。例如，标准的协演化在最大病毒性水平下运行，即 $v = 1.0$。定义域为 $f({\boldsymbol s})$，值域为 $f({\boldsymbol s}, 1.0)$ 的函数是单射（即一对一映射），因此能够鼓励个体最大化其适应度值。在 $v = 0.75$ 的情况下，函数是满射（即多对一映射），这意味着一些高值的 $f({\boldsymbol s})$ 被映射到较低的 $f({\boldsymbol s}, 0.75)$ 值，从而对种群中高表现个体进行惩罚，避免搜索集中于最大化相对适应度值的解方案。最后，我们引入了基于多目标进化算法框架的帕累托协演化方法 [27, 32, 52]。 

在竞争协演化中，基于简单加权求和形式的相对适应度值，例如 $$f({\boldsymbol s}) = \sum_{{\boldsymbol t} \in {\boldsymbol T}'} M({\boldsymbol s}, {\boldsymbol t})$$ 的使用可能受到质疑。这是因为针对某些测试案例的结果可能更加重要，因而需要非均匀权重。更重要的是，由协演化生成的不同种群配置可能会改变权重的集合。帕累托协演化直接解决了这一问题，它采用测试案例来表征问题的不同潜在目标，当中解决方案个体必须应对这些目标。选择过程基于种群中竞争个体相对于一组测试案例的帕累托支配关系。

对于一个由解决方案 $ {\boldsymbol s}_i, \ i = 1,\ldots,S$ 和测试案例 ${\boldsymbol t}_k, \ k = 1,\ldots,T$ 组成的两种群竞争协演化问题来说 [27]，帕累托支配关系可以表示为：
$$
{\boldsymbol s}_i \succ {\boldsymbol s}_j \iff \forall {\boldsymbol t}_k : g({\boldsymbol s}_i,{\boldsymbol t}_k) \geq g({\boldsymbol s}_j,{\boldsymbol t}_k) \land \exists {\boldsymbol t}_k : g({\boldsymbol s}_i,{\boldsymbol t}_k) > g({\boldsymbol s}_j,{\boldsymbol t}_k)。
$$
通俗地讲，当且仅当 ${\boldsymbol s}_i$ 在 $T$ 个目标中的每一个（每个由测试案例 ${\boldsymbol t}_k$ 表示）上至少具有与 ${\boldsymbol s}_j$ 相同的游戏收益，并且至少存在一个目标上 ${\boldsymbol s}_i$ 的游戏收益优于 ${\boldsymbol s}_j$ 时，${\boldsymbol s}_i$ 帕累托支配 ${\boldsymbol s}_j$。

如果 ${\boldsymbol s}_i$ 在至少一个目标上优于 ${\boldsymbol s}_j$，且 ${\boldsymbol s}_j$ 在至少另一个目标上优于 ${\boldsymbol s}_i$，并且在其余目标上两者是等价的，则两者互为不支配关系。一般而言，该机制会试图根据从帕累托支配关系中获得的帕累托层对个体进行排序。属于帕累托前沿的互为不支配个体是最优先选择的。可以通过逐步移除处于较高帕累托层中的互为不支配个体来系统性地获得后续帕累托层 [32]。

需注意的是，在单种群竞争环境中，此类帕累托支配排序将生成一个强连通组件的无环支配链，这些组件对应于帕累托层 [19]。
