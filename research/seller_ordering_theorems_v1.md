# 有限 K 件商品的卖家排序：强起始记忆下的均衡转移与收入定理

日期：2026-09-20。研究推导稿；不是已完成新颖性审查或同行评审的论文。

本稿接续 seller_ordering_general_n_K.md。优化目标始终是 seller 总预期收入。新的关键是：可以利用不同顺序之间的博弈等价关系比较收入，不必先为每个顺序分别求出全部 bidder 出价的闭式。

## 0. 结果状态

| 结果 | 条件 | 状态 |
|---|---|---|
| 第一轮赢家在所有后续历史下都保持估值优势 | 强起始记忆条件 | 已证明 |
| 后续逐轮当前价值出价是一组序贯最优响应 | 同上，无预算等限制 | 已证明；未称其为任意环境中的占优策略 |
| 多轮游戏归约到与后续顺序无关的辅助拍卖 | 一般有界私人估值，允许各商品独立 | 已证明 |
| 同一辅助均衡可转移到所有起始商品在前的顺序 | 辅助均衡存在，统一选择后续真实出价 | 已证明条件定理 |
| R(pi)=C+alpha_W W(pi)，从而比较顺序无需知道 C | 同上 | 已证明条件定理 |
| 辅助均衡存在性 | 有限类型和有限报价网格版本 | 由有限 Bayesian 游戏的混合均衡存在性得到 |
| 任意 n,K 的显式均衡与收入 | 单参数私人值 v_ik=a_k+lambda_k theta_i | 已证明；同一人的商品估值相关 |
| 强起始商品排在第一位严格优于放在后面 | 进一步加强起始边权，相关条件见第 6 节 | 已证明比较界 |
| NP-hardness、1/2 近似、子集 DP | 上述收入定理适用的均衡选择/离散版本 | 已证明归约与保证 |
| 第一轮少收、总收入更高 | 单参数显式五商品例子 | 已证明并核验 |

未解决：连续多维私人类型辅助拍卖的一般均衡存在性；没有强起始商品的一般图；弱记忆下多次更换赢家；任意均衡选择下统一收入比较。

## 1. 主模型：仍然是 n 个 bidder、K 件商品

n>=2、K>=2 均为任意有限整数。每件商品只拍一次，所有 bidder 参加全部拍卖，可以获得多件商品；没有 entry cost、预算或单位需求约束。每轮第二价格拍卖，非负出价，平局规则预先固定。效用按商品取得时价值减支付加总。

先取 rho=1、alpha_W>alpha_L>=0，记 delta=alpha_W-alpha_L>0。基础估值 x_ik in [0,1]，每个人事先知道自己的向量，其他人不知道。下面的均衡转移论证允许 x_ik 在 bidder 和商品之间完全独立，例如全部独立 U[0,1]；也允许其他已给定先验。

从 K 件商品中指定一件 s，作为候选起始商品；没有额外增加一件商品。其余集合记 J，|J|=K-1。M 非负。记普通商品之间的总边权

\[
G_{tot}=\sum_{j,k\in J:\,j\ne k}M_{jk},
\qquad B=\min_{k\in J}M_{sk}.
\]

本稿首先研究满足

\[
\boxed{\delta B>1+\delta G_{tot}}
\tag{A}
\]

的图。它表示：赢得 s 对后续商品提供的差异记忆，足以覆盖初始估值差距及其余商品可能给竞争者积累的差异记忆。

这是一个强条件，用于构造可解理论基准。它不是任意记忆图自动满足的性质，不应被写成原始一般模型已经解决。

## 2. 引理：起始赢家的估值优势在所有后续历史中保持

设 s 第一轮出售，赢家为 w。考虑任意之后的赢家历史，不限于均衡路径。对每件尚未出售的 k in J，以及 i!=w，

\[
\begin{aligned}
v_{wk}^t-v_{ik}^t
={}&x_{wk}-x_{ik}+\delta M_{sk}\\
&+\delta\sum_{r=2}^{t-1}
\bigl(\mathbf1\{w_r=w\}-\mathbf1\{w_r=i\}\bigr)M_{\pi_rk}\\
\ge{}&-1+\delta B-\delta G_{tot}>0.
\end{aligned}
\tag{1}
\]

共同的 alpha_L 位移在差值里抵消。因此第一轮赢家即使在一些后续轮次假想地让其他人获胜，仍然对每件剩余商品有最高当前估值。

## 3. 引理：后续真实出价策略在每个历史下都序贯理性

给定第一轮赢家 w，规定从第二轮起每个人报实际当前估值。由 (1)，w 赢得所有剩余商品。

这个策略组合可以通过反向归纳验证：

- 对 i!=w，如果本轮不赢，此后在该策略组合下继续不赢，净效用为零。如果偏离并赢本轮，需要支付至少 w 的当前出价，高于自己的当前价值；之后仍无法改变 (1) 的估值优势，未来净效用仍为零。因此偏离获胜不利。
- 对 w，按当前价值出价已能获胜，并获得严格正的当期净收益。保持获胜的其他出价不改变当前支付。若主动输掉本轮，相对于获胜，w 对每个未来商品的估值弱下降，获胜对手的未来出价弱上升；其他人的对应变化不使 w 获益。未来继续由 w 获胜，其未来净效用不增加，再加上失去当前正收益，所以主动输不利。

有限轮反向归纳排除跨多轮的有利偏离。本论证对所有类型实现、所有后续历史都成立，因而不要求特定的后验信念。

边界：这里证明的是指定策略组合下的最佳响应。没有声称真实出价是任意其他出价策略下的动态占优策略，也没有证明所有第二价格均衡都有同样收入。

## 4. 核心归约：后续顺序只改变第一轮的一个共同常数

对 pi=(s,pi_2,...,pi_K)，定义

\[
W(\pi)=\sum_{r<t}M_{\pi_r\pi_t},
\quad S_i=\sum_{k=1}^Kx_{ik},
\quad L_i(x_{-i})=\sum_{k\in J}\max_{j\ne i}x_{jk}.
\]

S_i 是 i 的全部商品基础估值之和；L_i 是在 i 获得所有后续商品时，后续商品的基础竞争价格之和。不同商品的最大值可以来自不同对手，不能把它误写成某个对手的总估值。

若 i 赢第一轮，第一轮支付 p，则后续他赢得全部商品。他的总毛价值为

\[
S_i+\alpha_W W(\pi),
\]

后续支付之和为

\[
L_i(x_{-i})+\alpha_L W(\pi).
\]

因此赢第一轮的总净效用为

\[
\boxed{U_i^{win}=S_i-L_i(x_{-i})+\delta W(\pi)-p.}
\tag{2}
\]

输第一轮时净效用为零。定义

\[
\phi_i(x)=S_i-L_i(x_{-i}).
\]

phi_i 不依赖后续商品的顺序，但依赖对手的私人值；所以辅助拍卖是相互依赖价值的拍卖，不能把 phi_i 当成 bidder 已知的私人价值并直接报真。

### 4.1 与顺序无关的辅助拍卖

构造一个数学辅助游戏：每人只观察原本自己的估值向量，报 q_i；最高者获胜、支付第二高 q，赢者效用 phi_i-q_(2)，输者为零。

由于 x_ik in [0,1]，

\[
-(K-1)\le\phi_i\le K.
\]

可把辅助出价限制在 [-(K-1),K]：低于全局价值下界的出价被该下界弱支配；高于全局价值上界的出价被该上界弱支配。这里用的是赢者 phi_i-q、输者零的辅助游戏效用，需逐对手类型实现比较，论证不要求 phi_i 被本人观察到。

辅助出价和支付可以为负，但这只是数学变换，并不向真实拍卖引入补贴。

### 4.2 均衡转移定理（连续版本有存在性前提）

假设辅助拍卖存在 Bayesian Nash 均衡 tau，允许混合策略。对每个以 s 开头的排列，在第一轮使用

\[
\boxed{b_i^\pi=q_i^\tau+\delta W(\pi),}
\tag{3}
\]

后续使用第 3 节的真实出价策略，则得到原多轮拍卖的一组序贯理性策略。

证明：所有出价加同一个常数，赢家不变，第一轮支付相应增加该常数。由 (2)，每个 bidder 在任何类型/随机化实现上的净效用恰好等于辅助游戏的净效用。

合法性：W(pi)>=(K-1)B，所以最小实际出价至少

\[
\delta(K-1)B-(K-1)>0.
\]

所有实际支付仍非负。原拍卖允许的任何第一轮偏离，减去 delta W 后对应辅助出价；超出辅助区间的偏离被区间端点弱支配。因此没有遗漏低价或高价偏离。后续序贯理性由第 3 节保证；信念可按 Bayes 规则在路径上更新，后续最优响应不依赖信念的选择。

这不是对一般连续辅助拍卖存在性的证明。连续多维情形仍须单独证明辅助均衡存在，或者明确作为条件定理使用。

## 5. 卖家收入定理：不必求出 C，也可以优化后续顺序

设辅助均衡 tau 的赢家为 w，第二高出价为 q_(2)。第一轮支付为 q_(2)+delta W，后续支付为 L_w+alpha_L W。因此逐条耦合路径上都有

\[
\boxed{\sum_tp_t=q_{(2)}+L_w+\alpha_WW(\pi).}
\tag{4}
\]

令

\[
C_s^\tau=\mathbb E_\tau[q_{(2)}+L_w],
\]

则

\[
\boxed{R^\tau(\pi)=C_s^\tau+\alpha_WW(\pi).}
\tag{5}
\]

C_s^tau 可能难算，但它不依赖后续顺序。这正是卖家优化所需要的性质。对任意两个以 s 开头的顺序，

\[
\boxed{R^\tau(\pi)-R^\tau(\pi')
=\alpha_W[W(\pi)-W(\pi')].}
\tag{6}
\]

相比上一份笔记中尚未知的 Q(pi)，现在在这个结构类中得到 Q(pi)=C_s^tau+delta W(pi)。战略竞价使赢家专属记忆的那部分作用也进入卖家收入。

### 均衡选择的范围

必须在不同顺序下使用同一个辅助均衡 tau，或统一使用该辅助均衡集合上的乐观/悲观收入标准。在“后续真实出价、第一轮辅助均衡”这一均衡族内，所有顺序的收入集合只是同一个集合平移 alpha_W W(pi)，因此其 supremum/infimum 也按相同公式比较（集合非空时）。

没有声称任意挑选不同均衡后仍保持同样的收入排序，也没有把本均衡族等同于所有可能的多轮第二价格均衡。

## 6. 起始商品是否应该真的放第一位？

第 5 节先解决了 s 已放第一位后的排列问题。下面给出足够条件，保证把它放后面更差。

为明确构造，取 alpha_L=1、alpha_W=2、delta=1，所有 x_ik in [0,1]。设 M_sk=B（k in J）、M_js=0（j in J），J 内任意非负矩阵 G，总权重 G_tot。选择

\[
\boxed{B=G_{tot}+K+1.}
\tag{7}
\]

它满足强起始条件 (A)。辅助出价下界与 L_w>=0 给出 C_s^tau>=-(K-1)，所以任何本均衡族中以 s 开头的顺序满足

\[
R\ge2(K-1)B-(K-1).
\tag{8}
\]

若 s 不在第一位，至少损失一条权重 B 的起始前向边。由于其余 G 的前向贡献至多 G_tot，

\[
W(\pi)\le(K-2)B+G_{tot}.
\]

在任意有定义的均衡中，bidder 都能选择全程零出价，保证自己不支付正价格、总效用非负。故每个人的 ex ante 均衡效用非负，seller 收入不超过期望总毛价值。每件基础值至多 1，每个记忆项至多按 alpha_W=2 计入，因此

\[
R(\pi)\le K+2[(K-2)B+G_{tot}].
\tag{9}
\]

(8) 减 (9) 等于 2(B-G_tot)-(2K-1)=3>0。因此，在存在本构造均衡的前提下，s 在前的一组可实现收入严格高于所有 s 不在前的均衡收入上界。

这一结论证明了一个明确类中的 seller 最优首件选择；不是把“强起始商品优先”未经证明地当成启发式。

## 7. 排序复杂性与算法

s 第一位后，它的所有出边贡献是常数。seller 的剩余问题为

\[
\max_{\sigma\in\mathfrak S_J}\sum_{j\prec_\sigma k}G_{jk}.
\tag{10}
\]

任意非负整数加权有向图都可嵌入 G，且 (7) 的 B 具有多项式表示长度。因此，在第 5 节指定均衡族/共同辅助均衡选择下，准确优化收入至少和最大加权前向边排序一样困难。辅助基线 C 无需计算，因为它对所有后续顺序相同。

离散版本见第 8 节，其辅助均衡存在性可以补齐；连续多维版本保留前述存在性前提。另有第 9 节的连续单参数显式子类，也直接提供非空的均衡构造。

### 1/2 近似

取任意后续顺序 sigma 和它的逆序 sigma_rev，二者前向边总权重之和等于 G_tot。选较好者，则

\[
W_G(\widehat\sigma)\ge\tfrac12G_{tot}\ge\tfrac12W_G(\sigma^*).
\]

因此对可优化的顺序收入项 alpha_W W_G 提供 1/2 保证。这里特意把固定的强起始收益剥离，避免大常数掩盖算法质量。

### 精确与可解特例

\[
D(\varnothing)=0,\qquad
D(S)=\max_{k\in S}\left[D(S\setminus\{k\})+
\sum_{j\in S\setminus\{k\}}G_{jk}\right].
\]

预计算子集入边和后为 O(K 2^K) 时间。每对商品较优方向形成 DAG 时，拓扑排序最优；无向支撑为森林时可以逐边选较优方向；对称边权则顺序项不变。这些结论只在已获得 (5) 的结构类中转化为 seller 收入保证。

## 8. 保留独立私人估值的有限网格版本：补齐均衡存在性

这是一项明确的技术变体，不假装等同于原连续报价机制：

- 每个基础估值取有限网格，例如 iid {0,1} 或 {0,1/m,...,1}；不同 bidder、不同商品相互独立。
- 出价允许的单位为 Delta=1/m，记忆参数及实际价值位移都在此网格上。
- 可选一个足够大的有限报价上限，覆盖本构造的第一轮区间和后续真实值；这不是预算约束，但确实是报价空间的技术限定。
- 取 alpha_L=1、alpha_W=2、整数 G、B 如 (7)，即可保证所有平移都是报价网格的整数倍。

辅助拍卖于是是有限类型、有限动作的 Bayesian 游戏，其有限战略式表示存在混合 Nash 均衡。因此第 4 节的 tau 非空，转移构造得到原多轮游戏中的序贯理性策略。

后续价值落在网格上，所以第 3 节真实出价策略合法。第一轮超出辅助区间的网格偏离被端点弱支配。于是第 5--7 节在这个离散私人值版本中不再以未证明的辅助均衡存在性为前提。

这给出一个任意 n、任意有限 K、非退化独立私人估值、赢家/输家不对称记忆下的完整理论基准。它仍然依赖强起始记忆，不覆盖一般 M。

## 9. 连续单参数子类：把策略与收入都显式写出来

为了提供一个完全可读的闭式例子，另设

\[
v^0_{ik}=a_k+\lambda_k\theta_i,
\quad a_k\ge0,\quad\lambda_k>0,\quad
\theta_i\overset{iid}{\sim}F\text{ on }[0,1].
\tag{11}
\]

每人只有一个私人参数 theta_i，但商品的 a_k、lambda_k 可不同。同一人的商品估值相关，不能把本节误当成各商品独立的模型。

记 A=sum a_k、Lambda=sum lambda_k。把强起始条件改成 delta B>max_k lambda_k+delta G_tot，后续优势与真实出价证明不变。

给定排列 pi，以 s 开头，第一轮存在显式均衡

\[
\boxed{b_{is}=a_s+\lambda_s\theta_i+\delta W(\pi).}
\tag{12}
\]

证明：对手最高类型为 T=max_{j ne i}theta_j，其最高出价为 a_s+lambda_s T+delta W。i 赢第一轮后，毛价值 A+Lambda theta_i+alpha_W W；后续支付为 sum_(k!=s)a_k+(Lambda-lambda_s)T+alpha_L W。减去第一轮支付后得到

\[
\boxed{U_i^{win}=\Lambda(\theta_i-T),\qquad U_i^{lose}=0.}
\tag{13}
\]

因此当且仅当 theta_i>T 时应获胜，(12) 实现这个门槛，对所有对手类型实现都成立。获胜者是最高 theta，随后获得所有商品。

总收入逐类型实现为

\[
\boxed{\sum_tp_t=A+\Lambda\theta_{(n-1:n)}+\alpha_WW(\pi).}
\tag{14}
\]

若 F=U[0,1]，则

\[
\boxed{R(\pi)=A+\Lambda\frac{n-1}{n+1}+\alpha_WW(\pi).}
\tag{15}
\]

在这个子类中，改变后续顺序的记忆收益完全进入 seller 收入，bidder 的净收益不变：赢家获得 Lambda(theta_max-theta_second)，其他人零。

注意：这也揭示子类的局限。各商品的基础排序本来就由同一个 theta 决定，强记忆又锁定了第一轮赢家，所以不能用它证明一般模型中记忆会怎样使赢家反复更替。

## 10. 一个精确的“第一轮少收、总收入更高”例子

取 n=3、K=5，商品为 s,H,B,C,D；lambda_k=1，theta_i iid U[0,1]，只有 a_H=100，其他 a_k=0。alpha_L=1、alpha_W=2。

普通商品图只有 B->C 权重 3、C->D 权重 2、D->B 权重 1；H 没有出边，也没有来自其他普通商品的入边。G_tot=6。令 s 到 H,B,C,D 的每条边权重为 B_seed=12，所有指向 s 的边为零。

比较两种事前承诺顺序：

- pi=(s,H,B,C,D)，W(pi)=4*12+3+2=53。
- pi'=(H,s,B,C,D)，W(pi')=3*12+3+2=41。

对 pi，使用上文均衡。对 pi'，H 的胜负不改变后续记忆，之后的起始拍卖有对所有类型实现最优的显式策略，因此第一轮 H 的真实出价最优，其信号不改变后续指定策略的最优性；之后套用剩余四商品版本。

由于 E theta_(2:3)=1/2，精确预期支付为：

| 顺序 | 第一轮收入 | 各轮预期支付 | 总收入 |
|---|---:|---|---:|
| s,H,B,C,D | 53.5 | 53.5, 112.5, 12.5, 15.5, 14.5 | 208.5 |
| H,s,B,C,D | 100.5 | 100.5, 41.5, 12.5, 15.5, 14.5 | 184.5 |

因此先卖 s：第一轮少收 47，但总收入多收 24。H 本身是当前价格更高的商品，但先经历 s 会增加从 H 获得的后续收入；赢家专属记忆的一部分还通过第一轮战略溢价被提前支付。

两种顺序均为均衡支付比较，不是把买家静态报真直接用于第一轮起始拍卖。例子的基础值结构属于第 9 节相关单参数子类。

## 11. 衰减扩展：可保留有限 K，但最优起始位置需另证

若 0<rho<=1，把强起始条件加强为

\[
\delta\rho^{K-2}B>1+\delta G_{tot}
\]

（单参数子类把 1 换成 max lambda_k），则起始赢家在所有后续历史下仍然领先。

定义

\[
W_\rho(\pi)=\sum_{r<t}\rho^{t-r-1}M_{\pi_r\pi_t}.
\]

第 4--5 节的归约与均衡转移仍成立，只需用 delta W_rho 平移第一轮出价，得到

\[
R^\tau(\pi)=C_s^\tau+\alpha_WW_\rho(\pi).
\]

这说明路径衰减可以进入模型，而不需要立即放弃收入结构。但是目标此时依赖距离，普通线性排序的算法不能无条件照搬；第 6 节关于 s 最先最优的比较界也需重新证明。rho 很小时维持同一强优势所需 B 可能极大，因此不应把本条件当作经验上普遍合理的假设。

## 12. 论文价值与下一步缺口

这次真正建立的桥梁是：特定记忆结构使多个销售顺序的竞价博弈通过一个共同出价平移相连，卖家可以在不知道每个 bidder 出价闭式、甚至不知道收入常数 C 的情况下，求出最优的后续顺序。

但这不等于完整论文已经完成：

1. 强起始记忆锁定后续赢家，是本定理的关键，也可能是审稿人认为过强的地方。最重要的推广是减弱它，让赢家会变化，同时尽量保留收入比较或近似界。
2. 原连续多维独立私人值版本的辅助均衡存在性仍待证明。不能拿一般私人价值第二价格拍卖的真实性来代替，因为辅助 phi_i 依赖对手信息。
3. 本稿收入比较使用统一的辅助均衡和后续真实出价族。若要覆盖所有序贯均衡，需更强的均衡选择/稳健结果。
4. 图排序困难性与基础 1/2 近似本身是已知算法知识；潜在贡献在于从战略记忆拍卖到这些目标的严格桥梁，而不是重新提出已有排序算法。
5. 必须与既有序贯互补估值、动态外部性、排序拍卖文献进行完整新颖性比对，才能判断投稿价值。

来源边界：最大无环子图的困难性可参见 [Austrin, Manokaran, Wenner, On the NP-Hardness of Approximating Ordering Constraint Satisfaction Problems](https://arxiv.org/abs/1307.5090)。关于不完全信息、不连续拍卖的均衡存在不能自动断言，可参见 [Jackson, Simon, Swinkels, Zame, Communication and Equilibrium in Discontinuous Games of Incomplete Information](https://authors.library.caltech.edu/records/pqb72-7qr31)。这些文献不等于对本稿定理新颖性的认证。

## 13. 核验记录

- 单参数例子：遍历 24 个后续排列，并检查 2880 个包括偏离路径的历史状态，验证起始赢家始终领先、后续当前价值出价的偏离不利、第一轮净效用等于 Lambda(theta_i-T)。
- 精确有理数计算验证总收入 208.5 与 184.5，以及第一轮支付 53.5 与 100.5。
- 一般多维类型：随机生成 30 组各商品分别取值的类型矩阵及任意辅助出价，对全部 24 个后续排列验证 720 个逐路径支付/效用恒等式。这里的任意辅助出价不被称为均衡；计算仅验证转移的代数恒等式。
- 子集 DP 与全排列最优普通图权重一致，均为 5。
- 数值检查不替代上述证明，也不证明连续辅助均衡存在。

下附核验代码，使用 Python 标准库和精确有理数。

```python
from fractions import Fraction as F
from itertools import permutations
import json
from random import Random

K, n = 5, 3
alpha_l, alpha_w = F(1), F(2)
delta = alpha_w-alpha_l
theta = [F(1,10), F(2,5), F(9,10)]
a = [F(0),F(100),F(0),F(0),F(0)]
lam = [F(1)]*K
M = [[F(0) for _ in range(K)] for _ in range(K)]
M[2][3],M[3][4],M[4][2]=F(3),F(2),F(1)
ordinary_total=sum(sum(row) for row in M)
B=ordinary_total+K+1
for k in range(1,K):M[0][k]=B

def W(pi):
    return sum(M[pi[s]][pi[t]] for s in range(K) for t in range(s+1,K))

def update(mu,j,w):
    return [[mu[i][k]+(alpha_w if i==w else alpha_l)*M[j][k] for k in range(K)] for i in range(n)]

def values(mu,j):
    return [a[j]+lam[j]*theta[i]+mu[i][j] for i in range(n)]

def truthful_tail(pi,t,mu):
    util=[F(0)]*n
    prices=[]
    winners=[]
    for j in pi[t:]:
        v=values(mu,j)
        w=max(range(n),key=lambda i:v[i])
        price=sorted(v)[-2]
        util[w]+=v[w]-price
        prices.append(price);winners.append(w)
        mu=update(mu,j,w)
    return util,prices,winners

state_checks=0
def verify_all_histories(pi,t,mu,leader):
    global state_checks
    if t==K:return
    j=pi[t];v=values(mu,j)
    assert all(v[leader]>v[i] for i in range(n) if i!=leader)
    util,_,winners=truthful_tail(pi,t,mu)
    assert all(w==leader for w in winners)
    for deviator in range(n):
        if deviator!=leader:
            after,_,_=truthful_tail(pi,t+1,update(mu,j,deviator))
            dev_util=v[deviator]-v[leader]+after[deviator]
            assert dev_util<util[deviator]
        else:
            for other in range(n):
                if other==leader:continue
                after,_,_=truthful_tail(pi,t+1,update(mu,j,other))
                assert after[leader]<util[leader]
    state_checks+=1
    for w in range(n):verify_all_histories(pi,t+1,update(mu,j,w),leader)

all_scores=[]
for suffix in permutations(range(1,K)):
    pi=(0,)+suffix
    total=W(pi)
    for leader in range(n):
        mu=[[F(0)]*K for _ in range(n)]
        mu=update(mu,0,leader)
        verify_all_histories(pi,1,mu,leader)
        tail_u,tail_p,tail_w=truthful_tail(pi,1,mu)
        T=max(theta[j] for j in range(n) if j!=leader)
        first_price=a[0]+lam[0]*T+delta*total
        winning_u=a[0]+lam[0]*theta[leader]-first_price+tail_u[leader]
        assert winning_u==sum(lam)*(theta[leader]-T)
    leader=max(range(n),key=lambda i:theta[i])
    T=sorted(theta)[-2]
    _,tail_p,_=truthful_tail(pi,1,update([[F(0)]*K for _ in range(n)],0,leader))
    first_price=a[0]+lam[0]*T+delta*total
    assert first_price+sum(tail_p)==sum(a)+sum(lam)*T+alpha_w*total
    all_scores.append((pi,total))

pi=(0,1,2,3,4)
rho_theta=F(n-1,n+1)
revenue=sum(a)+sum(lam)*rho_theta+alpha_w*W(pi)
other_revenue=revenue-alpha_w*B
first_anchor=a[0]+rho_theta+delta*W(pi)
first_terminal=a[1]+rho_theta
assert revenue==F(417,2)
assert other_revenue==F(369,2)

dp={0:F(0)}
items=list(range(1,K))
for mask in range(1,1<<len(items)):
    dp[mask]=max(dp[mask^(1<<q)]+sum(M[items[p]][items[q]] for p in range(len(items)) if p!=q and mask>>p&1) for q in range(len(items)) if mask>>q&1)
assert dp[(1<<len(items))-1]==max(score for _,score in all_scores)-(K-1)*B
best=max(all_scores,key=lambda z:z[1])
nonanchor_upper=sum(a)+sum(lam)+alpha_w*((K-2)*B+ordinary_total)
anchor_lower=sum(a)+sum(lam)*rho_theta+alpha_w*(K-1)*B
assert anchor_lower>nonanchor_upper

# This checks a pathwise identity for independent multidimensional types;
# arbitrary residual bids here are NOT claimed to constitute an equilibrium.
rng=Random(20260920)
general_checks=0
for _ in range(30):
    base=[[F(rng.randrange(101),100) for _ in range(K)] for _ in range(n)]
    residual=[F(rng.randrange(-400,501),100) for _ in range(n)]
    w=max(range(n),key=lambda i:residual[i])
    q_second=sorted(residual)[-2]
    competitor_cost=sum(max(base[i][k] for i in range(n) if i!=w) for k in range(1,K))
    phi=sum(base[w])-competitor_cost
    constant=q_second+competitor_cost
    for suffix in permutations(range(1,K)):
        order=(0,)+suffix
        first=q_second+delta*W(order)
        memory=update([[F(0)]*K for _ in range(n)],0,w)
        total_price=first
        utility=base[w][0]-first
        for j in suffix:
            v=[base[i][j]+memory[i][j] for i in range(n)]
            assert max(range(n),key=lambda i:v[i])==w
            price=sorted(v)[-2]
            total_price+=price
            utility+=v[w]-price
            memory=update(memory,j,w)
        assert total_price==constant+alpha_w*W(order)
        assert utility==phi-q_second
        general_checks+=1

print(json.dumps({
    'off_path_states_checked':state_checks,
    'suffix_permutations_checked':len(all_scores),
    'cue_strength':str(B),
    'best_suffix_gain':str(dp[(1<<len(items))-1]),
    'best_order':best[0],
    'best_W':str(best[1]),
    'anchor_first_expected_revenue':str(revenue),
    'terminal_first_expected_revenue':str(other_revenue),
    'anchor_first_round_expected_payment':str(first_anchor),
    'terminal_first_round_expected_payment':str(first_terminal),
    'anchor_any_order_revenue_lower_bound':str(anchor_lower),
    'nonanchor_any_equilibrium_revenue_upper_bound':str(nonanchor_upper),
    'independent_type_pathwise_identities_checked':general_checks,
},indent=2))
```
