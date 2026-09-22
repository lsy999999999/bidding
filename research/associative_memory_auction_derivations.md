# 联想记忆序贯拍卖：模型审计、均衡、排序与算法

日期：2026-09-19。性质：研究工作笔记；以下为本次推导及其证明，不代表已完成同行评审或新颖性认证。

## 0. 已取得的结果与适用边界

| 结果 | 信息/参数条件 | 状态 |
|---|---|---|
| 线性记忆展开式与充分状态 | 一般 K、固定顺序 | 已证明 |
| 固定顺序下的静态 bundle 表示 | 一般 rho、取得时价值计入效用 | 已证明 |
| 跨顺序存在统一 bundle valuation 的充要条件 | rho=1、所有分配历史 | 已证明 |
| 两商品一般私人信息的最优响应积分方程 | 一般 n、连续对手最高出价 | 已推导；不是一般均衡存在性定理 |
| 两商品、两人、公开估值的动态临界出价 | 临界出价非负 | 已证明 |
| 两商品、任意 n、提前知道自身全部估值的均衡与收入 | X,Y 独立 U[0,1]；d>=1 | 已证明 |
| 强记忆策略向阈值下方的近似均衡延伸 | 同上，0<=d<=1 | 已证明偏离收益不超过 1-d |
| 两商品、任意 n、未来估值延后揭示的均衡与收入 | iid 连续分布；c,d>=0 | 已证明；属于另一信息结构 |
| 信息揭示对收入影响的人数反转 | U[0,1]、d>=1；n=2 与 n=4 | 精确计算并证明 |
| 一般 K 的加权排序等价、NP-hardness | 公共位移边界，或公开异质强优势环境 | 已证明；未证明原始 iid 非退化非对称模型的同一等价 |
| 排序的 1/2 近似、DAG/森林特例、子集 DP | 非负加权排序子问题 | 已证明 |
| 权重误差、记忆衰减的加性鲁棒界 | 上述可分解环境 | 已证明 |
| 原始一般 K 模型的近似均衡 | 非负小差异记忆；固定顺序 | 已证明偏离收益界；不是精确均衡收入保证 |

原始“全部私人估值提前已知、任意强度记忆、一般 K”的完整均衡和最优排序仍未解决。不能把下述不同信息结构/特例的定理拼成一个覆盖原模型的定理。

## 1. 冻结效用与信息结构

有 n>=2 个风险中性 bidder、K 件商品，每件只拍卖一次。无参与费、预算或单位需求约束；每个人参加每轮，可以获得多件商品。seller 承诺一个排列 pi，每轮非负出价、零保留价、第二价格拍卖，平局规则预先固定。

默认 bidder i 开始就知道自己的完整基础估值向量 theta_i=(v^0_i1,...,v^0_iK)，不知道其他人的向量。不同 bidder 和商品的基础估值相互独立，商品 k 的分布为 F_k。商品、赢家身份公开，支付可公开；在两商品结果中，末轮真实出价对后验信念不敏感。

M 是公开固定的有向矩阵。基准非负结果取 M_jk>=0、alpha_W>=alpha_L>=0。定义

\[
\delta=\alpha_W-\alpha_L,\qquad
\mu_{ik}^{1}=0,\qquad
\mu_{ik}^{t+1}=\rho\mu_{ik}^{t}
+\bigl(\alpha_L+\delta\mathbf1\{w_t=i\}\bigr)M_{\pi_tk}.
\]

其中 rho in [0,1]，w_t 是本轮赢家，v^t_ik=v^0_ik+mu^t_ik。

本笔记明确采用“取得商品时”的价值：

\[
U_i=\sum_{t=1}^{K}\mathbf1\{w_t=i\}
\bigl(v^0_{i,\pi_t}+\mu_{i,\pi_t}^{t}-p_t\bigr).
\tag{1}
\]

已获得商品的效用不被后续记忆追溯重估，输掉商品本身没有额外直接心理效用。若改成终局重估所有已持有商品，这是不同模型，以下推导须重做。

seller 收入仅为 sum_t p_t。全文收入是指定均衡下的收入；第二价格博弈存在其他均衡的可能，不能省略均衡选择。

## 2. 记忆状态、路径依赖与静态表示

### 命题 1：记忆的显式展开

对任意历史，

\[
\mu_{ik}^{t}=
\sum_{r<t}\rho^{t-1-r}
\bigl(\alpha_L+\delta\mathbf1\{w_r=i\}\bigr)M_{\pi_r k}.
\tag{2}
\]

证明：对更新式递推展开即可。0^0 在相邻轮影响处按 1 处理。

rho=1 时，令 P_t 是已经拍卖的商品集合，S_i^t 是 i 已获得的商品集合，则

\[
\boxed{\mu_{ik}^{t}=\alpha_L\sum_{j\in P_t}M_{jk}
+\delta\sum_{j\in S_i^t}M_{jk}.}
\tag{3}
\]

因此，相同已售集合、相同本人持仓，给出相同本人记忆。相同所有人持仓也确定了已售集合，因此给出相同整个记忆状态。原先计划的“rho=1 下，相同完整持仓但不同拍卖顺序导致不同当前 mu”的证明不成立。

这不意味着顺序无关：后续价值可能相同，但在更早时刻已经取得的商品价值以及支付仍可能不同。

### 命题 2：固定顺序下存在 bundle 表示，衰减也不能消除此事实

设 pos_pi(j) 是 j 在 pi 中的位置，并定义 j 先于 k 时

\[
a_{jk}^{\pi}=\rho^{\operatorname{pos}_{\pi}(k)-\operatorname{pos}_{\pi}(j)-1}M_{jk}.
\]

对于固定 pi，i 的总毛效用恰为

\[
g_i^{\pi}(S)=\sum_{k\in S}
\left(v^0_{ik}+\alpha_L\sum_{j\prec_\pi k}a_{jk}^{\pi}\right)
+\delta\sum_{\substack{j,k\in S\\j\prec_\pi k}}a_{jk}^{\pi}.
\tag{4}
\]

证明：将 (2) 代入 (1)，对所有 i 赢得的商品加总，再把公共位移与本人此前赢得商品引起的差异项分开。每件此前商品总被某人获得，输家系数统一，因此其他人的具体身份不进入 i 的这项毛效用。

固定 pi 后，这就是普通的单项价值加成加两两互补/替代项。rho<1 只改变系数，不改变可表示性。因而“固定顺序下无法表示为任何静态 bundle valuation”是错误命题。

真正可研究的区别是：seller 改变 pi 时，g_i^pi 本身也改变；不是只改变一个预先固定的 bundle valuation 的销售程序。

### 命题 3：rho=1 时跨顺序统一表示的充要条件

要求有一个不依赖 pi 的 g_i(S)，对所有排列与所有分配历史都表示 i 的毛效用。当 n>=2 时，其充要条件为，对每对不同商品 j,k，

\[
\boxed{\alpha_LM_{jk}=0,\qquad
\alpha_W(M_{jk}-M_{kj})=0.}
\tag{5}
\]

必要性证明：

1. 把 j,k 放在前两轮，j 分给另一人，k 分给 i，i 不获得其他商品。交换这两轮而保持最终完整分配不变。i 的毛效用差为 alpha_L M_jk，所以它必须为零。
2. 把 j,k 都分给 i，再交换前两轮。毛效用差为 alpha_W(M_jk-M_kj)，因此也必须为零。

充分性证明：公共位移项由第一个条件消失。同一人获得的每对商品贡献 alpha_W M_jk，而第二个条件保证选择哪个方向都相同。于是

\[
g_i(S)=\sum_{k\in S}v^0_{ik}
+\sum_{\{j,k\}\subseteq S}\alpha_WM_{jk}
\]

成立，其中任取每个无序对的一个方向。

结论：alpha_L=0、rho=1、M 对称时，模型确实退化成普通两两互补估值；如果 alpha_L M_jk 非零，或 alpha_W 不为零且有方向不对称，则不能用一个跨顺序不变的 own-bundle valuation 表示。

### 命题 4：保留衰减可得到真正不同的当前记忆

取 A,B,C，i 在 A,B 两轮均获胜。比较 A->B->C 和 B->A->C，在 C 前有

\[
\mu_{iC}^{AB}=\alpha_W(\rho M_{AC}+M_{BC}),\qquad
\mu_{iC}^{BA}=\alpha_W(\rho M_{BC}+M_{AC}).
\]

因此

\[
\mu_{iC}^{AB}-\mu_{iC}^{BA}
=\alpha_W(\rho-1)(M_{AC}-M_{BC}).
\tag{6}
\]

例如 rho=1/2、alpha_W=1、M_AC=1、M_BC=0，两条历史持仓完全相同，但 C 的记忆增量分别为 1/2、1。这里比较的是不同拍卖顺序的历史；不能据此否认固定 pi 下的 (4)。

## 3. 动态出价究竟能写成什么

### 3.1 物理状态不等于博弈的信息状态

X_t=(mu_1^t,...,mu_n^t) 是物理状态。一般私人信息博弈还需要剩余商品/后续顺序、本人 theta_i，以及由公开历史得到的其他类型后验 B_t。继续效用应写为 V_i^t(theta_i,X_t,B_t;pi)，或直接写在完整信息历史上。

对一般 n，输给 j 与输给 k 会使不同竞争者获得 winner memory，因此不存在未经说明的单一“输的下一状态”。应使用 V_i^{L_j}。

只有在继续博弈已固定、没有依赖当前出价/价格的额外信号效应，且对手相关状态已知时，

\[
v_i+V_i^W-V_i^{L_j}
\]

才是相对输给 j 的确定性临界支付。这是一次胜负比较，不自动成为一般私人信息的出价策略。

### 3.2 两商品的一般最优响应：可以准确写下，但仍是固定点问题

先卖 A，再卖 B。记

\[
x_i=v^0_{iA},\quad y_i=v^0_{iB},\quad
c=\alpha_LM_{AB},\quad d=\delta M_{AB}.
\]

末轮选择真实出价均衡，赢家 A 的 B 估值为 y_i+c+d，其他人为 y_i+c。于是

\[
u_i^W(y)=\left(y_i+d-\max_{\ell\ne i}y_\ell\right)_+,
\]

\[
u_i^{L_j}(y)=\left(y_i-
\max\{y_j+d,\max_{\ell\ne i,j}y_\ell\}\right)_+.
\tag{7}
\]

共同位移 c 从净剩余中消去，但仍进入支付。

记 T_i 是对手最高第一轮出价，J_i 是该对手身份，D_iJ=u_i^W-u_i^{L_J}。在无原子出价环境下，i 出价 b 的期望总效用为

\[
\mathbb E[u_i^{L_{J_i}}\mid x_i,y_i]
+\mathbb E[\mathbf1\{T_i<b\}(x_i-T_i+D_{iJ_i})\mid x_i,y_i].
\tag{8}
\]

因末轮真实出价在任何后验下仍最优，这里不遗漏第一轮出价通过信号改变末轮策略的效应。

若最高对手出价具有正密度且候选最优出价在内部，必要条件为

\[
\boxed{b=x_i+\mathbb E[D_{iJ_i}\mid T_i=b,x_i,y_i].}
\tag{9}
\]

充分性还需条件期望边际收益在 b 附近及全域具有适当单交叉性。不能把 (9) 的必要条件误写成已经求出的均衡。

当 n=3，y_i=3、y_j=0、y_k=2、d=2 时，u_i^W=3，输给 j 时 u_i=1，输给 k 时 u_i=0。对应续期差分别为 2、3，直接展示了输家状态为何依赖赢家身份。

## 4. 两人公开估值：完整临界出价

本节是公开估值解析基准，不是私人信息模型。令 z=y_i-y_j，则

\[
D(z,d)=(z+d)_+-(z-d)_+.
\]

对 d>=0，

\[
D(z,d)=
\begin{cases}
0,&z\le-d,\\
z+d,&-d<z<d,\\
2d,&z\ge d.
\end{cases}
\tag{10}
\]

等价地 D=d+clip(z,-d,d)。因此

\[
\boxed{b_{iA}=x_i+D(y_i-y_j,d).}
\tag{11}
\]

证明：若面对对手出价 q，赢 A 的效用为 x_i-q+(y_i-y_j+d)_+，输 A 的效用为 (y_i-y_j-d)_+。前者减后者等于 x_i+D-q。因此在末轮真实出价的继续博弈中，临界出价就是 (11)。它给出一个子博弈精炼均衡，但不声明所有均衡相同。

bidder 若无论如何都输 B，其溢价为零；若无论如何都赢 B，溢价为 2d：获得 A 不仅提高自己的 B 价值 d，还避免对手提高 B 价值 d。这解释了两倍效应。

若 d<0，D 的符号反向，在出价和末轮估值仍非负的条件下可以得到压价。本文主要精确结果取 d>=0，避免混入负出价和保留价边界问题。

## 5. 核心私人信息结果：强记忆导致第一轮按两件商品估值之和竞争

### 定理 5：任意 n 的强记忆均衡

假设全部 X_i,Y_i 独立 U[0,1]；每人从开始知道自身两项值；c>=0、d>=1。定义 S_i=X_i+Y_i，其 CDF 为

\[
F_S(s)=\begin{cases}
s^2/2,&0\le s\le1,\\
1-(2-s)^2/2,&1<s\le2.
\end{cases}
\]

设 T=max_{j ne i}S_j，Z=max_{j ne i}Y_j，并令

\[
h_n(s)=\mathbb E[Z\mid T=s].
\]

存在一个对称序贯均衡：

\[
\boxed{b_{iA}=\beta_n(S_i),\qquad
\beta_n(s)=s+d-h_n(s),}
\tag{12}
\]

末轮按实际 B 估值真实出价。beta_n 严格递增，因此 A 分给 S_i 最大的人；该人也赢 B（边界平局为零概率）。

**证明，第一步：** d>=1 保证 A 的赢家在 B 中胜过任何输家，即使其 Y=0，而输家 Y=1。故赢 A 的总净效用为

\[
S_i+d-q-Z,
\]

输 A 的后续净效用为 0。

**第二步：** 若所有对手使用 (12)，面对对手最高 score T=t，支付为 beta_n(t)。条件赢得效用为

\[
S_i+d-\beta_n(t)-h_n(t)=S_i-t.
\tag{13}
\]

因此最优地赢过所有 t<S_i，输给所有 t>S_i，恰由出价 beta_n(S_i) 实现。这个论证验证了所有出价偏离，而不仅是局部一阶条件。

**第三步：证明 beta 严格递增。** 给定 T=s，最大 score 对手的 Y 分布为

\[
Y^{top}_s\sim U[\max(0,s-1),\min(1,s)],
\]

其他 n-2 个对手的 Y 独立服从 Y|X+Y<=s。记其 CDF 为 H_s(y)，top 的 CDF 为 L_s(y)，则

\[
h_n(s)=\int_0^1[1-L_s(y)H_s(y)^{n-2}]\,dy.
\tag{14}
\]

H_s(y) 在非平凡范围内为

\[
H_s(y)=\begin{cases}
(sy-y^2/2)/(s^2/2),&0<s\le1,\quad0\le y\le s,\\
y/F_S(s),&1<s<2,\quad0\le y\le s-1,\\
\bigl(sy-y^2/2-(s-1)^2/2\bigr)/F_S(s),&1<s<2,\quad s-1<y\le1.
\end{cases}
\]

top Y 的每个内部 quantile 随 s 单调且斜率严格小于 1。对 truncated Y：s<=1 时，其 quantile 是 s 乘一个 [0,1] 内的固定数；s>1 时由 H_s(q_s(p))=p 隐式求导可得 0<=dq_s/ds<1（内部 p）。在 q>s-1 的分段，令 r=s-1、A=(1+r)q-q^2/2-r^2/2、F=(1+2r-r^2)/2，有

\[
\frac{dq}{ds}=\frac{p(1-r)-(q-r)}{s-q}.
\]

非负性由

\[
A(1-r)-(q-r)F
=\tfrac12(1-q)[(1-r)q+r+r^2]\ge0
\]

保证，上界 1 则由 p(1-r)<=1 保证，内部严格。用相同独立 quantiles 耦合不同 s 的所有竞争者 Y，其最大值增长严格小于 s 的增长。故 h_n(s')-h_n(s)<s'-s，beta 严格递增。证毕。

注意 h_n 是由基础分布明确决定的一维积分，不是未解的均衡固定点。端点 h_n(0)=0、h_n(2)=1；n=2 时 h_2(s)=s/2。

### 推论 5.1：两人闭式出价

\[
\boxed{b_{iA}=d+\frac{x_i+y_i}{2}.}
\tag{15}
\]

例如 d=1，(x_1,y_1)=(0.2,0.9)、(x_2,y_2)=(0.7,0.1)，出价为 1.55、1.40；当前 A 估值更低的 bidder 1 赢 A。出价中出现 y_i 是严格的最优响应结果，而非把未来价值直接猜测性加到 x_i 上。

### 推论 5.2：任意 n 的总收入

记 S_(n-1:n) 为 S_i 的第二大值，则

\[
\boxed{R(A\to B)=c+d+\mathbb E[S_{(n-1:n)}].}
\tag{16}
\]

直接证明：令第一轮 runner-up score 为 s。第一轮支付 beta_n(s)=s+d-h_n(s)。条件于这个 runner-up score，所有输家的最高 Y 的期望恰为 h_n(s)，因此末轮支付的条件期望为 c+h_n(s)。两项相加为 c+d+s，再取期望即可。

等价地，赢家取得两件商品，毛效用为 S_max+c+d。其 interim utility 为 E[(S_i-T)_+]，与 score S 的普通第二价格拍卖相同，亦可由收入恒等式推出 (16)。

两人时，

\[
\mathbb E\min(S_1,S_2)=\int_0^2(1-F_S(s))^2ds=\frac{23}{30},
\]

故

\[
\boxed{R(A\to B)=\frac{23}{30}+c+d
=\frac{23}{30}+\alpha_WM_{AB}.}
\tag{17}
\]

如果反向也满足 delta M_BA>=1，两个商品均为 U[0,1]，则

\[
R(A\to B)-R(B\to A)=\alpha_W(M_{AB}-M_{BA}).
\tag{18}
\]

这说明强记忆区间的出价非常战略化，但对称分布下的顺序条件反而简单。不能为了故事而把它包装成复杂 threshold。

更一般地，若 iid pair (X,Y) 有界且 d 大到保证 A 赢家赢 B，且定义的 beta(s)=s+d-E[max Y_j|max(X_j+Y_j)=s] 严格递增，同一最优响应和收入证明仍成立。均匀情形已在上面验证单调性；任意分布则不能省略这个条件。

### 推论 5.3：在强记忆阈值下方仍有可证明的近似均衡

保留均匀分布与提前私人信息，取 0<=d<=1，仍使用 beta_n(s)=s+d-h_n(s)，末轮真实出价。则第一轮策略组成一个 epsilon-Bayes-Nash 均衡，其中 epsilon=1-d；末轮策略精确最优。

证明：定义一个仅用于比较的辅助收益：赢 A 时收益为 s+d-q-Z（暂不把末轮剩余截断为非负），输 A 时为 0。第 5 节的条件比较 (13) 对任意 d 仍成立，所以 beta 是这个辅助问题的最佳响应。在真实问题中，赢 A 时多出的收益为 (Z-y_i-d)_+，输 A 时的真实继续剩余也至多 1-d。因此，任何第一轮出价的真实收益与辅助收益之差都在 [0,1-d]，最优偏离相对 beta 至多提高 1-d。

该策略组合的实际收入还满足

\[
c+d+r_n(F_S)-(1-d)\le R_{\beta}(A\to B)
\le c+d+r_n(F_S).
\tag{18a}
\]

因为第一轮支付不变，末轮第二高价相对于“c+最高输家 Y”只会下降，且下降至多 1-d。这是明确的近似均衡与策略收入界，仍不能替代弱记忆区间的精确均衡定理。

## 6. 未来估值延后揭示：任意 n 的全强度精确基准

本节明确更换信息结构：第一轮时本人也尚不知道 Y_i，第二轮前才揭示；Y_i iid F_B 且与第一轮信息独立。取 c,d>=0。此条件不能暗中带回第 5 节的默认模型。

定义

\[
A_m(d)=\mathbb E\max(Y_1+d,Y_2,\ldots,Y_m),\qquad
G_m(d)=A_m(d)-A_m(0),\quad m\ge1.
\]

### 定理 6：动态溢价与总收入

\[
\boxed{b_{iA}=x_i+G_{n-1}(d).}
\tag{19}
\]

记 r_n(F)=E[第二大 iid F 样本]，则

\[
\boxed{R(A\to B)=r_n(F_A)+r_n(F_B)+c+H_n(d;F_B),}
\tag{20}
\]

其中

\[
\boxed{H_n(d;F)=nG_{n-1}(d)-(n-1)G_n(d).}
\tag{21}
\]

证明：A 赢家在 B 中的期望净剩余为 A_n(d)-A_(n-1)(0)；一个指定输家的期望净剩余为 A_n(d)-A_(n-1)(d)。两者之差恰为 G_(n-1)(d)，与所有第一轮出价和类型无关，因此 (19) 是全域最佳响应。

末轮赢家总价值为 c+A_n(d)，所有 bidder 的净剩余之和为

\[
nA_n(d)-A_{n-1}(0)-(n-1)A_{n-1}(d).
\]

由总价值=支付+总净剩余，末轮相对无记忆收入增量为 c+(n-1)[G_(n-1)-G_n]。再加第一轮溢价 G_(n-1)，得 (20)-(21)。

若 F 连续，则

\[
G_m'(d)=\int F(y+d)^{m-1}\,dF(y),
\]

\[
H_n'(d)=\int F(y+d)^{n-2}
[n-(n-1)F(y+d)]\,dF(y)\ge0.
\tag{22}
\]

特别地，G_m'(0)=1/m，

\[
H_n'(0)=\frac{2n-1}{n(n-1)}.
\tag{23}
\]

这也是“无参与费并不使 n 无关”的一个精确反例。

### 均匀分布的闭式

若 F_B=U[0,1]，则

\[
G_m(d)=\begin{cases}
\frac{d^2}{2}+\frac d m-\frac{d^{m+1}}{m(m+1)},&0\le d\le1,\\
d-\frac12+\frac1{m+1},&d\ge1.
\end{cases}
\tag{24}
\]

因为 G_m'(d)=d+(1-d^m)/m（0<=d<=1），对 d 积分即可。

两人时 G_1(d)=d，第一轮出价是 x_i+d。但总收入仍是非线性的：若两个商品都均匀，

\[
R(A\to B)=\frac23+c+g(d),\qquad
g(d)=\begin{cases}
\frac32d-\frac12d^2+\frac16d^3,&0\le d\le1,\\
d+\frac16,&d\ge1.
\end{cases}
\tag{25}
\]

直接核验末轮收入：E[min(Y_1+d,Y_2)] 在 d<=1 时为 1/3+d/2-d^2/2+d^3/6，d>=1 时为 1/2。

两种顺序的一般差为

\[
\boxed{R(A\to B)-R(B\to A)
=\alpha_L(M_{AB}-M_{BA})
+H_n(\delta M_{AB};F_B)-H_n(\delta M_{BA};F_A).}
\tag{26}
\]

这里出现了记忆方向、赢家输家差异和下一个商品估值分布的共同作用。若 F_A=F_B，且 alpha_L,delta>=0，H 单调意味着较大方向边更有利；不能在相同分布下硬造更复杂条件。

## 7. 信息时点的收入影响会随人数反转

比较第 5 节与第 6 节：X,Y 均为 U[0,1]，同样 c>=0、d>=1，仅改变 bidder 何时得知自身 Y。

提前知道 Y 的收入为

\[
R_{early}=c+d+r_n(F_S).
\]

延后得知 Y 的收入为

\[
R_{late}=c+d+\frac12+\frac{n-1}{n+1}.
\tag{27}
\]

因此差值是

\[
\boxed{R_{early}-R_{late}
=r_n(F_S)-\frac12-\frac{n-1}{n+1}.}
\tag{28}
\]

| n | r_n(F_S) | 提前减延后收入 |
|---|---:|---:|
| 2 | 23/30 | -1/15 |
| 3 | 1 | 0 |
| 4 | 947/840 | 23/840 |
| 5 | 611/504 | 23/504 |

这些为精确有理数，不是 simulation。仅 n=2 与 n=4 已足以证明“更多未来私人信息对收入的方向并不固定，并随竞争人数变化”。未在此声明对所有 n>=4 的符号定理。

计算方法：令 I_m=integral_0^2 F_S(s)^m ds，则

\[
I_m=\frac{1}{2^m(2m+1)}
+\sum_{j=0}^m\binom mj\frac{(-1)^j}{2^j(2j+1)},
\]

\[
r_n(F_S)=2-nI_{n-1}+(n-1)I_n.
\tag{29}
\]

两人时，提前信息也把期望毛总价值由 c+d+7/6 提高到 c+d+37/30，提高 1/15，同时收入降低 1/15，因此 bidder 的总期望剩余提高 2/15。该比较使用 (1) 的实现价值福利，不等于关于“心理福利”的外部规范判断。

## 8. 一般 K：哪些环境确实能化成加权排序

以下使用 rho=1。定义

\[
W(\pi)=\sum_{j\prec_\pi k}w_{jk}.
\]

### 定理 8A：公共位移边界

若 alpha_W=alpha_L=a>=0，记忆与谁赢无关。在每轮按实际当前价值真实出价的均衡中，

\[
R(\pi)=C+a\sum_{j\prec_\pi k}M_{jk},\qquad
C=\sum_k r_n(F_k).
\tag{30}
\]

证明：过去出价不改变任何未来人的价值；选择末轮起逐轮真实出价是一个序贯均衡。每轮所有出价共同增加 a sum_{j prior} M_jk，第二高出价增加同样的数。求和即可。

这是可用的控制组，但 delta=0 不满足研究中心的赢家输家不对称条件，不能作为核心模型已解决的证据。

### 定理 8B：保留赢家/输家不对称的一般 K 精确子类

使用两位 bidder、公开确定性基础估值：每件商品对 bidder 1 值 H，对 bidder 2 值 0。取 alpha_L=1、alpha_W=2（delta=1），M 非负，且

\[
H>\sum_{j\ne k}M_{jk}.
\tag{31}
\]

即使经历任意过去分配，bidder 1 的每件剩余商品当前价值仍高于 bidder 2。存在下面的子博弈精炼均衡：当前商品 j=pi_t，

\[
b_{2t}=v_{2j}^{t},\qquad
b_{1t}=v_{1j}^{t}+2\delta\sum_{k:\,j\prec_\pi k}M_{jk}.
\tag{32}
\]

bidder 1 赢得所有剩余商品，在均衡路径上 seller 收入为

\[
\boxed{R(\pi)=\alpha_L\sum_{j\prec_\pi k}M_{jk}.}
\tag{33}
\]

证明：由 H 的界，任何历史下 bidder 1 都有严格当前价值优势。对剩余轮数归纳。假设此后的所有商品都会由 bidder 1 获得：bidder 2 无论本轮输赢，未来净效用都是 0，所以临界出价是当前价值。bidder 1 本轮胜出，相对于让 bidder 2 胜出，会使每件未来商品的自身价值增加 delta M_jk、对手出价减少 delta M_jk，继续净效用增加 2 delta sum M_jk。因此其临界出价为 (32)，依旧严格高于 bidder 2。归纳成立。在实际路径上 bidder 2 每次都输，其出价是 alpha_L sum_{j prior}M_jk；加总为 (33)。

该结果保留非对称记忆，但放宽了“bidder 间 iid、非退化私人信息”的条件。不能宣称它证明了原始 iid 子类的 NP-hardness。

### 推论 8C：NP-hardness

给定一个非负整数加权有向图，把其边权设为 M_jk，构造 (31)-(33) 的环境，令 H=1+sum M_jk，其表示长度为多项式。每个排列的 seller 收入恰为该排列的前向边总权重。因此 maximum weighted acyclic subgraph / maximum linear ordering 的最优解与最优拍卖顺序一一对应。

所以一般的公开异质基础估值版本，即使 n=2、alpha_L=1、alpha_W=2 固定，最优排序仍 NP-hard。此归约没有假定一个未证明的 pairwise 收入表达式。

这里使用 (32) 明确规定的动态临界出价均衡选择；没有声称任意其他第二价格均衡都给出同样收入，也没有处理 seller 选择最有利均衡的另一种优化定义。

同样，公共位移版本通过 (30) 得到 hardness。但这两条都不能无条件覆盖全部最初的 iid 非对称私人信息设定。

## 9. 可分解子类中的算法

### 9.1 确定性 1/2 近似

对任意排列 pi 与其反向 pi^rev，

\[
W(\pi)+W(\pi^{rev})=\sum_{j\ne k}w_{jk}=W_{all}.
\]

选二者中较好者，就有

\[
\boxed{W(\widehat\pi)\ge\tfrac12W_{all}
\ge\tfrac12\max_\pi W(\pi).}
\tag{34}
\]

无需随机化；随机排列也满足 E W=W_all/2。若 R=C+W 且 C>=0，该算法同时保证总收入至少最优总收入的一半。对增量 W 给出保证更有解释力，避免大的常数基线掩盖排序质量。

### 9.2 反对称部分决定顺序

每个无序对 {j,k} 的常数为 (w_jk+w_kj)/2，剩余排序偏好只取决于 w_jk-w_kj。于是完全对称权重时所有顺序等价。

若只保留每对较大权重的方向后，得到的偏好图无环，则其任何拓扑序逐对取得 max(w_jk,w_kj)，全局最优。原始非负支持图为 DAG 是充分条件；更强的是“较优方向图”为 DAG 即可。

若非零权重的无向支撑为森林，每条边选较优方向不会形成有向环，因此可精确求解。这一结论适用于双向边不同权重的森林。

### 9.3 精确子集 DP

令 D(S) 是 S 内商品的最优前向边总权重。把 k 放最后：

\[
D(\varnothing)=0,\qquad
\boxed{D(S)=\max_{k\in S}\left[D(S\setminus\{k\})
+\sum_{j\in S\setminus\{k\}}w_{jk}\right].}
\tag{35}
\]

预计算子集入边和后，时间 O(K 2^K)，保守空间 O(K 2^K)；可用额外实现技巧节省内存。无需预计算时直接计算内层和得到 O(K^2 2^K)。这是小中等 K 的 exact benchmark。

## 10. 鲁棒性：已能给出的严格界

### 10.1 估计权重误差

设 E=sum_{j ne k}|w_jk-hat w_jk|。每个排列满足 |W(pi)-hat W(pi)|<=E。若 hat pi 精确最优化估计目标，则

\[
\boxed{W(\pi^*)-W(\widehat\pi)\le2E.}
\tag{36}
\]

若算法对估计目标为 gamma 近似，则

\[
W(\widehat\pi)\ge\gamma W(\pi^*)-(1+\gamma)E.
\tag{37}
\]

这只是已有估计误差下的优化稳定性，不提供 M 的识别或学习样本复杂度。在原始战略均衡可能跳变的模型中，不能直接套用这个 Lipschitz 界。

### 10.2 衰减

公共位移子类以及第 8B 节强优势子类，在 rho<1 时分别得到相应的

\[
R_\rho(\pi)=C+W_\rho(\pi),\qquad
W_\rho(\pi)=\sum_{s<t}\rho^{t-s-1}w_{\pi_s\pi_t}.
\tag{38}
\]

第 8B 节的 bidder 1 继续溢价相应变为 2 delta sum_{t'>t}rho^{t'-t-1}M_{pi_t,pi_t'}，同样归纳成立。

rho<1 时权重依赖两件商品的间隔，普通 pairwise linear ordering 结构不再精确。特别地 rho=0，只保留相邻项，成为加权有向 Hamiltonian path 目标。原先的 DAG/对称权重特例也不能照搬。

若 W_all=sum w_jk 且所有权重非负，Bernoulli 不等式给出

\[
0\le W_1(\pi)-W_\rho(\pi)
\le(1-\rho)(K-2)W_{all}=:\varepsilon_\rho.
\tag{39}
\]

若算法对 W_1 为 gamma 近似，则

\[
\boxed{W_\rho(\widehat\pi)
\ge\gamma\max_\pi W_\rho(\pi)-\varepsilon_\rho.}
\tag{40}
\]

因是单侧差异，此处只需减一个 epsilon，而不是随意使用两倍界。该界在 rho 接近 1 时有用；当误差项过大时它可能无信息。

### 10.3 原始一般 K 私人信息模型的近似均衡

对固定 pi，取 rho=1、M>=0、delta>=0。考虑所有人都只报“基础值+公共位移”的策略：

\[
b^0_{it}=v^0_{i,\pi_t}+\alpha_L\sum_{r<t}M_{\pi_r,\pi_t}.
\tag{41}
\]

这些策略故意忽略 winner-specific 记忆，不应称为真实出价。定义

\[
B_\pi=\sum_{j\prec_\pi k}M_{jk}.
\]

则该策略组合是一个 epsilon-Bayes-Nash 均衡，甚至偏离收益对固定他人类型也满足

\[
\boxed{\varepsilon=\delta B_\pi.}
\tag{42}
\]

证明：对手使用 (41) 时，其未来出价不依赖当前分配/信号。在 delta=0 的游戏中，每轮 (41) 是最佳响应；任何跨期偏离都不能提高总效用。delta>0 时，固定出价行动历史的新增毛效用为

\[
\delta\sum_{j\prec_\pi k}M_{jk}
\mathbf1\{i\text{ wins }j\}\mathbf1\{i\text{ wins }k\},
\]

它位于 [0,delta B_pi]。偏离最多获得这一幅度的额外收益，原有策略本身的额外收益非负，所以总偏离收益至多 delta B_pi。

该近似策略下的收入恰为 C+alpha_L B_pi。但这不是“真正精确均衡收入与 pairwise 目标相近”的证明，更不能直接声称对精确均衡最优收入有 1/2 近似。

## 11. 对论文故事的修正

已有序贯拍卖文献明确研究了未来拍卖机会引起的外部性与战略出价。因此“第二价格拍卖不再静态 truthful”“继续效用差进入当前策略”本身不能作为新颖性结论。

Paes Leme, Syrgkanis, Tardos, [*Sequential Auctions and Externalities*](https://arxiv.org/abs/1108.2452)。

Maximum Acyclic Subgraph 的已知困难性可引用：Austrin, Manokaran, Wenner, [*On the NP-Hardness of Approximating Ordering Constraint Satisfaction Problems*](https://arxiv.org/abs/1307.5090)。

上文的归约、均衡、表示性分类与信息时点比较是本笔记独立推导的待研究结果；已检索的这些文献不等于完成了全面新颖性核查。尤其强记忆下的“隐式打包”与固定顺序的 synergy 表示，应与既有 sequential synergy / bundling 文献详细比对后，才能声称是新定理。

当前更精确的候选论文主线是：

1. **模型表示性边界。** 哪些 memory 更新仅仅重写了互补估值，哪些会令 seller 的顺序改变偏好本身？
2. **战略行为。** 在明确的私人信息模型里，未来估值进入当前竞价；强记忆下出现可完全刻画的隐式打包。
3. **信息与竞争。** 未来私人估值的揭示时点对收入影响可随 bidder 数反转。
4. **算法桥梁。** 在什么条件下战略收入可分解为图排序？何时必须保留更复杂的状态依赖？

如果最终最强调“相同持仓、不同当前记忆”，应把 rho<1 放入主模型；但即使如此，固定顺序的 bundle 表示仍成立。需要明确比较跨顺序/自适应历史，而不是声称无法表示任何固定顺序游戏。

## 12. 尚未证明的具体缺口

- 全部私人估值提前已知、0<d<1 时，两商品的一般纯策略均衡存在性、唯一性与解析/可验证数值解。(9) 是起点，不是答案。
- 原始 iid、alpha_W!=alpha_L、一般 K 环境中，一个可证明的均衡选择、收入计算方法，以及排序近似保证。
- 从“公开异质强优势”hardness 升级到 iid 私人信息子类的保结构归约。
- rho<1 时的一般 K 正结果，特别是保留真实战略效应的结果，而非仅公共位移/强优势环境。
- 自适应 seller 的策略问题：bidder 若预见 seller 会根据出价/结果调整顺序，均衡也会变化，不能拿固定顺序的 continuation 直接拼成 Bellman policy。
- unknown M 的识别、估计与统计误差：目前只有给定权重误差时的优化界。
- 全面最近邻检索与正式新颖性判断。本笔记不声称已有可直接投稿的完整 theorem package。

## 13. 计算核验记录

计算是独立检查，不代替上文证明。

- 对 n=2,3,5，各生成 500,000 组独立均匀类型；设 c=0.3、d=1.2；第一轮按 (12) 出价，第二轮按真实价值竞价。
- 总收入理论值分别为 2.2666667、2.5、2.7123016；样本均值分别为 2.2675646、2.4996912、2.7122490，差异均在合理抽样误差内。
- 在 1001 点 score 网格上核查 beta 严格递增。
- 对 d=0,0.2,0.7,1,1.5，以数值积分核验 (25) 的末轮收入。
- 构造 5 件商品的非负整数图，遍历全部 120 个排列，并对每个排列递归检查所有二元赢家历史。动态临界出价与 (32) 一致，总收入与 (33) 一致。
- 子集 DP 与上述全排列最优值一致，均为 30。
- 使用有理数运算核验表 (28) 的信息时点差值。

以下附计算脚本，便于复查。依赖 Python、NumPy、SciPy；它不用于推断弱记忆区间的一般均衡。

```python
import itertools
import json
from fractions import Fraction as Q
import numpy as np
from scipy.integrate import quad

def fs(s):
    if s <= 0: return 0.
    if s <= 1: return s*s/2
    if s < 2: return 1-(2-s)**2/2
    return 1.

def h(s,n):
    if s <= 0: return 0.
    if s >= 2: return 1.
    lo,hi=max(0.,s-1),min(1.,s)
    def integrand(y):
        top=np.clip((y-lo)/(hi-lo),0,1)
        if s <= 1:
            trunc=1. if y>=s else (s*y-y*y/2)/fs(s)
        else:
            trunc=(y if y<=s-1 else s*y-y*y/2-(s-1)**2/2)/fs(s)
        return 1-top*trunc**(n-2)
    return quad(integrand,0,1,points=[lo,hi],epsabs=1e-10)[0]

rng=np.random.default_rng(20260919)
out={}
grid=np.linspace(0,2,1001)
for n in [2,3,5]:
    hv=np.array([h(s,n) for s in grid])
    bids=grid+1.2-hv
    assert np.min(np.diff(bids))>0
    x=rng.random((500000,n));y=rng.random((500000,n));s=x+y
    winner=np.argmax(s,axis=1)
    second=np.partition(s,-2,axis=1)[:,-2]
    pa=np.interp(second,grid,bids)
    y[np.arange(len(y)),winner]=-np.inf
    rev=pa+.3+np.max(y,axis=1)
    target=1.5+quad(lambda z:1-fs(z)**n-n*(1-fs(z))*fs(z)**(n-1),0,2,points=[1])[0]
    se=rev.std()/len(rev)**.5
    assert abs(rev.mean()-target)<6*se+1e-5
    out[f'strong_memory_n{n}']={'expected':target,'MC':float(rev.mean()),'SE':float(se),'min_bid_grid_increment':float(np.diff(bids).min())}

out['exact_n2_second_sum']=str(Q(23,30))
for d in [0.,.2,.7,1.,1.5]:
    numeric=quad(lambda z:(1-z)*(1-max(0,z-d)),0,1,points=[min(1,d)])[0]
    exact=1/3+d/2-d*d/2+d**3/6 if d<=1 else .5
    assert abs(numeric-exact)<1e-10
    out[f'delayed_shift_{d}']=exact

K=5
M=rng.integers(0,5,(K,K)).astype(float)
np.fill_diagonal(M,0)
delta=1.;alpha_l=1.;H=1+M.sum()
def solve(pi,mu,t=0):
    if t==K:return np.zeros(2),0.
    j=pi[t]
    nxt=[]
    for winner in [0,1]:
        updated=mu+alpha_l*M[j]
        updated=updated.copy();updated[winner]+=delta*M[j]
        nxt.append(solve(pi,updated,t+1))
    v=np.array([H,0.])+mu[:,j]
    b=np.array([v[0]+nxt[0][0][0]-nxt[1][0][0],v[1]+nxt[1][0][1]-nxt[0][0][1]])
    expected=np.array([v[0]+2*delta*sum(M[j,k] for k in pi[t+1:]),v[1]])
    assert np.allclose(b,expected)
    assert b[0]>b[1]>=0
    winner=int(np.argmax(b));price=b[1-winner]
    util=nxt[winner][0].copy();util[winner]+=v[winner]-price
    return util,price+nxt[winner][1]

brute=[]
for pi in itertools.permutations(range(K)):
    _,r=solve(pi,np.zeros((2,K)))
    target=sum(M[pi[s],pi[t]] for s in range(K) for t in range(s+1,K))
    assert abs(r-target)<1e-8
    brute.append(target)
dp={0:0.}
for mask in range(1,1<<K):
    dp[mask]=max(dp[mask^(1<<k)]+sum(M[j,k] for j in range(K) if j!=k and mask>>j&1) for k in range(K) if mask>>k&1)
assert dp[(1<<K)-1]==max(brute)
out['generalK_dominance']={'permutations_verified':len(brute),'optimum':max(brute),'DP':dp[(1<<K)-1]}
print(json.dumps(out,indent=2))
```
