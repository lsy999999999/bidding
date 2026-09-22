# 带正负边的记忆拍卖：收入分解、赢家更替均衡与有限 K 排序

日期：2026-09-20。推导稿；证明与数值核验不等于文献新颖性审查。

## 结果范围

- 任意有限 n、K、任意带符号 M：得到精确收入分解与条件加性收入保证。均衡存在性与选择仍须单独处理。
- n=2、K=2、独立均匀私人值：在强正/负影响区间得到显式 Bayesian 均衡；负影响时两件商品由不同 bidder 获得。
- n=2、任意有限 K、关联图为互不相连的商品对及孤立点：拼接上述均衡，精确计算两个修正项，获得所有排列上的收入最优算法，允许任意交错出售各组商品。
- 未解决：任意 n 的显式负影响均衡；一般连通有符号图；弱影响区间；一般图上由原始参数推出小的买家租金变化界；所有均衡统一的收入保证。

## 1. 模型与负边的含义

卖家事先承诺排列 pi，每件商品只出售一次。所有 bidder 每轮参加二价拍卖，非负出价、零底价、固定平局规则，可购买多件，无预算约束。bidder 事先知道自己的完整基础估值向量，其他 bidder 的向量未知。记忆、赢家和支付规则公知。

保持 rho=1，alpha_W>alpha_L>=0，delta=alpha_W-alpha_L。M 的非对角元允许任意符号。基础价值记 x_ik，取得时价值及总效用为

\[
v_{ik}^t=x_{ik}+\sum_{r<t}[\alpha_L+\delta\mathbf1\{w_r=i\}]M_{\pi_rk},
\qquad U_i=\sum_t\mathbf1\{w_t=i\}(v_{ik}^t-p_t),\quad k=\pi_t.
\]

负 M_jk 表示经历 j 后对 k 的估值下降；由于 alpha_W>alpha_L，赢 j 的人下降更多。alpha_L=0 时仅赢家受到该次影响；alpha_L>0 时输家也会下降。不能把所有负边自动称为满足经济学 formal gross-substitutes 条件的估值。

负边不必然使赢家更替：基础估值优势、其他记忆和前瞻策略都可能抵消它。下面给出确保更替的具体充分条件。

为避免负估值和零底价改变最后一轮支付表达，可取 x_ik=a_k+u_ik，u_ik>=0，并要求

\[
a_k>\alpha_W\sum_{j\ne k}(-M_{jk})_+.
\]

这保证所有历史下估值为正，保留线性模型；不是额外进行 max(v,0) 截断。截断、拒绝分配或保留价应作为另一机制重新分析。

## 2. 任意 n、K 的精确收入分解

令

\[
W(\pi)=\sum_{r<t}M_{\pi_r\pi_t},\quad
L(\pi,h)=\sum_{r<t}M_{\pi_r\pi_t}\mathbf1\{w_r\ne w_t\}.
\]

L 现在有正有负。沿每一条实际历史 h，

\[
\sum_t v_{w_t,\pi_t}^t
=\sum_t x_{w_t,\pi_t}+\alpha_W W(\pi)-\delta L(\pi,h).
\tag{1}
\]

证明：每一条先后边在最终赢家的价值中贡献 alpha_W M（两个端点同一赢家）或 alpha_L M（不同赢家），分别写成 alpha_W M-delta M 1{不同赢家} 后求和。

指定每个排列下所研究的均衡及一致的选择规则；假设预期效用和收入有限。设

\[
A(\pi)=\mathbb E\sum_t x_{w_t,\pi_t}-\mathbb E\sum_i U_i.
\]

则由总收入=总毛价值-总净效用得到

\[
R(\pi)=\alpha_W W(\pi)-\delta\mathbb E L(\pi,h)+A(\pi).
\tag{2}
\]

该恒等式甚至不要求策略是均衡；但用来比较经济预测时必须使用所研究的均衡。它不证明一般均衡存在，也不保证 A 与顺序无关。

## 3. 正确的带符号排序误差界

写 L=L^+-L^-，其中

\[
L^\pm=\sum_{r<t}(\pm M_{\pi_r\pi_t})_+\mathbf1\{w_r\ne w_t\}\ge0.
\]

令 pi* 最大化 R，hat pi 满足 W(hat pi)>=max_pi W(pi)-eta。设

\[
\beta=\max_\pi A(\pi)-\min_\pi A(\pi),\qquad
N_-=\sum_{j\ne k}(-M_{jk})_+.
\]

从 (2) 相减并舍去非正项得到

\[
\boxed{
R(\pi^*)-R(\widehat\pi)
\le\alpha_W\eta+\beta+
\delta[\mathbb E L^+(\widehat\pi)+\mathbb E L^-(\pi^*)]
\le\alpha_W\eta+\beta+
\delta[\mathbb E L^+(\widehat\pi)+N_-].
}
\tag{3}
\]

不能沿用非负图中 L(pi*)>=0 的推理。

另一对称版本：若对所有 pi 有 E|L(pi,h)|<=ell，则 regret<=alpha_W eta+beta+2 delta ell。若每条被考虑的先后边都有端点赢家不同的概率至多 epsilon，则 E|L|<=epsilon sum_{r<t}|M|，可代入；但该概率条件不是负边自动保证的。

这些是条件保证：一般图上必须进一步用已知的图与分布参数界定 beta，才能得到实际算法保证。仅把未知均衡量改名不是解决。

### 带符号 W 的图算法

有负边时，对原始 W 或 R 直接声称乘法 1/2 近似可能没有意义。对每对无序商品 {j,k}，取 c_jk=min(M_jk,M_kj)，d_jk=|M_jk-M_kj|。偏好的方向是系数较大的方向，则

\[
W(\pi)=\sum_{j<k}c_{jk}+\sum_{\text{偏好方向被满足}}d_{jk}.
\]

随机排序取得后一个非负增量目标的一半期望值；取任意排列及其反序中较好的一个，也至少取得全部 d 权重的一半。该保证不直接等于收入的 1/2 保证。精确子集 DP 为

\[
DP[S]=\max_{k\in S}\{DP[S\setminus\{k\}]+\sum_{j\in S\setminus\{k\}}M_{jk}\},
\]

预计算子集入边和后时间 O(K 2^K)。所得 W 的加性误差 eta 可代入 (3)。

## 4. 两件商品的显式强负影响均衡

这里开始限定 n=2。A 在 B 前，设

\[
x_{iA}=a+\xi_i,\quad x_{iB}=b+\zeta_i,\quad
\xi_1,\xi_2,\zeta_1,\zeta_2\ \text{iid }U[0,1].
\]

设 M_AB=-m<0，m>0，D=delta m>=1，a>=D，b>=alpha_W m（可取严格不等式避免边界）。M_BA 不影响此顺序。

如果 i 赢 A，则 B 时

\[
v_{iB}=b+\zeta_i-\alpha_Wm,\qquad
v_{jB}=b+\zeta_j-\alpha_Lm.
\]

由于 \zeta_i-\zeta_j-D<=0，j 在 B 中获胜，除零概率边界平局外。最后一轮报实际估值。因此

\[
U_i(\text{赢 A})=a+\xi_i-p_A,\quad
U_i(\text{输 A})=D+\zeta_i-\zeta_j.
\tag{4}
\]

后式是 bidder 之后赢得 B 时价值减支付的净效用，不是输掉 A 当下的惩罚或奖励。

定义私人得分 S_i=\xi_i-\zeta_i，取第一轮策略

\[
\boxed{b_{iA}=a-D+\frac{1+\xi_i-\zeta_i}{2}.}
\tag{5}
\]

### 全局最佳响应证明，不仅是一阶条件

对手使用严格递增 beta(s)=a-D+(1+s)/2。任意出价等价于选择对手得分阈值 t（支持外对应永赢/永输）。独立均匀分布给出

\[
\mathbb E[\zeta_j\mid S_j=t]=(1-t)/2,\quad -1<t<1.
\]

在边界 S_j=t，从输 A 改成赢 A 的条件期望净收益为

\[
a+\xi_i-\beta(t)-[D+\zeta_i-\mathbb E(\zeta_j\mid S_j=t)]
=S_i-t.
\]

所以若 G_i(t) 是选阈值 t 的期望总效用，f_S(t)=1-|t|，

\[
G_i(t)-G_i(S_i)=\int_{S_i}^{t}(S_i-z)f_S(z)\,dz\le0.
\tag{6}
\]

这验证所有第一轮出价偏离。结合最后一轮真实出价，得到该两轮私人信息拍卖的 Bayesian 均衡及序贯理性的最后一轮策略。

### 收入和租金

赢家是 S 最大者，记 w；另一人记 l。第一轮价格和第二轮价格分别是

\[
p_A=a-D+(1+S_l)/2,\quad p_B=b+\zeta_w-\alpha_Wm.
\]

E max S=7/30，E min S=-7/30，E \zeta_w=(1-E max S)/2=23/60。因此

\[
\boxed{R(A\to B)=a+b+\frac{23}{30}-(2\alpha_W-\alpha_L)m.}
\tag{7}
\]

此时 E 基础毛价值=a+b+37/30，实际记忆毛价值= -alpha_L m，故

\[
\mathbb E\sum_i U_i=\frac7{15}+2\delta m,
\qquad A(\pi)=a+b+\frac{23}{30}-2\delta m,
\qquad L=-m.
\tag{8}
\]

这精确算出了两个修正项。分给不同人避免了额外的记忆损失（-delta L 为正），但预期净效用增长 2 delta m，使卖家收入的净系数仍为负。不能仅分析总价值或 L 就宣称收入提高。

## 5. 与强正影响的统一表达

若 M_AB=m>0 且 delta m>=1，则第一轮赢家也赢 B；第一轮策略为

\[
b_{iA}=a+\delta m+(\xi_i+\zeta_i)/2.
\]

证明同理，得分 S=\xi+\zeta，E[\zeta_j|S_j=t]=t/2，阈值边际收益仍为 S_i-t。此时

\[
R=a+b+23/30+\alpha_Wm,\quad
\mathbb E\sum_i\xi_i=7/15,\quad A=a+b+23/30,\quad L=0.
\]

将有符号 m 直接作为 M_AB，定义仅在强影响区间使用的函数

\[
\psi(m)=
\begin{cases}
\alpha_Wm,&m\ge1/\delta,\\
(2\alpha_W-\alpha_L)m,&m\le-1/\delta.
\end{cases}
\]

则 R=a+b+23/30+psi(M_AB)。如果 M_AB、M_BA 都在该强影响区间，并且 a,b 足够大以保证两种顺序的价格和估值非负，

\[
R(A\to B)-R(B\to A)=\psi(M_{AB})-\psi(M_{BA}).
\]

psi 在这里覆盖的两个区间上严格递增，因此单个孤立商品对选择较大的 M 方向。这不表示一般有符号图上的问题也只需逐边选择；一般图会出现方向冲突、分配和租金的相互作用。

零边的两件独立商品收入为 a+b+2/3。不能把强影响公式中的 23/30 外推到 m=0 或 0<|delta m|<1。

## 6. 任意有限 K 的精确可解子类

仍取 n=2。图的弱连通分量为 q 对商品和 r 个孤立点，K=2q+r。不同商品对间 M 为零。每个商品对的两个方向均满足 |delta M_jk|>=1，所有基础估值为 x_ik=a_k+U_ik，U_ik 跨 bidder、跨商品全部独立 U[0,1]。对每对 {j,k}，令 a_j、a_k 都严格大于 alpha_W max(|M_jk|,|M_kj|)，从而两种方向的估值和首轮报价都为正。只有负入边界可以保证当前价值非负，却未必保证首轮动态报价非负，因此这里使用这个同时覆盖入边与出边的充分条件。

在每个商品对第一次出现时按第 4/5 节策略出价，第二次出现时报当前价值；孤立商品真实出价。这定义一个适用于任意排列（包括各对交错出现）的乘积均衡。

拼接理由：各组类型相互独立，记忆不跨组传播，无预算或单位需求约束，效用跨组相加。其他人只按组内类型和历史行动，组外历史在均衡路径上不提供该组未知类型的信息。一个 bidder 的跨组策略偏离只能为各组引入额外随机化或使用与该组类型无关的信息，不能改进该组已证明的最优响应；各组收益求和也不能改进。各组最后一次出售的真实出价对任意信念最优，首轮的同一阈值论证适用。因此无需限制组必须连续出售。这里选定这一乘积均衡族，不声称所有均衡收入相同。

令 m_l(pi) 是第 l 对中被排列实现的方向边权；令

\[
C=\sum_k a_k+q\frac{23}{30}+r\frac13.
\]

在上述均衡族中，两个修正项精确为

\[
L=\sum_{l:m_l(\pi)<0}m_l(\pi),\qquad
A(\pi)=C+2\delta\sum_{l:m_l(\pi)<0}m_l(\pi).
\]

于是

\[
\boxed{R(\pi)=C+\sum_{l=1}^q\psi(m_l(\pi)).}
\tag{9}
\]

凡实现负方向的一对，两个商品必由不同 bidder 获得。其他对可能由不同的人获胜，因此完全不要求第一个赢家包揽 K 件商品。

排序算法：逐对比较 psi(M_jk)、psi(M_kj)，将较大者对应的商品先出售；孤立商品任意放置。由于各对互不相交，不存在有向环冲突，可以直接逐对输出，再插入孤立商品。图以稀疏边表给定时，比较和输出为 O(K)。这是所有排列中、在上述选定均衡族下的精确收入最优算法。

这不是一般图 hardness 结果。图的非平凡性是符号可混合、私人类型多维、前瞻出价和赢家更替真实存在；局限是组间无相互作用，且 n=2。

## 7. 数值核验

脚本 verify_signed_memory.py 采用 200000 个独立类型样本，固定 seed=20260920，遍历四商品实例的全部 24 个排列：

- n=2，alpha_L=1，alpha_W=2；每个 x_ik=5+U_ik。
- M_AB=1，M_BA=-1，M_CD=-1.5，M_DC=-1；其他边零。
- 最优方向 A 在 B 前、D 在 C 前：精确收入 20+46/30+2-3=20.5333333333。
- 最差方向 B 在 A 前、C 在 D 前：精确收入 20+46/30-3-4.5=14.0333333333。
- 两个方向都负的 C/D 对，无论如何排序都会换赢家；A/B 对可通过顺序改变同一人赢两件或两人各赢一件。

代码分别核验：每条样本路径的 (2) 恒等式；正负边下赢家模式；24 个排列的理论收入与抽样均值；相同对内顺序、不同交错方式的逐路径收入不变；用有理数对 3362 个阈值组合精确核验 (6)。代码不是一般均衡求解器，Monte Carlo 也不代替第 4/5 节的证明。

运行：`python verify_signed_memory.py --samples 200000 --output signed_results.json`，需要 numpy。

## 8. 接下来真正要攻克的部分

1. 从匹配图推进到共享节点的最小连通结构，如 A->B->C，边可正可负；先精确描述最后两轮的私人信息 continuation，再分析第一轮，而不是把当期价值直接当报价。
2. 允许弱跨组边时，先检验候选均衡变化，证明策略/租金稳定后才能给收入误差界。对旧策略的 epsilon 最优响应保证不自动推出新均衡收入连续。
3. 任意 n 的 L、A 分解已经成立，但两 bidder 负边下“输 A 必赢 B”的简化在 n>=3 失效。此时输家之间仍竞争 B，必须重新计算 continuation。
4. 本次在可解子类中把修正项精确消去，证明目标有实现可能；任意连通图的近似保证仍是研究任务，不能把它写成已证定理。

## 9. 文献定位提醒

仅引入正负关联不足以声明创新。顺序拍卖中的 complements/substitutes 和 externalities 已有研究，包括：

- Paes Leme, Syrgkanis, Tardos, Sequential Auctions and Externalities, SODA 2012. https://arxiv.org/abs/1108.2452
- Donna and Espin-Sanchez, Complements and substitutes in sequential auctions: the case of water auctions, RAND Journal of Economics, 2018. https://doi.org/10.1111/1756-2171.12221

本稿只确立模型内的结果；需要进一步比较这些文献及其引用链，判断一般有向记忆图、输家状态变化和卖家顺序设计组合下哪些定理具有新颖性。
