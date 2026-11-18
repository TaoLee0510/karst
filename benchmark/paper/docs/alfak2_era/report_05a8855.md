# alfak2 05a8855（三种多时间点方法）与 alfakR 在 ALFA-K 文章检验上的对照报告

日期：2026-11-15。alfak2 版本 05a8855（`hpc/alfak2.lock`），alfakR 版本 alfakR_dev2 f93b2d2。

> 用语：输入的各列是不同时间取的细胞，本报告称"时间点"（图里是 "time points"），不称"样本"。
> 统计意义上的样本（bootstrap 等）和 Salehi 的来源系（SA609 等）不在此列。
> 图和表已按此重画（不重新拟合）；改动前的版本在 `results_05a8855/pre_time_points/`，除标签和 auto 的 Fig. 4–6 外数字相同。

- 结果在 `results_05a8855/`；前几轮的结果和报告都在原处：
  - 8327282：`results/`，`docs/report.md`
  - ff13cf4：`results_ff13cf4/`，`docs/report_ff13cf4.md`
  - 8512694：`results_8512694/`，`docs/report_8512694.md`
  - c7f809c：`results_c7f809c/`，只是代价测试（54 个任务），由 05a8855 取代
- 表在 `results_05a8855/analysis/{abm,salehi,summary,fig456}/`，图在 `results_05a8855/figures/`，文中写作 `figures/…`。
  总图是 `figures/summary_comparison.png`；Fig. 4–6 的图是 `figures/fig4_effects.png`、`fig5_wgd.png`、`fig6_dominance.png`。

## 0. 这一轮测的是什么

| commit | 改动 | 对本 benchmark 的影响 |
|---|---|---|
| c7f809c | `fit_sparse_joint(series = …)`：除默认的 `"anchored"` 外，加入 `"pairwise"`（相邻时间点逐对拟合）和 `"state_space"`（时间点之间的适应度可以漂移，`diagnostics$drift`） | 每个"全部时间点"输入多两个拟合；pairwise 要求每一对相邻时间点都有覆盖，Salehi 有 42 个扫荡输入因此不能拟合 |
| 05a8855 | pairwise / state_space 遇到没有覆盖的一对时，去掉时间点直到每一对都有覆盖（`sparse_joint_pair_left_out()`）；只剩首尾两个时退回两时间点拟合 | 42 个输入各去掉 1 个时间点，其中 4 个退回两时间点拟合；anchored、两时间点和完全覆盖的拟合与 c7f809c 逐位相同 |

- alfakR 的拟合与 8327282 那一轮相同，用符号链接复用。
- alfak2 的 5617 个任务全部重跑：
  - 两时间点（ends）与 anchored 全部时间点：3107 个，与 8512694 相同（第 6 节核对）；
  - pairwise 与 state_space 全部时间点：各 1255 个。

## 1. 结论

1. **三种 alfak2 全部时间点拟合几乎在每一项准确度检验上都优于 alfakR 全部时间点拟合**（第 3 节），
   与 8512694 的结论一致。唯一的例外是错分率实验 8 次传代的窗口上的 pairwise（见下）。
2. **在 alfak2 内部，三种方法各有所长，没有一种全面胜出：**
   - **state_space 景观最准**：ABM 配对中位 r 0.66 / 0.72（4 / 8 次传代的窗口），anchored 0.66 / 0.71；
     state_space 更高的输入占 65%。错分率实验上同样（0.77 / 0.83 对 0.76 / 0.82，67–71%）。
     但 ABM 预报略差：第 10 次传代胜过基线 0.84 / 0.89 对 0.86 / 0.91（Euclidean，p ≤ 5e-4）。
     Salehi 上 Wasserstein 预报也更差：0.49 对 0.54（p = 0.007）。
   - **pairwise 在 ABM 上明显最差**：r 0.57 / 0.63，比 anchored 低的输入占 68–70%；
     预报 0.76 / 0.78。错分率实验 8 次传代的窗口上甚至不如 alfakR 全部时间点（0.74 对 0.80）。
   - **pairwise 在 Salehi 新核型上最好**：适应度预测新核型出现的比例 0.66 对 anchored 0.52（p = 4e-9，351 对）；
     预报与 anchored 没有显著差别。
3. **扫荡子集**（anchored 去掉了时间点的 Salehi 轨迹，48 对预报、92 个新核型回归）：
   - 预报胜过亲代：Wasserstein pairwise 0.66、anchored 0.64、state_space 0.63；
     Euclidean 0.48、0.44、0.39。差别都不显著（p ≥ 0.29）。
   - 新核型：pairwise 0.58 对 anchored 0.47（p = 0.03），state_space 0.45。
   - alfakR 在这一子集上仍只有 0.11。
4. **Fig. 4–6 的生物学结论**（第 5 节）：
   - 文章的拟合用本仓库的代码重算，19 条结论全部复现；alfakR dev2 全部时间点复现 15 条，两时间点 14 条。
   - alfak2 各臂只复现 9–11 条（未采用的 auto 13 条）。最大的分歧是 **Fig. 4d：alfak2 给出的结论方向与文章相反**——
     文章说 PDX 中拷贝数改变的适应度效应比体外大（PDX_v +1.23），alfak2 各臂（包括 auto）都是负的（−0.38 到 −1.40），
     因为 alfak2 在体外的 |Δf| 是文章的 4–8 倍，PDX 中则相近。
   - Fig. 5c（WGD+ 与 WGD− 的 Δf 分布不同）：只有文章的拟合显著（p = 0.03）；alfak2 各臂 p 0.10–0.20，alfakR 0.38–0.60。
   - Fig. 6（主导核型随错分率改变）：anchored 0.80、state_space 0.82，与参照 0.80 相同；auto 0.77，pairwise 0.69，两时间点 0.60。
5. **代价**：state_space 最贵。全部时间点拟合的中位时间 ABM 804 秒（anchored 313，pairwise 419，alfakR 3525）；
   Salehi 最长 8.5 小时、13.6 GB（SA609_X10，10 次传代）。所有任务都完成，没有超时或内存不足。
6. **收敛检查不接受的比例**：state_space 在长窗口上最高（ABM w8_all 49%，Salehi all 38%），pairwise 最低（24%，14%）。
7. **auto（按规则在三种方法中选一个，阈值 0.3 / 0.5）没有通过事先登记的判据**（第 9 节），13 条中 3 条不过：
   - 错分率实验 8 次传代的窗口 r 0.804 对 0.815（允许差 0.01）；
   - Salehi Wasserstein 预报 0.506 对 0.541（p = 0.004）；
   - 扫荡 Euclidean 预报 0.417 对 0.438（48 个子代中差 1 个，允许差 0）。

   auto 因此没有被采用，alfak2 包删去了这个选项，默认仍是 anchored。
   auto 在 Salehi 新核型上比 anchored 好（0.607 对 0.521，p = 6e-5）。

## 2. 设置

| 拟合 | 输入 |
|---|---|
| alfak2 2 time points | 窗口或轨迹的首尾两个时间点 |
| alfak2 all time points | 全部时间点，默认 `series = "anchored"` |
| alfak2 pairwise all time points | 全部时间点，`series = "pairwise"` |
| alfak2 state space all time points | 全部时间点，`series = "state_space"` |
| alfak2 auto all time points | 每个输入从上面三个全部时间点拟合中按规则选一个，不另外拟合；评估过，未采用（第 9 节） |
| alfakR 2 time points / all time points | 首尾两个时间点 / 全部时间点（文章的做法） |
| paper alfakR | 文章自己的拟合，只作参照 |

- 数据、窗口、minobs 和参数与 8512694 相同（`docs/report_8512694.md` 第 2 节）。
- 只有两个时间点的输入（ABM 2 次传代的窗口、Salehi 2 次传代的轨迹）上，三种方法是同一个两时间点拟合。
- 准确度只在所有臂都有结果的同一输入上比较（配对）；失败、时间和内存按全部输入计。
- 扫荡子集：anchored 拟合去掉了中间时间点的 Salehi 轨迹（`summary_comparison.R` 里的 "all time points (sweeps)"）。

## 3. 准确度（`figures/summary_comparison.png`，`analysis/summary/summary_accuracy.tsv`）

全部时间点，配对中位数或胜出比例；括号里是比 anchored 好的配对占比（ABM 景观）或配对检验的 p。

| 检验 | anchored | pairwise | state_space | alfakR all |
|---|---|---|---|---|
| ABM 景观 r，4 次传代的窗口（1177 对） | 0.66 | 0.57（32%） | 0.66（65%） | 0.40 |
| ABM 景观 r，8 次传代的窗口（1190 对） | 0.71 | 0.63（30%） | 0.72（65%） | 0.45 |
| ABM rescaled R²，4 / 8 | 0.41 / 0.49 | 0.30 / 0.38 | 0.42 / 0.51 | 0.12 / 0.18 |
| 错分率实验 r，4 / 8（125 对） | 0.76 / 0.82 | 0.69 / 0.74（19% / 10%） | 0.77 / 0.83（67% / 71%） | 0.65 / 0.80 |
| ABM 第 10 次传代胜过基线，Euclidean，4 / 8 | 0.86 / 0.91 | 0.76 / 0.78 | 0.84 / 0.89 | 0.50 / 0.53 |
| ABM 第 10 次传代胜过基线，Wasserstein，4 / 8 | 0.87 / 0.92 | 0.72 / 0.79 | 0.85 / 0.91 | 0.45 / 0.49 |
| Salehi 预报胜过亲代，Wasserstein（283 对） | 0.54 | 0.52（p = 0.24） | 0.49（p = 0.007） | 0.27 |
| Salehi 预报胜过亲代，Euclidean | 0.42 | 0.45（p = 0.11） | 0.42（p = 0.9） | 0.26 |
| Salehi 新核型：适应度预测出现（351 对） | 0.52 | 0.66（p = 4e-9） | 0.55（p = 0.24） | 0.24 |
| 扫荡：预报胜过亲代，Wasserstein（48 对） | 0.64 | 0.66（p = 0.67） | 0.63（p = 0.86） | 0.11 |
| 扫荡：预报胜过亲代，Euclidean | 0.44 | 0.48（p = 0.29） | 0.39（p = 0.43） | 0.11 |
| 扫荡：新核型（92 对） | 0.47 | 0.58（p = 0.03） | 0.45（p = 0.62） | 0.11 |

- 三种 alfak2 对 alfakR 全部时间点的配对检验 p ≤ 0.006，方向都是 alfak2 更好，只有一处例外：
  错分率实验 8 次传代的窗口上，pairwise 比 alfakR 更好的只占 35%（p = 6e-5，alfakR 更好）。
- 两时间点的结果与 8512694 相同（拟合相同，第 6 节）。

## 4. 文章各项检验，按臂（两种设计）

### 4.1 ABM 与错分率实验（`analysis/abm/`）

**Supp. Fig. 2a**：景观 Pearson r 的中位数，minobs 10，按 λ = 0.2 / 0.4 / 0.8 / 1.6。

| 拟合 | 4 次传代的窗口 | 8 次传代的窗口 |
|---|---|---|
| alfak2 2 time points | 0.20 / 0.50 / 0.82 / 0.85 | 0.21 / 0.61 / 0.85 / 0.86 |
| alfak2 all（anchored） | 0.18 / 0.47 / 0.80 / 0.86 | 0.15 / 0.54 / 0.84 / 0.90 |
| alfak2 pairwise | 0.20 / 0.45 / 0.77 / 0.80 | 0.18 / 0.52 / 0.80 / 0.79 |
| alfak2 state space | 0.20 / 0.50 / 0.80 / 0.86 | 0.18 / 0.57 / 0.85 / 0.88 |
| alfakR 2 time points | 0.09 / 0.39 / 0.66 / 0.46 | 0.10 / 0.48 / 0.72 / 0.56 |
| alfakR all time points | 0.12 / 0.38 / 0.65 / 0.63 | 0.15 / 0.48 / 0.75 / 0.74 |
| paper alfakR | 0.09 / 0.38 / 0.66 / 0.59 | 0.11 / 0.46 / 0.76 / 0.73 |

- 频繁核型的 r（`rfq`）在 8 次传代的窗口上：state_space 0.98 / 0.96 / 0.95 / 0.92，与两时间点相近；
  anchored 在 λ = 0.2 时只有 0.85，pairwise 在 λ = 0.4 时 0.84。

**Supp. Fig. 2b、2f**（`supp2f_flows.tsv`，8 次传代的窗口）：频繁核型 R² > 0 的比例在 λ = 0.2 时
anchored 0.58、pairwise 0.81、state_space 0.86、alfakR 全部时间点 0.77、文章 0.83；λ ≥ 0.4 时 alfak2 各臂都 ≥ 0.98。

**Supp. Fig. 2c–e**（`supp2c_example.tsv`、`supp2de_example.tsv`，w_1p6_m_5e-05_rep_74）：
- 2c，8 次传代的窗口，minobs 10，频繁核型真假适应度的 r：anchored 0.955、state_space 0.957、pairwise 0.925、
  alfak2 两时间点 0.891、alfakR 全部时间点 0.831、文章 0.849。2 次传代的窗口：alfak2 0.46，alfakR −0.24，文章 −0.46。
- 2d–e，8 个最常见核型的频率轨迹，中位投影与观测的平均绝对误差：8 次传代的窗口上 alfakR 0.004、
  anchored 0.011、state_space 0.012、pairwise 0.026。alfak2 的投影误差比 alfakR 大。

**Supp. Note 3.2**（`note32_initial_frequencies.tsv`）：初始频率的估计在所有臂上相关系数中位 ≥ 0.9998，
RMSE 中位 0.002–0.003。

**Supp. Fig. 3**（`supp3_misseg_pearson.tsv`，8 次传代的窗口，错分率 5e-5 / 1e-4 / 2e-4 / 4e-4 / 8e-4）：

| 拟合 | r 中位数 |
|---|---|
| alfak2 2 time points | 0.84 / 0.88 / 0.88 / 0.85 / 0.83 |
| alfak2 all（anchored） | 0.83 / 0.85 / 0.84 / 0.77 / 0.70 |
| alfak2 pairwise | 0.80 / 0.77 / 0.77 / 0.68 / 0.62 |
| alfak2 state space | 0.85 / 0.85 / 0.85 / 0.79 / 0.77 |
| alfakR 2 time points | 0.68 / 0.76 / 0.83 / 0.80 / 0.81 |
| alfakR all time points | 0.72 / 0.80 / 0.81 / 0.79 / 0.81 |

  错分率高时，alfak2 的全部时间点拟合都比两时间点差；错分率 8e-4 时三种全部时间点拟合都不如 alfakR（0.70 / 0.62 / 0.77 对 0.81）。

**Supp. Fig. 4**（`supp4b_cv_flags_accuracy.tsv`，8 次传代的窗口）：CV > 0 时景观 R > 0 的比例，λ ≥ 0.4 时
alfak2 各臂都是 1.00；λ = 0.2 时 anchored 0.36、pairwise 0.30、state_space 0.41（CV < 0 时反而更高，0.47 / 0.41 / 0.54），
CV 在 λ = 0.2 时对 alfak2 没有诊断价值；alfakR 全部时间点 0.52 对 0.46，文章 0.10 对 0.21。

**Fig. 1g / Supp. Fig. 5b**（`fig1g_beat_baseline.tsv`，按 CV > 0 选择，第 10 次传代，8 次传代的窗口，λ = 0.4 / 0.8）：
anchored 0.83 / 0.98，state_space 0.85 / 0.96，pairwise 0.71 / 0.81，alfak2 两时间点 0.78 / 0.91，
alfakR 全部时间点 0.70 / 0.79，文章 0.63 / 0.73。

### 4.2 Salehi（`analysis/salehi/`）

| 检验（文章的选择规则） | alfak2 2 | anchored | pairwise | state_space | alfakR 2 | alfakR all | 文章 |
|---|---|---|---|---|---|---|---|
| Fig. 2b CV 分数中位数（只看各自） | 0.86 | 0.85 | 0.88 | 0.92 | 0.05 | 0.47 | 0.48 |
| Fig. 2c 预测夹角，各簇中位数的中位数（°） | 72 | 55 | 52 | 74 | 83 | 80 | 71 |
| Fig. 2d 预测与观测的 Wasserstein 距离，Spearman | 0.25 | 0.33 | 0.57 | 0.56 | 0.45 | 0.70 | 0.69 |
| Fig. 2f 胜过基线，Wasserstein | 0.49 | 0.56 | 0.62 | 0.61 | 0.25 | 0.30 | 0.38 |
| Fig. 2f 胜过基线，Euclidean | 0.36 | 0.50 | 0.53 | 0.51 | 0.25 | 0.24 | 0.36 |
| 拟合数（Fig. 2c–f） | 53 | 52 | 53 | 49 | 44 | 46 | 45 |
| Fig. 3e–f 适应度显著预测新核型，文章规则 | 0.50 | 0.49 | 0.55 | 0.50 | 0.23 | 0.15 | 0.32 |
| Fig. 3e–f 同上，全部拟合 | 0.57 | 0.53 | 0.63 | 0.55 | 0.21 | 0.23 | 0.38 |

- Fig. 2c：所有臂的夹角都显著小于 90°（簇自助法 p ≤ 0.011）；无关谱系 90.4°（p = 0.80），分叉的姐妹谱系 79°。
- Fig. 3 配对（`fig3_paired.tsv`，503 个回归）：pairwise 对 anchored，只有 pairwise 显著 59 个，只有 anchored 显著 10 个（McNemar p = 7.5e-9）；
  state_space 对 anchored 33 对 24（p = 0.29）。
- 按全部 minobs 的配对预报（`paired_forecast.tsv`，401–408 对）：三种 alfak2 对 alfakR 全部时间点 p ≤ 2e-7；
  pairwise 对 anchored 只有 overlap 指标显著（更差，p = 0.002）；state_space 对 anchored 只有 Wasserstein 显著
  （更差，0.47 对 0.51，p = 0.007），与第 3 节一致。

## 5. Fig. 4–6：能否得到文章的生物学结论（`analysis/fig456/conclusions.tsv`）

- **做法**：每个臂用自己的景观，按文章的规则选拟合（xv > 0，每条轨迹取最好的 minobs），重做文章的分析：
  - Fig. 4：200 次非重叠谱系的自助抽样，频繁核型的 44 个单染色体邻居的 Δf，glmmTMB 和 lmer 模型；
  - Fig. 5：SA906 的分段 lmer 和置换 KS 检验（1e4 次）；
  - Fig. 6：每个 clean_id 的最好景观，`find_steady_state()` 在 8 个错分率下的稳态。
- **参照**：
  - Fig. 4 是文章发表的自助抽样结果；
  - Fig. 5、6 是文章的拟合经本仓库同一代码重算（"paper alfakR, every input"）。
    Fig. 6 用重算而不用发表值，因为 alfakR dev2 的 `find_steady_state()` 与文章所用版本不同：
    同样的拟合、同样的核型，稳态最多差 0.18，35 个景观中主导核型在每个错分率下都相同的只有 26 个（`fig6_reproduction.tsv`）。
- **排除**：24 条不可分析的轨迹对所有臂都去掉（`FIG456_EXCLUDE=1`）；文章拟合在全部输入上的复现另算（`*_all_inputs`）。
- **"一致"**：与参照同号，且显著性相同。

| 臂 | 一致的结论（共 19 条） |
|---|---|
| paper alfakR（本代码重算） | 19 |
| paper alfakR，全部输入 | 19 |
| alfakR all time points | 15 |
| alfakR 2 time points | 14 |
| alfak2 auto（未采用） | 13 |
| alfak2 2 time points | 11 |
| alfak2 all（anchored） | 10 |
| alfak2 state space | 10 |
| alfak2 pairwise | 9 |

**Fig. 4d**（`figures/fig4_effects.png`；自助抽样中位数，* 为 95% 区间不含 0）：

| 系数 | 文章 | alfak2 2 | anchored | pairwise | state_space | auto | alfakR 2 | alfakR all |
|---|---|---|---|---|---|---|---|---|
| PDX_v（PDX 对体外，\|Δf\|） | +1.23* | −1.2* | −1.3* | −1.4* | −0.38* | −1.2* | +0.69* | +1.5* |
| PDX_kv（同上，方差） | +1.72* | −0.56 | −0.78* | −1.1* | +0.58 | −0.71* | +1.2* | +1.7* |
| cisplatin_v | +0.92* | +0.36* | +0.24 | +0.26 | +0.28* | +0.35* | +0.42 | +0.94* |

- 原因在体外（`fig4_deltaf.tsv`，|Δf| 中位数）：
  - 体外：alfak2 0.027–0.051（auto 0.024），文章 0.006，alfakR 0.004–0.009；
  - PDX 对照：alfak2 0.013–0.017，文章 0.012，alfakR 0.013–0.021。
  - alfak2 在体外给频繁核型的单染色体邻居很大的适应度差（而且多为正，邻居比频繁核型更差），PDX 中则与文章相当。
- 体外与 PDX 的数据差别：相邻传代的间隔中位 37.5 天对 15 天，每个时间点的细胞中位 602 对 316。
  所以不是"短间隔"的问题；原因待 alfak2 那边查。
- Fig. 4f（相似性模型）：alfak2 各臂在 `pair_typesame_traj` 及其交互项上与文章一致；
  `pair_typeparallel` 及其交互项在 alfak2 中显著，在文章中不显著，算作不一致。

**Fig. 5**（`fig5_tests.tsv`、`fig5_model.tsv`）：

| 臂 | 拟合数 | KS D | 置换 p |
|---|---|---|---|
| paper alfakR | 48 | 0.39 | **0.030** |
| alfak2 2 time points | 49 | 0.36 | 0.20 |
| alfak2 all（anchored） | 49 | 0.36 | 0.13 |
| alfak2 pairwise | 52 | 0.32 | 0.10 |
| alfak2 state space | 36 | 0.40 | 0.15 |
| alfak2 auto（SA906 全部选 anchored，与 anchored 相同） | 49 | 0.36 | 0.13 |
| alfakR 2 time points | 34 | 0.23 | 0.60 |
| alfakR all time points | 59 | 0.23 | 0.38 |

- alfak2 的 D 与文章相近，但都不显著；alfakR dev2 的 D 小得多。
- Fig. 5b 分段模型：各臂的 dk 与文章同号且显著，dk:cut3 也是（只有 alfakR 两时间点不显著）。
  `dk:wgd_strWGD+:cut3TRUE` 在 anchored、pairwise、state_space 和 alfakR 全部时间点中与文章反号；
  `dk:wgd_strWGD+` 在 pairwise 和 state_space 中反号。截距不能跨方法比较（适应度的零点不同）。

**Fig. 6**（`fig6_dominance.tsv`）：主导核型随错分率改变的景观比例：

| 臂 | 景观 | 改变的比例 |
|---|---|---|
| 参照（文章拟合重算） | 35 | 0.80 |
| 文章发表 | 35 | 0.77 |
| alfak2 all（anchored） | 39 | 0.80 |
| alfak2 state space | 38 | 0.82 |
| alfak2 pairwise | 39 | 0.69 |
| alfak2 auto | 39 | 0.77 |
| alfak2 2 time points | 40 | 0.60 |
| alfakR all / 2 time points | 36 / 35 | 0.75 / 0.77 |

  各臂都得到"多数景观的主导核型会改变"的结论。但逐个景观看，主导核型在全部错分率下都与参照相同的：
  alfak2 各臂 0%，alfakR 全部时间点 22%，两时间点 3%（`fig6_agreement.tsv`）。

## 6. 失败、时间、内存（`analysis/summary/summary_performance.tsv`）

| 数据 | 分组 | 方法 | 失败 | 时间中位，秒（90% 分位；最大） | 峰值内存中位（最大），MB |
|---|---|---|---|---|---|
| ABM | 全部时间点 | anchored | 0/2400 | 313（971；2441） | 1802（6582） |
| ABM | 全部时间点 | pairwise | 0/2400 | 419（1146；3257） | 1932（6261） |
| ABM | 全部时间点 | state_space | 0/2400 | 804（3816；18538） | 2371（8538） |
| ABM | 全部时间点 | alfakR | 33/2400 | 3525（8728；25142） | 12577（88804） |
| 错分率 | 全部时间点 | anchored | 0/250 | 310（1237；2589） | 2536（6714） |
| 错分率 | 全部时间点 | pairwise | 0/250 | 550（1355；2869） | 2694（7209） |
| 错分率 | 全部时间点 | state_space | 0/250 | 1045（7806；22514） | 3472（12104） |
| 错分率 | 全部时间点 | alfakR | 0/250 | 4426（18839；30045） | 25734（147156） |
| Salehi | 全部时间点 | anchored | 15/561 | 228（622；1508） | 1587（4624） |
| Salehi | 全部时间点 | pairwise | 15/561 | 244（741；1838） | 1725（5926） |
| Salehi | 全部时间点 | state_space | 15/561 | 663（4156；30720） | 2327（13622） |
| Salehi | 全部时间点 | alfakR | 15/561 | 1586（5111；12465） | 8003（41358） |

- 两时间点的行与 8512694 相同。
- 失败按输入 × minobs 计，已去掉被排除的输入；Salehi 的 15 个都是"minobs 20 下没有频繁核型"，四种方法相同。
- 5617 个任务全部 COMPLETED，没有 TIMEOUT 或 OUT_OF_MEMORY。最长的是 SA609_X10_l_10 的 state_space：8.5 小时，
  c7f809c 代价测试中 13 GB。
- 打分（Salehi 277 条轨迹）每条 ≤ 8 分钟。

## 7. 去掉的时间点（`analysis/summary/summary_dropped.tsv`）

| 方法 | 可分析的 Salehi 多次传代轨迹 | 去掉时间点的轨迹 | 去掉的时间点 | 退回两时间点 |
|---|---|---|---|---|
| anchored | 187 | 53 | 103 | 10 |
| pairwise | 187 | 42 | 42 | 4 |
| state_space | 187 | 42 | 42 | 4 |

- ABM 和错分率实验三种方法都没有去掉时间点。
- 排除的 24 条轨迹对所有方法相同（anchored 的规则）。其中 4 条按 pairwise 的规则其实能拟合，
  但为了所有臂用同样的输入，仍然排除。
- state_space 的 `drift`：ABM 中位 0（最大 0.0025）；Salehi 中位 0，90% 分位 0.010，最大 0.16。
  ABM 上 state_space 几乎不用漂移。

## 8. 与 8512694 相同（`analysis/compare_rounds/results_8512694_summary.tsv`）

05a8855 没有改 anchored 和两时间点的代码路径。`scripts/compare_rounds.R` 逐个核对 alfak2 的 8071 个拟合（两时间点和 anchored 全部时间点）：

- 状态完全相同；景观在数值舍入以内相同（均值最大差 2.2e-16，标准差 1.8e-15）。
- 交叉验证输出、history、预报的重复全部相同。
- 打分差最大 0.039，与上一轮（8512694 对 ff13cf4）的 0.038 一样，是 Spearman 相关里并列值的次序，不是拟合不同。

## 9. auto 臂（`analysis/summary/auto_criterion.tsv`，`tables/auto_picks.tsv`）

**规则**：由多时间点会话在看到任何 E4 结果之前定下。阈值是在 e1–e3 上选的，选择标准也在打分前登记：
过 e1 关和 e3 中性检验的规则里，取 e2 比值最低的，差 0.01 以内取较简单的。
暂定的 0.3 / 0.8 规则 e2 比值是 0.822，定下的 0.3 / 0.5 是 0.810。

- `ends_major`：中间各时间点的细胞落在首尾两个时间点主核型上的比例，取最小值；
- `weight`：anchored 拟合的 `intermediate_weight`；
- `ends_major` ≥ 0.3 → anchored；否则 `weight` 为 NA（中间没有剩下的时间点）或 ≥ 0.5 → pairwise；否则 state_space。

**做法**：auto 不另外拟合。

1. `scripts/auto_picks.R 0.3 0.5` 给每个多时间点输入选方法；
2. `scripts/compose_auto.R` 把选中的拟合目录链接为方法 `alfak2_auto`，共 1255 个链接，没有缺失；
3. 与其他臂一样打分、汇总，在同一批输入上和 anchored 配对比较。

**选择**：

| 数据 | anchored | pairwise | state_space |
|---|---|---|---|
| ABM（w4_all、w8_all） | 800 | 0 | 0 |
| 错分率实验 | 157 | 93 | 0 |
| Salehi | 80 | 102 | 23 |

- 28 个输入的 `weight` 为 NA。
- ABM 的中间时间点都落在首尾主核型上，所以全部选 anchored，ABM 上 auto 与 anchored 逐位相同。
- 36 条 SA906 轨迹也全部选 anchored，所以 Fig. 5 的 auto 与 anchored 相同。

**判据**（多时间点会话在 E4 之前登记）：

| 判据 | 配对 | auto | anchored | 差 | 允许 | 结果 |
|---|---|---|---|---|---|---|
| ABM 景观 r，4 次传代的窗口 | 1177 | 0.660 | 0.660 | 0 | −0.01 | 过 |
| ABM 景观 r，8 次传代的窗口 | 1190 | 0.709 | 0.709 | 0 | −0.01 | 过 |
| 错分率实验景观 r，4 次传代的窗口 | 125 | 0.754 | 0.764 | −0.010 | −0.01 | 过（正好在边界） |
| 错分率实验景观 r，8 次传代的窗口 | 125 | 0.804 | 0.815 | −0.011 | −0.01 | **不过** |
| ABM 第 10 次传代胜过基线（Euclidean / Wasserstein，4 / 8） | 1177 / 1190 | 与 anchored 相同 | | 0 | −0.02 | 过（4 条） |
| Salehi 预报胜过亲代，Euclidean | 283 | 0.410 | 0.422 | −0.012 | −0.02 | 过 |
| Salehi 预报胜过亲代，Wasserstein | 283 | 0.506 | 0.541 | −0.035 | −0.02 | **不过** |
| 扫荡：预报胜过亲代，Euclidean | 48 | 0.417 | 0.438 | −0.021 | 0 | **不过** |
| 扫荡：预报胜过亲代，Wasserstein | 48 | 0.639 | 0.639 | 0 | 0 | 过 |
| Salehi 新核型：适应度预测出现 | 351 | 0.607 | 0.521 | +0.086 | −0.02 | 过 |

**结论：auto 不通过**（13 条中 3 条不过），**评估过，未采用**。

- 用户的规则是：判据通过才实现 auto。多时间点会话因此从 alfak2 包中删去 `series = "auto"`，默认仍是 anchored：
  - 4ff8906 加入了 auto；
  - 94cef55（"Call the columns time points, and drop the automatic choice"）删去了它。

  从 94cef55 起，alfak2 没有 auto 选项。94cef55 的 anchored、pairwise、state_space 代码与 05a8855 相同，只改了文字。
- 规则和这里的结果作为负结果记在 alfak2 的 `validation/sparse_joint/series_methods/README.md` 和方法文档（tex）里。
- 不需要新的一轮：删去 auto 不改 anchored、pairwise、state_space 的拟合，改用语的 commit 只改文字。

- **错分率实验，8 次传代的窗口**：差 −0.011 是两位舍入后的中位数相减，只比允许值多 0.001。
  配对看是真实的差：auto 更差的输入占 32%，更好的占 6%（p = 4e-5）。原因是 93 个选 pairwise 的输入，
  pairwise 在错分率实验上本来就最弱（第 3 节）。
- **Salehi，Wasserstein**：更差的配对占 4.6%，更好的占 1.4%（p = 0.004），也是真实的差。
  auto 的 0.506 比 pairwise（0.52）和 state_space（0.49）单独用时都不占优。
- **扫荡，Euclidean**：差 −0.021 正好是 48 个子代中的 1 个（p = 0.85）。不过只是因为这一条的允许值是 0。
- **好的一面**：Salehi 新核型 auto 0.607，比 anchored 0.521 好（p = 6e-5），扫荡子集上 0.554 对 0.47（p = 0.08）。
  这来自 pairwise 的选择，与第 3 节 pairwise 在新核型上最好一致。

**Fig. 4–6 的 auto 臂**（第 5 节各表；重跑 21953562、21953563）：

- 复现文章 19 条结论中的 13 条，是 alfak2 各臂中最多的（其他 9–11 条）。
  比 anchored（10 条）多的 3 条是 Fig. 4d 的两条顺铂系数（cisplatin_v +0.35*、cisplatin_kv +0.62*，与文章同号且显著）
  和 Fig. 4f 的截距（不显著，与文章一致）。
- Fig. 4d 的 PDX 对体外仍与文章相反（PDX_v −1.22*），体外 |Δf| 中位 0.024，是文章的 4 倍。
- Fig. 5：与 anchored 相同。
- Fig. 6：主导核型改变的景观 0.77（30 / 39），参照 0.80；逐个景观与参照一致的仍是 0 / 31。
- 这些不在登记的判据里，不改变"不通过"的结论。

## 10. 收敛检查与 CV 分数

**收敛检查不接受的比例**（minobs 10，`tables/{abm,salehi}_runs.Rds`）：

| 数据 | anchored | pairwise | state_space | 两时间点 |
|---|---|---|---|---|
| ABM w4_all / w8_all | 21% / 41% | 20% / 24% | 23% / 49% | w4 20%，w8 36% |
| 错分率 w4_all / w8_all | 11% / 25% | 16% / 21% | 14% / 31% | w4 10%，w8 19% |
| Salehi all | 20% | 14% | 38% | 17% |

**CV 分数**（`analysis/summary/summary_cv.tsv`）：R²R 相对于各方法自己的全数据估计，不能跨方法比较，只报中位数和大于 0 的比例。

| 数据 | anchored | pairwise | state_space | alfakR all |
|---|---|---|---|---|
| ABM 全部时间点 | 0.56，72% | 0.12，58% | 0.51，73% | −0.01，43% |
| 错分率 全部时间点 | 0.89，96% | 0.52，93% | 0.84，98% | 0.07，78% |
| Salehi 全部时间点 | 0.80，77% | 0.84，82% | 0.86，77% | 0.04，56% |

## 11. 运行记录与决定

| 作业 | 内容 |
|---|---|
| 21945965、21945966 | c7f809c 编译和代价测试（54 个任务），由 05a8855 取代 |
| 21946021 | 编译 alfak2 05a8855 |
| 21946022–21946042（偶数） | ABM、错分率、Salehi 的各内存档，共 5617 个任务；各带一个加时监控（id + 1） |
| 21951672 | 编译 glmmTMB 1.1.11 到 `lib/extra_r`（镜像里没有） |
| 21951796、21952078 | Salehi 打分 277 条轨迹；SA609_X10 的 state_space 拟合在打分后才结束，单独重打 |
| 21952074 | ABM 汇总、Supp. Fig. 2c–e、k 集合核对、与 8512694 的比较 |
| 21952127 | Salehi 汇总、Fig. 5 |
| 21952128–21952131 | Fig. 4（7 臂 × 10 批，加文章拟合全部输入 10 批）、Fig. 6（7 臂加 1） |
| 21952227、21952262 | 总图、资源、CV；Fig. 4–6 汇总 |
| 21953043、21953044 | auto：`auto_picks.R 0.3 0.5`，`compose_auto.R` |
| 21953046 | auto 的 Salehi 打分（277 条轨迹），加时监控 21953049 |
| 21953047、21953048 | Fig. 4、6 的 auto 臂，第一次：读到的是还没有 auto 行的打分表（collect_salehi 在它们之后才跑），作废 |
| 21953051 | 全部汇总（含 auto_criterion、Fig. 5、Fig. 4–6 汇总、与 8512694 的比较） |
| 21953562、21953563 | Fig. 4、6 的 auto 臂重跑，加时监控 21953564 |
| 21953565 | 各轮的图和分析表按"时间点"用语重画，不重新拟合：8512694、05a8855 重画（旧的留在 `pre_time_points/`）；ff13cf4 用现在的代码重画后内容有变，原图恢复，重画版放在 `redrawn_2026-11-15/`；8327282 缺少输入覆盖表，没有重画 |
| 21953580 | Fig. 5 的表（`fig5_tests*.tsv`）按新标签重写，数字不变 |

决定（多时间点会话授权自动运行，按我的建议决定）：

- 记录轮是 05a8855；c7f809c 的 pair 检查会挡住 42 个 Salehi 扫荡输入，只做代价测试。
- 没有单独的测试轮：c7f809c 代价测试已覆盖两种新方法，05a8855 只改 pairwise / state_space 的时间点选择。
- 内存档：pairwise 16 GB、state_space 32 GB，每个任务 12 小时，加时监控每小时检查。
- Fig. 6 的参照改为文章拟合的重算（理由见第 5 节），发表值保留为一行。
- Fig. 4 选最好 minobs 时，xv 并列取最小的 minobs（文章的拟合没有并列，结果不变）。
- auto 的阈值由多时间点会话定（0.3 / 0.5），在给 auto 打分之前发来；在那之前没有向它透露任何 E4 准确度数字。
