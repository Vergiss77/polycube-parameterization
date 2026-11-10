# Proposal

## Why

上一个 change 得到的 polycube 中，每个面片的高度都是最终映射独立解出来的，几乎互不相等：CLEVIS 24 个面片有 24 个不同高度，TAILSITTER 59 个面片有 59 个。后续分块要沿面片所在平面切开 polycube，不同高度越多，切割面就越多，块也越多、越碎。Chen et al. 2019 §5 的有向图简化，正是在不改变 polycube 拓扑的前提下，把能对齐的同轴面片吸附到同一高度，减少切割面，为下一个 change（分块与 B 样条拟合）打基础。

## What Changes

- 新增高度简化（Chen 2019 §5.1–§5.4）：
  - 以面片高度组为顶点、平行于轴的棱线为有向边（低 → 高）建立有向图；
  - 只合并同轴、且彼此之间没有有向路径的顶点对，保证每条棱线两端的高度次序不变；
  - 每对候选的代价用棱线长度变化的弹簧畸变 `d = r̄/r + r/r̄` 衡量，并受最大畸变上限约束；
  - 贪心地每次合并代价最小的一对，直到没有可行的候选。
- 新增按简化后的高度重新求解体映射（Chen §5.5.1）：复用已有的固定高度 ARAP 最终映射，得到满足新高度的 polycube，并报告畸变和翻转四面体数。
- 流水线新增 `simplify` 阶段，位于结构提取之后，默认开启：
  - 可以关闭，也可以用 `--until structure` 停在简化之前；
  - 输出简化前后的高度、合并记录和简化后的 polycube 网格；
  - 统计报告新增简化相关字段。
- 不改变面片、角点、棱线的数量与连接关系。TAILSITTER 的 0.5 mm 台阶这类有棱线相连的上下层，按 Chen 的约束不能合并，本 change 不处理（需要改拓扑，视结果另开 change）。

## Capabilities

### New Capabilities
- `polycube-simplification`: polycube 结构的高度简化：有向图建立、候选对筛选、畸变代价与上限、贪心合并，按新高度重新求解体映射，以及该阶段在流水线中的开关、参数和输出。

### Modified Capabilities
<!-- 上一个 change（add-l1-polycube-pipeline）尚未归档，openspec/specs/ 下还没有主 spec。
     流水线新增阶段相关的行为写在 polycube-simplification 中，不修改 pipeline-cli。 -->

## Impact

- 新增代码：`src/pcp/structure/height_simplify.{h,cpp}`（有向图与贪心合并），单元测试和模型测试。
- 修改代码：
  - `src/pcp/pipeline/`：新增 `Stage::Simplify`、参数段 `simplify`、输出文件和统计字段；
  - `apps/pcp_cli/`：新增简化相关的命令行选项；
  - `final_map` 的固定高度模式会在这里第一次真正被流水线使用。
- 不新增第三方依赖。Chen 用 Ipopt 求解单对代价，这里每对代价只是一个一维凸问题，直接实现即可。
- 运行时间：图很小（最多几百个顶点），简化本身耗时可以忽略；主要开销是再做一次最终映射（TAILSITTER 约 2–3 分钟）。
- 不在本 change 范围内：
  - 改变拓扑的台阶去除；
  - 分块（切割、六面体网格或块提取，即 Chen §5.5.2）；
  - 曲面拟合与 Coons 体插值（Chen §6），以及三变量 B 样条拟合。这些放在下一个 change。
