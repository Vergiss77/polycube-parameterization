# polycube-simplification Specification

## Purpose

在不改变 polycube 拓扑的前提下，把能对齐的同轴面片吸附到同一高度，减少不同高度（切割面）的数量，并按新高度重新求解体映射，为后续分块提供更简洁的 polycube 结构（Chen et al. 2019 §5）。

## Requirements

### Requirement: Height groups and cut count
系统 SHALL 把同一坐标轴上高度相等（舍入误差以内）的面片视为一个高度组，并以各轴高度组的数量之和作为切割面数，在简化前后分别报告（总数及每个轴的数量）。

#### Scenario: Cut count before simplification
- **WHEN** 对 CLEVIS 的结构提取结果统计切割面数
- **THEN** 报告的总数等于三个轴各自的高度组数之和，且不大于面片数

### Requirement: Order-preserving merges
系统 SHALL 只合并同一坐标轴、且在"低 → 高"有向图中彼此之间没有有向路径的两个高度组。有向图的边有两类，都由较低一端的面片所在组指向较高一端：
- 棱线边：来自平行于该轴的棱线（Chen 2019）；
- 间隔边：来自隔着实体（板的两面）或隔着空隙（相对的两个面）沿该轴正对的两个同轴面片，这是对论文的扩展；
- 宽度边：在任一面片内部沿该轴方向，两侧分别以两个同轴面片为边界的条带，由较低的边界面片指向较高的边界面片，保证条带宽度为正，也是对论文的扩展。

合并后，每条边两端的高度 MUST 保持严格的"低 < 高"。

#### Scenario: Connected groups are not merged
- **GIVEN** 三个同轴高度组 A、B、C，棱线给出 A → B → C
- **WHEN** 运行简化
- **THEN** A 与 C 不会被合并，且结果中 A 的高度 < B 的高度 < C 的高度

#### Scenario: Disconnected groups are merged
- **GIVEN** 两个同轴高度组之间没有任何有向路径，合并后所有相关棱线的畸变都在上限之内
- **WHEN** 运行简化
- **THEN** 两组被合并为一个高度组

#### Scenario: Slab not collapsed
- **GIVEN** 一块板的底面和顶面之间没有任何棱线直接相连（例如两端各有一个凸台的 U 形底板）
- **WHEN** 运行简化
- **THEN** 底面与顶面不会被合并，板的厚度保持为正

#### Scenario: Gap not closed
- **GIVEN** 两个隔着空隙相对的同轴面片（例如 U 形两个凸台的内侧面）
- **WHEN** 运行简化
- **THEN** 两者不会被合并，空隙宽度保持为正

#### Scenario: Strip inside a patch not collapsed
- **GIVEN** 一个面片内部沿某轴的一条窄带，两侧分别以两个同向的同轴面片为边界（例如 CLEVIS 底板顶面上，耳片外侧面与底板端面之间的窄带）
- **WHEN** 运行简化
- **THEN** 两个边界面片不会被合并，窄带宽度保持为正

#### Scenario: Different axes never merged
- **WHEN** 对任一提供的模型运行简化
- **THEN** 每个合并后的高度组内所有面片都属于同一坐标轴

### Requirement: Distortion-bounded greedy simplification
每次合并 SHALL 选择使全部边（棱线边、间隔边和宽度边）弹簧畸变总和 Σ d(l) 增加最少的候选对，并为合并后的组取使该增量最小的高度。其中 d(l) = r̄/r + r/r̄，r̄ 为简化前、r 为当前的两端高度差。任何边的长度比 r/r̄ MUST 保持在 [1/s, s] 内（s 为可配置的最大拉伸比）。没有满足上限的候选对时，简化终止。

#### Scenario: Stretch bound respected
- **WHEN** 以最大拉伸比 s 对任一提供的模型完成简化
- **THEN** 每条边（棱线边、间隔边、宽度边）的高度差之比都位于 [1/s, s] 区间内

#### Scenario: Optimal merged height for a single pair
- **GIVEN** 两个无路径相连的同轴高度组，各自只有一条入边（其余高度固定）
- **WHEN** 计算这一对的合并代价
- **THEN** 得到的合并高度使两条棱线的畸变之和最小，与解析解之差不超过 1e-9 倍包围盒对角线

#### Scenario: Larger bound merges more
- **WHEN** 对 TAILSITTER 分别以较小和较大的最大拉伸比运行简化
- **THEN** 较大拉伸比结果的切割面数不多于较小拉伸比的结果

### Requirement: Inversion-safe rollback
贪心合并结束后，系统 SHALL 只保留按合并顺序的一个最长前缀。判据是：按该前缀的高度重新映射（不解缠）后，翻转四面体数不超过按简化前高度重新映射（不解缠）的翻转数。前缀长度用二分查找确定。统计报告 MUST 给出贪心得到的合并数、保留的合并数，以及两种映射的翻转数。

#### Scenario: Harmful merges are rolled back
- **GIVEN** 只用棱线边建图的 U 形 polycube，贪心会先合并两个顶面，再把底板的底面和顶面合并
- **WHEN** 运行带回退的简化
- **THEN** 只保留合并两个顶面的那一次，底板不被压扁，重新映射后的翻转数为 0

#### Scenario: No rollback when not needed
- **WHEN** 对 I-SPAR 运行简化
- **THEN** 保留的合并数等于贪心得到的合并数

### Requirement: Topology unchanged by simplification
简化 SHALL 只改变面片高度：面片、棱线、角点的数量和连接关系，以及面片的轴与方向，MUST 与简化前完全一致。

#### Scenario: Same structure after simplification
- **WHEN** 对任一提供的模型完成简化
- **THEN** 面片数、角点数、棱线数与简化前相同，拓扑检查仍然通过

#### Scenario: Cut count does not increase
- **WHEN** 对任一提供的模型完成简化
- **THEN** 简化后的切割面数不大于简化前；对 CLEVIS 和 TAILSITTER，切割面数严格减少

### Requirement: Volumetric map with simplified heights
系统 SHALL 在"每个面片的节点在其轴向上的坐标等于简化后的面片高度"的约束下重新求解全部节点位置，最小化相对于原始四面体网格的体积 ARAP 畸变，并报告平均与最大畸变以及翻转四面体数（解缠前后）。

#### Scenario: Heights imposed exactly
- **WHEN** 对任一提供的模型完成简化后的映射
- **THEN** 每个面片内所有节点的轴向坐标与该面片简化后高度之差不超过包围盒对角线的 1e-9 倍

#### Scenario: Inversions on provided models
- **WHEN** 对三个提供的模型完成简化后的映射
- **THEN** I-SPAR 与 CLEVIS 的翻转四面体数为 0；TAILSITTER 的翻转四面体数不超过四面体总数的 0.01%

### Requirement: Deterministic simplification
相同输入与参数下，简化结果（合并顺序、最终高度）SHALL 在同一平台上逐位一致。

#### Scenario: Repeated run identical
- **WHEN** 对同一结构提取结果连续运行两次简化
- **THEN** 两次的合并记录与最终高度完全相同

### Requirement: Simplification stage in the pipeline
流水线 SHALL 在结构提取之后运行简化阶段（默认开启），并支持：
- 通过参数关闭简化；
- 通过 `--until structure` 在简化之前停止；
- 通过命令行或参数文件设置最大拉伸比。

实际使用的简化参数 MUST 写入 `params.json`。

#### Scenario: Default run includes simplification
- **WHEN** 以默认参数完整运行 CLEVIS
- **THEN** 输出目录包含简化阶段的图文件和简化后的 polycube 网格，统计报告包含简化字段

#### Scenario: Simplification disabled
- **WHEN** 关闭简化运行 CLEVIS
- **THEN** 不产生简化阶段的输出，其余输出与启用前的流水线相同，且 `params.json` 中记录简化为关闭

#### Scenario: Stop before simplification
- **WHEN** 以 `--until structure` 运行
- **THEN** 不产生简化阶段的输出

### Requirement: Simplification outputs
简化阶段 SHALL 输出：
- 简化后的面片图（每个面片的简化前后高度、所属高度组，以及按顺序记录的合并：两组编号、合并高度、代价）；
- 简化后的 polycube 四面体网格，提供无损格式和带逐四面体畸变与翻转字段的可视化格式；
- 统计报告中的简化字段：简化前后每轴和总切割面数、贪心与保留的合并次数、各类边的数量、最大拉伸比、实际最大与最小长度比、简化后平均与最大畸变、翻转四面体数、耗时。

所有数值 MUST 合法（非 NaN）。

#### Scenario: Outputs complete
- **WHEN** 对 TAILSITTER 完整运行
- **THEN** 上述文件都存在，统计报告的简化字段齐全且均为合法数值
