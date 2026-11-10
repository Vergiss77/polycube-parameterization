# Spec Delta

## Purpose

沿切割面把轴对齐的 polycube 分成相互之间面对面完全贴合的长方体块，作为逐块参数化和 B 样条拟合的分区（Chen et al. 2019 §5.5.2）。

## ADDED Requirements

### Requirement: Cut planes from height groups
系统 SHALL 以 polycube 每个坐标轴上的高度组（同轴、等高的面片集合）的高度作为该轴的切割坐标。切割坐标按升序排列，相邻两个切割坐标之间是一个区间。简化启用时使用简化后的高度组，关闭时使用简化前的高度组。

#### Scenario: Cut count matches simplification
- **WHEN** 对 CLEVIS 完整运行到分块阶段
- **THEN** 每个轴的切割坐标数等于高度简化报告的该轴简化后切割面数

#### Scenario: Simplification disabled
- **WHEN** 关闭简化并运行到分块阶段
- **THEN** 每个轴的切割坐标数等于简化前该轴的高度组数

### Requirement: Blocks cover the polycube exactly
系统 SHALL 把三个轴的区间组成的网格中、位于 polycube 内部的每个格子作为一个块。所有块的并集 MUST 等于 polycube：块体积之和与 polycube 体积之差不超过 polycube 体积的 1e-9 倍，块之间没有重叠，每个块的三个边长都为正。

#### Scenario: L-shaped polycube
- **GIVEN** 由三个单位立方体组成的 L 形 polycube（每轴 3 个或 2 个切割坐标）
- **WHEN** 运行分块
- **THEN** 得到 3 个块，体积之和为 3

#### Scenario: Volume preserved on provided models
- **WHEN** 对三个提供的模型运行分块
- **THEN** 块体积之和与 polycube 四面体网格的总有向体积之差不超过 1e-9 倍

### Requirement: Conforming block faces
每个块的每个面 SHALL 恰好属于以下两种情况之一：
- 与另一个块的一个面完全重合（内部面）；
- 整个面位于 polycube 表面上（边界面），并属于一个确定的面片。

块之间 MUST 不存在 T 型连接（一个块的面只与另一个块的面的一部分重合）。

#### Scenario: Every face classified
- **WHEN** 对任一提供的模型运行分块
- **THEN** 每个块的 6 个面都被标记为内部面（并记录相邻块及其对应面）或边界面（并记录所属面片），内部面的相邻关系是对称的

#### Scenario: U-shaped polycube
- **GIVEN** 两端各有一个凸台的 U 形 polycube
- **WHEN** 运行分块
- **THEN** 两个凸台之间空隙处的格子不是块，凸台与底板之间的面是内部面

### Requirement: Deterministic block numbering
块的编号 SHALL 只由网格位置决定：按 (z, y, x) 区间编号的字典序排列。相同输入下，块编号、包围盒和相邻关系 MUST 逐位一致。

#### Scenario: Repeated run identical
- **WHEN** 对同一 polycube 连续运行两次分块
- **THEN** 两次的块列表（编号、包围盒、各面的类型与相邻关系）完全相同

### Requirement: Block decomposition outputs
分块阶段 SHALL 输出：
- 分块结构文件：每个轴的切割坐标（输出单位），以及每个块的网格位置、包围盒（polycube 坐标）、体积和 6 个面的类型与相邻块或所属面片；
- 带逐四面体块编号字段（按四面体在 polycube 上的重心所在的块；四面体网格不沿切割面剖分，跨越切割面的四面体只归入一个块）的四面体网格，分别在 polycube 上和原模型上各一份，可在 ParaView 中查看；
- 统计报告中的分块字段：块数、每轴切割坐标数、网格格子总数、内部面数和边界面数、耗时。

#### Scenario: Outputs complete
- **WHEN** 对 TAILSITTER 运行到分块阶段
- **THEN** 上述文件都存在，每个四面体都有块编号，统计报告的分块字段齐全且均为合法数值
