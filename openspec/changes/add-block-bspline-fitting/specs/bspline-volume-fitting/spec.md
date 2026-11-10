# Spec Delta

## Purpose

在 polycube 分块的基础上，为每个块拟合一个三变量 B 样条体，块间 C0 连续，逼近从 polycube 到原模型的体映射；报告拟合误差和缩放 Jacobian，并以 G+Smo MultiPatch XML 与 VTK 采样网格输出多片体（Chen et al. 2019 §6）。

## ADDED Requirements

### Requirement: Block parameterization
每个块的参数域 [0,1]³ SHALL 与该块在 polycube 上的包围盒仿射对应：参数 u、v、w 分别沿 polycube 坐标系的 +x、+y、+z 方向增加。参数点对应的物理位置 MUST 由 polycube 四面体网格到原模型四面体网格的分片线性映射给出（输出单位，包含全局旋转的逆变换，与输入模型位于同一坐标系）。

#### Scenario: Identity-like map
- **GIVEN** 一个本身就是长方体的输入（polycube 映射近似为恒等）
- **WHEN** 对其分块并拟合
- **THEN** 块角点参数 (0,0,0) 与 (1,1,1) 处的样条体位置与长方体对应角点之差不超过包围盒对角线的 1e-6 倍

#### Scenario: Unlocated sample points
- **WHEN** 某个采样点在 polycube 四面体网格中找不到包含它的四面体（例如落在翻转单元附近）
- **THEN** 系统用距离最近的四面体给出位置，不中断运行，并在统计报告中给出这类采样点的数量

### Requirement: Conforming B-spline spaces
每个块的样条体 SHALL 在三个参数方向上使用相同的次数 p（可配置，默认 3）和开放均匀节点向量。每个方向的单元数只由该块在这一方向上所占的 polycube 区间长度决定。因此，共享一个面的两个块在该面上具有相同的节点向量，并共享这个面上的控制点。

#### Scenario: C0 across shared faces
- **WHEN** 对任一提供的模型完成拟合
- **THEN** 对每个内部面，在面上的采样点处分别用两侧块的样条体求值，位置之差不超过包围盒对角线的 1e-12 倍

#### Scenario: Element count from interval length
- **GIVEN** 目标单元尺寸 h
- **WHEN** 一个区间的长度为 L
- **THEN** 该区间上的单元数为 max(1, round(L / h))，所有占用该区间的块在这一方向上的节点向量完全相同

### Requirement: Global least-squares fitting
系统 SHALL 以去重后的控制点为未知量，对所有块的采样点做一次全局最小二乘拟合，使样条体在每个采样参数点处逼近映射给出的物理位置。为保证方程组可解，加一个可配置的小权重控制网光顺项。拟合结果 MUST 报告：
- 全部采样点上的均方根误差和最大误差；
- 位于模型表面的采样点上的均方根误差和最大误差；
- 从样条体边界面采样到原模型表面的最大距离和平均距离。

以上误差都按原模型包围盒对角线归一化。

#### Scenario: Linear map reproduced
- **GIVEN** 一个仿射映射给出的采样（物理位置是参数的线性函数），光顺权重为 0
- **WHEN** 对其拟合
- **THEN** 最大误差不超过包围盒对角线的 1e-9 倍

#### Scenario: Fitting accuracy on provided models
- **WHEN** 以默认参数对三个提供的模型完成拟合
- **THEN** I-SPAR 的最大误差不超过包围盒对角线的 1e-3 倍；CLEVIS 与 TAILSITTER 表面采样点的均方根误差不超过包围盒对角线的 5e-3 倍

### Requirement: Jacobian quality report
系统 SHALL 在每个块的规则评估网格上计算样条体的缩放 Jacobian（Jacobian 行列式除以三列长度之积），报告全体的最小值、平均值、非正值点的个数与比例，以及含非正值点的块数。本阶段 MUST NOT 为消除负 Jacobian 而修改拟合结果。

#### Scenario: Valid on a simple model
- **WHEN** 以默认参数对 I-SPAR 完成拟合
- **THEN** 所有评估点的缩放 Jacobian 都为正

#### Scenario: Reported on complex models
- **WHEN** 对 CLEVIS 和 TAILSITTER 完成拟合
- **THEN** 统计报告给出缩放 Jacobian 的最小值、平均值和非正值比例，且都是合法数值

### Requirement: Deterministic fitting
相同输入与参数下，拟合得到的控制点 SHALL 在同一平台上逐位一致。

#### Scenario: Repeated run identical
- **WHEN** 对 I-SPAR 连续运行两次完整流水线
- **THEN** 两次输出的样条体文件完全相同

### Requirement: G+Smo multipatch output
系统 SHALL 把拟合结果写成 G+Smo 可读取的 MultiPatch XML：每个块一个三变量张量积 B 样条几何（次数、三个方向的节点向量、按 u 最快、w 最慢顺序排列的控制点，输出单位），以及整个多片体的块列表、块间接口（两侧块和对应面）和边界面（块和对应面）。数值 MUST 以足够的有效位数写出，使读回后的控制点与写出前之差不超过 1e-15 倍相对误差。

#### Scenario: Round trip
- **WHEN** 读回写出的 XML
- **THEN** 块数、次数、节点向量和控制点与内存中的结果一致，接口数等于分块的内部面数，边界面数等于分块的边界面数

### Requirement: Sampled volume output
系统 SHALL 输出每个块评估网格上的样条体采样，作为六面体单元组成的可视化网格，带逐单元的块编号和逐点的缩放 Jacobian 字段，可在 ParaView 中查看。

#### Scenario: Sample mesh complete
- **WHEN** 对 CLEVIS 完成拟合
- **THEN** 采样网格的单元数等于各块评估网格单元数之和，每个单元都有块编号，每个点都有缩放 Jacobian

### Requirement: Block and fitting stages in the pipeline
流水线 SHALL 在高度简化之后依次运行分块阶段和拟合阶段（默认运行到拟合），并支持：
- 通过 `--until simplify` 在分块之前停止，通过 `--until blocks` 在拟合之前停止；
- 通过命令行或参数文件设置次数、目标单元尺寸（相对于 polycube 包围盒对角线）、每个单元的采样数和光顺权重，非法取值（次数小于 1、尺寸或采样数不为正、权重为负）MUST 报错。

实际使用的拟合参数 MUST 写入 `params.json`。统计报告 MUST 包含拟合字段：次数、控制点总数（去重后）、采样点数、上述各项误差、缩放 Jacobian 统计、C0 最大偏差、找不到包含四面体的采样点数和耗时，所有数值合法（非 NaN）。

#### Scenario: Default run includes fitting
- **WHEN** 以默认参数完整运行 CLEVIS
- **THEN** 输出目录包含分块结构文件、样条体 XML 和采样网格，统计报告包含分块与拟合字段

#### Scenario: Stop before fitting
- **WHEN** 以 `--until blocks` 运行
- **THEN** 输出分块阶段的文件，不产生拟合阶段的输出

#### Scenario: Invalid parameter
- **WHEN** 以次数 0 运行
- **THEN** 工具以非零状态退出，并说明参数非法
