# Spec Delta

## Purpose

按 Huang 2014 的 ℓ1 方法把四面体网格变形为边界法向轴对齐的低畸变 polycube 形状，并提供复杂度、特征边与收敛方面的控制参数。

## ADDED Requirements

### Requirement: Polycube error metric
系统 SHALL 定义并报告归一化 polycube 误差：所有边界三角形 `面积 × (‖n‖₁ − 1)` 之和除以边界总面积，其中 n 为单位法向。该值对完全轴对齐的边界为 0。

#### Scenario: Already axis-aligned input
- **WHEN** 对未变形的 I-SPAR 计算 polycube 误差
- **THEN** 报告值小于 1e-6

#### Scenario: Curved input
- **WHEN** 对未变形的 TAILSITTER 计算 polycube 误差
- **THEN** 报告值大于 0.05

### Requirement: l1 polycube deformation
系统 SHALL 在保持边界总面积不变的约束下，联合最小化法向 ℓ1 能量与体积 ARAP 畸变能量，把四面体网格变形为 polycube 形状；运行结束时 MUST 报告最终 polycube 误差、平均 ARAP 畸变、迭代次数与耗时。

#### Scenario: Converge on provided models
- **WHEN** 以默认参数对 I-SPAR 和 CLEVIS 运行变形
- **THEN** 每个模型的最终 polycube 误差小于 1e-3，且边界总面积与输入的相对差小于 1e-6

#### Scenario: TAILSITTER converges within relaxed threshold
- **WHEN** 以默认参数对 TAILSITTER 运行变形
- **THEN** 最终 polycube 误差小于 5e-3，边界总面积与输入的相对差小于 1e-6；若未达到 1e-3 阈值，统计报告中 `reached_threshold` 为 false，日志给出警告（已知限制：0.5 mm 台阶和小孔附近在 20 个阶段内无法完全压平）

#### Scenario: Near-identity for aligned input
- **WHEN** 以默认参数对 I-SPAR 运行变形
- **THEN** 平均 ARAP 畸变小于 1e-3

### Requirement: Global orientation search
系统 SHALL 在变形过程中同时求解一个全局旋转，使结果与输入摆放方向无关；该功能 SHALL 可通过参数关闭，关闭时结果对齐世界坐标轴。

#### Scenario: Diagonal arms get aligned
- **WHEN** 以默认参数对 TAILSITTER（四根臂沿 45° 对角方向）运行变形
- **THEN** 报告的全局旋转绕 Z 轴的角度与 45° 的差（模 90°）小于 5°

#### Scenario: Orientation independence
- **WHEN** 先将 CLEVIS 绕任意轴旋转 30° 后运行变形
- **THEN** 最终 polycube 误差小于 1e-3，且面片数量与未旋转输入的结果相同

#### Scenario: Rotation disabled
- **WHEN** 关闭全局旋转对 CLEVIS 运行变形
- **THEN** 报告的全局旋转为单位矩阵

### Requirement: Complexity control
系统 SHALL 提供复杂度权重 β，惩罚相邻边界三角形的法向差；在其他参数相同时，较大的 β MUST 不产生更多的 polycube 角点。

#### Scenario: Larger beta simplifies
- **WHEN** 对 TAILSITTER 分别以较小与较大的 β 运行变形并提取结构
- **THEN** 较大 β 结果的 polycube 角点数不多于较小 β 的结果

### Requirement: Foldover prevention
变形过程中 SHALL 保持任意相邻边界三角形的单位法向点积严格大于 −1（不出现相邻面法向相反的翻折）；若因此无法继续下降，系统 MUST 记录该情况及涉及的三角形数量。

#### Scenario: No opposite adjacent normals
- **WHEN** 对任一提供的模型完成变形
- **THEN** 所有相邻边界三角形法向点积均大于 −1 + 1e-6

#### Scenario: Barrier stall is reported
- **WHEN** 优化因防翻折约束而步长持续过小
- **THEN** 日志中记录停滞事件及涉及的细长三角形数量

### Requirement: Optional feature-edge guidance
系统 SHALL 接受用户指定的边集合，将其引导为直角（polycube 棱）或平坦（同一面片内部），并以可调权重加入优化；未指定时该项不起作用。

#### Scenario: No guidance by default
- **WHEN** 不提供特征边运行变形
- **THEN** 结果与特征边权重为 0 时相同

### Requirement: Scheduled convergence
系统 SHALL 使用 ℓ1 权重 α 递增、平滑参数 ε 递减的分阶段调度；每阶段在 polycube 误差停止明显下降时结束，整体在误差低于阈值或达到最大阶段数时停止。阈值、初值、倍率与最大阶段数 SHALL 可配置，默认值为 Huang 2014 的取值。

#### Scenario: Stop on threshold
- **WHEN** 某阶段结束时 polycube 误差低于 1e-3
- **THEN** 优化停止，不再进入下一阶段

#### Scenario: Stop on stage limit
- **WHEN** 达到最大阶段数（默认 20）仍未低于阈值
- **THEN** 优化停止并报告未达阈值的警告及最终误差

### Requirement: Deterministic results
在相同输入、参数与平台下，变形结果 SHALL 完全可重复（不使用随机数或依赖线程调度的归约）。跨平台时浮点差异可能累积，因此只要求两平台结果都满足收敛阈值且提取出的面片数相同。

#### Scenario: Repeat run
- **WHEN** 在同一平台上对同一模型以相同参数运行两次
- **THEN** 两次输出的节点坐标逐位相同

#### Scenario: Cross-platform consistency
- **WHEN** 分别在 Windows 与 macOS 上以默认参数处理同一提供的模型
- **THEN** 两者最终 polycube 误差均小于 1e-3，且提取出的面片数相同
