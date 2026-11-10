# polycube-structure Specification

## Purpose

从 ℓ1 变形结果中提取拓扑合法的 polycube 结构（面片标号、面片图、棱线与角点、面片高度），并计算与该结构严格一致的体内 polycube 映射，供后续结构简化与分块使用。

## Requirements

### Requirement: Axis labeling of boundary triangles
系统 SHALL 为每个边界三角形赋予六个方向标号之一（±X、±Y、±Z），初始标号取变形后（含全局旋转）单位法向的绝对值最大分量及其符号。

#### Scenario: Aligned input labels
- **WHEN** 对 I-SPAR 的变形结果进行标号
- **THEN** 每个三角形的标号与其原始法向所指的坐标轴方向一致

### Requirement: Label cleanup
系统 SHALL 清理标号中的拓扑缺陷：消除面片边界上的锯齿三角形，并迭代地将度为 1 或 2 的面片并入相邻面片（度为 2 时并入共享边界最长的邻居），直到不存在此类面片。

#### Scenario: Zigzag removed
- **WHEN** 某三角形的标号与其三个邻居中的两个不同、且那两个邻居标号相同
- **THEN** 清理后该三角形取那两个邻居的标号

#### Scenario: No low-degree patches remain
- **WHEN** 对任一提供的模型完成清理
- **THEN** 每个面片相邻的其他面片数都不少于 3

### Requirement: Valid polycube topology
清理后的面片划分 SHALL 满足：每个面片是边连通的三角形集合；相邻面片的标号属于不同坐标轴；面片划分后的表面欧拉示性数与输入一致。不满足时系统 MUST 报告违反项而不是继续输出。

#### Scenario: Topology checks pass
- **WHEN** 对三个提供的模型提取结构
- **THEN** 三项检查均通过

#### Scenario: Known patch count
- **WHEN** 对 I-SPAR 提取结构
- **THEN** 得到 14 个面片（工字形截面 12 个侧面加 2 个端面）

### Requirement: Patch graph export
系统 SHALL 输出机器可读的面片图：每个面片的编号、轴向、符号、三角形列表、高度（该轴坐标）与面积；面片间的邻接关系；棱线（相邻面片之间的公共边链，含其平行轴与两端连接的面片）；角点（三条及以上棱线的交点）。

#### Scenario: Patch graph completeness
- **WHEN** 导出 I-SPAR 的面片图
- **THEN** 文件中含 14 个面片，每条棱线都记录了平行轴及两端相连的同轴面片，角点数为 24

### Requirement: Final polycube volumetric map
系统 SHALL 在"同一面片的边界节点在其轴向上坐标相等"的线性约束下，以低畸变为目标重新求解全部节点位置，得到严格的 polycube 形状及原模型到它的逐四面体线性映射，并报告平均与最大 ARAP 畸变及翻转（有向体积非正）的四面体数量。

#### Scenario: Exact planarity
- **WHEN** 对任一提供的模型计算最终映射
- **THEN** 每个面片内所有节点的轴向坐标与该面片高度之差不超过包围盒对角线的 1e-9 倍

#### Scenario: No inverted tetrahedra on provided models
- **WHEN** 对 I-SPAR 和 CLEVIS 计算最终映射
- **THEN** 报告的翻转四面体数量为 0

#### Scenario: TAILSITTER inversions bounded
- **WHEN** 对 TAILSITTER 计算最终映射
- **THEN** 报告的翻转四面体数量不超过四面体总数的 0.01%，并在统计报告中给出解缠前后的翻转数（已知限制：凹角孔洞附近的局部解缠不保证完全消除翻转）

### Requirement: Map evaluation between domains
系统 SHALL 支持对 polycube 空间中的任意点，找到包含它的四面体并用重心坐标映射回原空间；点落在 polycube 外部时 MUST 明确返回"不在域内"。

#### Scenario: Corner round trip
- **WHEN** 把原模型任一节点经映射得到的 polycube 坐标再映射回原空间
- **THEN** 得到的位置与原节点坐标之差不超过包围盒对角线的 1e-9 倍

#### Scenario: Outside query
- **WHEN** 查询 polycube 包围盒之外的点
- **THEN** 返回不在域内的结果
