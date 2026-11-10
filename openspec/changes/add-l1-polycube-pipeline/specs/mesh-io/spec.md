# Spec Delta

## Purpose

读取并校验输入三角表面网格，导出表面网格与四面体网格（含逐单元/逐面字段）以供查看和后续阶段使用。

## ADDED Requirements

### Requirement: Read OFF triangle meshes
系统 SHALL 读取 OFF 格式的三角表面网格，得到顶点坐标与三角面索引；遇到非三角面、索引越界或文件格式错误时 MUST 报错并说明原因，而不是静默继续。

#### Scenario: Read provided models
- **WHEN** 读取 `models/` 下的 CLEVIS、I-SPAR、TAILSITTER 三个 OFF 文件
- **THEN** 得到的顶点数与面数分别与文件头一致（20840/41684、14767/29530、28326/56668）

#### Scenario: Malformed file
- **WHEN** 读取一个面索引越界或含四边形面的 OFF 文件
- **THEN** 读取失败并返回指明问题所在的错误信息

### Requirement: Validate input surface
系统 SHALL 在流水线开始前校验输入表面为闭合、边流形、无退化三角形，并报告连通分量数、欧拉示性数与亏格；校验失败时 MUST 终止流水线并给出失败项。

#### Scenario: Valid models report topology
- **WHEN** 校验三个提供的模型
- **THEN** 三者均通过，且报告亏格分别为 CLEVIS 2、I-SPAR 0、TAILSITTER 5，连通分量数均为 1

#### Scenario: Open surface rejected
- **WHEN** 校验一个含边界边的表面
- **THEN** 校验失败并报告边界边数量

### Requirement: Export meshes with fields
系统 SHALL 将四面体网格导出为 ParaView 可直接打开的 `.vtu` 文件，可附带逐单元标量字段；SHALL 将三角表面导出为 OFF 或 OBJ，并可附带逐面整数标号字段（以 `.vtu` 表面形式）。

#### Scenario: Tet mesh with per-cell field
- **WHEN** 导出带逐四面体畸变值的四面体网格
- **THEN** 生成的 `.vtu` 在 ParaView 中可打开，单元数与网格一致，并含同名标量字段

#### Scenario: Surface with labels
- **WHEN** 导出带逐面轴向标号的表面
- **THEN** 生成的文件中每个三角面都有对应的整数标号值

### Requirement: Round-trip fidelity
导出再读取的 OFF 表面 SHALL 保持顶点坐标（在双精度十进制往返误差内）与面连接关系不变。

#### Scenario: OFF round trip
- **WHEN** 将读入的模型导出为 OFF 后再次读取
- **THEN** 顶点数、面数、面索引完全相同，坐标差不超过 1e-12 相对误差
