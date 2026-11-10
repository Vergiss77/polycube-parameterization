# tet-meshing Specification

## Purpose

把经过校验的闭合三角表面转为四面体网格，作为 polycube 体变形的计算域，并保证边界与输入表面一致、单元质量可控。

## Requirements

### Requirement: Boundary-preserving tetrahedralization
系统 SHALL 从闭合三角表面生成四面体网格，其边界三角形集合 MUST 与输入表面完全一致（不插入边界 Steiner 点），并提供每个边界三角形到输入面索引的对应关系。

#### Scenario: Boundary matches input
- **WHEN** 对任一提供的模型进行四面体化
- **THEN** 四面体网格边界面数等于输入面数，且每个边界面的三个顶点坐标与对应输入面一致

### Requirement: Valid element orientation
生成的每个四面体 SHALL 具有严格正的有向体积，所有四面体体积之和 MUST 等于输入表面围成的体积（相对误差不超过 1e-9）。

#### Scenario: Positive volumes
- **WHEN** 对三个提供的模型四面体化
- **THEN** 不存在体积小于等于零的四面体，且总体积与表面围成体积一致

### Requirement: Controllable element quality and size
系统 SHALL 接受单元质量（半径-边比上限）与最大单元体积参数，未指定时使用由输入平均边长推导的默认值，使内部单元尺寸与表面边长同量级。

#### Scenario: Default parameters
- **WHEN** 不指定质量与尺寸参数对 CLEVIS 四面体化
- **THEN** 成功生成网格，并报告四面体数、节点数与最小二面角

#### Scenario: Explicit size limit
- **WHEN** 指定更小的最大单元体积
- **THEN** 生成的四面体数多于默认参数的结果

### Requirement: Tetrahedralization failure reporting
四面体化失败时，系统 MUST 终止流水线并报告失败原因，不得输出部分结果作为成功结果。

#### Scenario: Self-intersecting input
- **WHEN** 输入表面存在自相交导致四面体化失败
- **THEN** 流水线以非零状态退出并报告四面体化错误
