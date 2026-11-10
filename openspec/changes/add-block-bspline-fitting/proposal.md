# Proposal

## Why

前两个 change 得到了高度简化后的 polycube，以及从 polycube 到原模型的四面体映射，但还没有项目最终要的东西：三变量 B 样条多片体。Chen et al. 2019 的后半部分（§5.5.2 与 §6）沿切割面把 polycube 分成长方体块，每块对应参数域 [0,1]³，再拟合样条体。本 change 完成这一步，让流水线第一次端到端输出可用于等几何分析的多片体。

## What Changes

- **分块**：
  - 用每个轴上简化后的高度组作为切割面，在 polycube 的包围盒上构造网格，保留落在 polycube 内部的格子作为块；
  - 全局网格保证相邻块的公共面完全重合（不存在 T 型连接）；
  - 记录每块的包围盒、相邻关系（共享面）和位于模型表面的面。
- **逐块参数化与采样**：
  - 每块的参数域 [0,1]³ 与块的 polycube 包围盒仿射对应；
  - 在参数域上取规则采样点，经 polycube 四面体网格定位，再按重心坐标插值得到原模型上的物理位置。
- **B 样条拟合**：
  - 每块一个三变量 B 样条体，默认 3 次，节点向量为开放均匀；
  - 每个轴向区间的单元数只由区间长度决定，同一区间上所有块的节点向量相同，公共面上的控制点直接共享，因此块间 C0 连续；
  - 所有块一起做全局最小二乘，未知量是去重后的控制点；加一个小权重的控制网光顺项，保证方程组正定。
- **质量与误差报告**：
  - 在采样网格上计算缩放 Jacobian，报告最小值和负值比例；
  - 报告拟合误差（全部采样点和表面采样点），以及公共面两侧的 C0 偏差；
  - 本 change 不修复负 Jacobian（控制点质量优化留到下一个 change）。
- **输出**：
  - G+Smo MultiPatch XML（各块的样条体、块间接口、边界面）；
  - 带缩放 Jacobian 与块编号的 VTK 采样网格；
  - 分块结构 JSON 与带块编号的四面体网格 VTK。
- **流水线**：
  - 新增 `blocks` 与 `fit` 两个阶段，位于 `simplify` 之后，默认运行到 `fit`；
  - 关闭简化时，直接用简化前的 polycube 分块；
  - 用 `--until simplify` 可以在分块之前停止；
  - `pcp_cli` 新增拟合参数选项（次数、目标单元尺寸）。

## Capabilities

### New Capabilities
- `polycube-block-decomposition`：沿切割面把 polycube 分成长方体块，提供块的包围盒、相邻关系和表面面，以及块编号的可视化输出。
- `bspline-volume-fitting`：逐块参数化与采样、共享控制点的全局最小二乘拟合、缩放 Jacobian 与误差报告、G+Smo XML 与 VTK 采样输出，以及 `blocks` / `fit` 两个流水线阶段的开关、参数和统计字段。

### Modified Capabilities
<!-- 流水线新增阶段相关的行为写在 bspline-volume-fitting 中，沿用上一个 change 的做法，不修改 pipeline-cli。 -->

## Impact

- **新增代码**：
  - `src/pcp/blocks/`：分块与 B 样条拟合（基函数、全局最小二乘、Jacobian 评估）；
  - `src/pcp/io/`：G+Smo XML 写出与读回、VTU 六面体单元；
  - 对应的单元测试和模型测试。
- **修改代码**：
  - `src/pcp/pipeline/`：新增 `Stage::Blocks`、`Stage::Fit`、参数段 `fit`、输出文件和统计字段；
  - `apps/pcp_cli/`：新增 `--until blocks|fit` 和拟合参数选项；
  - README 与 `docs/experiments.md`。
- **依赖**：不新增第三方依赖。稀疏最小二乘用 Eigen 自带的求解器（不用 SuiteSparse）；G+Smo 只是输出格式，不链接它的库。
- **分支**：本 change 依赖尚未合并到 main 的 `add-chen-height-simplification`，在其基础上的分支 `add-block-bspline-fitting` 开发。
- **不在本 change 范围内**：
  - 控制点质量优化（消除负 Jacobian）；
  - 合并相邻块以减少块数；
  - 改变拓扑的台阶去除；
  - 高阶连续（C1）拼接；
  - 曲面拟合加 Coons 体插值的做法（Chen §6 的另一条路线）。
