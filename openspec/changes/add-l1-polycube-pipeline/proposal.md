# Proposal

## Why

本项目的目标是参照 Huang 2014（ℓ1 polycube）与 Chen 2019（有向图简化 + 多片 B 样条），从三角网格自动构造基于 polycube 的三变量 B 样条多片体参数化。整条流水线的第一块基石是"从任意闭合三角网格得到低畸变 polycube 映射"，而当前仓库还没有任何代码。本 change 搭建跨平台工程骨架并完整实现 Huang 2014 的 ℓ1 polycube 构造，使后续的结构简化（Chen §5）、分块与样条拟合可以直接在其输出上开展。

## What Changes

- 新建 C++17 / CMake 工程，依赖全部通过 FetchContent 锁定版本获取（libigl 及其锁定的 Eigen 版本、TetGen，另有若干仅头文件的小工具库），在 Windows (MSVC) 与 macOS (Apple clang) 上均可构建；不引入 SuiteSparse、Ipopt、OpenMP 等平台负担较重的依赖。
- 网格 I/O：读取 OFF 三角网格并校验（闭合、2-流形、无退化面），导出用于查看的 `.vtu`（四面体网格及逐单元/逐面字段）与 OFF/OBJ 表面。
- 四面体化：基于 TetGen 将闭合表面转为保边界的四面体网格。
- ℓ1 polycube 变形（Huang 2014 §2–§4）：面积约束下最小化 `α·E1 + Eδ(ARAP) + Eη(β 复杂度) [+ Ee 可选]`，含全局旋转求解、光滑 ℓ1 近似、α/ε 调度、防翻折 log barrier、Schur 补 KKT 求解。
- Polycube 结构提取（Huang §4.1）：按法向标号、锯齿清理、度 1/2 面片清理、面片邻接图与高度提取，以及最终"同面片共面"线性约束下的 ARAP 体映射。
- 命令行工具：对单个模型运行上述阶段，输出 polycube 网格、标号、面片图、统计指标（polycube 误差、ARAP 畸变、角点数、耗时）。
- 以 `models/` 中三个模型（I-SPAR、CLEVIS、TAILSITTER）作为验收用例。

## Capabilities

### New Capabilities
- `build-platform`: 跨平台构建与依赖管理约束（支持的平台/编译器、依赖获取方式、禁止的依赖）。
- `mesh-io`: 三角表面网格读取与校验，四面体/表面网格及字段的导出。
- `tet-meshing`: 由闭合三角表面生成保边界的四面体网格。
- `polycube-deformation`: 基于 ℓ1 法向范数的 polycube 变形优化（全局旋转、正则化、复杂度控制、调度与收敛判据）。
- `polycube-structure`: 从变形结果提取合法的 polycube 拓扑结构（标号、清理、面片图、面片高度）并计算最终体内 polycube 映射。
- `pipeline-cli`: 运行流水线阶段并输出中间结果与统计信息的命令行入口。

### Modified Capabilities
<!-- 仓库当前没有任何既有 spec。 -->

## Impact

- 新增代码：`CMakeLists.txt`、`cmake/`、`src/`（核心库）、`apps/`（CLI）、`tests/`。
- 新增第三方依赖（构建期拉取，不需系统安装）：Eigen、libigl、TetGen（AGPLv3，项目仅用于研究，已确认可接受）；polyscope 作为可选调试查看器。
- 新增 CI：GitHub Actions 在 `windows-latest` 与 `macos-latest` 上构建并运行单元测试；模型集成测试手动触发。
- 不在本 change 范围内：Chen 的有向图简化与按新高度重变形、分块、B 样条拟合与质量优化（后续 change）；细长三角形塌缩（Huang §4 中 barrier 卡住时的重网格化）本次只做检测与日志。
