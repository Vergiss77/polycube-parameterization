# Tasks

## 1. 分块（spec: polycube-block-decomposition；design D1、D2）

- [x] 1.1 新建 `src/pcp/blocks/block_layout.{h,cpp}` 并接入 CMake：
  - 由面片轴向与高度组得到每轴升序切割坐标；
  - 构造全局网格，按格子中心的广义绕数判定块；
  - 按 (z, y, x) 字典序编号；
  - 给每个块的 6 个面分类：内部面记录相邻块及对应面；边界面用 AABB 找所属面片，并检查轴向一致，不一致时报错。

  验证：单元测试覆盖
  - L 形三立方体得到 3 个块、体积和为 3；
  - U 形 polycube（复用 `make_voxels`）空隙格子不是块，凸台与底板之间为内部面；
  - 每个面都被分类，相邻关系对称；
  - 块体积和与四面体有向体积相差 ≤ 1e-9；
  - 重复运行结果逐位一致。
- [x] 1.2 流水线新增 `Stage::Blocks`（`--until blocks`）：
  - 简化启用时用简化后的 polycube 与高度组，关闭时用简化前的；
  - 写出 `06_blocks.json`（切割坐标为输出单位）、`06_blocks_polycube.vtu`、`06_blocks_input.vtu`（逐四面体 `block`，按重心）；
  - `stats.json` 中新增 `blocks` 字段；
  - README 补充新阶段与输出文件。

  验证：I-SPAR 上分别以 `--until blocks` 和关闭简化运行，输出文件齐全，每轴切割坐标数等于简化报告的切割面数（关闭时为简化前的高度组数），统计字段无 NaN。

## 2. B 样条基础（design D4、D9）

- [x] 2.1 新建 `src/pcp/blocks/bspline.{h,cpp}`：
  - 开放均匀节点向量；
  - 基函数及一阶导数（Cox–de Boor），含 u = 1 端点的处理；
  - 三变量张量积求值及 ∂/∂u、∂/∂v、∂/∂w。

  验证：单元测试覆盖
  - 次数 1–4 的单位分解（误差 ≤ 1e-14）；
  - 导数与有限差分一致（≤ 1e-7）；
  - 控制点取 Greville 点时线性函数被精确重现；
  - 端点插值角控制点。

## 3. 拟合（spec: bspline-volume-fitting；design D3–D6）

- [x] 3.1 在 `src/pcp/blocks/spline_fit.{h,cpp}` 中实现共享控制点的样条空间（D4）：
  - 每个区间的单元数 max(1, round(L/h))；
  - 全局一维偏移 O_{k+1} = O_k + n_k + p − 1；
  - 块局部控制点映射到全局格点，只为被引用的格点编号未知量。

  验证：单元测试中 L 形分块的每个内部面两侧引用的全局控制点集合完全相同，未知量数与手工计数一致，单元数公式在 L/h 取 0.3、1.4、2.6 时分别给出 1、1、3。
- [x] 3.2 实现采样与目标映射（D3）：
  - 每块每方向 n·s + 1 个参数值；
  - 经 `TetLocator(polycube X, T)` 映射到 `tets.X`；
  - 定位失败时向块中心收缩 1e-10 × 对角线后重试，再退回最近四面体，并计数；
  - 标记表面采样点。

  验证：单元测试中 L 形 polycube（目标网格为其仿射像）的全部采样点映射值与仿射像之差 ≤ 1e-12，定位失败数为 0，表面采样点数与手工计数一致。
- [x] 3.3 实现全局最小二乘（D5）：
  - 组装 A 与二阶差分光顺矩阵 S；
  - κ = tr(AᵀA)/tr(S)；
  - 用 `SimplicialLDLT` 一次分解求三个分量，分解失败时报错并提示 λ；
  - 串行、按固定顺序组装。

  验证：单元测试覆盖
  - 仿射目标、λ = 0 时最大误差 ≤ 1e-9 × 对角线；
  - 带一个空格子的分块在 λ > 0 时可解；
  - 重复运行控制点逐位一致。
- [x] 3.4 实现误差、Jacobian 与 C0 自检（D6）：
  - 全部与表面采样点的 RMS 和最大误差；
  - 边界评估点到原表面的距离（`point_mesh_squared_distance`）；
  - 每块每方向 n·4 + 1 个评估点上的缩放 Jacobian 统计（最小、平均、非正点数与比例、含非正点的块数）；
  - 每个内部面两侧求值的最大差。

  验证：单元测试覆盖
  - 仿射目标时缩放 Jacobian 恒等于其解析值、到表面距离 ≤ 1e-9；
  - 人为翻折一个控制点后报告的非正点数 > 0；
  - C0 偏差 ≤ 1e-12 × 对角线。

## 4. 输出格式（spec: G+Smo multipatch output / Sampled volume output；design D7、D8）

- [x] 4.1 新建 `src/pcp/io/gismo_xml.{h,cpp}`：
  - 写出 Geometry（TensorBSpline3）、节点向量、控制点（u 最快，17 位有效数字）；
  - 写出 MultiPatch 的 patches、interfaces、boundary，面编号按 boxSide 1–6；
  - 对照 G+Smo 源码中读取接口的代码确定接口行的写法，并在文件注释和 README 中注明出处；
  - 实现只解析本项目输出的读取器。

  验证：单元测试中 L 形拟合结果写出后读回，块数、次数、节点向量、控制点逐位一致，接口数和边界面数分别等于内部面数和边界面数。
- [x] 4.2 扩展 `io/vtu`：支持六面体单元（VTK_HEXAHEDRON）和逐点标量字段；生成每块评估网格的六面体采样网格，带逐单元 `block` 和逐点 `scaled_jacobian`。

  验证：单元测试中写出的 L 形采样网格单元数等于各块评估网格单元数之和，文件可被现有的 VTU 读回检查（或解析 XML 计数）；在 ParaView 中打开 I-SPAR 的结果，确认块和 Jacobian 着色正常。

## 5. 流水线与命令行接入（spec: Block and fitting stages in the pipeline；design D9）

- [x] 5.1 新增 `Stage::Fit`（默认运行到此阶段）和 `PipelineParams.fit`（`FitParams{degree, element_size, samples_per_element, smoothing}`，含 to_json / update_from_json，未知键和非法值报错）：
  - 写出 `07_bspline.xml`、`07_bspline_samples.vtu`；
  - `stats.json` 新增 `fit` 字段：次数、控制点数、采样点数、各项误差、Jacobian 统计、C0 偏差、定位失败数、耗时。

  验证：
  - I-SPAR 上以默认参数和 `--until blocks` 运行，检查输出文件有无、`params.json` 中的 `fit` 段，以及 `00_`–`05_` 输出和顶层统计字段与改动前一致；
  - 参数 JSON 往返与非法值报错有单元测试。
- [x] 5.2 `pcp_cli` 新增 `--degree`、`--element-size`，`--until` 支持 `blocks|fit`；更新 README（流水线图、阶段列表、选项表、输出文件、ParaView 查看步骤、用 G+Smo 读取 XML 的方法）。

  验证：
  - `pcp_cli --help` 列出新选项；
  - `--degree 0` 以非零状态退出并说明原因；
  - 按 README 的命令在 I-SPAR 上运行，结果与描述一致。
- [x] 5.3 在 `test_profile.cpp` 中增加默认跳过的 `profile_fit`：从已有输出目录读取 `01_tet.tmsh`、`05_polycube_simplified.tmsh`（或 `04_polycube.tmsh`）和分组，只运行分块与拟合，参数从环境变量读取。

  验证：对 CLEVIS 的已有输出运行该用例，结果与完整流水线的分块、拟合统计一致。

## 6. 模型验收与实验记录

- [x] 6.1 在 `model_ispar` 和 `model_clevis` 中增加分块与拟合检查：
  - 块体积和、面分类与相邻对称；
  - C0 偏差 ≤ 1e-12；
  - XML 读回一致；
  - 字段无 NaN；
  - I-SPAR 最大误差 ≤ 1e-3 且缩放 Jacobian 全为正；
  - CLEVIS 表面 RMS ≤ 5e-3；
  - I-SPAR 重复运行的 `07_bspline.xml` 完全相同。

  验证：`ctest -R '^model_(ispar|clevis)$'` 通过。
- [x] 6.2 在 `model_tailsitter` 中增加相同的检查（表面 RMS ≤ 5e-3），实测并记录块数、最薄块尺寸、Jacobian 非正比例。若误差超标，先分析原因，放宽前与用户确认，并记为已知限制。

  验证：`ctest -R '^model_tailsitter$'` 通过。
- [x] 6.3 用 `profile_fit` 在三个模型上对比单元尺寸（0.1、0.05、0.025）和光顺权重（1e-4、1e-6、1e-8）对误差、Jacobian 非正比例、控制点数和耗时的影响，确定默认值；在 `docs/experiments.md` 新增"分块与 B 样条拟合"一节（含全局网格分块与 Chen 局部切割的差异、块数统计、Jacobian 分布与已知限制、ParaView 截图）。

  验证：文档中三个模型、各组参数的数据完整，默认值与代码一致。
- [ ] 6.4 在本地按 CI 命令序列复现（全新构建目录：构建、`ctest -L unit`、各模型测试），更新 README 中的模型测试时长；推送后，CI 的构建加单元测试在两个平台通过，并手动触发一次带模型测试的运行。

  验证：本地全部通过，GitHub Actions 两个平台的任务通过。

## Workflow follow-up

- 本分支基于 `add-chen-height-simplification`。该 change 的模型测试通过并合并到 main 后，先归档它（生成 `polycube-simplification` 主 spec），再把本分支变基到 main。
- 本 change 完成后归档，生成 `polycube-block-decomposition` 和 `bspline-volume-fitting` 主 spec。
- 下一个 change：控制点质量优化（消除负 Jacobian），视块数实测结果再考虑合并块。
