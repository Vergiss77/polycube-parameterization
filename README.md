# polycube-parameterization

基于 polycube 的多块体参数化实验代码。当前流水线：

```
OFF 表面 -> 校验 -> TetGen 四面体化 -> l1 polycube 变形 (Huang et al. 2014)
         -> 标号清理 / 面片图 -> 带约束 ARAP 的最终 polycube 映射
         -> 有向图高度简化 (Chen et al. 2019) -> 按新高度重新映射
         -> 沿切割面分块 -> 三变量 B 样条多片体拟合（块间 C0）
```

后续计划：控制点质量优化（消除负 Jacobian）。

## 构建

依赖（libigl 及其锁定的 Eigen 3.4、TetGen、nlohmann/json、CLI11、doctest）在 CMake 配置时
自动下载，只需要编译器、CMake（>= 3.20）和 git。

macOS（Apple clang）：

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

Windows（Visual Studio 2022，在 "x64 Native Tools" 命令行中）：

```bat
cmake -S . -B build
cmake --build build --config Release --parallel
```

可选的 polyscope 查看器默认不编译（`-DPCP_BUILD_VIEWER=ON` 打开）。

## 测试

```sh
ctest --test-dir build -C Release -L unit --output-on-failure   # 单元测试（含 Web 查看器的 Node 测试），约 15 秒
ctest --test-dir build -C Release -L models --output-on-failure # 三个模型的集成测试，较慢
```

### CI（GitHub Actions）

- 每次 push 或 PR：在 Windows 和 macOS 上构建并运行单元测试，几分钟内结束。
- 模型集成测试：只在手动触发时运行（本地约 80 分钟，云端更慢）。在仓库的 **Actions** 页面选择 **CI** → **Run workflow**，勾选 *Run the model integration tests* 后运行。

## 命令行用法

```sh
build/apps/pcp_cli/pcp_cli models/CLEVIS-D64EEBFA.off -o out/clevis
```

常用选项（完整列表见 `pcp_cli --help`）：

| 选项 | 说明 |
|---|---|
| `--until validate\|tet\|l1\|structure\|simplify\|blocks\|fit` | 只运行到指定阶段（默认 `fit`，即全部阶段） |
| `--resume-from out/clevis/01_tet.tmsh` | 用保存的四面体网格续跑，跳过四面体化 |
| `--params p.json` | 参数文件，格式同输出目录中的 `params.json`；命令行选项会覆盖它 |
| `--beta 2` | 复杂度权重 β（越大 polycube 越简单） |
| `--no-rotation` | 关闭全局朝向搜索，结果对齐世界坐标轴 |
| `--feature-edges e.txt --omega 100` | 特征边引导，文件每行 `顶点a 顶点b corner\|flat` |
| `--linear-solver auto\|direct\|cg` | Newton 线性系统求解方式（默认按规模自动选择） |
| `--tet-size-factor 2` | 四面体边长 / 表面平均边长 |
| `--no-simplify` | 跳过高度简化（等价于旧版流水线的输出） |
| `--max-stretch 2` | 高度简化中每条棱线/间隔的最大拉伸比 s（≥ 1，越大合并越多、畸变越大） |
| `--degree 3` | B 样条体的次数（≥ 1） |
| `--element-size 0.05` | B 样条目标单元尺寸，相对于 polycube 包围盒对角线；每个切割区间分 max(1, round(长度 / 尺寸)) 个单元 |
| `-q` / `-v` | 安静 / 详细输出 |

退出码：成功为 0，任何错误为 1（错误信息输出到 stderr，并写入 `stats.json` 的 `error` 字段）。

## 输出文件

| 文件 | 内容 |
|---|---|
| `params.json` | 本次运行实际使用的全部参数 |
| `stats.json` | 统计：输入规模与亏格、四面体数、全局旋转、polycube 误差、ARAP 畸变、面片/角点/棱线数、翻转四面体数、阶段与迭代数、各阶段耗时；`simplify` 字段为高度简化的结果（简化前后切割面数、合并次数、长度比范围、简化后畸变与翻转数）；`blocks` 字段为分块结果（块数、每轴切割坐标数、网格格子数、内部面与边界面数、体积相对差、最薄块的边长与对角线之比）；`fit` 字段为拟合结果（控制点数、采样点数、拟合误差与到原表面的距离（均除以输入包围盒对角线）、缩放 Jacobian 的最小/平均值与非正点数、C0 最大偏差、找不到所在四面体的采样点数）。顶层的畸变和翻转数指 `04` 的结果 |
| `log.txt` | 完整日志（含每个调度阶段的 α、ε、误差、畸变） |
| `00_validate.json` | 输入表面校验结果 |
| `01_tet.tmsh` / `01_tet.vtu` | 四面体网格（`.tmsh` 为无损文本格式，用于续跑） |
| `02_l1.tmsh` / `02_l1.vtu` | ℓ1 变形结果，`.vtu` 带逐四面体 `arap` 字段 |
| `03_labels_input.vtu` | 原始表面，带逐面 `label`（0..5 = +X,-X,+Y,-Y,+Z,-Z）和 `patch` 字段 |
| `03_labels_polycube.vtu` | 最终 polycube 表面，字段同上 |
| `03_patch_graph.json` | 面片图：面片（轴向、符号、高度、面积、三角形）、棱线（平行轴、两侧面片、两端同轴面片）、角点 |
| `04_polycube.tmsh` / `04_polycube.vtu` | 最终 polycube 四面体网格（简化前），带 `arap` 与 `inverted` 字段 |
| `05_simplify.json` | 高度简化：每个面片的简化前后高度和所属高度组、按顺序的合并记录（两组、合并高度、代价）、每轴切割面数、简化用的有向图（归一化单位） |
| `05_polycube_simplified.tmsh` / `.vtu` | 按简化后高度重新求解的 polycube，带 `arap` 与 `inverted` 字段 |
| `06_blocks.json` | 分块结构：每轴切割坐标（输出单位）；每个块的网格位置、包围盒与体积（归一化的 polycube 坐标）、6 个面（`x_lo`…`z_hi`）的相邻块或所属面片 |
| `06_blocks_polycube.vtu` / `06_blocks_input.vtu` | polycube 与原模型上的四面体网格，带逐四面体 `block` 字段（按四面体在 polycube 上的重心所在的块） |
| `07_bspline.xml` | B 样条多片体，G+Smo MultiPatch XML：每块一个 `TensorBSpline3`（id 与块编号相同，控制点按 u 最快、w 最慢排列，输入模型的长度单位与坐标系），以及块间接口和边界面 |
| `07_bspline_samples.vtu` | 每块评估网格上的样条体采样（六面体），带逐单元 `block` 和逐点 `scaled_jacobian` 字段 |

所有几何量都使用输入模型的长度单位；polycube 结果包含全局旋转。

## 在 ParaView 中查看

1. `File > Open`，选择 `03_labels_polycube.vtu`，点 `Apply`，着色选 `label` 或 `patch`。
2. 同时打开 `03_labels_input.vtu` 对照原模型上的分区。
3. 打开 `04_polycube.vtu`，着色选 `arap` 查看畸变分布；用 `Filters > Threshold`
   （字段 `inverted`，范围 1~1）检查是否有翻转单元；用 `Clip` 查看内部。
4. 打开 `05_polycube_simplified.vtu` 与 `04_polycube.vtu` 对比，查看高度简化把哪些面对齐到了同一平面。
5. 打开 `06_blocks_input.vtu`，着色选 `block`，查看分块落在原模型上的位置；`06_blocks_polycube.vtu` 是同样的分块在 polycube 上的样子。
6. 打开 `07_bspline_samples.vtu`，着色选 `scaled_jacobian`（点数据）查看样条体质量；用 `Threshold`（`scaled_jacobian` ≤ 0）找出 Jacobian 非正的位置。

## 在浏览器中查看样条

`apps/web_viewer/` 是一个不需要构建、可以离线使用的 Web 查看器。它在浏览器里直接读取 `07_bspline.xml`，用 JavaScript 求值 B 样条，结果和 C++ 一致：位置误差 ≤ 1e-12 × 对角线，Jacobian 统计逐项相同。

启动：在仓库根目录起一个静态服务，再用浏览器打开页面。

```sh
python3 -m http.server 8765
# 浏览器打开 http://localhost:8765/apps/web_viewer/
```

端口被占用时换一个即可（例如 `python3 -m http.server 9876`），地址中的端口号随之修改。页面不能直接从 `file://` 打开，因为浏览器不允许从本地文件加载 ES 模块。

加载文件有三种方式：
- 用工具栏的"文件 A / 文件 B"选择文件；
- 把文件拖到视图里；
- 用地址参数（路径相对于页面），例如：
  - `http://localhost:8765/apps/web_viewer/?a=../../out/clevis/07_bspline.xml`
  - `...?a=../../out/clevis/07_bspline.xml&b=../../out/clevis_h0025/07_bspline.xml`：左右并排对比两个结果

功能：

| 控件 | 作用 |
|---|---|
| 细分 | 每个单元每个方向的细分数（1–32） |
| 单元线 / 细分线 | 显示样条的单元划分，或全部离散网格线 |
| 着色 | 按块（相邻块颜色不同）、缩放 Jacobian（viridis 色标，≤ 0 为红色）、纯色 |
| 非正点 | 标出评估网格上缩放 Jacobian ≤ 0 的点，包括体内部的点，与 `stats.json` 中 `fit.nonpositive_jacobian` 一致 |
| 接口面 | 半透明显示块与块之间的内部接口面 |
| 剖切 | 沿 x/y/z 剖切，滑块调位置，"翻转"切换保留的一侧 |
| 控制网 | 显示选中块或全部块的控制点与控制网 |
| 点击模型 | 选中一个块，信息面板显示块编号、单元数、控制点数和该块的 Jacobian |
| 相机同步 / 重置视角 | 两个视图同步旋转缩放；恢复默认视角 |

左上角的信息面板显示块数、接口数、控制点数（去重后与 `stats.json` 中 `fit.control_points` 一致）、Jacobian 统计以及各步耗时。

three.js 0.170.0（MIT）放在 `apps/web_viewer/vendor/three/`，版本和校验和见其中的 `VERSION`。计算部分的模块（`gismo.js`、`bspline.js`、`quality.js`、`tessellate.js`、`controlnet.js`、`colors.js`）不依赖浏览器，Node 测试在 `tests/web/`，通过 ctest 的 `unit` 标签运行（需要安装 Node，没有时自动跳过）。

## 用 G+Smo 读取样条体

`07_bspline.xml` 按 G+Smo（stable 分支）读取 MultiPatch 的代码写出：`src/gsIO/gsXml.cpp` 的 `getInterfaces` / `getBoundaries`，`src/gsIO/gsXmlUtils.hpp` 的 `gsXml<gsMultiPatch>::get`，面编号见 `src/gsCore/gsBoundary.h`（1–6 = u=0, u=1, v=0, v=1, w=0, w=1）。每个接口一行：`块1 面1 块2 面2`，接着是方向映射（恒等 `0 1 2`）和方向符号（法向为 0，切向为 1），例如 `0 2 1 1 0 1 2 0 1 1`。

本项目不依赖 G+Smo，CI 也不检查这一步。装有 G+Smo 时可以这样验证：

```cpp
#include <gismo.h>
int main() {
  gismo::gsMultiPatch<> mp;
  gismo::gsReadFile<>("out/clevis/07_bspline.xml", mp);
  gismo::gsInfo << mp.nPatches() << " patches, " << mp.interfaces().size() << " interfaces, "
                << mp.boundaries().size() << " boundary sides\n";
  // 与 stats.json 中 blocks.blocks / interior_faces / boundary_faces 一致即说明拓扑读入正确
}
```

## 参考文献

- J. Huang, T. Jiang, Z. Shi, Y. Tong, H. Bao, M. Desbrun. ℓ1-Based Construction of Polycube Maps
  from Complex Shapes. ACM TOG 33(3), 2014.
- L. Chen, G. Xu, S. Wang, Z. Shi, J. Huang. Constructing volumetric parameterization based on
  directed graph simplification of ℓ1 polycube structure from complex shapes. CMAME 351, 2019.

TetGen 采用 AGPLv3 许可证，本项目仅用于研究。
