# polycube-parameterization

基于 polycube 的多块体参数化实验代码。当前已实现流水线的前半段：

```
OFF 表面 -> 校验 -> TetGen 四面体化 -> l1 polycube 变形 (Huang et al. 2014)
         -> 标号清理 / 面片图 -> 带约束 ARAP 的最终 polycube 映射
         -> 有向图高度简化 (Chen et al. 2019) -> 按新高度重新映射
```

后续计划：分块、三变量 B 样条拟合与质量优化。

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
ctest --test-dir build -C Release -L unit --output-on-failure   # 单元测试，约 15 秒
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
| `--until validate\|tet\|l1\|structure\|simplify` | 只运行到指定阶段（默认 `simplify`，即全部阶段） |
| `--resume-from out/clevis/01_tet.tmsh` | 用保存的四面体网格续跑，跳过四面体化 |
| `--params p.json` | 参数文件，格式同输出目录中的 `params.json`；命令行选项会覆盖它 |
| `--beta 2` | 复杂度权重 β（越大 polycube 越简单） |
| `--no-rotation` | 关闭全局朝向搜索，结果对齐世界坐标轴 |
| `--feature-edges e.txt --omega 100` | 特征边引导，文件每行 `顶点a 顶点b corner\|flat` |
| `--linear-solver auto\|direct\|cg` | Newton 线性系统求解方式（默认按规模自动选择） |
| `--tet-size-factor 2` | 四面体边长 / 表面平均边长 |
| `--no-simplify` | 跳过高度简化（等价于旧版流水线的输出） |
| `--max-stretch 2` | 高度简化中每条棱线/间隔的最大拉伸比 s（≥ 1，越大合并越多、畸变越大） |
| `-q` / `-v` | 安静 / 详细输出 |

退出码：成功为 0，任何错误为 1（错误信息输出到 stderr，并写入 `stats.json` 的 `error` 字段）。

## 输出文件

| 文件 | 内容 |
|---|---|
| `params.json` | 本次运行实际使用的全部参数 |
| `stats.json` | 统计：输入规模与亏格、四面体数、全局旋转、polycube 误差、ARAP 畸变、面片/角点/棱线数、翻转四面体数、阶段与迭代数、各阶段耗时；`simplify` 字段为高度简化的结果（简化前后切割面数、合并次数、长度比范围、简化后畸变与翻转数）。顶层的畸变和翻转数指 `04` 的结果 |
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

所有几何量都使用输入模型的长度单位；polycube 结果包含全局旋转。

## 在 ParaView 中查看

1. `File > Open`，选择 `03_labels_polycube.vtu`，点 `Apply`，着色选 `label` 或 `patch`。
2. 同时打开 `03_labels_input.vtu` 对照原模型上的分区。
3. 打开 `04_polycube.vtu`，着色选 `arap` 查看畸变分布；用 `Filters > Threshold`
   （字段 `inverted`，范围 1~1）检查是否有翻转单元；用 `Clip` 查看内部。
4. 打开 `05_polycube_simplified.vtu` 与 `04_polycube.vtu` 对比，查看高度简化把哪些面对齐到了同一平面。

## 参考文献

- J. Huang, T. Jiang, Z. Shi, Y. Tong, H. Bao, M. Desbrun. ℓ1-Based Construction of Polycube Maps
  from Complex Shapes. ACM TOG 33(3), 2014.
- L. Chen, G. Xu, S. Wang, Z. Shi, J. Huang. Constructing volumetric parameterization based on
  directed graph simplification of ℓ1 polycube structure from complex shapes. CMAME 351, 2019.

TetGen 采用 AGPLv3 许可证，本项目仅用于研究。
