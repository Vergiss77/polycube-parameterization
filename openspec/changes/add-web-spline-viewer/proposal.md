# Proposal

## Why

流水线现在输出 B 样条多片体（`07_bspline.xml`），但没有合适的方式查看它：
- `07_bspline_samples.vtu` 要用 ParaView 打开，开发机上没有装；
- 文档里的截图是 matplotlib 静态渲染，看不清细节，也没法交互。

下一步要做控制点质量优化（消除负 Jacobian），需要反复对比优化前后的样条体、查看负 Jacobian 的位置和控制网。一个直接读样条文件、可交互的查看器，是那项工作的前提。

## What Changes

- **新增 Web 查看器** `apps/web_viewer/`：
  - 纯 HTML + ES 模块 + three.js，不需要构建步骤；
  - three.js 用固定版本，放进仓库，离线可用；
  - 用本地静态服务（如 `python3 -m http.server`）打开。
- **前端直接读样条文件**：
  - 在浏览器里解析流水线输出的 G+Smo MultiPatch XML；
  - 用 JavaScript 求值三变量 B 样条，和 C++ 是同一套公式；
  - 不修改 C++ 流水线和它的输出。
- **显示与交互**：
  - 边界面的离散网格：每个单元的细分数可调，可显示网格线；
  - 着色：按块、按缩放 Jacobian（带色标）、纯色；
  - 标出缩放 Jacobian 非正的点；
  - 信息面板：块数、接口数、控制点数，以及与流水线统计一致的 Jacobian 统计；
  - 可显示块间内部接口面，以及沿坐标轴的剖切平面；
  - 可显示控制点与控制网；
  - 点选一个块，显示它的编号和质量；
  - 同时加载两个结果并排对比，两边相机同步。
- **测试**：
  - 用 Node 测试 JS 的解析与求值；
  - 对同一个小算例，JS 求值与 C++ 输出的 `07_bspline_samples.vtu` 逐点一致，Jacobian 统计与 C++ 一致；
  - 测试接入 ctest 的 `unit` 标签，CI 中自动运行。
- **文档**：README 新增查看器的使用说明。

## Capabilities

### New Capabilities
- `web-spline-viewer`：在浏览器中加载并显示 B 样条多片体，包括文件加载、求值与离散、着色与 Jacobian 标记、内部接口面与剖切、控制网显示、两结果对比，以及离线运行和与 C++ 求值一致的要求。

### Modified Capabilities
<!-- 无：查看器只读取已有的 07_bspline.xml，不改变流水线的行为与输出。 -->

## Impact

- **新增代码**：
  - `apps/web_viewer/`：页面与 JS 模块；
  - `apps/web_viewer/vendor/three/`：three.js 模块文件及其 MIT 许可证；
  - `tests/web/`：Node 测试。
- **修改代码**：
  - `tests/CMakeLists.txt`：找到 Node 时注册 JS 测试，并注册生成测试数据的前置测试；
  - 一个 C++ 单元测试，用来生成小算例的 XML 与采样网格；
  - README。
- **依赖**：
  - three.js（MIT）放进仓库；
  - 运行测试需要 Node，找不到 Node 时跳过 JS 测试并给出提示。GitHub Actions 的 Windows 和 macOS 环境都自带 Node；
  - 不引入 npm 依赖和打包工具。
- **分支**：基于尚未合并的 `add-block-bspline-fitting`，在分支 `add-web-spline-viewer` 上开发。
- **不在本 change 范围内**：
  - 在查看器里编辑控制点；
  - 读取 G+Smo 的其他几何类型（如 NURBS、二维片）；
  - 显示四面体网格等其他流水线输出；
  - 部署成公开网站。
