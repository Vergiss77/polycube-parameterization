# Design

## Context

动机见 proposal.md，行为要求见 `specs/web-spline-viewer/spec.md`。

查看器依赖的现有内容：

- **`07_bspline.xml`**：由 `src/pcp/io/gismo_xml.cpp` 写出。
  - 每块一个 `TensorBSpline3`：开放均匀节点向量，三个方向次数相同，控制点按 u 最快、w 最慢排列，17 位有效数字；
  - `<interfaces>` 每行 10 个整数，`<boundary>` 每行"块 面"，面编号 1–6 对应 u=0、u=1、v=0、v=1、w=0、w=1。
- **C++ 求值**：`src/pcp/blocks/bspline.cpp`。用 Piegl & Tiller A2.2/A2.3 算基函数及一阶导数，`u = 1` 归入最后一个非空节点区间。
- **C++ 质量统计**：`evaluate_fit` 在每块每方向 `n·4 + 1` 个评估点上计算缩放 Jacobian，`sample_mesh` 把这些点写成 `07_bspline_samples.vtu`（点按块依次排列，块内 u 最快）。
- **工具链**：开发机上有 Node v26 和 python3；GitHub Actions 的 Windows 和 macOS 环境都自带 Node。项目目前没有任何前端代码，也没有 npm 依赖。

## Goals / Non-Goals

**Goals:**
- 计算部分（解析、求值、离散、统计）不依赖 three.js 和 DOM，在浏览器和 Node 中运行同一份代码，能用 Node 直接测试。
- 不需要任何构建步骤：克隆仓库后起一个静态服务就能用。
- 与 C++ 的求值和统计逐点一致，查看器显示的就是真实的样条。

**Non-Goals:**
- 不做通用的 G+Smo 读取器（只支持本项目写出的 `TensorBSpline3` 多片体）。
- 不追求移动端或旧浏览器支持，只保证新版 Chrome、Safari、Firefox。
- 不在浏览器端修改或保存样条。

## Decisions

### D1 模块划分：计算层与显示层分离

```
apps/web_viewer/
  index.html            页面骨架、importmap、控件
  style.css
  src/
    gismo.js            解析 XML 字符串 -> { patches, interfaces, boundary }，严格校验
    bspline.js          findSpan、基函数及导数、体求值与 Jacobian、缩放 Jacobian
    tessellate.js       边界面/接口面 -> 位置、法向、块编号、缩放 Jacobian、索引、网格线
    quality.js          评估网格上的 Jacobian 统计与非正点列表（与 C++ 同网格）
    controlnet.js       控制点与控制网线
    view.js             three.js 场景、相机、拾取、剖切（一个视图一个实例）
    app.js              加载、界面状态、两个视图的同步
  vendor/three/         three.module.js、OrbitControls.js、LICENSE、VERSION
tests/web/
  *.test.mjs            Node 内置测试框架（node:test）
```

- 前五个模块只用纯 JavaScript 和类型化数组，可以在 Node 中直接导入。
- `view.js` 和 `app.js` 依赖 three.js 和 DOM，不在 Node 中测试，靠浏览器手动验证。
- **备选**：用 TypeScript 或打包工具（Vite 等）。类型更安全，但会引入 npm 依赖和构建步骤，和"克隆即用"冲突，所以不用。

### D2 XML 解析：基于字符串，不用 DOMParser

Node 没有 `DOMParser`，而计算层需要在 Node 中测试，所以用字符串查找来提取 `<Geometry>`、`<KnotVector degree=...>`、`<coefs>`、`<MultiPatch>`、`<patches>`、`<interfaces>`、`<boundary>`。格式与 C++ 的 `read_gismo_xml` 对应。

解析后严格校验以下各项，任何一项不满足就抛出带块编号和原因的错误：
- Geometry 的 id 为 0..N−1；
- 每个方向的节点向量单调不减，两端重复 p+1 次；
- 控制点数等于三个方向基函数数之积；
- 接口与边界引用的块存在，面编号在 1–6；
- 每个块的 6 个面恰好出现在接口或边界中一次。

界面层捕获错误并显示，不渲染任何几何（spec: Invalid file）。

### D3 求值：逐行对应 C++ 实现

- `bspline.js` 按 `bspline.cpp` 的算法和循环顺序移植：`findSpan`、`basisFunctions`（带一阶导数）、`evalPoint`、`evalJacobian`。全部用 JS 的双精度数，结果与 C++ 只差舍入量级。spec 要求 ≤ 1e-12 × 对角线。
- 每块预先算好三个方向的节点向量和"单元边界参数值"（去重后的节点），离散和统计都用它们。
- **性能**：给定一组参数值，基函数按方向预先计算后在面上复用（张量积结构），不在每个点上重复计算三次。

### D4 边界面与接口面的离散

- **采样点**：每个面（块 b、面 f）固定一个参数（f 为奇数时取 1，偶数时取 0），另外两个方向各取 `n_e · s + 1` 个均匀参数值（n_e 为该方向单元数，s 为细分数），全部包含单元边界。
- **无裂缝**：两个相邻块在公共棱上用相同的参数值求值，而公共棱上的控制点是共享的，因此位置一致到舍入量级。同一块相邻两个面在棱上用的是同一个体和同一组参数，结果逐位相同。
- **顶点属性**：位置、法向（取两个切向 Jacobian 列的叉积，按面的朝外方向定号）、块编号、该点的缩放 Jacobian（由同一次 `evalJacobian` 得到，不额外计算）。
- **几何体组织**：边界面合并成一个三角网格，内部接口面合并成另一个（只取每个接口的一侧，避免重复）。切换着色方式时只重写颜色数组，不重新求值；改细分数时才重新离散。
- **网格线**：单元边界线和全部细分线分别生成线段集合，按开关显示。

### D5 Jacobian 统计与非正点标记

- `quality.js` 在每块每方向 `n·4 + 1` 个评估点上计算缩放 Jacobian，网格和 C++ 的 `evaluate_fit` 完全相同，得到评估点数、非正点数、含非正点的块数、最小值、平均值，以及每块的最小值。
- 非正点作为一组点精灵显示，渲染时不做深度测试，所以位于体内部的点也能看到。
- 统计在加载时算一次，与细分数无关。
- TAILSITTER 共 67 750 个评估点，计算量远小于离散。

### D6 着色

- **按块**：用块之间的邻接关系（接口）做贪心图着色，在 12 色调色板里给相邻块分配不同颜色。
- **按缩放 Jacobian**：用连续色标（viridis），范围默认为 [min, 1]，也可以手动设定；非正值用醒目的单独颜色。旁边显示色标图例。
- **纯色**：配合网格线查看离散。

### D7 剖切、拾取与控制网

- **剖切**：用 three.js 材质的 `clippingPlanes`（启用 `renderer.localClippingEnabled`）。一个可选轴向的平面，位置由滑块在模型包围盒范围内调节，可以翻转保留的一侧，对边界面、接口面、控制网和非正点标记同时生效。
- **拾取**：对边界面网格做射线检测，根据命中三角形的块编号属性高亮该块（单独叠加一层半透明网格），并在信息面板显示该块的信息（spec: Block picking）。
- **控制网**：每块的控制点及沿 u、v、w 三个方向相邻控制点的连线，可以显示全部块或只显示选中的块。TAILSITTER 全部块约 3 万个控制点，用一个线段几何体即可。

### D8 两个结果的并排对比

- 页面分成左右两个视图，各自一个 `<canvas>`、一个 `View` 实例（独立的场景和相机）和一个信息面板。只加载一个文件时只显示一个视图。
- **相机同步**：一个视图的 OrbitControls 触发 `change` 时，把相机位置、朝向和目标点复制到另一个视图（用标志位防止相互递归触发）。可以取消同步。
- 显示设置保存在一个共享的状态对象里，改动时同时作用于两个视图。
- **备选**：一个 renderer 加 scissor 分屏。可以省一个 WebGL 上下文，但拾取和尺寸处理更复杂。两个视图的开销可以接受，选择简单的做法。

### D9 three.js 的引入与页面运行方式

- **版本**：固定一个版本（实现时取 0.170.0，若该版本的模块拆分方式不便，就换邻近版本，并在 `VERSION` 中记录），用 `npm pack three@<版本>` 取得官方发布包。复制以下文件到 `vendor/three/`：
  - `build/three.module.js`（以及它导入的同目录文件，如有）；
  - `examples/jsm/controls/OrbitControls.js`；
  - `LICENSE`。
  在 `VERSION` 中写明版本号和文件的 SHA-256。
- **importmap**：`index.html` 中把 `"three"` 映射到 `./vendor/three/three.module.js`，把 `"three/addons/"` 映射到 `./vendor/three/addons/`，这样 OrbitControls 的 `import ... from 'three'` 无需修改即可工作。
- **运行方式**：ES 模块不能从 `file://` 加载，必须通过 HTTP 访问。推荐在仓库根目录运行 `python3 -m http.server 8765`，打开 `http://localhost:8765/apps/web_viewer/`。
  - 地址参数 `?a=<相对路径>&b=<相对路径>` 直接加载一个或两个文件，例如 `?a=../../out/clevis/07_bspline.xml`；
  - 也可以用文件选择或拖放，这两种方式不受服务根目录限制。
- 页面不引用任何外部地址（spec: Offline operation）。

### D10 测试与 CI 接入

- **测试数据**：新增一个 C++ 测试用例（doctest，属于 `unit`）。它对 U 形 polycube 做一次小规模拟合，在环境变量 `PCP_WEB_FIXTURE_DIR` 给出的目录里写出 `07_bspline.xml`、`07_bspline_samples.vtu` 和一份 JSON，其中是 C++ 的 Jacobian 统计。环境变量未设置时，这个用例只做内存中的检查，不写文件。
  - U 形 polycube 有多个块、内部接口和边界面，而且样条不是仿射的，求值对比有意义。
- **ctest 注册**：`tests/CMakeLists.txt` 中：
  - `web_fixture`：运行上述用例，设置环境变量，属性 `FIXTURES_SETUP web`；
  - `web_viewer_js`：`node --test tests/web/`，属性 `FIXTURES_REQUIRED web`，标签 `unit`；
  - 用 `find_program(NODE_EXECUTABLE node)` 查找 Node，找不到时不注册 JS 测试，只打印提示。
  - 现有 CI 的 `ctest -L unit` 会自动带上这些测试。
- **Node 测试内容**：
  - XML 解析正常与各类错误；
  - 求值与 `07_bspline_samples.vtu` 逐点对比（位置 ≤ 1e-12 × 对角线，缩放 Jacobian ≤ 1e-9）；
  - Jacobian 统计与 C++ 的 JSON 一致；
  - 离散结果无裂缝（相邻块公共棱上的点位置一致）；
  - 细分数改变时三角形数符合公式；
  - 控制网的点数与线段数符合公式。
- **界面验证**：界面部分在浏览器中手动验证。验证时用三个模型的真实输出，检查各功能并截图写进文档。若本机可以用浏览器自动化工具，就用它截图；否则由用户确认。

## Risks / Trade-offs

- **[three.js 版本的模块结构变化]** 新版本把 `three.module.js` 拆成多个文件，importmap 的映射需要对应调整。
  → 实现时以实际发布包为准，复制它导入的全部文件，并在 `VERSION` 中记录。页面加载失败会直接在控制台报错，容易发现。
- **[仓库体积]** three.js 的模块文件约 1 MB 多。
  → 只复制用到的文件（主模块、OrbitControls、许可证），不复制整个包。
- **[JS 与 C++ 不是逐位相同]** 浮点求和顺序的细微差别可能导致最后几位不同。
  → spec 只要求 1e-12 × 对角线的一致性；求值代码按 C++ 的循环顺序编写，实际差异预计远小于这个量级。统计中的非正点判定（≤ 0）在 TAILSITTER 上不处于临界值附近（最小值为 −0.49）。
- **[大文件解析性能]** TAILSITTER 的 XML 约 3 万行控制点。
  → 字符串解析加 `parseFloat` 的开销在毫秒到百毫秒量级。信息面板显示加载和离散的耗时，便于发现问题。
- **[Windows CI 上的 Node 路径]** Windows 运行器上 `node` 可能叫 `node.exe`。
  → `find_program` 会自动处理扩展名。若 CI 环境异常，测试会被跳过而不是失败，CI 日志里能看到提示。
- **[界面无法自动测试]** 拾取、剖切、同步等交互没有自动测试。
  → 交互逻辑尽量放进可测的纯函数（如剖切范围计算、图着色），三维交互部分靠手动验证清单覆盖。

## Migration Plan

- 只新增文件，不改变流水线的行为和输出，没有迁移问题。
- 回退：删除 `apps/web_viewer/`、`tests/web/` 和 `tests/CMakeLists.txt` 中新增的测试即可。
