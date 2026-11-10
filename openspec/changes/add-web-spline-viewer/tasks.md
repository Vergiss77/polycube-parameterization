# Tasks

## 1. 页面骨架与 three.js（design D1、D9）

- [x] 1.1 用 `npm pack three@0.170.0`（必要时换邻近版本）取得发布包，把 `three.module.js`（及其导入的同目录文件）、`OrbitControls.js`、`LICENSE` 复制到 `apps/web_viewer/vendor/three/`，写 `VERSION`（版本与 SHA-256）；新建 `index.html`（importmap、左右两个视图容器、控件栏）、`style.css` 和 `src/app.js` 的最小版本（显示一个坐标轴和空场景）。

  验证：在仓库根目录运行 `python3 -m http.server 8765`，打开 `http://localhost:8765/apps/web_viewer/`，页面显示空场景并可旋转；浏览器控制台无错误，网络面板中没有外部请求。

## 2. 计算层：解析、求值、离散（spec: Load a multipatch file / Evaluation consistent with the pipeline / Boundary surface mesh / Control net display；design D2–D5、D10）

- [x] 2.1 实现 `src/gismo.js`：基于字符串的解析，以及 D2 列出的全部校验，错误信息带块编号和原因。

  验证：`tests/web/gismo.test.mjs` 覆盖：
  - 解析 C++ 写出的文件，块数、节点向量、控制点、接口和边界与文件一致；
  - 缺少 MultiPatch、控制点数不足、接口引用不存在的块、面编号越界、某个面既不在接口也不在边界中，都抛出说明原因的错误。

  `node --test tests/web/` 通过。
- [x] 2.2 实现 `src/bspline.js`（按 `bspline.cpp` 的算法和循环顺序移植）和 `src/quality.js`（与 C++ 相同的评估网格和统计）。新增 C++ 用例：在 `PCP_WEB_FIXTURE_DIR` 中写出 U 形 polycube 的 `07_bspline.xml`、`07_bspline_samples.vtu` 和 Jacobian 统计 JSON。在 `tests/CMakeLists.txt` 注册 `web_fixture`（FIXTURES_SETUP）和 `web_viewer_js`（FIXTURES_REQUIRED，标签 `unit`）；找不到 Node 时跳过并提示。

  验证：`ctest -L unit` 运行 JS 测试，其中：
  - 求值与采样网格逐点一致：位置 ≤ 1e-12 × 对角线，缩放 Jacobian ≤ 1e-9；
  - 统计与 C++ 的 JSON 一致：评估点数、非正点数、含非正点块数相同，最小值与平均值差 ≤ 1e-9；
  - 单位分解与导数的有限差分检查通过。
- [x] 2.3 实现 `src/tessellate.js`（边界面与接口面的离散，生成顶点属性、索引、单元边界线与细分线）和 `src/controlnet.js`。

  验证：Node 测试覆盖：
  - 细分数 s 时每个面的顶点数为 (n_u·s+1)(n_v·s+1)，三角形数为 2·n_u·n_v·s²；
  - 相邻块公共棱上的顶点位置差 ≤ 1e-12 × 对角线（无裂缝）；
  - 法向朝外（与面的朝外方向点积为正）；
  - 接口面只取一侧；
  - 控制网点数为各块控制点数之和，线段数符合公式。

## 3. 显示与交互（spec: Load / Coloring and Jacobian markers / Interior interfaces and clipping / Control net display / Block picking / Side-by-side comparison；design D4–D8）

- [x] 3.1 实现 `src/view.js` 和 `src/app.js` 的单视图功能：
  - 文件选择、拖放、地址参数 `?a=` 加载，错误信息显示；
  - 边界面网格与细分数滑块（1–32，默认 8），单元边界线和细分线开关；
  - 三种着色方式与色标；
  - 非正点标记开关；
  - 信息面板：块数、接口数、边界面数、控制点数、次数、Jacobian 统计、加载与离散耗时。

  验证：浏览器中加载三个模型的 `07_bspline.xml`，信息面板数值与各自 `stats.json` 一致（TAILSITTER：67 750 / 29 / 15）；切换着色和细分数正常；加载损坏的文件时显示错误信息。
- [x] 3.2 实现内部接口面（半透明，可开关）和剖切平面（轴向选择、位置滑块、翻转），剖切对所有图层生效；按块着色使用接口邻接关系做图着色（纯函数，加 Node 测试：相邻块颜色不同）。

  验证：
  - CLEVIS 打开接口面后显示 62 个内部面；
  - 沿 z 轴剖切到中部可以看到内部的块结构；
  - 图着色的 Node 测试通过。
- [x] 3.3 实现块拾取（高亮，并在信息面板显示块编号、次数、单元数、控制点数、块内 Jacobian 最小值与非正点数），以及控制网显示（全部块 / 只显示选中块）。

  验证：
  - 点击 I-SPAR 的块，显示其信息，Jacobian 最小值约为 1；
  - 只显示选中块的控制网时，控制点数等于三个方向之积；
  - 点击 TAILSITTER 含非正点的块时，非正点数与统计一致。
- [x] 3.4 实现两个结果并排对比：`?a=&b=` 或分别选择文件，两个视图各自一个信息面板；相机同步（可取消）；显示设置同时作用于两边。

  验证：左边加载默认参数的 CLEVIS 结果、右边加载 `--element-size 0.025` 的结果，两边同步旋转，各自显示控制点数（5 160 / 13 488）。

## 4. 文档与整体验证

- [x] 4.1 README 新增"在浏览器中查看样条"一节：启动命令、地址参数、各功能说明、不能用 `file://` 打开的原因。用三个模型做一遍完整的手动验证清单：
  - TAILSITTER 以默认细分数在 10 秒内显示完成，信息面板给出耗时；
  - 断网后功能正常，网络面板无外部请求；
  - 对比视图正常。

  截取 CLEVIS 和 TAILSITTER（按 Jacobian 着色、显示非正点）的截图放进 `docs/images/`，在 `docs/experiments.md` 中引用。

  验证：按 README 的命令操作，结果与描述一致，清单各项通过。
- [ ] 4.2 推送分支后，CI 的构建加单元测试（包含 `web_viewer_js`）在 Windows 和 macOS 两个平台通过。

  验证：GitHub Actions 两个平台的任务通过，日志中 `web_viewer_js` 被运行而不是被跳过。

## Workflow follow-up

- 本分支基于 `add-block-bspline-fitting`。该 change 的模型测试通过并合并到 main 后，把本分支变基到 main 再合并。
- 本 change 完成后归档，生成 `web-spline-viewer` 主 spec。
- 之后开始负 Jacobian 的诊断与控制点质量优化（可以用本查看器对比优化前后的结果）。
