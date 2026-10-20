<div align="center">

# SDU-TCS · XCPC Template

**算法模板 · 专题笔记 · 赛场参考**

日常训练中积累的 C++ 板子，以及按专题整理、可构建为 PDF 的模板手册。

[![PDF Build](https://github.com/wtz2333/SDU-TCS_XCPC_template/actions/workflows/build.yml/badge.svg?branch=main)](https://github.com/wtz2333/SDU-TCS_XCPC_template/actions/workflows/build.yml)
[![C++](https://img.shields.io/badge/Language-C%2B%2B-00599C?logo=cplusplus&logoColor=white)](cpp/)
[![PDF](https://img.shields.io/badge/Read-PDF-B42318?logo=readthedocs&logoColor=white)](https://github.com/wtz2333/SDU-TCS_XCPC_template/blob/gh-pages/template.pdf)

[**阅读 PDF**](https://github.com/wtz2333/SDU-TCS_XCPC_template/blob/gh-pages/template.pdf) · [**下载 PDF**](https://raw.githubusercontent.com/wtz2333/SDU-TCS_XCPC_template/gh-pages/template.pdf) · [**C++ 模板**](cpp/) · [**电竞指南**](电竞指南.md)

维护：**wtz2333 · yifusuyi · VictorYuan**

</div>

---

## 从这里开始

| 你想做什么 | 入口 |
| :--- | :--- |
| 阅读或打印模板手册 | [PDF 预览](https://github.com/wtz2333/SDU-TCS_XCPC_template/blob/gh-pages/template.pdf) / [下载](https://raw.githubusercontent.com/wtz2333/SDU-TCS_XCPC_template/gh-pages/template.pdf) |
| 按专题查算法、公式和实现 | 下方专题导航 |
| 快速取用独立代码片段 | [cpp/](cpp/) · 52 份 C++ 文件 |
| 调整字体、封面和打印版排版 | [构建说明](pandoc/README.md) · [LaTeX 模板](pandoc/algo.latex) |
| 看训练、配合与比赛策略 | [电竞指南](电竞指南.md) |

## 专题导航

| 专题 | 内容速览 |
| :--- | :--- |
| [01 · 数据结构](1-数据结构.md) | ST 表、线段树、主席树、莫队、树状数组、并查集、Splay、LCT、李超树 |
| [02 · 图论](2-图论.md) | 树链剖分、LCA、虚树、网络流、图匹配、SCC、点双与边双、圆方树、2-SAT |
| [03 · 字符串](3-字符串.md) | Manacher、后缀数组、回文自动机、后缀自动机、AC 自动机、KMP、Z 函数、哈希 |
| [04 · 数学](4-数学.md) | 同余与逆元、CRT、BSGS、数论筛法、反演、组合计数、线性代数、博弈 |
| [05 · 多项式](5-多项式.md) | 拉格朗日插值、FFT、FWT、生成函数、集合幂级数、多项式运算 |
| [06 · 计算几何](6-计算几何.md) | 二维几何、动态凸包、最小圆覆盖、多种几何模板 |
| [07 · 杂项](7-杂.md) | 位运算、bitset、PBDS、对拍、Python、Bash、常数优化与注意事项 |

## C++ 模板索引

独立代码统一存放在 [`cpp/`](cpp/)，保留原文件名和不同实现版本，按需选取。

<details>
<summary><strong>展开全部 52 份 C++ 模板</strong></summary>

| 分类 | 文件 |
| :--- | :--- |
| 起手模板与工具 | [template](cpp/template.cpp) · [cf_template](cpp/cf_template.cpp) · [custom_hash](cpp/custom_hash.cpp) · [bitset](cpp/bitset.cpp) |
| 基础数据结构 | [DSU](cpp/DSU.cpp) · [Fenwick](cpp/Fenwick.cpp) · [ST](cpp/ST.cpp) |
| 线段树及扩展 | [SegmentTree](cpp/SegmentTree.cpp) · [SegTree](cpp/SegTree.cpp) · [lazysegment](cpp/lazysegment.cpp) · [Persistent_SegmentTree](cpp/Persistent_SegmentTree.cpp) · [LCTree（李超树）](cpp/LCTree.cpp) · [DynamicConvexShell](cpp/DynamicConvexShell.cpp) |
| 树上算法 | [HLD](cpp/HLD.cpp) · [LCA](cpp/LCA.cpp) · [Virtual_tree](cpp/Virtual_tree.cpp) · [TreeHash](cpp/TreeHash.cpp) · [DynamicDiameter](cpp/DynamicDiameter.cpp) |
| 连通性、流与匹配 | [BlockCutTree](cpp/BlockCutTree.cpp) · [EBCC](cpp/EBCC.cpp) · [euler_circuit](cpp/euler_circuit.cpp) · [Flow](cpp/Flow.cpp) · [MCF](cpp/MCF.cpp) · [FlowerTree](cpp/FlowerTree.cpp) |
| 字符串 | [ACAM](cpp/ACAM.cpp) · [KMP](cpp/KMP.cpp) · [Hash](cpp/Hash.cpp) · [manacher](cpp/manacher.cpp) · [PAM](cpp/PAM.cpp) · [SAM](cpp/SAM.cpp) · [Z_function](cpp/Z_function.cpp) |
| 数论与代数 | [Exgcd](cpp/Exgcd.cpp) · [QuickPow](cpp/QuickPow.cpp) · [dusum](cpp/dusum.cpp) · [Min25](cpp/Min25.cpp) · [NumberSqrt](cpp/NumberSqrt.cpp) · [Gauss](cpp/Gauss.cpp) · [LinearBasis](cpp/LinearBasis.cpp) · [Frac](cpp/Frac.cpp) · [Z](cpp/Z.cpp) |
| 多项式与插值 | [FFT](cpp/FFT.cpp) · [Poly](cpp/Poly.cpp) · [Poly_](cpp/Poly_.cpp) · [Poly_easyt](cpp/Poly_easyt.cpp) · [jiangly_poly](cpp/jiangly_poly.cpp) · [Language（拉格朗日插值）](cpp/Language.cpp) |
| 计算几何 | [ComputationalGeometry](cpp/ComputationalGeometry.cpp) · [CG_](cpp/CG_.cpp) · [CG_cjj](cpp/CG_cjj.cpp) · [CG_x](cpp/CG_x.cpp) · [dls_C](cpp/dls_C.cpp) · [DynamicConvexHull](cpp/DynamicConvexHull.cpp) |

</details>

部分文件是代码片段，需要补充头文件、类型别名、常量或配套模板。使用前请确认下标约定、数值范围与初始化方式；不同实现之间的接口不一定相同。

## PDF 构建

### 在线获取

推送到 `main` 后，[GitHub Actions](https://github.com/wtz2333/SDU-TCS_XCPC_template/actions/workflows/build.yml) 自动构建 PDF，并将结果发布到 `gh-pages` 分支。上方 PDF 入口指向最近一次成功发布的版本；最新提交是否已完成构建，可查看构建状态。

### 本地生成

在 Linux 或 WSL 中准备 Pandoc、XeLaTeX、latexmk，以及 Python 包 `pandocfilters`、`Pygments`。排版使用 Source Han Serif CN、Source Han Sans CN 和 Source Code Pro 字体，完整环境安装步骤可参考 [构建工作流](.github/workflows/build.yml)。

在仓库根目录运行：

```bash
mkdir -p ~/.pandoc/templates
cp pandoc/algo.latex ~/.pandoc/templates/
bash pandoc/gen.sh
```

生成结果为根目录的 `template.pdf`。字体配置在 [`pandoc/gen.sh`](pandoc/gen.sh)，封面与排版配置在 [`pandoc/algo.latex`](pandoc/algo.latex)。

> 打印版由根目录的 `0-一切的开始.md` 至 `7-杂.md` 构建。`cpp/` 中的独立代码和《电竞指南》不会自动编入 PDF。

## 仓库结构

```text
.
├── 0-一切的开始.md       # 手册起始页
├── 1-数据结构.md
├── 2-图论.md
├── 3-字符串.md
├── 4-数学.md
├── 5-多项式.md
├── 6-计算几何.md
├── 7-杂.md
├── cpp/                 # 独立 C++ 模板与代码片段
├── assets/              # 文档插图
├── pandoc/              # PDF 构建脚本与 LaTeX 排版模板
├── 电竞指南.md           # 训练、团队配合与比赛策略
└── .github/workflows/   # 自动构建与发布
```

## 来源与维护

本仓库在 [ECNU F0RE1GNERS 模板](https://github.com/xcpcio/template-Markdown-ECNU-F0RE1GNERS) 的基础上持续整理，结合日常训练补充代码与笔记，保留文件中已有的作者和来源说明。

发现错误或有改进建议，欢迎提交 [Issue](https://github.com/wtz2333/SDU-TCS_XCPC_template/issues) 或 Pull Request。算法修正请尽量附上适用条件、复现用例或题目链接，方便核对。
