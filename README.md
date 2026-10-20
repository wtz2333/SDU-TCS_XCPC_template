# SDU-TCS XCPC template

Author:wtz2333,yifusuyi,VictorYuan

## 内容

- `0-一切的开始.md` 至 `7-杂.md`：按专题整理的打印版模板，PDF 构建入口为 [pandoc/gen.sh](pandoc/gen.sh)。
- [cpp/](cpp/)：日常使用的独立 C++ 模板与代码片段，保留不同实现版本。部分文件需要补充头文件、类型别名、常量或配套模板，不能直接作为独立程序编译。
- [DP 优化](dp优化.md)、[Trick](Trick.md)、[电竞指南](电竞指南.md)：补充笔记，不自动编入打印版 PDF。
- [dp.bat](dp.bat)：Windows 对拍脚本，需在当前目录准备 `gen.exe`、`std.exe`、`brute.exe`，其中生成器负责写入 `data.in`。
- `assets/`：模板插图；`pandoc/`：PDF 排版模板与构建工具。

PDF 构建环境及使用方式见 [pandoc/README.md](pandoc/README.md)。推送到 `main` 后，GitHub Actions 会构建 PDF 并发布到 `gh-pages` 分支。
