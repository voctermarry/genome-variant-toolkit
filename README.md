## 用途

本项目是「基因组变异分析工具链」的代码仓库，用于逐步实现该方向的序列处理、比对与变异分析能力。

当前实现了 FASTA/FASTQ 序列的读取、校验与规范化。

## 环境与安装

- Python 3.11 及以上

```bash
python -m pip install -e .
```

## 测试

```bash
python -m pytest
```

## 命令行入口

安装后提供 `genome-variant-toolkit` 命令：

```bash
genome-variant-toolkit version                 # 打印版本号
genome-variant-toolkit --help                  # 打印用法
genome-variant-toolkit normalize-sequences \
    INPUT [--input-format {fasta,fastq,auto}] \
    [--output-format {fasta,fastq}] [--output OUTPUT] [--line-width N]
```

`normalize-sequences` 将输入序列规范化后重新写出：

- `INPUT` 为文本文件路径，或 `-` 表示标准输入；`--output` 缺省或为 `-` 时写到标准输出。
- `--input-format` 缺省为 `auto`，依据首个非空行的 `>` / `@` 判定。
- `--output-format` 缺省时保持原格式（FASTA 读出即写 FASTA，FASTQ 读出即写 FASTQ）。
- FASTA 序列按 `--line-width`（正整数，默认 80）换行；该选项不影响 FASTQ，FASTQ 始终输出固定四行，质量值按 Phred+33 还原。
- FASTA 输入（无质量值）写成 FASTQ 属于用法错误。
- 输出统一使用 `\n`，末尾恰有一个换行符；相同输入与参数产生逐字节相同的结果。
- 退出码：成功 `0`；参数或格式错误 `2`（单行错误信息写到标准错误）；文件读写错误 `1`。

## 现有公开接口

- 命令行程序 `genome-variant-toolkit`
- Python 包 `genome_variant`，其 `__version__` 为当前版本号
- `genome_variant.sequence_io`：
  - `SequenceRecord`：记录（`id`、`sequence`、`description`、`quality`）
  - `read_sequences(source, format="auto")`：接受文本路径或文本流，惰性按序产生记录
  - `write_sequences(records, output, format="fasta", *, line_width=80)`
  - `SequenceFormatError`、`SequenceValidationError`

序列解析规则：

- 标识取标题行第一个非空字段，其余部分为描述。
- 序列去除行内空白后转大写，只允许 IUPAC DNA 符号 `ACGTRYSWKMBDHVN`。
- FASTA 支持多行序列；FASTQ 支持序列与质量值跨行，`+` 分隔行可只写加号或原样重复标题（重复时必须与起始标题完全一致）。
- FASTQ 质量按 Phred+33（ASCII 33-126）解码，长度必须等于序列长度。
- 错误消息包含来源名与一基行号；非法碱基错误指出记录标识、符号与一基位置。

## 限制

- 比对、变异检测与注释等后续功能尚未实现。
