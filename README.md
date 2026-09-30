## 用途

本项目是「基因组变异分析工具链」的代码仓库，用于逐步实现该方向的序列处理、比对与变异分析能力。

当前版本提供 FASTA/FASTQ 序列的读取、校验与规范化功能。

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
genome-variant-toolkit version       # 打印版本号
genome-variant-toolkit --help        # 打印用法（无子命令时同样打印帮助并返回 0）
```

### normalize-sequences

读取 FASTA 或 FASTQ 并输出规范化结果：

```bash
genome-variant-toolkit normalize-sequences INPUT \
    [--input-format fasta|fastq|auto] \
    [--output-format fasta|fastq|auto] \
    [--output OUTPUT] \
    [--line-width N]
```

- `INPUT` 与 `--output` 为 `-` 时分别表示标准输入、标准输出；默认输出到标准输出，格式自动识别并保持原格式。
- FASTA 序列按 `--line-width`（正整数，默认 60）换行；FASTQ 始终输出每条记录四行，质量按 Phred+33 还原，`--line-width` 不影响 FASTQ 布局。
- 将没有质量值的 FASTA 记录写成 FASTQ 属于用法错误。
- 相同输入与参数产生逐字节相同的输出，文本统一使用 `\n`，文件末尾恰有一个换行符。
- 返回码：成功 `0`；参数或格式错误 `2`（标准错误单行消息）；文件读写错误 `1`。

## Python 公开接口

模块 `genome_variant.sequence_io` 提供：

- `SequenceRecord(identifier, sequence, description="", quality=None)`：记录标识、大写序列、标题余部（描述），FASTQ 记录的 `quality` 为 Phred 质量整数元组，FASTA 为 `None`。
- `read_sequences(source, format="auto")`：接受文本路径或文本流，`format` 为 `fasta`、`fastq` 或 `auto`，按输入顺序惰性产生记录。
- `write_sequences(records, output, format="fasta", line_width=60)`：写 FASTA 或 FASTQ。

序列只允许 IUPAC DNA 符号 `ACGTRYSWKMBDHVN`（小写输入会转大写）；FASTQ 质量字符范围为 ASCII 33–126，质量长度必须等于序列长度。

异常：

- `SequenceFormatError`：空输入、无法判定格式、空标识、空序列、标题前出现内容、FASTQ 结构不完整或长度不等；消息包含来源名与一基行号。
- `SequenceValidationError`：非法碱基；消息指出记录标识、符号与一基位置。
- 打开路径失败保留 `OSError`，迭代期间的其他异常不会被吞掉。

## 限制

- 仅支持文本 FASTA/FASTQ 与 IUPAC DNA 符号，尚不支持比对、变异调用等后续流程。
