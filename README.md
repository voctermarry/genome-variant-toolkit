## 用途

本项目是「基因组变异分析工具链」的代码仓库，用于逐步实现该方向的序列处理、比对与变异分析能力。

当前版本提供 FASTA/FASTQ 序列的读取、校验与规范化、读段质量过滤，以及可复用的 k-mer 出现位置索引。

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

子命令 `normalize-sequences`、`filter-reads`、`kmer-index` 分别见下文。

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

### filter-reads

只接受 FASTQ 输入，按质量修剪末端并过滤读段：

```bash
genome-variant-toolkit filter-reads INPUT \
    [--min-end-quality N] \
    [--min-mean-quality N] \
    [--min-length N] \
    [--output OUTPUT]
```

- `--min-end-quality`（默认 20）：从左右两端连续删除 Phred 分值**严格低于**该阈值的碱基；等于阈值的碱基保留，内部低质量碱基不切分读段。
- `--min-mean-quality`（默认 20）：修剪后质量分总和小于「阈值 × 修剪后长度」的记录被丢弃。
- `--min-length`（默认 30，至少 1）：修剪后长度短于该值的记录被丢弃；全部碱基被修剪同样丢弃。
- 质量阈值只接受 0 至 93 的非布尔整数；非法阈值、FASTA 输入、格式或碱基错误以及记录质量错误均返回 `2`（标准错误单行消息）；文件读写失败返回 `1`。
- `INPUT` 与 `--output` 为 `-` 时沿用标准流约定，未指定 `--output` 时写标准输出；记录按输入顺序写出，四行 Phred+33 布局、`\n` 换行、文件末尾恰有一个换行符；全部记录被过滤时输出为空，成功时标准错误为空。
- 输出到文件时先写同目录临时文件，成功后原子替换；任一读取、校验、处理或写入错误都不会留下部分结果，已有目标内容保持不变。
- 相同输入与参数产生逐字节相同的输出。

### kmer-index

读取 FASTA 或 FASTQ，构建可复用的 k-mer 出现位置索引，以 JSON Lines 输出（每个 k-mer 一行）：

```bash
genome-variant-toolkit kmer-index INPUT --k N \
    [--input-format fasta|fastq|auto] \
    [--no-canonical] \
    [--output OUTPUT]
```

- `--k`（必填）：k-mer 长度，只接受至少为 1 的非布尔整数。
- 默认为规范化（canonical）模式：每个窗口取其正向片段与反向互补序列中字典序较小者作键，方向分别记为 `+` 或 `-`；二者相等（回文）时固定为 `+`。指定 `--no-canonical` 时仅用正向片段作键，方向恒为 `+`。
- 含 `A、C、G、T` 以外符号（包括 IUPAC 简并符号 `N` 等）的窗口跳过；长度短于 `k` 的记录不产生条目；重复记录标识不合并，位置按输入顺序追加，并带零基记录序号、记录标识、零基序列起点与方向。
- 每行固定包含 `kmer`、`count`、`occurrences`；位置对象固定包含 `record`、`id`、`position`、`strand`；`count` 等于位置数，行按 `kmer` 字典序排列。JSON 为紧凑形式（键值与元素间无空格），保留非 ASCII 标识（不转义）。
- 统一使用 `\n`；非空结果末尾恰有一个换行，空索引输出为空。
- `INPUT` 与 `--output` 沿用标准流约定（`-` 分别表示标准输入、标准输出）；k-mer 处理不按质量值过滤记录。
- 输出到文件时先完成全部读取与校验、再写同目录临时文件并原子替换；任一参数、格式、序列或写入错误都返回非零码且不改变目标、不留下部分结果。
- 返回码：成功 `0`；参数、格式或序列错误 `2`（标准错误单行消息且不输出数据）；文件读写错误 `1`。

## Python 公开接口

模块 `genome_variant.sequence_io` 提供：

- `SequenceRecord(identifier, sequence, description="", quality=None)`：记录标识、大写序列、标题余部（描述），FASTQ 记录的 `quality` 为 Phred 质量整数元组，FASTA 为 `None`。
- `read_sequences(source, format="auto")`：接受文本路径或文本流，`format` 为 `fasta`、`fastq` 或 `auto`，按输入顺序惰性产生记录。
- `write_sequences(records, output, format="fasta", line_width=60)`：写 FASTA 或 FASTQ。

模块 `genome_variant.quality` 提供：

- `filter_reads(records, min_end_quality=20, min_mean_quality=20, min_length=30)`：惰性的读段修剪与过滤入口，接收 `SequenceRecord` 迭代器并按输入顺序产生新记录。每条记录先从左右两端连续删除 Phred 分值严格低于末端阈值的碱基（等于阈值保留，内部低质量碱基不切分）；全部碱基被删除、修剪后长度小于最短长度、或质量分总和小于平均质量阈值乘以修剪后长度时丢弃该记录；保留记录的标识与描述不变，序列与质量同步切片。三个阈值只接受非布尔整数，质量阈值限于 0 至 93，最短长度至少为 1，非法值在调用时抛出 `ValueError`；记录错误在迭代到该记录时惰性出现。

模块 `genome_variant.kmer` 提供：

- `build_kmer_index(records, k, canonical=True)`：接收 `SequenceRecord` 可迭代对象，返回 `KmerIndex`。`k` 只接受至少为 1 的非布尔整数，`canonical` 只接受布尔值；非法参数在消费 `records` 前抛出 `ValueError`。记录按输入顺序处理：重复标识不合并，短于 `k` 的记录无条目；含非 `ACGT` 符号的窗口跳过（IUPAC 简并符号不报错），手工构造记录中含 IUPAC 字母表外符号时，在处理到该记录时抛出 `SequenceValidationError`，消息指出记录标识、符号与一基位置。规范化模式下键为正向片段与其反向互补序列中字典序较小者（回文固定 `+`），方向为 `+` 或 `-`；非规范化模式仅用正向片段且恒为 `+`。同一记录流无论上游如何分批，键与位置顺序都相同。
- `KmerIndex`：按 `kmer` 字典序迭代（`kmers()`、`items()`、`__iter__`），`occurrences(kmer)` 与 `count(kmer)` 返回该键的位置（输入顺序）。
- `KmerOccurrence`：不可变位置对象，字段为零基 `record`（记录序号）、`id`（记录标识）、零基 `position`（序列起点）与 `strand`（`+`/`-`）。

序列只允许 IUPAC DNA 符号 `ACGTRYSWKMBDHVN`（小写输入会转大写）；FASTQ 质量字符范围为 ASCII 33–126，质量长度必须等于序列长度。

异常：

- `SequenceFormatError`：空输入、无法判定格式、空标识、空序列、标题前出现内容、FASTQ 结构不完整或长度不等；消息包含来源名与一基行号。
- `SequenceValidationError`：非法碱基；消息指出记录标识、符号与一基位置。
- `ReadQualityError`（`ValueError` 子类，位于 `genome_variant.quality`）：记录没有质量值、质量数量与序列长度不等或质量值超出 0 至 93；消息包含记录标识。
- 打开路径失败保留 `OSError`，迭代期间的其他异常不会被吞掉。

## 限制

- 仅支持文本 FASTA/FASTQ 与 IUPAC DNA 符号，尚不支持比对、变异调用等后续流程。
