## 用途

本项目是「基因组变异分析工具链」的代码仓库，用于逐步实现该方向的序列处理、比对与变异分析能力。

当前版本提供 FASTA/FASTQ 序列的读取、校验与规范化、读段质量过滤、k-mer 索引、成对序列比对以及参考序列感知的 VCF 规范化功能。

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

读取 FASTA 或 FASTQ，构建 k-mer 索引并以 JSON Lines 输出：

```bash
genome-variant-toolkit kmer-index INPUT --k N \
    [--input-format fasta|fastq|auto] \
    [--no-canonical] \
    [--output OUTPUT]
```

- `--k` 为必填的正整数 k-mer 长度；`--no-canonical` 关闭正反向互补归并（默认以正向片段与其反向互补中字典序较小者为键，反向互补来源标记为 `-`，回文固定为 `+`）。
- 每个 k-mer 输出一行紧凑 JSON，固定包含 `kmer`、`count`、`occurrences`；位置对象固定包含 `record`（零基记录序号）、`id`、`position`（零基起点）、`strand`；行按 kmer 字典序排列，非 ASCII 标识原样保留。
- 含非 `A`、`C`、`G`、`T` 符号的窗口被跳过，短于 k 的记录无条目；重复标识不合并；不按质量值过滤记录。
- 输出统一使用 `\n`，非空结果末尾恰有一个换行符，空索引输出为空；`INPUT` 与 `--output` 为 `-` 时沿用标准流约定。
- 输出到文件时先完成全部读取与校验，再经同目录临时文件原子替换；失败不改变已有目标、不留下部分结果。
- 返回码：成功 `0`；参数、格式或序列错误 `2`（标准错误单行消息，且不输出数据）；文件读写错误 `1`。

### align-pair

对两个各含一条记录的 FASTA/FASTQ 输入进行成对序列比对，输出单行紧凑 JSON；FASTQ 质量值不参与打分：

```bash
genome-variant-toolkit align-pair REFERENCE QUERY \
    [--reference-format fasta|fastq|auto] \
    [--query-format fasta|fastq|auto] \
    [--mode global|local] \
    [--match N] [--mismatch N] [--gap-open N] [--gap-extend N] \
    [--output OUTPUT]
```

- `REFERENCE`、`QUERY` 与 `--output` 为 `-` 时沿用标准流约定（输出默认标准输出），但两个输入不能同时为标准输入；输入格式默认自动识别，各自必须恰好包含一条记录。
- `--match`（别名 `--match-score`，默认 2）为正整数匹配分；`--mismatch`（别名 `--mismatch-penalty`，默认 3）、`--gap-open`（默认 5）、`--gap-extend`（默认 2）为非负整数罚分；长度为 `L` 的连续缺口扣 `gap_open + gap_extend × (L-1)`。仅完全相同的字符算匹配，其他 IUPAC 组合（如 `R` 对 `A`）均算错配。
- `--mode global`（默认）覆盖两条完整序列；`--mode local` 只返回正分最佳片段，无正分片段时返回 0 分、四个坐标均为 0、空 CIGAR 与空对齐串。同分时依次选择 `reference_start`、`reference_end`、`query_start`、`query_end` 较小者，仍相同则选择 CIGAR 字典序较小者。
- 输出按固定顺序包含 `reference`、`query`（记录标识）、`mode`、`score`、`reference_start`、`reference_end`、`query_start`、`query_end`、`cigar`、`aligned_reference`、`aligned_query`；坐标零基、右端不含；CIGAR 使用 `=`、`X`、`I`（只消耗查询）、`D`（只消耗参考），相邻同类操作合并；对齐串以 `-` 表示缺口。成功输出末尾恰有一个 `\n`。
- 输出到文件时仅在全部校验与比对成功后经同目录临时文件原子替换；失败保留已有内容且无部分结果。相同输入与参数产生逐字节相同的输出。
- 返回码：成功 `0`；参数错误、空输入、多条记录、两边同时使用标准输入、格式或序列错误均为 `2`（标准错误单行，标准输出为空）；文件读写错误 `1`。

### normalize-vcf

读取 VCF，按参考 FASTA 校验并规范化后写出：

```bash
genome-variant-toolkit normalize-vcf INPUT --reference REFERENCE [--output OUTPUT]
```

- `INPUT`、`--reference` 与 `--output` 为 `-` 时沿用标准流约定（输出默认标准输出），`--reference` 为必填的参考 FASTA；VCF 输入与参考不能同时为标准输入。
- `##` 元信息行按原有内容与顺序保留；要求恰有一条合法的 `#CHROM` 列头，记录的固定列、FORMAT 与样本列数必须与列头一致；记录按输入顺序处理与写出，顺序不变。
- 普通序列型 REF/ALT 统一为大写并限于 IUPAC DNA 符号；符号 ALT（`<DEL>` 等）、断点 ALT、`*` 与 `.` 原样保留。
- 规范化按 CHROM 在参考中查找唯一序列，先验证 POS 对应的 REF 片段，再把一个记录内全部序列型等位基因作为整体移除共同的最长后缀与前缀（每个等位基因始终保留至少一个碱基，POS 同步前移）；存在长度差时继续左移到表示相同单倍型的最小 POS，每次移动后重新最简化。含任一符号、断点、`*` 或 `.` ALT 的记录不剪裁也不左移，但仍校验 REF。ALT 次序以及 ID、QUAL、FILTER、INFO、FORMAT 和样本文本保持不变。
- 缺少或重复列头、列数不符、非正整数 POS、空或非法 REF/ALT 以及越界记录为 VCF 格式错误；参考中缺少或重复 CHROM、或 REF 与参考不一致为参考不匹配错误。
- 输出使用制表符分隔与 `\n` 换行，末尾恰有一个换行符；空记录 VCF 仍写完整头部；相同输入与参考产生逐字节相同的输出。
- 输出到文件时在全部记录验证成功后经同目录临时文件原子替换；失败保留已有内容且无部分结果。
- 返回码：成功 `0`；参数错误、两边同时使用标准输入、VCF 格式错误、参考不匹配或参考序列格式错误均为 `2`（标准错误单行，标准输出为空）；文件读写错误 `1`。

## Python 公开接口

模块 `genome_variant.sequence_io` 提供：

- `SequenceRecord(identifier, sequence, description="", quality=None)`：记录标识、大写序列、标题余部（描述），FASTQ 记录的 `quality` 为 Phred 质量整数元组，FASTA 为 `None`。
- `read_sequences(source, format="auto")`：接受文本路径或文本流，`format` 为 `fasta`、`fastq` 或 `auto`，按输入顺序惰性产生记录。
- `write_sequences(records, output, format="fasta", line_width=60)`：写 FASTA 或 FASTQ。

模块 `genome_variant.quality` 提供：

- `filter_reads(records, min_end_quality=20, min_mean_quality=20, min_length=30)`：惰性的读段修剪与过滤入口，接收 `SequenceRecord` 迭代器并按输入顺序产生新记录。每条记录先从左右两端连续删除 Phred 分值严格低于末端阈值的碱基（等于阈值保留，内部低质量碱基不切分）；全部碱基被删除、修剪后长度小于最短长度、或质量分总和小于平均质量阈值乘以修剪后长度时丢弃该记录；保留记录的标识与描述不变，序列与质量同步切片。三个阈值只接受非布尔整数，质量阈值限于 0 至 93，最短长度至少为 1，非法值在调用时抛出 `ValueError`；记录错误在迭代到该记录时惰性出现。

模块 `genome_variant.kmer` 提供：

- `KmerOccurrence(record, id, position, strand)`：一次 k-mer 出现，含零基记录序号、记录标识、零基序列起点与方向（`+` 或 `-`）。
- `build_kmer_index(records, k, canonical=True)`：接收 `SequenceRecord` 可迭代对象，返回以 k-mer 为键、按键字典序排列的映射，值为按输入顺序排列的 `KmerOccurrence` 列表。`canonical` 为真时以正向片段及其反向互补中字典序较小者为键（反向互补来源标 `-`，回文固定 `+`），为假时仅索引正向片段并标 `+`；含非 `ACGT` 符号的窗口跳过，短于 k 的记录无条目，重复标识不合并。`k` 只接受至少为 1 的非布尔整数，`canonical` 只接受布尔值，非法参数在消费 `records` 前抛出 `ValueError`；记录含非法碱基时在迭代到该记录时抛出 `SequenceValidationError`。

模块 `genome_variant.alignment` 提供：

- `PairwiseAlignment`：不可变结果对象，字段依次为 `reference`、`query`（记录标识）、`mode`、`score`、`reference_start`、`reference_end`、`query_start`、`query_end`、`cigar`、`aligned_reference`、`aligned_query`；坐标零基、右端不含，对齐串等长且以 `-` 表示缺口。
- `align_pair(reference, query, mode="global", match_score=2, mismatch_penalty=3, gap_open=5, gap_extend=2)`：比对两个 `SequenceRecord`。`mode` 为 `global`（覆盖两条完整序列）或 `local`（只取正分最佳片段；无正分片段时返回 0 分、四个坐标均为 0、空 CIGAR 与空对齐串）。长度为 `L` 的连续缺口扣 `gap_open + gap_extend × (L-1)`；仅完全相同字符算匹配，其他 IUPAC 组合均算错配；质量值不参与打分。同分时依次选择四个坐标较小者，仍相同则选择 CIGAR 字典序较小者。`match_score` 必须为非布尔正整数，三个罚分必须为非布尔非负整数，`mode` 非法时在处理记录前抛出 `ValueError`；非法碱基抛出 `SequenceValidationError`。

模块 `genome_variant.vcf` 提供：

- `VcfHeader(meta_lines, columns)`：VCF 头部；`meta_lines` 为按原有内容与顺序保留的 `##` 行（含 `##` 前缀），`columns` 为 `#CHROM` 列头的全部字段（八个固定列，可选 FORMAT 与样本列）；`sample_names` 属性返回样本列名。
- `VcfRecord(chrom, pos, id=".", ref="", alts=(), qual=".", filter=".", info=".", format=None, samples=(), line_number=0)`：一条 VCF 记录；`ref` 与序列型 `alts` 为大写 IUPAC 序列，符号、断点、`*`、`.` 等位基因原样保留，`format` 为 `None` 表示无 FORMAT 列，`line_number` 为一基来源行号（写出时忽略）。
- `read_vcf(source)`：接受文本路径或文本流，返回 `(VcfHeader, list[VcfRecord])`，记录按输入顺序排列；`##` 元信息原样保留，要求恰有一条合法的 `#CHROM` 列头，记录列数必须与列头一致。
- `normalize_vcf(records, reference, source="<stream>")`：按参考校验并规范化记录，返回同序的新记录列表。`reference` 为「CHROM → 序列」映射或 `SequenceRecord` 可迭代对象；重复 CHROM 视为参考不匹配。先验证 POS 对应的 REF 片段，再把记录内全部序列型等位基因作为整体移除共同最长后缀与前缀（每个等位基因至少保留一个碱基，POS 同步调整），存在长度差时继续左移到表示相同单倍型的最小 POS 并在每次移动后重新最简化；含符号、断点、`*` 或 `.` ALT 的记录只校验不修改。`source` 用于越界错误消息中的来源名。
- `write_vcf(header, records, output)`：向文本路径或文本流写 VCF；`##` 行原样写出，记录按制表符分隔、`\n` 换行，末尾恰有一个换行符，空记录仍写完整头部。

序列只允许 IUPAC DNA 符号 `ACGTRYSWKMBDHVN`（小写输入会转大写）；FASTQ 质量字符范围为 ASCII 33–126，质量长度必须等于序列长度。

异常：

- `SequenceFormatError`：空输入、无法判定格式、空标识、空序列、标题前出现内容、FASTQ 结构不完整或长度不等；消息包含来源名与一基行号。
- `SequenceValidationError`：非法碱基；消息指出记录标识、符号与一基位置。
- `ReadQualityError`（`ValueError` 子类，位于 `genome_variant.quality`）：记录没有质量值、质量数量与序列长度不等或质量值超出 0 至 93；消息包含记录标识。
- `VcfFormatError`（位于 `genome_variant.vcf`）：缺少或重复列头、列数不符、非正整数 POS、空或非法 REF/ALT 以及越界记录；消息包含来源与一基行号。
- `ReferenceMismatchError`（位于 `genome_variant.vcf`）：参考中缺少或重复 CHROM，或 REF 与参考不一致；消息包含 CHROM、POS 和可判定的期望值与实际值，并以同名属性携带。
- 打开路径失败保留 `OSError`，迭代期间的其他异常不会被吞掉。

## 限制

- 仅支持文本 FASTA/FASTQ 与 IUPAC DNA 符号，尚不支持变异调用等后续流程。
