## 用途

本项目是「基因组变异分析工具链」的代码仓库，用于逐步实现该方向的序列处理、比对与变异分析能力。

当前版本提供 FASTA/FASTQ 序列的读取、校验与规范化、读段质量过滤、读段级确定性去重、k-mer 索引、成对序列比对、多参考多读段的确定性映射、单样本 SNV（可选短插入/短缺失）变异调用、参考序列感知的 VCF 规范化、基于 GFF3 CDS 的 VCF 变异注释与影响分级，以及把多样本单样本 VCF 清单汇总为逐变异 JSON Lines 的功能。

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

### deduplicate-reads

在映射与变异调用前对 FASTA 或 FASTQ 读段做确定性去重，输出格式与输入一致：

```bash
genome-variant-toolkit deduplicate-reads INPUT \
    [--input-format fasta|fastq|auto] \
    [--canonical] \
    [--output OUTPUT]
```

- 默认仅把完整大写序列完全相同的记录归为一组；标识相同但序列不同的记录不合并。
- `--canonical`：以序列与其完整 IUPAC 反向互补序列中字典序较小者为分组键，因此两种方向也视为重复；不启用时正反向序列保持独立。
- 每组只输出一条代表记录：FASTA 组保留最早出现的记录；FASTQ 组保留 Phred 质量总和最高的记录，分数相同时保留最早出现者。代表记录的标识、描述、原始方向、序列与质量值均原样保留，不把序列改写成分组键。
- 结果按各分组首次出现的输入位置排序，输入迭代器采用不同分块方式时结果不变；空的 Python 迭代器返回空元组，命令行空文件仍沿用现有读取规则视为格式错误。
- `INPUT` 与 `--output` 为 `-` 时沿用标准流约定，`--input-format` 默认自动识别；FASTQ 始终输出每条记录四行，FASTA 使用默认行宽（60）；统一使用 `\n`，非空结果末尾恰有一个换行符；空文件为格式错误而非空输出。
- 输出到文件时仅在全部记录完成校验和去重后经同目录临时文件原子替换；任何参数、格式、序列、质量错误或读写失败都保留已有目标且不留下部分结果。相同输入与参数产生逐字节相同的输出。
- 返回码：成功 `0`；参数、格式、序列或质量错误，以及同一输入混有带质量与不带质量记录均为 `2`（标准错误单行、标准输出为空）；文件读写错误 `1`。

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

### map-reads

把多条 FASTA/FASTQ 读段确定性地映射到多记录 FASTA 参考上，按读段输入顺序输出 JSON Lines，每条一行；FASTQ 质量值不参与打分：

```bash
genome-variant-toolkit map-reads REFERENCE READS \
    [--reads-format fasta|fastq|auto] \
    [--match N] [--mismatch N] [--gap-open N] [--gap-extend N] \
    [--min-score N] \
    [--output OUTPUT]
```

- `REFERENCE` 按 FASTA 解析，可含多条参考记录；`READS` 为 FASTA 或 FASTQ，默认自动识别（可用 `--reads-format` 指定）。两者与 `--output` 为 `-` 时沿用标准流约定（输出默认标准输出），但两个输入不能同时为标准输入。
- 每条读段分别以原序列（正链）与完整 IUPAC 反向互补序列（负链），与每条参考记录做局部比对，计分参数与 `align-pair` 的局部比对一致（默认 `--match 2`、`--mismatch 3`、`--gap-open 5`、`--gap-extend 2`）；`--min-score` 为默认值 1 的非布尔正整数，只保留分数不低于该值的候选。
- 候选先按分数从高到低选择；同分时依次取参考输入序号、`reference_start`、`reference_end`、正链优先、换算回原读段后的 `query_start`、`query_end`、CIGAR 字典序较小者。负链的查询坐标为换算到原读段后的零基半开区间，CIGAR 仍按参考与反向互补读段的方向解释；重复标识不合并。
- 无达标候选的读段仍输出一条未映射记录。每行固定按顺序包含 `record`（零基读段序号）、`id`、`mapped`、`reference_record`（零基参考序号）、`reference`、`reference_start`、`reference_end`、`query_start`、`query_end`、`strand`（`+` 或 `-`）、`score`、`cigar`；未映射记录仅 `record`、`id` 与 `mapped:false` 有值，其余字段均为 `null`。
- 输出统一使用 `\n`，非空结果末尾恰有一个换行符（每条读段恰好对应一行，故成功输出不会为空）；写文件时仅在全部读段成功处理后经同目录临时文件原子替换，任一读取、校验或处理失败都保留已有目标且不留部分结果。相同输入与参数无论如何分批都产生逐字节相同的结果。
- 返回码：成功 `0`；参数错误、空参考/空读段输入、两边同时使用标准输入、格式或序列错误均为 `2`（标准错误单行、标准输出为空）；文件读写错误 `1`。

### call-variants

从多记录参考 FASTA 与 FASTQ 读段调用单样本 SNV，输出可被现有 VCF 接口读回的 VCFv4.2 文档；映射沿用 `map-reads` 的局部比对计分、正反链搜索与候选裁决，质量不参与映射打分、只过滤碱基证据：

```bash
genome-variant-toolkit call-variants REFERENCE READS \
    [--reads-format fasta|fastq|auto] \
    [--min-base-quality N] [--min-alt-count N] \
    [--min-alt-fraction F] [--homozygous-fraction F] \
    [--sample-name NAME] \
    [--call-indels] [--max-indel-length N] \
    [--output OUTPUT]
```

- `REFERENCE` 按 FASTA 解析且可含多条参考记录（标识不得重复，集合不得为空）；`READS` 默认自动识别（可用 `--reads-format` 指定），必须携带质量值。两者与 `--output` 为 `-` 时沿用标准流约定（输出默认标准输出），但两个输入不能同时为标准输入。
- 默认仅调用 `A`、`C`、`G`、`T` 之间的替换；插入、缺失、未对齐部分、未映射读段与歧义碱基（参考或观测碱基非 `ACGT`）均忽略。每条已映射读段在同一参考位置最多贡献一次观测；负链证据按参考正链方向报告观测碱基，质量值对应回原读段位置。
- `--call-indels`（默认关闭）：额外调用短插入与短缺失。只从胜出比对的连续 `I`、`D` CIGAR 片段提取事件；插入以左侧参考碱基加插入序列表示，缺失以左锚点加被删除参考片段表示。无可用左锚点、缺少左右任一已对齐读段碱基、任一相关参考或插入碱基非 `ACGT`、或长度超过 `--max-indel-length` 的事件一律忽略。插入证据要求两侧及全部插入碱基质量达标，缺失证据要求两侧碱基质量达标（负链质量换算回原读段坐标）。候选先按参考感知规则最简化并左对齐，再按 `CHROM`、`POS`、`REF`、`ALT` 合并等价事件；事件的 `DP` 为同时跨过其左右边界、两侧质量达标且未在该边界产生其他插入或缺失的已映射读段数，`AC` 为其中支持该事件的读段数，其余阈值与 SNV 相同，`AD` 为 `DP-AC,AC`。
- `--max-indel-length`（默认 50，正整数）：可调用的最长插入/缺失碱基数；必须与 `--call-indels` 同用，单独给出或数值非法均为参数错误。
- `--min-base-quality`（默认 20，0 至 93 的非布尔整数）：只接纳 Phred 不低于该值的碱基；`DP` 为通过门槛的 `ACGT` 观测总数。
- `--min-alt-count`（默认 2，正整数）与 `--min-alt-fraction`（默认 0.2，0 至 1 的有限数）：非参考碱基中计数最高者为唯一 ALT，并列时按 `A`、`C`、`G`、`T` 选择；ALT 计数达到门槛且 `AC/DP` 不低于最小比例才输出。
- `--homozygous-fraction`（默认 0.8，0 至 1 的有限数，且不得低于 `--min-alt-fraction`）：ALT 比例不低于该值时 GT 为 `1/1`，否则为 `0/1`。
- 输出记录按参考输入次序、再按 `POS`、`REF`、`ALT` 升序排列（同一位置的多个合格 indel 分别成行）；`REF`、`ALT` 大写，`ID` 为 `.`，`QUAL` 为点，`FILTER` 为 `PASS`；`INFO` 为 `DP`、`AC`、`AF`（`AF` 固定六位小数），`FORMAT` 为 `GT:DP:AD`，`AD` 为参考与所选 ALT 的计数，样本名默认 `SAMPLE`（`--sample-name` 覆盖，不得为空、纯空白或含制表符）。
- 无候选时仍输出完整头部（7 行 `##` 元信息加 `#CHROM` 列头）。写文件时仅在全部读取与调用成功后经同目录临时文件原子替换，失败保留原目标；读段分批不影响逐字节结果。
- 返回码：成功 `0`；参数错误、非法阈值（在消费输入前报错）、空参考、重复参考标识、非法样本名、缺少质量或质量长度错误、FASTQ 与序列错误均为 `2`（标准错误单行、标准输出为空）；文件读写错误 `1`。

### batch-call-variants

按 JSON Lines 清单对多个样本的 FASTQ 读段逐一调用变异并直接合并汇总，不产生中间 VCF：

```bash
genome-variant-toolkit batch-call-variants MANIFEST \
    --reference REFERENCE \
    [--min-base-quality N] [--min-alt-count N] \
    [--min-alt-fraction F] [--homozygous-fraction F] \
    [--call-indels] [--max-indel-length N] \
    [--output OUTPUT]
```

- 清单为 JSON Lines：每个非空行是仅含字符串字段 `sample` 与 `reads` 的 JSON 对象；`sample` 不得为空、纯空白或含制表符，且在清单内唯一；`reads` 为 FASTQ 路径且不得为 `-`。文件清单中的相对 `reads` 按清单所在目录解析，流清单（`MANIFEST` 为 `-`）按当前工作目录解析。
- 每个清单项按原顺序读取读段，以 `sample` 为样本名，按 `call-variants` 的质量门槛、等位计数、等位比例、基因型及可选短插入缺失语义独立调用；各阈值选项与 `--call-indels`、`--max-indel-length` 的默认值、约束和结果与逐样本调用一致。
- 各样本结果按 `summarize-variants` 的参考规范化规则合并，输出与其相同的紧凑 JSON Lines：字段、统计口径与排序一致，`samples` 保持清单顺序，没有变异时输出为空。相同参考、清单、读段与参数跨运行产生逐字节相同的结果。
- `MANIFEST`、`--reference` 与 `--output` 沿用标准流约定，但清单与参考不能同时来自标准输入。
- 输出到文件时仅在全部样本调用与汇总成功后经同目录临时文件原子替换；失败保留已有目标且不留临时结果。
- 返回码：成功 `0`；参数、清单、序列、质量或参考数据错误 `2`（标准错误单行、标准输出为空）；文件访问或写入失败 `1`。

### coverage-report

统计多记录参考 FASTA 的覆盖度与一致性，按参考输入顺序输出 JSON Lines，每条参考记录一行；映射沿用 `map-reads` 的局部比对计分、正反链搜索与候选裁决，FASTQ 质量值不参与统计：

```bash
genome-variant-toolkit coverage-report REFERENCE READS \
    [--reads-format fasta|fastq|auto] \
    [--match N] [--mismatch N] [--gap-open N] [--gap-extend N] \
    [--min-score N] \
    [--output OUTPUT]
```

- 输入、映射选项（`--match`、`--mismatch`、`--gap-open`、`--gap-extend`、`--min-score`）、`--output` 与标准流约定均同 `map-reads`；两个输入不能同时为标准输入。
- 每条读段只采用唯一胜出的映射，未映射读段不计入。CIGAR 的 `=` 与 `X` 各为对应参考位置增加一次深度，`D` 不增加，`I` 无参考位置；一致仅指观测字符与参考字符逐字相同（包括 IUPAC 字符，如 `N` 对 `N` 一致、`R` 对 `A` 不一致）。
- 每条参考记录输出一行，无命中也输出，重复标识不合并、以零基 `reference_record` 区分。每行固定按顺序包含 `reference_record`、`reference`、`length`、`mapped_reads`、`observed_bases`（深度总和）、`covered_bases`（深度非零的位置数）、`coverage_fraction`、`mean_depth`、`concordant_bases`、`concordance_fraction`；三个小数分别按 `covered_bases/length`、`observed_bases/length`、`concordant_bases/observed_bases` 计算，分母为零取 `0`，均为固定六位小数字符串。
- 成功输出为 UTF-8 与 `\n`，末尾恰有一个换行；写文件时仅在全部校验与统计完成后经同目录临时文件原子替换，失败保留已有目标。相同输入与参数无论如何分批都产生逐字节相同的结果。
- 返回码：成功 `0`；参数错误、空参考/空读段输入、两边同时使用标准输入、格式或序列错误均为 `2`（标准错误单行、标准输出为空）；文件读写错误 `1`。

### normalize-vcf

读取 VCF 并结合必需的参考 FASTA 进行参考序列感知的规范化（等位基因最简化与左对齐）：

```bash
genome-variant-toolkit normalize-vcf INPUT --reference REFERENCE \
    [--output OUTPUT]
```

- `INPUT`、`--reference` 与 `--output` 为 `-` 时沿用标准流约定（输出默认标准输出），但 VCF 与参考不能同时使用标准输入；参考按 FASTA 解析。
- 读取时保留全部 `##` 元信息行的原有内容与顺序，要求恰有一条合法的 `#CHROM` 列头，记录按输入顺序产生；固定列、FORMAT 与样本列数必须与列头一致。
- 普通序列型 REF/ALT 统一为大写并只允许 IUPAC DNA 符号 `ACGTRYSWKMBDHVN`。
- 规范化时按 CHROM 在参考中查找唯一记录：先校验 POS 处的 REF 片段，再把一个记录内全部序列型等位基因作为整体移除共同的最长后缀与前缀（始终为每个等位基因保留至少一个碱基，POS 随前缀裁剪同步前移）；等位基因长度不同时继续向左移动到在参考上表示相同单倍型的最小 POS，每次移动后重新最简化。ALT 次序以及 ID、QUAL、FILTER、INFO、FORMAT 与样本文本保持不变。
- 含任一符号型 ALT（`<...>`）、断点 ALT（含 `[`、`]` 或 `.` 连接）、`*` 或缺失值 `.` 的记录不剪裁也不左移，但仍校验 REF；记录顺序不变。
- 缺少或重复 `#CHROM` 列头、列数不符、非正整数 POS、空或非法 REF/ALT、记录越界统一为 `VcfFormatError`（消息含来源与一基行号）；参考中缺少或重复 CHROM、REF 与参考不一致为 `ReferenceMismatchError`（消息含 CHROM、POS 及可判定的期望值与实际值）。
- 成功输出使用制表符与 `\n`，末尾恰有一个换行；无记录的 VCF 仍写出完整头部。写文件时先完成全部记录的读取与校验，再经同目录临时文件原子替换，失败保留原文件。相同输入与参考产生逐字节一致的结果。
- 返回码：成功 `0`；参数或格式错误、`VcfFormatError`、`ReferenceMismatchError` 以及参考 FASTA 的序列格式错误均为 `2`（标准错误单行、标准输出为空）；文件读写错误 `1`。

### annotate-vcf

读取现有 VCF、参考 FASTA 与 GFF3 注释特征，在 INFO 列追加 `GVANN` 注释后输出 VCF：

```bash
genome-variant-toolkit annotate-vcf INPUT --reference REFERENCE --features FEATURES \
    [--output OUTPUT]
```

- `INPUT`、`--reference`、`--features` 与 `--output` 为 `-` 时沿用标准流约定（输出默认标准输出），但三个输入至多一个来自标准输入；参考按 FASTA 解析，特征按 GFF3 解析。
- 仅以 GFF3 的 `CDS` 特征为注释范围；同一 `Parent` 的多个 CDS 片段按坐标与链组合为一条转录本，片段顺序与 GFF3 行序无关。转录本命中按 `Parent` 字典序排列。
- 依据链与每段 `phase` 建立阅读框：正链按坐标升序读取，负链反向读取并对密码子取反向互补；某片段声明的 `phase` 必须与已拼装的阅读框一致，否则为格式错误。
- 仅判定 `A`、`C`、`G`、`T` 单碱基替换，后果限定为 `START_LOST`、`STOP_GAINED`、`STOP_LOST`、`SYNONYMOUS`、`MISSENSE`，影响依次为 `HIGH`、`HIGH`、`HIGH`、`LOW`、`MODERATE`（使用标准遗传密码，起始密码子按 `ATG`/甲硫氨酸判定，终止密码子记为 `*`）。未命中任何 CDS 的 SNV 输出 `NON_CODING`/`MODIFIER`；非 SNV、歧义（IUPAC 兼并码）或符号/缺失等位基因输出 `UNSUPPORTED`/`MODIFIER`，原记录一律保留。
- 元信息末尾追加唯一的 `##INFO=<ID=GVANN,...>` 定义；INFO 项目格式为 `ALT|CONSEQUENCE|IMPACT|TRANSCRIPT|CDS_POS|CODON_CHANGE|AA_CHANGE`，同一 ALT 命中多条转录本时以逗号分隔，多个 ALT 各自产生一项并保持 ALT 次序；`CDS_POS` 为编码序列一基坐标，密码子按编码链给出（如 `ATG>ACG`），氨基酸用单字母（终止为 `*`，如 `K>*`）；无命中时转录本与后续字段为 `.`。INFO 原本缺失（`.`）时以 `GVANN` 开始，否则追加到末尾。
- 记录顺序、ALT 顺序、样本列与既有字段保持不变；注释前先校验每条记录的 REF。CDS 行列字段数、坐标或 `phase` 非法、缺少 `Parent`、链非 `+`/`-`、同一 `Parent` 的片段跨序列或阅读框矛盾、CDS 越过序列末端均抛出 `AnnotationFormatError`（消息含来源与一基行号）；VCF 已声明或已使用 `GVANN` 同样抛出该异常。参考中缺少或重复 CHROM、REF 不匹配为 `ReferenceMismatchError`，记录越界为 `VcfFormatError`。
- 成功输出使用制表符与 `\n`，末尾恰有一个换行；无记录的 VCF 仍写出含 `GVANN` 定义的完整头部。写文件时先完成全部读取、校验与注释，再经同目录临时文件原子替换，失败保留原文件。相同输入跨批次逐字节一致。
- 返回码：成功 `0`；参数或数据错误（`AnnotationFormatError`、`VcfFormatError`、`ReferenceMismatchError`、参考 FASTA 序列错误）均为 `2`（标准错误单行、标准输出为空）；文件读写错误 `1`。

### summarize-variants

读取 JSON Lines 清单中列出的多个单样本 VCF，结合参考 FASTA 规范化并按变异合并，输出逐变异 JSON Lines：

```bash
genome-variant-toolkit summarize-variants MANIFEST --reference REFERENCE \
    [--output OUTPUT]
```

- 清单每个非空行为一个 JSON 对象，且仅含字符串字段 `sample` 与 `vcf`；`sample` 非空、不含制表符且在清单内全局唯一，`vcf` 指向单样本 VCF。空行忽略。文件清单中的相对 `vcf` 路径按清单文件所在目录解析，标准输入清单中的相对路径按当前工作目录解析；绝对路径原样使用。
- `MANIFEST`、`--reference` 与 `--output` 为 `-` 时沿用标准流约定（输出默认标准输出），但清单与参考不能同时使用标准输入；参考按 FASTA 解析。
- 每个 VCF 恰有一个样本且样本名与对应清单 `sample` 相同。每条记录仅接受一个普通序列型 ALT（拒绝符号型、断点、`*` 与缺失值 `.`，以及多等位），FILTER 只能为 `PASS` 或 `.`，FORMAT 必须包含 `GT`、`DP`、`AD`（允许额外字段与任意字段次序，但字段不得重复，样本列字段数须与 FORMAT 一致）。
- GT 仅接受 `0/1` 或 `1/1`；DP 为非负整数；AD 恰含两个非负整数，且二者之和等于 DP。记录先按现有规则做参考校验、最简化与左对齐，再以规范化后的 `CHROM`、`POS`、`REF`、`ALT` 合并等价调用；同一个 VCF 规范化后键重复视为数据错误。
- 每行固定按顺序包含 `chrom`、`pos`、`ref`、`alt`、`sample_count`、`allele_count`、`depth`、`samples`；`samples` 仅收录实际调用该变异的样本并保持清单次序，每项含 `sample`、`gt`、`dp`、`ad`（两个整数的数组）。`sample_count` 为样本项数，`allele_count` 对 `0/1` 计 1、对 `1/1` 计 2，`depth` 为各样本 DP 之和。
- 结果按参考记录次序，再按 `POS`、`REF`、`ALT` 升序排列；清单为空或无任何变异时输出为空。成功结果为 UTF-8 紧凑 JSON，非 ASCII 原样保留，统一使用 `\n`，非空输出末尾恰有一个换行符。
- 写文件时先完成清单、全部 VCF 与参考的读取校验及合并，再经同目录临时文件原子替换；任一失败都保留已有目标且不产生部分标准输出。相同输入跨运行逐字节一致。
- 返回码：成功 `0`；清单 JSON/字段/样本唯一性错误、VCF 结构或基因型字段不合规（消息含清单行号或 VCF 来源与一基行号）、`VcfFormatError`、`ReferenceMismatchError`、参考 FASTA 错误均为 `2`（标准错误单行、标准输出为空）；文件读写失败 `1`。

## Python 公开接口

模块 `genome_variant.sequence_io` 提供：

- `SequenceRecord(identifier, sequence, description="", quality=None)`：记录标识、大写序列、标题余部（描述），FASTQ 记录的 `quality` 为 Phred 质量整数元组，FASTA 为 `None`。
- `read_sequences(source, format="auto")`：接受文本路径或文本流，`format` 为 `fasta`、`fastq` 或 `auto`，按输入顺序惰性产生记录。
- `write_sequences(records, output, format="fasta", line_width=60)`：写 FASTA 或 FASTQ。

模块 `genome_variant.quality` 提供：

- `filter_reads(records, min_end_quality=20, min_mean_quality=20, min_length=30)`：惰性的读段修剪与过滤入口，接收 `SequenceRecord` 迭代器并按输入顺序产生新记录。每条记录先从左右两端连续删除 Phred 分值严格低于末端阈值的碱基（等于阈值保留，内部低质量碱基不切分）；全部碱基被删除、修剪后长度小于最短长度、或质量分总和小于平均质量阈值乘以修剪后长度时丢弃该记录；保留记录的标识与描述不变，序列与质量同步切片。三个阈值只接受非布尔整数，质量阈值限于 0 至 93，最短长度至少为 1，非法值在调用时抛出 `ValueError`；记录错误在迭代到该记录时惰性出现。

模块 `genome_variant.deduplication` 提供：

- `deduplicate_reads(records, canonical=False)`：接收 `SequenceRecord` 可迭代对象，确定性去重后按各分组首次出现顺序返回代表记录元组。默认仅完整大写序列完全相同的记录归为一组，标识相同但序列不同不合并；`canonical` 为真时以序列与其完整 IUPAC 反向互补序列中字典序较小者为分组键，两种方向视为重复。每个分组只输出一条代表记录：不带质量（FASTA）的组保留最早出现的记录，带质量（FASTQ）的组保留 Phred 质量总和最高者、同分时保留最早者；代表记录的标识、描述、原始方向、序列与质量值原样保留，不把序列改写成分组键。结果与迭代器分块方式无关；空迭代器返回空元组。`canonical` 不是布尔值时在消费输入前抛出 `ValueError`；非法 IUPAC 碱基抛出 `SequenceValidationError`；带质量记录的质量缺失、数量不符或数值越界抛出 `ReadQualityError`；同一迭代器混有带质量与不带质量记录时抛出 `ValueError`。

模块 `genome_variant.kmer` 提供：

- `KmerOccurrence(record, id, position, strand)`：一次 k-mer 出现，含零基记录序号、记录标识、零基序列起点与方向（`+` 或 `-`）。
- `build_kmer_index(records, k, canonical=True)`：接收 `SequenceRecord` 可迭代对象，返回以 k-mer 为键、按键字典序排列的映射，值为按输入顺序排列的 `KmerOccurrence` 列表。`canonical` 为真时以正向片段及其反向互补中字典序较小者为键（反向互补来源标 `-`，回文固定 `+`），为假时仅索引正向片段并标 `+`；含非 `ACGT` 符号的窗口跳过，短于 k 的记录无条目，重复标识不合并。`k` 只接受至少为 1 的非布尔整数，`canonical` 只接受布尔值，非法参数在消费 `records` 前抛出 `ValueError`；记录含非法碱基时在迭代到该记录时抛出 `SequenceValidationError`。

模块 `genome_variant.alignment` 提供：

- `PairwiseAlignment`：不可变结果对象，字段依次为 `reference`、`query`（记录标识）、`mode`、`score`、`reference_start`、`reference_end`、`query_start`、`query_end`、`cigar`、`aligned_reference`、`aligned_query`；坐标零基、右端不含，对齐串等长且以 `-` 表示缺口。
- `align_pair(reference, query, mode="global", match_score=2, mismatch_penalty=3, gap_open=5, gap_extend=2)`：比对两个 `SequenceRecord`。`mode` 为 `global`（覆盖两条完整序列）或 `local`（只取正分最佳片段；无正分片段时返回 0 分、四个坐标均为 0、空 CIGAR 与空对齐串）。长度为 `L` 的连续缺口扣 `gap_open + gap_extend × (L-1)`；仅完全相同字符算匹配，其他 IUPAC 组合均算错配；质量值不参与打分。同分时依次选择四个坐标较小者，仍相同则选择 CIGAR 字典序较小者。`match_score` 必须为非布尔正整数，三个罚分必须为非布尔非负整数，`mode` 非法时在处理记录前抛出 `ValueError`；非法碱基抛出 `SequenceValidationError`。

模块 `genome_variant.mapping` 提供：

- `MappingReferenceError`（`ValueError` 子类）：参考记录集合为空时抛出。
- `ReadMapping`：不可变结果对象，字段依次为 `record`（零基读段序号）、`id`、`mapped`、`reference_record`（零基参考输入序号）、`reference`（参考标识）、`reference_start`、`reference_end`、`query_start`、`query_end`、`strand`（`+` 或 `-`）、`score`、`cigar`；坐标零基、右端不含；未映射时除 `record`、`id`、`mapped` 外其余字段均为 `None`。
- `map_reads(references, reads, match_score=2, mismatch_penalty=3, gap_open=5, gap_extend=2, min_score=1)`：把读段惰性映射到多条参考，按读段输入顺序产生 `ReadMapping`。每条读段分别以原序列与完整 IUPAC 反向互补序列与每条参考做局部比对（计分语义同 `align_pair` 的 `local`），仅保留分数不低于 `min_score` 的候选；同分依次按参考输入序号、`reference_start`、`reference_end`、正链优先、原读段上的 `query_start`、`query_end`、CIGAR 字典序裁决。负链查询坐标换算回原读段，CIGAR 仍按参考与反向互补读段解释；无达标候选时仍产生未映射结果，重复标识不合并。`min_score` 必须为非布尔正整数，计分参数约束同 `align_pair`，非法参数在消费任一输入前抛出 `ValueError`；参考集合为空在消费 `reads` 前抛出 `MappingReferenceError`；非法碱基在迭代到该记录时抛出 `SequenceValidationError`。

模块 `genome_variant.calling` 提供：

- `VariantCallingError`（`ValueError` 子类）：参考集合为空、参考标识重复或样本名为空、纯空白、含制表符（或非字符串）时抛出。
- `call_variants(references, reads, *, min_base_quality=20, min_alt_count=2, min_alt_fraction=0.2, homozygous_fraction=0.8, sample_name="SAMPLE", call_indels=False, max_indel_length=50)`：返回可由 `read_vcf` 读回的 `VcfFile`（VCFv4.2，单样本）。映射与 `map_reads` 使用相同的局部比对计分（固定默认值）、正反链搜索与候选裁决；质量不参与映射，只过滤碱基证据。仅统计参考碱基与观测碱基都属于 `ACGT` 的比对列，插入、缺失、未对齐部分、未映射读段与歧义碱基忽略；每条已映射读段在同一参考位置至多贡献一次观测，负链观测碱基按参考正链报告而质量取自原读段对应位置。`DP` 为通过质量门槛的 `ACGT` 观测总数，唯一 ALT 为非参考碱基中计数最高者（并列按 `A`、`C`、`G`、`T`），ALT 计数至少 `min_alt_count` 且 `AC/DP` 至少 `min_alt_fraction` 才输出记录；比例不低于 `homozygous_fraction` 时 GT 为 `1/1`，否则 `0/1`。`call_indels` 为真时额外从胜出比对的连续 `I`、`D` CIGAR 片段调用短插入/短缺失：插入以左侧参考碱基加插入序列、缺失以左锚点加被删除参考片段表示；无左锚点、缺少任一侧已对齐读段碱基、相关参考或插入碱基非 `ACGT`、或长度超过 `max_indel_length` 的事件忽略；插入证据要求两侧及全部插入碱基质量达标，缺失证据要求两侧碱基质量达标。候选经参考感知的最简化与左对齐后按 `CHROM`、`POS`、`REF`、`ALT` 合并；事件 `DP` 为同时跨过左右边界、两侧质量达标且未在该边界产生其他插入或缺失的已映射读段数，`AC` 为其中支持该事件的读段数，`AD` 为 `DP-AC,AC`。记录按参考输入次序、`POS`、`REF`、`ALT` 排列，`INFO` 为 `DP`、`AC`、`AF`（六位小数），`FORMAT` 为 `GT:DP:AD`。`min_base_quality` 限于 0 至 93 的非布尔整数，`min_alt_count` 与 `max_indel_length` 为非布尔正整数，两个比例为 0 至 1 的有限数且纯合比例不得低于最小 ALT 比例；非法阈值在消费任一输入前抛出 `ValueError`，样本名与参考错误抛出 `VariantCallingError`，读段缺少质量、质量长度错误或质量越界在迭代到该读段时抛出 `ReadQualityError`。

模块 `genome_variant.coverage` 提供：

- `ReferenceCoverage`：不可变结果对象，字段依次为 `reference_record`（零基参考输入序号）、`reference`、`length`、`mapped_reads`、`observed_bases`、`covered_bases`、`coverage_fraction`、`mean_depth`、`concordant_bases`、`concordance_fraction`；三个小数字段为固定六位小数字符串，分母为零时为 `"0.000000"`。
- `coverage_report(references, reads, match_score=2, mismatch_penalty=3, gap_open=5, gap_extend=2, min_score=1)`：按 `map_reads` 的映射语义（相同的局部比对计分、正反链搜索与候选裁决）聚合每条参考记录的覆盖度与一致性，按参考输入顺序产生 `ReferenceCoverage`，无命中的参考也产生结果，重复标识不合并。每条读段只采用唯一胜出的映射；`=`/`X` 列各为对应参考位置增加一次深度，`D` 不增加、`I` 无参考位置；一致仅指字符逐字相同（含 IUPAC 字符）；质量值不参与。参数校验与异常时序同 `map_reads`：非法参数在消费任一输入前抛出 `ValueError`，空参考在消费 `reads` 前抛出 `MappingReferenceError`，非法碱基在迭代到该记录时抛出 `SequenceValidationError`。

模块 `genome_variant.vcf` 提供：

- `VcfHeader(meta_lines=(), samples=())`：`meta_lines` 为按输入顺序保留的 `##` 元信息行（不含行终止符），`samples` 为 `#CHROM` 列头声明的样本名元组；属性 `column_names` 给出完整列头。
- `VcfRecord`：不可变数据对象，字段为 `chrom`、`pos`（一基）、`id`、`ref`、`alt`（ALT 等位基因字符串元组，保持次序）、`qual`、`filter`、`info`，带样本时另有 `format_text` 与逐样本的 `sample_text`（均为原文），`line_number` 记录来源行号。
- `VcfFile(header, records=(), source="")`：头部与按输入顺序排列的记录元组，以及来源标签。
- `read_vcf(source)`：接受文本路径或文本流，读取并校验整个 VCF，返回 `VcfFile`；普通序列型 REF/ALT 转大写并限于 IUPAC DNA 符号，符号型（`<...>`）、断点（`[`、`]`、`.` 连接）、`*` 与缺失值 `.` ALT 原样保留。
- `normalize_record(record, reference)` 与 `normalize_vcf(document, reference)`：`reference` 为 CHROM 到参考序列的映射或 `SequenceRecord` 可迭代对象。先校验 POS 处 REF；再将全部序列型等位基因作为整体移除共同最长后缀与前缀（每等位至少留一个碱基，POS 随前缀调整），长度不同时左移到表示相同单倍型的最小 POS，每次移动后重新最简化。特殊 ALT 记录仅校验 REF，不剪裁不左移；ALT 次序及 ID、QUAL、FILTER、INFO、FORMAT、样本文本与记录顺序不变。
- `render_vcf(document)` 与 `write_vcf(document, output)`：使用制表符与 `\n` 序列化，末尾恰有一个换行；无记录时仍写出完整头部。`write_vcf` 接受文本路径或文本流。

模块 `genome_variant.annotation` 提供：

- `AnnotationFormatError`（`ValueError` 子类）：GFF3 的 CDS 行列字段数、坐标或 `phase` 非法、缺少 `Parent`、链非 `+`/`-`、同一 `Parent` 的片段跨序列或阅读框矛盾、CDS 越过序列末端，或 VCF 已声明/已使用 `GVANN` 时抛出；消息含来源名与一基行号。
- `annotate_vcf(document, reference, features)`：`document` 为 `VcfFile`，`reference` 为 CHROM 到序列的映射或 `SequenceRecord` 可迭代对象，`features` 为 GFF3 文本路径或文本流。以 GFF3 的 CDS 为注释范围，按 `Parent` 组合转录本片段、依据链与 `phase` 建立阅读框，用标准遗传密码判定 `A/C/G/T` 单碱基替换。每个 ALT 依次注释，多转录本命中按 `Parent` 字典序；后果为 `START_LOST`、`STOP_GAINED`、`STOP_LOST`（均 `HIGH`）、`SYNONYMOUS`（`LOW`）、`MISSENSE`（`MODERATE`），未命中 CDS 的 SNV 为 `NON_CODING`/`MODIFIER`，非 SNV、歧义或符号等位基因为 `UNSUPPORTED`/`MODIFIER`。返回追加唯一 `##INFO=<ID=GVANN,...>` 定义、并在每条记录 INFO 末尾（INFO 缺失时以其开始）追加 `GVANN=ALT|CONSEQUENCE|IMPACT|TRANSCRIPT|CDS_POS|CODON_CHANGE|AA_CHANGE` 的新 `VcfFile`；记录顺序、ALT 次序、样本列与既有字段不变。参考缺少/重复 CHROM 或 REF 不匹配抛 `ReferenceMismatchError`，记录越界抛 `VcfFormatError`。

模块 `genome_variant.summary` 提供：

- `ManifestError`（`ValueError` 子类）：清单非合法 JSON、行不是对象、字段缺失/多余/类型非法、`sample` 为空或含制表符、样本名重复，或样本 VCF 结构与基因型字段不合规（样本数不为一、样本名不符、ALT 非单一普通序列型、FILTER 非 `PASS`/`.`、FORMAT 缺 `GT`/`DP`/`AD`、GT 非 `0/1`/`1/1`、DP/AD 非法或 AD 之和不等于 DP、同一 VCF 规范化后键重复）时抛出；消息含清单文件与行号，或 VCF 来源与一基记录行号。
- `ManifestEntry(line_number, sample, vcf)`：一条清单项，含一基清单行号、唯一样本名与已解析的 VCF 路径。
- `VariantSummary`：一条合并结果，字段为 `ref_index`（CHROM 在参考中的零基序号）、`chrom`、`pos`、`ref`、`alt`、`samples`；属性 `sample_count`、`allele_count`（`0/1` 计 1、`1/1` 计 2）、`depth`（DP 之和）与排序键 `ref_order_key`。`SampleCall(sample, gt, dp, ad)` 为其中的单个样本调用。
- `read_manifest(source)`：接受清单文本路径或文本流，忽略空行，返回按清单次序的 `ManifestEntry` 元组；文件清单的相对 `vcf` 按其目录解析，流清单按当前工作目录解析。
- `summarize_variants(manifest, reference)`：`manifest` 为清单路径/流或 `ManifestEntry` 序列，`reference` 为 FASTA 路径/流或 `SequenceRecord` 可迭代对象。逐个读取、校验并以参考规范化每个单样本 VCF，按规范化后的 `CHROM`/`POS`/`REF`/`ALT` 跨样本合并等价调用；返回按参考次序、`POS`、`REF`、`ALT` 排序的 `VariantSummary` 元组，每项的 `samples` 保持清单次序。参考问题抛 `ReferenceMismatchError`/`VcfFormatError`，文件访问失败抛 `OSError`。
- `render_summaries(summaries)`：序列化为紧凑 UTF-8 JSON Lines（非 ASCII 原样、`\n` 换行），非空结果末尾恰有一个换行，空结果为空字符串。

模块 `genome_variant.batch` 提供：

- `BatchManifestError`（`ValueError` 子类）：清单非合法 JSON、行不是对象、字段缺失/多余/类型非法、`sample` 为空、纯空白或含制表符、样本名重复，或 `reads` 为空、非字符串或为 `-` 时抛出；消息含清单来源与一基行号。
- `BatchManifestEntry(line_number, sample, reads)`：一条清单项，含一基清单行号、唯一样本名与已解析的 FASTQ 路径。
- `read_batch_manifest(source)`：接受清单文本路径或文本流，忽略空行，返回按清单次序的 `BatchManifestEntry` 元组；文件清单的相对 `reads` 按其目录解析，流清单按当前工作目录解析。
- `batch_call_variants(manifest, reference, *, min_base_quality=20, min_alt_count=2, min_alt_fraction=0.2, homozygous_fraction=0.8, call_indels=False, max_indel_length=50)`：`manifest` 为清单路径/流或 `BatchManifestEntry` 序列，`reference` 为 FASTA 路径/流或 `SequenceRecord` 可迭代对象。按清单顺序对每个样本的 FASTQ 读段以 `call_variants` 语义独立调用（阈值参数与约束相同），样本名取清单 `sample`；各样本结果按 `summarize_variants` 的参考规范化规则合并，不产生中间 VCF。返回按参考次序、`POS`、`REF`、`ALT` 排序的 `VariantSummary` 元组，每项的 `samples` 保持清单次序；无变异时返回空元组。非法阈值在消费任何读段前抛 `ValueError`，空参考或重复参考标识抛 `VariantCallingError`，清单问题抛 `BatchManifestError`，读段格式、碱基与质量问题抛对应现有异常，文件访问失败抛 `OSError`。

序列只允许 IUPAC DNA 符号 `ACGTRYSWKMBDHVN`（小写输入会转大写）；FASTQ 质量字符范围为 ASCII 33–126，质量长度必须等于序列长度。

异常：

- `SequenceFormatError`：空输入、无法判定格式、空标识、空序列、标题前出现内容、FASTQ 结构不完整或长度不等；消息包含来源名与一基行号。
- `VcfFormatError`（位于 `genome_variant.vcf`）：缺少或重复 `#CHROM` 列头、列头不合法、固定列/FORMAT/样本列数与列头不符、非正整数 POS、空或非法 REF/ALT、以及规范化时记录越界；消息包含来源与一基行号。
- `ReferenceMismatchError`（位于 `genome_variant.vcf`）：参考中缺少或重复 CHROM，或 REF 与参考不一致；消息包含 CHROM、POS 以及可判定的期望值与实际值。
- `AnnotationFormatError`（位于 `genome_variant.annotation`）：GFF3 CDS 行列字段数、坐标或 `phase` 非法、缺少 `Parent`、链非 `+`/`-`、片段跨序列或阅读框矛盾、CDS 越过序列末端，或 VCF 已声明/已使用 `GVANN`；消息包含来源与一基行号。
- `ManifestError`（位于 `genome_variant.summary`）：清单结构/字段/样本唯一性错误，或样本 VCF 结构与基因型字段不合规、同一 VCF 规范化后键重复；消息含清单文件与行号，或 VCF 来源与一基行号。
- `BatchManifestError`（位于 `genome_variant.batch`）：批量清单结构/字段/样本唯一性错误，或 `reads` 为 `-`；消息含清单来源与一基行号。
- `SequenceValidationError`：非法碱基；消息指出记录标识、符号与一基位置。
- `ReadQualityError`（`ValueError` 子类，位于 `genome_variant.quality`）：记录没有质量值、质量数量与序列长度不等或质量值超出 0 至 93；消息包含记录标识。
- `VariantCallingError`（`ValueError` 子类，位于 `genome_variant.calling`）：变异调用时参考集合为空、参考标识重复或样本名非法。
- 打开路径失败保留 `OSError`，迭代期间的其他异常不会被吞掉。

## 限制

- 仅支持文本 FASTA/FASTQ 与 IUPAC DNA 符号；VCF 支持记录读取、参考校验、序列型等位基因的最简化及左对齐，以及基于 GFF3 CDS 的单碱基替换注释（标准遗传密码），尚不支持插入/缺失与非标准遗传密码等后续注释流程。
- 变异调用目前仅支持单样本的 `A`、`C`、`G`、`T` 单碱基替换，不调用插入、缺失、多等位位点或歧义碱基；映射计分参数固定为 `map-reads` 的默认值。
