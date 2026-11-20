> 更新：项目已获用户同意改为探索性组织组比较。下文关于逐样本数据的缺口仍成立，但不再阻止当前描述性项目；参见 docs/project-status.zh.md。

# 这次上传文件的检查结果

已直接用 gzip 读取原始上传文件，无需手动解压。上传文件与项目已有 GEO 下载文件的 SHA-256 完全一致。可重复运行：`python3 scripts/prepare_uploaded_retina.py`。

## 第一步：认识表格

原表有 57,815 行记录、25 列。前 9 列是注释：tracking_id、class_code、nearest_ref_id、gene_id、gene_short_name、tss_id、locus、length、coverage。其余 16 列属于四个组织组，每组含 FPKM、conf_lo、conf_hi、status 四个字段。它们不是 16 个样本。

四个组是 Macular_RPE_choroid_sclera、Macular_Retina、Peripheral_RPE_choroid_sclera、Peripheral_Retina。FPKM 是表达量；conf_lo/conf_hi 是原文件给出的估计区间上下界，不能当成不同 donor；status 是原始估计状态。

## 第二步：确认样本身份

通过项目保存的 GEO Series Matrix 样本标题和 GSM 编号，重新生成 data/metadata.csv。这里每行是一个样本，不是一个基因。sample_id 是 GEO 编号；eye_id 是配对眼标识；region 是 macular 或 peripheral；geo_title 保留来源标题。八个 eye_id 分别为 11-1516、11-1556、11-1614、11-1624、11-1648、11-1833、11-1875、11-2043，每个恰好有一个黄斑和一个周边样本。共 16 个 retinal samples，排除了另 16 个支持组织样本。

样本身份来自 GEO 元数据，不是从上传表达表推断出的独立表达列。保留 eye_id 这个原始含义，不从它额外推断人口学信息或双眼关系。

## 第三步：提取视网膜表达

data/retinal_group_expression.csv 保留全部 57,815 行的 tracking_id、gene_id、gene_short_name、locus，以及两个 Retina 组的 FPKM、区间、状态。RPE/choroid/sclera 的表达列已移除。没有改写表达量、填补数据、合并基因名或做统计筛选。tracking_id 全部唯一，但不据此假定 gene_short_name 也唯一。

两个 Retina FPKM 列没有缺失值，全部是有限的非负数；两个状态列均为 OK。OK 只说明原表的估计状态，不代表完成了样本层面的质量控制。

## 第四步：理解目前能做什么

这份 retinal-only 表可以用于查看组织组表达值、识别候选基因；它不是 proposal 要求的 16 样本输入。配对分析需要对每个基因算出同一 eye_id 的黄斑减周边差值，共 8 个差值。目前每个基因只有两个组织组估计，无法算出这 8 个差值。不能复制组值来凑样本，也不能用区间上下界替代样本。

项目继续保持 processed FPKM、8 对样本、预定义 rod/cone gene sets 的范围。下一项真正缺少的材料是带 GSM 对应关系的逐样本 gene-level FPKM，优先向原研究作者或导师获取。未擅自转为原始测序重分析，也未运行配对检验、样本 PCA 或报告研究结论。

results/upload_inspection.json 保存完整列名、文件指纹与检查数字；docs/data-access.md 说明所缺输入的格式。
