# QIIME2 2021.2分析流程
	
	# 注：本文以QIIME 2 2021.2版本为例，访问https://docs.qiime2.org/或查询最新版
	
## 0. 软件安装
	# 仅限在Linux/Mac系统中运行
	# Windows用户可选使用服务器、内置Linux子系列Ubuntu(支持右键粘贴)或VirtualBox虚拟机中的Linux运行
	# 详细教程参阅官网https://docs.qiime2.org/2021.2/，或中文教程见宏基因组公众号 https://mp.weixin.qq.com/s/5jQspEvH5_4Xmart22gjMA QIIME2教程部分
	# 安装代码见附录

## 1. 准备工作
	# 设置工作目录，如服务器为~/qiime2，Win子系统C盘中为
	# wd=/home/htt-2022/文档/gmu/dai/nc/
	wd=/home/htt-2022/文档/gmu/dai/case
	# 进入工作目录
	# mkdir -p ${wd}
	cd ${wd}
	
	#linux下运行代码自动保存
  script -a dai_record.txt
  #脚本已启动，文件为 typescript
  #结束后退出代码保存
  # exit
	
	#激活conda后激活qiime2工作环境
	# # 2021.2
	# source /data/root/miniconda3/bin/activate
	# conda activate qiime2-2021.2
	# # 2022.2
	# source /data/root/miniconda3/bin/activate
	# conda activate qiime2-2022.2
	#2022.8（最新版，容易不兼容，报错较多）/2022.2
	source /data/root/miniconda3/bin/activate
  conda activate qiime2-2022.8
  # ##序列质量检查
  
	# 准备样本元信息metadata、原始数据seq/*.fastq和manifest
  # 检查文件大小，确定是否下载完整或正常
  ls -lsh seq 
	# 根据metadata生成manifest文件
	awk 'NR==1{print "sample-id\tforward-absolute-filepath\treverse-absolute-filepath"} \
	  NR>1{print $1"\t$PWD/seq/"$1"_1.fq.gz\t$PWD/seq/"$1"_2.fq.gz"}' \
	  metadata.txt > manifest

	# # 数据导入qiime2，格式为双端33格式
	qiime tools import \
	  --type 'SampleData[PairedEndSequencesWithQuality]' \
	  --input-path manifest \
	  --output-path paired-end-demux.qza \
	  --input-format PairedEndFastqManifestPhred33V2
	# fq文件1G用时7m，fq.gz压缩格式仅34s

  #可视化查看质量
  time qiime demux summarize \
    --i-data paired-end-demux.qza \
    --o-visualization paired-end-demux.qzv

## 2. 生成特征表和代表序列

### 方法1. DADA2

	# # 支持多线程加速，服务器0/96p, 34m；24p, 44m；8p, 77m；1p, 462m
	# # Win10下12线程65m，机时538分
	# trim-left-f 上游引物+barcode切除,trim-left-r 下游引物切除
	# trunc-len-f 序列最短长度；trunc-len-r 序列最大长度
	# threads 运行线程
	time qiime dada2 denoise-paired \
    --i-demultiplexed-seqs paired-end-demux.qza \
    --p-n-threads 30 \
    --p-trim-left-f 19 --p-trim-left-r 20 \
    --p-trunc-len-f 215 --p-trunc-len-r 0 \
    --o-table dada2-table.qza \
    --o-representative-sequences dada2-rep-seqs.qza \
    --o-denoising-stats denoising-stats.qza

	# # 确定使用dada2结果并导入主流程
	cp dada2-table.qza table.qza
	cp dada2-rep-seqs.qza rep-seqs.qza

 # 统计代表序列。
  qiime feature-table tabulate-seqs \
    --i-data dada2-rep-seqs.qza \
    --o-visualization dada2-rep-seqs.qzv

 # 生成统计结果：
  qiime metadata tabulate \
    --m-input-file denoising-stats.qza \
    --o-visualization denoising-stats.qzv

 #按特征表的数据量过滤，只有特征序列总测序量大于10条以上的才保留
 # 详者注：(实验中会有大量低丰度的特征/OTU，
 # 它们会增加计算工作量和影响高丰度结果差异比较时FDR校正Pvalue，
 # 导致假阴性，通常需要选择一定的阈值进行筛选，
 # 常用的有相对丰度千分之五、千分之一、万分之五、万分之一；
 # 也可根据测序总量，reads频率的筛选阈值常用2、5、10等，
 # 大项目样本多数据量大，有时甚至超过100，推荐最小丰度为百万分之一)
 #   # 过滤低丰度，< 10
  qiime feature-table filter-features \
    --i-table table.qza \
    --p-min-frequency 2 \
    --o-filtered-table feature-frequency-filtered-table.qza

   # 统计特征表。
  qiime feature-table summarize \
    --i-table  feature-frequency-filtered-table.qza \
    --o-visualization feature-frequency-filtered-table.qzv \
    --m-sample-metadata-file metadata.txt

## 4. 物种组成分析
  #训练器路径
  # 2021.2 /home/dell/文档/classifier/
  #2022.8 /data/database/silva/silva-138-99-nb-classifier.qza(全长)
	### 物种注释
  # 4m classifier_gg_13_8_99_V5-V7.qza 是我用V5-V7训练的文件，详见附录或官方教程
  time qiime feature-classifier classify-sklearn \
    --i-classifier /data/database/silva/silva-138-99-nb-classifier.qza \
    --i-reads rep-seqs.qza \
    --o-classification taxonomy.qza

	#过滤线粒体、叶绿体
  time qiime taxa filter-table \
    --i-table feature-frequency-filtered-table.qza \
    --i-taxonomy taxonomy.qza \
    --p-exclude mitochondria,chloroplast \
    --o-filtered-table table-no-mitochondria-no-chlorplast.qza

  #过滤不能分类
  time qiime taxa filter-table  \
    --i-table table-no-mitochondria-no-chlorplast.qza  \
    --i-taxonomy taxonomy.qza \
    --p-exclude Unassigned \
    --o-filtered-table table-no-m-c-u.qza

  #过滤真核细胞
  time qiime taxa filter-table  \
    --i-table table-no-mitochondria-no-chlorplast.qza  \
    --i-taxonomy taxonomy.qza \
    --p-exclude d__Eukaryota \
    --o-filtered-table table-no-m-c-u-e.qza

  #保留d__开头的物种，即为只筛选比较确定属于细菌或古菌门的序列，可以有效去除宿主非特异扩增污染和其它未知来源的序列
  qiime taxa filter-table \
    --i-table table-no-m-c-u-e.qza \
    --i-taxonomy taxonomy.qza \
    --p-include d__ \
    --o-filtered-table table-no-m-c-u-e+d.qza

	# 可视化物种注释
	qiime metadata tabulate \
    --m-input-file taxonomy.qza \
	  --o-visualization taxonomy.qzv

	# 统计特征表。
  qiime feature-table summarize \
    --i-table table-no-m-c-u-e+d.qza \
    --o-visualization table-no-m-c-u-e+d.qzv \
    --m-sample-metadata-file metadata.txt

	# 堆叠柱状图展示
	qiime taxa barplot \
	  --i-table table-no-m-c-u-e+d.qza \
	  --i-taxonomy taxonomy.qza \
	  --m-metadata-file metadata.txt \
	  --o-visualization taxa-bar-plots.qzv
  
  mv usearch usearch_v1
  
  # 导出代表序列
  qiime tools export \
    --input-path rep-seqs.qza \
    --output-path usearch/result/raw

  # 导出table
  qiime tools export \
    --input-path table-no-m-c-u-e+d.qza \
    --output-path usearch/result
  # 转换为tsv格式
  biom convert -i usearch/result/feature-table.biom \
    -o usearch/result/table.tsv \
    -to-tsv
  # 删除注释行(可选)
  sed -i '/# Const/d' usearch/result/table.tsv
  cp -i usearch/result/table.tsv usearch/result/otutab.txt
  cp -i usearch/result/raw/dna-sequences.fasta usearch/result/raw/otus.fa
  # 根据筛选后的otutab,生成筛选后的代表序列。
  cut -f 1 usearch/result/otutab.txt | tail -n+2 > usearch/result/otutab.id
  usearch10 -fastx_getseqs usearch/result/raw/otus.fa -labels usearch/result/otutab.id -fastaout usearch/result/otus.fa
  wc -l usearch/result/otutab.txt
  less usearch/result/otus.fa|grep '>' -c
  #
# 3. Alpha和beta多样性分析

# 构建进化树用于多样性分析 53s
  qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences rep-seqs.qza \
    --o-alignment aligned-rep-seqs.qza \
    --o-masked-alignment masked-aligned-rep-seqs.qza \
    --o-tree unrooted-tree.qza \
    --o-rooted-tree rooted-tree.qza


### 方法2. 外部导入特征表和代表序列(常用)，详见附录


## 计算核心多样性
	# 13s，采样深度通常选择最小值，来自table.qzv##	feature-frequency-filtered-
	  qiime diversity core-metrics-phylogenetic \
	    --i-phylogeny rooted-tree.qza \
	    --i-table table-no-m-c-u-e+d.qza \
	    --p-sampling-depth 7391 \
	    --m-metadata-file metadata.txt \
	    --output-dir core-metrics-results


### Alpha多样性组间显著性分析和可视化
	# 7s, 可选的alpha指数有 faith_pd、shannon、observed_features、evenness
	  index=shannon
	  qiime diversity alpha-group-significance \
	    --i-alpha-diversity core-metrics-results/${index}_vector.qza \
	    --m-metadata-file metadata.txt \
	    --o-visualization core-metrics-results/${index}-group-significance.qzv
  

### Beta多样性组间显著性分析和可视化
	# 可选的beta指数有 unweighted_unifrac、bray_curtis、weighted_unifrac和jaccard
	# 7s, 指定分组是减少计算量，置换检验较耗时

	  distance=weighted_unifrac
	  column=Group
	  qiime diversity beta-group-significance \
	    --i-distance-matrix core-metrics-results/${distance}_distance_matrix.qza \
	    --m-metadata-file metadata.txt \
	    --m-metadata-column ${column} \
	    --o-visualization core-metrics-results/${distance}-${column}-significance.qzv \
	    --p-pairwise

    distance=unweighted_unifrac
	  column=Group
	  qiime diversity beta-group-significance \
	    --i-distance-matrix core-metrics-results/${distance}_distance_matrix.qza \
	    --m-metadata-file metadata.txt \
	    --m-metadata-column ${column} \
	    --o-visualization core-metrics-results/${distance}-${column}-significance.qzv \
	    --p-pairwise

    distance=bray_curtis
	  column=Group
	  qiime diversity beta-group-significance \
	    --i-distance-matrix core-metrics-results/${distance}_distance_matrix.qza \
	    --m-metadata-file metadata.txt \
	    --m-metadata-column ${column} \
	    --o-visualization core-metrics-results/${distance}-${column}-significance.qzv \
	    --p-pairwise

    #导出table
    qiime tools export \
      --input-path core-metrics-results/weighted_unifrac_distance_matrix.qza \
      --output-path usearch/result/beta
  
    #导出table
    qiime tools export \
      --input-path core-metrics-results/weighted_unifrac_distance_matrix.qza \
      --output-path usearch/result/beta
  
#usearch&R alpha beta
 
  # 设置流程EasyAmplicon(ea)和工作目录(work directory, wd)，添加环境变量并进入wd
    # **每次打开linux必须运行下面2行**
    ##软件安装后需要配置权限 chmod 775 file 777读写执行；775读执行
    ##环境配置
    ea=/data/software
    PATH=$PATH:${ea}
    ea=/home/htt-2022/文档/gmu/EasyAmplicon/
    PATH=$PATH:${ea}
    
 #可选统计方法：OTU表简单统计 Summary OTUs table
    usearch10 -otutab_stats usearch/result/otutab.txt \
      -output usearch/result/otutab.stat
    cat usearch/result/otutab.stat
    #注意最小值、分位数，或查看result/raw/otutab_nonBac.stat中样本详细数据量，用于重采样

### 等量抽样标准化 normlize by subsample

    #使用vegan包进行等量重抽样，输入reads count格式Feature表result/otutab.txt
    #可指定输入文件、抽样量和随机数，输出抽平表result/otutab_rare.txt和多样性alpha/vegan.txt
    cp -i metadata.txt usearch/metadata.txt
    cd usearch
    mkdir -p temp/
    mkdir -p result/alpha/
    Rscript ${ea}/script/otutab_rare.R --input result/otutab.txt \
      --depth 7505 --seed 123 \
      --normalize result/otutab_rare.txt \
      --output result/alpha/vegan.txt
    usearch10 -otutab_stats result/otutab_rare.txt \
      -output result/otutab_rare.stat
    cat result/otutab_rare.stat
    
    cut -f 1-2 metadata.txt > temp/group.txt

## Alpha多样性 Alpha diversity

### 计算多样性指数 Calculate alpha diversity index
    #Calculate all alpha diversity index(Chao1有错误勿用)
    #details in http://www.drive5.com/usearch/manual/alpha_metrics.html
    usearch10 -alpha_div result/otutab_rare.txt \
      -output result/alpha/alpha_rare.txt
    usearch10 -alpha_div result/otutab.txt \
      -output result/alpha/alpha.txt
## Beta多样性 Beta diversity
    
    #结果有多个文件，需要目录
    mkdir -p result/beta/
    ##qiime2 生成的weighted_unifrac改名，由于缺少weighted_unifrac，写入txt中
    cp -i result/beta/distance-matrix.tsv result/beta/wunifrac.txt
    
    # #基于OTU构建进化树 Make OTU tree, 30s
    usearch10 -cluster_agg result/otus.fa -treeout result/otus.tree
    #生成5种距离矩阵：bray_curtis, euclidean, jaccard, manhatten, unifrac, 3s
    usearch10 -beta_div result/otutab_rare.txt -tree result/otus.tree \
      -filename_prefix result/beta/ # 1s
    # usearch10 -beta_div result/otutab_rare.txt -filename_prefix result/beta/

    ##R生成α多样性
    R
    ## 退出R
    q()
    y
    cd ..
