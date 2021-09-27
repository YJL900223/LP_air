**LEfSe**
==============

LEfSe (Linear discriminant analysis Effect Size) determines the features
(organisms, clades, operational taxonomic units, genes, or functions)
most likely to explain differences between classes by coupling standard
tests for statistical significance with additional tests encoding
biological consistency and effect relevance.

LEfSe is available as a [Galaxy module](http://huttenhower.org/galaxy/),
a Conda formula, a Docker image, and included in bioBakery (VM and
cloud). For additional information, please refer to the [LEfSe
paper](http://www.ncbi.nlm.nih.gov/pubmed/21702898).

## Installation

LEfSe can be installed with Conda or run from a Docker image. Please
note, if you are using bioBakery (Vagrant VM or cloud) you do not need
to install LEfSe because the tool and its dependencies are already
installed.

Install with Conda: `$ conda install -c biobakery lefse`

Install with Docker: `$ docker run -it biobakery/lefse bash`

We provide support for LEfSe users. Please join our [bioBakery Support Forum](https://forum.biobakery.org/c/Downstream-analysis-and-statistics/LEfSe) designated specifically for LEfSe users. 

```sh
conda create -n lefse python=2.7
conda activate lefse
conda install -c biobakery lefse
format_input.py tax_lefse.txt input.in -c 1 -o 1000000
run_lefse.py input.in input.res
plot_cladogram.py input.res cladogram.pdf --format pdf
plot_res.py input.res res.pdf --format pdf
```

```sh
# 23、LEfSe、STAMP统计分析和可视化

## LEfSe windows + 网页分析

# 1. 23STAMP目录中准备otutab.txt, design.txt, taxonomy.txt三个文件；
# 2. Rstudio打开taxonomy_summary.Rmd并Knit生成输入文件和可重复计算网页；
# 3. 打开tax_lefse.txt并在线提交 http://www.ehbio.com/ImageGP/index.php/Home/Index/LEFSe.html


## LEfSe Linux分析(选学)

mkdir -p ../23STAMP
cd ../23STAMP
# 上传 tax_0LEfSe.txt 文件

# 格式转换为lefse内部格式
lefse-format_input.py tax_0LEfSe.txt input.in -c 1 -o 1000000
# 运行lefse
run_lefse.py input.in input.res
# 绘制物种树注释差异
lefse-plot_cladogram.py input.res cladogram.pdf --format pdf
# 绘制所有差异features柱状图
lefse-plot_res.py input.res res.pdf --format pdf
# 绘制单个features柱状图(同STAMP中barplot)
head input.res # 查看差异features列表
lefse-plot_features.py -f one --feature_name "Bacteria.Firmicutes.Bacilli.Bacillales.Planococcaceae.Paenisporosarcina" \
   --format pdf input.in input.res Bacilli.pdf 
# 批量绘制所有差异features柱状图，慎用(几百张差异结果柱状图阅读也很困难)
mkdir -p features
lefse-plot_features.py -f diff --archive none --format pdf \
  input.in input.res features/
```

