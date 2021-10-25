# EpiGenToolKit
Is a small library created to deal with data from `EPIC / 450K` microarrays. The tool allows to:

a) Simply visualize methylation levels of specific CpG or genomic region.

b) Perform enrichment analysis of a selected subset of CpG against the whole array. In this type of analysis expected frequency [%] (based on mynorm) of genomic regions is compared to observed (based on provided cpgs set), results are comapred using chi-square test.

# How to start?

a) using env


```
python -m venv env
source env/bin/activate # Windows: env\Scripts\activate
pip install epitoolkit
```

b) using poetry

```
poetry new .
poetry add epitoolkit
```

c) or just clone the repository:


```
git clone https://github.com/ClinicalEpigeneticsLaboratory/EpiGenToolKit.git
cd EpiGenToolKit && poetry install
```

# How to use?

## Visualization

To visualize single **CpG** site or specific genomic region initialize **Visualise** object:

```
from epitoolkit.tools import Visualize

viz = Visualize(manifest=<path_to_array_manifest>, # path to manifest file
                mynorm=<path_to_mynorm_file>, # path to mynorm file
                poi=<path_to_poi_file>, # path to poi file
                poi_col=<column_name> # name of column containing sample phenotype
                skiprows=0) # many manifest contains headers, set skiprows argument to ignore them.
```
all files must have *.csv extension, mynorm must contain sample names as `columns` and cpgs as `rows`, the proper
EPIC manifest may be downloaded from [here](https://emea.support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html),
poi file must contain sample names `rows` (only samples overlapped between poi and mynorm will be used)
and POI (phenotype of interest) column containing names of phenotype e.g. Control and Case.

To visualize single CpG:
```
viz.plot_CpG("cg07881041", # cpg ID
    static=False, # plot type static / interactive [default]
    height=400, # plot size [default]
    width=700, # plot size [default]
    title="", # plot title [default]
    legend_title="", # legend title [default]
    font_size=22, # font size [default]
    show_legend=True, # False to hide legedn [default]
    x_axis_label="CpG", # x axsis label [default]
    category_order=["Cohort 1", "Cohort 2], # box order [default]
    y_axis_label="beta-values") # y axis label [default]
```
> NOTE: most of those arguments are default! So you don't need to specify most of them!

![CpGPlot](https://github.com/ClinicalEpigeneticsLaboratory/EpiGenToolKit/blob/main/Plots/Plot1.png?raw=true)


To visualize specific genomic region:
```
vis.plot_Range(chr=17, start=5999, end=7000)
```

> NOTE: please note that all arguments available in `viz.plot_CpG` are also in `plot_Range`

![CpGPlot](https://github.com/ClinicalEpigeneticsLaboratory/EpiGenToolKit/blob/main/Plots/Plot2.png?raw=true)


To visualize specific CpGs in genomic order, instead of whole region, just pass collection of CpGs:
```
viz.plot_Range(cpgs=["cg04594855", "cg19812938", "cg05451842"]
```

![CpGPlot](https://github.com/ClinicalEpigeneticsLaboratory/EpiGenToolKit/blob/main/Plots/Plot3.png?raw=true)


To save plots use *export* argument, for instance:
```
viz.plot_Range(chr=17, start=5999, end=6770, export="plot.html") # if static = False only html format is supported if static = True, use png extension.
```

### Enrichment analysis

To perform enrichment analysis against any type of genomic region specified in the manifest file, the user needs to initialize **EnrichemntAnalysis** object.
```
from src.epitoolkit.tools import EnrichmentAnalysis

ea = EnrichmentAnalysis(manifest=<path_to_array_manifest>,
        mynorm=<path_to_mynorm_file>)
```

or if `Visualize` object already exists use `load` method (this approach makes you not have to load the data again):
```
ea = EnrichmentAnalysis.load(<Visualize_object_name>)
```
To start analysis:

```
ea.enrichmentAnalysis(categories_to_analyse=["UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island"],  # list of categories to analyse
                cpgs=cpgs) # list of cpgs to analyse against background
```

![examplePlot](https://github.com/ClinicalEpigeneticsLaboratory/EpiGenToolKit/blob/main/Plots/Plot4.png?raw=true)
