# EpiGenToolKit
Is a small library created to deal with data from EPIC / 450K microarrays. The tool allows to:

a) Simply visualize methylation levels of specific CpG or genomic region. 

b) Perform enrichment analysis of a selected subset of CpG against the whole array. In this type of analysis expected frequency [%] (based on mynorm) of genomic regions is compared to observed (based on provided cpgs set), results are comapred using chi-square test.

## How to start?

a) using env

        python -m venv env 
        source env/bin/activate # Windows: env\Scripts\activate
        pip install epitoolkit

b) using poetry

        poetry new .
        poetry add epitoolkit

c) or just clone the repository:

        git clone https://github.com/ClinicalEpigeneticsLaboratory/EpiGenToolKit.git
        cd EpiGenToolKit && poetry install

## How to use?


### Visualization

To visualize single CpG site or specific genomic region initialize **Visualise** object:
    
        from epitoolkit.tools import Visualize
        viz = Visualize(manifest=<path_to_array_manifest>,
                        mynorm=<path_to_mynorm_file>,
                        poi=<path_to_poi_file>,
                        skiprows=7) # many manifest contains headers, set skiprows argument to ignore few first rows.
    
all files must have *.csv extension, mynorm must contain samples [**columns**] and cpgs [**rows**], the proper EPIC manifest may be downloaded from the [link](https://emea.support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html), poi file must contain sample names [**rows**] (the same as in mynorm file) and POI (phenotype of interest) column containing names of phenotype e.g. Control and Case.

To visualize single CpG:

    viz.plot_CpG("cg07881041", # cpg ID 
        static=False, # plot type static / interactive
        height=400, # plot size
        width=700, # plot size
        title="", # plot title
        legend_title="", # legend title
        font_size=22, # font size
        show_legend=True, # False to hide legedn
        x_axis_label="CpG", # x axsis label
        y_axis_label="beta-values") # y axis label



![CpGPlot](https://github.com/ClinicalEpigeneticsLaboratory/EpiGenToolKit/blob/main/Plots/Plot1.png?raw=true)


To visualize specific genomic region:

        vis.plot_Range(chr=17, start=5999, end=6770) # All arguments describing plot details such as height or width (described above) are also available in this method.

![CpGPlot](https://github.com/ClinicalEpigeneticsLaboratory/EpiGenToolKit/blob/main/Plots/Plot2.png?raw=true)


To visualize specific CpGs in genomic order, instead of whole region pass collection of CpGs:
        
        viz.plot_Range(cpgs=["cg04594855", "cg09002677"]) 


![CpGPlot](https://github.com/ClinicalEpigeneticsLaboratory/EpiGenToolKit/blob/main/Plots/Plot3.png?raw=true)


To save plots use *export* argument for instance:

        viz.plot_Range(chr=17, start=5999, end=6770, export="plot.html") # if static = False only html format is supported if static = True, use png extension.


### Enrichment analysis

To perform enrichment analysis against any type of genomic region specified in the manifest file, the user needs to initialize **EnrichemntAnalysis** object.

        from src.epitoolkit.tools import EnrichemntAnalysis
        ea = EnrichemntAnalysis(manifest=<path_to_array_manifest>,
                mynorm=<path_to_mynorm_file>)

or if Analysis object already exists:

        ea = EnrichemntAnalysis(manifest=viz.manifest,
                mynorm=viz.mynorm)

To start analysis:

        ea.enrichemntAnalysis(categories_to_analyse=["UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island"],  # list of categories to analyse
                        cpgs=cpgs) # list of cpgs to analyse against backgorund

![examplePlot](https://github.com/ClinicalEpigeneticsLaboratory/EpiGenToolKit/blob/main/Plots/Plot4.png?raw=true)
