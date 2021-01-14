---
layout: tutorial_hands_on
title: Tutorial Básico de RNA-Seq
zenodo_link: "https://zenodo.org/record/4249555"
contributors:
- mblue9
- bphipson
- hdashnow
- ccano
---

# Tutorial Básico de RNA-Seq

Las tecnologías de secuenciación de ARN (RNA-Sequencing o RNA-Seq) tienen como objetivo identificar qué loci se expresan en una población celular en un instante dado, es decir, identificar secuencias de ARN en una población celular y cuantificar su abundancia. O sea, caracterizar el transcriptoma celular. Estas tecnologías permiten cuantificar la expresión de genes, descubrir nuevas secuencias transcritas a partir de ADN, identificar genes con splicing alternativo o detectar expresión específica de alelo, entre otros. Además, estas tecnologías han permitido caracterizar no sólo RNA mensajero (mRNA), sino también otros tipos de RNAs como los RNAs que no codifican proteínas (los llamados RNAs no codificantes o non-coding RNAs, ncRNAs) que incluyen los lncRNAs y los miRNAs, entre otros.
Puedes consultar más detalles sobre la secuenciación de ARN en los siguientes enlaces: 

http://cshprotocols.cshlp.org/content/early/2015/04/11/pdb.top084970.abstract

https://www.nature.com/articles/nrg2484

https://galaxyproject.org/tutorials/rb_rnaseq/

Hay muchos pasos involucrados en el análisis de datos de RNA-Seq. Típicamente, este proceso comienza con el procesamiento de lecturas (reads) de un fichero FASTQ que se alinean contra un genoma de referencia para cuantificar el número de secuencias de ARN asociadas a cada loci genético. Esto produce una matriz de números sobre la que podemos realizar análisis estadísticos y computacionales para identificar genes y pathways diferencialmente expresados entre dos grupos de muestras. 

Este tutorial se divide en varias secciones. Primero, se enseña cómo alinear reads contra un genoma de referencia para cuantificar la abundancia de las distintas secuencias de ARN en la muestra. Después, se aplican análisis estadísticos para identificar genes diferencialmente expresados. 

**Este tutorial se ha adaptado a partir de los tutoriales siguientes para Galaxy y R, simplificando numerosos detalles para permitir un primer acercamiento al análisis de datos de RNA-Seq**: 

[Galaxy RNA-Seq reads to counts](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/rna-seq-reads-to-counts/tutorial.md).

[Galaxy RNA-seq counts to genes](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/rna-seq-counts-to-genes/tutorial.md). 

[COMBINE R RNAseq workshop](http://combine-australia.github.io/RNAseq-R/07-rnaseq-day2.html)

[Taller de análisis de expresión de RNA en R](https://sites.google.com/view/taller-r-analisis-rna-seq/análisis-de-expresión-de-rnas?authuser=2)

Te invito a que explores los dos primeros tutoriales de Galaxy cuando estés más familiarizado con la herramienta y los análisis bioinformáticos para RNA-Seq para conocer más detalles sobre estos análisis. 

## El problema

En este tutorial abordaremos el problema propuesto en un paper de Nature Cell Biology: [EGF-mediated induction of Mcl-1 at the switch to lactation is essential for alveolar cell survival](https://www.ncbi.nlm.nih.gov/pubmed/25730472) en 2015. Este trabajo analiza los transcriptomas de células basales y luminales en las glándulas mamarias de ratonas embarazadas, con descendencia lactante y sin descendencia. Se presentan, por lo tanto, seis grupos de estudio, uno para cada combinación de tipo celular y condición.  Para cada grupo de estudio se toman dos muestras para el análisis. Por tanto, el número total de muestras del estudio es 12, como muestra el siguiente esquema. 

![](https://i.imgur.com/IM6SPmS.png)


Los investigadores de este estudio han puesto los datos a disposición de la comunidad, tanto los datos en bruto (raw data) con las lecturas generadas por el secuenciador, como los datos ya preprocesados (matrices numéricas). Todos los datos estan disponibles en el siguiente enlace a la base de datos Gene Expression Omnibus database (GEO) [GSE60450](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450).



# 1. Preparación de las lecturas

## 1. A. Cargar los datos

Existen numerosos algoritmos para alinear lecturas contra un genoma de referencia: HISAT2, STAR y Subread. Estos algoritmos requieren recursos computacionales elevados, y aunque Galaxy nos provee de estos recursos de forma gratuita ([usegalaxy.eu](https://usegalaxy.eu)), **vamos a continuar este tutorial utilizando una pequeña parte de las lecturas originales para aligerar la carga computacional**. Este subconjunto de lecturas está disponible en el siguiente enlace de [Zenodo](https://zenodo.org/record/4249555). Además, **utilizaremos sólo una muestra por cada grupo de estudio, para reducir así el número de muestras a 6**. 

Para facilitar la tarea, a continuación se indica cómo cargar los  ficheros FASTQ comprimidos del enlace de Zenodo como una colección de ficheros (permite asignar cómodamente un nombre de muestra a cada fichero)
```
MCL1-DK-basallactate	https://zenodo.org/record/4249555/files/SRR1552454.fastq.gz
MCL1-DI-basalpregnant	https://zenodo.org/record/4249555/files/SRR1552452.fastq.gz
MCL1-DG-basalvirgin	https://zenodo.org/record/4249555/files/SRR1552450.fastq.gz
MCL1-LE-luminallactate	https://zenodo.org/record/4249555/files/SRR1552448.fastq.gz
MCL1-LC-luminalpregnant	https://zenodo.org/record/4249555/files/SRR1552446.fastq.gz
MCL1-LA-luminalvirgin	https://zenodo.org/record/4249555/files/SRR1552444.fastq.gz
```


Para cargar estos ficheros en Galaxy, haremos lo siguiente:

> ### Carga de datos
>
> 1. Crea una nueva historia para este tutorial e.j. `RNA-seq`
>
> 2. Carga los ficheros con galaxy: 
>    - Abre el *Upload Manager* de Galaxy
>    - Haz click en la pestaña **Rule-based**
>        - *"Upload data as"*: `Collection(s)`
>        - *"Load tabular data from"*: `Pasted Table`
>    - Pega el contenido de la tabla gris de arriba
>    - Click **Build**
>

>    - En el `rules editor`:
>
>        - **Define las columnas con el identificador y URL**. `Add / Modify Column Definitions`
>            - Click `Add Definition` 
>                - *"List Identifier(s)"*: `A`
>            - Click `Add Definition` 
>                - *"URL"*: `B`
>            - Click `Apply`
>
>        - **Nombra la colección**. *"Name"*: `fastqs` 
>        - Click `Upload`
>
>
>        Debes tener una colección (lista) llamada `fastqs` en tu historia que contiene 6 FASTQs.
>
>

## 1. B. Control de Calidad de las lecturas

El control de calidad es un paso esencial en el análisis de datos de secuenciación. En [este tutorial de Galaxy](https://bit.ly/3q7XJ6i) vimos como utilizar [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) para realizar este control de calidad.  

> ### **FastQC**
>
> 1. **FastQC**:
>    - *"Short read data from your current history"*: Seleccionar el tercer icono: `dataset collection` y `fastqs` (Nombre que dimos a la colección de datos de entrada). Click en `Execute`. 
>    
> 2. Revisa las `Webpage` de **FastQC** para las distintas muestras de nuestro conjunto de datos haciendo click en el icono del ojo de cada muestra. 

Como la colección de datos tiene 6 ficheros FASTQ, obtenemos 6 informes de FASTQC. Existe una herramienta llamada MultiQC que agrega y resume la información de FASTQC sobre un conjunto de ficheros FASTQ. 

> ### Utilizar **MultiQC** para agregar informes de FASTQ
>
> 1. **MultiQC**
>      - *"Results"*
>        - *"Which tool was used generate logs?"*: `FastQC`
>        - *"FastQC output"*
>           - *"Type of FastQC output?"*: `Raw data`
>           - *"FastQC output"*: `RawData`



El tutorial de Galaxy [Quality Control tutorial](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/rna-seq-reads-to-counts/tutorial.html) tiene más detalle sobre FASTQC.

> ### ¿Qué opinas sobre la calidad de las secuencias?
>
> En general, las muestras tienen buena pinta y la calidad general es alta. Se han encontrado algunos adaptadores de Illumina y secuencias duplicadas que debieran ser eliminadas, aunque este paso excede los objetivos de este tutorial. Para más detalle, consultar el tutorial de Galaxy: https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/rna-seq-reads-to-counts/tutorial.html). 

# 2. Mapeo contra genoma de referencia

El siguiente paso es alinear las lecturas de las 12 muestras contra el genoma de referencia. Vamos a utilizar la versión del genoma de referencia del ratón `mm10`. Utilizaremos el software [**HISAT2**](https://ccb.jhu.edu/software/hisat2/index.shtml) para el alineamiento. El tutorial de Galaxy [RNA-seq ref-based tutorial](/topics/transcriptomics/tutorials/ref-based/tutorial.html#mapping) contiene más información sobre el proceso de alineamiento / mapeo para RNASeq. Las muestras que estamos estudiando contienen lecturas de tipo single-end y sin especificar hebra (unstranded), así que necesitaremos especificar estas características en los parámetros de HISAT2. Como solo estamos ejecutando el alineamiento sobre un conjunto de 1000 lecturas de cada fichero, este proceso será rápido (si utilizas todo el dataset original, este proceso tomará más tiempo).


> ### **HISAT2**
>
> 1. Ejecuta **HISAT2** con argumentos:
>    -  *"Source for the reference genome"*: `Use a built-in genome`
>        -  *"Select a reference genome"*: `mm10`
>    -  *"Is this a single or paired library?"*: `Single-end`
>        -  *"FASTA/Q file"*: 3er icono `Dataset Collection`: `fastqs` >    - *"Specify strand information"*: `unstranded`
>    - *"Summary Options"*:
>        -  *"Output alignment summary in a more machine-friendly style."*: `Yes`
>        -  *"Print alignment summary to a file."*: `Yes`
> 2. **MultiQC**  with the following parameters to aggregate the HISAT2 summary files
>    - In *"Results"*
>        -  *"Which tool was used generate logs?"*: `HISAT2`
>        -  *"Output of HISAT2"*: `Mapping summary` (output of **HISAT2** )
> 3. Add a tag `#hisat` to the `Webpage` output from MultiQC and inspect the webpage


Es importante comprobar si el porcentaje de lecturas mapeadas al genoma de referencia es alto. En este caso, se obtiene un porcentaje de lecturas mapeadas con éxito superior al 90% y la mayoría de lecturas se han mapeado a una única posición (poca ambigüedad, lo deseable). 

**HISAT2** genera un fichero con el alineamiento en formato BAM.


# 3. Crear la matriz de conteos

Ahora que el alinemiento indica de qué posición del genoma proviene cada lectura, podemos cuantificar la abundancia de lecturas transcritas de cada gen. Para ello, utilizamos un software que se llama `featureCounts`. Las lecturas que mapean a posiciones de exones de los genes se suman para contar cuan abundante es la transcripción de cada gen. La salida es una matriz de conteos para cada gen de Entrez (los identificadores de genes de Entrez son códigos como  `100008567`). 


## 3.A. Cuantificar lecturas asociadas a genes

> **featureCounts**
>
> 1. **featureCounts** :
>    - *"Alignment file"*: `HISAT2 on aligned reads (BAM)`
>    - *"Gene annotation file"*: `featureCounts built-in`
>        - *"Select built-in genome"*: `mm10`
>
> 2. **MultiQC** se utiliza para agregar resultados, no solo de FASTQC sino también de otro sofware como featureCounts. Para obtener útiles tablas resumen a partir de los resultados de featureCounts, ejecuta **MultiQC**:
>    - *"Which tool was used generate logs?"*: `featureCounts`
>        - *"Output of FeatureCounts"*: `featureCounts summary` 

> ### Pregunta
>
> ¿Qué % de lecturas se asignan a exones?
>
> > ###
> > ~60-70% de las lecturas mapean a exones. 


Los conteos de RNAs para cada muestra están en un fichero distinto. Abre uno de estos ficheros (icono del ojo) para echarle un vistazo. La primera columna contiene el Entrez ID del gen y la segunda el número de moléculas de RNA que han mapeado contra alguna posición de ese gen. 

## 3. B. Crear la matriz de conteos

En lugar de disponer de los conteos de cada muestra en distintos ficheros, vamos a agregar los resultados de los conteos en una única matriz, con un gen por cada fila y una columna por cada muestra. 

> ### **Column Join on Collection**
>
> **Column Join on Collection** :
>    - *"Tabular files"*: `Counts` (salida de **featureCounts** )
>    - *"Identifier column"*: `1`
>    - *"Number of header lines in each input file"*: `1`
>    - *"Add column name to header"*: `No`
>

Si echas un vistazo a esta matriz de conteos (utilizando el icono del ojo) comprobarás que el número de moléculas de RNA asociadas a cada gen en las muestras es 0 para casi todos los genes y muestras. Esto se debe a que hemos utilizado solo 1000 lecturas de cada muestra. Si exploras la matriz encontrarás que algunos genes si tienen cierta presencia (1, 2, 3 copias) en algunas muestras. En experimentos reales, utilizando todas las lecturas, estas matrices muestran valores mucho mayores. 


# 4. Análisis de la matriz de conteos

La obtención de la matriz de conteos ha requerido los pasos 1, 2 y 3 de este tutorial. Llegados a este punto, y por si alguno de estos pasos ha supuesto un problema, se proporcionan las matrices de conteos ya generadas para que os las podáis descargar y continuar siguiendo el resto del tutorial desde aquí. Estas matrices de conteos contienen los datos de todas las lecturas para todas las muestras del estudio original, así que el resultado del análisis será más preciso si utilizáis estas matrices en lugar de las que habéis generado en la parte anterior del tutorial. 


## Descarga estos datos para empezar por aquí. 

Proporcionamos tres datos para empezar por este punto el análisis de las matrices de conteos: 

 * **Matriz de conteos** (un gen por fila y una muestra por columna)
 * **Información de las muestras** (id muestra, grupo)



> ### Datos
>
> 1. Crea una nueva historia para este ejercicio `RNA-seq tutorial 2`
>
> 2. Importa los ficheros:

>     - Cuadro Carga de datos - **Paste/Fetch**:
>
>     ```
>     https://zenodo.org/record/4273218/files/countdata.tsv
>     https://zenodo.org/record/4273218/files/factordata_fixed.tsv
>     ```
>
> 2. Cambia el nombre a los ficheros para que se llamen simplemente `countdata` (matriz de conteos) y`factordata` (información de las muestras) utilizando el icono del lápiz. 
> 3. Comprueba en la pestaña `Datatypes` que el tipo de fichero es `tabular`. Si no lo fuera, cambia el tipo a `tabular`.

Vamos a explorar los datos. El fichero `countdata` contiene un gen por fila (el ID es el Entrez ID) y una columna por muestra, indicando con números el número de lecturas que han mapeado contra ese gen en cada muestra del estudio. 

El fichero `factordata` contiene el tipo de cada una de las muestras. Observa que el id de la muestra coincide con el que aparecía en el fichero `countdata`.


# 4.1 Expresión diferencial 

## Filtrado de genes planos

Un paso habitual para simplificar el análisis y mejorar la potencia estadística de los test que potencialmente se apliquen sobre los datos es eliminar los genes que apenas se expresan en ninguna muestra (los genes que tienen bajo número de copias en las muestras).   Para ello, se suele emplear una escala de medida de número de copias de un gen (loci) por millón de lecturas (counts-per-million , CPM), de manera que sólo se seleccionan los loci que presentan un número de lecturas mayor que el umbral en, al menos, un cierto número de muestras. En nuestro caso, seleccionamos los genes con, al menos, 0.5 CPM en, al menos, dos muestras del estudio. Por las características de este experimento, un CPM de 0.5  significa que seleccionaremos aquellos genes para los que hayamos detectado, al menos, 10-15 lecturas (copias) en al menos dos de las doce muestras del estudio. Consideramos que los genes que no lleguen a ese mínimo apenas se expresan y serán descartados del análisis. 
Para realizar este filtrado necesitamos calcular el número de CPMs para cada gen. Para ello, tendríamos que tener en cuenta, además,  que cada experimento de secuenciación produce un número total de lecturas distinto y distinto volumen de lecturas según la región genómica de que se trate (esto se denomina profundidad de lectura o sequencing depth). 

## Definir los grupos a contrastar

Para los estudios de expresión diferencial debemos definir dos grupos de muestras entre los que vamos a realizar los tests estadísticos. Por ejemplo, si estamos interesados en conocer genes diferencialmente expresados entre los grupos pregnant y lactating en las células basales, especificaremos como *Contrast of Interest* los grupos: `basalpregnant-basallactate`. Los nombres de los grupos deben coincidir escrupulosamente con los nombres que aparecían en el fichero `factordata`. Se pueden evaluar distintas comparaciones entre grupos utilizando el botón `Insert Contrast`. Primero, vamos a observar las diferencias entre `basalpregnant-basallactate`.

> ### Expresión diferencial con limma-voom
>
> 1. **limma** {% icon tool %}:
>      - *"Differential Expression Method"*: `limma-voom`
>      -  *"Count Files or Matrix?*": `Single Count Matrix`
>          -  *"Count Matrix"*: Select `countdata`
>      -  *"Input factor information from file?"*: `Yes`
>          -  *"Factor File"*: Select `factordata`
>      -  *"Use Gene Annotations?"*: `Yes`
>          -  *"Factor File"*: Select `annodata`
>      -  *"Contrast of Interest"*: `basalpregnant-basallactate`
>      -  *"Filter lowly expressed genes?"*: `Yes`
>          -  *"Filter on CPM or Count values?"*: `CPM`
>          -  *"Minimum CPM"*: `0.5`
>          -  *"Minimum Samples"*: `2`
>      - **Output Options**
>          - {% icon param-check %} *"Additional Plots"* tick:
>              - `Glimma Interactive Plots`
>              - `Density Plots (if filtering)`
>              - `CpmsVsCounts Plots (if filtering on cpms)`
>              - `Box Plots (if normalising)`
>              - `MDS Extra (Dims 2vs3 and 3vs4)`
>              - `MD Plots for individual samples`
>              - `Heatmaps (top DE genes)`
>              - `Stripcharts (top DE genes)`


> 2. Una vez ejecutado, inspecciona el `Report` haciendo click en el icono del ojo. 



# QC

Antes de comprobar los genes diferencialmente expresados, podemos echar un vistazo al informe generado por `limma` (icono del ojo) que contiene distintas visualizaciones de los datos. Revisamos algunas de ellas a continuación. 

## Multidimensional scaling plot

Una de las visualizaciones más interesantes para analizar datos de RNA-Seq son los diagramas multidimensionales (MDS plots). Este tipo de análisis está íntimamente ligado al Análisis de Componentes Principales (PCA), y gracias a los diagramas MDS podemos visualizar los datos respecto a las variables que introducen una mayor fuente de variación. Este tipo de gráficos también permite detectar outliers en los datos.
Si el problema es sencillo, con un número de muestras reducido, puede observarse a simple vista si la mayor fuente de variación en los datos se debe al grupo/tratamiento de la muestra. Por ejemplo, en el MDS plot generado se observa que bastan las dos primeras componentes principales de los datos de expresión para discriminar cada uno de los 6 tipos de muestras.


> ### Preguntas
>
> Comprueba los nombres de las muestras en el MDS plot. ¿Cuál es la principal fuente de variabilidad en los datos? (es decir, ¿qué significa la dimensión 1 del plot?)
> ¿Cuál es la segunda fuente de variabilidad de los datos?
>
> ¿Cuántos genes se han eliminado del análisis por ser genes planos?
>


## Density plots y Box plots

Los plots de densidad o Density plots muestran la distribución de los counts de las distintas muestras antes y después del filtrado de genes planos. Los diagramas de cajas o box plots también permiten visualizar estas distribuciones. En estos plots podemos apreciar que las distribuciones de counts no son idénticas entre las muestras, pero tampoco resulta demasiado diferentes. Si una muestra estuviera significativamente desplazada hacia arriba o abajo de la línea horizontal azul (que marca la media de las distribuciones de las muestras), sería necesario analizarla en mayor detalle o retirarla del análisis. 


## Heatmaps

Haz click en el enlace `Heatmap_basalpregnant-basallactate.pdf`, que muestra la expresión de los 10 genes más diferencialmente expresados (según p-valor ajustado) entre los dos conjuntos de muestras bajo estudio: `basalpregnant-basallactate`. 

> ### Tarea
>
> Comprueba si en el heatmap se observa una expresión diferenciada de los genes entre los grupos `basalpregnant-basallactate`. 


# Conclusion

Este tutorial pretende aportar un primer contacto con el análisis de datos de RNA-Seq. Te recomiendo que ahora eches un vistazo a los tutoriales de [Galaxy Training](https://training.galaxyproject.org/training-material/), en particular, a los [tutoriales completos de RNA-Seq](https://training.galaxyproject.org/training-material/topics/transcriptomics/) a partir de los cuales se ha obtenido, a modo de resumen, este tutorial. Otra forma natural de continuar el trabajo es explorando el siguiente paso del pipeline de análisis: [la anotación funcional de los genes diferencialmente expresados para facilitar su interpretación biológica](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/rna-seq-genes-to-pathways/tutorial.html).   

