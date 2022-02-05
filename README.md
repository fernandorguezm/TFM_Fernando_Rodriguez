# TFM_Fernando_Rodriguez

Scripts utilizados para el análisis de los distintos datos necesarios para la realización de mi Trabajo de Fin de Máster (TFM) 

Las carpetas se basan en:
- Procesamiento de RNA-seq: scripts utilizados para el procesamiento de los datos brutos del RNA-seq desde el formato fastq hasta archivos gtf. Se incluyen la generación de los índices del genoma de referencia (piperna_pairend.sh), el análisis de calidad, el mapeado de las lecturas sobre el genoma de referencia, ensamblado y cuantificación de transcritos (rna_seq_processing.sh). Estos scripts pertecen a Francisco Romero Campero (https://github.com/fran-romero-campero), los cuales he modificado para las necesidades de este trabajo en concreto (procesamientos de muestra pair-end).
- Análisis del RNA-seq: esta carpeta consta de diferentes scripts:
  1. DESeq2_analysis: análisis de la expresión génica diferencial con el paquete DESeq2. Análisis de PCA.
  2. DEGs: establecimiento de los criterios para la selección de los DEGs. Generación de tablas con los DEGs de cada una de las condiciones y los mutantes analizados. Generación de Volcano plots para los mutantes.
  3. Hierarchical_clustering: análisis de agrupamiento jerárquico para los DEGs de cada especie. Generación y representación de los 30 clústeres y selección y representación de los 10 clústeres con más DEGs. Anotación de los diferentes clústeres para la posterior representación en la red de co-expresión génica. 
  4. Co-expression_network: generación de la red de co-expresión génica. Análisis del ajuste de dicha red a una red libre de escala.
  5. Orthogroups: análisis de ortogrupos para obtener aquellos ortogrupos comunes entre los DEGs de las diferentes especies regulados por luz, temperatura o ambos. 
  6. Venn_diagrams: generación de diagramas de Venn de los genes regulados por luz y temperatura de todas las especies.
  7. Heatmap_mutants: análisis de la expresión de los DEGs de los mutantes frente a la cepa silvestre. Generación de los diagramas de Venn de los mutantes de las plantas terrestres.
  8. Hierarchical_clustering_mutants: análisis de los clústeres seleccionados en el script 3.Hierarchical_clustering.R en los mutantes.
  9. Orthogroups_mutants: análisis de los ortogrupos comunes entre los DEGs de los mutantes de las plantas terrestres.
- Análisis fenotípico: script utilizado para la generación de las figuras y análisis estadístico correspondiente para el análisis fenotípico.
- Análisis filogenético: script utilizado para la generación del dotplot del análisis filogenético.
