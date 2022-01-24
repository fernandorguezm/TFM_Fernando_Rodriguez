# TFM_Fernando_Rodriguez

Scripts utilizados para el análisis de los distintos datos necesarios para la realización de mi Trabajo de Fin de Máster (TFM) 

Las carpetas se basan en:
- Procesamiento de RNA-seq: scripts utilizados para el procesamiento de los datos brutos del RNA-seq desde el formato fastq hasta archivos gtf. Se incluyen la generación de los índices del genoma de referencia (piperna_pairend.sh), el análisis de calidad, el mapeado de las lecturas sobre el genoma de referencia, ensamblado y cuantificación de transcritos (rna_seq_processing.sh). Estos scripts pertecen a Francisco Romero Campero (https://github.com/fran-romero-campero), los cuales he modificado para las necesidades de este trabajo en concreto (procesamientos de muestra pair-end).
- Análisis de RNA-seq: esta carpeta consta de diferentes scripts, desde el análisis con DESeq2 y la extracción de genes expresados diferencialmente hasta el clustering jerárquico y la generación de redes de co-expresión génica.
- Análisis fenotípico: script utilizado para la generación de las figuras y análisis estadístico correspondiente para el análisis fenotípico.
- Análisis filogenético: script utilizado para la generación del dotplot del análisis filogenético.
