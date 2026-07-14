```bash
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir <genomeDir> \
--genomeFastaFiles <genome.fa> --sjdbGTFfile <genomeannotation.gtf> --sjdbOverhang <readlength - 1> 
```

snRNAseq 10xv3
```bash
STAR --genomeDir <genomDir> --readFilesIn <R2.fastz.gz> <R1.fastz.gz> --soloType CB_UMI_Simple --soloCBwhitelist 3M-february-2018.txt --soloUMIlen 12 --clipAdapterType CellRanger4 \
--outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
--soloCellFilter EmptyDrops_CR 10000 0.99 10 45000 90000 500 0.01 20000 0.01 10000 --soloFeatures GeneFull --soloMultiMappers EM \
--outSAMtype None --runThreadN 12 --outFilterMultimapNmax 20 --readFilesCommand zcat --outFileNamePrefix <outDir>
```
