#' Fusioncatcher data
#'
#' Documentation for the Fusioncatcher example data.
#'
#' @name raw_fusioncatcher
#'
#' @section fusioncatcher_833ke_final-list-candidate-fusion-genes.txt:
#'
#' This file has the results from a run of Fusioncatcher-0.99.3e on the 833ke
#' cell line. The program was ran with the standard configuration file and with
#' the parameters "-p 8 -z --keep-preliminary".
NULL

#' InFusion data
#'
#' Documentation for the InFusion example data.
#'
#' @name raw_infusion
#'
#' @section infusion_fusions.txt:
#'
#' This is example data from the InFusion getting started page located at
#' https://bitbucket.org/kokonech/infusion/wiki/Getting%20Started .
NULL

#' PRADA data
#'
#' Documentation for the PRADA example data.
#'
#' @name raw_prada
#'
#' @section PRADA.acc.fusion.fq.TAF.tsv:
#'
#' This is example data thankfully provided by PRADA authors Siyuan Zheng and
#' Roeland Verhaak.
NULL

#' JAFFA data
#'
#' Documentation for the JAFFA example data.
#'
#' @name raw_jaffa
#'
#' @section jaffa_results.csv:
#'
#' This is example data from the described JAFFA example run documented at
#' https://github.com/Oshlack/JAFFA/wiki/Example .
NULL

#' EricScript data
#'
#' Documentation for the EricScript example data.
#'
#' @name raw_ericscript
#'
#' @section ericscript_SRR1657556.results.total.tsv:
#'
#' This is example data thankfully provided by EricScript author Matteo Benelli.
NULL

#' FusionMap data
#'
#' Documentation for the FusionMap example data.
#'
#' @name raw_fusionmap
#'
#' @section FusionMap_01_TestDataset_InputFastq.FusionReport.txt:
#'
#' This is example data provided with the FusionMap version released 2015-03-31.
NULL

#' STAR-Fusion data
#'
#' Documentation for the STAR-Fusion example data.
#'
#' @name raw_starfusion
#'
#' @section star-fusion.fusion_candidates.final.abridged.txt:
#'
#' This example data was retrieved from the STAR-Fusion github page April 13.th
#' 2016.
NULL

#' SOAPfuse data
#'
#' Documentation for the SOAPfuse example data.
#'
#' @name raw_soapfuse
#'
#' @section soapfuse_833ke_final.Fusion.specific.for.genes:
#'
#' This file has the results from a run of soapfuse-1.26 on the 833ke cell line.
#' The program was ran with the standard configuration file.
NULL

#' deFuse data
#'
#' Documentation for the deFuse example data.
#'
#' @name raw_defuse
#'
#' @section defuse_833ke_results.filtered.tsv:
#'
#' This file has the results from a run of deFuse-0.7.0 on the 833ke cell line.
#' The program was ran with the standard configuration, but with the parameter
#' span_count_threshold=5 instead of the standard 3. The resulting
#' results.filtered.tsv file was then manually filtered to only include 17
#' fusion events in the interest of saving computing time for tests and
#' examples. The original results contained 171 fusion events.
#'
#' @section reads_supporting_defuse_fusion_5267.*.fq:
#'
#' These two files, reads_supporting_defuse_fusion_5267.1.fq and
#' reads_supporting_defuse_fusion_5267.2.fq, contains the reads that support the
#' fusion event with cluster_id 5267.
#'
#' @section 5267readsAligned.bam:
#'
#' The bamfile 5267readsAligned.bam and the 5267readsAligned.bam.bai index file
#' contains the reads supporting the fusion event with cluster_id 5267 aligned
#' to the fusion sequence. It is used with plotFusionReads().
#'
NULL

#' Homo_sapiens.GRCh37.74_subset.gtf
#'
#' @name raw_Homo_sapiens.GRCh37.74
#'
#' @section Homo_sapiens.GRCh37.74_subset.gtf:
#'
#' The Homo_sapiens.GRCh37.74.gtf file is a subset version of the Ensembl
#' Homo_sapiens.GRCh37.74.gtf file, located here:
#' \url{ftp://ftp.ensembl.org/pub/release-74/gtf/homo_sapiens}. This gtf file
#' contains transcripts for the partner genes in two of the fusion transcripts
#' from the deFuse example data provided with this package: The fusion
#' transcript with cluster_id=5267, and the fusion transcript with
#' cluster_id=11759.
#'
#' The file is the result of running this command:
#'
#' # grep "ENST00000373831\|ENST00000373832\|ENST00000373833\|ENST00000398958\|ENST00000411533\|ENST00000419074\|ENST00000427469\|ENST00000429051\|ENST00000430407\|ENST00000434290\|ENST00000478232\|ENST00000486790\|ENST00000370031\|ENST00000370032\|ENST00000402983\|ENST00000420055\|ENST00000483729\|ENST00000493676\|ENST00000382073\|ENST00000299665\|ENST00000382064" Homo_sapiens.GRCh37.74.gtf > Homo_sapiens.GRCh37.74_subset.gtf
#'
#' The transcript names given in the command above are all transcripts available
#' for the genes CLEC6A, CLEC4D, HENMT1, and RCC1 in Ensembl version 74.
#'
#' @section Homo_sapiens.GRCh37.74.sqlite:
#'
#' The Homo_sapiens.GRCh37.74.sqlite file is the sqlite database that the
#' Ensembldb package creates from the corresponding gtf file. It was created
#' using this command:
#'
#' # ensDbFromGtf(
#' #   gtf = "Homo_sapiens.GRCh37.74_subset.gtf",
#' #   organism = "Homo_sapiens",
#' #   genomeVersion = "GRCh37",
#' #   version = 74)
NULL

#' Fusion5267and11759 bamfile
#'
#' Documentation for the fusion5267and11759reads.bam file containing reads
#' mapped to the region where the genes in the fusions with cluster_id=5267 and
#' cluster_id=11759 is.
#'
#' @name raw_fusion5267reads
#'
#' @section fusion5267and11759reads.bam:
#'
#' This file is the result of running these commands:
#'
#' samtools view -b original_bamfile.bam "1:28831455-28866812" "1:109189912-109205148" "12:8608225-8677832" > fusion5267and11759reads.bam
#' samtools index fusion5267and11759reads.bam fusion5267and11759reads.bam.bai
#'
#' where we extract the reads mapping to the region where we know the fusions
#' with cluster_id=5267 and cluster_id=11759 from the deFuse example data is.
#'
#' The original_bamfile.bam is from a study of the 833KE cell line by Andreas M. Hoff et al., documented in the paper [Identification of Novel Fusion Genes in Testicular Germ Cell Tumors](http://cancerres.aacrjournals.org/content/76/1/108.full).
NULL

#' protein_domains_5267 bed file
#'
#' Documentation for the protein_domains_5267.bed file containing protein
#' domains for the genes in the fusion with cluster_id=5267.
#'
#' @name raw_fusion5267proteindomains
#'
#' @section protein_domains_5267.bed:
#'
#' This file is an excerpt from a larger file that we created by:
#' - downloading domain name annotation from Pfam database (PfamA version 31)
#'   and domain region annotation from Ensembl database through BioMart API
#' - switching the domain coordinates in the protein level to these in
#'   transcript level.
NULL

#' Fusion5267and11759 bedGraph file
#'
#' Documentation for the fusion5267and11759reads.bedGraph file containing read
#' count data from the regions of the fusion event with cluster_id=5267.
#'
#' @name raw_fusion5267readsBedGraph
#'
#' @section fusion5267and11759reads.bedGraph:
#'
#' This file is the result of running this command:
#'
#' bedtools genomecov -ibam fusion5267and11759reads.bam -bg > fusion5267and11759reads.bam.bedGraph
#'
#' fusion5267and11759reads.bam has its own documentation entry for how it was
#' created.
NULL

#' Cytoband information HG19
#'
#' Cytoband information for the HG19 assembly from UCSC. Downloaded from
#' http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/
#'
#' @name raw_cytobandhg19
#'
#' @section UCSC.HG19.Human.CytoBandIdeogram.txt:
#'
#' This data is used with RCircos in plotCircle().
NULL

#' Cytoband information HG138
#'
#' Cytoband information for the HG38 assembly from UCSC. Downloaded from
#' http://hgdownload.cse.ucsc.edu/goldenpath/hg38database/
#'
#' All _alt or _random entries has been manually removed, as has the chrM entry.
#'
#' @name raw_cytobandhg38
#'
#' @section UCSC.HG38.Human.CytoBandIdeogram.txt:
#'
#' This data is used with RCircos in plotCircle().
NULL
