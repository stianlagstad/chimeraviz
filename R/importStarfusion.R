#' read the results from runnning FusionInspector
#'
#' @export
importFusionInspector <- function (filename='FusionInspector-inspect/finspector.igv.FusionJuncSpan',
                                  limit=Inf) {
    ## Try to read the FusionInspector report.
    ## Within one 'scaffold', only keep the leftmost and rightmost coordinate
    report <- withCallingHandlers({
        col_types_fusioninspector = readr::cols_only(
          "#scaffold" = col_character(),
          "fusion_break_name" = col_skip(),
          "break_left" = col_integer(),
          "break_right" = col_integer(),
          "num_junction_reads" = col_skip(),
          "num_spanning_frags" = col_skip(),
          "spanning_frag_coords" = col_skip()
          )
        if (missing(limit)) {
            readr::read_tsv(
              file = filename,
              col_types = col_types_fusioninspector)
        } else {
            readr::read_tsv(
              file = filename,
              col_types = col_types_fusioninspector,
              n_max = limit)
        }
    },
      error = function(cond) {
          message(paste0("Reading ", filename, " caused an error: ", cond[[1]]))
          stop(cond)
      },
      warning = function(cond) {
          message(paste0("Reading ", filename, " caused a warning: ", cond[[1]]))
          warning(cond)
      })
    colnames(report)[1] <- 'id'
    d <-
      data.frame(id=with(report, id[ !duplicated(id) ]),
                 starts=with(report, tapply(start, id, min)),
                 ends=with(report, tapply(end, id, max)))
    rownames(d) <- d$id
    d
}                                        #importFusionInspector

#' Import results from a STAR-Fusion run into a list of Fusion objects.
#'
#' A function that imports the results from a STAR-Fusion run, typically from
#' a star-fusion.fusion_candidates.final.abridged file, into a list of Fusion
#' objects.
#'
#' @param filename Filename for the STAR-Fusion
#' star-fusion.fusion_candidates.final.abridged results file.
#' @param genomeVersion Which genome was used in mapping (hg19, hg38, etc.).
#' @param limit A limit on how many lines to read.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' starfusionData <- system.file(
#'   "extdata",
#'   "star-fusion.fusion_candidates.final.abridged.txt",
#'   package = "chimeraviz")
#' fusions <- importStarfusion(starfusionData, "hg19", 3)
#' # This should import a list of 3 fusions described in Fusion objects.
#'
#' @export
importStarfusion <- function (filename, genomeVersion, limit=Inf,
                              useFusionInspector=FALSE) {
  # Is the genome version valid?
  validGenomes <- c("hg19", "hg38", "mm10")
  if (is.na(match(tolower(genomeVersion), tolower(validGenomes)))) {
    stop("Invalid genome version given")
  }

  # If the limit is set, is the value valid?
  if (missing(limit) == FALSE) {
    if (is.numeric(limit) == FALSE || limit <= 0) {
      stop("limit must be a numeric value bigger than 0")
    }
  }

  # Try to read the fusion report
  report <- withCallingHandlers(
    {
      col_types_starfusion = readr::cols_only(
##        "#FusionName" = col_skip(),
        "#FusionName" = col_character(), 
        "JunctionReadCount" = col_integer(),
        "SpanningFragCount" = col_integer(),
        "SpliceType" = col_skip(),
        "LeftGene" = col_character(),
        "LeftBreakpoint" = col_character(),
        "RightGene" = col_character(),
        "RightBreakpoint" = col_character(),
        "LargeAnchorSupport" = col_character(),
        "LeftBreakDinuc" = col_character(),
        "LeftBreakEntropy" = col_number(),
        "RightBreakDinuc" = col_character(),
        "RightBreakEntropy" = col_number(),
        "FFPM" = col_number(),
        "PROT_FUSION_TYPE" = col_character()
      )
      if (missing(limit)) {
        # Read all lines
        readr::read_tsv(
          file = filename,
          col_types = col_types_starfusion
        )
      } else {
        # Only read up to the limit
        readr::read_tsv(
          file = filename,
          col_types = col_types_starfusion,
          n_max = limit
        )
      }
    },
    error = function(cond) {
      message(paste0("Reading ", filename, " caused an error: ", cond[[1]]))
      stop(cond)
    },
    warning = function(cond) {
      message(paste0("Reading ", filename, " caused a warning: ", cond[[1]]))
      warning(cond)
    }
  )

  if(useFusionInspector) {
      l <- limit+10
      fi.table <- importFusionInspector(limit=l)
  }

  
  # Set variables
  id                    <- NA
  inframe               <- NA
  fusionTool            <- ifelse(useFusionInspector, "starfusion+fusioninspector", "starfusion")
  spanningReadsCount    <- NA
  splitReadsCount       <- NA
  junctionSequence      <- NA
  fusionReadsAlignment  <- NA

  # List to hold all Fusion objects
  fusionList <- vector("list", dim(report)[1])

  # Iterate through each line in the .tsv file
  for (i in 1:dim(report)[1]) {

    # Import starfusion-specific fields
    fusionToolSpecificData <- list()
    fusionToolSpecificData[["FusionName"]] = report[[i, "#FusionName"]]
    fusionToolSpecificData[["LargeAnchorSupport"]] = report[[i, "LargeAnchorSupport"]]
    fusionToolSpecificData[["LeftBreakDinuc"]] = report[[i, "LeftBreakDinuc"]]
    fusionToolSpecificData[["LeftBreakEntropy"]] = report[[i, "LeftBreakEntropy"]]
    fusionToolSpecificData[["RightBreakDinuc"]] = report[[i, "RightBreakDinuc"]]
    fusionToolSpecificData[["RightBreakEntropy"]] = report[[i, "RightBreakEntropy"]]
    fusionToolSpecificData[["FFPM"]] = report[[i, "FFPM"]]

    if(!is.null(report$PROT_FUSION_TYPE)) {
        ## only available when loading from coding_effect file (i.e.
        ## star-fusion.fusion_predictions.abridged.annotated.coding_effect.tsv)
        ft <- report[[i, "PROT_FUSION_TYPE"]]
        fusionToolSpecificData[["inframe"]] <- ft
        if (ft == 'INFRAME')
          inframe <- TRUE
        if (ft == 'FRAMESHIFT')
          inframe <- FALSE
    }
    
    # id for this fusion
    id <- as.character(i)

    leftBreakPoint <- unlist(strsplit(report[[i, "LeftBreakpoint"]], split=":"))
    rightBreakPoint <- unlist(strsplit(report[[i, "RightBreakpoint"]], split=":"))

    # Chromosome names
    chromosomeA <- leftBreakPoint[1]
    chromosomeB <- rightBreakPoint[1]

    # Breakpoints
    breakpointA <- as.numeric(leftBreakPoint[2])
    breakpointB <- as.numeric(rightBreakPoint[2])

    # Strand
    strandA <- leftBreakPoint[3]
    strandB <- rightBreakPoint[3]

    # Number of supporting reads
    splitReadsCount <- report[[i, "JunctionReadCount"]]
    spanningReadsCount <- report[[i, "SpanningFragCount"]]

    # STAR-Fusion doesn't provide the fusion sequence
    junctionSequenceA <- Biostrings::DNAString()
    junctionSequenceB <- Biostrings::DNAString()

    # Gene names
    geneNames1 <- unlist(strsplit(report[[i, "LeftGene"]], split="\\^"))
    geneNames2 <- unlist(strsplit(report[[i, "RightGene"]], split="\\^"))
    nameA <- geneNames1[1]
    nameB <- geneNames2[1]

    # Ensembl ids
    ensemblIdA <- geneNames1[2]
    ensemblIdB <- geneNames2[2]

    if (useFusionInspector) {
        ## override things to go with the FI's virtual 'mini genome', but
        ## first save them in fusionToolSpecificData:
        fusionToolSpecificData[['orig_chromosomeA']] <- chromosomeA
        fusionToolSpecificData[['orig_breakpointA']] <- breakpointA
        fusionToolSpecificData[['orig_strandA']] <- strandA
        fusionToolSpecificData[['orig_chromosomeB']] <- chromosomeB
        fusionToolSpecificData[['orig_breakpointB']] <- breakpointB
        fusionToolSpecificData[['orig_strandB']] <- strandB

        fusion.name <- report[[i, "#FusionName"]]
        if(fusion.name %in% fi.table$id) {
            chromosomeA <- fusion.name
            chromosomeB <- fusion.name
            strandA <- '+'
            strandB <- '+'
            breakpointA <- fi.table[fusion.name, 'start']
            breakpointB <- fi.table[fusion.name, 'end']
        } else {
            warning(sprintf("Could not find location for fusion %s
on virtual contig, ignoring it (line %d)", fusion.name, i+1))
            next
        }
    }

    # PartnerGene objects
    geneA <- new(Class = "PartnerGene",
                 name = nameA,
                 ensemblId = ensemblIdA,
                 chromosome = chromosomeA,
                 breakpoint = breakpointA,
                 strand = strandA,
                 junctionSequence = junctionSequenceA,
                 transcripts = GenomicRanges::GRangesList())
    geneB <- new(Class = "PartnerGene",
                 name = nameB,
                 ensemblId = ensemblIdB,
                 chromosome = chromosomeB,
                 breakpoint = breakpointB,
                 strand = strandB,
                 junctionSequence = junctionSequenceB,
                 transcripts = GenomicRanges::GRangesList())

    fusionList[[i]] <- new(Class = "Fusion",
                           id = id,
                           fusionTool = fusionTool,
                           genomeVersion = genomeVersion,
                           spanningReadsCount = spanningReadsCount,
                           splitReadsCount = splitReadsCount,
                           fusionReadsAlignment = Gviz::AlignmentsTrack(),
                           geneA = geneA,
                           geneB = geneB,
                           inframe = inframe,
                           fusionToolSpecificData = fusionToolSpecificData)
  }                                     #for i

  # Return the list of Fusion objects
  fusionList
}
