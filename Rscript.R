#' Calculate PWM enrichment upstream of transcription start sites (TSS)
#' This function calculates the enrichment of a position weight matrix (PWM)
#' in promoter regions around transcription start sites (TSS), with optional 
#' @param pwm A `PWMatrix` object.
#' @param genome Path to the genome FASTA file.
#' @param upstream Integer, number of bp upstream of TSS (default: 1000).
#' @param downstream Integer, number of bp downstream of TSS (default: 500).
#' @param Txdb A `TxDb` object containing transcript annotations.
#' @param bedfile Optional path to a BED file of peaks for filtering promoter regions.
enrich.in.TSS <- function(pwm,
                          genome,
                          upstream = 1000,
                          downstream = 500,
                          Txdb,
                          bedfile = NULL) {
  seqlevels(Txdb) <- c(paste0("Chr", 1:5), "Mt", "Pt")
  seqlevels(Txdb, pruning.mode = "coarse") <- seqlevels(Txdb)[
    stringr::str_detect(seqlevels(Txdb), "[cC][hH][rR][0-9]+")
  ]

  genome <- Biostrings::readDNAStringSet(genome)

  promoter <- GenomicFeatures::promoters(Txdb,
                                         upstream = upstream,
                                         downstream = downstream) %>%
    GenomicRanges::trim()

  promoter <- promoter[
    promoter$tx_name %>%
      gsub(".[2-9]$", "", .) %>%
      grep(pattern = ".1$", .)
  ]
  promoter$tx_name <- gsub("\\.[1-9]$", "", promoter$tx_name)
  if (!is.null(bedfile)) {
    peakfile <- ChIPseeker::readPeakFile(bedfile)
    seqlevels(peakfile) <- c(paste0("Chr", 1:5), "Mt", "Pt")
    seqlevels(peakfile, pruning.mode = "coarse") <- seqlevels(Txdb)[
      stringr::str_detect(seqlevels(Txdb), "[cC][hH][rR][0-9]+")
    ]
    
    peakanno <- ChIPseeker::annotatePeak(peakfile,
                                         tssRegion = c(-upstream, downstream),
                                         TxDb = Txdb,
                                         annoDb = "org.At.tair.db")
    data <- as.data.frame(peakanno)
    data$annotation <- gsub("Promoter.*$", "Promoter", data$annotation)
    geneid <- data[data$annotation == "Promoter", "geneId"]
    promoter <- promoter[promoter$tx_name %in% geneid]
  }

  seqs <- genome[promoter]
  enrich <- pattern.enrich(seq = seqs, pwm = pwm, type = "pwm")
  
  return(enrich)
}

