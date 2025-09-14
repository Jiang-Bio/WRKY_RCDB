#读取narrowpeak文件并取peak交集做下游分析
mergePeakFiles <- function(path, pattern = "narrowPeak$") {
  peakfiles <- list.files(path = path, pattern = pattern, full.names = TRUE)
  if (length(peakfiles) == 0) {
    stop("No peak files found in the given path.")
  }

  g <- rtracklayer::import(peakfiles[1])

  if (length(peakfiles) > 1) {
    for (i in 2:length(peakfiles)) {
      g2 <- rtracklayer::import(peakfiles[i])
      a <- findOverlaps(g, g2)
      g2 <- g2[-unique(a@to)]
      g <- c(g, g2)
    }
  }
  return(g)
}
merged_publicdap <- mergePeakFiles("2024_WRKY/5_CALLPEAK/Public_DAP_ampDAP/DAP/callpeak/")
merged_dap <- mergePeakFiles("2024_WRKY/5_CALLPEAK/DAP_callpeak_minus_pixHalo/")
  a <- findOverlapsOfPeaks(unique(merged_dap),unique(merged_publicdap))
  dat<-c("this_study" = a$peaklist$unique.merged_dap %>% length(), 
         "previous_study" = a$peaklist$unique.merged_publicdap %>% length(), 
         "this_study&previous_study" =a$peaklist$`unique.merged_dap///unique.merged_publicdap` %>% length()
  )
  p1 <- plot(euler(dat),
             quantities = T,
             fills=list(fill=c("#FF6E72","#F7BC66","#61CA9D")),
             edges = list(col=c("black"),lty=1,lwd=2),main = "DAP")

merged_publicampdap <- mergePeakFiles("2024_WRKY/5_CALLPEAK/Public_DAP_ampDAP/ampDAP/callpeak/")
merged_ampdap <- mergePeakFiles("2024_WRKY/5_CALLPEAK/AMPdap_callpeak_minus_pixHalo")
  a <- findOverlapsOfPeaks(unique(merged_ampdap),unique(merged_publicampdap))
  dat<-c("this_study" = a$peaklist$merged_ampdap %>% length(), 
         "previous_study" = a$peaklist$unique.merged_publicampdap %>% length(), 
         "this_study&previous_study" =a$peaklist$`unique.merged_ampdap///unique.merged_publicampdap` %>% length()
  )
  p2 <- plot(euler(dat),
             quantities = T,
             fills=list(fill=c("#FF6E72","#F7BC66","#61CA9D")),
             edges = list(col=c("black"),lty=1,lwd=2),main = "DAP")

