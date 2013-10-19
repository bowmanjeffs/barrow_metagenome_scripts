setwd('~/deming_lab/Barrow 10.2/barrow_metagenome')

mapped <- read.table('mapped_genomes_regions.txt')

l <- 'NC_007205'

mapped_ratio <- c()

for(l in levels(mapped$V1)){
  
  temp_map <- mapped[which(mapped$V1 == l),]
  breaks <- ceiling(max(temp_map$V2) / 5000)
  temp_hist <- hist(temp_map$V2,
                       breaks = breaks,
                       plot = F)
  r  <- length(which(temp_hist$counts == 0)) / breaks
  
  pdf(paste(l, '.pdf', sep = ''),
      width = 6,
      height = 4
      )
  
  cov <- (temp_hist$counts * 96) / 5000
  
  plot(cov ~ temp_hist$mids,
       ylim = c(1, max(cov) + 5),
#       log = 'y',
       type = 'h',
       xlab = expression(paste('position x 10'^{6})),
       ylab = 'coverage per 5000 bp bin',
       xaxt = 'n',
       )
  
  axis(side = 1,
       labels = seq(0,1,1),
       at = c(0,1000000)
  )

  temp_gaps <- temp_hist$mids[which(temp_hist$counts == 0)]

  points(temp_gaps, c(rep(max(cov) + 5, length(temp_gaps))),
       pch = '|',
       cex = 0.6,
       col = 'red')
  title(main = l)
  title(sub = paste('gap ratio = ', round(r,2)))
  dev.off()
}

ratio <- read.table('top_reads_v_length.txt', sep = '\t')

t.test(ratio[,5] ~ ratio[,3])