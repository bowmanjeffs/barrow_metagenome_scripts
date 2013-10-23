setwd('~/deming_lab/Barrow 10.2/barrow_metagenome')

mapped <- read.table('mapped_genomes_regions.txt.gz')
refs <- read.table('mapped_reads_length.txt', sep = '\t', header = F)

##### coverage plots #####

## for testing
l <- 'NC_003047'

## generate coverage plots for all genomes above the coverage cutoff,
## as output from read_length_recruit_from_sam.py

mapped_ratio <- c()

for(l in levels(mapped$V1)){
  
  strain <- refs[which(refs$V3 == l),]
  
  temp_map <- mapped[which(mapped$V1 == l),]
  temp_map <- temp_map[order(temp_map$V2),]
  
  start <- NULL
  end <- NULL
  
  for(i in seq(1:(length(temp_map$V2) - 1))){
    if(temp_map$V2[i + 1] - temp_map$V2[i] > 5000){
      start <- append(start, temp_map$V2[i])
      end <- append(end, temp_map$V2[i + 1])
    }
  }
  
  ## count gaps at beginning and end if present
  
  if(temp_map$V2[1] > 5000){
    start <- append(start, 0)
    end <- append(end, temp_map$V2[1])
  }
  
  if(max(temp_map$V4) - temp_map$V2[length(temp_map$V2)] > 5000){
    start <- append(start, temp_map$V2[length(temp_map$V2)])
    end <- append(end, max(temp_map$V4))
  }
  
  pdf(paste('coverage_plots/', l, '.pdf', sep = ''),
      width = 6,
      height = 4
      )
  
  plot(temp_map$V2,
       temp_map$V3,
       ylim = c(1, max(temp_map$V3) + 5),
       type = 'h',
       xlab = expression(paste('position x 10'^{3})),
       ylab = 'reads mapped',
       xaxt = 'n',
       xlim = c(0,max(temp_map$V4))
       )

  if(max(temp_map$V4) <= 10000){
    axis(side = 1,
         labels = seq(0,10,1),
         at = c(0,1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
    )
  }
  
  else if(max(temp_map$V4) <= 100000){
    axis(side = 1,
    labels = seq(0,100,10),
    at = c(0,10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000)
         )
  }
  
  else if(max(temp_map$V4) <= 1000000){
    axis(side = 1,
         labels = seq(0,1000,100),
         at = c(0,100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000)
    )
  } 
  
  else if(max(temp_map$V4) > 1000000){
    axis(side = 1,
         labels = seq(0,10000,1000),
         at = c(0,1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000)
    )
  } 
  
  for(i in seq(1:length(start))){
    lines(c(start[i], end[i]),
          c(max(temp_map$V3), max(temp_map$V3)),
          col = 'red',
          lwd = 4)
  }
  
  r = round(length(temp_map$V2) / max(temp_map$V4), 3)
  
  title(main = paste(l, '\n', strain$V2[1], strain$V8[1], sep = ' '),
        sub = paste('fraction mapped = ', r, sep = ''),
        cex.main = 0.8)
  
  dev.off()
}

##### compare coverage bewteen genetic elements #####

rhizobiales <- refs[which(refs$V1 == "Rhizobiales"),]

plot(refs$V5 ~ refs$V4,
     type = 'n',
     xlab = expression(paste('length x 10'^{6})),
     ylab = expression(paste('reads mapped x 10'^{3})),
     xaxt = 'n',
     yaxt = 'n')

axis(side = 1,
     labels = seq(0,14,2),
     at = seq(0,14000000,2000000))

axis(side = 2,
     labels = seq(0,50,10),
     at = seq(0,50000,10000))

points(refs$V5[which(refs$V1 != "Rhizobiales" & refs$V2 != "Candidatus Pelagibacter ubique")] ~ refs$V4[which(refs$V1 != "Rhizobiales" & refs$V2 != "Candidatus Pelagibacter ubique")],
       pch = 19,
       cex = 0.6,
       col = 'grey')

points(refs$V5[which(refs$V1 == "Rhizobiales")] ~ refs$V4[which(refs$V1 == "Rhizobiales")],
       pch = 19,
       cex = 0.6,
       col = 'black')

points(refs$V5[which(refs$V2 == "Candidatus Pelagibacter ubique")] ~ refs$V4[which(refs$V2 == "Candidatus Pelagibacter ubique")],
       pch = 3,
       cex = 0.6,
       col = 'black')

legend('topright',
       legend = c('non-Rhizobiales', 'Rhizobiales', 'Pelagibacter ubique'),
       col = c('grey', 'black', 'black'),
       pch = c(19,19,3),
       cex = 0.8)

## just Rhizobiales, plasmids vs chromosomes

plot(refs$V5 ~ refs$V4,
     type = 'n',
     xlab = expression(paste('length x 10'^{6})),
     ylab = expression(paste('reads mapped x 10'^{3})),
     xaxt = 'n',
     yaxt = 'n')

axis(side = 1,
     labels = seq(0,14,2),
     at = seq(0,14000000,2000000))

axis(side = 2,
     labels = seq(0,50,10),
     at = seq(0,50000,10000))

points(refs$V5[which(refs$V1 == "Rhizobiales" & refs$V8 == 'chromosome')] ~ refs$V4[which(refs$V1 == "Rhizobiales" & refs$V8 == 'chromosome')],
       pch = 19,
       cex = 0.8,
       col = 'black')

points(refs$V5[which(refs$V1 == "Rhizobiales" & refs$V8 == 'plasmid')] ~ refs$V4[which(refs$V1 == "Rhizobiales" & refs$V8 == 'plasmid')],
       pch = 1,
       cex = 0.8,
       col = 'black')

legend('topright',
       legend = c('chromosome', 'plasmid'),
       col = c('black', 'black'),
       pch = c(19,1),
       cex = 0.8)

## just rhizobiales plasmids

plot(refs$V5[which(refs$V1 == "Rhizobiales" & refs$V8 == 'plasmid')] ~ refs$V4[which(refs$V1 == "Rhizobiales" & refs$V8 == 'plasmid')],
     type = 'n',
     xlab = expression(paste('length x 10'^{6})),
     ylab = expression(paste('reads mapped x 10'^{3})),
     xaxt = 'n',
     yaxt = 'n')

axis(side = 1,
     labels = seq(0,3,0.5),
     at = seq(0,3000000,500000))

axis(side = 2,
     labels = seq(0,10,2),
     at = seq(0,10000,2000))

points(refs$V5[which(refs$V1 == "Rhizobiales" & refs$V8 == 'plasmid')] ~ refs$V4[which(refs$V1 == "Rhizobiales" & refs$V8 == 'plasmid')],
       pch = 1,
       cex = 0.8,
       col = 'black')


       