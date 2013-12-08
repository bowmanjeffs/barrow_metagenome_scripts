setwd('~/deming_lab/Barrow 10.2/barrow_metagenome')

mapped <- read.table('mapped_genomes_regions.txt.gz')
refs <- read.table('mapped_reads_length.txt', sep = '\t', header = T)

## find top 5 % covered

t <- ceiling(0.05 * length(refs$cov))

refs_top_cov <- refs[order(refs$cov, decreasing = T)[1:t],]

##### coverage plots #####

## for testing
l <- 'NC_011986'

## generate coverage plots for all genomes above the coverage cutoff,
## as output from read_length_recruit_from_sam.py

mapped_ratio <- c()

for(l in refs_top_cov$id){
  
  strain <- refs[which(refs$id == l),]
  
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
  
  if(strain$length[1] - temp_map$V2[length(temp_map$V2)] > 5000){
    start <- append(start, temp_map$V2[length(temp_map$V2)])
    end <- append(end, strain$length[1])
  }
  
  pdf(paste('coverage_plots/', l, '.pdf', sep = ''),
      width = 6,
      height = 4
      )
  
  plot(temp_map$V2,
       temp_map$V3,
       ylim = c(1, max(temp_map$V3) + (max(temp_map$V3) * 0.07)),
       type = 'h',
       #log = 'y',
       #ylog = TRUE,
       xlab = expression(paste('position x 10'^{3})),
       ylab = 'reads mapped',
       xaxt = 'n',
       xlim = c(0,strain$length[1])
       )

  if(strain$length[1] <= 10000){
    axis(side = 1,
         labels = seq(0,10,1),
         at = c(0,1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
    )
  }
  
  else if(strain$length[1] <= 100000){
    axis(side = 1,
    labels = seq(0,100,10),
    at = c(0,10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000)
         )
  }
  
  else if(strain$length[1] <= 1000000){
    axis(side = 1,
         labels = seq(0,1000,100),
         at = c(0,100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000)
    )
  } 
  
  else if(strain$length[1] > 1000000){
    axis(side = 1,
         labels = seq(0,10000,1000),
         at = c(0,1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000)
    )
  } 
  
  for(i in seq(1:length(start))){
    rect(start[i],
         max(temp_map$V3),
         end[i],
         max(temp_map$V3) + (0.05 * max(temp_map$V3)),
         col = 'red',
         border = NA)
  }
  
  r = round(length(temp_map$V2) / strain$length[1], 3)
  
  title(main = paste(l, '\n', strain$strain[1], strain$type[1], sep = ' '),
        sub = paste('fraction mapped = ', r, sep = ''),
        cex.main = 0.8)
  
  dev.off()
}

##### compare coverage bewteen genetic elements #####

## coverage v breadth

col = colorRampPalette(c('blue','red'))(2)

rhizobiales <- refs[which(refs$order == "Rhizobiales"),]

plot(refs$cov ~ refs$breadth,
     type = 'n',
     xlab = expression(paste('breadth')),
     ylab = expression(paste('coverage'))
)

points(refs$cov[which(refs$order != "Rhizobiales" & refs$strain != "Candidatus Pelagibacter ubique")] ~ refs$breadth[which(refs$order != "Rhizobiales" & refs$strain != "Candidatus Pelagibacter ubique")],
       pch = 1,
       cex = 0.6,
       col = col[refs$type])

points(refs$cov[which(refs$order == "Rhizobiales")] ~ refs$breadth[which(refs$order == "Rhizobiales")],
       pch = 19,
       cex = 0.6,
       col = col[refs$type])

points(refs$cov[which(refs$strain == "Candidatus Pelagibacter ubique")] ~ refs$breadth[which(refs$strain == "Candidatus Pelagibacter ubique")],
       pch = 3,
       cex = 0.6,
       col = col[refs$type])

legend('topright',
       legend = c('non-Rhizobiales', 'Rhizobiales', 'Pelagibacter ubique'),
       col = c('grey', 'black', 'black'),
       pch = c(19,19,3),
       cex = 0.8)

## coverage v length

opar <- par()

pdf(file = 'fig_1A.pdf',
    width = 7,
    height = 5)

par(mar = c(5, 6, 4, 1) + 0.1)

rhizobiales <- refs[which(refs$order == "Rhizobiales"),]

plot(refs$mapped ~ refs$length,
     type = 'n',
     xlab = expression(paste('length x 10'^{6})),
     ylab = expression(paste('reads mapped x 10'^{3})),
     xaxt = 'n',
     yaxt = 'n',
     cex.axis = 1.2)

axis(side = 1,
     labels = seq(0,14,2),
     at = seq(0,14000000,2000000))

axis(side = 2,
     labels = seq(0,50,10),
     at = seq(0,50000,10000))

points(refs$mapped[which(refs$order != "Rhizobiales" & refs$strain != "Candidatus Pelagibacter ubique")] ~ refs$length[which(refs$order != "Rhizobiales" & refs$strain != "Candidatus Pelagibacter ubique")],
       pch = 19,
       cex = 0.6,
       col = 'blue')

points(refs$mapped[which(refs$order == "Rhizobiales")] ~ refs$length[which(refs$order == "Rhizobiales")],
       pch = 19,
       cex = 0.6,
       col = 'red')

points(refs$mapped[which(refs$strain == "Candidatus Pelagibacter ubique")] ~ refs$length[which(refs$strain == "Candidatus Pelagibacter ubique")],
       pch = 3,
       cex = 0.6,
       col = 'blue')

legend('topright',
       legend = c('non-Rhizobiales', 'Rhizobiales', 'Pelagibacter ubique'),
       col = c('blue', 'red', 'blue'),
       pch = c(19,19,3),
       cex = 1)

dev.off()

## just Rhizobiales, plasmids vs chromosomes

pdf(file = 'fig_1B.pdf',
    width = 7,
    height = 5)

par(mar = c(5, 6, 4, 1) + 0.1)

plot(refs$mapped ~ refs$length,
     type = 'n',
     xlab = expression(paste('length x 10'^{6})),
     ylab = expression(paste('reads mapped x 10'^{3})),
     xaxt = 'n',
     yaxt = 'n',
     cex.axis = 1)

axis(side = 1,
     labels = seq(0,14,2),
     at = seq(0,14000000,2000000))

axis(side = 2,
     labels = seq(0,50,10),
     at = seq(0,50000,10000))

points(refs$mapped[which(refs$order == "Rhizobiales" & refs$type == 'chromosome')] ~ refs$length[which(refs$order == "Rhizobiales" & refs$type == 'chromosome')],
       pch = 19,
       cex = 0.8,
       col = 'red')

points(refs$mapped[which(refs$order == "Rhizobiales" & refs$type == 'plasmid')] ~ refs$length[which(refs$order == "Rhizobiales" & refs$type == 'plasmid')],
       pch = 19,
       cex = 0.8,
       col = 'blue')

legend('topright',
       legend = c('chromosome', 'plasmid'),
       col = c('red', 'blue'),
       pch = c(19,19),
       cex = 1)

dev.off()

## t-test between these two groups ##

t.test(refs$cov[which(refs$order == "Rhizobiales" & refs$type == 'chromosome')],
       refs$cov[which(refs$order == "Rhizobiales" & refs$type == 'plasmid')])

t.test(refs$breadth[which(refs$order == "Rhizobiales" & refs$type == 'chromosome')],
       refs$breadth[which(refs$order == "Rhizobiales" & refs$type == 'plasmid')])

## just rhizobiales plasmids

plot(refs$mapped[which(refs$order == "Rhizobiales" & refs$type == 'plasmid')] ~ refs$length[which(refs$order == "Rhizobiales" & refs$type == 'plasmid')],
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

points(refs$mapped[which(refs$order == "Rhizobiales" & refs$type == 'plasmid')] ~ refs$length[which(refs$order == "Rhizobiales" & refs$type == 'plasmid')],
       pch = 1,
       cex = 0.8,
       col = 'black')

##### look at plasmids corresponding to best covered chromosomes

t <- ceiling(0.01 * length(refs$cov))

refs_0.01_top_cov <- refs[order(refs$cov, decreasing = T)[1:t],]

ids <- refs_0.01_top_cov$taxid[which(refs_0.01_top_cov$type == 'chromosome' & refs_0.01_top_cov$order == 'Rhizobiales')]

select <- refs[which(refs$taxid %in% ids),]

select <- select[-grep('Brucella|Bartonella|Candidatus|Liberibacter', select$strain),]

select <- select[order(select$strain, select$type),]

opar <- par()

par(mar = c(10, 7, 2, 2) + 0.1)

color <- NULL

for(c in select$type){
  if(c == 'chromosome'){
    color <- append(color, 'red')
  }else(color <- append(color, 'blue'))
}

labels <- select$strain[-which(duplicated(select$strain))]

text_x <- c(2,5.5,10.5,17) * 1.2

barplot(select$cov,
        col = color,
        ylab = 'Coverage',
        ylim = c(0,max(select$cov) * 1.1))

box()

legend('topleft',
       legend = c('chromosome', 'plasmid'),
       fill = c('red', 'blue'))

text(x = text_x,
     y = -0.3,
     labels = labels,
     xpd = TRUE,
     srt = 45,
     adj = 1,
     cex = 1
     )

par(opar)

##### barplot of top chromosomes #####

species_bp <- barplot(refs_top_cov$cov[which(refs_top_cov$type == 'chromosome')][1:10],
        ylab = 'Coverage',
        xaxt = 'none',
        cex.axis = 0.9)

labels <- refs_top_cov$strain[which(refs_top_cov$type == 'chromosome')][1:10]

text(species_bp[,1],
     y = -max(refs_top_cov$cov) * 0.1,
     labels = labels,
     xpd=TRUE,
     srt=45,
     adj = 1,
     cex = 0.9
)

##### boxplot of chromosome vs plasmids #####

pdf(file = 'fig_1C.pdf',
    width = 6,
    height = 5)

par(mar = c(5, 2, 4, 5) + 0.1)

boxplot(rhizobiales$breadth ~ rhizobiales$type,
        notch = T,
        yaxt = 'n')

axis(side = 4)

mtext('Breadth',
      side = 4,
      line = 2.5)

dev.off()

pdf(file = 'fig_1D.pdf',
    width = 6,
    height = 5)

par(mar = c(5, 2, 4, 5) + 0.1)

boxplot(rhizobiales$cov ~ rhizobiales$type,
        notch = T,
        yaxt = 'n')

axis(side = 4)

mtext('Coverage',
      side = 4,
      line = 2.5)

dev.off()

par(opar)

##### look at plasmids and chromosomes for strains with plasmids not fitting pattern ######

rhizobiales <- rhizobiales[order(rhizobiales$cov, decreasing = T),]

top_plasmids <- rhizobiales[which(rhizobiales$type == 'plasmid')[1:4],]

top_plasmids_w_chrom <- rhizobiales[which(rhizobiales$taxid %in% top_plasmids$taxid),]

top_plasmids_w_chrom <- top_plasmids_w_chrom[order(top_plasmids_w_chrom$taxid, top_plasmids_w_chrom$type),]

pdf(file = 'fig_2.pdf',
    width = 8,
    height = 6)

par(mar = c(10, 7, 2, 2) + 0.1)

color <- NULL

for(c in top_plasmids_w_chrom$type){
  if(c == 'chromosome'){
    color <- append(color, 'red')
  }else(color <- append(color, 'blue'))
}

labels <- top_plasmids_w_chrom$strain[-which(duplicated(top_plasmids_w_chrom$strain))]

text_x <- c(3.5,11.5,19,27.5) * 1.2

barplot(top_plasmids_w_chrom$cov,
        col = color,
        ylab = 'Coverage',
        ylim = c(0,max(top_plasmids_w_chrom$cov) * 1.1))

box()

legend('topright',
       legend = c('chromosome', 'plasmid'),
       fill = c('red', 'blue'))

text(x = text_x,
     y = -0.3,
     labels = labels,
     xpd = TRUE,
     srt = 45,
     adj = 1,
     cex = 1
)

dev.off()

##### get best covered Rhizobiales strains, all elements ##

pos_mapped <- cbind(rhizobiales$taxid, rhizobiales$length, rhizobiales$length * rhizobiales$breadth)

summed_pos <- data.frame(tapply(pos_mapped[,3], pos_mapped[,1], sum))

summed_length <- data.frame(tapply(pos_mapped[,2], pos_mapped[,1], sum))

summed <- cbind(summed_length$row.names, summed_pos[,1] / summed_length[,1])

summed <- as.matrix(summed[order(summed[,1], decreasing = T),])

top_five <- rhizobiales[which(rhizobiales$taxid %in% row.names(summed)[1:5]),]

top_five <- top_five[order(top_five$taxid, top_five$type),]

write.table(top_five$id, "ncbi_id_top_five_rhizobiales.txt", quote = F, col.names = F, row.names = F)