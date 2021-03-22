setwd("~/age/era/motfi")
library(MultiAssayExperiment)
library(ELMER)
met <- read.csv("~/age/era/data/meth.csv", row.names=1)

exp <- read.csv("~/age/era/data/exp.csv", row.names=1)
colnames(exp) <- colnames(met)

assay <- c(rep("DNA methylation", ncol(met)),
           rep("Gene expression", ncol(exp)))

primary <- c(colnames(met),colnames(exp))

colname <- c(colnames(met),colnames(exp))

sampleMap <- data.frame(assay,primary,colname)

distal.probes <- get.feature.probe(genome = "hg38", 
                                   met.platform = "EPIC")
                                   #rm.chr = paste0("chr",c(1:21,"X","Y")))

colData <- data.frame(sample = colnames(met))
colnames(colData)='primary'

rownames(colData) <- colnames(met)

colData$Group=c('n','t','n','t','n','t','n','t')

mae <- createMAE(exp = exp, 
                 met = met,
                 save = TRUE,
                 filter.probes = distal.probes,
                 colData = colData,
                 sampleMap = sampleMap,
                 linearize.exp = TRUE,
                 save.filename = "mae.rda",
                 met.platform = "EPIC",
                 genome = "hg38",
                 TCGA = FALSE)

sig.diff <- get.diff.meth(data = mae, 
              group.col = "Group",
              group1 =  "t",
              group2 = "n",
              minSubgroupFrac = 1, # if supervised mode set to 1
              sig.dif = .3,
              diff.dir = "hypo", # hyper = Search for hypomethylated probes in group 1
              cores = 4, 
              dir.out ="result", 
              pvalue = 0.9)

enriched.motif <- get.enriched.motif(data = mae,
                                     probes = sig.diff$probe, 
                                     dir.out = "result", 
                                     label = "hyper",
                                     min.incidence = 10,
                                     lower.OR = 4)
TF <- get.TFs(data = mae, 
              group.col = "Group",
              group1 =  "n",
              group2 = "t",
              mode = "unsupervised",
              enriched.motif = enriched.motif,
              dir.out = "result", 
              cores = 1, 
              label = "hyper")


nearGenes <- GetNearGenes(data = mae, 
                          probes = sig.diff$probe, 
                          numFlankingGenes = 20) # 10 upstream and 10 dowstream genes

Hyper.pair <- get.pair(data = mae,
                      group.col = "Group",
                      group1 =  "n",
                      group2 = "t",
                      nearGenes = nearGenes,
                      mode = "unsupervised",
                      permu.dir = "result/permu",
                      permu.size = 100, # Please set to 100000 to get significant results
                      raw.pvalue = 0.1,   
                      Pe = 0.01, # Please set to 0.001 to get significant results
                      filter.probes = TRUE, # See preAssociationProbeFiltering function
                      filter.percentage = 0.05,
                      filter.portion = 0.3,
                      dir.out = "result",
                      cores = 1,
                      label = "hyper")
cpg='cg08878436'
scatter.plot(data = mae,
             byProbe = list(probe = c(cpg), numFlankingGenes = 20), 
             category = "Group", 
             lm = TRUE, # Draw linear regression curve
             save = FALSE) 

load("result/getTF.hypo.TFs.with.motif.pvalue.rda")
motif <- colnames(TF.meth.cor)[1]
TF.rank.plot(motif.pvalue = TF.meth.cor, 
             motif = motif,
             save = FALSE) 

