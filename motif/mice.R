setwd("~/age/era/motif")
library(MultiAssayExperiment)
library(ELMER)
library(gridExtra)
library(dplyr)
met <- read.csv("~/age/era/data/vein_meth.csv", row.names=1)

exp <- read.csv("~/age/era/data/expf.csv", row.names=1)
#exp=exp[colnames(met)]
colnames(exp)=colnames(met) # change

assay <- c(rep("DNA methylation", ncol(met)),
           rep("Gene expression", ncol(exp)))

primary <- c(colnames(met),colnames(exp))

colname <- c(colnames(met),colnames(exp))

sampleMap <- data.frame(assay,primary,colname)

distal.probes <- get.feature.probe(genome = "hg38", 
                                   met.platform = "EPIC",
                                   rm.chr = paste0("chr",c("X","Y")))

colData <- data.frame(sample = colnames(met))
colnames(colData)='primary'

rownames(colData) <- colnames(met)

colData$Group=c('Aged','Treated','Aged','Treated','Aged','Treated','Aged','Treated')

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
              group1 =  "Treated",
              group2 = "Aged",
              minSubgroupFrac = 1, # if supervised mode set to 1
              sig.dif = .3,
              diff.dir = "hypo", # hyper = Search for hypomethylated probes in group 1
              cores = 4, 
              dir.out ="result", 
              pvalue = 0.9)

enriched.motif <- get.enriched.motif(data = mae,
                                     probes = sig.diff$probe, 
                                     dir.out = "result", 
                                     label = "hypo",
                                     min.incidence = 10,
                                     lower.OR = 2)
TF <- get.TFs(data = mae, 
              group.col = "Group",
              group1 =  "Treated",
              group2 = "Aged",
              mode = "unsupervised",
              enriched.motif = enriched.motif,
              dir.out = "result", 
              cores = 1, 
              label = "hypo")


nearGenes <- GetNearGenes(data = mae, 
                          probes = sig.diff$probe, 
                          numFlankingGenes = 5) # 10 upstream and 10 dowstream genes

Hypo.pair <- get.pair(data = mae,
                      group.col = "Group",
                      group1 =  "Treated",
                      group2 = "Aged",
                      nearGenes = nearGenes,
                      raw.pvalue = .1,
                      permu.dir = "result/permu1",
                      permu.size = 100, # Please set to 100000 to get significant results
                      mode = "unsupervised",
                      Pe = .1, # Please set to 0.001 to get significant results
                      filter.probes = FALSE, # See preAssociationProbeFiltering function
                      #filter.percentage = 0.03,
                      #filter.portion = 0.3,
                      dir.out = "result",
                      cores = 1,
                      label = "hypo")
cpg='cg23468453'
scatter.plot(data = mae,
             byProbe = list(probe = c(cpg), numFlankingGenes = 10), 
             category = "Treated", 
             lm = TRUE, # Draw linear regression curve
             save = FALSE) 
m=(Hypo.pair %>% filter(startsWith(Symbol,'RP')))$Probe
g=(Hypo.pair %>% filter(startsWith(Symbol,'RP')))$GeneID
for (i in c(1:length(m))){
p=scatter.plot(data = mae,
             byPair = list(probe = m[i], 
                           gene = g[i]), 
             category = "Group", save = TRUE, lm_line = TRUE) 
p
}
pair <- read.csv("result/getPair.hypo.pairs.significant.csv")
schematic.plot(pair = pair, 
               data = mae,
               group.col = "Group",
               byProbe = "cg10844712",
               save = FALSE)
load("result/getTF.hypo.TFs.with.motif.pvalue.rda")
motif <- colnames(TF.meth.cor)[1]
TF.rank.plot(motif.pvalue = TF.meth.cor, 
             motif = motif,
             save = FALSE,)


