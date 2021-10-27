library(Biostrings)
library(reshape2)
library(ggplot2)
library(ggpubr)
catsSims<- readRDS("C:/DANE/R/pixy/Cat simulaitons.RDS")
barbusSims <- readRDS("C:/DANE/R/pixy/Barbus simulaitons.RDS")
PIs <- function(data){
  dnas <- as.data.frame(as.character(data))
  freqs <- apply(dnas,1,table,exclude = "n")
  nDiff <- sapply(freqs,function(z){ 
    tabs <- outer(z,z)
    sum(tabs[lower.tri(tabs)])})
  nComp <- sapply(freqs,sum)
  nComp <- (nComp*(nComp-1))/2
  WeightedPi <- mean(nDiff/nComp,na.rm = TRUE)
  PEGASres <- pegas::nuc.div(data)
  PIXYres <- sum(nDiff)/sum(nComp)
  return(data.frame(Ns = 0, weighted = WeightedPi, 
                    pixy = PIXYres))
}
ftPIs <- function(data){
  fasta <- t(data.frame(as.character(data)))
  gathered <- unlist(fasta) 
  Ns <- sample(1:(0.9*length(gathered)),1)
  gathered[sample(1:length(gathered),Ns)] <- "n"
  dnas <- matrix(gathered,nrow = dim(fasta)[1],ncol = dim(fasta)[2],byrow = FALSE,
                 dimnames = list(names(data),c()))
  freqs <- apply(dnas,2,table,exclude = "n")
  emptySites <- sapply(unlist(freqs),dim)
  if (!is.list(emptySites)) dnas <- dnas[,emptySites!=0]
  nDiff <- sapply(freqs,function(z){ 
    tabs <- outer(z,z)
    sum(tabs[lower.tri(tabs)])})
  nComp <- sapply(freqs,sum)
  nComp <- (nComp*(nComp-1))/2
  pixyres <- sum(nDiff)/sum(nComp)
  WeightedPi <- mean(nDiff/nComp,na.rm = TRUE)
  return(data.frame(Ns = Ns/length(gathered), 
                    weighted = WeightedPi, 
                    pixy = pixyres))
}
# cl <- makeCluster(detectCores())
# clusterExport(cl,c("ftPIs"))
# clusterEvalQ(cl, {library(ape)})
# 
barbs <- ape::read.FASTA("./fastas/barb.fasta",type = "DNA")
# clusterExport(cl,"barbs")
# barbusSims <- parSapply(cl,1:100000,function(z){
#   return(ftPIs(barbs))})
# barbusSims <- apply(barbusSims,1,unlist)

barbNeiPi <- PIs(barbs)[1,3]

cats <- ape::read.FASTA("./fastas/cat.fasta",type = "DNA")
# clusterExport(cl,"cats")
# catsSims <- parSapply(cl,1:100000,function(z){
#   return(ftPIs(cats))})

# stopCluster(cl)

# catsSims <- apply(catsSims,1,unlist)
catsNeiPi <- PIs(cats)[1,3]


colours <- (c(weighted = "grey65", pixy = "grey30"))

##### Analysis of simulations of feline mtDNA CR ####
nClass <- cut(catsSims[,1],breaks = seq(0,.9,by = 0.1))
levels(nClass) <- sub("[(]","",levels(nClass))
levels(nClass) <- sub("[]]","",levels(nClass))
levels(nClass) <- sub("[,]","-",levels(nClass))
catmelt <- melt(catsSims[,2:3]) 
catmelt <- cbind(Ns = sapply(catmelt$Var1, function (z) nClass[z]),catmelt)
cats_ggList <- list()
for (NsScale in levels(catmelt$Ns)){
  cats_ggList[[NsScale]] <- ggplot(subset(catmelt,Ns == NsScale)) +
    geom_density(aes(x=value, fill=as.factor(Var2),colour = as.factor(Var2)), 
                 alpha = 0.2,size = 0.7) +
    geom_vline(xintercept = catsNeiPi,size = 0.4,linetype = 3)+
    theme_classic2() +
    xlim(0,0.01)+
    scale_fill_manual(name = "", values=colours,
                      labels = c(expression(italic(pi[W])),
                                 expression(italic(pi[pixy])))) +
    scale_colour_manual(name = "", values=colours,
                        labels = c(expression(italic(pi[W])),
                                   expression(italic(pi[pixy])))) +
    xlab(expression(italic(pi)))+
    labs(legend_title = "")+
    
    if (NsScale==levels(catmelt$Ns)[1]) {
      theme(legend.position = "top")}
  else
    theme(legend.position = "none")
}

pdf("cats_densities_3_plots.pdf", width = 4, height = 6,onefile = FALSE)
ggarrange(plotlist =cats_ggList[c(1,5,9)],labels = "AUTO",
          ncol = 1,common.legend = TRUE,align = "hv")
dev.off()


catdevs <- data.frame(Ns = nClass,abs((catsSims[,2:3] - catsNeiPi)/catsNeiPi))
catDevMelt <- melt(t(catdevs[,2:3]),value.name = "RAE")
catDevMelt <- cbind(Ns = sapply(catDevMelt$Var2,function(z)catdevs[z,1]),catDevMelt)
quants <- data.frame(piW = sapply(levels(catdevs$Ns),function(lev)
  round(quantile(catdevs$weighted[catdevs$Ns==lev],c(.95),na.rm = TRUE),4)))
quants <- cbind(quants,pixy = sapply(levels(catdevs$Ns),function(lev)
  round(quantile(catdevs$pixy[catdevs$Ns==lev],c(.95),na.rm = TRUE),4)))
write.table(quants,file = "catsQuants.xls",quote = FALSE,col.names = NA,sep = "\t")

cat_relErr <-  ggplot(catDevMelt,aes(x = Ns,y=RAE,fill = Var1) )+ geom_boxplot() +
  scale_fill_manual(name = "", values=colours,
                    labels = c(expression(italic(pi[W])),
                               expression(italic(pi[pixy])))) +
  geom_hline(yintercept = 0)+ 
  xlab("proportion of missing data") + 
  ylab("absolute relative error") +
  theme_classic2()+
  theme(legend.position = "none")


pdf("cats_errors_2_plots.pdf", width = 4, height = 6,onefile = FALSE)
cat_relErr
dev.off()

####### Analysis of Barbel DNA simulations
nClass <- cut(barbusSims[,1],breaks = seq(0,.9,by = 0.1))
levels(nClass) <- sub("[(]","",levels(nClass))
levels(nClass) <- sub("[]]","",levels(nClass))
levels(nClass) <- sub("[,]","-",levels(nClass))

barbmelt <- melt(barbusSims[,2:3]) 
barbmelt <- cbind(Ns = sapply(barbmelt$Var1, function (z) nClass[z]),barbmelt)
barb_ggList <- list()
for (NsScale in levels(barbmelt$Ns)){
  barb_ggList[[NsScale]] <- ggplot(subset(barbmelt,Ns == NsScale)) +
    geom_density(aes(x=value, fill=as.factor(Var2),colour = as.factor(Var2)), 
                 alpha = 0.2, size = 0.7) +
    geom_vline(xintercept = barbNeiPi,size = 0.4,linetype = 3)+
    theme_classic2() +
    xlim(0,0.001)+
    scale_fill_manual(name = "", values=colours,
                      labels = c(expression(italic(pi[W])),
                                 expression(italic(pi[pixy])))) +
    scale_colour_manual(name = "", values=colours,
                        labels = c(expression(italic(pi[W])),
                                   expression(italic(pi[pixy])))) +
    xlab(expression(italic(pi)))+
    labs(legend_title = "")+
    
    if (NsScale==levels(barbmelt$Ns)[1]) {
      theme(legend.position = "top")}
  else
    theme(legend.position = "none")
}

pdf("barb_densities_3_plots.pdf", width = 8, height = 15,onefile = FALSE)
ggarrange(plotlist =barb_ggList[c(1,5,9)],labels = "AUTO",
          ncol = 1,common.legend = TRUE,align = "hv")
dev.off()

barbdevs <- data.frame(Ns = nClass,abs((barbusSims[,2:3] - barbNeiPi)/barbNeiPi))
barbDevMelt <- melt(t(barbdevs[,2:3]),value.name = "RAE")
barbDevMelt <- cbind(Ns = sapply(barbDevMelt$Var2,function(z)barbdevs[z,1]),barbDevMelt)
quants <- data.frame(piW = sapply(levels(barbdevs$Ns),function(lev)
  round(quantile(barbdevs$weighted[barbdevs$Ns==lev],c(.95),na.rm = TRUE),4)))
quants <- cbind(quants,pixy = sapply(levels(barbdevs$Ns),function(lev)
  round(quantile(barbdevs$pixy[barbdevs$Ns==lev],c(.95),na.rm = TRUE),4)))
write.table(quants,file = "barbQuants.xls",quote = FALSE,col.names = NA,sep = "\t")

barb_relErr <- ggplot(barbDevMelt,aes(x = Ns,y=RAE,fill = Var1) )+ geom_boxplot() +
  scale_fill_manual(name = "", values=colours,
                    labels = c(expression(italic(pi[W])),
                               expression(italic(pi[pixy])))) +
  geom_hline(yintercept = 0)+ 
  xlab("proportion of missing data") + 
  ylab("absolute relative error") +
  theme_classic2()+
  theme(legend.position = "none")


pdf("barb_errors_2_plots.pdf", width = 8, height = 12,onefile = FALSE)
barb_relErr
dev.off()

# save.image("sims.RData")


####### common plots ############
pdf("errors_2_plots.pdf", width = 6, height = 6,onefile = FALSE)
ggarrange(cat_relErr,barb_relErr, 
          ncol = 1,nrow = 2,common.legend = TRUE,align = "hv",labels = c("A","B"))
dev.off()
pdf("densities_3_plots.pdf", width = 8, height = 8,onefile = FALSE)
ggarrange(plotlist =c(cats_ggList[1],barb_ggList[1],
                      cats_ggList[5],barb_ggList[5],
                      cats_ggList[9],barb_ggList[9]),labels = c("A","D","B","E","C","F"),
          nrow = 3,ncol = 2,common.legend = TRUE,align = "hv")
dev.off()

sum(catdevs$weighted==0)
sum(catdevs$pixy==0)
sum(barbdevs$weighted==0)
sum(barbdevs$pixy==0)


# desc_statby(barbmelt,measure.var= "value",grps = c("Var2","Ns"))
# desc_statby(catDevMelt,measure.var= "RAE",grps = c("Var1","Ns"))
# head(barbDevMelt)

