library(limma)

#shared otus

shared_otu_table = read.table("shared_otus.txt", row.names=1,sep="\t", header=TRUE)
shared_otu_table = read.table("ITS_shared_otus.txt", row.names=1,sep="\t", header=TRUE)
ps_bray = phyloseq::distance(ps, method = "bray")
ps_jaccard = phyloseq::distance(ps, method = "jaccard")
ps_unifrac = phyloseq::distance(ps, method = "unifrac")

#transform to presence/absence OTU table
ps.pa <- transform_sample_counts(ps,function(x)1*(x>0))
#transform to relative abundance table
ps.ra <- transform_sample_counts(ps,function(x)x/sum(x))

ps.df <- data.frame(sample_data(ps))
ps.df <- subset(ps.df,select=c(Year, Experiment,Seed_source, Disinfection, Inoculation,Tissue, Farm, Experiment))
View(ps.df)


osu2013 <- subset(ps.df,Experiment==1 & Year==2013)

ofrf2014 <- subset(ps.df,Experiment==1 & Year==2014)
ps_2012.ss <- subset(ps.df, Experiment==2 & Year==2012)

ps_2012.ss1 <- subset(ps_2012.ss, Seed_source==1)
ps_2012.ss2 <- subset(ps_2012.ss, Seed_source==2)



ps_2013.ss <- subset(ps.df, Experiment==2 & Year==2013)
ps_2013.ss3 <- subset(ps_2013.ss, Seed_source==3)
ps_2013.ss4 <- subset(ps_2013.ss, Seed_source==4)


#Shared OTUs between seed sources
sot.2013 <- shared_otu_table[rownames(ps_2012.ss1),rownames(ps_2012.ss2)]

#Shared OTUs between 2013 crown and seed samples and the 2012 seed source #1

sot.2013 <- shared_otu_table[rownames(osu2013),rownames(ps_2012.ss2)]

sot.2013 <-cbind(sot.2013,meanshared=rowMeans(sot.2013),perc_otus_shared=100*rowMeans(sot.2013)/sample_sums(ps)[rownames(sot.2013)])

sot.2013 <- cbind(osu2013,meanshared=sot.2013$meanshared)
fit <-aov(meanshared ~ Seed_source*Tissue,data=sot.2013)
summary(fit)

plot(sot.2013$meanshared ~ sot.2013$Tissue,xlab="Tissue",
          ylab="Shared OTUs", main="Mean Plot\nwith 95% CI")

interaction.plot(sot.2013$Seed_source, sot.2013$Tissue, sot.2013$meanshared, type="b", col=c(1:3),
                 leg.bty="o", leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Seed_source",
                 ylab="Shared OTUs",
                 main="Interaction Plot")

#Shared OTUs between 2013 crown and seed samples and the 2012 seed source #2

sot.2013 <- shared_otu_table[sample_names(ps.2013),sample_names(ps_2012.ss2)]
sot.2013 <-cbind(sot.2013,meanshared=rowMeans(sot.2013),perc.shared=100*rowMeans(sot.2013)/sample_sums(ps.2013.pa))
ps.2013.df <- data.frame(sample_data(ps.2013))
sot.2013 <- cbind(sot.2013,Seed_source=ps.2013.df$Seed_source,Tissue=ps.2013.df$Tissue)
fit <-aov(meanshared ~ Seed_source*Tissue,data=sot.2013)
summary(fit)

attatch(sot.2013)
plotmeans(meanshared ~ Tissue,xlab="Tissue",
          ylab="Shared OTUs", main="Mean Plot\nwith 95% CI")
interaction.plot(Seed_source, Tissue, meanshared, type="b", col=c(1:3),
                 leg.bty="o", leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Seed_source",
                 ylab="Shared OTUs",
                 main="Interaction Plot")

#2014

sot.2014 <- shared_otu_table[sample_names(ps.2014),sample_names(ps_2013.ss4)]
sot.2014 <-cbind(sot.2014,meanshared=rowMeans(sot.2014),perc.shared=100*rowMeans(sot.2014)/sample_sums(ps.2014.pa))
ps.2014.df <- data.frame(sample_data(ps.2014))
sot.2014 <- cbind(sot.2014,Farm=ps.2014.df$Farm,Seed_source=ps.2014.df$Seed_source,Tissue=ps.2014.df$Tissue)
fit <-aov(meanshared ~ Seed_source+Farm+Tissue+Seed_source:Farm+Seed_source:Tissue+Farm:Tissue,data=sot.2014)
fit <-aov(meanshared ~ Seed_source,data=sot.2014)
summary(fit)

attach(sot.2014)
plotmeans(meanshared ~ Seed_source,xlab="Seed_source",
          ylab="Shared OTUs", main="Mean Plot\nwith 95% CI")
interaction.plot(Seed_source, Farm, meanshared, type="b", col=c(1:3),
                 leg.bty="o", leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Seed_source",
                 ylab="Shared OTUs",
                 main="Interaction Plot")

