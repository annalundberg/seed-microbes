#### Illustrate inoculation effect, OSU 2014
OFRF2014seeds.osu_bray <- phyloseq::distance(OFRF2014seeds.osu, method='bray')
View(sample_data(OFRF2014seeds.osu))

ord = ordinate(OFRF2014seeds.osu, method = "NMDS", distance = OFRF2014seeds.osu_bray)
plot_ordination(OFRF2014seeds.osu, ord, color = "Inoculation") + ggtitle("NMDS: Bray-Curtis")
ord = ordinate(OFRF2014seeds.osu, method = "NMDS", distance = OFRF2014seeds.osu_unifrac)
plot_ordination(OFRF2014seeds.osu, ord, color = "Inoculation") + ggtitle("NMDS: Unifrac Distance")



#### 

#### Illustrate inoculation effect, OSU 2014
F.Full.ps <- subset_samples(ps,Experiment ==)



bray <- phyloseq::distance(ps_ITS, method='bray')


ord = ordinate(ps_ITS, method = "NMDS", distance = bray)
plot_ordination(ps, ord, color = "Tissue") + ggtitle("NMDS: Bray-Curtis")


jaccard <- phyloseq::distance(ps_ITS, method='jaccard')


ord = ordinate(ps_ITS, method = "NMDS", distance = jaccard)
plot_ordination(ps, ord, color = "Tissue") + ggtitle("NMDS: Jaccard")