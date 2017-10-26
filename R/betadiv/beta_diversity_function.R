
Field_model

OFRF2014_model <- Farm+Tissue+Disinfection+Inoculation+Farm*Tissue+Farm*Disinfection+Farm*Inoculation+Tissue*Disinfection
+Tissue*Inoculation+Disinfection*Inoculation

OFRF2014


# function for taking a phyloseq object and model and running it through a range of tests
field_model <- function(ps, design) {
  
ps_bray = phyloseq::distance(ps, method = "bray")
ps_unifrac = phyloseq::distance(OFRF2014, method = "unifrac")
OFRF2014df = data.frame(sample_data(OFRF2014))
adonis(OFRF2014_bray ~ Farm+Tissue+Disinfection+Inoculation
       +Farm*Tissue+Farm*Disinfection+Farm*Inoculation+Tissue*Disinfection
       +Tissue*Inoculation+Disinfection*Inoculation,data=OFRF2014df)
adonis(OFRF2014_unifrac ~ Farm+Tissue+Disinfection+Inoculation
       +Farm*Tissue+Farm*Disinfection+Farm*Inoculation+Tissue*Disinfection
       +Tissue*Inoculation+Disinfection*Inoculation,data=OFRF2014df)


plot_ordination(Field, ord, color = "Tissue") + ggtitle("NMDS: Bray-Curtis")

plot_ordination(Field, ord, color = "Seed_source") + ggtitle("NMDS: Bray-Curtis")

}


