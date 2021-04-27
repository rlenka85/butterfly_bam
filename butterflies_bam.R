#Required packages
install.packages("vegan")
library(vegan)
install.packages("permute")
install.packages("lattice")
library(lattice)
library(permute)

#upload the file
butterfly <- read.csv('butterflybu.csv', row.names = 1)
butterfly
View(butterfly)
#Total number of species in each site
speciesrichness <- apply(butterfly[,-1]>0,1,sum)
speciesrichness
View(speciesrichness)
write.table(speciesrichness, file = "speciesrichness.csv", sep = ",")


#Diveristy indices

#1. Mechinik's index
n <-apply(butterfly[,-1]>0,1,sum)
N <- apply(butterfly[,-1],1,sum)
mechinik <- n/sqrt(N)
mechinik
write.table(mechinik, file = "mechinik.csv", sep = ",")
View(mechinik)

#2. Margalef's Index
n<-apply(butterfly[,-1]>0,1,sum)
N <- apply(butterfly[,-1],1,sum)
Meg <- (n-1)/log(N)
Meg
write.table(Meg, file = "margalef.csv", sep = ",")
View(Meg)

#3. Shannon-Wiener Index (H')
shannon <- diversity(butterfly[-1], index="shannon")
shannon
write.table(shannon, file = "shannon.csv", sep = ",")
View(shannon)


#5. Simpson's Index (λ)
simpson <- diversity(butterfly[-1], index="simpson")
simpson
write.table(simpson, file = "simpson.csv", sep = ",")
View(simpson)

hist(shannon)
#alpha, beta and gamma diversity
#alpha is nothing but shannon and simpson
#but beta diveristy indices are of 24 types for presence-absence data. 
#For more details, refer to Koleff, P., Gaston, K.J. and Lennon, J.J. (2003) 
#Measuring beta diversity for presence-absence data. Journal of Animal Ecology. 72: 367-382.

#Species abundance and density (Number of individuals)
speciesabundance <- apply(butterfly[,-1],1,sum)
speciesabundance
write.table(speciesabundance, file = "speciesabundance.csv", sep = ",")
View(speciesabundance)

#Rarefaction
rarefaction <- rarefy(butterfly[-1], sample=10, MARGIN=1)
rarefaction
write.table(rarefaction, file = "rarefraction.csv", sep = ",")
View(rarefaction)

#Evenness
#1. Pilou's Evenness (J)
library(vegan)
S <- apply(butterfly[,-1]>0,1,sum)
pilou <- diversity(butterfly[-1], index="simpson")/log(S)
pilou
write.table(pilou, file = "pilou.csv", sep = ",")
View(pilou)

#2. Hill's ratios (Ea:b)
library(vegan)
S <- apply(butterfly[,-1]>0,1,sum)
hills <- exp(diversity(butterfly[-1], index="simpson"))/S
write.table(hills, file = "hills.csv", sep = ",")
View(hills)



#species accumulation curve
sac <- specaccum(butterfly) #sac-species accumulation curve
plot(sac, ci.type="polygon", ci.col="yellow")

#Sorensen index of dissimilarity
beta <- vegdist(butterfly, binary=TRUE)
mean(beta)

#For beta diversity, no specific functions are needed, but
#this index can be easily found with the help of vegan function specnumber
ncol(butterfly)/mean(specnumber(birds)) - 1

bray = vegdist(butterfly, "bray") 
gower = vegdist(butterfly, "gower")
hist(bray, xlim = range(0.0,1.0))
hist(gower, xlim = range(0.0,1.0))

#Species Abundance
spAbund <- rowSums(butterfly)  #gives the number of individuals found in each plot
spAbund # view observations per plot (number of individuals per plot)

#Rarefaction
raremin <- min(rowSums(butterfly))  #rarefaction uses the smallest number of observations per sample to extrapolate the expected number if all other samples only had that number of observations
raremin # view smallest # of obs (site 17)
sRare <- rarefy(butterfly, raremin) # now use function rarefy
sRare #gives an "expected"rarefied" number of species (not obs) if only 15 individuals were present
rarecurve(butterfly, col = "blue")



#MDS, NMDS
# Our community-by-species matrix
NMDS=metaMDS(butterfly,k=2)
plot(NMDS)
ordiplot(butterfly,type="n") #Ordination plot function especially for congested plots
orditorp(butterfly,display="species",col="red",air=0.01) #The function adds text or points to ordination plots
orditorp(butterfly,display="sites",cex=1.25,air=0.01)


#The phyloseq package (McMurdie and Holmes (2013)) can be used to quickly 
#plot a variety of alpha diversity indexes per sample using the plot_richness function. 
install.packages("phyloseq") #NOT FOR THIS VERSION, TRY SOMEWHERE ELSE
library(phyloseq)

