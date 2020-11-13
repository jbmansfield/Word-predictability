# extract-compounds.R
library(plyr)
library(dplyr)
library(forcats)
library(ggplot2)
library(entropy)

# # # # # # # # #
# SELECT language

# choose Murrinhpatha input
mp <- ct
mp <- ct.inclSimples
mp <- ct.samp.gram
mp <- ct.samp.lex
mp$m.classifier.separable <- mp$m.classifier.bound # dummy variable to make same structure as DE
compounds <- as.data.frame(table(mp$w.coverb.stem, mp$m.classifier.base, mp$m.classifier.bound, mp$w.coverb.bound, mp$m.classifier.separable))

# do reformat line below....

# choose German input
de <- gt
de <- gt.inclSimples
de <- gt.samp
de <- gt.samp.gram
de <- gt.samp.lex

compounds <- as.data.frame(table(de$Stem, de$Preverb, de$Preverb.bound, de$Stem.bound, de$Preverb.separable))

# # # # # 
# Reformat (both langs)
colnames(compounds) <- c("Lex", "Gram", "Gram.bound", "Lex.bound", "Gram.separable", "Freq")
compounds <- compounds[compounds$Freq > 0,]

# hacky labelling :(
compounds$Gram.bound.lab <- as.factor(compounds$Gram.bound)
# german
levels(compounds$Gram.bound.lab)[levels(compounds$Gram.bound.lab)=="FALSE"] <- "Phrase"
levels(compounds$Gram.bound.lab)[levels(compounds$Gram.bound.lab)=="TRUE"] <- "Word"
# murrinhpatha
levels(compounds$Gram.bound.lab)[levels(compounds$Gram.bound.lab)=="FALSE"] <- "Free finite stem"
levels(compounds$Gram.bound.lab)[levels(compounds$Gram.bound.lab)=="TRUE"] <- "Classifier stem"

# get graph probs and gram entropies
grams <- ddply(compounds, .(Gram, Gram.bound, Gram.bound.lab), summarise, Gram.freq = sum(Freq))

# graph the grams
free.gram.count <- sum(subset(grams, Gram.bound==FALSE)$Gram.freq)
bound.gram.count <- sum(subset(grams, Gram.bound==TRUE)$Gram.freq)
grams$prop <- rep(0, nrow(grams))
grams$prop[grams$Gram.bound==FALSE] <- grams$Gram.freq[grams$Gram.bound==FALSE] / free.gram.count
grams$prop[grams$Gram.bound==TRUE] <- grams$Gram.freq[grams$Gram.bound==TRUE] / bound.gram.count

grams$Gram <- reorder(grams$Gram, -grams$Gram.freq)

ggplot(grams, aes(Gram)) + geom_bar(aes(weight=prop), width=0.8) + facet_grid(.~grams$Gram.bound.lab, scales="free_x", space = "free_x") + theme_bw() + theme(text = element_text(size = 24), axis.text=element_text(size=11, angle = 45), axis.title.x=element_blank()) + ylab("Probability")

# Entropy of the gram distributions
ent.grams.free <- entropy.ChaoShen(subset(grams, Gram.bound==FALSE)$Gram.freq, unit = "log2")
ent.grams.bound <- entropy.ChaoShen(subset(grams, Gram.bound==TRUE)$Gram.freq, unit = "log2")


# SimLang input
compounds <- as.data.frame(table(sim.tokens$Lex, sim.tokens$Gram, sim.tokens$Gram.bound, sim.tokens$Lex.bound))
