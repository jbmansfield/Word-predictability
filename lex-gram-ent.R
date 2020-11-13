# lex-gram-ent.R

# The compounds dataframe needs to be built first, in extract-compounds.R
library(plyr)
library(dplyr)
library(forcats)
library(ggplot2)
library(entropy)
library(VennDiagram)

# divide compounds dataset into free-lex and bound-lex portions
compounds.freelex <- subset(compounds, Lex.bound==FALSE)
compounds.boundlex <- subset(compounds, Lex.bound==TRUE)
compounds.freelex <- droplevels(compounds.freelex)
compounds.boundlex <- droplevels(compounds.boundlex)

### Now select which type to analyse
# (freelex is the main study for German, boundlex for Murrinhpatha)
cmps <- compounds.freelex
cmps <- compounds.boundlex

# with the restricted dataset calculate independent weighted average entropy of grams, free and bound
grams <- ddply(cmps, .(Gram, Gram.bound), summarise, Gram.freq = sum(Freq))
ent.grams.free <- entropy.ChaoShen(subset(grams, Gram.bound==FALSE)$Gram.freq, unit = "log2")
ent.grams.bound <- entropy.ChaoShen(subset(grams, Gram.bound==TRUE)$Gram.freq, unit = "log2")

# calculate independent weighted average entropy of lexes, in free and bound complexes
lexes <- ddply(cmps, .(Lex, Gram.bound), summarise, Lex.freq = sum(Freq))
ent.lexes.freeGram <- entropy.ChaoShen(subset(lexes, Gram.bound==FALSE)$Lex.freq, unit = "log2")
ent.lexes.boundGram <- entropy.ChaoShen(subset(lexes, Gram.bound==TRUE)$Lex.freq, unit = "log2")

# calculate conditional entropy of grams for each lex
lex.ent <- ddply(cmps, .(Lex, Gram.bound, Gram.bound.lab), summarise, Lex.freq = sum(Freq), Type.freq = length(Lex), Ent = entropy.empirical(Freq, unit = "log2"), Ent.CS = entropy.ChaoShen(Freq, unit = "log2"))

# proportional predictability
lex.ent$MI.CS.prop <- rep(NA, nrow(lex.ent))
lex.ent$MI.CS.prop[lex.ent$Gram.bound==FALSE] <- (ent.grams.free - lex.ent$Ent.CS[lex.ent$Gram.bound==FALSE]) / ent.grams.free 
lex.ent$MI.CS.prop[lex.ent$Gram.bound==TRUE] <- (ent.grams.bound - lex.ent$Ent.CS[lex.ent$Gram.bound==TRUE]) / ent.grams.bound
# a handful come out with negative MI.CS.prop values, because CS may overestimate the conditional entropy at low counts to be more than the independent entropy. For these cases, we switch to using empirical entropy instead
lex.ent$MI.CS.prop[lex.ent$MI.CS.prop<0 & lex.ent$Gram.bound==FALSE] <- (ent.grams.free - lex.ent$Ent[lex.ent$MI.CS.prop<0 & lex.ent$Gram.bound==FALSE]) / ent.grams.free 
lex.ent$MI.CS.prop[lex.ent$MI.CS.prop<0 & lex.ent$Gram.bound==TRUE] <- (ent.grams.bound - lex.ent$Ent[lex.ent$MI.CS.prop<0 & lex.ent$Gram.bound==TRUE]) / ent.grams.bound

# calculate weighted average conditional entropy of grams, i.e. H(G|L), separating free and bound
lex.ent$weighted.ent <- rep(NA, nrow(lex.ent))

lex.ent$weighted.ent[lex.ent$Gram.bound==FALSE] <- lex.ent$Ent.CS[lex.ent$Gram.bound==FALSE] * (lex.ent$Lex.freq[lex.ent$Gram.bound==FALSE] / sum(lex.ent$Lex.freq[lex.ent$Gram.bound==FALSE]))

lex.ent$weighted.ent[lex.ent$Gram.bound==TRUE] <- lex.ent$Ent.CS[lex.ent$Gram.bound==TRUE] * (lex.ent$Lex.freq[lex.ent$Gram.bound==TRUE] / sum(lex.ent$Lex.freq[lex.ent$Gram.bound==TRUE]))


# histogram by lex 
# change threshold count depending on whether you are doing main data or minority data

ggplot(subset(lex.ent, Lex.freq > 9), aes(MI.CS.prop)) + geom_histogram(binwidth = 0.1, boundary = 0) + theme_bw() + theme(text = element_text(size = 24)) + facet_wrap(~Gram.bound.lab) + scale_x_continuous(limits = c(0,1)) + ylab("Stem types") + xlab("Internal predictability")
# See also relabelled code below

# counting how many lexes fit into groups (change bound and greater/lesser values)
nrow(subset(lex.ent, (Lex.freq > 9 & Gram.bound==TRUE & MI.CS.prop > 0.9)))

# " The Mann-Whitney test ranked all the values from low to high, and then compared the means ranks."
# aka "Mann-Whitney-Wilcoxon"
# https://www.graphpad.com/support/faq/the-mann-whitney-test-doesnt-really-compare-medians/
# https://data.library.virginia.edu/the-wilcoxon-rank-sum-test/
wilcox.test(MI.CS.prop ~ Gram.bound, data = subset(lex.ent, Lex.freq > 9),  paired=FALSE)
library(rstatix)
subset(lex.ent, Lex.freq > 9) %>% wilcox_effsize(MI.CS.prop ~ Gram.bound)

# Deprecated: Kruskal-Wallis test of specific internal predictability for free versus bound lexes
#kruskal.test(MI.CS.prop ~ Gram.bound, data = subset(lex.ent, Lex.freq > 9))

#CE = weighted average conditional entropy of grams, i.e. H(G|L)
lex.gram.CE.freeGram <- sum(lex.ent$weighted.ent[lex.ent$Gram.bound==FALSE])
#MI = mutual information, i.e. H(G) - H(G|L)
lex.gram.MI.freeGram <- ent.grams.free - lex.gram.CE.freeGram
lex.gram.MI.prop.freeGram <- lex.gram.MI.freeGram / ent.grams.free
lex.gram.MI.prop.lex.freeGram <- lex.gram.MI.freeGram / ent.lexes.freeGram

# Venn diagram
grid.newpage()
draw.pairwise.venn(area1 = ent.grams.free, area2 = ent.lexes.freeGram, cross.area = lex.gram.MI.freeGram, category = c(paste0("H(G)=",ent.grams.free), paste0("H(L)=",ent.lexes.freeGram) ), inverted = TRUE, cat.pos=0)

# and for bound complexes
lex.gram.CE.boundGram <- sum(lex.ent$weighted.ent[lex.ent$Gram.bound==TRUE])
lex.gram.MI.boundGram <- ent.grams.bound - lex.gram.CE.boundGram
lex.gram.MI.prop.boundGram <- lex.gram.MI.boundGram / ent.grams.bound
lex.gram.MI.prop.lex.boundGram <- lex.gram.MI.boundGram / ent.lexes.boundGram

grid.newpage()
draw.pairwise.venn(area1 = ent.grams.bound, area2 = ent.lexes.boundGram, cross.area = lex.gram.MI.boundGram, category = c(paste0("H(G)=",ent.grams.bound), paste0("H(L)=",ent.lexes.boundGram) ), inverted = TRUE, cat.pos=0)

# Relablling required for German bound stems graph
lex.ent$Gram.bound.lab <- as.character(lex.ent$Gram.bound.lab)
lex.ent$Gram.bound.lab[lex.ent$Gram.bound.lab=="Phrase"] <- "Recursive word"
lex.ent$Gram.bound.lab[lex.ent$Gram.bound.lab=="Word"] <- "Minimal word"
lex.ent$Gram.bound.lab <- as.factor(lex.ent$Gram.bound.lab)
lex.ent$Gram.bound.lab <- factor(lex.ent$Gram.bound.lab, levels = c("Recursive word", "Minimal word"))

# THIS IS NO LONGER REPORTED IN THE STUDY
# Now  divide the whole set into free and bound lexes
#

# labelling again :(
compounds$Lex.bound.lab <- as.factor(compounds$Lex.bound)
# german
levels(compounds$Lex.bound.lab)[levels(compounds$Lex.bound.lab)=="FALSE"] <- "Free stem"
levels(compounds$Lex.bound.lab)[levels(compounds$Lex.bound.lab)=="TRUE"] <- "Bound stem"

lexes <- ddply(compounds, .(Lex, Lex.bound), summarise, Lex.freq = sum(Freq))
ent.lexes.free <- entropy.ChaoShen(subset(lexes, Lex.bound==FALSE)$Lex.freq, unit = "log2")
ent.lexes.bound <- entropy.ChaoShen(subset(lexes, Lex.bound==TRUE)$Lex.freq, unit = "log2")

# calculate independent weighted average entropy of grams, in free and bound lex complexes
grams <- ddply(compounds, .(Gram, Lex.bound), summarise, Gram.freq = sum(Freq))
ent.grams.freeLex <- entropy.ChaoShen(subset(grams, Lex.bound==FALSE)$Gram.freq, unit = "log2")
ent.grams.boundLex <- entropy.ChaoShen(subset(grams, Lex.bound==TRUE)$Gram.freq, unit = "log2")

lex.ent <- ddply(compounds, .(Lex, Lex.bound, Lex.bound.lab), summarise, Lex.freq = sum(Freq), Type.freq = length(Lex), Ent = entropy.empirical(Freq, unit = "log2"), Ent.CS = entropy.ChaoShen(Freq, unit = "log2"))

# proportional predictability: proportion of gram entropy that is consumed by knowing the lex
lex.ent$MI.CS.prop <- rep(NA, nrow(lex.ent))
lex.ent$MI.CS.prop[lex.ent$Lex.bound==FALSE] <- (ent.grams.freeLex - lex.ent$Ent.CS[lex.ent$Lex.bound==FALSE]) / ent.grams.freeLex 
lex.ent$MI.CS.prop[lex.ent$Lex.bound==TRUE] <- (ent.grams.boundLex - lex.ent$Ent.CS[lex.ent$Lex.bound==TRUE]) / ent.grams.boundLex

# calculate weighted average conditional entropy of lexes, i.e. H(L|G) * P(L), separating free and bound lexes
lex.ent$weighted.ent <- rep(NA, nrow(lex.ent))

lex.ent$weighted.ent[lex.ent$Lex.bound==FALSE] <- lex.ent$Ent.CS[lex.ent$Lex.bound==FALSE] * (lex.ent$Lex.freq[lex.ent$Lex.bound==FALSE] / sum(lex.ent$Lex.freq[lex.ent$Lex.bound==FALSE]))

lex.ent$weighted.ent[lex.ent$Lex.bound==TRUE] <- lex.ent$Ent.CS[lex.ent$Lex.bound==TRUE] * (lex.ent$Lex.freq[lex.ent$Lex.bound==TRUE] / sum(lex.ent$Lex.freq[lex.ent$Lex.bound==TRUE]))

# histogram by lex
ggplot(subset(lex.ent, Lex.freq > 4), aes(MI.CS.prop)) + geom_histogram(binwidth = 0.05, boundary = 0) + theme_bw() + theme(text = element_text(size = 24)) + facet_wrap(~Lex.bound, labeller = label_both)  + scale_x_continuous(limits = c(0, 1))

#CE = weighted average conditional entropy of grams, i.e. H(G|L)
lex.gram.CE.freeLex <- sum(lex.ent$weighted.ent[lex.ent$Lex.bound==FALSE])
#MI = mutual information, i.e. H(G) - H(G|L)
lex.gram.MI.freeLex <- ent.grams.free - lex.gram.CE.freeLex
lex.gram.MI.prop.freeLex <- lex.gram.MI.freeLex / ent.grams.freeLex
lex.gram.MI.prop.lex.freeLex <- lex.gram.MI.freeLex / ent.lexes.free

# Venn diagram
grid.newpage()
draw.pairwise.venn(area1 = ent.lexes.free, area2 = ent.grams.freeLex, cross.area = lex.gram.MI.freeLex, category = c("Lex", "Gram"), inverted = TRUE, cat.pos=0) # something wrong here

# and for bound complexes
lex.gram.CE.boundLex <- sum(lex.ent$weighted.ent[lex.ent$Lex.bound==TRUE])
lex.gram.MI.boundLex <- ent.grams.boundLex - lex.gram.CE.boundLex
lex.gram.MI.prop.boundLex <- lex.gram.MI.boundLex / ent.grams.boundLex
lex.gram.MI.prop.lex.boundLex <- lex.gram.MI.boundLex / ent.lexes.bound

grid.newpage()
draw.pairwise.venn(area1 = ent.grams.boundLex, area2 = ent.lexes.bound, cross.area = lex.gram.MI.boundLex, category = c("Lex", "Gram"), inverted = TRUE, cat.pos=0)


