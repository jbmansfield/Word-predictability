# exploring possible factors to explain some more variation

# trying to understand among DE phrasal verbs, why some fall into a minority group with high internal predictability.
# could be to do with distribution of verb stem among phrasal, complex word, and simple verb types???
# hypothesis: it will be more internally predictable within the phrasal verb construction type, when it has low freq as a simple verb. But freq in complex word subtype will not make much difference

library(tidyr)

# This is used for German
compounds.inclSimples <- compounds
compounds.inclSimples$Gram.bound <- as.character(compounds.inclSimples$Gram.bound)
compounds.inclSimples$Gram.bound[compounds.inclSimples$Gram=="NONE"] <- "NONE"
compounds.inclSimples$Gram.bound <- as.factor(compounds.inclSimples$Gram.bound)

cxtypes <- ddply(compounds.inclSimples, .(Lex, Gram.bound), summarise, Freq = sum(Freq))
cxtypes.wide <- spread(cxtypes, Gram.bound, Freq, fill = 0)
names(cxtypes.wide )[names(cxtypes.wide ) == 'FALSE'] <- 'Freq.Freegram'
names(cxtypes.wide )[names(cxtypes.wide ) == 'TRUE'] <- 'Freq.Boundgram'
names(cxtypes.wide )[names(cxtypes.wide ) == 'NONE'] <- 'Freq.Nogram'

lex.ent.cxtypes <- inner_join(lex.ent, cxtypes.wide, by="Lex")

# explore relationship of internal predictability in phrasal verbs with the frequency as a simple verb
ggplot(subset(lex.ent.cxtypes, Gram.bound=="FALSE"), aes(MI.CS.prop, Freq.Nogram)) + geom_jitter(width = 0.05) + scale_y_sqrt(breaks =c(100, 500, 1000, 2000, 5000))

lex.ent.cxtypes$rareAsSimple <- rep(FALSE, nrow(lex.ent.cxtypes))
lex.ent.cxtypes$rareAsSimple[lex.ent.cxtypes$Freq.Nogram < 20] <- TRUE
# hacky convert to factor for ordering
lex.ent.cxtypes$rareAsSimple <- factor(lex.ent.cxtypes$rareAsSimple, levels = c("TRUE", "FALSE"))
# histo

ggplot(subset(lex.ent.cxtypes, (Lex.freq > 9)), aes(MI.CS.prop, fill=rareAsSimple)) + geom_histogram(binwidth = 0.1, boundary = 0) + theme_bw() + theme(text = element_text(size = 24)) + facet_wrap(.~Gram.bound.lab) + scale_x_continuous(limits = c(0,1))   + scale_fill_manual(name="Simple freq.", values=c("grey60", "grey30"), labels=c("Low", "High"), guide = guide_legend(reverse=TRUE) ) + ylab("Stem types") + xlab("Internal predictability")

# Examples
rareAsSimple <- subset(lex.ent.cxtypes, (Lex.freq > 9 & rareAsSimple==TRUE & Gram.bound=="FALSE"))

# # #
# Combinatoric predictability of an element within a construction type is highly correlated with its frequency outside that construction type
# This encompasses "bound" elements, that just have zero freq outside the construction type
# # #

# This is for Murrinhpatha
# Want to know the weighted average simpleFreq of the finite stems with which coverbs combine

# get gram freqs
gramfreqs.simple <- subset(cmps, Lex=="none")
gramfreqs.simple <- subset(gramfreqs.simple, select=c(Gram, Freq))
names(gramfreqs.simple)[names(gramfreqs.simple) == 'Freq'] <- 'Gram.freq'
# merge gram freqs in to compounds
compounds.gramFreqs <- left_join(cmps, gramfreqs.simple, by="Gram", fill=0)
compounds.gramFreqs[is.na(compounds.gramFreqs)] <- 0

# calculate conditional entropy of grams for each lex
lex.ent <- ddply(compounds.gramFreqs, .(Lex, Gram.bound), summarise, Lex.freq = sum(Freq), Type.freq = length(Lex), Ent = entropy.empirical(Freq, unit = "log2"), Ent.CS = entropy.ChaoShen(Freq, unit = "log2"), GramFreq.av = mean(Gram.freq))

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

# exploratory scatter plot
ggplot(subset(lex.ent, Gram.bound=="FALSE"), aes(MI.CS.prop, GramFreq.av)) + geom_jitter(width = 0.1) 

# histogram by lex 
lex.ent$highGramFreq <- rep(FALSE, nrow(lex.ent))
lex.ent$highGramFreq[lex.ent$GramFreq.av > 1500] <- TRUE
# hacky convert to factor for ordering
lex.ent$highGramFreq <- factor(lex.ent$highGramFreq, levels = c("TRUE", "FALSE"))

ggplot(subset(lex.ent, Lex.freq > 4), aes(MI.CS.prop, fill=highGramFreq)) + geom_histogram(binwidth = 0.1, boundary = 0) + theme_bw() + theme(text = element_text(size = 24)) + facet_wrap(.~Gram.bound, labeller = label_both) + scale_x_continuous(limits = c(0,1))  + scale_fill_manual( values=c("grey60", "grey30") ) 


