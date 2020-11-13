# sampling corpus sizes
library(directlabels)
library(ggplot2)
library(entropy)
library(forcats)
library(plyr)
library(dplyr)

# Sample sizes for individual stem conditional entropy

ent.samples <- data.frame(Lex=character(), 
                          Sample.N=integer(),
                          Ent=numeric(),
                          Method=character(),
                          stringsAsFactors=FALSE) 

lexes.freegram <- c("arbeiten", "streichen", "werfen")
lexes.boundgram <- c("kaufen", "kommen", "reichen")

# sample sizes from 2^3 up to 2^12
# Need to change values each time for "lexes.freegram/boundgram" and "Preverb.bound==TRUE/FALSE
#for (e in 2:12) {
#  n <- 2^e
  
for (n in seq(5, 300, 5)) {
  for (lex in lexes.freegram) {
    lex.tokens <- subset(de, (Stem==lex & Preverb.bound==FALSE))
    samples <- split(lex.tokens, (seq(nrow(lex.tokens))-1) %/% n) # creates a list
    l <- length(samples)
    samples[l] <- NULL # remove the last one, because this one is the residue after all samples of size n have been taken
    for (sm in samples) {
      gram.freqs <- as.data.frame(table(sm$Preverb))
      colnames(gram.freqs) <- c("Gram", "Freq")
      gram.freqs <- gram.freqs[gram.freqs$Freq > 0,]
      Ent.emp <- entropy.empirical(gram.freqs$Freq, unit = "log2")
      Ent.CS <- entropy.ChaoShen(gram.freqs$Freq, unit = "log2")
      #write rows
      ent.samples[nrow(ent.samples)+1,] = list(lex, n, Ent.emp, "Empirical")
      ent.samples[nrow(ent.samples)+1,] = list(lex, n, Ent.CS, "Chao-Shen")
    }
  }
}
ent.samples$Method <- as.factor(ent.samples$Method)

# Plotting
# Change Y axis label: Particle, Prefix; Only add geom_smooth for bound type; and scale breaks
ggplot(ent.samples, aes(Sample.N, Ent)) + geom_jitter(aes(shape=Method), size=2, width=0, height=0.1) + scale_shape_manual(values=c(16, 4)) + facet_wrap(.~Lex) + scale_x_log10(breaks=c(5, 10, 20, 50, 100, 300)) + theme_bw() + theme(text = element_text(size = 24))  + xlab("Sample size") + ylab("Particle entropy") 
+ geom_smooth(aes(group=Method, linetype=Method), method="lm", formula = y ~ x + I(x^2), level=0.95, se=FALSE, fullrange=FALSE, colour="gray50")
#  + scale_size_manual(values=c(2,1))



# Sample sizes for overall IP measure
IP.samples <- data.frame(Gram.bound=logical(), 
                          Sample.N=integer(),
                          IP=numeric(),
                          Method=character(),
                          stringsAsFactors=FALSE) 

# sample sizes from 2^3 up to 2^12
#for (e in 10:16) {
#  n <- 2^e
for (n in seq(1000, 77000, 1000)) {
  #shuffle corpus rows so the split will be different next time
  #set.seed(27)
  de <- de[sample(nrow(de)), ]
  
  # split corpus into random samples of this size
  samples <- split(de, (seq(nrow(de))-1) %/% n) # creates a list of dfs
  l <- length(samples)
  samples[l] <- NULL # remove the last one, because this one is the residue after all samples of size n have been taken
  
  for (sm in samples) {
    # for each sample, calculate the IP score of free and bound construction types; add it to IP samples
    compounds <- as.data.frame(table(sm$Stem, sm$Preverb, sm$Preverb.bound, sm$Stem.bound, sm$Preverb.separable))
    colnames(compounds) <- c("Lex", "Gram", "Gram.bound", "Lex.bound", "Gram.separable", "Freq")
    compounds <- compounds[compounds$Freq > 0,]
    compounds.freelex <- subset(compounds, Lex.bound==FALSE)
    compounds.freelex <- droplevels(compounds.freelex)
    cmps <- compounds.freelex
    
    # with the restricted dataset calculate independent weighted average entropy of grams, free and bound
    grams <- ddply(cmps, .(Gram, Gram.bound), summarise, Gram.freq = sum(Freq))
    ent.grams.free.CS <- entropy.ChaoShen(subset(grams, Gram.bound==FALSE)$Gram.freq, unit = "log2")
    ent.grams.bound.CS <- entropy.ChaoShen(subset(grams, Gram.bound==TRUE)$Gram.freq, unit = "log2")
    ent.grams.free <- entropy.empirical(subset(grams, Gram.bound==FALSE)$Gram.freq, unit = "log2")
    ent.grams.bound <- entropy.empirical(subset(grams, Gram.bound==TRUE)$Gram.freq, unit = "log2")
    
    # calculate independent weighted average entropy of lexes, in free and bound complexes
    lexes <- ddply(cmps, .(Lex, Gram.bound), summarise, Lex.freq = sum(Freq))
    ent.lexes.freeGram.CS <- entropy.ChaoShen(subset(lexes, Gram.bound==FALSE)$Lex.freq, unit = "log2")
    ent.lexes.boundGram.CS <- entropy.ChaoShen(subset(lexes, Gram.bound==TRUE)$Lex.freq, unit = "log2")
    ent.lexes.freeGram <- entropy.empirical(subset(lexes, Gram.bound==FALSE)$Lex.freq, unit = "log2")
    ent.lexes.boundGram <- entropy.empirical(subset(lexes, Gram.bound==TRUE)$Lex.freq, unit = "log2")
    
    # calculate conditional entropy of grams for each lex
    lex.ent <- ddply(cmps, .(Lex, Gram.bound), summarise, Lex.freq = sum(Freq), Type.freq = length(Lex), Ent = entropy.empirical(Freq, unit = "log2"), Ent.CS = entropy.ChaoShen(Freq, unit = "log2"))
    
    # proportional predictability
    lex.ent$MI.CS.prop <- rep(NA, nrow(lex.ent))
    lex.ent$MI.CS.prop[lex.ent$Gram.bound==FALSE] <- (ent.grams.free.CS - lex.ent$Ent.CS[lex.ent$Gram.bound==FALSE]) / ent.grams.free.CS 
    lex.ent$MI.CS.prop[lex.ent$Gram.bound==TRUE] <- (ent.grams.bound.CS - lex.ent$Ent.CS[lex.ent$Gram.bound==TRUE]) / ent.grams.bound.CS
    
    lex.ent$MI.prop <- rep(NA, nrow(lex.ent))
    lex.ent$MI.prop[lex.ent$Gram.bound==FALSE] <- (ent.grams.free - lex.ent$Ent[lex.ent$Gram.bound==FALSE]) / ent.grams.free 
    lex.ent$MI.prop[lex.ent$Gram.bound==TRUE] <- (ent.grams.bound - lex.ent$Ent[lex.ent$Gram.bound==TRUE]) / ent.grams.bound
    
    # calculate weighted average conditional entropy of grams, i.e. H(G|L), separating free and bound
    lex.ent$weighted.ent.CS <- rep(NA, nrow(lex.ent))
    lex.ent$weighted.ent.CS[lex.ent$Gram.bound==FALSE] <- lex.ent$Ent.CS[lex.ent$Gram.bound==FALSE] * (lex.ent$Lex.freq[lex.ent$Gram.bound==FALSE] / sum(lex.ent$Lex.freq[lex.ent$Gram.bound==FALSE]))
    lex.ent$weighted.ent.CS[lex.ent$Gram.bound==TRUE] <- lex.ent$Ent.CS[lex.ent$Gram.bound==TRUE] * (lex.ent$Lex.freq[lex.ent$Gram.bound==TRUE] / sum(lex.ent$Lex.freq[lex.ent$Gram.bound==TRUE]))
    
    lex.gram.CE.CS.freeGram <- sum(lex.ent$weighted.ent.CS[lex.ent$Gram.bound==FALSE])
    #MI = mutual information, i.e. H(G) - H(G|L)
    lex.gram.MI.CS.freeGram <- ent.grams.free.CS - lex.gram.CE.CS.freeGram
    lex.gram.MI.CS.prop.freeGram <- lex.gram.MI.CS.freeGram / ent.grams.free.CS
    lex.gram.MI.CS.prop.lex.freeGram <- lex.gram.MI.CS.freeGram / ent.lexes.freeGram.CS
    
    IP.CS.free <- lex.gram.MI.CS.freeGram / ent.grams.free.CS
    
    lex.gram.CE.CS.boundGram <- sum(lex.ent$weighted.ent.CS[lex.ent$Gram.bound==TRUE])
    lex.gram.MI.CS.boundGram <- ent.grams.bound.CS - lex.gram.CE.CS.boundGram
    lex.gram.MI.CS.prop.boundGram <- lex.gram.MI.CS.boundGram / ent.grams.bound.CS
    lex.gram.MI.CS.prop.lex.boundGram <- lex.gram.MI.CS.boundGram / ent.lexes.boundGram.CS
    
    IP.CS.bound <- lex.gram.MI.CS.boundGram / ent.grams.bound.CS
    
    #write rows
    IP.samples[nrow(IP.samples)+1,] = list(FALSE, n, IP.CS.free, "Chao-Shen")
    IP.samples[nrow(IP.samples)+1,] = list(TRUE, n, IP.CS.bound, "Chao-Shen")
    
    # Repeat the last bits, but now with empirical entropy
    lex.ent$weighted.ent <- rep(NA, nrow(lex.ent))
    lex.ent$weighted.ent[lex.ent$Gram.bound==FALSE] <- lex.ent$Ent[lex.ent$Gram.bound==FALSE] * (lex.ent$Lex.freq[lex.ent$Gram.bound==FALSE] / sum(lex.ent$Lex.freq[lex.ent$Gram.bound==FALSE]))
    lex.ent$weighted.ent[lex.ent$Gram.bound==TRUE] <- lex.ent$Ent[lex.ent$Gram.bound==TRUE] * (lex.ent$Lex.freq[lex.ent$Gram.bound==TRUE] / sum(lex.ent$Lex.freq[lex.ent$Gram.bound==TRUE]))
    
    lex.gram.CE.freeGram <- sum(lex.ent$weighted.ent[lex.ent$Gram.bound==FALSE])
    #MI = mutual information, i.e. H(G) - H(G|L)
    lex.gram.MI.freeGram <- ent.grams.free - lex.gram.CE.freeGram
    lex.gram.MI.prop.freeGram <- lex.gram.MI.freeGram / ent.grams.free
    lex.gram.MI.prop.lex.freeGram <- lex.gram.MI.freeGram / ent.lexes.freeGram
    
    IP.free <- lex.gram.MI.freeGram / ent.grams.free
    
    lex.gram.CE.boundGram <- sum(lex.ent$weighted.ent[lex.ent$Gram.bound==TRUE])
    lex.gram.MI.boundGram <- ent.grams.bound - lex.gram.CE.boundGram
    lex.gram.MI.prop.boundGram <- lex.gram.MI.boundGram / ent.grams.bound
    lex.gram.MI.prop.lex.boundGram <- lex.gram.MI.boundGram / ent.lexes.boundGram
    
    IP.bound <- lex.gram.MI.boundGram / ent.grams.bound
    
    #write rows
    IP.samples[nrow(IP.samples)+1,] = list(FALSE, n, IP.free, "Empirical")
    IP.samples[nrow(IP.samples)+1,] = list(TRUE, n, IP.bound, "Empirical")
  }
}

#asymptotes
final.free <- IP.samples$IP[IP.samples$Sample.N==77000 & IP.samples$Gram.bound==FALSE & IP.samples$Method=="Chao-Shen"]
final.bound <- IP.samples$IP[IP.samples$Sample.N==77000 & IP.samples$Gram.bound==TRUE & IP.samples$Method=="Chao-Shen"]

# Plotting with boundedness colouring
ggplot(IP.samples, aes(Sample.N, IP)) + geom_jitter(aes(shape=Method), size=2, width=0, height=0) + scale_shape_manual(values=c(16, 4)) +  scale_y_continuous(limits = c(0.4, 0.9)) + geom_hline(yintercept=final.free, linetype="dashed", color = "black")  + geom_hline(yintercept=final.bound, linetype="dashed", color = "black") + annotate(geom="text", x=n-5000, y=final.free+0.04, label="Phrase", size = 10) + annotate(geom="text", x=n-5000, y=final.bound+0.04, label="Word", size = 10) + scale_x_continuous(breaks = c(1000, 5000, 10000, 20000, 40000, 60000, 80000)) + theme_bw() + theme(text = element_text(size = 24))  + xlab("Sample size") + ylab("Internal predictability")
#+ stat_function(fun=function(x) -log2(1/2^x), colour="black")
#  + geom_smooth(aes(group=Gram.bound, colour=Gram.bound), method="lm", formula = y ~ x + I(x^2), level=0.95, se=FALSE)

#geom_jitter(aes(shape=Method), size=2, width=0, height=0.1) + scale_shape_manual(values=c(16, 4)) + 

    