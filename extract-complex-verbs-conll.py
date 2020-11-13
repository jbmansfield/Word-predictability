# -*- coding: utf-8 -*-
# This tells python that parts of the script are utf-8

"""
extract-complex-verbs-conll.py
"""

import sys, re, os
import pyconll

def usage():
	print("""
#################################################
# USAGE INSTRUCTIONS: extract-complex-verbs-conll.py
#################################################

Must use python3

See here for doc of CONLL format:
https://universaldependencies.org/format.html

This extracts complex verbs from CONLL format treebanks.
It uses the pyconll package, see documentation for details.

""")

outfile = open('hamburger-verbs.txt', 'w')
outfile.write('Stem'+'\t'+'Preverb'+'\t'+'Preverb.bound'+'\t'+'Preverb.separable'+'\t'+'Stem.bound'+'\t'+'Separated'+'\n')

# Lists taken from Dodd et al, pages 142ff. See also Whittle et al 2011: 109, 166
prvBound = ['be', 'emp', 'ent', 'er', 'ge', 'miss', 'ver', 'zer']

# beyond those in Dodd, I also add "bei, fest, fort, frei, hoch, los"
prvFree = ['ab', 'an', 'auf', 'aus', 'bei', 'ein', 'entgegen', 'fern', 'fest', 'fort', 'frei', 'her', 'hin', 'hoch', 'los', 'mit', 'nach', 'vor', 'weg', 'zu', 'zurück', 'zusammen']
prvMixed = ['durch', 'über', 'um', 'unter', 'voll', 'wieder']

prvInsep = ['be', 'emp', 'ent', 'er', 'ge', 'ver', 'zer', 'hinter', 'miss', 'wider']


stopWords = ['herrschen', 'geistern', 'zugen', 'gen', 'beugen', 'gegnen', 'begnügen', 'geben', 'gehen', 'gelten', 'missen', 'mitteln', 'hindern', 'einigen', 'bergen', 'fernen', 'bessern', 'behren', 'antworten', 'wegen', 'analysieren', 'hinken', 'einen', 'freien', 'animieren', 'absolvieren', 'festigen', 'ankern', 'abonnieren', 'absorbieren', 'losen', 'widern', 'angeln', 'annullieren', 'anonymisieren', 'nachten', 'beginnen']
stopForms = ['herr', 'einander', 'ander', 'einheit', 'hintan'] # don't parse a preverb if the string starts like this
derivs = ['recht'] # when checking whether a stem is bound, first remove any other preverbs, and these "derivs"

# See also:
# http://www.dartmouth.edu/~deutsch/Grammatik/Wortbildung/Separables.html
# http://www.dartmouth.edu/~deutsch/Grammatik/Wortbildung/Inseparables.html

contxtCounter = 0

wordsFile = open('derewo-v-ww-bll-320000g-2012-12-31-1.0/derewo-v-ww-bll-320000g-2012-12-31-1.0.txt', 'r')
wordList = {}
lines = wordsFile.readlines()
for line in lines:
    #print(line)
    #skip empty
    if not re.search(r'\w', line):
        continue
    if re.search('#', line):
        continue
    line = re.sub('[\n\r]*', '', line) #chomp
    line = re.sub('"', '', line) #clean out any quotes, because R will choke on them

    if re.search(r'.+? .+? .+?', line):
        (lemma, freq, pos) = line.split(' ')
    else:
        (lemma, freq) = line.split(' ')
    
    if lemma in wordList.keys():
        wordList[lemma] += float(freq)
    else:
         wordList[lemma] = float(freq)

def main():
    # need to build this index of known words from a separate database

    # And here is the treebank corpus, using just the manually annotated section
    tbank = pyconll.load_from_file('hamburg-dependency-treebank-conll/part_A.conll')
 
    # each s may contain one or more compounding relations, because some compounds are recursive like "ab-be-kommen"
    for s in tbank:
        verbs = [] # list of tuples like (verbStem, particle, preverb, separated)
        verbStem = ''
        particle = ''
        verbIndex = ''
        preverb = ''
        separated = False
        #outfile.write(' '.join(w.form for w in s))
        #outfile.write('\n')

        for w in s:
            if w.upos=='V':
                verbStem = w.lemma
                verbIndex = w.id
                verbs.append((verbStem, particle, preverb, separated))

            elif (w.upos=='PTKVZ' and w.head==verbIndex):
                particle = w.lemma
                # update last verb with this particle
                verbs.pop()
                verbs.append((verbStem, particle, preverb, separated))

            elif (w.lemma==',' and w.deprel=='ROOT'):
                # every time you hit a comma root, reset the sentence
                # this may miss a few datapoints, since list commas seem to also be annotated as ROOTs in the corpus 
                for v in verbs:
                    process_verb(v)
                verbs = []
                verbStem = ''
                particle = ''
                verbIndex = ''
                preverb = ''
                separated = False

            
            """# debugging
            if re.search('verabsch', w.form):
                outfile.write('EXX: ')
                outfile.write(' '.join(w.form for w in s))
                outfile.write('\n')"""

        for v in verbs:
            process_verb(v)

        #outfile.write('\n')
        """
        # First we check if the final particle is a preverb, and if so we already have a pair that can be written, flagged as separated
        if (particle in (prvFree + prvMixed)):
            preverb = particle
            separated = True
            write_verb(verbStem, preverb, separated)

        # Now we parse the input verb stem, recursively until it can't parse no more
        parse_de_verb(verbStem)"""


def process_verb (v):
    (verbStem, particle, preverb, separated) = v
    # First we check if the final particle is a preverb, and if so we already have a pair that can be written, flagged as separated
    if (particle in (prvFree + prvMixed)):
        preverb = particle
        separated = True
        write_verb(verbStem, preverb, separated)
    # Now we parse the input verb stem, recursively until it can't parse no more
    parse_de_verb(verbStem, preverb)

def parse_de_verb (verbStem, preverb):
    global contxtCounter
    # first check for stopwords, these ones look like they have preverbs but they don't, so skip to next sentence when you reach one of these
    if verbStem in stopWords:
        # if no preverb has already been identified, then write this as a simple verbs
        if preverb=='':
            #write_verb(verbStem, "NONE"+str(contxtCounter), False)
            write_verb(verbStem, "NONE", False)
            contxtCounter += 1
    else:
        # cycle through stopForms, starting from the longest
        stopped = False
        for sf in sorted(stopForms, key = len, reverse = True):
            rgx = r"^"+sf
            #print(rgx)
            if re.search(rgx, verbStem):
                #write_verb(verbStem, "NONE"+str(contxtCounter), False)
                contxtCounter += 1
                stopped = True
                break
        # If you didn't hit a stopform, then check for preverbs
        if stopped==False:
            # cycle through possible preverbs, starting from the longest
            for prv in sorted((prvInsep + prvMixed + prvFree), key = len, reverse = True):
                #rgx = r"^" + re.escape(prv) + r"([.+])"
                rgx = r"^"+prv
                #print(rgx)
                if re.search(rgx, verbStem):
                    verbStem = re.sub(rgx, '', verbStem)
                    #matchObj = re.match(rgx, mainVerb)
                    #verbStem = matchObj.group(1)
                    preverb = prv
                    write_verb(verbStem, preverb, False)
                    parse_de_verb(verbStem, preverb)
                    break
            if preverb=='':
                #write_verb(verbStem, "NONE"+str(contxtCounter), False)
                write_verb(verbStem, "NONE", False)
            contxtCounter += 1

def write_verb (verbStem, preverb, separated):
    if (verbStem!='' and preverb!=''):
        boundStem = True
        if verbStem in wordList.keys():
            boundStem = False
        for fm in sorted((prvInsep + prvMixed + prvFree + derivs), key = len, reverse = True):
            rgx = r"^"+fm
            if re.search(rgx, verbStem):
                boundStem = False


        boundPrv = False
        if preverb in prvBound:
            boundPrv = True

        separablePrv = 'Separable'
        if preverb in prvInsep:
            separablePrv = 'Attached'
        elif preverb in prvMixed:
            separablePrv = 'Mixed'
        outfile.write(verbStem+'\t'+preverb+'\t'+str(boundPrv).upper()+'\t'+separablePrv+'\t'+str(boundStem).upper()+'\t'+str(separated).upper()+'\n')

if __name__ == "__main__":
    main()
