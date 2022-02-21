from typing import Counter
from Bio.SubsMat import MatrixInfo as mf
import numpy as np 

#Get score form blossom
def blossom62Score(row , col):
        Blossom62 = mf.blosum62
        try:
            val = Blossom62[row , col]
        except KeyError:
            val = Blossom62[col , row]
        return val

#Read dataset
def read_data(path):
        seqData= {}
        with open(path , "r") as f:
                for line in f:
                        (key, seq) = line.split()
                        seqData[int(key)] = seq
        return seqData

#Step 1: Remove Low complexity regions
def remove_lowComplx(qSeq):
        rep = []
        #getting all scores
        for i in range(len(qSeq)-1):
                fletter = qSeq[i]
                sletter = qSeq[i+1]
                val = blossom62Score(fletter , sletter)
                rep.append([val , i, i+1])

        #finding similar neighborhing aa
        cc=0
        qlowSeq=""
        i = 0 
        while cc < len(qSeq):
                if i < len(rep):
                        if rep[i][0] > 1 and rep[i+1][0] > 1:
                                qlowSeq = qlowSeq +"XX"
                                cc+=2
                        else: 
                                qlowSeq += qSeq[cc]
                                cc+=1
                if i>= len(rep):
                        qlowSeq += qSeq[cc]
                        cc+=1
                i+=1

        # valid_rep = np.array(rep)
        # indexes = valid_rep>=1
        # Tcount =[]
        # counter =0
        # for i in indexes:
        #         if i == False and len(Tcount)>=3:
        #                 continue
        #         elif i == True:
        #                 Tcount.append(counter)
        #                 counter = counter+1
        #                 continue
        #         counter = counter+1
        #         # else:
        #         #         Tcount.clear()
        #         #         counter = counter+1
        # #if found, exclude 'em and save what they are
        # if len(Tcount) >= 3:
        #         qSeqLowCom = ""
        #         xseq =""
        #         for i in range(len(qSeq)):
        #                 if i in Tcount:
        #                         qSeqLowCom += ""
        #                         xseq += qSeq[i]
        #                 else:
        #                         qSeqLowCom+=(qSeq[i])
        #         return (qSeqLowCom, xseq)
        # else:
        #         return qSeq, "NONE"

        return qlowSeq
        
#Step 2: Make a W-letter word list of the query sequence
def Listing_Words(seq):
    wordList = []
    wordListindex = []
    for i in range(len(seq) - 2):
         wordList.append(seq[i]+seq[i+1]+seq[i+2])
         wordListindex.append([i , i+1 ,i+2])
    return (wordList , wordListindex)

#Step 3&4: 
def Finding_seeds(wordlist, wordListindex , t):
        aa = ["L","I","V","G","A","P","Q","N",
           "M","T","S","C","E","D","K","R","Y","F","W","H"]
        neighborhoodWords = []
        newListofINdex = []

        for i in range(len(wordlist)):
                neighborhoodWords.append([])
                newListofINdex.append([])
                w = wordlist[i]
                neighborhoodWords[i].append(w)
                newListofINdex[i].append(wordListindex[i])
                for j in range(len(aa)):
                        if aa[j] != w[0]:
                                if w[0] == "X":
                                       if aa[j] != qSeq[wordListindex[i][0]]:
                                                neighborhoodWords[i].append(aa[j]+w[1]+w[2])
                                                newListofINdex[i].append(wordListindex[i])
                                else:
                                        neighborhoodWords[i].append(aa[j]+w[1]+w[2])
                                        newListofINdex[i].append(wordListindex[i])
                        if aa[j] != w[1]:
                                if w[1] == "X":
                                       if aa[j] != qSeq[wordListindex[i][1]]:
                                                neighborhoodWords[i].append(w[0]+aa[j]+w[2])
                                                newListofINdex[i].append(wordListindex[i])
                                else:
                                        neighborhoodWords[i].append(w[0]+aa[j]+w[2])
                                        newListofINdex[i].append(wordListindex[i])
                        if aa[j] != w[2]:
                                if w[2] == "X":
                                       if aa[j] != qSeq[wordListindex[i][2]]:
                                                neighborhoodWords[i].append(w[0]+w[1]+aa[j])
                                                newListofINdex[i].append(wordListindex[i])
                                else:
                                        neighborhoodWords[i].append(w[0]+w[1]+aa[j])
                                        newListofINdex[i].append(wordListindex[i])

        #calculating their scores
        listscores = []
        for i in range(len(wordlist)):
                parole = wordlist[i]
                listscores.append([])
                for j in range(len(neighborhoodWords[i])):
                        word = neighborhoodWords[i][j]
                        s0 = blossom62Score(parole[0],word[0])
                        s1 = blossom62Score(parole[1],word[1])
                        s2 = blossom62Score(parole[2],word[2])
                        listscores[i].append(s0+s1+s2)
        
        #collect the seeds
        # Lscores = np.array(listscores)
        # NeighborhoodWords = np.array(neighborhoodWords)
        # Lindices = np.array(newListofINdex)
        # cond = Lscores >= t
        # Lscores = Lscores[cond]
        # NeighborhoodWords = NeighborhoodWords[cond]
        # Lindices = Lindices[cond]

        # Lscores = Lscores.tolist()
        # Lindices = Lindices.tolist()
        # NeighborhoodWords = NeighborhoodWords.tolist()
        Lscores = []
        Lindices = []
        seeds = []
        
        for i in range(len(neighborhoodWords)):
                for j in range(len(neighborhoodWords[i])):
                        if listscores[i][j]>=t:
                                Lscores.append(listscores[i][j])
                                Lindices.append(newListofINdex[i][j])
                                seeds.append(neighborhoodWords[i][j])


        return Lscores,Lindices,seeds

#Step 5
def finding_hits(seqData, Lseeds):
        hits = []
        for i in range(len(seqData)):
                hits.append([])
                s =seqData[i]
                for j in range(len(Lseeds)):
                        num = s.find(Lseeds[j])
                        hits[i].append(num)
        return hits

#Step 6
def extension(xthreshold,qSeq, seqData, Lseeds,seedsScore, hits, wordlistIndex):
        lHSPs = []
        for i in range(len(seqData)):
                s = seqData[i]
                for j in range(len(Lseeds)):
                        HSP= Lseeds[j]
                        score = seedsScore[j]
                        #if seq has no hits
                        if hits[i][j] == -1:
                                continue
                        #position of seed in relation to qSeq
                        leftofSeed = wordlistIndex[j][0]
                        rightofSeed = wordlistIndex[j][-1]
                        fcounter = hits[i][j] + 3
                        bcounter = hits[i][j]-1
                        # if seq has hits
                        for k in range(len(qSeq)):
                                # backward counter
                                bscore=0
                                # forward counter
                                fscore=0
                                #right alignment
                                condition1 = fcounter <len(s) and rightofSeed+1+k <len(qSeq)
                                #left alignment
                                condition2 = bcounter >= 0 and leftofSeed-1-k>=0
                                if condition1 == True:
                                        HSP += qSeq[rightofSeed+1+k]
                                        fscore = blossom62Score(qSeq[rightofSeed+1+k] , s[fcounter])
                                        fcounter+=1
                                if condition2 == True:
                                        HSP = qSeq[leftofSeed-1-k] + HSP
                                        bscore = blossom62Score(qSeq[leftofSeed-1-k] , s[bcounter])
                                        bcounter-=1
                                score += fscore + bscore

                                if score >= xthreshold:
                                        if condition1 == True:
                                              HSP = HSP[:len(HSP)-1]
                                              score -= fscore
                                        if condition2 == True:
                                                HSP =HSP[1:]
                                        lHSPs.append([i, HSP, [bcounter+1 , fcounter-1], score])
                                        break
                                ## if seq in database est fini
                                if condition1 == False and condition2 == False:
                                        lHSPs.append([i, HSP, [bcounter+1 , fcounter-1], score])
                                        break
        return lHSPs

#step 7
def displayOutputs(HSPs, seqData, ThresholdHSP):
        x = 0
        rev = 0
        lrev = []
        for i in range(len(HSPs)):
                if HSPs[i][3]<=ThresholdHSP:
                        remvoed = HSPs[i]
                        lrev.append(remvoed)
                        
        for i in range(len(lrev)):
                HSPs.remove(lrev[i])

        #remove HSPs with same places but only with high score are left out 
        HSPs = sorted(HSPs, key = lambda x: (x[0], x[2]))
        for l in range(len(HSPs)):
                for m in range(len(HSPs)):
                        if HSPs[l][0]== HSPs[m][0] and HSPs[l][2]== HSPs[m][2]:
                                if HSPs[l][3] > HSPs[m][3]:
                                        HSPs[m] =[0,0,0,0]
                                        rev+=1
                                elif HSPs[l][3] < HSPs[m][3]:
                                        HSPs[l]=[0,0,0,0]
                                        rev+=1
        for m in range(rev):
                HSPs.remove([0,0,0,0])
        
        print(HSPs)

        while x<len(HSPs):
                seq = ""
                # num of seq in database
                seqNum = HSPs[x][0]
                print(seqData[seqNum])
                # start and end of HSPs to seq
                start = HSPs[x][2][0]
                end = HSPs[x][2][1]
                totalScore = 0
                flag = True
                #is used for overlap HSPs
                split = 0
                c=0
                while c < len(seqData[seqNum]):
                        if  c>= start and c<=end:
                                seq += HSPs[x][1][split:] + " "
                                c+=len(HSPs[x][1][split:])
                        else:
                                seq += " "
                
                        if c > end:
                                #list of HSP est fini
                                if x+1==len(HSPs):
                                        flag = False
                                #when one seq has more than one HSP
                                elif HSPs[x+1][0] == seqNum:
                                        x+=1
                                        totalScore += HSPs[x][3]
                                        #when they overlap
                                        if (HSPs[x][2][0]< end):
                                                split = end - HSPs[x][2][0] + 1 
                                                start = end + 1
                                        else:
                                                start = HSPs[x][2][0]
                                                split = 0
                                        end = HSPs[x][2][1]
                                        flag = True
                                else:
                                        
                                        flag = False
                                
                        if flag == False:
                                totalScore += HSPs[x][3]
                                print(seq)
                                print("Total score = " + str(totalScore))
                                print("")
                                x+=1
                                break
                        c+=1

def BLAST(qSeq, wordThreshold,ThresholdHSP):
        seqData = read_data("D:\FCIS\Fall 2021\IBT\Labs\project\dataset.txt")
        qSeqLowCom, xseq = remove_lowComplx(qSeq)
        wordList , wordListindex = Listing_Words(qSeqLowCom)
        Lscores,Lindices,NeighborhoodWords = Finding_seeds(wordList,wordListindex, wordThreshold)
        hits = finding_hits(seqData,NeighborhoodWords)
        flag = False
        for i in range(len(hits)):
                for j in range(len(hits[i])):
                        if hits[i][j] != -1:
                                flag = True
        
        if flag == False:
                print("Query sequnce is not found in our dataset")
        else:
                lHSPs = extension(ThresholdHSP,qSeq, seqData,NeighborhoodWords,Lscores,hits,Lindices)
                displayOutputs(lHSPs,seqData, ThresholdHSP)



seqData = read_data("D:\FCIS\Fall 2021\IBT\Labs\project\dataset.txt")

qSeq = "PQGEFG"
qSeqLowCom = remove_lowComplx(qSeq)
# print(qSeqLowCom)
# print(len(qSeqLowCom))
# print(len(qSeq))
wordList , wordListindex = Listing_Words(qSeqLowCom)
Lscores,Lindices,NeighborhoodWords = Finding_seeds(wordList,wordListindex, 15)

hits = finding_hits(seqData,NeighborhoodWords)
lHSPs = extension(10,qSeq, seqData,NeighborhoodWords,Lscores,hits,Lindices)
print(lHSPs)
#displayOutputs(lHSPs,seqData, 16)




