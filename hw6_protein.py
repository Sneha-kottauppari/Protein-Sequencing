"""
Protein Sequencing Project
Name: Sneha.K
Roll Number:2021501022
"""

from json import load
from os import read
import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    fp=open(filename,'r')
    text=fp.read()
    temp=text.split('\n')
    file_content=""
    for each in temp:
        file_content=file_content+each
    fp.close()
    return file_content


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    dna=dna.replace('T','U')
    length=len(dna)
    temp_list=[]
    codon_list=[]
    stop_codons=['UAA','UAG','UGA']
    for index in range(startIndex,length,3):
        temp_list.append(dna[index:index+3])
    for each in temp_list:
        if each in stop_codons:
            codon_list.append(each)
            break
        else:
            codon_list.append(each)
    return codon_list


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    fj=open(filename,'r')
    result_dict={}
    json_dict=load(fj)
    for k,v in json_dict.items():
        for each in v:
            each=each.replace('T','U')
            result_dict[each]=k
    fj.close()
    return result_dict


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    protein_list=[]
    stop_codons=['UAA','UAG','UGA']
    for i in range(len(codons)):
        if codons[i]=='AUG'and i==0:
            protein_list.append('Start')
        elif codons[i] in stop_codons or i==len(codons):
            protein_list.append('Stop')
        else:
            protein_list.append(codonD[codons[i]])
    return protein_list


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    dna_in_file=readFile(dnaFilename)
    codon_dict= makeCodonDictionary(codonFilename)
    count=0
    final_protein_list=[]
    i=0
    codons_list=[]
    while i<len(dna_in_file):
        if dna_in_file[i:i+3] == 'ATG':
            start_index=i
            count+=1
            codons_list=dnaToRna(dna_in_file,start_index)
            protein_seq = generateProtein(codons_list,codon_dict)
            final_protein_list.append(protein_seq)
            i=i+len(codons_list)*3
        else: i=i+1
    return final_protein_list


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    unique_list=[]
    for each_list in proteinList1:
        if each_list in proteinList2:
            if each_list not in unique_list:
                unique_list.append(each_list)
    return unique_list


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    collapse_list=[]
    for each_list in proteinList:
        for i in each_list:
            collapse_list.append(i)
    return sorted(collapse_list)


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    aa_count={}
    for each_acid in aaList:
        if each_acid not in aa_count:
            aa_count[each_acid]=1
        else:
            aa_count[each_acid]+=1
    return aa_count


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    result_AA=[]
    ##proteinlist2 list and dict##
    aa_list1=combineProteins(proteinList1)
    aa_dict1=aminoAcidDictionary(aa_list1)
   
    aa_list2=combineProteins(proteinList2)
    aa_dict2=aminoAcidDictionary(aa_list2)
    
    ## differences##
    notinlist1=list(set(aa_list2)-set(aa_list1))
    notinlist2=list(set(aa_list1)-set(aa_list2))

    ##updating dict##
    for each_acid_a in notinlist1:
        aa_dict1[each_acid_a]=0
    for each_acid_b in notinlist2:
        aa_dict2[each_acid_b]=0
    
    ##frequencies##
    length1=len(aa_list1)
    aa_freq1={}
    for each_acid_c in aa_dict1:
        aa_freq1[each_acid_c]=aa_dict1[each_acid_c]/length1

    length2=len(aa_list2)
    aa_freq2={}
    for each_acid_d in aa_dict2:
        aa_freq2[each_acid_d]=aa_dict2[each_acid_d]/length2

    ##result list##
    for k,v in aa_dict1.items():
        inner_list=[]
        if k!='Start' and k!='Stop':
            if k in aa_dict2.keys():
                diff_freq = abs(aa_freq1[k]-aa_freq2[k])
                if diff_freq>cutoff:
                    inner_list.append(k)
                    inner_list.append(aa_freq1[k])
                    inner_list.append(aa_freq2[k])
                    result_AA.append(inner_list)
    return result_AA



'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    ##display for commoalities##
    for each_list in sorted(commonalities):
        each_list=sorted(each_list)
        display_string=""
        for i in range(len(each_list)):
            if each_list[i]!="Start" and each_list[i]!="Stop":
                display_string = display_string + "-"+ str(each_list[i])
        print(display_string.strip('-'))
    print('^'*10)
    ##display for differences##
    for each_diff in sorted(differences):
        freq1=round(each_diff[1]*100,2)
        freq2=round(each_diff[2],2)
        print(each_diff[0], ":", freq1, "% in seq1,", freq2, "% in seq2")
    return


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    return


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    return


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    return


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    # print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    # test.week1Tests()
    # print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    # runWeek1()

    ## Uncomment these for Week 2 ##

    # print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    # test.week2Tests()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()


    ## Uncomment these for Week 3 ##
    """
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    """
    # test.testReadFile()
    # test.testDnaToRna()
    # test.testMakeCodonDictionary()
    # test.testGenerateProtein()
    # test.testSynthesizeProteins()
    # test.testCommonProteins()
    # test.testCombineProteins()
    # test.testAminoAcidDictionary()
    # test.testFindAminoAcidDifferences()
