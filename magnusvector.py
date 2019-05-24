import numpy as np
from numpy import empty
np.set_printoptions(threshold=np.inf)




def mainmagnus(stri):
    # Initial DNA Sequence (can be upper or lower case)
    inputfilename='testb'+stri+'.txt'
    file = open(inputfilename, "r")
    dna = file.read()

    #Output file
    outputfilename='bacoutput'+stri+'.txt'
    output = open(outputfilename, "w")

    #testing
    #dna='accccgc'

    #remove numbers and spaces, and linebreaks
    dna = dna.replace(" ","")
    dna = dna.replace('\n', '').replace('\r', '')

    dna = ''.join([i for i in dna if not i.isdigit()])
    dna = dna.upper()

    output.write('DNA Sequence:\n')
    output.write(dna+'\n')

    #Length of dna
    dnalength = len(dna)

    output.write('DNA Length:\n')
    output.write(str(dnalength)+'\n')

    #Length of each window (i.e. k for k-mer)
    N=5
    output.write('Window Length:\n')
    output.write(str(N)+'\n')

    #Length of Magnus Vector
    mvl = (4**(N+1)-4)/3

    output.write('Length of Magnus Vector (for each window):\n')
    output.write(str(mvl)+'\n')

    #Number of Windows
    numwin=dnalength//N
    output.write('Number of Windows:\n')
    output.write(str(numwin)+'\n')

    #Initialize Array of Windows
    winarray = empty([numwin,mvl])
    #print(winarray)



    #print(dnabase4)

    # bottom-up, dynamic programming solution using a single array
    #Reference: https://stackoverflow.com/questions/6877249/find-the-number-of-occurrences-of-a-subsequence-in-a-string

    def num_subsequences(seq, sub):
        m, n = len(seq), len(sub)
        table = [0] * n
        for i in xrange(m):
            previous = 1
            for j in xrange(n):
                current = table[j]
                if seq[i] == sub[j]:
                    table[j] += previous
                previous = current
        return table[n-1] if n else 1

    #Change base
    #Reference: https://stackoverflow.com/questions/2267362/how-to-convert-an-integer-in-any-base-to-a-string
    def numberToBase(n, b):
        if n == 0:
            return '0'
        digits = ''
        while n:
            digits+=str(int(n % b))
            n //= b
        return digits[::-1]

    for globalcounter in range(0,numwin):
        startindex=globalcounter*N
        endindex=startindex+N
        #Base 4 DNA
        dnabase4 = ''

        #Magnus Vector
        magnusvec = [0] * mvl

        for i in range(startindex,endindex):
            if (dna[i]=='A'):
                dnabase4+='0'
            elif (dna[i]=='C'):
                dnabase4+='1'
            elif(dna[i]=='G'):
                dnabase4+='2'
            elif(dna[i]=='T' or dna[i]=='U'):
                dnabase4+='3'
            else:
                print('DNA Sequence contains unallowed letters.')
                break

        counter = 0
        maxds = 4**N
        while counter<mvl:
            for ds in range(0,maxds):
                s = numberToBase(ds,4)
                while N-len(s)>=0:
                    #print(s)
                    ps = (4**len(s)-4)/3 +ds + 1
                    #print(ps)
                    alphas = num_subsequences(dnabase4,s)
                    magnusvec[ps-1]=alphas
                    s = '0'+s
                    counter+=1



        winarray[globalcounter]=magnusvec
        #end loop

    #print("Magnus Vectors for R=Z:")
    #print(winarray)


    summagnus=[0]*mvl
    for i in range(0,len(winarray)):
        summagnus=summagnus+winarray[i]
    meanmagnus=summagnus/numwin
    output.write('Mean Magnus Vector\n')
    output.write(str(meanmagnus)+'\n')

    output.close()

    npyname='mean'+stri
    np.save(npyname,meanmagnus)






for j in range(1,2):
    mainmagnus(str(j))


######
#Short Magnus Window
'''
shortmagnusvec=[mvl]
for i in range(0,len(magnusvec)):
    if magnusvec[i]!=0:
        shortmagnusvec.append(i+1)

print('Short Magnus Vector for R=Z:')
print(shortmagnusvec)
'''

#Mod 2 Magnus
'''
def reduceMod2(l):
    return [i % 2 for i in l]


magnusvecmod2=reduceMod2(magnusvec)
print("Magnus Vector for R=Z/2:")
print(magnusvecmod2)

shortmagnusvecmod2=[mvl]
for i in range(0,len(magnusvecmod2)):
    if magnusvecmod2[i]!=0:
        shortmagnusvecmod2.append(i+1)

print('Short Magnus Vector for R=Z/2:')
print(shortmagnusvecmod2)
'''