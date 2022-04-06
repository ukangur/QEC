from unittest import result
import numpy as np
from numpy.linalg import matrix_rank
from collections import OrderedDict
import time

def parseToBinary(inputlist):
    fullmatrix = []
    for i in inputlist:
        outputleft = []
        outputright = []
        for j in i:
            binarytuple = pauliToBinary(j)
            outputleft.append(binarytuple[0])
            outputright.append(binarytuple[1])
        outputrow = outputleft + outputright
        fullmatrix.append(outputrow)
    return fullmatrix


def pauliToBinary(pchar):
    if pchar == 'I':
        return (0,0)
    elif pchar == 'X':
        return (1,0)
    elif pchar == 'Z':
        return (0,1)
    elif pchar == 'Y':
        return (1,1)

def binaryToPauli(left,right):
    if (left,right) == (0,0):
        return 'I'
    elif (left,right) == (1,0):
        return 'X'
    elif (left,right) == (0,1):
        return 'Z'
    elif (left,right) == (1,1):
        return 'Y'

def parseToPauli(left,right):
    operator = ''
    for i in range(len(left)):
        operator += binaryToPauli(left[i],right[i])
    return operator

def ge(M,right):

    qswaps = []
    ops = []

    m,n = M.shape
    r = matrix_rank(M)

    if(r == n-1 and right):
        return [M, qswaps, ops]

    oldrank = r

    i=0
    j=0

    if(right):
        i = r
        j = r
        r = n-1

    while i < r and j < r:

        #if pivot is zero then
        if(M[i,i] == 0):

            #find first index, where column element is 1
            k = np.argmax(M[i:, j]) +i

            # first try to swap rows
            temp = np.copy(M[k])
            M[k] = M[i]
            M[i] = temp
            ops.append(["r",i,k])

        #if pivot still zero, then
        if(M[i,i] == 0):

            #swap columns
            k = np.argmax(M[i,oldrank:]) + oldrank
            temp = np.copy(M[:,k])
            qswaps.append([i,k])
            ops.append(["c",i,k])
            M[:,k] = M[:,i]
            M[:,i] = temp

        if(i == r-1):
            #backsubstitution to get upper triangle into 0-s
            if(right):
                r = oldrank
            else:
                r = 0

            while i > r and j > r:
                aijn = M[i, j:]
                col = np.copy(M[:, j])

                col[i] = 0
                ops.append(["op",i,j,col])    
                flip = np.outer(col, aijn)

                M[:, j:] = M[:, j:] ^ flip

                i -= 1
                j -= 1        
            break
        
        #print("M: ",M)
        aijn = M[i, j:]
        #print("ajin: ",aijn)
        col = np.copy(M[:, j])
        #print("col: ",col)
        col[i] = 0

        ops.append(["op",i,j,col])

        flip = np.outer(col, aijn)
        #print("flip: ",flip)

        #print("M[:, j:]: ",M[:, j:])
        M[:, j:] = M[:, j:] ^ flip

        #print("M[:, j:]: ",M[:, j:])

        i += 1
        j += 1 
    return [M, qswaps, ops]

def applyoperationsmatrix(M,ops):
    for op in ops:
        i = op[1]
        j = op[2]

        if(op[0] == "r"):
            temp = np.copy(M[j])
            M[j] = M[i]
            M[i] = temp

        elif(op[0] == "c"):
            temp = np.copy(M[:,j])
            M[:,j] = M[:,i]
            M[:,i] = temp

        elif(op[0] == "op"):
            aijn = M[i, j:]
            col = op[3]
            flip = np.outer(col, aijn)
            M[:, j:] = M[:, j:] ^ flip

    return M

def applyswapsvector(M,qswaps):
    for swap in qswaps:
        i = swap[0]
        k = swap[1]
        temp = M[k]
        M[k] = M[i]
        M[i] = temp
    return M

def findcommutinglist(elemlist):

    commutinglists = []
    commdict = {}

    n = len(elemlist[0])

    lambdamatrix = np.zeros(shape=(n,n))
    x = 0
    for i in range(int(n/2),n):
        lambdamatrix[i,x] = 1
        lambdamatrix[x,i] = 1
        x += 1

    for indi,i in enumerate(elemlist):
        subcommutelist = []
        for indj,j in enumerate(elemlist):
            x,y = len(i),len(j)
            if(x != y):
                print("Pauli group elements in your input must be of same length. Program will exit now.")
                raise SystemExit(0)
            if(np.array(i).dot(lambdamatrix).dot(np.array(j)) % 2 == 0):
                subcommutelist.append(indj)
        commdict[indi] = subcommutelist
        commutinglists.append(subcommutelist)

    commdictordered = OrderedDict(sorted(commdict.items(), key=lambda x: len(x[1]), reverse=True))
    commutinglists = sorted(commutinglists, key=len)
    output = []

    counter = 0
    keysordered = list(commdictordered.keys())

    while counter < len(commdictordered):
        smallcounter = 0
        nextlist = commdictordered[keysordered[counter]]
        for i in nextlist:
            for j in nextlist:
                if(i in commdictordered[j]):
                    smallcounter += 1
                else:
                    break
            else:
                continue
            break
        if smallcounter == len(nextlist)**2:
            for i in nextlist:
                output.append(elemlist[i])
            break
        counter += 1
            
    return output      

def slicematrix(fullM):
    leftM = []
    rightM = []
    n = int(len(fullM[0])/2)

    for i in fullM:
        leftM.append(i[:n])
        rightM.append(i[n:])

    return leftM, rightM

def getNiceLogicalPauliOps(commutinglist,logicallist):

    logicaloptions = set()

    for i in commutinglist:
        a = tuple(x^y for x,y in zip(logicallist,i))
        logicaloptions.add(a)
        for j in commutinglist:
            b = tuple(x^y for x,y in zip(a,j))
            logicaloptions.add(b)
            for p in commutinglist:
                c = tuple(x^y for x,y in zip(b,p))
                logicaloptions.add(c)

    logicaloptionsinPauli = []

    for i in logicaloptions:
        op = parseToPauli(i[:n],i[n:])
        if('X' in op and 'Z' in op or 'X' in op and 'Y' in op or 'Z' in op and 'Y' in op):
            continue
        else:
            logicaloptionsinPauli.append(parseToPauli(i[:n],i[n:]))
    
    return logicaloptionsinPauli

# -- HERE STARTS THE MAIN PROGRAM --

first_times = []
second_times = []
third_times = []
fourth_times = []

# FIRST I READ THE INPUT FILE AND GET THE DATA
begin = time.perf_counter()
input = []

f = open("input.txt", "r")
for line in f.readlines():
    input.append(line.strip())
f.close()

n = len(input[0])

# THEN I TRANSFORM MY INPUT INTO ITS BINARY REPRESENTATION
# I SEPARATE FULL, LEFT AND RIGHT MATRICES (BUT I ONLY USE LEFT AND RIGHT)

full = parseToBinary(input)
end = time.perf_counter()
first_times.append(end-begin)
# WE TAKE THE INPUT PAULI GROUP ELEMENTS AND FIND THE LARGEST SUBSET OF COMMUTING PAULI GROUP ELEMENTS. IF THE NUMBER IS LESS THAN N-1, THEN THE PROGRAM ENDS.
begin = time.perf_counter()
commutinglist = findcommutinglist(full)
if(len(commutinglist) < n-1):
    print("The input does not include at least n-1 commuting pauli group elements, thus we can not calculate single logical X and Z operations")
    raise SystemExit(0)
end = time.perf_counter()
second_times.append(end-begin)
# WE SLICE THE FULL MATRIX INTO HALF (LEFT AND RIGHT MATRIX)
begin = time.perf_counter()
left, right = slicematrix(commutinglist)
leftmatrix = np.matrix(left)
rightmatrix = np.matrix(right)
rleft = matrix_rank(leftmatrix)
rright = matrix_rank(rightmatrix)

# THEN I DO THE GAUSSIAN ELIMINATION ON THE LEFT MATRIX UPPER SUBMATRIX
resultleft = ge(leftmatrix,False)
changedleftmatrix = resultleft[0]
qswapsleft = resultleft[1]
opsleft = resultleft[2]

# I APPLY THE SAME OPERATIONS TO THE RIGHT MATRIX

changedrightmatrix = applyoperationsmatrix(rightmatrix,opsleft)

# THEN I DO THE GAUSSIAN ELIMINATION ON THE RIGHT MATRIX BOTTOM-MIDDLE SUBMATRIX

if(n-1 != rright):
    resultright = ge(changedrightmatrix,True)
    changedrightmatrix = resultright[0]
    qswapsright = resultright[1]
    opsright = resultright[2]

else:
    changedrightmatrix = rightmatrix
    qswapsright = []
    opsright = []

# I APPLY THE SAME OPERATIONS TO THE LEFT MATRIX

changedleftmatrix = applyoperationsmatrix(changedleftmatrix,opsright)

# I APPEND THE RIGHT AND LEFT MATRIX QUBIT SWITCHES FOR LATER USE

qswapstotal = []
qswapstotal.extend(qswapsleft)
qswapstotal.extend(qswapsright)

#I GET THE Z COLUMN PART A2 OF LEFT MATRIX
z_arrA2 = changedleftmatrix[:rleft,n-1]

#I USE THE LOGICAL Z EQUATION WITH MY LEFT MATRIX COLUMN

leftz = []
for i in range(n):
    leftz.append(0)

rightz = []
for i in range(len(z_arrA2)):
    rightz.append(z_arrA2[i,0])
for i in range(n-1-rleft):
    rightz.append(0)
rightz.append(1)

#I DO THE QUBIT SWITCHES ON THE LOGICAL Z EQUATION

leftz = applyswapsvector(leftz,qswapstotal)
rightz = applyswapsvector(rightz,qswapstotal)

#I TRANSFORM THE LOGICAL QUBIT BACK INTO ITS ORIGINAL FORM

listz = leftz + rightz
finalz = parseToPauli(leftz,rightz)
print(f"THE INITIAL LOGICAL Z FOR {len(finalz)} QUBITS: {finalz}")

#I GET THE X COLUMN PARTS C and E OF RIGHT MATRIX
x_arrC = changedrightmatrix[:rright,n-1]

if(rright != n-1):
    x_arrE = changedrightmatrix[rright:,n-1]
else:
    x_arrE = []


#I USE THE LOGICAL X EQUATION WITH MY RIGHT MATRIX COLUMNS

leftx = []
for i in range(rright):
    leftx.append(0)
for i in range(len(x_arrE)):
    leftx.append(x_arrE[i,0])
leftx.append(1)

rightx = []
for i in range(len(x_arrC)):
    rightx.append(x_arrC[i,0])
for i in range(n-rright):
    rightx.append(0)

#I DO THE QUBIT SWITCHES ON THE LOGICAL X EQUATION

leftx = applyswapsvector(leftx,qswapstotal)
rightx = applyswapsvector(rightx,qswapstotal)

#I TRANSFORM THE LOGICAL QUBIT BACK INTO ITS ORIGINAL FORM

listx = leftx + rightx
finalx = parseToPauli(leftx,rightx)
print(f"THE INITIAL LOGICAL X FOR {len(finalx)} QUBITS: {finalx}")

#I FIND ALL NICE (NON Z AND X MIXED) LOGICAL OPERATORS

zchoices = getNiceLogicalPauliOps(commutinglist,listz)
xchoices = getNiceLogicalPauliOps(commutinglist,listx)

print(f"THE SET OF ALL NICE LOGICAL Z CHOICES FOR {len(finalz)} QUBITS: {zchoices}")
print(f"THE SET OF ALL NICE LOGICAL X CHOICES FOR {len(finalx)} QUBITS: {xchoices}")

#I WRITE RESULT TO FILE

file1 = open("outputX.txt", "w")
for i in xchoices:
    file1.write(i + '\n')
file1.close()

file2 = open("outputZ.txt", "w")
for i in zchoices:
    file2.write(i + '\n')
file2.close()

