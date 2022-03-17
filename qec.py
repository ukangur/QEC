import numpy as np
from numpy.linalg import matrix_rank

def parseToBinary(inputlist):
    fullmatrix = []
    leftmatrix = []
    rightmatrix = []
    for i in inputlist:
        outputleft = []
        outputright = []
        for j in i:
            binarytuple = pauliToBinary(j)
            outputleft.append(binarytuple[0])
            outputright.append(binarytuple[1])
        leftmatrix.append(outputleft)
        rightmatrix.append(outputright)
        outputrow = outputright.copy().append(outputleft.copy())
        fullmatrix.append(outputrow)
    return fullmatrix, leftmatrix, rightmatrix


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

    m,n = M.shape
    r = matrix_rank(M)
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

        #if pivot still zero, then
        if(M[i,i] == 0):

            #swap columns
            k = np.argmax(M[i,oldrank:]) + oldrank
            temp = np.copy(M[:,k])
            qswaps.append([i,k])
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

                flip = np.outer(col, aijn)

                M[:, j:] = M[:, j:] ^ flip

                i -= 1
                j -= 1        
            break

        aijn = M[i, j:]
        col = np.copy(M[:, j])

        col[i] = 0

        flip = np.outer(col, aijn)

        M[:, j:] = M[:, j:] ^ flip

        i += 1
        j += 1
    return [M, qswaps]

def applyswapsmatrix(M,qswaps):
    for swap in qswaps:
        i = swap[0]
        k = swap[1]
        temp = np.copy(M[:,k])
        M[:,k] = M[:,i]
        M[:,i] = temp
    return M

def applyswapsvector(M,qswaps):
    for swap in qswaps:
        i = swap[0]
        k = swap[1]
        temp = M[k]
        M[k] = M[i]
        M[i] = temp
    return M

# -- HERE STARTS THE MAIN PROGRAM --

# FIRST I READ THE INPUT FILE AND GET THE DATA

input = []

f = open("input.txt", "r")
for line in f.readlines():
    input.append(line.strip())
f.close()

n = len(input)

# THEN I TRANSFORM MY INPUT INTO ITS BINARY REPRESENTATION
# I SEPARATE FULL, LEFT AND RIGHT MATRICES (BUT I ONLY USE LEFT AND RIGHT)

full, left, right = parseToBinary(input)
leftmatrix = np.matrix(left)
rightmatrix = np.matrix(right)

# THEN I DO THE GAUSSIAN ELIMINATION ON THE LEFT MATRIX UPPER SUBMATRIX

resultleft = ge(leftmatrix,False)
changedleftmatrix = resultleft[0]
qswapsleft = resultleft[1]

# I APPLY THE SAME QUBIT SWAPS TO THE RIGHT MATRIX

changedrightmatrix = applyswapsmatrix(rightmatrix,qswapsleft)

# THEN I DO THE GAUSSIAN ELIMINATION ON THE RIGHT MATRIX BOTTOM-MIDDLE SUBMATRIX

resultright = ge(changedrightmatrix,True)
changedrightmatrix = resultright[0]
qswapsright = resultright[1]

# I APPLY THE SAME QUBIT SWAPS TO THE LEFT MATRIX

changedleftmatrix = applyswapsmatrix(changedleftmatrix,qswapsright)

# I APPEND THE RIGHT AND LEFT MATRIX QUBIT SWITCHES FOR LATER USE

qswapstotal = []
qswapstotal.extend(qswapsleft)
qswapstotal.extend(qswapsright)

#I GET THE Z COLUMN PART A2 OF LEFT MATRIX
r = matrix_rank(changedleftmatrix)
z_arrA2 = changedleftmatrix[:r,n]

#I USE THE LOGICAL Z EQUATION WITH MY LEFT MATRIX COLUMN

leftz = []
for i in range(n+1):
    leftz.append(0)

rightz = []
for i in range(len(z_arrA2)):
    rightz.append(z_arrA2[i,0])
for i in range(n-r):
    rightz.append(0)
rightz.append(1)

#I DO THE QUBIT SWITCHES ON THE LOGICAL Z EQUATION

leftz = applyswapsvector(leftz,qswapstotal)
rightz = applyswapsvector(rightz,qswapstotal)

#I TRANSFORM THE LOGICAL QUBIT BACK INTO ITS ORIGINAL FORM

finalz = parseToPauli(leftz,rightz)
print(f"THE FINAL LOGICAL Z FOR 7 QUBITS: {finalz}")

#I GET THE X COLUMN PARTS C and E OF RIGHT MATRIX
r = matrix_rank(changedleftmatrix)
x_arrC = changedleftmatrix[r:,n]
x_arrE = changedleftmatrix[:r,n]

#I USE THE LOGICAL X EQUATION WITH MY RIGHT MATRIX COLUMNS

leftx = []
for i in range(r):
    leftx.append(0)
for i in range(len(x_arrE)):
    leftx.append(x_arrE[i,0])
leftx.append(1)

rightx = []
for i in range(len(x_arrC)):
    rightx.append(x_arrC[i,0])
for i in range(n-r+1):
    rightx.append(0)

#I DO THE QUBIT SWITCHES ON THE LOGICAL X EQUATION

leftx = applyswapsvector(leftx,qswapstotal)
rightx = applyswapsvector(rightx,qswapstotal)

#I TRANSFORM THE LOGICAL QUBIT BACK INTO ITS ORIGINAL FORM

finalx = parseToPauli(leftx,rightx)
print(f"THE FINAL LOGICAL X FOR 7 QUBITS: {finalx}")