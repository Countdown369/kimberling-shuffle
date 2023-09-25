# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import pandas as pd
import math
import time

pd.set_option('display.max_columns', 2000)

def genNextRILI(inrow, stage):
    nextrow = []
    cap = 2 * (stage - 1)
    
    if stage == 1:
        cap = 1
    
    for x in range(cap):
        if (x+1) % 2 == 1:
            nextrow.append(inrow[(stage-1) + math.ceil((x+1)/2)])
        else:
            nextrow.append(inrow[(stage-1) - math.ceil((x+1)/2)])
    for y in inrow[(stage+math.ceil(cap/2)):]:
        nextrow.append(y)
    return(nextrow)

def genNextROLO(inrow, stage):
    nextrow = []
    cap = 2 * (stage - 1)
    
    if stage == 1:
        for z in inrow[1:]:
            nextrow.append(z)
        return nextrow
    
    else:
        for x in reversed(range(cap)):
            if (x+1) % 2 == 0:
                nextrow.append(inrow[(stage-1) + int(round((x/2)+.75))])
            else:
                nextrow.append(inrow[(stage-1) - int(round((x/2)+.75))])
        for y in inrow[(stage+math.ceil(cap/2)):]:
            nextrow.append(y)
        return(nextrow)

def genNextLIRI(inrow, stage):
    nextrow = []
    cap = 2 * (stage - 1)
    
    if stage == 1:
        cap = 0
    
    for x in range(cap):
        if (x+1) % 2 == 1:
            nextrow.append(inrow[(stage-1) - math.ceil((x+1)/2)])
        else:
            nextrow.append(inrow[(stage-1) + math.ceil((x+1)/2)])
    for y in inrow[(stage+math.ceil(cap/2)):]:
        nextrow.append(y)
    return(nextrow)

def genNextLORO(inrow, stage):
    nextrow = []
    cap = 2 * (stage - 1)
    
    if stage == 1:
        for z in inrow[1:]:
            nextrow.append(z)
        return nextrow
    
    else:
        for x in reversed(range(cap)):
            if (x+1) % 2 == 0:
                nextrow.append(inrow[(stage-1) - int(round((x/2)+.75))])
            else:
                nextrow.append(inrow[(stage-1) + int(round((x/2)+.75))])
        for y in inrow[(stage+math.ceil(cap/2)):]:
            nextrow.append(y)
        return(nextrow)
    
def genNextROLI(inrow, stage):
    nextrow = []
    cap = stage - 1
    
    if stage == 1:
        for z in inrow[1:]:
            nextrow.append(z)
        return nextrow
    
    else:
        for x in range(cap):
            nextrow.append(inrow[(stage-1) + (cap - x)])
            nextrow.append(inrow[(stage-1) - (x+1)])
        for y in inrow[(stage+cap):]:
            nextrow.append(y)
        return(nextrow)
            
def genNextLIRO(inrow, stage):
    nextrow = []
    cap = stage - 1
    
    if stage == 1:
        for z in inrow[1:]:
            nextrow.append(z)
        return nextrow
    
    else:
        for x in range(cap):
            nextrow.append(inrow[(stage-1) - (x+1)])
            nextrow.append(inrow[(stage-1) + (cap - x)])
        for y in inrow[(stage+cap):]:
            nextrow.append(y)
        return(nextrow)
    
def genNextLORI(inrow, stage):
    nextrow = []
    cap = stage - 1
    
    if stage == 1:
        for z in inrow[1:]:
            nextrow.append(z)
        return nextrow
    
    else:
        for x in range(cap):
            nextrow.append(inrow[(stage-1) - (cap - x)])
            nextrow.append(inrow[(stage-1) + (x+1)])
        for y in inrow[(stage+cap):]:
            nextrow.append(y)
        return(nextrow)

def genNextRILO(inrow, stage):
    nextrow = []
    cap = stage - 1
    
    if stage == 1:
        for z in inrow[1:]:
            nextrow.append(z)
        return nextrow
    
    else:
        for x in range(cap):
            nextrow.append(inrow[(stage-1) + (x+1)])
            nextrow.append(inrow[(stage-1) - (cap - x)])
        for y in inrow[(stage+cap):]:
            nextrow.append(y)
        return(nextrow)

# Maximum stage ~ 0.333347(numlim)+0.528058
def generate(numlim, stagelim, form):
    rows = []
    rows.append(list(np.linspace(1, numlim, num=numlim)))
    if form == 'RILI':
        for x in range(stagelim):
            rows.append(genNextRILI(rows[x], x+1))
    if form =='ROLO':
        for x in range(stagelim):
            rows.append(genNextROLO(rows[x], x+1))
    if form == 'LIRI':
        for x in range(stagelim):
            rows.append(genNextLIRI(rows[x], x+1))
    if form == 'LORO':
        for x in range(stagelim):
            rows.append(genNextLORO(rows[x], x+1))
    if form == 'ROLI':
        for x in range(stagelim):
            rows.append(genNextROLI(rows[x], x+1))
    if form == 'LIRO':
        for x in range(stagelim):
            rows.append(genNextLIRO(rows[x], x+1))
    if form == 'LORI':
        for x in range(stagelim):
            rows.append(genNextLORI(rows[x], x+1))
    if form == 'RILO':
        for x in range(stagelim):
            rows.append(genNextRILO(rows[x], x+1))
    return rows

def orbit(kimberling, value):
    places = []
    for i in kimberling:
        try:
            places.append(i.index(value) + 1)
        except ValueError:
            places.append(None)
    return places

def diagonal(kimberling):
    return np.diag(kimberling)

# from https://oeis.org/A035486
def K(i, j):
    if j >= 2*i-3: return i+j-1
    q, r = divmod(j+1, 2)
    return K(i-1, i-1+(1-2*r)*q)

# from https://oeis.org/wiki/User:Enrique_P%C3%A9rez_Herrero/Kimberling
def L(n):
    i = math.floor((n+4)/3)
    j = math.floor((2*n + 1)/3)
    while (i!=j):
        j = max(2*(i - j), 2*(j - i) - 1)
        i = i + 1
    return i

#------------------------------------------------------------------------------
print(K(25, 25))
print(L(2))


# start = time.time()
# rows = generate(10000,3334, 'RILI')
# middle= time.time()
# print(middle - start)
# df = pd.DataFrame(rows)
# df.columns += 1
# df.rows += 1
# df.to_csv('RILI_10000.csv')
# end = time.time()
# print(end - middle)
    
# df = pd.DataFrame(rows)
# df.columns += 1
# print(diagonal(df))
# print()
# print(orbit(rows, 17))
