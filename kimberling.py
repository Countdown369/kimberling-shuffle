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
        if x % 2 == 0:
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
            if x % 2 == 1:
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
        if x % 2 == 0:
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
            if x % 2 == 1:
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
def generate(stagelim, form):
    rows = []
    rows.append(list(np.linspace(1, stagelim*3, num=stagelim*3)))
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
    places2 = []
    for k in places:
        if k != None:
            places2.append(k)
    return places2

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

# print(K(25, 25))
# for i in range(1,10):
#     print(L(i))
#     print(orbit(generate(25, 'RILI'),i))

# RED BROWN YELLOW!
# for i in range(1,10):
#     print([3*i - 4, 3*i - 3, 3*i - 2])


"""

# print(orbit(generate(10000, 3334, 'RILI'),20))
# print(orbit(generate(10000, 3334, 'RILI'),53))

# These numbers, when do they get expelled? They have a long run of decreasing
# by 1, and then very soon after (how soon after?), they are expelled.

# print(orbit(generate(10000, 3334, 'RILI'),122))
# print(generate(30, 5, 'LORO'))
"""
# start=time.time()
# thething = orbit(generate(49600, 'RILI'),19)
# print(max(thething))
# print(thething.index(max(thething)))
# end=time.time()
# print(end - start)

kimberling1 = generate(3100, 'RILI')

fifteen_orbit = orbit(kimberling1, 15)
fifteen_res = [fifteen_orbit[i + 1] - fifteen_orbit[i] for i in range(len(fifteen_orbit)-1)]
fortytwo_orbit = orbit(kimberling1, 42)
fortytwo_res = [fortytwo_orbit[i + 1] - fortytwo_orbit[i] for i in range(len(fortytwo_orbit)-1)]
ninetynine_orbit = orbit(kimberling1, 99)
ninetynine_res = [ninetynine_orbit[i + 1] - ninetynine_orbit[i] for i in range(len(ninetynine_orbit)-1)]
twosixteen_orbit = orbit(kimberling1, 216)
twosixteen_res = [twosixteen_orbit[i + 1] - twosixteen_orbit[i] for i in range(len(twosixteen_orbit)-1)]
fourfiftythree_orbit = orbit(kimberling1, 453)
fourfiftythree_res = [fourfiftythree_orbit[i + 1] - fourfiftythree_orbit[i] for i in range(len(fourfiftythree_orbit)-1)]
print(15)
print(fifteen_res)
print()
print(42)
print(fortytwo_res)
print()
print(99)
print(ninetynine_res)
print()
print(216)
print(twosixteen_res)
print()
print(453)
print(fourfiftythree_res)
print()

#UNCOIMMENT THIS AFTER TALKING TO DAVID
#FAMILY 2: 15*2^k - 3k - 12
for i in range(456,1002, 3):
    print(i)
    epic = orbit(kimberling1, i)
    epic_res = [epic[i + 1] - epic[i] for i in range(len(epic)-1)]
    if epic_res[len(epic_res)-1] == -158:
        epic2 = orbit(kimberling1, i+3)
        epic_res2 = [epic2[i + 1] - epic2[i] for i in range(len(epic2)-1)]
        print(epic_res)
        print()
        print(i+3)
        print(epic_res2)
        
        break
    else:
        print("Failed")
        print()


def perfectlyExpelledFamily():
    # Family 9*2^k - 3k - 10
    
    for k in range(1, 10):
        print(9 * (2**k) - (3*k) - 10)
        # print(orbit(generate(1000, 334, 'RILI'),9 * (2**k) - (3*k) - 10))
        # print(orbit(generate(1000, 334, 'RILI'),k))
    
    # How many times does the number go down by only one row?
    [2, 7, 18, 41, 88, 183]
    # A[095151 | 077802]!!!
    
    # Outside of the family? Starts at 1, like
  # [1, 2, 3, 4, 5, 6, 7, 8, 9, ...]
    [0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, ...]
    # floor((1/3)*(row + 1.5))
    
    # A[095151 | 077802] follows this ^ function, where row is 9*2^k - 3k - 10.
    
    # Minimum row expelled, according to Herrero? (Starts at 1)
    [1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, ...]
    # So, the above list plus one. Makes sense!
    # Make sure to be able to explain why, but this makes "intuitive" sense:
    # After decreasing by one, getting attracted towards the expelled value,
    # if "expelled perfectly", it will be expelled on row
    # floor((1/3)*(row + 1.5)) + 1. And Herrero's formula is equivalent!


# print(generate(4, 'RILI'))

# ceil to (x+2)/2

# start = time.time()
# rows = generate(5, 'RILI')
# middle= time.time()
# for i in range(len(rows)):
#     print(rows[i][i])
"""
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
"""

# epic = generate(10, 4, 'LORO')
# print(epic)
# perfectlyExpelledFamily()

# for i in range(len(epic)):
#     print(epic[i][0])

# in class - bound by triangular numbers (quadratic) 
# bound? 2n + 1 + n for expelled number for row n
