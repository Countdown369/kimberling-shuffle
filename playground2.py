#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 15:52:33 2023

@author: cabrown802
"""

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

# Maximum stage ~ 0.333347(numlim)+0.528058
def generate(stagelim, form):
    rows = []
    rows.append(list(np.linspace(1, stagelim*3, num=stagelim*3)))
    if form == 'RILI':
        for x in range(stagelim):
            rows.append(genNextRILI(rows[x], x+1))
    else:
        pass
    return rows

def orbit(kimberling, value):
    places = []
    for i in kimberling:
        try:
            places.append(i.index(value) + 1)
        except ValueError:
            pass
    return places

def efficientOrbit(value, ultraSneaky):
    if ultraSneaky:
        start = time.time()
    orbit = []
    a_n = value
    n = 1
    while a_n != n:
        if ultraSneaky:
            if time.time() > start + 0.2:
                return [69,0]
        orbit.append(a_n)
        if a_n < n:
            a_n = 2*(n - a_n)
            n = n+1
            continue
        elif a_n > 2*n:
            a_n = a_n - 1
            n = n+1
            continue
        else:
            a_n = 2*(a_n - n) - 1
            n = n+1
            continue
    orbit.append(a_n)
    return orbit

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

def makeshiftFamilyChecker(kim, start, end, skipAnnoy):
    # for families seven and eight
    for i in range(start, end+1):
        if skipAnnoy:
            if i in [2, 6, 11, 12, 17, 19, 27, 29, 30, 32, 40, 43, 49, 50]:
                print(i)
                continue
        bro = orbit(kim, i)
        reltodiag = []
        for k in range(len(bro)):
            reltodiag.append(bro[k] - (k+1))
            
        signs = [math.copysign(1, j) for j in reltodiag]
        
        if signs[-3:] == [-1, -1, 1] and -1 not in signs[:(len(signs)-3)]:
            print(str(i) + '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
            print(bro)
            break
            
        else:
            print(i)
            
def nineFamilyChecker(kim, start, end, skipAnnoy):
    facts = []
    # for yet discovered family?
    for i in range(start, end+1):
        if skipAnnoy:
            if i in [2, 6, 11, 12, 17, 19, 27, 29, 30, 32, 40, 43, 49, 50]:
                print(i)
                continue
        bro = orbit(kim, i)
        reltodiag = []
        for k in range(len(bro)):
            reltodiag.append(bro[k] - (k+1))
            
        signs = [math.copysign(1, j) for j in reltodiag]
        
        if signs[-5:] == [1, -1, 1, -1, 1] and -1 not in signs[:(len(signs)-5)]:
            print(str(i) + '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
            print(bro)
            facts.append(i)
            
        else:
            print(i)
            
            
def fastFamilyChecker(start, end):
    # for families seven and eight
    for i in range(start, end+1):
        bro = efficientOrbit(i, True)
        reltodiag = []
        for k in range(len(bro)):
                reltodiag.append(bro[k] - (k+1))
                
        signs = [math.copysign(1, j) for j in reltodiag]
            
        if signs[-3:] == [-1, -1, 1] and -1 not in signs[:(len(signs)-3)]:
            print(str(i) + '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
            print(bro)
            return 0
                
        else:
            print(i)
            
def fastNineFamilyChecker(start, end):
    facts = []
    # for yet discovered family?
    for i in range(start, end+1):
        bro = efficientOrbit(i, True)
        reltodiag = []
        for k in range(len(bro)):
            reltodiag.append(bro[k] - (k+1))
            
        signs = [math.copysign(1, j) for j in reltodiag]
        
        if signs[-5:] == [1, -1, 1, -1, 1] and -1 not in signs[:(len(signs)-5)]:
            print(str(i) + '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
            print(bro)
            facts.append(i)
            print(facts)
            
        else:
            print(i)
            
        if i % 20 == 0:
            print(facts)
        
    print(facts)
    return facts

def twoFastNineFamilyChecker(start, middle, end):
    facts = []
    # for yet discovered family?
    for i in range(start, end+1):
        # I changed this to False  
        bro = efficientOrbit(i, True)
        reltodiag = []
        for k in range(len(bro)):
            reltodiag.append(bro[k] - (k+1))
            
        signs = [math.copysign(1, j) for j in reltodiag]
        
        if signs[-4:] == [-1, -1, -1, 1] and -1 not in signs[:(len(signs)-4)]:
            print(str(i) + '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
            print(bro)
            facts.append(i)
            print(facts)
            
        else:
            print(i)
            
        if i % 100 == 0:
            print(facts)
        
    print(facts)
    return facts

def volcano(start, middle, end):
    volcano_list = [middle]
    the_range = min(end-middle, middle-start) // 3
    for i in range(1, the_range + 1):
        volcano_list.append(middle-(i*3))
        volcano_list.append(middle+(i*3))
    for i in range(volcano_list[-1]+3, end+1, 3):
        volcano_list.append(i)
    return(volcano_list)
            
def diffChecker(arr):
  result = []
  for i in range(len(arr)-1):
    result.append(arr[i+1]-arr[i])
  return result

def fracCheck():
    for i in range(1, 100):
        for j in range (1, i):
            if (j/i) > 0.6923 and (j/i) < 0.69232:
                print(j/i)
                print(str(j) + '/' + str(i))


'''
LOOK HERE!
fastNineFamilyChecker(96210, 1000000)

myboy = generate(7800, 'RILI')

[-2, -1, 4]  (7 = 1 (mod 3)) {1}
[-2, -6, -14, -30, -62, -25, 52] (283 = 1 (mod 3)) {5}
[-2, -6, -14, -30, -62, -126, -254, -510, -1022, -409, 820] (4879 = 1 (mod 3)) {9}
Hypothetical next:
[-2, -6, -14, -30, -62, -126, -254, -510, -1022, ...] (????? = 1 (mod 3)) {13}

YES!
[-2, -6, -14, -30, -62, -126, -254, -510, -1022, -2046, -4094, -8190, -16382, -6553, 13108] (78595 = 1 (mod 3)) {13}

[-4, -10, -22, -9, 20] (95 = 2 (mod 3)) {3}
[-4, -10, -22, -46, -94, -190, -382, -153, 308] (1811 = 2 (mod 3)) {7}
Hypothetical next:
[-4, -10, -22, -46, -94, -190, -382, ...] 29447 = 2 (mod 3)) {11}
'''
'''
# [-4, 7, -5, 5] (8 = 2 (mod 3)) {1}
# [oops] (90692 = 2 (mod 3)) {13}
# k inc. by 12???
# YEP! (371679776 = 2 (mod 3)) {25}

# [-2,-6, -14, -30, 35, -29, 21] (97 = 1 (mod 3)) {3}
# Next at {15}?
# Yepp...483901 = 1 (mod 3) {15}

# [-3, -8, -18, -38, -78, -158, -318, -638, -1278, 1379, -1181, 789] (4689 = 0 (mod 3)) {8}
# Next at {20}?
'''

# [28, 402, 2698, 27264]
# 1 mod 3
# 28 {1}
# 2698 {7}
#
# 0 mod 3
# 402 {4}
# 27264 {10}

# [382, 6514] 26667
# {5}, {9}

# start = time.time()
# efficientOrbit(1747566, False)
# end = time.time()
# print(end-start)