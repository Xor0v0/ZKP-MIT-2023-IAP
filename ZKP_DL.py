#!/usr/bin/env python
from Crypto.Util.number import getPrime
from random import random
import numpy as np
from functools import reduce
from math import floor

ROUND = 100

def qpow(x, exp):
    if exp == 0:
        return 1
    if exp % 2 == 0:
        half = qpow(x, exp // 2)
        return half * half
    else:
        return x * qpow(x, exp - 1)

# 生成100个随机bit
def pseudoRandomBits(h, g, p):
    seed = np.sum(h)
    hash = qpow(g, seed) % p

    b = np.array([])
    for i in range(ROUND):
        temp = hash & i
        bit = 0
        if temp > 0:
            bit = 1
        b = np.append(b, bit)
    return b

# x: secret, g: generator, p: big prime
def dlogProof(x, g, p):
    y = qpow(g, x) % p # 计算g^x (mod p)

    # 生成100个随机数r，满足0 <= r < p-1
    r, h= np.array([]), np.array([])
    for _ in range(ROUND):
        r = np.append(r, floor(random() * (p - 1)))
    for _ in r:
        h = np.append(h, qpow(g, _) % p)
    b = pseudoRandomBits(h, g, p) # fiat-shamir
    s = (r + b * x) % (p - 1)
    pf = (h, s)
    return (y, pf)

def verify(y, g, p, pf):
    (h, s) = pf
    b = pseudoRandomBits(h, g, p)
    lhs = (g ** s) % p
    rhs = (h * (y ** b)) % p
    res = reduce(lambda l,r: l and r, lhs == rhs)
    return res

if __name__ == '__main__':
    p = getPrime(4)
    g = getPrime(4)
    x = getPrime(4)

    (y, pf) = dlogProof(x, g, p)
    res = verify(y, g, p, pf)
    print("The proof is: ", res)
    res = verify(y - 1, g, p, pf)
    print("The proof is: ", res)




