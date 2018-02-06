from gf2_mul import *
def eval_poly(points,coefs,mul = gf832_mul):
    ret = [0]*len(points)
    for i in reversed(coefs):
        for j in range(len(points)):
            ret[j] = mul(ret[j], points[j])
            ret[j] ^= i
    return ret 

def square(points,mul):
    return [mul(p,p) for p in points]
    
def exp(points,n,mul):
    ret = [1]*len(points)
    for i in range(n):
        ret = [mul(ret[i],points[i]) for i in range(len(points))]
    return ret
        
def square_add(points,mul):
    return [mul(p,p)^p for p in points]

def int_to_bitvec(a, n):
  r = [((a>>i)&1) for i in range(n)]
  return r

def bitvec_to_int(a):
  r=0
  for i in range(len(a)):
    r+=(a[i]<<i)
  return r

#def transform_matrix(points, deg, gf_size, mul):
#    ret = []
#    exp = [1]*len(points)
#    for i in range(deg):
#        l=[] 
#        for i in exp:
#            l += int_to_bitvec(i, gf_size)
#        for i,v in enumerate(points):
#            exp[i] = mul(exp[i],v)
#        ret.append(l)
#    return ret

def transform_matrix(points, deg, gf_size = None, mul = gf832_mul):
    ret = [[] for i in range(deg)]
    for p in points:
        tmp = 1
        if gf_size is None:
            max_number_of_bits = 0
            for i in range(deg):
                max_number_of_bits = max(max_number_of_bits, tmp.bit_length())
                tmp = mul(tmp,p)
        else:
            max_number_of_bits = gf_size
        
        tmp = 1
        for i in range(deg):
            ret[i] += int_to_bitvec(tmp, max_number_of_bits)
            tmp = mul(tmp,p)

    return ret


def gf2_gauss_elim(origin_mat):
    mat = [list(vec) for vec in origin_mat]
    m = len(mat)
    n = len(mat[0])
    r = 0
    for i in range(n):
        if mat[r][i] == 0:
            for j in range(r+1,m):
                if mat[j][i] == 1:
                    mat[r],mat[j] = mat[j],mat[r]
        if mat[r][i] == 1:
            for j in range(m):
                if mat[j][i] == 1 and j!=r:
                    for k in range(i,n):
                        mat[j][k] ^= mat[r][k]
            r+=1
            if r==m:
                break
    return mat

def gf2_mat_show(mat):
    for v in mat:
        s = ""
        for i in v:
            s+=str(i)
        print(s)

def gf2_rank(mat):
    rank = 0
    mat = gf2_gauss_elim(mat)
    for v in mat:
        if all([i==0 for i in v]) != True:
            rank+=1
    return rank
    

def ran(s,n):
    return [s+i for i in range(n)]

def pseq(n):
    if n==0:
        return [0]
    if n==1:
        return [1]
    s = [1]
    j=0
    for i in range(2,n+1):
        if i-1==2**(j):
            j+=1
            s = [v*2 for v in s]
        else:
            s = [v*2 + a for v in s for a in [0,1]]
    return s

def fbasis(n):
    if n==0:
        return [0]
    if n==1:
        return [1]
    s = [1]
    j=0
    for i in range(2,n+1):
        if i-1==2**(j):
            j+=1
            s = [v*2 for v in s]
        else:
            s = [1] +  [v*2 for v in s ]
    return s


import random
def gf2_rand(n):
    return [random.randint(0,1) for i in range(n)]

def gf256_rand(n):
    return [random.randint(0,256) for i in range(n)]

def gf2_poly_mul(a,b):
    l = len(a)+len(b)-1
    ret = [0]*l
    for i in range(len(a)):
        if a[i]==1 :
            for j in range(len(b)):
                ret[i+j] = ret[i+j]^b[j]

    return ret

def gf2_poly_div(a,b):
    a=strip_zero(a)
    la = len(a)
    b=strip_zero(b)
    lb = len(b)
    assert lb>0
    count = 0
    rm = list(a)
    ret = []
    for i in range(la-1,lb-2,-1):
        if rm[i] == 1:
            for j in range(lb):
                rm[i-j] = rm[i-j]^b[lb-1-j]
            ret.append(1)
        else:
            ret.append(0)
    ret.reverse() 
    rm = strip_zero(rm)
    return ret,rm

def gf2_poly_div_count(a,b):
    a=strip_zero(a)
    la = len(a)
    b=strip_zero(b)
    lb = len(b)
    assert lb>0
    count = 0
    rm = list(a)
    ret = []
    for i in range(la-1,lb-2,-1):
        if rm[i] == 1:
            for j in range(lb):
                if b[lb-1-j]==1:
                    rm[i-j] = rm[i-j]^b[lb-1-j]
                    count+=1
            ret.append(1)
        else:
            ret.append(0)
    ret.reverse() 
    rm = strip_zero(rm)
    return ret,rm,count

def strip_zero(l):
    ret = list(l)
    while len(ret)>0 and ret[-1]==0:
        ret.pop(-1)
    return ret

def gf2_poly_add(a,b):
    l=max(len(a),len(b))
    ret = [0]*l
    for i in range(len(a)):
        ret[i] = ret[i] ^ a[i]
    for i in range(len(b)):
        ret[i] = ret[i] ^ b[i]
    return ret

def gf2_poly_test():
    a = gf2_rand(12)
    b = gf2_rand(8)
    print(a)
    print(b)
    q,r = gf2_poly_div(a,b)
    print(q,r)
    tmp = (strip_zero(gf2_poly_add(gf2_poly_mul(q,b),r)))
    return tmp == strip_zero(a)


def extend_zero(l,n):
    return l+[0]*(n-len(l))


def cantor_square(x):
    return x^(x>>1)

def cantor_square_n(x,n):
    """
    Example:
    [(bin(cantor_square_n(1<<63,i))) for i in range(64)]
    """
    for i in range(n):
        x=cantor_square(x)
    return x

