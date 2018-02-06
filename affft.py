from gf2_poly import *

def s(k,x,mul=gf832_mul):
    if k==0:
        return x
    tmp = s(k-1,x,mul)
    return mul(tmp,tmp)^ tmp


def fft_constants(k,a):
    if k==0:
        return [[s(k,a)]]
    l1 = fft_constants(k-1,a)
    l2 = fft_constants(k-1,a^(1<<k))
    return [[s(k,a)]] + [ l1[i] + l2[i] for i in range(k)]

def fft(coeffs, k, a, mul = gf832_mul):
    if k==-1:
        return [c for c in coeffs ]
    assert len(coeffs) == 2**(k+1)
    w = s(k,a,mul)

    p0 = list(coeffs[0:len(coeffs)//2])
    p1 = list(coeffs[len(coeffs)//2::])
    h0 = [ c0 ^ mul(c1,w) for c0,c1 in zip(p0,p1)]
    h1 = [ c0 ^ c1 for c0,c1 in zip(h0,p1)]

    return fft(h0,k-1,a,mul) + fft(h1,k-1,a^(1<<k),mul)


"""
h0 = p0+w p1
h1 = p0+(w+1)p1

p1 = h0 + h1
p0 = h0 + w p1
"""

def ifft(points, k, a, mul = gf832_mul):
    if k==-1:
        return [c for c in points]
    w = s(k,a,mul)

    h0 = ifft(list(points[0:len(points)//2]),k-1,a,mul)
    h1 = ifft(list(points[len(points)//2::]),k-1,a^(1<<k),mul)
    p1 = [ c0 ^ c1 for c0,c1 in zip(h0,h1)]
    p0 = [ c0 ^ mul(c1,w) for c0,c1 in zip(h0,p1)]

    return p0 + p1



def affft(coeffs, k, a, mul = gf832_mul, n=0):
    if k==-1:
        return [c for c in coeffs ]
    w = s(k,a,mul)

    p0 = list(coeffs[0:len(coeffs)//2])
    p1 = list(coeffs[len(coeffs)//2::])
    h0 = [ c0 ^ mul(c1,w) for c0,c1 in zip(p0,p1)]
    h1 = [ c0 ^ c1 for c0,c1 in zip(h0,p1)]

    if n==0:
        return [affft(h0,k-1,a,mul,n) , affft(h1,k-1,a^(1<<(k)),mul,n+1)]

    #assert len(coeffs) == 2**(k+1)
    #do encode
    if 2**((n).bit_length()-1) == n:
        return affft(h0,k-1,a,mul,n+1)

    return affft(h0,k-1,a,mul,n+1) + affft(h1,k-1,a^(1<<k),mul,n+1)

def iaffft(points, k, a, mul = gf832_mul, inv = gf832_inv, n=0):
    if k==-1:
        return points
        return [c for c in points ]

    w = s(k,a,mul)

    if n==0:
        h0 = iaffft(points[0],k-1,a,mul,inv,n)
        h1 = iaffft(points[1],k-1,a^(1<<k),mul,inv,n+1)
        p1 = [ c0 ^ c1 for c0,c1 in zip(h0,h1)]
        p0 = [ c0 ^ mul(c1,w) for c0,c1 in zip(h0,p1)]
        return p0+p1

    #do decode
    if 2**(n.bit_length()-1) == n:
        h1 = iaffft(points,k-1,a,mul,inv,n+1)

        c1 = (w) >> n
        c0 = (w) & (2**n-1)
        
        """
        aa = a + b*(c1 x + c0) = (a+b*c0) + (b*c1) *x
                                    u         v
        b = v * c1^(-1)
        a = u + b*c0
        """
        p0=[]
        p1=[]
        for aa in h1:
            u = aa & (2**n-1)
            v = aa >> n
            #b = mul(v,inv(c1))
            b = v
            a = mul(b,c0) ^ u
            p0.append(a)
            p1.append(b)
            
        return p0+p1

    h0 = iaffft(points[0:len(points)//2],k-1,a,mul,inv,n+1)
    h1 = iaffft(points[len(points)//2::],k-1,a^(1<<k),mul,inv,n+1)
    p1 = [ c0 ^ c1 for c0,c1 in zip(h0,h1)]
    p0 = [ c0 ^ mul(c1,w) for c0,c1 in zip(h0,p1)]

    return p0 + p1





def affft2(coeffs, k, a, mul = gf832_mul, n=0):
    if k==-1:
        return [c for c in coeffs ]
    w = s(k,a,mul)

    if n==0:
        return affft2(coeffs,k,a^(1<<(k+1)),mul,n+1)

    assert len(coeffs) == 2**(k+1)

    p0 = list(coeffs[0:len(coeffs)//2])
    p1 = list(coeffs[len(coeffs)//2::])
    h0 = [ c0 ^ mul(c1,w) for c0,c1 in zip(p0,p1)]
    h1 = [ c0 ^ c1 for c0,c1 in zip(h0,p1)]

    #do encode
    if 2**((n).bit_length()-1) == n:
        return affft2(h1,k-1,a^(1<<k),mul,n+1)

    return affft2(h0,k-1,a,mul,n+1) + affft2(h1,k-1,a^(1<<k),mul,n+1)

def iaffft2(points, k, a, mul = gf832_mul, inv = gf832_inv, n=0):
    if k==-1:
        return [c for c in points ]

    w = s(k,a,mul)

    if n==0:
        return iaffft2(points,k,a^(1<<(k+1)),mul,inv,n+1)

    #do decode
    if 2**(n.bit_length()-1) == n:
        h1 = iaffft2(points,k-1,a^(1<<k),mul,inv,n+1)

        c1 = (w^1) >> n
        c0 = (w^1) & (2**n-1)
        
        """
        aa = a + b*(c1 x + c0) = (a+b*c0) + (b*c1) *x
                                    u         v
        b = v * c1^(-1)
        a = u + b*c0
        """
        p0=[]
        p1=[]
        for aa in h1:
            u = aa & (2**n-1)
            v = aa >> n
            #b = mul(v,inv(c1))
            b = v
            a = mul(b,c0) ^ u
            p0.append(a)
            p1.append(b)
            
        return p0+p1

    h0 = iaffft2(list(points[0:len(points)//2]),k-1,a,mul,inv,n+1)
    h1 = iaffft2(list(points[len(points)//2::]),k-1,a^(1<<k),mul,inv,n+1)
    p1 = [ c0 ^ c1 for c0,c1 in zip(h0,h1)]
    p0 = [ c0 ^ mul(c1,w) for c0,c1 in zip(h0,p1)]

    return p0 + p1




def s_poly(k):
    if k==0:
        return [0,1]
    else:
        p = s_poly(k-1)
        return gf2_poly_add(gf2_poly_mul(p,p),p)

def sparse_s_poly(k):
    if k==0:
        return [1]
    else:
        p = sparse_s_poly(k-1)
        #sparse mul
        ret = set()
        q=list(p)
        if 0 in q:
            q.remove(0)
        else:
            q.insert(0,0)
        
        for i in p:
            for j in q:
                if i+j not in ret:
                    ret.add(i+j)
                else:
                    ret.remove(i+j)
                    
        return sorted(list(ret))

def basis_conversion(coeffs,k):
    if k==0:
        return coeffs
    p1,p0 = gf2_poly_div(coeffs,s_poly(k))
    return basis_conversion(extend_zero(p0,2**k),k-1) + basis_conversion(extend_zero(p1,2**k),k-1)

def basis_conversion_count(coeffs,k):
    if k==0:
        return coeffs,0
    p1,p0,count = gf2_poly_div_count(coeffs,s_poly(k))
    b0 = basis_conversion_count(extend_zero(p0,2**k),k-1)
    b1 = basis_conversion_count(extend_zero(p1,2**k),k-1)
    return b0[0]+b1[0], b0[1]+b1[1]+count

def ibasis_conversion(coeffs,k):
    if k==0:
        return coeffs
    p0 = ibasis_conversion(list(coeffs[0:len(coeffs)//2]),k-1)
    p1 = ibasis_conversion(list(coeffs[len(coeffs)//2::]),k-1)
    return gf2_poly_add(gf2_poly_mul(p1,s_poly(k)),p0)
    
def test_basis_conversion(k):
    n=2**(k+1)
    a=gf2_rand(n)
    b=ibasis_conversion(basis_conversion(a,k),k)
    print(a)
    print(b)

def test_fft(k):
    n=2**(k+1)
    a=gf256_rand(n)
    fa = fft(a,k,0)
    ifa = ifft(fa,k,0)
    print(a)
    print(ifa)

def test_affft(k):
    n=2**(k+1)
    a=gf2_rand(n)
    fa = affft(a,k,0)
    print(fa)
    ifa = iaffft(fa,k,0)
    print(a)
    print(ifa)


def poly_mul_by_fft(a,b,k,mul=gf832_mul):
    aa = basis_conversion(a,k)
    bb = basis_conversion(b,k)
    faa = fft(aa,k,0,mul)
    fbb = fft(bb,k,0,mul)
    faabb = point_mul(faa,fbb,mul)
    aabb = ifft(faabb,k,0,mul)
    ab = ibasis_conversion(aabb,k)
    return ab

def poly_mul_by_affft(a,b,k,mul=gf832_mul):
    aa = basis_conversion(a,k)
    bb = basis_conversion(b,k)
    faa = affft(aa,k,0,mul)
    fbb = affft(bb,k,0,mul)
    faabb = point_mul(faa,fbb,mul)
    print("point mul result",faabb)
    aabb = iaffft(faabb,k,0,mul)
    ab = ibasis_conversion(aabb,k)
    return ab

def poly_mul_by_affft2(a,b,k,mul=gf832_mul):
    aa = basis_conversion(a,k)
    bb = basis_conversion(b,k)
    faa = affft2(aa,k,0,mul)
    fbb = affft2(bb,k,0,mul)
    faabb = point_mul(faa,fbb,mul)
    print("point mul result",faabb)
    aabb = iaffft2(faabb,k,0,mul)
    ab = ibasis_conversion(aabb,k)
    return ab
    
def point_mul(a,b,mul=gf832_mul):
    if isinstance(a, list):
        return [point_mul(u,v,mul) for u,v in zip(a,b)]
    else:
        return mul(a,b)
        

def test_poly_mul_fft(k):
    n=2**k
    a=extend_zero(gf2_rand(n),2*n)
    b=extend_zero(gf2_rand(n),2*n)
    print(a,b)
    result = poly_mul_by_fft(a,b,k)
    print (result)
    ref_result = gf2_poly_mul(a,b)
    ref_result = ref_result[0:len(ref_result)//2+1]
    print (ref_result)



def test_poly_mul_affft(k):
    n=2**k
    a=extend_zero(gf2_rand(n),2*n)
    b=extend_zero(gf2_rand(n),2*n)
    print(a,b)
    result = poly_mul_by_affft(a,b,k)
    print ("    ",result)
    ref_result = gf2_poly_mul(a,b)
    ref_result = ref_result[0:len(ref_result)//2+1]
    #print("ref pmul",affft(basis_conversion(ref_result,k),k,0))
    print ("ref ",ref_result)
    return ref_result == result


def test_poly_mul_affft2(k):
    n=2**k
    a=extend_zero(gf2_rand(n),2*n)
    b=extend_zero(gf2_rand(n),2*n)
    print(a,b)
    result = poly_mul_by_affft2(a,b,k)
    print ("    ",result)
    ref_result = gf2_poly_mul(a,b)
    ref_result = ref_result[0:len(ref_result)//2+1]
    #print("ref pmul",affft(basis_conversion(ref_result,k),k,0))
    print ("ref ",ref_result)
    return ref_result == result

def ft(coeffs,k):
    return eval_poly([i for i in range(2**(k+1))],coeffs)


def test_fft_bc(k):
    a = gf2_rand(2**(k+1))
    print(a)
    aa = basis_conversion(a,k)
    faa = fft(aa,k,0)
    fa = ft(a,k)
    print(faa)
    print(fa)
    print(ibasis_conversion((ifft(faa,k,0)),k))



def test_affft_bc(k):
    a = gf2_rand(2**(k+1))
    print(a)
    aa = basis_conversion(a,k)
    faa = affft(aa,k,0)
    print(faa)
    print(ibasis_conversion((iaffft(faa,k,0)),k))

def test_affft_fft(k):
    a = gf2_rand(2**(k+1))
    print(a)
    aa = basis_conversion(a,k)
    fa = ft(a,k)
    faa = fft(aa,k,0)
    ffaa = affft(aa,k,0)
    print(fa)
    print(faa)
    print(ffaa)


dp = dict()
# 2**(n+1) points
def rec(n):
    if n in dp:
        return dp[n]
    if n==0:
        dp[n] = 0
        return 0
    if n==1:
        dp[n] = 1
        return 1
    m = 0
    while 2**(m+1) <=n:
        m=m+1
    ret =  rec(n-2**m)+rec(2**m-1)+n-2**m+1
    dp[n] = ret
    return ret


def nbc(n):
    if n==1:
        return 0
    if n==0:
        return 0
    d = 0
    while 2**(d+1) <n:
        d=d+1
    m = 2**d
    """
    x^{2^m}+x
    """
    ret =  nbc(n-m)+nbc(m)+n-m
    return ret



