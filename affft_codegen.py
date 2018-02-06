from symbol import *
import sys

def affft_codegen(coeffs, k, a, eqs, mul = gf832_mul, n=0):
    if k==-1:
        return [c for c in coeffs ]
    w = s(k,a,mul)

    p0 = list(coeffs[0:len(coeffs)//2])
    p1 = list(coeffs[len(coeffs)//2::])
    """
    h0=p0 ^ p1*w
    h1=h0 ^ p1
    """

    bit_length = 2**((n).bit_length())
    h0=[]
    h1=[]
    for c0,c1 in zip(p0,p1):
        if 2**((n).bit_length()-1) == n:
            mul_exp = Mul_Const(c1,w)
            mul_var = Var(bit_length)
            eqs.append((mul_var,mul_exp))

            add0 = Add(c0,mul_var)
            add_var0 =  Var(bit_length)
            eqs.append((add_var0,add0))
            h0.append(add_var0)
            continue
        butterfly_out = Var(bit_length*2)

        eqs.append((butterfly_out,Butterfly(c0,c1,w)))
        hh0 = Var(bit_length)
        hh1 = Var(bit_length)

        h0.append(hh0)
        h1.append(hh1)
        eqs.append((hh0,LowerHalf(butterfly_out)))
        eqs.append((hh1,UpperHalf(butterfly_out)))
        
    #do encode
    if 2**((n).bit_length()-1) == n:
        return affft_codegen(h0,k-1,a,eqs,mul,n+1)
    
    if n==0:
        return [affft_codegen(h0,k-1,a,eqs,mul,n), affft_codegen(h1,k-1,a^(1<<(k)),eqs,mul,n+1)]

    #assert len(coeffs) == 2**(k+1)

    return affft_codegen(h0,k-1,a,eqs,mul,n+1) + affft_codegen(h1,k-1,a^(1<<k),eqs,mul,n+1)


def iaffft_codegen(points, k, a,eqs, mul = gf832_mul, inv = gf832_inv, n=0):
    if k==-1:
        return points
        return [c for c in points ]

    w = s(k,a,mul)

    p1 = []
    p0 = []
    if n==0:
        h0 = iaffft_codegen(points[0],k-1,a,eqs,mul,inv,n)
        h1 = iaffft_codegen(points[1],k-1,a^(1<<k),eqs,mul,inv,n+1)
        for c0,c1 in zip(h0,h1):
            add0 = Add(c0,c1)
            add_var0 =  Var(max(len(c0.bits),len(c1.bits)))
            eqs.append((add_var0,add0))
            p1.append(add_var0)

            mul_exp = Mul_Const(add_var0,w)
            mul_var = Var(len(add_var0.bits))
            eqs.append((mul_var,mul_exp))

            add1 = Add(c0,mul_var)
            add_var1 = Var(max(len(c0.bits),len(c1.bits)))
            p0.append(add_var1)
            eqs.append((add_var1,add1))
        return p0+p1

    #do decode
    if 2**(n.bit_length()-1) == n:
        h1 = iaffft_codegen(points,k-1,a,eqs,mul,inv,n+1)

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
            #u = aa & (2**n-1)
            var_u = Var(len(aa.bits)//2)
            eqs.append((var_u,LowerHalf(aa)))
            #v = aa >> n
            var_v = Var(len(aa.bits)//2)
            eqs.append((var_v,UpperHalf(aa)))
            #b = mul(v,inv(c1))
            #b = v
            mul_var = Var(len(var_v.bits))
            eqs.append((mul_var,Mul_Const(var_v,c0)))
            var_a = Var(len(var_u.bits))
            eqs.append((var_a,Add(mul_var,var_u)))
            #a = mul(b,c0) ^ u
            p0.append(var_a)
            p1.append(var_v)
            
        return p0+p1

    h0 = iaffft_codegen(points[0:len(points)//2],k-1,a,eqs,mul,inv,n+1)
    h1 = iaffft_codegen(points[len(points)//2::],k-1,a^(1<<k),eqs,mul,inv,n+1)
    for c0,c1 in zip(h0,h1):
        add0 = Add(c0,c1)
        add_var0 =  Var(max(len(c0.bits),len(c1.bits)))
        eqs.append((add_var0,add0))
        p1.append(add_var0)

        mul_exp = Mul_Const(add_var0,w)
        mul_var = Var(len(add_var0.bits))
        eqs.append((mul_var,mul_exp))

        add1 = Add(c0,mul_var)
        add_var1 = Var(max(len(c0.bits),len(c1.bits)))
        p0.append(add_var1)
        eqs.append((add_var1,add1))
    return p0 + p1




counter = 0
def xors(v, eqs, start, length, direction, offset, size, distance):
    for i in range(start, start+length*direction, direction):
        for j in range(size):
            #print(offset, size, distance)
            if eqs == None:
                v[offset + size*(i-distance)+j] ^= v[offset + size*i + j]
                global counter
                counter+=1
            else:
                yy = v[offset + size*(i-distance)+j] 
                xx = v[offset + size*i + j]
                eqs.append( (yy,Add(xx,yy)) )
    

def substitute(v, eqs, n, offset, size, m):
    """
    x^{2^m}+x
    """
    a = n-m-1
    if a<0:
        return
    dist = (2**m-1)*2**a
    xors(v, eqs, 2**n-1, 2**(n-1), -1, offset, size, dist)
    substitute(v, eqs, n-1, offset, size, m)
    substitute(v, eqs, n-1, offset+2**(n-1)*size, size, m)


def bc(v, eqs, n, offset=0, size=1):
    if n==0:
        return v
    if n==1:
        return v
    d = 0
    while 2**(d+1) <n:
        d=d+1
    m = 2**d
    """
    x^{2^m}+x
    """
    substitute(v, eqs, n, offset, size, m)
    bc(v, eqs, n-m, offset, (size*2**m))
    for i in range(2**(n-m)):
        bc(v, eqs, m, offset+i*size*2**m, size)

    return v

def isubstitute(v, eqs, n, offset, size, m):
    """
    x^{2^m}+x
    """
    a = n-m-1
    if a<0:
        return
    dist = (2**m-1)*2**a
    isubstitute(v, eqs, n-1, offset+2**(n-1)*size, size, m)
    isubstitute(v, eqs, n-1, offset, size, m)
    xors(v, eqs, 2**(n-1), 2**(n-1), +1, offset, size, dist)



def ibc(v, eqs, n, offset=0, size=1):
    if n==0:
        return v
    if n==1:
        return v
    d = 0
    while 2**(d+1) <n:
        d=d+1
    m = 2**d
    """
    x^{2^m}+x
    """
    for i in range(2**(n-m)):
        ibc(v, eqs, m, offset+i*size*2**m, size)
    ibc(v, eqs, n-m, offset, (size*2**m))
    isubstitute(v, eqs, n, offset, size, m)

    return v


def dfs_var_mul(a,b,eqs):
    ret = []
    for i,j in zip(a,b):
        if isinstance(i,list):
            ret.append(dfs_var_mul(i,j,eqs))
        else:
            v = Var(len(i.bits))
            eqs.append((v,Mul(i,j)))
            ret.append(v)
    return ret

code = ""

#print(sys.argv[1])
k=int(sys.argv[1])

var_in0 = [Var(1) for i in range(2**(k+1))]
var_in1 = [Var(1) for i in range(2**(k+1))]
eqs_bc0 = []
bc(var_in0, eqs_bc0, k)
eqs0 = []
var_out0=(affft_codegen(var_in0,k,0,eqs0))

eqs_bc1 = []
bc(var_in1, eqs_bc0, k)
eqs1 = []
var_out1=(affft_codegen(var_in1,k,0,eqs1))

#TODO Point multiplication here
var_point_mul_out = var_out0
eqs2=[]
var_point_mul_out = dfs_var_mul(var_out0,var_out1,eqs2)

eqs3 = []
var_out2=(iaffft_codegen(var_point_mul_out,k,0,eqs3))
eqs_ibc = []
ibc(var_out2, eqs_ibc, k+1)
#print(k+1)

field_in0 = gf2_rand(2**(k))
field_in1 = gf2_rand(2**(k))

#init input value

code_in = "#input1\n"
code_zero = ""

for i in range(2**(k)):
    code_in += str(var_in0[i].bits[0]) + " = " +  str(field_in0[i]) + "\n"
code_in += "#input2\n"
for i in range(2**(k)):
    code_in += str(var_in1[i].bits[0]) + " = " +  str(field_in1[i]) + "\n"
for i in range(2**(k)):
    code_zero += str(var_in1[i+2**k].bits[0]) + " = " +  "0;\n"
for i in range(2**(k)):
    code_zero += str(var_in0[i+2**k].bits[0]) + " = " +  "0;\n"

code += eq_to_str(eqs_bc0 + eqs0+ eqs_bc1 + eqs1+eqs2+eqs3 + eqs_ibc)

stat={1<<i:0 for i in range(k+1)}
def dfs_mul(a,b):
    ret = []
    for i,j in zip(a,b):
        if isinstance(i,list):
            ret.append(dfs_mul(i,j))
        else:
            ret.append(gf832_mul(i,j))
    return ret



    
def dfs_sub(vs):
    ret = []
    for v in vs:
        if isinstance(v,list):
            ret.append(dfs_sub(v))
        else:
            stat[len(v.bits)]+=1
            ret.append( "+".join(["{}*(1<<{})".format(b,i) for i,b in enumerate(v.bits)]))
    return "["+",".join(ret)+"]"


out0 = affft(basis_conversion(field_in0+[0]*2**k, k),k,0)
out1 = affft(basis_conversion(field_in1+[0]*2**k, k),k,0)

point_mul = dfs_mul(out0,out1)
point_mul_iaffft_out = iaffft(point_mul,k,0)

result_via_affft_val = list(point_mul_iaffft_out)
ibc(result_via_affft_val,None,k+1)
result_via_affft_val = ibasis_conversion(iaffft(point_mul,k,0),k)

code_out=""

#code_out+=("print({})\n".format(dfs_sub(var_out0)))
#code_out+=("print({})\n".format( str(out0)))
#
#code_out+=("print({})\n".format(dfs_sub(var_out1)))
#code_out+=("print({})\n".format( str(out1)))
#
#code_out+=("print('========')\n")
#code_out+=("print({})\n".format(dfs_sub(var_point_mul_out)))
#code_out+=("print({})\n".format(point_mul))
#code+=("print({})\n".format( str(affft(field_in1+[0]*2**k,k,0)) ))

#code_out+=("print({})\n".format(dfs_sub(var_out2)))
#code_out+=("print({})\n".format(result_via_affft_val ))

result_via_direct_mul = gf2_poly_mul(field_in0,field_in1)+ [0]
#code_out+=("print({})\n".format(result_via_direct_mul ))
#code_out+=("print({} == {})\n".format(dfs_sub(var_out2),( result_via_affft_val )))

code_out+="#output\n"
code_out+=("{}\n".format([i.bits[0] for i in var_out2]))

#code_out+=("print({} == {})\n".format([i.bits[0] for i in var_out2], result_via_direct_mul))
#code_out+=("print({} == {})\n".format(result_via_affft_val, result_via_direct_mul))



code_prune = ""
code = code.replace(";","")
code = code_zero + code
zeros=set()
for l in code.split("\n"):
    l=l.strip()
    t = l.split()
    if len(t) == 0 :
        continue
    if len(t)!=0 and t[-1] == "0;":
        zeros.add(t[0])
    t = l.split(" = ")
    if len(t)==1:
        print(l)
        continue
    tmp = t[1].split(" ^ ")
    tmp2 = []
    for i in tmp:
        if i not in zeros and i!="0":
            tmp2.append(i)
    if len(tmp2)==0:
        code_prune+=("{} = 0\n".format(t[0]))
        zeros.add(t[0])
        continue
        
    if len(tmp2)==1 and tmp2[0] in zeros:
        zeros.add(t[0])
    code_prune+=("{} = {}\n".format(t[0], " ^ ".join(tmp2)))

print_code = False
print_code = True

if(print_code):
    #print(code_in)
    print("#input1")
    print([i.bits[0] for i in var_in0[0:2**k]])
    print("#input2")
    print([i.bits[0] for i in var_in1[0:2**k]])
    print(code_prune)
    print(code_out)
#code_prune = code

#exec(code_in+code_prune+code_out)
#print(2**(k))
bc_bt_pm_xor = code_prune.count("^")
pm_and = code_prune.count("&")
#print("# bc+butterfly+point_mul ^",bc_bt_pm_xor)
#print("# point_mul &",pm_and)

#bc_xor = rec(k-1)*2**k/2 *2 + rec(k)*2**(k+1)/2
#print("# basis conversion ^", bc_xor)
#print("# total",bt_pm_xor+pm_and+bc_xor)
#print("# total",bc_bt_pm_xor+pm_and)

    
