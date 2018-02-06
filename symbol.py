from affft import *
from gen_xor_equation import *

class Bit:
    counter = 0

    def __init__(self):
        self.id = self.new_id()

    @classmethod
    def new_id(cls):
        id = cls.counter
        cls.counter += 1
        return id

    def __repr__(self):
        return "Bit_"+str(self.id)

class Var:
    counter = 0

    def __init__(self, num_bit):
        self.bits = [Bit() for i in range(num_bit)]
        self.id = self.new_id()
    def __repr__(self):
        return "Var_"+str(self.id)

    @classmethod
    def new_id(cls):
        id = cls.counter
        cls.counter += 1
        return id

class Exp:
    def __init__(self):
        self.traverse = None
    
class Mul(Exp):
    def __init__(self, var0, var1):
        super().__init__()
        self.var0 = var0
        self.var1 = var1
    
    def __repr__(self):
        return str(self.var0) + " * " + str(self.var1)

class Mul_Const(Exp):
    def __init__(self, var0, var1):
        super().__init__()
        self.var0 = var0
        self.var1 = var1
    
    def __repr__(self):
        return str(self.var0) + " * " + str(self.var1)

class Add(Exp):
    def __init__(self, var0, var1):
        super().__init__()
        self.var0 = var0
        self.var1 = var1

    def __repr__(self):
        return str(self.var0) + " + " + str(self.var1)

class LowerHalf(Exp):
    def __init__(self, var0):
        super().__init__()
        self.var0 = var0

    def __repr__(self):
        return "LowerHalf(" + str(self.var0)+ ")"

class UpperHalf(Exp):
    def __init__(self, var0):
        super().__init__()
        self.var0 = var0

    def __repr__(self):
        return "UpperHalf(" + str(self.var0)+ ")"

class MergeHalf(Exp):
    def __init__(self, var0, var1):
        super().__init__()
        self.var0 = var0
        self.var1 = var1
    def __repr__(self):
        return "MergeHalf(" + str(self.var0)+ "," + str(self.var1) + ")"

class Identity(Exp):
    def __init__(self, var0):
        super().__init__()
        self.var0 = var0
    def __repr__(self):
        return "Identity(" + str(self.var0)+ ")"

class Butterfly(Exp):
    def __init__(self, var0, var1, var2):
        super().__init__()
        self.var0 = var0
        self.var1 = var1
        self.var2 = var2
    def __repr__(self):
        return "Butterfly(" + str(self.var0)+ "," + str(self.var1) + "," + str(self.var2) + ")"
    

def cse(in_length, out_length, linear_deps):
    ncol = in_length
    nrow = out_length
    while(True):
        order = list(range(ncol))
        #random.shuffle(order)
        hmax = 0
        for ii in range(ncol):
           for jj in range(ii+1,ncol):
                i = order[ii]
                j = order[jj]
                hw = 0
                for s in linear_deps:
                    if i in s and j in s:
                        hw += 1
                if hw > hmax: #  or (hw==hmax and random.randint(0,1)==0):
                    hmax = hw
                    maxi = i
                    maxj = j
        if hmax > 1:
            linear_deps.append({maxi,maxj})
            for i in range(nrow):
                s = linear_deps[i]
                if maxi in s and maxj in s:
                    s.remove(maxi)
                    s.remove(maxj)
                    s.add(ncol)
            ncol+=1
        else:
            break
    return (ncol-in_length,linear_deps)

def gen_assign(out_var_name,in_var_name,xor_list):
	assign_str = ""
	assign_str += "{} = ".format(out_var_name)
	equation = []
	for j in sorted(list(xor_list)):
		if j == -1:
			equation.append("1")
		else:
			equation.append( str(in_var_name[j]))
	assign_str += " ^ ".join(equation)
	if len(equation) == 0 :
		assign_str += "0"
	assign_str += ";\n"
	return assign_str

def gen_code(in_var_name ,in_length, out_var_name, out_length, linear_func, prefix=""):
    ret = ""
    linear_deps = get_linear_deps(in_length,out_length,linear_func)
    ncs,linear_deps = cse(in_length,out_length,linear_deps)
    cs_var = Var(ncs)
    out_var_bits = out_var_name.bits + cs_var.bits
    for i in list(range(out_length,out_length+ncs))+ list(range(out_length)):
        ret += gen_assign(prefix + str(out_var_bits[i]),in_var_name + cs_var.bits,
            linear_deps[i])
    return ret


def exp_graph_to_eqs(exp,eqs):
        if exp.__class__ == Var:
            return exp
        if exp.traverse is not None:
            return exp.traverse
        if exp.__class__ == Mul_Const:
            a = exp_graph_to_eqs(exp.var0,eqs)
            b = Var(len(a.bits))
            eqs.append((b,Mul_Const(a,exp.var1)))
            #TODO bits length maybe wrong because of mul constant
            exp.traverse = b
            return b
        elif exp.__class__ == Add:
            a = exp_graph_to_eqs(exp.var0,eqs)
            b = exp_graph_to_eqs(exp.var1,eqs)
            c = Var(len(a.bits))
            eqs.append((c,Add(a,b)))
            exp.traverse = c
            return c
        elif exp.__class__ == UpperHalf:
            a = exp_graph_to_eqs(exp.var0,eqs)
            b = Var(len(a.bits)//2)
            eqs.append((b,UpperHalf(a)))
            exp.traverse = b
            return b
        elif exp.__class__ == LowerHalf:
            a = exp_graph_to_eqs(exp.var0,eqs)
            b = Var(len(a.bits)//2)
            eqs.append((b,LowerHalf(a)))
            exp.traverse = b
            return b
        elif exp.__class__ == Mul:
            a = exp_graph_to_eqs(exp.var0,eqs)
            b = exp_graph_to_eqs(exp.var1,eqs)
            c = Var(len(a.bits))
            eqs.append((c,Mul(a,b)))
            exp.traverse = c
            return c
        elif exp.__class__ == MergeHalf:
            a = exp_graph_to_eqs(exp.var0,eqs)
            b = exp_graph_to_eqs(exp.var1,eqs)
            if len(a.bits)!=len(b.bits):
                raise "NOT MERGE HALF!"
            c = Var(len(a.bits)+len(b.bits))
            eqs.append((c,MergeHalf(a,b)))
            exp.traverse = c
            return c
        else:
            raise
    
    

def eq_to_str(eqs):
    code = ""
    for lhs,rhs in eqs:
        #print(lhs,rhs)
        if rhs.__class__ == Mul_Const:
            #code+="{} = mul({}, {})\n".format(lhs,rhs.var0,rhs.var1)
            tmp = [i for i in rhs.var0.bits] 
            tmp += [0]*(len(lhs.bits)-len(tmp))
            #print(tmp)
            if len(lhs.bits)==2:
                if rhs.var1 == 0:
                    code+="{} = 0;\n".format(lhs.bits[0])
                    code+="{} = 0;\n".format(lhs.bits[1])
                if rhs.var1 == 1:
                    code+="{} = {}\n".format(lhs.bits[0],tmp[0])
                    code+="{} = {}\n".format(lhs.bits[1],tmp[1])
                if rhs.var1 == 2:
                    code+="{} = {}\n".format(lhs.bits[0],tmp[1])
                    code+="{} = {} ^ {}\n".format(lhs.bits[1],tmp[0],tmp[1])
                if rhs.var1 == 3:
                    code+="{} = {} ^ {}\n".format(lhs.bits[0],tmp[0],tmp[1])
                    code+="{} = {}\n".format(lhs.bits[1],tmp[0])
            else:
                code+=gen_code(rhs.var0.bits,len(rhs.var0.bits),lhs,len(lhs.bits),
                    lambda x:gf832_mul(x,rhs.var1))
        elif rhs.__class__ == Butterfly:
            #code+="{} = mul({}, {})\n".format(lhs,rhs.var0,rhs.var1)
            #print(tmp)

            c0 = rhs.var0
            c1 = rhs.var1
            hh0 = Add(c0,Mul_Const(c1,rhs.var2))
            tmp = []
            r = exp_graph_to_eqs(hh0,tmp)

            hh1 = Var(len(c1.bits))
            tmp.append((hh1,Add(r,c1)))
            tmp.append((lhs,MergeHalf(r,hh1)))
            tmp_code1 = eq_to_str(tmp)


            c0 = rhs.var0
            c1 = rhs.var1
            hh1 = Add(c0,Mul_Const(c1,rhs.var2 ^ 1))
            tmp = []
            r = exp_graph_to_eqs(hh1,tmp)

            hh0 = Var(len(c0.bits))
            tmp.append((hh0,Add(r,c1)))
            tmp.append((lhs,MergeHalf(hh0,r)))
            tmp_code2 = eq_to_str(tmp)


            l = len(lhs.bits)//2
            def butterfly(x,w):
                c0 = x&(2**l - 1)
                c1 = x>>l
                hh0 = gf832_mul(c1,w) ^ c0
                hh1 = gf832_mul(c1,(w^1)) ^ c0
                return (hh1<<l)|hh0

            tmp_code3 = gen_code(rhs.var0.bits + rhs.var1.bits ,len(rhs.var0.bits)*2,lhs,l*2,
                lambda x:butterfly(x,rhs.var2))

            tmp_code1 = tmp_code1.replace(";","")
            tmp_code2 = tmp_code2.replace(";","")
            tmp_code3 = tmp_code3.replace(";","")

            if tmp_code2.count("^")<tmp_code1.count("^"):
                if tmp_code3.count("^")<tmp_code2.count("^"):
                    code+=tmp_code3
                    #code+=tmp_code2
                else:
                    code+=tmp_code2
            else:
                if tmp_code3.count("^")<tmp_code1.count("^"):
                    code+=tmp_code3
                    #print(">>>>>>")
                    #print(tmp_code3)
                    #print("=====")
                    #print(tmp_code1)
                    #print("<<<<<<")
                    #print(tmp_code3)
                    #print(tmp_code3.count("^"),tmp_code1.count("^"),tmp_code2.count("^"))
                    #code+=tmp_code1
                else:
                    code+=tmp_code1

            #[ print(i) for i in tmp]

        elif rhs.__class__ == Add:
            #code+="{} = {} ^ {}\n".format(lhs,rhs.var0,rhs.var1)
            for i,v in enumerate(lhs.bits):
                if i<len(rhs.var0.bits):
                    if i<len(rhs.var1.bits):
                        code+="{} = {} ^ {}\n".format(lhs.bits[i],rhs.var0.bits[i],rhs.var1.bits[i])
                    else:
                        code+="{} = {}\n".format(lhs.bits[i],rhs.var0.bits[i])
                else:
                    if i<len(rhs.var1.bits):
                        code+="{} = {}\n".format(lhs.bits[i],rhs.var1.bits[i])
                    else:
                        code+="{} = 0\n".format(lhs.bits[i])
            #if len(rhs.var0.bits) != len(rhs.var1.bits):
            #    print(len(rhs.var0.bits), len(rhs.var1.bits))
        elif rhs.__class__ == UpperHalf:
            for i,v in enumerate(lhs.bits):
                    code+="{} = {}\n".format(lhs.bits[i],rhs.var0.bits[i+len(lhs.bits)])
        elif rhs.__class__ == LowerHalf:
            for i,v in enumerate(lhs.bits):
                    code+="{} = {}\n".format(lhs.bits[i],rhs.var0.bits[i])
        elif rhs.__class__ == Mul:
            var0 = rhs.var0
            var1 = rhs.var1
            if len(lhs.bits) == 1:
                code+="{} = {} & {}\n".format(lhs.bits[0],var0.bits[0],var1.bits[0])
            if len(lhs.bits) == 2:
                a = rhs.var0
                a0 = LowerHalf(a)
                a1 = UpperHalf(a)
                b = rhs.var1
                b0 = LowerHalf(b)
                b1 = UpperHalf(b)
                ab0=Mul(a0,b0)
                ab1 = Add(Mul(a1,b0),Mul(a0,b1))
                ab2 = Mul(a1,b1)
                ab0 = Add(ab0,ab2)
                ab1 = Add(ab1,ab2)
                tmp = []
                r = exp_graph_to_eqs(MergeHalf(ab0,ab1),tmp)
                #[ print(i) for i in tmp]
                tmp.append((lhs,Identity(r)))
                code+=eq_to_str(tmp)
            if len(lhs.bits) == 4:
                a = rhs.var0
                a0 = LowerHalf(a)
                a1 = UpperHalf(a)
                b = rhs.var1
                b0 = LowerHalf(b)
                b1 = UpperHalf(b)


                a0b0=Mul(a0,b0)
                a1b1=Mul(a1,b1)
                a0a1xb0b1_a0b0 = Add(Mul(Add(a0,a1),Add(b0,b1)),a0b0)
                rd0=Mul_Const(a1b1,2)
                a0b0 = Add(a0b0,rd0)
                merge = MergeHalf(a0b0,a0a1xb0b1_a0b0)
                tmp = []
                r = exp_graph_to_eqs(merge,tmp)
                #[ print(i) for i in tmp]
                tmp.append((lhs,Identity(r)))
                code+=eq_to_str(tmp)
            if len(lhs.bits) == 8:
                a = rhs.var0
                a0 = LowerHalf(a)
                a1 = UpperHalf(a)
                b = rhs.var1
                b0 = LowerHalf(b)
                b1 = UpperHalf(b)
                a0b0 = Mul(a0,b0)
                ab = Mul(Add(a0,a1),Add(b0,b1))
                a1b1 = Mul(a1,b1)
                
                a0b1_a1b0_a1b1 = Add(ab,a0b0)
                a1b10x8 = Mul_Const(a1b1,8)
                merge = MergeHalf(Add(a1b10x8,a0b0),a0b1_a1b0_a1b1)
                tmp = []
                r = exp_graph_to_eqs(merge,tmp)
                tmp.append((lhs,Identity(r)))
                code+=eq_to_str(tmp)
            if len(lhs.bits) == 16:
                a = rhs.var0
                a0 = LowerHalf(a)
                a1 = UpperHalf(a)
                b = rhs.var1
                b0 = LowerHalf(b)
                b1 = UpperHalf(b)
                a0b0 = Mul(a0,b0)
                ab = Mul(Add(a0,a1),Add(b0,b1))
                a1b1 = Mul(a1,b1)
                
                a0b1_a1b0_a1b1 = Add(ab,a0b0)
                a1b10x8 = Mul_Const(a1b1,0x80)
                merge = MergeHalf(Add(a1b10x8,a0b0),a0b1_a1b0_a1b1)
                tmp = []
                r = exp_graph_to_eqs(merge,tmp)
                tmp.append((lhs,Identity(r)))
                code+=eq_to_str(tmp)
                

        elif rhs.__class__ == Identity:
            for i,v in enumerate(lhs.bits):
                    code+="{} = {}\n".format(lhs.bits[i],rhs.var0.bits[i])
        elif rhs.__class__ == MergeHalf:
            for i,v in enumerate(rhs.var0.bits):
                    code+="{} = {}\n".format(lhs.bits[i],rhs.var0.bits[i])
            for i,v in enumerate(rhs.var0.bits):
                    code+="{} = {}\n".format(lhs.bits[i+len(rhs.var0.bits)],rhs.var1.bits[i])
            
        else:
            print(rhs.__class__)
            raise
    return code


