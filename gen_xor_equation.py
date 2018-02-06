def get_linear_deps(in_length,out_length,linear_func):
	constant = linear_func(0)
	ret = [set() for i in range(out_length)]
	for i in range(out_length):
		if ((constant>>i)&1):
			ret[i].add(-1)
	for i in range(in_length):
		out = linear_func(1<<i) ^ constant
		for j in range(out_length):
			if ((out>>j)&1) != 0:
				ret[j].add(i)
	return ret

def gen_assign(out_var_name,in_var_name,xor_list):
	assign_str = ""
	assign_str += "{} = ".format(out_var_name)
	equation = []
	for j in sorted(list(xor_list)):
		if j == -1:
			equation.append("1")
		else:
			equation.append( "{}[{}]".format(in_var_name,j))
	assign_str += " ^ ".join(equation)
	if len(equation) == 0 :
		assign_str += "0"
	assign_str += ";\n"
	return assign_str

def gen_code(in_var_name ,in_length, out_var_name, out_length, linear_func, prefix=""):
	ret = ""
	for i in range(out_length):
		ret += gen_assign(prefix + "{}[{}]".format(out_var_name,i),in_var_name,
			get_linear_deps(in_length,out_length,linear_func)[i])
	return ret
			

if __name__ == "__main__":
	#print gen_code("x",5,"y",2,lambda x: (x^(x>>2)^2))
	from gf2_tower import *
	for c in range(0,16):
	    print(gen_code("x",4,"y",4,lambda x:gf16_mul(x,c)))
