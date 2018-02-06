# Copyright (C) 2017 Ming-Shing Chen


def gf2_mul( a , b ):
  return a&b

# gf4 := gf2[x]/x^2+x+1
# 4 and , 3 xor
def gf4_mul( a , b ):
  a0 = a&1
  a1 = (a>>1)&1
  b0 = b&1
  b1 = (b>>1)&1
  ab0 = a0&b0
  ab1 = (a1&b0)^(a0&b1)
  ab2 = a1&b1
  ab0 ^= ab2
  ab1 ^= ab2
  ab0 ^= (ab1<<1)
  return ab0

# gf16 := gf4[y]/y^2+y+x
# gf16 mul: xor: 18  ,and: 12
def gf16_mul( a , b ):
  a0 = a&3
  a1 = (a>>2)&3
  b0 = b&3
  b1 = (b>>2)&3
  a0b0 = gf4_mul( a0 , b0 )
  a1b1 = gf4_mul( a1 , b1 )
  a0a1xb0b1_a0b0 = gf4_mul( a0^a1 , b0^b1 ) ^ a0b0
  rd0 = gf4_mul( 2 , a1b1 )
  a0b0 ^= rd0
  return a0b0|(a0a1xb0b1_a0b0<<2)


# gf256 := gf16[x]/x^2 + x + 0x8
def gf256_mul( a , b ):
  a0 = a&15
  a1 = (a>>4)&15
  b0 = b&15
  b1 = (b>>4)&15
  ab0 = gf16_mul( a0 , b0 )
  ab1 = gf16_mul( a1 , b0 ) ^ gf16_mul( a0 , b1 )
  ab2 = gf16_mul( a1 , b1 )
  ab0 ^= gf16_mul( ab2 , 8 )
  ab1 ^= ab2
  ab0 ^= (ab1<<4)
  return ab0

#382 bit operations
def gf216_mul( a , b ):
  a0 = a&0xff
  a1 = (a>>8)&0xff
  b0 = b&0xff
  b1 = (b>>8)&0xff
  a0b0 = gf256_mul( a0 , b0 )
  a1b1 = gf256_mul( a1 , b1 )
  #a0b1_a1b0 = gf16_mul( a0^a1 , b0^b1 ) ^ a0b0 ^ a1b1     ^a1b1
  a0b1_a1b0 = gf256_mul( a0^a1 , b0^b1 ) ^ a0b0
  rd0 = gf256_mul( a1b1 , 0x80 )
  return (a0b1_a1b0<<8)|(rd0^a0b0)

#gf65536 := gf256[x]/x^2 + x + 0x80
def gf65536_mul( a , b ):
  a0 = a&0xff;
  a1 = (a>>8)&0xff;
  b0 = b&0xff;
  b1 = (b>>8)&0xff;
  ab0 = gf256_mul( a0 , b0 );
  ab2 = gf256_mul( a1 , b1 );
  ab1 = gf256_mul( a0^a1 , b0^b1 )^ab0;
  return (ab1<<8)^(ab0^gf256_mul(ab2,0x80));


#gf832 := gf65536[x]/x^2 + x + 0x8000
def gf832_mul( a , b ):
  a0 = a&0xffff;
  a1 = (a>>16)&0xffff;
  b0 = b&0xffff;
  b1 = (b>>16)&0xffff;
  ab0 = gf65536_mul( a0 , b0 );
  ab2 = gf65536_mul( a1 , b1 );
  ab1 = gf65536_mul( a0^a1 , b0^b1 )^ab0;
  return (ab1<<16)^(ab0^gf65536_mul(ab2,0x8000));

def gf832_inv(a) :
  r = a
  for i in range(2,32):
    r = gf832_mul(r,r)
    r = gf832_mul(r,a)
  return gf832_mul(r,r)
