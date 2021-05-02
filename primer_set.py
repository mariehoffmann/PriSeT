import sys

f = "./primer_set.txt"

KAPPA_MIN = 10
KAPPA_MAX = 11

PRIMER_MIN_TM = 50.0
PRIMER_MAX_TM = 64.0

CG_MIN_CONTENT = .4
CG_MAX_CONTENT = .6

SUFFIX_10 = (1 << 10) - 1
SUFFIX_10 = (1 << 20) - 1

# return 0 code and elongated versions are not suitable
# return 1 if code not suitable, but elongated one might be
# return 2 if code is suitable, and maybe elongations
def filter(code):
    CG = 0
    AT = 0
    while (code != 1):
        if (code & 3) | ((code & 3) + 1) == 3:
            CG += 1
        else:
            AT += 1
        # check consecutive runs
        if code > (1 << 10):
            code_suffix = code & SUFFIX_10
            if (code_suffix == 0 or code_suffix == 0b0101010101 or code_suffix == 0b1111111111 or code_suffix == 0b1010101010):
                return 0
        if code > (1 << 20):
            code_suffix = code & SUFFIX_20
            if (code_suffix == 0b00010001000100010001 or code_suffix == 0b00100010001000100010 or code_suffix == 0b00110011001100110011 or \
                code_suffix == 0b01000100010001000100 or code_suffix == 0b01100110011001100110 or code_suffix == 0b01110111011101110111 or \
                code_suffix == 0b10001000100010001000 or code_suffix == 0b10011001100110011001 or code_suffix == 0b10111011101110111011 or \
                code_suffix == 0b11001100110011001100 or code_suffix == 0b11011101110111011101 or code_suffix == 0b11101110111011101110 \
                ):
                return 0

    if CG/(CG + AT) < CG_MIN_CONTENT or CG/(CG + AT) > CG_MAX_CONTENT:
        return 1
    if ((AT << 1) + (CG << 2) < PRIMER_MIN_TM) or ((AT << 1) + (CG << 2) > PRIMER_MAX_TM):
        return 1
    return 2

# check for CG clamps on both sides (one side is ok for forward on 5' and reverse on 3')
def filter_CG_clamp(code):
    CG_suffix = 0
    CG_prefix
    suffix = 0
    while code != 1:
        if suffix < 5 and ((code & 3) | (code & 3 + 1)) == 3:
            CG_suffix += 1
            suffix += 1
        if code < (1 << 10) and ((code & 3) | (code & 3 + 1)) == 3:
            CG_prefix += 1
        code >>= 2
    return False if (CG_suffix > 3 or CG_prefix > 3) else True


def compute_primers(offset = 0):
    ctr = 0
    heap = [1]  # valid prefixes
    with open(f, 'a') as f_hdlr:
        while len(heap) > 0:
            code = heap.pop()
            if code > (1 << (2 * KAPPA_MIN)): # add code if large enough
                if offset < ctr and filter_CG_clamps(code):
                    f.write(dna_decoder(code) + "\n")
                ctr += 1
                if ctr > 10000:
                    print(ctr)
            code <<= 2
            if code < (1 << (2 * KAPPA_MAX + 1)):
                if filter(code):
                    heap.append(code)
                if filter(code + 1):
                    heap.append(code + 1)
                if filter(code + 2):
                    heap.append(code + 2)
                if filter(code + 3):
                    heap.append(code + 3)
    print(ctr, " primers found")


if __name__ == "__main__":
    #sys.argc == 2
    compute_primers() #int(sys.argv[1]))
