"""Put an OpenSWMM 6 .out header into EPA SWMM 5.2.4's spelling, so the two files can be
compared byte for byte on results alone.

OpenSWMM 6 differs from EPA 5.2.4 in exactly seven header bytes on every file:
  * the version word: 60000 (OpenSWMM 6.0) against 52004 (EPA 5.2.4);
  * the input-property codes. OpenSWMM's enum (legacy output.c) adds OUTPUT_INVERT,
    OUTPUT_MAX_DEPTH and OUTPUT_OFFSET, so node codes read 0,2,4 (EPA 0,2,3) and link
    codes 0,6,7,4,8 (EPA 0,4,4,3,5, where EPA writes INPUT_OFFSET for both ends).
The values that follow the codes are the same quantities in the same order.
"""
import struct, sys, hashlib

def normalize(b: bytes) -> bytes:
    b = bytearray(b); o = 0
    def i():
        nonlocal o; v = struct.unpack_from('<i', b, o)[0]; o += 4; return v
    magic, ver, fu, ns, nn, nl, npol = [i() for _ in range(7)]
    assert magic == 516114522, 'not a SWMM binary output file'
    struct.pack_into('<i', b, 4, 52004)
    for n in (ns, nn, nl, npol):
        for _ in range(n): o += 4 + i()
    o += 4 * npol                              # pollutant units
    k = i(); o += 4 * k + 4 * k * ns           # subcatchment: count, codes, values
    k = i(); assert k == 3
    struct.pack_into('<3i', b, o, 0, 2, 3); o += 12 + 12 * nn
    k = i(); assert k == 5
    struct.pack_into('<5i', b, o, 0, 4, 4, 3, 5)
    return bytes(b)

if __name__ == '__main__':
    for p in sys.argv[1:]:
        print(hashlib.sha256(normalize(open(p, 'rb').read())).hexdigest(), p)
