#!/usr/bin/env pypy3
import sys,os
import gzip
import struct
from format import *

infile = sys.argv[1]
outfile = sys.argv[2]

ifp = gzip.open(infile, mode='rb')
ofp = gzip.open(outfile,mode='wb',compresslevel=1)

def ol(line):
    # print(line,end='')
    ofp.write(line.encode('ascii'))

all_binary = ifp.read()
ifp.close()
binary_index = 0

def rd(fmt):
    global all_binary
    global binary_index
    size = struct.calcsize(fmt)
    tmp = struct.unpack_from(fmt,all_binary,binary_index)
    binary_index += size 
    return tmp

while True:
    if True: # M lines
        fmt = fshort 
        try:
            unpackdata = rd(fmt)
        except struct.error:
            break
        smilelen, = unpackdata
        
        fmt = fM
        unpackdata = rd(fmt)
        zinc_index, codex, atoms,bonds,nxyz,nconf,nset,nrigid,nmline,nclu, charge, polar, apolar, tsol, surf,m4,m5 = unpackdata

        fmt = f"{smilelen}" + fstring
        smile, = rd(fmt)
        smile = smile.decode('ascii')

        if not codex:
            codex = 'none'
        if not nmline:
            nmline = 5

        zinc = "ZINC%012d"%zinc_index
        fline = "M %16s %9s %3d %3d %6d %6d %6d %6d %6d %6d\n"%(zinc,codex,atoms,bonds,nxyz,nconf,nset,nrigid,nmline,nclu)
        ol(fline)
        sline = "M %+9.4f %+10.3f %+10.3f %+10.3f %9.3f\n"%(charge,polar,apolar,tsol,surf)
        ol(sline)
        tline = "M %-76s\n"%smile
        ol(tline)
        nline = 'M %-76s\n'%"NO_LONG_NAME"
        ol(nline)
        nline = 'M  +999.9990\n'
        ol(nline)

    if True: # A lines 
        index = 0
        for _ in range(atoms):
            index += 1
            fmt = fA
            unpackdata = rd(fmt)
            atomname,atomtype,dt,co,charge,polar,apolar,total,surf = unpackdata
            charge = charge * 0.01
            polar = polar * 0.01
            apolar = apolar * 0.01
            total = total * 0.01
            surf = surf * 0.01
            atomname = atomname.decode('ascii')
            atomtype = atomtype.decode('ascii')
            
            aline = "A %3d %-4s %-5s %2d %2d %+9.4f %+10.3f %+10.3f %+10.3f %9.3f\n"%(index,atomname,atomtype,dt,co,charge,polar,apolar,total,surf)
            ol(aline)

    if True: # B lines
        index = 0
        for _ in range(bonds):
            index += 1
            fmt = fB
            unpackdata = rd(fmt)
            atom1,atom2,btype = unpackdata
            btype = fbtype[btype]
            bline = "B %3d %3d %3d %-2s\n"%(index,atom1,atom2,btype)
            ol(bline)

    if True: # X 
        index = 0
        for _ in range(nxyz):
            index += 1
            fmt = fX
            atom,confnum,x,y,z = rd(fmt)
            xline = "X %9d %3d %6d %+9.4f %+9.4f %+9.4f\n"%(index,atom,confnum,x,y,z)
            ol(xline)

    if True: # R
        index = 0
        for _ in range(nrigid):
            index += 1
            fmt = fR
            co,x,y,z = rd(fmt)
            rline = "R %6d %2d %+9.4f %+9.4f %+9.4f\n"%(index,co,x,y,z)
            ol(rline)

    if True: # C
        index = 0
        for _ in range(nconf):
            index += 1
            fmt = fC
            start,stop = rd(fmt)
            cline = "C %6d %9d %9d\n"%(index,start,stop)
            ol(cline)

    if True: # S
        index = 0 
        for _ in range(nset):
            index += 1
            fmt = fS1
            nconfs,broken,hydrogens,omega_energy = rd(fmt)
            nlines = (nconfs+7)//8
            shline = "S %6d %6d %3d %1d %1d %+11.3f\n"%(index,nlines,nconfs,broken,hydrogens,omega_energy)
            ol(shline)

            fmt = fshort*nconfs
            all_confs = rd(fmt)
            slineindex = 0
            confindex = 0
            while True:
                slineindex += 1
                if confindex == nconfs:
                    break
                if nconfs - confindex > 8:
                    c_confs = all_confs[confindex:confindex+8]
                    confindex += 8
                else:
                    c_confs = all_confs[confindex:]
                    confindex = nconfs

                c_len = len(c_confs)
                csline = ("S %6d %6d %1d"+" %6d"*c_len +"\n")%( (index, slineindex, c_len) + tuple(c_confs) )
                ol(csline)
            
    if True: # D
        index = 0
        lineindex = 0
        for _ in range(nclu):
            index += 1 
            fmt = fD1
            nlines,sstart,send,matchstart,matchend = rd(fmt)

            cline = "D %6d %6d %6d %3d %3d %3d\n"%(index,sstart,send,matchstart,matchend,nlines)
            ol(cline)

            for tmp in range(matchstart):
                lineindex += 1
                fmt = fD2
                color,x,y,z = rd(fmt)
                cline = "D %6d %2d %+9.4f %+9.4f %+9.4f\n"%(lineindex,color,x,y,z)
                ol(cline)

    if True: # E
        eline = "E\n"
        ol(eline)

ofp.close()

