#!/usr/bin/env python3
import sys,os
import gzip
import struct

infile = sys.argv[1]
outfile = sys.argv[2]

ofp = open(outfile,'wb') 

# for line in gzip.open(infile):
#     if line[0] == 69: 
#         print( line )

def next_fami( infile ):
    lines = list()
    for line in gzip.open(infile):
        if b"ZINC" in line and line[0] == 77 :
            lines = list()
            lines.append(line)
        elif line[0] == 69:
            lines.append(line)
            yield lines
        else:
            lines.append(line)
        
if True:
    for lines in next_fami( infile ):
        if True: # M lines
            # decode first line
            line = lines[0].decode('ascii')
            zinc = line[2:21]
            zinc_index = int(zinc[4:])
            codex = line[21:28].strip()
            assert codex == 'none'
            codex = None
            atoms = int( line[29:32] )
            bonds = int( line[33:36] )
            nxyz  = int( line[37:43] )
            nconf = int( line[44:50] )
            nset  = int( line[51:57] )
            nrigid= int( line[58:64] )
            nmline= int( line[65:72] )
            if nmline == 5:
                nmline = None
            nclu  = int( line[73:]   )

            #012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
            #M ZINC000008615381      none  40_ 42_    40_     1_     1_    21_     5_     1
            #M %16s %9s %3d %3d %6d %6d %6d %6d &6d %6d\n

            # decode second line
            line = lines[1].decode('ascii')
            charge = float( line[2:11]  )
            polar  = float( line[12:22] )
            apolar = float( line[23:33] )
            tsol   = float( line[34:44] )
            surf   = float( line[45:]   )

            #012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
            #M_  -1.0000_   -45.900_    +6.920_   -38.980   300.990
            #M %+9.4f %+10.3f %+10.3f %+10.3f %9.3f\n

            # decode third line (smiles) 
            line = lines[2].decode('ascii')
            smile = line[2:].strip().encode('ascii')
            smilelen = len(smile)

            #M CC1=CC2=C(C=C1C)C(C)(C)C1=C(N=C(N)[N-]C1=O)N2C
            #M %77s\n

            # decode fourth line ( NO_LONG_NAME )
            assert b'NO_LONG_NAME' in lines[3]
            m4 = None

            #M %77s\n

            # decode fifth line ( +999.9990 )
            assert b'M  +999.9990\n' in lines[4]
            m5 = None
            
            #M %77s\n

            
            buffer = struct.pack(f"I?BBHHHH?Bfffff{smilelen}s??", zinc_index, codex, atoms,bonds,nxyz,nconf,nset,nrigid,nmline,nclu, charge, polar, apolar, tsol, surf,smile,m4,m5 )

            ofp.write(buffer)

        if True: # A
            readindex = 4
            for _ in range(atoms):
                readindex += 1
                line = lines[readindex] 

                #012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
                #A   1_C1  _C.3  _ 5_ 7_  -0.1100_    -0.550_    +0.190_    -0.360_    9.600
                #A %3d %-4s %-5s %2d %2d %+9.4f %+10.3f %+10.3f %+10.3f %9.3f\n
                index = None
                atomname = line[6:10].decode("ascii").encode("ascii")
                atomtype = line[11:16].decode("ascii").encode("ascii")
                dt     = int(   line[17:19] )
                co     = int(   line[20:22] )
                charge = float( line[23:32] )
                polar  = float( line[33:43] )
                apolar = float( line[44:54] )
                total  = float( line[55:65] )
                surf   = float( line[66:] )
                buffer = struct.pack(f"4s5sHHfffff",atomname,atomtype,dt,co,charge,polar,apolar,total,surf)
                ofp.write(buffer)

        if True: # B

            for _ in range(bonds):
                readindex += 1
                line = lines[readindex] 

                #012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
                #B   5_  2_  7_ar
                #B %3d %3d %3d %-2s\n
                index = None
                atom1 = int(line[6:9])
                atom2 = int(line[10:13])
                btype = line[14:16].decode("ascii").encode("ascii")
                buffer = struct.pack(f"HH2s",atom1,atom2,btype)
                ofp.write(buffer)

            pass
        if True: # X
            for _ in range(nxyz):
                readindex += 1
                line = lines[readindex] 


                #012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
                #X %9d %3d %6d %+9.4f %+9.4f %+9.4f\n
                #X         1   1_     1_  -0.2954   +1.4879   +0.2588
                index = None
                atom = int( line[12:15] )
                confnum =  int( line[16:22] )
                x  = float( line[23:32] )
                y = float( line[33:42] )
                z = float( line[43:] )
                buffer = struct.pack(f"HHfff",atom,confnum,x,y,z)
                ofp.write(buffer)

            pass

        if True: # R
            for _ in range(nrigid):
                readindex += 1
                line = lines[readindex] 

                #012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
                #R %3d %2d %+9.4f %+9.4f %+9.4f\n
                #R      1  7   -0.2954   +1.4879   +0.2588
                index = None
                co = int( line[9:11] )
                x = float( line[12:21] )
                y = float( line[22:31] )
                z = float( line[32:] )

                buffer = struct.pack(f"Hfff",co,x,y,z)
                ofp.write(buffer)

            pass

        if True: # C
            for _ in range(nconf):
                readindex += 1
                line = lines[readindex] 


                #012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
                #C %6d %9d %9d\n
                #C      1         1        31

                index = None
                start = int(line[9:18])
                stop  = int(line[19:])
                buffer = struct.pack(f"HH",start,stop)
                ofp.write(buffer)

        if True: # def S function
            def next_set_lines( setlines ):
                flag_header = True
                lines = list()
                tmp_count = 0
                for line in setlines:
                    if flag_header:
                        line_num = float(line[9:15])
                        flag_header = False
                        lines = list()
                        lines.append(line)
                    else:
                        tmp_count += 1
                        if tmp_count == line_num:
                            lines.append(line)
                            yield lines
                            flag_header = True
                            tmp_count = 0
                            continue
                        else:
                            lines.append(line)
        if True: # test
            tmp = readindex
            while True:
                tmp += 1
                if lines[tmp][0] == 83:
                    pass
                else:
                    tmp -= 1
                    break
            setlines = lines[ readindex +1:tmp+1 ]
            readindex = tmp

            for this in next_set_lines( setlines ):
                fline = this[0]
                slines = this[1:]
                setnum = None 
                nlines = int( fline[9:15] )
                nconfs = int( fline[16:19] )
                broken = int( fline[20:21] )
                # assert broken == 0
                # broken = False
                hydrogens = int( fline[22:23] )
                # assert hydrogens == 0
                # hydrogens = False
                omega_energy = float( fline[24:35] )

                buffer = struct.pack(f"H??f",nconfs,broken,hydrogens,omega_energy)
                ofp.write(buffer)

                for line in slines:
                    nconf = int(line[15:17])
                    # print("nconf",nconf)
                    for _ in range(nconf):
                        conf_number = int(line[7*_+18:7*_+24])
                        buffer = struct.pack(f"H",conf_number)
                        ofp.write(buffer)

                #012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
                #S      1      1   7 0 0      +0.000
                #S      1      1 7      1     34     48     49    202    203    204
                #S %6d %6d %3d %1d %1d %+11.3f\n
                #S %6d %6d %1d %6d %6d %6d %6d %6d %6d %6d %6d\n 
            pass

        if True: # D
            readindex += 1
            dheader_line = lines[readindex]
            

            #012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
            #D %6d %6d %6d %3d %3d %3d\n
            #D %3d %2d %+9.4f %+9.4f %+9.4f\n
            #D      1      1      3   5   1   5
            #D      1  7   -2.3383   -0.3063   +0.8216
            #D      2  7   -3.2398   -1.0436   +0.6238
            #D      3  7   -1.4261   -0.4810   -0.1858
            #D      4  7   -2.9312   -1.7105   -0.5142
            #D      5  7   -1.8143   -1.3433   -0.9957 
            pass
        if True: # E
            pass

ofp.close()


