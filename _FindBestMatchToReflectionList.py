tThets = [ 32.7, 42.1, 48.3, 64.77, 76.7 ] # in degrees
wavelength = 2.41 # in Angrstrom

#----- Instructions -----#
# Above is a list, tThets
# For your diffraction pattern, enter into tThets a list of 2-theta values, in degrees, where there are reflections
# Above is a variable, wavelength
# change wavelength to the wavelength your experiment was done at
# Put this file into a folder with all the .cif files that you want to search through
# Run the file from the anaconda prompt/command prompt/terminal (it is assumed the user knows how to do this)
# A file, 2ThetaHits.csv should appear in the folder with all the .cif files that you want to search through
# Have a look at 2ThetaHits.csv

import os
from os import listdir
import re
from numpy import array, sin, cos, pi, sqrt, dot, vstack, arcsin, savetxt, abs
from numpy import mod, array_equal, flip, zeros, ones
from numpy.linalg import inv, norm
from matplotlib import pyplot as plt
from cmath import exp as cexp

#thank you to John Machin, https://stackoverflow.com/questions/4703390/how-to-extract-a-floating-number-from-a-string
numeric_const_pattern = r'[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?'
rx = re.compile(numeric_const_pattern, re.VERBOSE)

cwd = os.getcwd()

files = [ f for f in os.listdir() if '.cif' in f and not 'Shortcut' in f ]
#print(files)
    
def getabc( file ):
    with open(file, encoding='utf8') as f:
        for line in f:
            if '_cell_length_a' in line:
                a = float( rx.findall( line )[0]  )
            elif '_cell_length_b' in line:
                b = float( rx.findall( line )[0]  )
            elif '_cell_length_c' in line:
                c = float( rx.findall( line )[0]  )
            elif '_cell_angle_alpha' in line:
                alpha = float( rx.findall( line )[0]  ) * pi / 180
            elif '_cell_angle_beta' in line:
                beta = float( rx.findall( line )[0]  ) * pi / 180
            elif '_cell_angle_gamma' in line:
                gamma = float( rx.findall( line )[0]  ) * pi / 180
                
        ahat = array( [ 1, 0, 0 ] )
        bhat = array( [ cos( gamma ), sin( gamma ), 0 ] )
        
        x = cos( beta )
        y = ( cos( alpha ) - cos( beta ) * cos(gamma) ) / sin( gamma )
        z = sqrt( 1 - x**2 - y**2 )
        
        chat = array( [x,y,z] )
                
        B = vstack( ( a*ahat, b*bhat, c*chat ) )#.transpose()
        
        return(B)
#B = getabc( 'Ce_C12_m1_CollCode2284.cif' )

def getAtoms( file ):
    with open(file, encoding='utf8') as f:
        look = False
        coords = []
        for line in f:
            if look:
                if 'loop_' in line or '#End' in line:
                    break
                else:
                    dat = line.split(' ')[4:6+1]
                    for i in range( len( dat ) ):
                        dat[i] = float( rx.findall( dat[i] )[0] )
                    coords.append(dat)
            if '_atom_site_occupancy' in line:
                look = True
    return array( coords )

def getHKLs( file, lambd, writeCSV = False ):
    knum = 2*pi / lambd
    B = getabc( file )
    R = 2 * pi * inv(B)
    h,k,l, tThet, d = ( [], [], [], [], [] )
    for i in range(-10,10):
        for j in range(-10,10):
            for kk in range(-10,10):
                q = norm( dot( R, array( [i,j,kk] ) ) )
                for n in range(1,3+1):
                    if( abs( n * q / 2 / knum ) <= 1 and q>0 ):
                        h.append( n*i )
                        k.append( n*j )
                        l.append( n*kk )
                        tThet.append( 2 * arcsin( n * q / 2 / knum ) * 180 / pi )
                        d.append( 2 * pi / n / q )
    sortx = array(tThet).argsort()
    hs, ks, ls, tTs, ds = ( array( h )[sortx], array( k )[sortx], array( l )[sortx], array( tThet )[sortx], array( d )[sortx] )
    
    keep = [True]
    for i in range( 1, len(hs) ):
        if abs( tTs[i] - tTs[i-1] ) < 0.001:
            keep.append(False)
        else:
            keep.append(True)
    data = vstack( ( hs[keep], ks[keep], ls[keep], tTs[keep], ds[keep] ) )#.transpose()
    if writeCSV:
        savetxt( file[:file.index('.')] + '.csv', data.transpose(), fmt='%i,%i,%i,%.4f,%.4f', header = 'h,k,l,2theta,d', delimiter = ',', comments ='' )
    return data

#getHKLs( 'Ce_C12_m1_CollCode2284.cif', 2.41 )

def getSymOps( file ):
    with open( file, encoding='utf8' ) as f:
        lines = f.readlines()
        ops = []
        count = 1
        lookEnd = False
        for line in lines:
            if lookEnd:
                if not 'loop_' in line:
                    id1 = line.index("'")
                    id2 = line.index("'", id1 + 1)
                    opsta = line[id1+1: id2].split(",")
                    op=[]
                    for opst in opsta:
                        if '-x' in opst:
                            sgn, i = ( -1, 0 )
                        elif 'x' in opst:
                            sgn, i = ( 1, 0 )
                        elif '-y' in opst:
                            sgn, i = ( -1, 1 )
                        elif 'y' in opst:
                            sgn, i = ( 1, 1 )
                        elif '-z' in opst:
                            sgn, i = ( -1, 2 )
                        elif 'z' in opst:
                            sgn, i = ( 1, 2 )
                        fracsta = re.findall( r"[-+]?\d/\d", opst )
                        if len(fracsta) > 0:
                            fracst = fracsta[0]
                            if '-' in fracst:
                                sgnfrac = -1
                            else:
                                sgnfrac = 1
                            frac = sgnfrac * float( re.findall(r"\d", fracst)[0] ) / float( re.findall(r"\d", fracst)[1] )
                        else:
                            frac = 0
                        op.append([sgn, i, frac])
                    ops.append(op)
                else:
                    break
            if '_space_group_symop_operation_xyz' in line:
                lookEnd = True
    return( ops )

def symOp( v, symOps ):
    vs = [v]
    for op in symOps:
        x=[]
        for oi in op:
            x.append( oi[0] * v[oi[1]] + oi[2] )
        vs.append( x )
    vs = array(vs)
    vs = mod(vs,1)
    uvs = array( [ vs[0] ] )
    for vec in vs:
        includeIt = True
        for uvec in uvs:
            if array_equal( vec, uvec ):
                includeIt = False
        if includeIt:
            uvs = vstack( ( uvs, vec ) )
    return(uvs)

def getHKLs2( file, lambd, writeCSV = True ):
    data = getHKLs( file, lambd, writeCSV = False )
    ops = getSymOps( file )
    atoms = getAtoms( file )
    allAtoms = []
    for atom in atoms:
        oppedAtoms = symOp( atom, ops )
        for oppedAtom in oppedAtoms:
            allAtoms.append( oppedAtom )
    keep = []
    Fc = []
    for i in range( len( data[0] ) ):
        S = 0.0 + 0.0J
        hkl = array( [ data[0][i], data[1][i], data[2][i] ] )
        for a in allAtoms:
            S += cexp( 2j*pi * dot( a, hkl ) )
        if abs(S) > 0.001:
            keep.append(True)
            Fc.append( abs(S) )
        else:
            keep.append(False)
    newData = []
    for i in range( len( data ) ):
        newData.append( data[i][keep] )
    newData.append(Fc)
    if writeCSV:
        savetxt( file[:file.index('.')] + '.csv', array(newData).transpose(), fmt='%i,%i,%i,%.4f,%.4f,%.4f', header = 'h,k,l,2theta,d,Fc', delimiter = ',', comments ='' )
    return newData
#getHKLs2( 'CeGaGe_CollCode621125.cif', 2.41 )

def compareHKLs(tThetas, files, lambd, writeCSV = False, name=''):

    hitMarkers = [0, 0.25, 0.5, 1, 2]
    allHits = []
    totalRefs=[]
    for marker in hitMarkers:
        if marker != 0:
            allHits.append( [] )
    for file in files:
        tTs = getHKLs2( file, lambd, writeCSV )[3]
        totalRefs.append( len( tTs ) )
        hits = zeros( len( allHits ) )
        for tThet in tThetas:
            looks = ones( len( allHits ) )
            for tT in tTs:
                for i in range( len(allHits) ):
                    if  hitMarkers[i] <= abs( tT - tThet ) < hitMarkers[i+1] and looks[i]:
                        hits[i] += 1
                        looks[i] = 0 # or False
        for i in range( len(allHits) ):
            allHits[i].append( hits[i] )
    allHits = array( allHits )
    totalRefs = array( totalRefs )
    sortx = flip( (allHits[0]**2 / totalRefs ).argsort() )# The "goodness" of the match goes like the square of the number of matches, divided by the total number of matches since there are always massive heterostructures with reflections at nearly every half degree in two theta
    files = array(files)[sortx]
    totalRefs = totalRefs[sortx]
    for i in range( len( allHits ) ):
        allHits[i] = allHits[i][sortx]
    
    with open( '2ThetaHits'+name+'.csv', 'w+' ) as f:
        print( 'lambda, %.3f' % lambd, file=f )
        f.write('2Thetas, ')
        for tThet in tThetas:
            f.write( '%.3f,' % tThet)
        f.write('\ncif,')
        for i in range( len(allHits) ):
            f.write('hits within %.2f degree,' % ( hitMarkers[i+1] ))
        f.write('reflections total,\n')
        for j in range( len( allHits[0] ) ):
            f.write( '%s,' % ( files[j] ) )
            for i in range( len( allHits ) ):
                f.write('%i,' % ( allHits[i][j] ) )
            f.write('%i,\n' %( totalRefs[j] ))
        
tThets_flux10K = [ 32.75, 47.9, 53, 56, 58.5, 63.8, 65, 75.9, 88.5, 90.6, 92.1, 110.8, 127 ]
tThets_zr10K = [ 32.7, 42.1, 48.3, 64.77, 76.7 ]

compareHKLs( tThets, files, wavelength )
#for val in tThets_flux10K:
#    compareHKLs( [val], files, 2.41, False, name = 'flux' + str(val) )

'''
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.scatter( vs[:,0], vs[:,1], vs[:,2], c='blue' )
ax.scatter( [0,0,0,0,1,1,1,1], [0,1,0,1,0,1,0,1], [0,0,1,1,0,0,1,1], c='black' )
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()
'''
'''
count = 0
for file in files:
    with open( file ) as f:
        lines = f.readlines()
        for line in lines:
            if '_space_group_IT_number' in line:
                count += 1
                print(file, line)
                break

print( "files: %i, sg#s: %i" % (len(files), count) )
'''