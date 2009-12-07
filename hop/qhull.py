# $Id$
# Hop --- a framework to analyze solvation dynamics from MD simulations
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
qhull
=====

Interface to some functions of the `qhull`_ and `qconvex`_ programs. These must
be installed separately (see links).

The main functionality is to define a region in space within the convex hull of
a protein. The hull is typically defined by a selection of atoms and written as
a "density" file for use in :mod:`hop`.

.. _qhull: http://www.qhull.org/html/index.htm
.. _qconvex: http://www.qhull.org/html/qconvex.htm

"""
from __future__ import with_statement

from subprocess import Popen
import numpy

#: Comparisons of distances less than EPSILON yield equal.
EPSILON = 1e-6

def make_ca_points(psf="1ifc_xtal.psf", pdb="1ifc_xtal.pdb", coordinates="ca.dat"):
    from MDAnalysis import Universe

    u = Universe(psf, pdbfilename=pdb)
    CA = u.selectAtoms('name CA').coordinates()

    with open(coordinates, 'w') as data:
        data.write('%d\n' % CA.shape[1])  # dimension
        data.write('%d\n' % CA.shape[0])  # number of points
        for point in CA:
            data.write("%f %f %f\n" % tuple(point))

    print "Wrote %(coordinates)r." % vars()


class ConvexHull(object):
    """The convex hull of a set of points.

    The convex hull is calculated with the `qhull`_ program.

    .. _qhull: http://www.qhull.org/
    """

    def __init__(self, coordinates, planes="planes.dat", vertices="vertices.dat"):
        """Compute convex hull and populate data structures.

        :Arguments:
        - coordinates: input suitable for qconvex
        - planes: filename for normal form plane file (``qconvex n``)
        - vertices: filename for vertices coordinates file (``qconvex p``)
        """

        self.files = {'coordinates': coordinates,
                      'planes': planes,
                      'vertices': vertices,}

        args = ['n', 'TO', "'"+self.files['planes']+"'"]
        rc = self.qconvex(args)
        print "Wrote planes %(planes)r [%(rc)r]" % vars()

        args = ['p', 'TO', "'"+self.files['vertices']+"'"]
        rc = self.qconvex(args)    
        print "Wrote vertices %(vertices)r [%(rc)r]" % vars()

        self.planes = self.read_planes()
        self.vertices = self.read_vertices()


    def qconvex(self, args):
        with open(self.coordinates) as coord: 
            Q = Popen(['qconvex']+args, stdin=coord)
            rc = Q.wait()
        return rc

    def read_vertices(self):
        """Read vertices from qconvex p file.

        Numpy array of points [[x,y,z], ...]
        """
        with open(self.files['vertices']) as vertices:
            vertices.readline()
            npoints = vertices.readline()
            a = []
            for line in vertices:
                a.append(map(float, line.strip().split()))
        if len(a) != npoints:
            raise IOError("Wrong number of datapoints %d, should be %d" % (len(a), npoints))
        return numpy.array(a)

    def read_planes(self):
        """Read planes from qconvex n file.

        Numpy array [[n1,n2,n3,-p], ...] for planes n*x = -p.

        Planes are oriented and point outwards.
        """
        with open(self.files['planes']) as planes:
            planes.readline()  # dimension (3+1)
            npoints = planes.readline()
            a = []
            for line in planes:
                a.append(map(float, line.strip().split()))
        if len(a) != npoints:
            raise IOError("Wrong number of datapoints %d, should be %d" % (len(a), npoints))
        return numpy.array(a)

    def _data_reader(self, name):
        """Read simple data structures from qhull files.

        :Arguments:
        - name: keyword in self.files

        File format::
           dimension
           npoints
           x1 x2 x3 x4 ...
        """
        with open(self.files[name]) as data:
            data.readline()  # dimension (3+1)
            npoints = data.readline()
            a = []
            for line in data:
                a.append(map(float, line.strip().split()))
        if len(a) != npoints:
            raise IOError("Wrong number of datapoints %d, should be %d" % (len(a), npoints))
        return numpy.array(a)
        

    def point_inside(self, point, planes):
        """Check if point [x,y,z] is inside the polyhedron defined by plains.

        Iff for all i: plain[i](x,y,z) = ax + by +cz + d > 0 <==> (x,y,z) outside
        """
        # crappy implementation, I am sure one can do this better with broadcasts
        # or a better algorithm, eg 
        # http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
        # 1. shoot semi-inifinite ray. 2. count how many faces F it crosses
        # crosses iff (1) P under plane of F (2) projection of P on plane is inside F (2D problem)
        # (Problems if ray hits a vertex.)

        # check ALL faces.... :-p
        return numpy.all(numpy.dot(planes[:,:3], point) + planes[:,3] <= EPSILON)

    def points_inside(self, points, planes):
        """Return bool array for all points:

        True: inside
        False: outside

        :Arguments:
        - points = [[x1,y1,z1], ...]
        - planes: normal forms of planes

        :Returns:
        [True, False, True, ...]
        """
        # crappy code, should optimize for numpy
        return numpy.array([point_inside(point, planes) for point in points])

    def write_vertices_pdb(filename, pdb="vertices.pdb"):
        points = read_vertices(filename)
        ppw = _PrimitivePDBWriter(pdb)
        ppw.write(points)
        print "Wrote points to pdb file %(pdb)r." % vars()


class _PrimitivePDBWriter(object):
    """PDB writer that implements a subset of the PDB 3.2 standard.
    http://www.wwpdb.org/documentation/format32/v3.2.html
    """
    #          1         2         3         4         5         6         7         8
    # 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.
    # ATOM__seria nameAres CressI   xxxxxxxxyyyyyyyyzzzzzzzzOCCUPAtempft          elCH
    # ATOM  %5d   %-4s %-3s %4d %1s %8.3f   %8.3f   %8.3f   %6.2f %6.2f           %2s
    #                 %1s  %1s                                                      %2d
    #            =        =      ===                                    ==========
    # ATOM  %5d %-4s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2d
    # ATOM  %(serial)5d %(name)-4s%(altLoc)1s%(resName)-3s %(chainID)1s%(resSeq)4d%(iCode)1s   %(x)8.3f%(y)8.3f%(z)8.3f%(occupancy)6.2f%(tempFactor)6.2f          %(element)2s%(charge)2d

    fmt = {'ATOM':   "ATOM  %(serial)5d %(name)-4s%(altLoc)1s%(resName)-3s %(chainID)1s%(resSeq)4d%(iCode)1s   %(x)8.3f%(y)8.3f%(z)8.3f%(occupancy)6.2f%(tempFactor)6.2f          %(element)2s%(charge)2d\n",
           'REMARK': "REMARK     %s\n",
           'TITLE':  "TITLE    %s\n",
           }

    def __init__(self,filename):
        self.filename = filename
        self.pdb = open(self.filename,'w')

    def close(self):
        self.pdb.close()

    def write(self,coordinates,name="CA",resname="VRT",resid=1):
        """Write coordinates as CA.

        write(points)
        """
        
        self.TITLE("points as CA")
        for i, atom in enumerate(coordinates):
            self.ATOM(serial=i+1, name=name.strip(), resName=resname.strip(), resSeq=resid,
                      x=coordinates[i,0], y=coordinates[i,1], z=coordinates[i,2])
        # get bfactor, too?
        self.close()

    def TITLE(self,*title):
        """Write TITLE record.
        http://www.wwpdb.org/documentation/format32/sect2.html
        """        
        line = " ".join(title)    # should do continuation automatically
        self.pdb.write(self.fmt['TITLE'] % line)

    def REMARK(self,*remark):
        """Write generic REMARK record (without number).
        http://www.wwpdb.org/documentation/format32/remarks1.html
        http://www.wwpdb.org/documentation/format32/remarks2.html
        """
        line = " ".join(remark)
        self.pdb.write(self.fmt['REMARK'] % line)
        
    def ATOM(self,serial=None,name=None,altLoc=None,resName=None,chainID=None,
             resSeq=None,iCode=None,x=None,y=None,z=None,occupancy=1.0,tempFactor=0.0,
             element=None,charge=0):
        """Write ATOM record. 
        http://www.wwpdb.org/documentation/format32/sect9.html
        Only some keword args are optional (altLoc, iCode, chainID), for some defaults are set.

        All inputs are cut to the maximum allowed length. For integer
        numbers the highest-value digits are chopped (so that the
        serial and reSeq wrap); for strings the trailing characters
        are chopped.

        Note: Floats are not checked and can potentially screw up the format.
        """
        for arg in ('serial','name','resName','resSeq','x','y','z',
                    'occupancy','tempFactor','charge'):
            if locals()[arg] is None:
                raise ValueError('parameter '+arg+' must be defined.')
        serial = int(str(serial)[-5:])  # check for overflow here?
        name = name[:4]
        if len(name) < 4:
            name = " "+name   # customary to start in column 14
        altLoc = altLoc or " "
        altLoc= altLoc[:1]
        resName = resName[:3]
        chainID = chainID or ""   # or should we provide a chainID such as 'A'?
        chainId = chainID[:1]
        resSeq = int(str(resSeq)[-4:]) # check for overflow here?
        iCode = iCode or ""
        iCode = iCode[:1]
        element = element or name.strip()[0]  # could have a proper dict here...
        element = element[:2]
        self.pdb.write(self.fmt['ATOM'] % vars())        
        
    def __del__(self):
        self.close()