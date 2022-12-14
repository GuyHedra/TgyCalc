import math

""" simpleTools.py is a simplified version of polyTools intended for use within TgyCalc"""


class Point(tuple):
    """ Point is used to store and operate on points and free vectors in 3-space.
    Free vectors are vectors whose initial point is the origin """
    length = 3

    def __init__(self, coords):
        tuple.__init__(self)
        if (isinstance(coords, tuple) and len(coords) == self.length
                and False not in [isinstance(coordinate, numbers.Number) for coordinate in coords]):
            # Can't just assign coords to self.coords because coords could be of type Point and coords should be a tuple
            self.coords = tuple(coords)
        else:
            raise TypeError('coords argument must be a 3 element tuple of Number')

    @property
    def x(self):
        return self.coords[0]

    @property
    def y(self):
        return self.coords[1]

    @property
    def z(self):
        return self.coords[2]

    @property
    def cyl_mag(self):
        return (self.x ** 2 + self.y ** 2) ** 0.5

    @property
    def cyl_theta(self):
        theta = math.acos(self.x / self.cyl_mag)
        if self.y < 0:
            theta = 2 * math.pi - theta
        return theta

    def translate(self, x_offset, y_offset, z_offset):
        return Point((self.coords[0] + x_offset, self.coords[1] + y_offset, self.coords[2] + z_offset))

    def rotate(self, center, axis, angle):
        """ Rotates self about center and axis by angle radians.
        Translates self such that center is at the origin, rotates self by angle radians about axis and then
        translates self such that center is back at its original location"""
        if not (isinstance(center, tuple) and len(center) == self.length):
            raise TypeError('center must be a tuple of length 3')
        if not (isinstance(axis, tuple) and len(axis) == self.length):
            raise TypeError('axis must be a tuple of length 3')
        if not isinstance(angle, numbers.Number):
            raise TypeError('angle must be a number')
        point = self - center
        quat = axisangle_to_q(axis, angle)
        return quat * point + center


    def __add__(self, other):
        """ adds self.coords to other element by element. Returns a Point object"""
        if (isinstance(other, tuple) and len(other) == self.length and
                False not in [isinstance(coord, numbers.Number) for coord in other]):
            return Point(tuple([p_coordinate + o_coordinate
                                for p_coordinate, o_coordinate in zip(self.coords, other)]))
        else:
            raise TypeError('other must have a length of 3 and all elements of other must be of type Number.')

    def __radd__(self, other):
        """ adds self.coords to other element by element. Invoked when Point object is 2nd operand.
        Returns a Point object """
        return self + other

    def __sub__(self, other):
        """ subtracts other from self.coords element by element. Returns a Point object """
        if (isinstance(other, tuple) and len(other) == 3 and
                False not in [isinstance(coord, numbers.Number) for coord in other]):
            return Point(tuple([p_coordinate - o_coordinate
                                # for p_coordinate, o_coordinate in zip(self.coords, other)]))
                                for p_coordinate, o_coordinate in zip(self, other)]))
        else:
            raise TypeError('other must have a length of 3 and all elements of other must be of type Number.')

    def __rsub__(self, other):
        """ subtracts self.coords from  other element by element. Invoked when Point object is 2nd operand.
        Returns a Point object """
        if (isinstance(other, tuple) and len(other) == 3 and
                False not in [isinstance(coord, numbers.Number) for coord in other]):
            return Point(tuple([o_coordinate - p_coordinate
                                for p_coordinate, o_coordinate in zip(self.coords, other)]))
        else:
            raise TypeError('other must have a length of 3 and all elements of other must be of type Number.')

    def __mul__(self, other):
        """ multiplies self.coords by other element by element. Other must be a scalar Number.
        Returns a Point object """
        if isinstance(other, numbers.Number):
            return Point(tuple([p_coordinate * other for p_coordinate in self.coords]))
        else:
            raise TypeError('other must be of type Number.')

    def __rmul__(self, other):
        """ multiplies self.coords by other element by element. Other must be a scalar Number.
        Invoked when Point object is 2nd operand. Returns a Point object """
        return self * other

    def __truediv__(self, other):
        """ divides self.coords by other element by element. Other must be a scalar Number.
        Returns a Point object """
        if isinstance(other, numbers.Number):
            return Point(tuple([p_coordinate / other for p_coordinate in self.coords]))
        else:
            raise TypeError('other must be of type Number.')

    def __neg__(self):
        """ negates each coord"""
        return Point((-self.coords[0], -self.coords[1], -self.coords[2]))

    def __eq__(self, other, verbose=False):
        if isinstance(other, tuple) and len(other) == self.length:
            if False not in [abs(s_coord - o_coord) < tolerance
                             for s_coord, o_coord in zip(self, other)]:
                return True
            else:
                if verbose:
                    print('Point.__eq__ info:', self, other, ' are not equal')
                return False
        elif other is None:
            return False
        else:
            raise TypeError('bad operand for Point.__eq__', other.__class__)

    @property
    def magnitude(self):
        return sum([coord ** 2 for coord in self]) ** 0.5

    @property
    def normalize(self):
        """ return a unit vector parallel to the vector defined by the origin and self"""
        mag = sum(coord * coord for coord in self) ** 0.5
        if mag > tolerance:
            return Point(tuple(coord / mag for coord in self))
        else:
            raise ValueError('Point magnitude too close to zero')
