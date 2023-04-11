""" TgyPar creates stable prism towers from a set of input parameters"""
import math
import numpy as np
import matplotlib.pyplot as plt
import olof_1 as olof

len_tol = 0.1

tendon_type_list = ['bot_top', 'bot_bot', 'top_top', 'waist']


def rotate_list(li, x):
    return li[-x % len(li):] + li[:-x % len(li)]


def distance(p0, p1):
    return sum([(c0 - c1) ** 2 for c0, c1 in zip(p0, p1)]) ** 0.5


def angle(vec0, vec1, acute=True, debug=False):
    """ returns the acute angle between two free vectors """
    # if not isinstance(vec0, Point):
    #     raise TypeError('vec0 must be of type Point')
    # if not isinstance(vec1, Point):
    #     raise TypeError('vec1 must be of type Point')
    # value = math.acos(dotproduct(vec0, vec1) / (vec0.magnitude * vec1.magnitude))
    p0 = np.array(vec0)
    p1 = np.array(vec1)
    cross = np.cross(p0, p1)
    normal = cross / np.linalg.norm(cross)
    # dot = np.dot(p0, p1)
    # norm0 = np.linalg.norm(p0)
    # norm1 = np.linalg.norm(p1)
    # value = math.acos(np.dot(p0, p1) / (np.linalg.norm(p0) * np.linalg.norm(p1)))
    value = math.atan2(np.dot(cross, normal), np.dot(p0, p1))
    if debug:
        print('value from angle', value)
    # value = math.acos(dotproduct(vec0, vec1) / (point0.magnitude * vec1.magnitude))
    # if the angle is greater than pi / 2 then it is the obtuse angle, and we want the acute angle
    if acute and value > math.pi / 2:
        value = math.pi - value
    return value


def midpoint(p0, p1):
    return [(c0 + c1) / 2 for c0, c1 in zip(p0, p1)]


def trilateration(P1, P2, P3, r1, r2, r3):
    p1 = np.array([0, 0, 0])
    p2 = np.array([P2[0] - P1[0], P2[1] - P1[1], P2[2] - P1[2]])
    p3 = np.array([P3[0] - P1[0], P3[1] - P1[1], P3[2] - P1[2]])
    v1 = p2 - p1
    v2 = p3 - p1
    Xn = v1 / np.linalg.norm(v1)
    tmp = np.cross(v1, v2)
    Zn = tmp / np.linalg.norm(tmp)
    Yn = np.cross(Xn, Zn)
    i = np.dot(Xn, v2)
    d = np.dot(Xn, v1)
    j = np.dot(Yn, v2)
    X = ((r1 ** 2) - (r2 ** 2) + (d ** 2)) / (2 * d)
    Y = (((r1 ** 2) - (r3 ** 2) + (i ** 2) + (j ** 2)) / (2 * j)) - ((i / j) * X)
    Z1 = np.sqrt(max(0, r1 ** 2 - X ** 2 - Y ** 2))
    Z2 = -Z1
    K1 = P1 + X * Xn + Y * Yn + Z1 * Zn
    K2 = P1 + X * Xn + Y * Yn + Z2 * Zn
    return K1, K2


def prism_height(strut_theta, strut_length, radii):
    """ Return the height of a prism tensegrity"""
    # Define a triangle whose hypotenuse is the strut, height is the height of the prism and base is the distance
    # between the strut vertex_list projection onto the xy plane
    strut_theta = abs(strut_theta)
    strut_bottom = (radii[0], 0, 0)
    strut_top_on_xy = (radii[1] * math.cos(strut_theta), radii[1] * math.sin(strut_theta), 0)
    base = sum([(p0_coord - p1_coord) ** 2 for p0_coord, p1_coord in zip(strut_bottom, strut_top_on_xy)]) ** 0.5
    return (strut_length ** 2 - base ** 2) ** 0.5


class Hub:
    def __init__(self, length, radius, stitch=None):
        self.length = length
        self.radius = radius
        self.stitch = stitch  # wire allowance for attaching to hub. Includes both ends!


class Vertex:
    def __init__(self, coords, level, hub=None):
        self.coords = coords
        self.level = level
        self.hub = hub
        self.strut = None
        self.tendons = []

    def set_coordinates(self, coords):
        self.coords = coords

    def set_strut(self, strut):
        self.strut = strut

    def add_tendon(self, tendon):
        self.tendons.append(tendon)

    def tendon_radial_vec(self, tendon, debug=False):
        strut_vector = np.array(self.strut.position_vec(self.strut.other_vtx_index(self)))
        tendon_vector = np.array(tendon.position_vec(tendon.other_vtx_index(self)))
        tendon_strut_component = ((np.dot(tendon_vector, strut_vector) / np.linalg.norm(strut_vector) ** 2)
                                  * strut_vector)
        vector = tendon_vector - tendon_strut_component
        if debug:
            return vector, strut_vector, tendon_vector, tendon_strut_component
        else:
            return vector

    def tendon_radial_angle(self, tendon, debug=True):
        """ returns the angle between tendon[0] and tendon after projecting both tendon_list onto the plane that
        passes through this vertex and is perpendicular to strut. This is useful for the design of tendon angles
        on a hub"""
        # normal vector to the plane is the strut
        if tendon is self.tendons[0]:
            # return 0
            normal_component = [0, 0, 0]
            value = 0
        else:
            ref_vector = self.tendon_radial_vec(self.tendons[0])
            vector = self.tendon_radial_vec(tendon)
            value = angle(ref_vector, vector)
            # strut_vector = np.array(self.strut.position_vec(self.strut.other_vertex(self)))
            # ref_tendon_vector = np.array(self.tendon_list[0].position_vec(self))
            # ref_normal_component = (np.dot(ref_tendon_vector, strut_vector) /
            #                         np.linalg.norm(strut_vector) ** 2) * strut_vector
            # ref_vector = ref_tendon_vector - ref_normal_component
            # tendon_vector = np.array(tendon.position_vec(self))
            # normal_component = (np.dot(tendon_vector, strut_vector) / np.linalg.norm(strut_vector) ** 2) * strut_vector
            # vector = tendon_vector - normal_component
            # # if debug:
            # #     print('cross ', np.cross(ref_vector, vector), 'dot cross', np.dot())
            # value = angle(ref_vector, vector, acute=False)
        return value

    @property
    def cyl_coordinates_deg(self):
        """ returns [rho(radii), phi(azimuth), z] """
        x = self.coords[0]
        y = self.coords[1]
        z = self.coords[2]
        r = (x ** 2 + y ** 2) ** 0.5
        theta = math.degrees(math.atan2(y, x))
        return [r, theta, z]

    def distance(self, vertex):
        return (sum([(c0 - c1) ** 2 for c0, c1 in zip(self.coords, vertex.coords)])) ** 0.5


class Strut:
    def __init__(self, vertices, level, targ_length=None):
        self.vertices = vertices
        self.level = level
        self.force = None
        if not targ_length:
            self.targ_length = self.vertices[0].distance(self.vertices[1])

    def set_targ_length(self, targ_length):
        self.targ_length = targ_length

    def set_force(self, force):
        self.force = force

    def other_vertex(self, vertex):
        return self.vertices[(self.vertices.index(vertex) + 1) % 2]

    # @staticmethod
    def position_vec(self, vertex):
        """ returns the position vector (position vectors start at the origin) that points at the vertex arg"""
        return [c0 - c1 for c0, c1 in zip(vertex.coords, self.other_vertex(vertex).coords)]

    @property
    def top_vertex(self):
        return self.vertices[1]

    @property
    def bot_vertex(self):
        return self.vertices[0]

    @property
    def curr_length(self):
        return distance(self.vertices[0].coords, self.vertices[1].coords)

    @property
    def build_length(self):
        return self.curr_length - self.vertices[0].hub.length - self.vertices[1].hub.length


class Tendon:
    def __init__(self, vertices, level=None, tendon_type=None, targ_length=None):
        self.vertices = vertices
        self.level = level
        if tendon_type in tendon_type_list:
            self.tendon_type = tendon_type
        else:
            raise Exception('tendon_type must be in tendon_type_list')
        self.force = None
        if not targ_length:
            self.targ_length = self.vertices[0].distance(self.vertices[1])

    def set_targ_length(self, targ_length):
        self.targ_length = targ_length

    def other_vertex(self, vertex):
        return self.vertices[(self.vertices.index(vertex) + 1) % 2]

    # @staticmethod
    def position_vec(self, vertex):
        """ returns the position vector (position vectors start at the origin) that points at the vertex arg"""
        return [c0 - c1 for c0, c1 in zip(vertex.coords, self.other_vertex(vertex).coords)]

    @property
    def curr_length(self):
        return distance(self.vertices[0].coords, self.vertices[1].coords)

    @property
    def build_length(self):
        angle0 = angle(self.position_vec(self.vertices[0]), self.vertices[0].strut.position_vec(self.vertices[0]))
        offset0 = self.vertices[0].hub.radius * math.cos(angle0)
        angle1 = angle(self.position_vec(self.vertices[1]), self.vertices[1].strut.position_vec(self.vertices[1]))
        offset1 = self.vertices[1].hub.radius * math.cos(angle1)
        return self.curr_length - offset0 - offset1

    @property
    def loose(self):
        return self.curr_length + len_tol < self.targ_length

    def set_force(self, force):
        self.force = force

    def contract(self, d, vertex, constraint0, constraint1):
        """ contract self by d and find the new vertex location while maintaining the constraint distances.
        the constraints can be any member (tendon or strut) """


class Tensegrity:
    def __init__(self, vertices, struts, tendons):
        self.vertices = vertices
        self.struts = struts
        self.tendons = tendons
        self.populate_members()
        self.debug_points = []

    def add_debug_points(self, coords, color='blue'):
        self.debug_points.extend([[coords, color]])

    def populate_members(self):
        """ Populate each vertex's list of tendon_list and strut (members) """
        for strut in self.struts:
            # for vertex in strut.vertex_list:
            for vertex in strut.vtx_coords:
                vertex.set_strut(strut)
        for tendon in self.tendons:
            # for vertex in tendon.vertex_list:
            for vertex in tendon.vtx_coords:
                vertex.add_tendon(tendon)

    def get_vertices(self):
        """ returns a list of all vertex coordinates """
        return [vertex.coords for vertex in self.vertices]

    def get_struts(self):
        """ returns a list of vertex index pairs that define all struts """
        return [[self.vertices.index(strut.vtx_coords[0]),
                 self.vertices.index(strut.vtx_coords[1])] for strut in self.struts]

    def get_tendons(self):
        """ returns a list of vertex index pairs that define all tendons """
        return [[self.vertices.index(tendon.vtx_coords[0]),
                 self.vertices.index(tendon.vtx_coords[1])] for tendon in self.tendons]

    def get_forces(self):
        struts = [[self.vertices.index(strut.vtx_coords[0]),
                   self.vertices.index(strut.vtx_coords[1])] for strut in self.struts]
        tendons = [[self.vertices.index(tendon.vtx_coords[0]),
                    self.vertices.index(tendon.vtx_coords[1])] for tendon in self.tendons]
        vertices = [vertex.coords for vertex in self.vertices]
        forces, termination_type = olof.solve_tensegrity_tensions(np.array(struts), np.array(tendons),
                                                                  np.array(vertices), verbose=False)
        if termination_type:
            for member, force in zip(self.struts + self.tendons, forces):
                member.set_force(force)
        print('termination type', termination_type)

    def plot(self, struts=True, tendons=True, lateral_f=False, vertex_f=False, axes=False, debug=False):
        """ plots tendon_list and strut_list using matplotlib """

        def plot_cross(size=0.2):
            offset = size / 2
            ax.plot3D([point[0] - offset, point[0] + offset], [point[1], point[1]], [point[1], point[1]],
                      color, linewidth=1)
            ax.plot3D([point[0], point[0]], [point[1] - offset, point[1] + offset], [point[1], point[1]],
                      color, linewidth=1)
            ax.plot3D([point[0], point[0]], [point[1], point[1]], [point[1] - offset, point[1] + offset],
                      color, linewidth=1)
            # ax.plot3D(x, y, z, 'blue', linewidth=1)

        # fig = plt.figure(figsize=(4, 4))
        ax = plt.axes(projection='3d')
        # ax.set_xlim(-5, 5)
        # ax.set_ylim(-5, 5)
        # ax.set_zlim(0, 10)
        x_span = 2*max(vertex.coords[0] for vertex in self.vertices) - 2*min(vertex.coords[0] for vertex in self.vertices)
        y_span = 2*max(vertex.coords[1] for vertex in self.vertices) - 2*min(vertex.coords[1] for vertex in self.vertices)
        z_span = 2*max(vertex.coords[2] for vertex in self.vertices) - 2*min(vertex.coords[2] for vertex in self.vertices)
        span = max(x_span, y_span, z_span)
        ax.set_xlim(-span / 2, span / 2)
        ax.set_ylim(-span / 2, span / 2)
        ax.set_zlim(0, span)
        if not axes:
            ax.set_axis_off()
        if struts:
            for strut in self.struts:
                x = [vertex.coords[0] for vertex in strut.vtx_coords]
                y = [vertex.coords[1] for vertex in strut.vtx_coords]
                z = [vertex.coords[2] for vertex in strut.vtx_coords]
                ax.plot3D(x, y, z, 'red', linewidth=3)
        if tendons:
            for tendon in self.tendons:
                x = [vertex.coords[0] for vertex in tendon.vtx_coords]
                y = [vertex.coords[1] for vertex in tendon.vtx_coords]
                z = [vertex.coords[2] for vertex in tendon.vtx_coords]
                ax.plot3D(x, y, z, 'grey', linewidth=1)
        if debug:
            # for vertex in self.vtx_coords:
            strut = self.struts[0]
            vertex = strut.vtx_coords[0]
            tendon = vertex.tendon_list[0]
            strut_vector = np.array(strut.position_vec(strut.other_vtx_index(vertex)))
            # ax.quiver(*vertex.coords, *strut_vector, color='orange', linewidth=3)
            coords = [val for pair in zip(vertex.coords, strut_vector) for val in pair]
            x = [vertex.coords[0] for vertex in [vertex.coords, strut_vector]]
            y = [vertex.coords[1] for vertex in [vertex.coords, strut_vector]]
            z = [vertex.coords[2] for vertex in [vertex.coords, strut_vector]]
            # y = [vertex.coords[1] for vertex in tendon.vtx_coords]
            # z = [vertex.coords[2] for vertex in tendon.vtx_coords]
            ax.plot3D(x, y, z, 'orange', linewidth=5)
                      # *self.strut_list[0].position_vec(self.strut_list[0].vtx_coords[1]), 'green', linewidth=3)
            tendon_vector = np.array(tendon.position_vec(tendon.other_vtx_index(vertex)))
            # ax.quiver(*vertex.coords, *tendon_vector, color='orange', linewidth=3)
            tendon_strut_component = ((np.dot(tendon_vector, strut_vector) / np.linalg.norm(strut_vector) ** 2)
                                      * strut_vector)
            # ax.quiver(*vertex.coords, *tendon_strut_component, color='green', linewidth=3)
            tendon_lateral_component = tendon_vector - tendon_strut_component
            # ax.quiver(*vertex.coords, *tendon_lateral_component, color='blue', linewidth=3)

            for vertex in self.vertices[0:1]:
                for tendon in vertex.tendon_list[0:1]:
                    radial_vec, strut_vec, tendon_vec, tendon_strut_component = vertex.tendon_radial_vec(tendon, debug=True)
                    # ax.quiver(0, 0, 0, 10, 10, 10, 'green', linewidth=2)
                    # ax.plot3D([vertex.coords[0], vertex.tendon_radial_vec(tendon)[0]],
                    #           [vertex.coords[1], vertex.tendon_radial_vec(tendon)[1]],
                    #           [vertex.coords[2], vertex.tendon_radial_vec(tendon)[2]],
                    #           'blue', linewidth=1)
                    # ax.plot3D([vertex.coords[0], radial_vec[0]],
                    #           [vertex.coords[1], radial_vec[1]],
                    #           [vertex.coords[2], radial_vec[2]],
                    #           'blue', linewidth=1)
                    # ax.plot3D([vertex.coords[0], tendon_vec[0]],
                    #           [vertex.coords[1], tendon_vec[1]],
                    #           [vertex.coords[2], tendon_vec[2]],
                    #           'orange', linewidth=1)
                    # ax.plot3D([vertex.coords[0], tendon_strut_component[0]],
                    #           [vertex.coords[1], tendon_strut_component[1]],
                    #           [vertex.coords[2], tendon_strut_component[2]],
                    #           'yellow', linewidth=1)
                    # ax.plot3D([vertex.coords[0], strut_vec[0]],
                    #           [vertex.coords[1], strut_vec[1]],
                    #           [vertex.coords[2], strut_vec[2]],
                    #           'orange', linewidth=1)

            # offset = 0.1
            # for debug_point in self.debug_points:
            #     point, color = debug_point
            #     ax.plot3D([point[0] - offset, point[0] + offset], [point[1], point[1]], [point[2], point[2]],
            #               color, linewidth=1)
            #     ax.plot3D([point[0], point[0]], [point[1] - offset, point[1] + offset], [point[2], point[2]],
            #               color, linewidth=1)
            #     ax.plot3D([point[0], point[0]], [point[1], point[1]], [point[2] - offset, point[2] + offset],
            #               color, linewidth=1)
            #     # plot_cross()
        plt.show()

    def print_build(self):
        precision = 3
        np.set_printoptions(formatter={'float': '{: 10.3f}'.format})
        label_width = len("Element")
        number_width = len("Finish Length")
        # number_width = 11
        print('*** Prism Tower Build Data ***')
        print('n =', self.n)
        print('levels =', self.levels)
        print('radii =', np.array(self.radii))
        print('heights =', np.array(self.heights))
        print(f'{"Element": <{label_width}}',
              f'{"Level": <{number_width}}',
              f'{"Index": <{number_width}}',
              f'{"Finish Length": <{number_width}}',
              f'{"Build Length": <{number_width}}')
        for index, strut in enumerate(self.struts):
            print(f'{"Strut": <{label_width}}',
                  f'{strut.level: <{number_width}}',
                  # f'{index % (self.levels - 1): <{number_width}}',
                  f'{index % self.n: <{number_width}}',
                  f'{str(round(strut.curr_length, precision)): <{number_width}}',
                  f'{str(round(strut.build_length, precision)): <{number_width}}')
        print(f'{"Element": <{label_width}}',
              f'{"Type": <{label_width}}',
              f'{"Level": <{number_width}}',
              f'{"Index": <{number_width}}',
              f'{"Finish Length": <{number_width}}',
              f'{"Build Length": <{number_width}}')
        for index, tendon in enumerate(self.tendons):
            print(f'{"Tendon": <{label_width}}',
                  f'{tendon.tendon_type: <{label_width}}',
                  f'{tendon.level: <{number_width}}',
                  # todo the line below incorrectly labels the index for waist tendon_list that are not at the very top or
                  # bottom
                  f'{index % self.n: <{number_width}}',
                  f'{str(round(tendon.curr_length, precision)): <{number_width}}',
                  f'{str(round(tendon.build_length, precision)): <{number_width}}')
        print(f'{"Element": <{label_width}}',
              f'{"Level": <{number_width}}',
              f'{"Index": <{number_width}}',
              f'{"Tendon Angles": <{number_width}}')
        for index, vertex in enumerate(self.vertices):
            t_angles = []
            for tendon in vertex.tendon_list:
                # tangle, normal_component = vertex.tendon_radial_angle(tendon)
                tangle = vertex.tendon_radial_angle(tendon)
                t_angles.append(round(math.degrees(tangle)))
            print(f'{"Vertex": <{label_width}}',
                  f'{vertex.level: <{number_width}}',
                  f'{index % self.n: <{number_width}}',
                  f'{t_angles}'
                  )

    def print_cyl(self, vertices=True, struts=True, tendons=True):
        np.set_printoptions(formatter={'float': '{: 10.3f}'.format})
        # type_width = len('vertex')
        coord_width = 10  # make this equal to the width specifier in set_printoptions call above
        array_3_width = 3 * coord_width + 4  # + 4 to account for two space and two brackets
        number_width = 11
        precision = 3
        label_vertex = 'Vertex'
        label_tendon = 'Tendon'
        label_strut = 'Strut'
        label_element = 'Element'
        label_width = max([len(e)
                           for e in [label_vertex, label_tendon, label_strut, label_element] + tendon_type_list]) + 1
        heading_coord = 'Coordinates'
        heading_v0 = f'V0 {heading_coord}'
        heading_v1 = f'V1 {heading_coord}'
        if vertices:
            print(f'{"Element": <{label_width}}',
                  f'{heading_coord: <{array_3_width}}',
                  f'{"Force": <{array_3_width}}')
            print(f'{" ": <{label_width}}',
                  f'{"radii   ": >{coord_width}}',
                  f'{"theta": >{coord_width}}',
                  f'{"z   ": >{coord_width}}',
                  " ",
                  f'{"radii   ": >{coord_width}}',
                  f'{"theta": >{coord_width}}',
                  f'{"z   ": >{coord_width}}')
            for vertex in self.vertices:
                print(f'{label_vertex: <{label_width}}',
                      f'{str(np.array(vertex.cyl_coordinates_deg)): <{array_3_width}}',
                      # f'{str(np.around(np.array(vertex.f_vector), precision)): <{array_3_width}}'
                      )
        if tendons:
            print(f'{"Element": <{label_width}}',
                  f'{heading_v0: <{array_3_width}}',
                  f'{heading_v1: <{array_3_width}}',
                  f'{"type": <{label_width}}',
                  f'{"targ len": <{number_width}}',
                  f'{"curr len": <{number_width}}',
                  # f'{heading_v1: <{array_3_width}}',
                  # f'{"Force": <{array_3_width}}'
                  )
            print(f'{" ": <{label_width}}',
                  f'{"radii   ": >{coord_width}}',
                  f'{"theta": >{coord_width}}',
                  f'{"z   ": >{coord_width}}',
                  " ",
                  f'{"radii   ": >{coord_width}}',
                  f'{"theta": >{coord_width}}',
                  f'{"z   ": >{coord_width}}')
            for tendon in self.tendons:
                vertex0 = tendon.vtx_coords[0]
                vertex1 = tendon.vtx_coords[1]
                print(f'{"Tendon": <{label_width}}',
                      f'{str(np.array(vertex0.cyl_coordinates_deg)): <{array_3_width}}',
                      f'{str(np.array(vertex1.cyl_coordinates_deg)): <{array_3_width}}',
                      f'{tendon.tendon_type: <{label_width}}'
                      f'{str(round(tendon.targ_length, precision)): <{number_width}}'
                      f'{str(round(tendon.curr_length, precision)): <{number_width}}'
                      # f'{str(round(tendon.spring_f_mag, precision)): <{number_width}}'
                      )
        if struts:
            print(f'{"Element": <{label_width}}',
                  f'{heading_v0: <{array_3_width}}',
                  f'{heading_v1: <{array_3_width}}',
                  # f'{"Twist": <{number_width}}',
                  # f'{"Force Vector": <{array_3_width}}',
                  )
            print(f'{" ": <{label_width}}',
                  f'{"radii   ": >{coord_width}}',
                  f'{"theta": >{coord_width}}',
                  f'{"z   ": >{coord_width}}',
                  " ",
                  f'{"radii   ": >{coord_width}}',
                  f'{"theta": >{coord_width}}',
                  f'{"z   ": >{coord_width}}')
            for strut in self.struts:
                vertex0_cyl = strut.vtx_coords[0].cyl_coordinates_deg
                vertex1_cyl = strut.vtx_coords[1].cyl_coordinates_deg
                print(f'{"Strut": <{label_width}}',
                      f'{str(np.array(vertex0_cyl)): <{array_3_width}}',
                      f'{str(np.array(vertex1_cyl)): <{array_3_width}}',
                      # f'{str(round(vertex0_cyl[1] - vertex1_cyl[1], precision)): <{number_width}}',
                      # f'{str(np.array(strut.spring_f_vec)): <{array_3_width}}'
                      )


# class PrismTower(Tensegrity):
#     def __init__(self, n=3, levels=2, radii=[3, 3], heights=[10, 10], overlaps=0.06, verbose=0):
#         """ return 3 numpy.arrays that describe a stable tensegrity structure:
#         vertex coordinates, tendon vertex_list, strut vertex_list """
#         # Assume that the bottom vertex_list lie on the x,y plane
#         # bottom vertex_list are fully described by the input parameters n and radii[0]
#         self.n = n
#         self.levels = levels
#         temp0 = [heights[0:i] for i in range(1, len(heights) - 1)]
#         print('temp0', temp0)
#         # temp = [[sum([heights[0:i]] for i in range(len(heights))]
#         # bot_z_list = [0] + [[sum([heights[0:i]])] for i in range(len(heights))]
#         # for level, height in enumerate(heights):
#         for level in range(levels):
#             # bot_z = bot_z_list[level]
#             bot_z = 0
#             top_z = heights[0]
#             # top_z = bot_z + height
#             theta_offset = 0
#             # twist_sign = 1 will make a left-handed twist (right-handed = strut_list look like right-handed threads)
#             twist_sign = -1
#             twist = twist_sign * (math.pi / 2 - math.pi / self.n)
#             theta_step = 2 * math.pi / self.n
#             # todo create multi level element lists
#             self.bot_vertices = [Vertex(np.array([radii[0] * math.cos(i * theta_step),
#                                                   radii[0] * math.sin(i * theta_step),
#                                                   bot_z], level)) for i in range(self.n)]
#             self.top_vertices = [Vertex(np.array([radii[1] * math.cos(i * theta_step - twist),
#                                                   radii[1] * math.sin(i * theta_step - twist),
#                                                   top_z], level))]
#             self.bot_tendons = [Tendon([v0, v1]) for v0, v1 in zip(self.bot_vertices, rotate_list(self.bot_vertices, 1))]
#             self.top_tendons = [Tendon([v0, v1]) for v0, v1 in zip(self.top_vertices, rotate_list(self.top_vertices, 1))]
#             self.vert_tendons = [Tendon([v0, v1]) for v0, v1 in zip(self.bot_vertices, self.top_vertices)]
#
#         strut_list = [Strut([v0, v1]) for v0, v1 in zip(self.bot_vertices, rotate_list(self.top_vertices, twist_sign))]
#         vertex_list = self.bot_vertices + self.top_vertices
#         tendon_list = self.bot_tendons + self.top_tendons + self.vert_tendons
#         Tensegrity.__init__(self, vertex_list, strut_list, tendon_list)

class PrismTower(Tensegrity):
    def __init__(self, n=3, levels=2, radii=[4, 3, 4], heights=[10, 10], strut_lengths=None, overlaps=[0.2],
                 hub_length=0, hub_radius=8, verbose=0):
        """ return 3 numpy.arrays that describe a stable tensegrity structure:
        vertex coordinates, tendon vertex_list, strut vertex_list
        if strut_lengths are given, then heights are calculated from the strut_lengths and heights are ignored."""
        if len(radii) < levels + 1:
            raise Exception('the length of radii must be greater than levels + 1')
        if len(heights) < levels and not strut_lengths:
            raise Exception('levels must be equal to or less than the length of heights')
        if len(overlaps) < levels - 1:
            raise Exception('the length of overlaps must 1 less than the length of heights or greater')
        if strut_lengths and len(strut_lengths) < levels:
            raise Exception('levels must be equal to or less than the length of strut_lengths')
        if strut_lengths and len(strut_lengths) + 1 > len(radii):
            raise Exception('The length of strut_lengths must be one greater or more than the length of radii')
        # redundant_tendons = True
        redundant_tendons = False
        self.n = n
        self.levels = levels
        self.radii = radii
        self.overlaps = overlaps
        hub = Hub(hub_length, hub_radius)
        vertices = []
        struts = []
        tendons = []
        # twist_sign = 1 will make a left-handed twist (right-handed = strut_list look like right-handed threads)
        twist_sign = -1
        twist = twist_sign * (math.pi / 2 - math.pi / n)
        if strut_lengths:
            self.heights = [prism_height(abs(twist) + 2 * math.pi / n, length, [self.radii[i], self.radii[i+1]])
                            for i, length in enumerate(strut_lengths)]
        else:
            self.heights = heights
        theta_step = 2 * math.pi / n
        theta_offset = 0
        if levels > 1:
            bot_z_list = [0] + [sum(self.heights[0:i]) * (1 - self.overlaps[i - 1]) for i in range(1, levels)]
        else:
            bot_z_list = [0]
        for level, radius, height in zip(range(levels), self.radii, self.heights):
            # First we create vertex_list and strut_list, they fit neatly into a layer
            bot_z = bot_z_list[level]
            bot_radius = self.radii[level]
            top_radius = self.radii[level + 1]
            top_z = bot_z + height
            bot_vertices = [Vertex(np.array([bot_radius * math.cos(i * theta_step + theta_offset),
                                             bot_radius * math.sin(i * theta_step + theta_offset),
                                             bot_z]), level=level, hub=hub)
                            for i in range(n)]
            top_vertices = [Vertex(np.array([top_radius * math.cos(i * theta_step - twist + theta_offset),
                                             top_radius * math.sin(i * theta_step - twist + theta_offset),
                                             top_z]), level=level, hub=hub)
                            for i in range(n)]
            vertices.extend(bot_vertices + top_vertices)
            # Always put bot vertex in strut.vertex_list[0]
            struts.extend([Strut([v0, v1], level=level) for v0, v1 in
                           zip(bot_vertices, rotate_list(top_vertices, twist_sign))])
            if redundant_tendons:
                tendons.extend([Tendon([v0, v1], level=level, tendon_type='bot_top') for v0, v1 in
                               zip(bot_vertices, rotate_list(top_vertices, -twist_sign))])
                tendons.extend([Tendon([v0, v1], level=max(v0.level, v1.level), tendon_type='waist') for v0, v1 in
                                zip(bot_vertices, rotate_list(bot_vertices, 1))])
                tendons.extend([Tendon([v0, v1], level=max(v0.level, v1.level), tendon_type='waist') for v0, v1 in
                                zip(top_vertices, rotate_list(top_vertices, 1))])
            theta_offset += math.pi / n - twist
        # Tendons: now let's create the tendon_list, they can bridge layers
        # bottom tendon_list
        if not redundant_tendons:
            vertex_list = [strut.bot_vertex for strut in [strut for strut in struts if strut.level == 0]]
            tendons.extend([Tendon([v0, v1], level=max(v0.level, v1.level), tendon_type='waist') for
                            v0, v1 in zip(vertex_list, rotate_list(vertex_list, 1))])
            # top tendon_list
            vertex_list = [strut.top_vertex for strut in [strut for strut in struts if strut.level == self.levels - 1]]
            tendons.extend([Tendon([v0, v1], level=self.levels, tendon_type='waist') for
                            v0, v1 in zip(vertex_list, rotate_list(vertex_list, 1))])
        for level in range(levels - 1):
            # mid tendon_list (waist tendon_list that are not bottom or top, they connect the top of each layer to bottom of
            # the layer above)
            lower_vertex_list = [strut.top_vertex for strut in [strut for strut in struts if strut.level == level]]
            upper_vertex_list = [strut.bot_vertex for strut in [strut for strut in struts if strut.level == level + 1]]
            tendons.extend([Tendon([v0, v1], level=max(v0.level, v1.level), tendon_type='waist')
                            for v0, v1 in zip(lower_vertex_list, rotate_list(upper_vertex_list, -1))])
            tendons.extend([Tendon([v0, v1], level=max(v0.level, v1.level), tendon_type='waist')
                            for v0, v1 in zip(rotate_list(lower_vertex_list, -1), rotate_list(upper_vertex_list, -1))])
            # vertical tendon_list top to top
            lower_vertex_list = [strut.top_vertex for strut in [strut for strut in struts if strut.level == level]]
            upper_vertex_list = [strut.top_vertex for strut in [strut for strut in struts if strut.level == level + 1]]
            tendons.extend([Tendon([v0, v1], level=level, tendon_type='top_top')
                            for v0, v1 in zip(lower_vertex_list, rotate_list(upper_vertex_list, 1))])
            # vertical tendon_list bottom to bottom
            # if level > 1:
            if redundant_tendons:
                lower_vertex_list = [strut.bot_vertex for strut in [strut for strut in struts if strut.level == level]]
                upper_vertex_list = [strut.bot_vertex
                                     for strut in [strut for strut in struts if strut.level == level + 1]]
                tendons.extend([Tendon([v0, v1], level=level, tendon_type='bot_bot')
                                for v0, v1 in zip(lower_vertex_list, rotate_list(upper_vertex_list, 1))])
        # vertical tendon_list bottom to top
        lower_vertex_list = [strut.bot_vertex for strut in [strut for strut in struts if strut.level == 0]]
        upper_vertex_list = [strut.top_vertex for strut in [strut for strut in struts if strut.level == 0]]
        tendons.extend([Tendon([v0, v1], level=0, tendon_type='bot_top')
                        for v0, v1 in zip(lower_vertex_list, rotate_list(upper_vertex_list, 1))])
        self.vertices = vertices
        self.tendons = tendons
        self.struts = struts
        Tensegrity.__init__(self, vertices, struts, tendons)
        # ** End __init__ for PrismTower

    # def strut_list(self, level=None):
    #     if level:
    #         return [strut for strut in self.strut_list if strut.level == level]
    #     else:
    #         return self.strut_list

    def set_overlap(self, level=0):
        """ set overlaps such that waist tendon_list lie on same plane as strut
        Assumes symmetrical (not bent) tower prism
        level references the lowest of the two interfacing levels """
        # for strut in self.strut_list:
        #     print('strut level', strut.level)
        #     if strut.level == level + 1:
        #         waist_tendons = [tendon for tendon in strut.bot_vertex.tendon_list if tendon.tendon_type == 'waist']
        #         waist_vectors = [tendon.position_vec(strut.bot_vertex) for tendon in waist_tendons]
        #         waist_cross = np.cross(np.array(waist_vectors[0]), np.array(waist_vectors[1]))
        #         overlap_cos = (np.dot(waist_cross, np.array(strut.position_vec(strut.bot_vertex))) /
        #                        (np.linalg.norm(waist_cross) * np.linalg.norm(strut.position_vec(strut.bot_vertex))))
        #         overlap_angle = math.acos(overlap_cos)
        #         print('overlap angle', math.degrees(overlap_angle))
        for strut in [strut for strut in self.struts if strut.level == level + 1]:
            waist_tendons = [tendon for tendon in strut.bot_vertex.tendon_list if tendon.tendon_type == 'waist']
            waist_vectors = [tendon.position_vec(strut.bot_vertex) for tendon in waist_tendons]
            waist_cross = np.cross(np.array(waist_vectors[0]), np.array(waist_vectors[1]))
            overlap_cos = (np.dot(waist_cross, np.array(strut.position_vec(strut.bot_vertex))) /
                           (np.linalg.norm(waist_cross) * np.linalg.norm(strut.position_vec(strut.bot_vertex))))
            overlap_angle = math.acos(overlap_cos)
            print('overlap angle', math.degrees(overlap_angle))

    def stabilize(self):
        if self.levels == 1:
            pass
            """ Assumes we start with vertical strut_list """
            # Confirm that there are two intersections b


class KitePar(Tensegrity):
    """ Demonstrate trilateration method of constructing a stable tensegrity"""
    # def __init__(self, t_len=2**0.5/2, s_len=1):
    def __init__(self, t_len=1.5, s_len=1):
        bot_vertices = [Vertex(c, level=1) for c in [[-1, -1, 0], [1, -1, 0]]]
        top_vertices = [Vertex(c, level=1) for c in [[0.7, 0.7, 0], [-0.7, 0.7, 0]]]
        struts = [Strut([v0, v1], level=1) for v0, v1 in [[bot_vertices[0], top_vertices[0]],
                                                          [bot_vertices[1], top_vertices[1]]]]
        vertices = bot_vertices + top_vertices
        tendons = [Tendon([v0, v1]) for v0, v1 in zip(vertices, rotate_list(vertices, 1))]
        for tendon in tendons:
            tendon.set_targ_length(t_len)
        for strut in struts:
            strut.set_targ_length(s_len)
        Tensegrity.__init__(self, vertices, struts, tendons)

    @property
    def loose(self):
        return True in [tendon.loose for tendon in self.tendons]

    def stretch_struts(self, max_steps=3):
        """ lengthen the strut_list until tendon_list are completely tight """
        d_tol = len_tol / 10
        step_count = 0
        strut_step = 0.2
        done = False
        # while self.loose and step_count < max_steps:
        while not done and step_count < max_steps:
            for strut in self.struts:
            # for strut in self.strut_list[:1]:
                vertex = strut.top_vertex  # this is the vertex that is moving
                tendons = vertex.tendon_list
                #  we use tendon target lengths and strut current lengths to support our algorithmic strategy
                trial_positions = trilateration(strut.other_vtx_index(vertex).coords,
                                                tendons[0].other_vtx_index(vertex).coords,
                                                tendons[1].other_vtx_index(vertex).coords,
                                                strut.curr_length, tendons[0].targ_length, tendons[1].targ_length
                                                )
                vertex.coords = midpoint(trial_positions[0], trial_positions[1])
                # if distance(trial_positions[0], trial_positions[1]) < d_tol:
                #     done = True
                # print('trial positions', trial_positions)
                self.add_debug_points(trial_positions[0])
                self.add_debug_points(trial_positions[1])
                self.add_debug_points(vertex.coords, color='black')
                # print('strut.other_vertex(vertex).coords', strut.other_vertex(vertex).coords, strut.length)
                # self.add_debug_points(vertex.coords, color='green')
                # print('vertex.coords', vertex.coords)
                self.add_debug_points(tendons[0].other_vtx_index(vertex).coords, color='green')
                # print('tendon_list[0].other_vertex(vertex).coords', tendon_list[0].other_vertex(vertex).coords, tendon_list[0].length)
                self.add_debug_points(tendons[1].other_vtx_index(vertex).coords, color='green')
                # print('tendon_list[1].other_vertex(vertex).coords', tendon_list[1].other_vertex(vertex).coords, tendon_list[1].length)
            for strut in self.struts:
                strut.set_targ_length(strut.curr_length + strut_step)
            step_count += 1
        print('step count', step_count)


if __name__ == '__main__':
    # prism_tower = True
    prism_tower = False
    # prism_1_s_len = False
    prism_1_s_len = True
    # prism_1_tower = True
    prism_1_tower = False
    # prism_2_tower = True
    prism_2_tower = False
    # bojum_tower = True
    bojum_tower = False
    # kite_par = True
    kite_par = False
    if kite_par:
        # Note: KitePar() is under construction and not ready for show time!
        kite = KitePar()
        # for strut in kite.strut_list:
        #     strut.set_length(strut.length * 0.9)
        print('loose is', kite.loose)
        kite.print_cyl(vertices=False)
        kite.stretch_struts()
        kite.print_cyl(vertices=False)
        kite.plot(debug=True)
    if prism_1_tower:
        """ For levels = 1 PrismTower creates a Tensegrity that Olof's alglib code can balance 
            For levels > 1 PrismTower creates 'reasonable' structures, but the interlayer overlaps is just a guess
            provided by the user, so multi-level PrismTowers can't be balanced by Olof's alglib code"""
        tower = PrismTower(n=3, levels=1, radii=[4, 4], heights=[8], verbose=True)
        tower.get_forces()
        tower.print_build()
        tower.set_overlap(0)
        tower.plot()
    if prism_1_s_len:  # specify strut_length instead of height
        """ For levels = 1 PrismTower creates a Tensegrity that Olof's alglib code can balance 
            For levels > 1 PrismTower creates 'reasonable' structures, but the interlayer overlaps is just a guess
            provided by the user, so multi-level PrismTowers can't be balanced by Olof's alglib code"""
        hub_len = 10
        finish_strut_len = 17.5 * 25.4 + 2 * hub_len
        tower = PrismTower(n=3, levels=1, strut_lengths=[finish_strut_len],
                           radii=[finish_strut_len/2.78, finish_strut_len/2.78], heights=[999, 999],
                           overlaps=[0.26], hub_length=10, hub_radius=4.25, verbose=True)
        tower.get_forces()
        tower.print_build()
        tower.set_overlap(0)
        tower.plot(struts=False, debug=True)
    if prism_2_tower:
        """ For levels = 1 PrismTower creates a Tensegrity that Olof's alglib code can balance 
            For levels > 1 PrismTower creates 'reasonable' structures, but the interlayer overlaps is just a guess
            provided by the user, so multi-level PrismTowers can't be balanced by Olof's alglib code"""
        tower = PrismTower(n=3, levels=2, radii=[4, 4, 4], heights=[8, 8], overlaps=[0.26], verbose=True)
        tower.get_forces()
        tower.print_cyl()
        tower.print_build()
        tower.set_overlap(0)
        tower.plot()
    if prism_tower:
        """ For levels = 1 PrismTower creates a Tensegrity that Olof's alglib code can balance 
            For levels > 1 PrismTower creates 'reasonable' structures, but the interlayer overlaps is just a guess
            provided by the user, so multi-level PrismTowers can't be balanced by Olof's alglib code"""
        # tower = PrismTower(n=3, levels=3, radii=[4, 4, 4, 4, 3], heights=[8, 8, 8, 8, 8], overlaps=[0.268, 0.268],
        #                    verbose=True)
        tower = PrismTower(n=4, levels=4, radii=[4, 4, 4, 4, 3], heights=[8, 8, 8, 8, 8],
                           overlaps=[0.172, 0.172, 0.1785], verbose=True)
        tower = PrismTower(n=4, levels=6, radii=[4, 4, 4, 4, 4, 4, 4], heights=[8, 8, 8, 8, 8, 8, 8],
                           overlaps=[0.172, 0.172, 0.1785, 0.1785, 0.1785, 0.1785], verbose=True)
        tower.get_forces()
        tower.print_cyl()
        tower.set_overlap(0)
        tower.set_overlap(1)
        tower.set_overlap(2)
        tower.plot()
    if bojum_tower:
        tower = PrismTower(n=4, levels=5, radii=[6, 5, 2, 2, 2, 5], heights=[8, 8, 8, 8, 8], overlaps=[0, 0, 0, 0], verbose=True)
        # tower.print_cyl()
        tower.print_build()
        tower.plot(axes=False)

