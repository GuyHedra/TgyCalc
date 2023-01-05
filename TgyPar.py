""" TgyPar creates stable prism towers from a set of input parameters"""
import math
import numpy as np
import matplotlib.pyplot as plt
import olof_1 as olof


def rotate_list(li, x):
    return li[-x % len(li):] + li[:-x % len(li)]


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


class Vertex:
    def __init__(self, coords, level):
        self.coords = coords
        self.level = level

    @property
    def cyl_coordinates_deg(self):
        """ returns [rho(r), phi(azimuth), z] """
        x = self.coords[0]
        y = self.coords[1]
        z = self.coords[2]
        r = (x ** 2 + y ** 2) ** 0.5
        theta = math.degrees(math.atan2(y, x))
        return [r, theta, z]


class Strut:
    def __init__(self, vertices, level):
        self.vertices = vertices
        self.level = level
        self.force = None

    def set_force(self, force):
        self.force = force


class Tendon:
    def __init__(self, vertices):
        self.vertices = vertices
        self.force = None

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

    def get_forces(self):
        struts = [[self.vertices.index(strut.vertices[0]),
                   self.vertices.index(strut.vertices[1])] for strut in self.struts]
        tendons = [[self.vertices.index(tendon.vertices[0]),
                   self.vertices.index(tendon.vertices[1])] for tendon in self.tendons]
        vertices = [vertex.coords for vertex in self.vertices]
        forces, termination_type = olof.solve_tensegrity_tensions(np.array(struts), np.array(tendons),
                                                                  np.array(vertices), verbose=False)
        if termination_type:
            for member, force in zip(self.struts + self.tendons, forces):
                member.set_force(force)
        print('termination type', termination_type)

    def plot(self, struts=True, tendons=True, lateral_f=False, vertex_f=False):
        """ plots tendons and struts using matplotlib """
        ax = plt.axes(projection='3d')
        # ax.set_xlim(-5, 5)
        # ax.set_ylim(-5, 5)
        # ax.set_zlim(0, 10)
        x_span = max(vertex.coords[0] for vertex in self.vertices) - min(vertex.coords[0] for vertex in self.vertices)
        y_span = max(vertex.coords[1] for vertex in self.vertices) - min(vertex.coords[1] for vertex in self.vertices)
        z_span = max(vertex.coords[2] for vertex in self.vertices) - min(vertex.coords[2] for vertex in self.vertices)
        span = max(x_span, y_span, z_span)
        ax.set_xlim(-span / 2, span / 2)
        ax.set_ylim(-span / 2, span / 2)
        ax.set_zlim(0, span)
        ax.set_axis_off()
        if struts:
            for strut in self.struts:
                x = [vertex.coords[0] for vertex in strut.vertices]
                y = [vertex.coords[1] for vertex in strut.vertices]
                z = [vertex.coords[2] for vertex in strut.vertices]
                ax.plot3D(x, y, z, 'red', linewidth=3)
        if tendons:
            for tendon in self.tendons:
                x = [vertex.coords[0] for vertex in tendon.vertices]
                y = [vertex.coords[1] for vertex in tendon.vertices]
                z = [vertex.coords[2] for vertex in tendon.vertices]
                ax.plot3D(x, y, z, 'grey', linewidth=1)
        plt.show()

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
        label_width = max([len(e) for e in [label_vertex, label_tendon, label_strut, label_element]])
        heading_coord = 'Coordinates'
        heading_v0 = f'V0 {heading_coord}'
        heading_v1 = f'V1 {heading_coord}'
        if vertices:
            print(f'{"Element": <{label_width}}',
                  f'{heading_coord: <{array_3_width}}',
                  f'{"Force": <{array_3_width}}')
            print(f'{" ": <{label_width}}',
                  f'{"r   ": >{coord_width}}',
                  f'{"theta": >{coord_width}}',
                  f'{"z   ": >{coord_width}}',
                  " ",
                  f'{"r   ": >{coord_width}}',
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
                  # f'{"Force": <{array_3_width}}'
                  )
            print(f'{" ": <{label_width}}',
                  f'{"r   ": >{coord_width}}',
                  f'{"theta": >{coord_width}}',
                  f'{"z   ": >{coord_width}}',
                  " ",
                  f'{"r   ": >{coord_width}}',
                  f'{"theta": >{coord_width}}',
                  f'{"z   ": >{coord_width}}')
            for tendon in self.tendons:
                vertex0 = tendon.vertices[0]
                vertex1 = tendon.vertices[1]
                print(f'{"Tendon": <{label_width}}',
                      f'{str(np.array(vertex0.cyl_coordinates_deg)): <{array_3_width}}',
                      f'{str(np.array(vertex1.cyl_coordinates_deg)): <{array_3_width}}',
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
                  f'{"r   ": >{coord_width}}',
                  f'{"theta": >{coord_width}}',
                  f'{"z   ": >{coord_width}}',
                  " ",
                  f'{"r   ": >{coord_width}}',
                  f'{"theta": >{coord_width}}',
                  f'{"z   ": >{coord_width}}')
            for strut in self.struts:
                vertex0_cyl = strut.vertices[0].cyl_coordinates_deg
                vertex1_cyl = strut.vertices[1].cyl_coordinates_deg
                print(f'{"Strut": <{label_width}}',
                      f'{str(np.array(vertex0_cyl)): <{array_3_width}}',
                      f'{str(np.array(vertex1_cyl)): <{array_3_width}}',
                      # f'{str(round(vertex0_cyl[1] - vertex1_cyl[1], precision)): <{number_width}}',
                      # f'{str(np.array(strut.spring_f_vec)): <{array_3_width}}'
                      )


# class PrismTower(Tensegrity):
#     def __init__(self, n=3, levels=2, radii=[3, 3], heights=[10, 10], overlap=0.06, verbose=0):
#         """ return 3 numpy.arrays that describe a stable tensegrity structure:
#         vertex coordinates, tendon vertices, strut vertices """
#         # Assume that the bottom vertices lie on the x,y plane
#         # bottom vertices are fully described by the input parameters n and radii[0]
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
#             # twist_sign = 1 will make a left-handed twist (right-handed = struts look like right-handed threads)
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
#         struts = [Strut([v0, v1]) for v0, v1 in zip(self.bot_vertices, rotate_list(self.top_vertices, twist_sign))]
#         vertices = self.bot_vertices + self.top_vertices
#         tendons = self.bot_tendons + self.top_tendons + self.vert_tendons
#         Tensegrity.__init__(self, vertices, struts, tendons)

class PrismTower(Tensegrity):
    def __init__(self, n=3, levels=2, radii=[3, 3], heights=[10, 10], overlap=0.2, verbose=0):
        """ return 3 numpy.arrays that describe a stable tensegrity structure:
        vertex coordinates, tendon vertices, strut vertices """
        # Assume that the bottom vertices lie on the x,y plane
        # bottom vertices are fully described by the input parameters n and radii[0]
        # redundant_tendons = True
        redundant_tendons = False
        self.n = n
        self.levels = levels
        vertices = []
        struts = []
        tendons = []
        # twist_sign = 1 will make a left-handed twist (right-handed = struts look like right-handed threads)
        twist_sign = -1
        twist = twist_sign * (math.pi / 2 - math.pi / n)
        theta_step = 2 * math.pi / n
        theta_offset = 0
        bot_z_list = [0] + [sum(heights[0:i]) * (1 - overlap) for i in range(1, len(heights))]
        for level in range(levels):
            # first we create vertices and struts, they fit neatly into a layer
            bot_z = bot_z_list[level]
            top_z = bot_z + heights[level]
            bot_vertices = [Vertex(np.array([radii[0] * math.cos(i * theta_step + theta_offset),
                                             radii[0] * math.sin(i * theta_step + theta_offset),
                                             bot_z]), level=level)
                            for i in range(n)]
            top_vertices = [Vertex(np.array([radii[1] * math.cos(i * theta_step - twist + theta_offset),
                                             radii[1] * math.sin(i * theta_step - twist + theta_offset),
                                             top_z]), level=level)
                            for i in range(n)]
            vertices.extend(bot_vertices + top_vertices)
            # always put bot vertex in strut.vertices[0]
            struts.extend([Strut([v0, v1], level=level) for v0, v1 in
                           zip(bot_vertices, rotate_list(top_vertices, twist_sign))])
            if redundant_tendons:
                tendons.extend([Tendon([v0, v1]) for v0, v1 in
                               zip(bot_vertices, rotate_list(top_vertices, -twist_sign))])
                tendons.extend([Tendon([v0, v1]) for v0, v1 in zip(bot_vertices, rotate_list(bot_vertices, 1))])
                tendons.extend([Tendon([v0, v1]) for v0, v1 in zip(top_vertices, rotate_list(top_vertices, 1))])
            theta_offset += math.pi / n - twist
        # Tendons: now let's create the tendons, they can bridge layers
        # bottom tendons
        if not redundant_tendons:
            vertex_list = [strut.vertices[0] for strut in [strut for strut in struts if strut.level == 0]]
            tendons.extend([Tendon([v0, v1]) for v0, v1 in zip(vertex_list, rotate_list(vertex_list, 1))])
            # top tendons
            vertex_list = [strut.vertices[1] for strut in [strut for strut in struts if strut.level == self.levels - 1]]
            tendons.extend([Tendon([v0, v1]) for v0, v1 in zip(vertex_list, rotate_list(vertex_list, 1))])
        for level in range(levels - 1):
            # mid tendons (waist tendons that are not bottom or top, they connect the top of each layer to bottom of
            # the layer above)
            lower_vertex_list = [strut.vertices[1] for strut in [strut for strut in struts if strut.level == level]]
            upper_vertex_list = [strut.vertices[0] for strut in [strut for strut in struts if strut.level == level + 1]]
            tendons.extend([Tendon([v0, v1]) for v0, v1 in zip(lower_vertex_list, rotate_list(upper_vertex_list, -1))])
            tendons.extend([Tendon([v0, v1]) for v0, v1 in zip(rotate_list(lower_vertex_list, -1),
                                                               rotate_list(upper_vertex_list, -1))])
            # vertical tendons top to top
            lower_vertex_list = [strut.vertices[1] for strut in [strut for strut in struts if strut.level == level]]
            upper_vertex_list = [strut.vertices[1] for strut in [strut for strut in struts if strut.level == level + 1]]
            tendons.extend([Tendon([v0, v1]) for v0, v1 in zip(lower_vertex_list, rotate_list(upper_vertex_list, 1))])
            # vertical tendons bottom to bottom
            lower_vertex_list = [strut.vertices[0] for strut in [strut for strut in struts if strut.level == level]]
            upper_vertex_list = [strut.vertices[0] for strut in [strut for strut in struts if strut.level == level + 1]]
            tendons.extend([Tendon([v0, v1]) for v0, v1 in zip(lower_vertex_list, rotate_list(upper_vertex_list, 1))])
        # vertical tendons bottom to top
        lower_vertex_list = [strut.vertices[0] for strut in [strut for strut in struts if strut.level == 0]]
        upper_vertex_list = [strut.vertices[1] for strut in [strut for strut in struts if strut.level == 0]]
        tendons.extend([Tendon([v0, v1]) for v0, v1 in zip(lower_vertex_list, rotate_list(upper_vertex_list, 1))])

        Tensegrity.__init__(self, vertices, struts, tendons)

    def strut_list(self, level=None):
        if level:
            return [strut for strut in self.struts if strut.level == level]
        else:
            return self.struts

    def stabilize(self):
        if self.levels == 1:
            pass
            """ Assumes we start with vertical struts """
            # Confirm that there are two intersections b


if __name__ == '__main__':
    tower = PrismTower(n=3, levels=1, radii=[3, 2, 2, 3, 3], heights=[8, 8, 8, 8, 8], verbose=True)
    tower.get_forces()
    tower.print_cyl()
    tower.plot()
