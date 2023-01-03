""" TgyPar creates stable prism towers from a set of input parameters"""
import math
import numpy as np
import matplotlib.pyplot as plt


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
    def __init__(self, coords):
        self.coords = coords

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
    def __init__(self, vertices):
        self.vertices = vertices


class Tendon:
    def __init__(self, vertices):
        self.vertices = vertices


class PrismTower:
    def __init__(self, n=3, levels=2, radii=[3, 3], heights=[10, 10], verbose=0):
        """ return 3 numpy.arrays that describe a stable tensegrity structure:
        vertex coordinates, tendon vertices, strut vertices """
        # Assume that the bottom vertices lie on the x,y plane
        # bottom vertices are fully described by the input parameters n and radii[0]
        twist = 1  # twist = -1 will make a left (or is it right?) handed twist
        theta_step = 2 * math.pi / n
        self.bot_vertices = [Vertex(np.array([radii[0] * math.cos(i * theta_step),
                                              radii[0] * math.sin(i * theta_step),
                                              0]))
                             for i in range(n)]
        self.top_vertices = [Vertex(np.array([radii[1] * math.cos(i * theta_step),
                                              radii[1] * math.sin(i * theta_step),
                                              heights[0]]))
                             for i in range(n)]
        self.struts = [Strut([v0, v1]) for v0, v1 in zip(self.bot_vertices, self.top_vertices)]
        self.bot_tendons = [Tendon([v0, v1]) for v0, v1 in zip(self.bot_vertices, rotate_list(self.bot_vertices, 1))]
        self.top_tendons = [Tendon([v0, v1]) for v0, v1 in zip(self.top_vertices, rotate_list(self.top_vertices, 1))]
        self.vert_tendons = [Tendon([v0, v1]) for v0, v1 in zip(self.bot_vertices, rotate_list(self.top_vertices, 1))]
        self.vertices = self.bot_vertices + self.top_vertices
        self.tendons = self.bot_tendons + self.top_tendons + self.vert_tendons
        # for i in range(n):
        #     Vertex(np.array(math.cos(i * theta_step), math.sin(i * theta_step), 0))

    def plot(self, struts=True, tendons=True, lateral_f=False, vertex_f=False):
        """ plots tendons and struts using matplotlib """
        ax = plt.axes(projection='3d')
        ax.set_xlim(-5, 5)
        ax.set_ylim(-5, 5)
        ax.set_zlim(0, 10)
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
                  f'{"Twist": <{number_width}}',
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
                      f'{str(round(vertex0_cyl[1] - vertex1_cyl[1], precision)): <{number_width}}',
                      # f'{str(np.array(strut.spring_f_vec)): <{array_3_width}}'
                      )



if __name__ == '__main__':
    tower = PrismTower(n=3, levels=2, radii=[3, 3], heights=[10, 10], verbose=True)
    tower.print_cyl()
    tower.plot()
