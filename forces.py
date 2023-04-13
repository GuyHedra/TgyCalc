""" forces.p implements a set of tools to calculate the forces acting on a tensegrity structure"""

import numpy as np
import math
import copy
import matplotlib.pyplot as plt


def other_vtx_index(member, vertex_index):
    if member[0] == vertex_index:
        return member[1]
    elif member[1] == vertex_index:
        return member[0]
    else:
        raise Exception('vertex_index not found in member')


# def prism1(n=3, hr_ratio=1, twist=math.pi / 2 - math.pi / 3, verbose=False):
#     """ Implements a single level prism and associated functions for exploring stability space. """
#     h = 1
#     radius = h / hr_ratio
#     vertices = []
#     # create vtx_coords
#     theta_step = 2 * math.pi / n
#     for h, twist in zip([0, h], [0, twist]):
#         vertices.extend([[radius * math.cos(i * theta_step + twist),
#                           radius * math.sin(i * theta_step + twist), h] for i in range(n)])
#     vertices = np.array(vertices)
#     # struts
#     # struts = [[v0, v1] for v0, v1 in zip(range(n), tp.rotate_list(list(range(n, 2*n)), -1))]
#     struts = [[v0, v1] for v0, v1 in zip(range(n), np.roll(list(range(n, 2*n)), -1))]
#     # bottom waist tendons
#     tendons = [[v0, v1] for v0, v1 in zip(range(n), np.roll(list(range(n)), 1))]
#     # top waist tendons
#     tendons.extend([[v0, v1] for v0, v1 in zip(range(n, 2*n), np.roll(list(range(n, 2*n)), 1))])
#     # vertical tendons
#     tendons.extend([[v0, v1] for v0, v1 in zip(range(n), range(n, 2*n))])
#     if verbose:
#         print('vtx_coords', vertices)
#         print('struts', struts)
#         print('tendons', tendons)
#     return vertices, struts, tendons

# todo levels=1 crashes when optimized
# def prism_tower(n=3, levels=2, height=[1, 1], radius=[[3, 3], [3, 3]], level_twist=[math.pi*(1/2-1/3), math.pi*(1/2-1/3)],
#                 interface_twist=[math.pi/3], interface_overlap=[0.2], f_strut=[1.0, 1.0], f_interlayer_v_tendon=[0.4],
#                 verbose=False):
def prism_tower(n, levels, height, radius, level_twist, interface_twist, interface_overlap, f_strut,
                f_interlayer_v_tendon, verbose=False):
    """ Returns a numpy array of vertex coordinates, a list of struts, a list of tendons, a list of pre-determined
        strut forces, a list of pre-determined tendon forces and a list of tendon_types
        len(hr_ratio) and len(level_twist) must both equal levels
        len(interface_twist) len(interface_overlap) must equal levels - 1
        len( f_interlayer_v_tendon) must equal ceil((level_count - 1) / 2) since every other layer interface needs it
    """
    # todo should prism_tower be a class that inherits Tensegrity???
    # todo add choice of chirality [left, right, alternating starting with left, alternating starting with right]
    # todo add choice of skipping struts with inter-layer vertical struts for high n values to get a stiffer structure
    if not isinstance(height, np.ndarray):
        height = np.array(height)
    if not isinstance(radius, np.ndarray):
        radius = np.array(radius)
    # Build arrays containing z, theta_start so that we can build the array of vtx_coords
    # theta_start_array and z_array have shape (levels, 2) since each level contains [bottom, top]
    theta_start_array = [[0, level_twist[0]]]
    z_array = [[0, h[0]]]
    if levels > 1:
        # todo put for loop below into two comprehensions
        for level in range(1, levels):
            bottom_theta = theta_start_array[level - 1][1] + interface_twist[level - 1]
            top_theta = bottom_theta + level_twist[level]
            theta_start_array.extend([[bottom_theta, top_theta]])
            bottom_z = z_array[level - 1][1] - interface_overlap[level - 1]
            top_z = bottom_z + height[level]
            z_array.extend([[bottom_z, top_z]])
    # Build theta array
    theta_step = 2 * math.pi / n
    theta_array = []
    # todo put for loop below in list comprehension
    for level in range(levels):
        theta_array.extend([[[theta_start_array[level][layer] + i * theta_step for i in range(n)] for layer in [0, 1]]])
    # vtx_coords.shape is (2 * level * n, 3) where 3 is [x, y, z]
    vertices = np.array([[radius[2 * level + layer] * math.cos(theta_array[level][layer][i]),
                          radius[2 * level + layer] * math.sin(theta_array[level][layer][i]),
                          z_array[level][layer]] for level in range(levels) for layer in [0, 1] for i in range(n)])
    """ Some vtx_coords index offsets:
    1:e vertex in layer=0 (bottom layer) of level: 2 * level * n
    1:e vertex in layer=1 (top layer) of level: 2 * level * n + n
    1: vertex in layer=1 (top layer) of top level: 2 * levels * n - n
    """
    # # struts shape is (levels, n, 2)
    # struts = np.array([[[v0, (v1 + 1) % n] for v0, v1 in zip(range(2 * level * n, 2 * level * n + n),
    # struts = np.array([[[v0, v1] for v0, v1 in zip(range(2 * level * n, 2 * level * n + n),
    #                                                range(2 * level * n + n, 2 * level * n + 2 * n))]
    #                    for level in range(levels)], dtype=int)
    struts = np.array([[[v0, v1] for v0, v1 in zip(range(2 * level * n, 2 * level * n + n),
                                                   np.roll(np.arange(2 * level * n + n, 2 * level * n + 2 * n,
                                                                     dtype=int), -1))]
                       for level in range(levels)], dtype=int)
    f_struts = np.array([f_strut[level] for level in range(levels) for i in range(n)])
    # todo use struts as pointer to simplify the vertex index generation of the tendon
    cap_waist_tendons = np.empty([2, n, 2], dtype=int)  # shape = (len[top, bottom], n, len[v0, v1])
    cap_waist_ten_type = np.empty([2, n], dtype=object)
    # bottom cap waist tendons
    cap_waist_tendons[0] = [[v0, v1] for v0, v1 in zip(range(n), np.roll(list(range(n)), 1))]
    cap_waist_ten_type[0].fill('cap waist bottom')
    # top cap waist tendons
    cap_waist_tendons[1] = [[v0, v1] for v0, v1 in
                            zip(range(2 * levels * n - n, 2 * levels * n),
                                np.roll(list(range(2 * levels * n - n, 2 * levels * n)), 1))]
    cap_waist_ten_type[1].fill('cap waist top')
    f_cap_waist_tendons = np.full(shape=[2 * n], fill_value=np.nan)
    # intra layer vertical tendons
    lower_vertices = np.array([v for v in range(n)], dtype=int)
    upper_vertices = np.array([v for v in range(n, 2 * n)], dtype=int)
    # intra_layer_vertical_tendons = np.array([[v0, v1] for v0, v1 in zip(lower_vertices, np.roll(upper_vertices, 1))])
    intra_layer_vertical_tendons = np.array([[v0, v1] for v0, v1 in zip(lower_vertices, upper_vertices)])
    intra_layer_vertical_ten_type = np.full(shape=[n], fill_value='intra layer vertical: layer 0', dtype=object)
    # f_intra_layer_vertical_tendons = np.full(shape=[n * (levels - 1)], fill_value=np.nan)
    f_intra_layer_vertical_tendons = np.full(shape=[n], fill_value=np.nan)
    if levels > 1:
        # interface waist tendons
        interface_waist_tendons = np.empty([levels - 1, 2 * n, 2], dtype=int)
        interface_waist_ten_type = np.empty([levels - 1, 2 * n], dtype=object)
        for interface_layer in range(levels - 1):
            # interface_layer points to the level below the interface
            # create list of lower level and list of upper level vtx_coords
            upper_vertices = np.array([v for v in range(2 * interface_layer * n + n, 2 * interface_layer * n + 2 * n)],
                                      dtype=int)
            lower_vertices = np.array([v for v in range(2 * (interface_layer + 1) * n,
                                                        2 * (interface_layer + 1) * n + n)], dtype=int)
            interleaved_vertices = np.empty(upper_vertices.size + lower_vertices.size, dtype=int)
            interleaved_vertices[0::2] = upper_vertices
            interleaved_vertices[1::2] = lower_vertices
            interface_waist_tendons[interface_layer] = [[v0, v1] for v0, v1 in
                                                        zip(interleaved_vertices, np.roll(interleaved_vertices, 1))]
            interface_waist_ten_type[interface_layer].fill('interface waist: iface layer ' + str(interface_layer))
        f_interface_waist_tendons = np.full(shape=[2 * n * (levels - 1)], fill_value=np.nan)
        # interlayer vertical tendons, top to top
        inter_layer_vertical_tendons = np.empty([levels - 1, n, 2], dtype=int)
        inter_layer_vertical_ten_type = np.empty(shape=[levels - 1, n], dtype=object)
        f_inter_layer_vertical_tendons = np.full([levels - 1, n], fill_value=np.nan)
        for interface_layer in range(levels - 1):
            lower_vertices = np.array([v for v in range(2 * interface_layer * n + n, 2 * interface_layer * n + 2 * n)],
                                      dtype=int)
            upper_vertices = np.array([v for v in
                                       range(2 * (interface_layer + 1) * n + n, 2 * (interface_layer + 1) * n + 2 * n)],
                                      dtype=int)
            inter_layer_vertical_tendons[interface_layer] = [[v0, v1] for v0, v1 in zip(lower_vertices,
                                                                                        np.roll(upper_vertices, 1))]
            inter_layer_vertical_ten_type.fill('inter layer vertical: iface layer ' + str(interface_layer))
            if interface_layer % 2 == 0:  # only even interface layers get pre-assigned forces
                f_inter_layer_vertical_tendons[interface_layer] = n * [f_interlayer_v_tendon[interface_layer]][0]
        tendons = np.concatenate((cap_waist_tendons.reshape(-1, 2), intra_layer_vertical_tendons,
                                  interface_waist_tendons.reshape(-1, 2), inter_layer_vertical_tendons.reshape(-1, 2)),
                                 axis=0)
        tendon_types = np.concatenate((cap_waist_ten_type.reshape(-1), intra_layer_vertical_ten_type,
                                       interface_waist_ten_type.reshape(-1),
                                       inter_layer_vertical_ten_type.reshape(-1)))
        f_tendons = np.concatenate((f_cap_waist_tendons, f_intra_layer_vertical_tendons,
                                    f_interface_waist_tendons,
                                    f_inter_layer_vertical_tendons.reshape(-1)))
    else:
        tendons = np.concatenate((cap_waist_tendons.reshape(-1, 2), intra_layer_vertical_tendons), axis=0)
        tendon_types = np.concatenate((cap_waist_ten_type.reshape(-1, 2), intra_layer_vertical_ten_type), axis=0)
        f_tendons = np.concatenate((f_cap_waist_tendons.reshape(-1, 2), f_intra_layer_vertical_tendons), axis=0)

    if verbose:
        # print('vtx_coords', vtx_coords)
        print('struts', struts)
        print('tendons', tendons)
    # todo add name as a return value (or better yet, make prism_tower a class that inherits Tensegrity)
    return n, vertices, struts.reshape(-1, 2), tendons.reshape(-1, 2), f_struts, f_tendons, tendon_types


class Tensegrity:
    # todo  remove n argument and also remove from return values from prism_tower
    def __init__(self, n, vertices, struts, tendons, f_struts, f_tendons, tendon_types, name):
        self.n = n
        self.vtx_coords = vertices
        self.struts = struts
        self.tendons = tendons
        self.f_struts = f_struts
        self.f_tendons = f_tendons
        self.tendon_types = tendon_types
        self.name = name
        self.tendon_force = 1
        self.vtx_f_unit_vectors = {}  # key is vertex index: value is a list unit force vectors [strut, tendon0, tendon1,..]
        """ Build dictionaries that let us find the struts and tendons that are connected to a particular vertex"""
        self.vtx_struts = {}  # key is vertex index: value is strut index
        self.vtx_tendons = {}  # key is vertex index: values is list of tendon indices
        for i, strut in enumerate(self.struts):
            if strut[0] in self.vtx_struts or strut[1] in self.vtx_struts:
                raise Exception('vertex has more than one associated strut')
            else:
                self.vtx_struts[strut[0]] = i
                self.vtx_struts[strut[1]] = i
        for i, tendon in enumerate(self.tendons):
            if tendon[0] in self.vtx_tendons:
                self.vtx_tendons[tendon[0]].extend([i])
            else:
                self.vtx_tendons[tendon[0]] = [i]
            if tendon[1] in self.vtx_tendons:
                self.vtx_tendons[tendon[1]].extend([i])
            else:
                self.vtx_tendons[tendon[1]] = [i]
        """ We will use the convention that strut forces are directed towards the vertex and tendon forces are directed 
        away from the vertex:
        For each vertex the vector sum of all forces is equal to 0 for a balanced tensegrity """
        """ Find the strut force vectors """
        for strut in self.struts:
            for i, vertex_index in enumerate(strut):
                if vertex_index in self.vtx_f_unit_vectors:
                    raise Exception('multiple struts are associated with vertex', vertex_index)
                # todo use other_vtx_index function in line below and for tendon code below
                vector = self.vtx_coords[vertex_index] - self.vtx_coords[strut[(i + 1) % 2]]
                self.vtx_f_unit_vectors[vertex_index] = [vector / np.linalg.norm(vector)]
        """ Find the tendon force vectors """
        for tendon in self.tendons:
            for i, vertex_index in enumerate(tendon):
                vector = self.vtx_coords[tendon[(i + 1) % 2]] - self.vtx_coords[vertex_index]
                if np.linalg.norm(vector) <= 1e-15:
                    print('*+*+ Warning tendon', tendon, 'has zero magnitude *+*+')
                    self.vtx_f_unit_vectors[vertex_index].extend([vector])
                else:
                    self.vtx_f_unit_vectors[vertex_index].extend([vector / np.linalg.norm(vector)])

    def cost(self, verbose=False):
        """ Pre Determined Forces: Each vertex has one or more pre-determined forces such that a square full rank
            matrix can represent all the forces in 3-space. Struts are always assigned a pre-determined force,
            and when needed vertical interlayer tendons are also assigned a pre-determined force in order to produce
            matrices that can be solved.
            Virtual Tendons: Vertices that have only two tendons are assigned a virtual tendon. The virtual tendon
            is the normalized cross product of the two real tendons. The square of the force of the virtual tendon
            is added to the cost so that the optimization process drives the virtual tendon force to zero.
            Force Differences: Once the unknown member force magnitudes have been found for all the vertices,
            each real tendon will have two force magnitudes, one from each of the members
            vertices. The force difference is the difference between these two forces.
            Returns: cost: sum of squares of all real tendon force differences plus the sum of the squares of all the
            virtual tendons.
            f_diff_squared: a list of the square of the force magnitude differences for each tendon
            f_mean: a list of the mean of the two force magnitudes for each tendon
            f_virtual: a list of the force magnitudes for each virtual tendon
            vtf_success: a boolean flag, True if all the force matrices described above were successfully solved
        """
        # Solve for the unknown forces at each vertex
        vtx_tendon_forces = []
        f_virtual = []
        for vtx_index in range(len(self.vtx_coords)):
            forces, vtf_success = self.vtx_tendon_forces(vtx_index)
            if vtf_success:
                vtx_tendon_forces.extend([forces])
                if len(self.vtx_tendons[vtx_index]) == 2:
                    f_virtual.extend([forces[2]])
            else:
                break
        # Find the force differences and means
        if vtf_success:
            f_diff_squared = []
            f_mean = []
            # once for each associated tendon with this implementation
            for tendon_index, tendon in enumerate(self.tendons):
                # prepare the index that selects the correct tendon from the vertex's list of tendons
                vtx0_to_tendon_index = self.vtx_tendons[tendon[0]].index(tendon_index)
                # get the vertex 0 force for this tendon
                vtx0_force = vtx_tendon_forces[tendon[0]][vtx0_to_tendon_index]
                # prepare the index that selects the correct tendon from the vertex's list of tendons
                vtx1_to_tendon_index = self.vtx_tendons[tendon[1]].index(tendon_index)
                # get the vertex 1 force for this tendon
                vtx1_force = vtx_tendon_forces[tendon[1]][vtx1_to_tendon_index]
                # todo add a parameter 'cost_only' so we don't compute f_mean or return f_diff_square and f_mean while
                # todo in an optimization loop unless desired
                # f_diff_squared.extend(abs(vtx0_force - vtx1_force))
                f_diff_squared.extend((vtx0_force - vtx1_force) ** 2)
                f_mean.extend((vtx0_force + vtx1_force) / 2)
            cost = np.sum(np.array(f_diff_squared)) + np.sum(np.array(f_virtual) ** 2)
            return cost, f_diff_squared, f_mean, f_virtual, vtf_success
        else:
            return -1, -1, -1, -1, False

    def vtx_tendon_forces(self, vtx_index, verbose=False):
        """ Return the tendon forces on the vertex at vtx_index """
        a = np.column_stack([np.transpose(self.vtx_f_unit_vectors[vtx_index])])
        member_count = len(self.vtx_f_unit_vectors[vtx_index])
        if self.name == 'kite':
            a[2] = [1, 0, 0]
            b = np.array([(member_count - 1) * [0], [1]])
        elif self.name == 'prism':
            if member_count == 3:  # 1 strut, 2 tendons
                """ using the normal approach yields a singular matrix for two-tendon vertices. We will
                add a virtual tendon perpendicular to the plane of the two real tendons in order to create a
                useful matrix for the solver. We will add the magnitude of the force on said virtual tendon to
                the cost (error) total so that the optimizer will drive it to zero. When the virtual tendon force
                is zero, then the two real tendons lie in a plane with the strut as desired """
                virtual_tendon_vector = np.cross(self.vtx_f_unit_vectors[vtx_index][1],
                                                 self.vtx_f_unit_vectors[vtx_index][2])
                virtual_tendon_unit_vec = virtual_tendon_vector / np.linalg.norm(virtual_tendon_vector)
                a = np.column_stack((a, virtual_tendon_unit_vec))
                a = np.vstack((a, [1] + member_count * [0]))
                b = np.array(member_count * [[0.0]] + [[self.f_struts[self.vtx_struts[vtx_index]]]])
            elif member_count == 4:  # 1 strut, 3 tendons
                a = np.vstack((a, [1] + (member_count - 1) * [0]))
                b = np.array((member_count - 1) * [[0.0]] + [[self.f_struts[self.vtx_struts[vtx_index]]]])
            elif member_count == 5:  # 1 strut, 4 tendons
                a = np.vstack((a, [1] + (member_count - 1) * [0]))
                a = np.vstack((a, (member_count - 1) * [0.0] + [1]))
                interlayer_tendon_index = [t_index for t_index in self.vtx_tendons[vtx_index]
                                           if not np.isnan(self.f_tendons[t_index])]
                if len(interlayer_tendon_index) > 1:
                    raise Exception('Multiple tendons connected to the same vertex with pre-assigned forces '
                                    'not supported')
                b = np.array((member_count - 2) * [[0.0]] + [[self.f_struts[self.vtx_struts[vtx_index]]]] +
                             [[self.f_tendons[interlayer_tendon_index[0]]]])
            """ return None if matrix is singular """
            rank = np.linalg.matrix_rank(a)
            if rank == a.shape[0]:
                """ Solve the force matrix"""
                x = np.linalg.solve(a, b)
                if verbose:
                    f_vectors = (a.T * x)[0:member_count, 0:3]
                    sum_f_vectors = np.sum(f_vectors, axis=0)
                    print('vtx_index', vtx_index)
                    print('a', a)
                    print('b', b)
                    print('x', x)
                    print('f_vectors', f_vectors)
                    print('sum of force vectors', sum_f_vectors)
                    # print('total_tendon_forces', total_tendon_forces)
                return x[1:], True  # skip the strut forces at x[0]
            else:
                return None, False

    def plot(self, struts=True, tendons=True, lateral_f=False, vertex_f=False, axes=False, debug=False):
        """ plots tendons and struts using matplotlib """

        ax = plt.axes(projection='3d')
        x_span = 2 * max(vertex[0] for vertex in self.vtx_coords) - 2 * min(vertex[0] for vertex in self.vtx_coords)
        y_span = 2 * max(vertex[1] for vertex in self.vtx_coords) - 2 * min(vertex[1] for vertex in self.vtx_coords)
        z_span = 2 * max(vertex[2] for vertex in self.vtx_coords) - 2 * min(vertex[2] for vertex in self.vtx_coords)
        span = max(x_span, y_span, z_span)
        ax.set_xlim(-span / 2, span / 2)
        ax.set_ylim(-span / 2, span / 2)
        ax.set_zlim(0, span)
        if not axes:
            ax.set_axis_off()
        if struts:
            for strut in self.struts:
                # x = [vertex[0] for vertex in strut.vtx_coords]
                x = [self.vtx_coords[vertex][0] for vertex in strut]
                y = [self.vtx_coords[vertex][1] for vertex in strut]
                z = [self.vtx_coords[vertex][2] for vertex in strut]
                # y = [vertex[1] for vertex in strut.vtx_coords]
                # z = [vertex[2] for vertex in strut.vtx_coords]
                ax.plot3D(x, y, z, 'red', linewidth=3)
        if tendons:
            for tendon in self.tendons:
                # x = [vertex[0] for vertex in tendon.vtx_coords]
                # y = [vertex[1] for vertex in tendon.vtx_coords]
                # z = [vertex[2] for vertex in tendon.vtx_coords]
                x = [self.vtx_coords[vertex][0] for vertex in tendon]
                y = [self.vtx_coords[vertex][1] for vertex in tendon]
                z = [self.vtx_coords[vertex][2] for vertex in tendon]
                ax.plot3D(x, y, z, 'grey', linewidth=1)
        if debug:
            # for vertex in self.vtx_coords:
            strut = self.struts[0]
            vertex = strut.vtx_coords[0]
            tendon = vertex.tendon_list[0]
            strut_vector = np.array(strut.position_vec(strut.other_vtx_index(vertex)))
            # ax.quiver(*vertex, *strut_vector, color='orange', linewidth=3)
            coords = [val for pair in zip(vertex, strut_vector) for val in pair]
            x = [vertex[0] for vertex in [vertex, strut_vector]]
            y = [vertex[1] for vertex in [vertex, strut_vector]]
            z = [vertex[2] for vertex in [vertex, strut_vector]]
            ax.plot3D(x, y, z, 'orange', linewidth=5)
            # *self.strut_list[0].position_vec(self.strut_list[0].vtx_coords[1]), 'green', linewidth=3)
            tendon_vector = np.array(tendon.position_vec(tendon.other_vtx_index(vertex)))
            # ax.quiver(*vertex, *tendon_vector, color='orange', linewidth=3)
            tendon_strut_component = ((np.dot(tendon_vector, strut_vector) / np.linalg.norm(strut_vector) ** 2)
                                      * strut_vector)
        plt.show()

    @staticmethod
    def print_build(self):
        """ Print the information needed to build the tensegrity """
        # todo add strut to ground angle
        # todo add tendon lengths (vtx to vtx, hub compensated and hook compensated)
        vv_strut_lengths = [np.linalg.norm(self.vtx_coords[v0] - self.vtx_coords[v1])
                            for v0, v1 in self.struts[::self.n]]
        for length in vv_strut_lengths:
            print('strut length ', length)

    @staticmethod
    def print_tendon_forces(self, f_tendon):
        print('*** tendon forces by type ***')
        for ten_type, ten_force in zip(self.tendon_types[::self.n], f_tendon[::self.n]):
            print(ten_type, 'force: ', str(ten_force))

    @staticmethod
    def plot_gradient_descent(self, param_history, cost_history, f_tendon_diff_history):

        def legend_without_duplicate_labels(figure):
            handles, labels = plt.gca().get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            figure.legend(by_label.values(), by_label.keys(), loc='lower left')

        plt.style.use('seaborn-whitegrid')
        markers = ['.', ',', 'o', 'v', '^', '<', '>', '1', '2', '3', '4', '8', 's', 'p', 'P',
                   '*', '.', ',', 'o', 'v', '^', '<', '>', '1', '2', '3', '4', '8', 's', 'p', 'P', '*',
                   '*', '.', ',', 'o', 'v', '^', '<', '>', '1', '2', '3', '4', '8', 's', 'p', 'P', '*']
        colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black',
                  'blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black',
                  'blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black',
                  'blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black']
        # Tendon Force Differences
        diff_fig = plt.figure(1)
        diff_ax = plt.axes()
        # param_ax = param_fig.add_axes()
        levels = param_history[0].levels
        diff_fig.suptitle("Tendon Force Difference")
        diff_ax.set_xlabel('Gradient Descent Step')
        plot_tendons = self.tendons[::self.n]  # For radially symmetric prism towers only plot every nth tendon
        tendon_labels = [str(tendon) for tendon in plot_tendons]
        f_diff_array = np.array([f_tendon_diff_of_t[::self.n] for f_tendon_diff_of_t in f_tendon_diff_history])
        for tendon_index, tendon_label in enumerate(tendon_labels):
            diff_ax.plot(np.arange(len(f_diff_array[:, tendon_index])), f_diff_array[:, tendon_index],
                         color=colors[tendon_index], marker=markers[tendon_index], label=tendon_label)
        legend_without_duplicate_labels(diff_fig)
        # Tune Params
        param_fig = plt.figure(2)
        param_ax = plt.axes()
        # param_ax = param_fig.add_axes()
        levels = param_history[0].levels
        param_fig.suptitle("Parameters")
        param_ax.set_xlabel('Gradient Descent Step')
        p_history_of_t = np.array([params.tune_param_array for params in param_history])
        for p_index, p_label in enumerate(param_history[0].tune_param_labels):
            param_ax.plot(np.arange(len(p_history_of_t[:, p_index])), p_history_of_t[:, p_index], color=colors[p_index],
                          marker=markers[p_index], label=p_label)
        legend_without_duplicate_labels(param_fig)
        # Cost
        cost_fig = plt.figure(3)
        cost_ax = plt.axes()
        cost_fig.suptitle("Cost History")
        cost_ax.set_xlabel('Gradient Descent Step')
        cost_ax.plot(np.arange(len(cost_history)), cost_history, color=colors[p_index],
                     marker=markers[p_index], label=p_label)
        legend_without_duplicate_labels(cost_fig)
        plt.show()


class TowerTuneParams:
    def __init__(self, n, levels, height, radius, level_twist, interface_twist, interface_overlap, f_strut,
                 f_interlayer_v_tendon, tune_param_list):
        """ Stores prism tower parameters and provides several means of setting and accessing the parameters,
        including support for stability optimization """
        # todo add check for correct argument lengths and types
        self.n = n
        self.levels = levels
        # self.strut_length = strut_length  todo add support for strut_length
        self.height = height
        self.radius = radius
        self.level_twist = level_twist
        self.interface_twist = interface_twist
        self.interface_overlap = interface_overlap
        self.f_strut = f_strut
        self.f_interlayer_v_tendon = f_interlayer_v_tendon
        self.p_names = ['n', 'levels', 'height', 'radius', 'overlap radius', 'level twist', 'interface twist',
                        'interface overlap', 'strut force', 'interlayer tendon force']
        # exclude n, levels, height and radius with self.p_names[4:]
        if False not in [p_name in self.p_names[4:] for p_name in tune_param_list]:
            self.tune_param_list = tune_param_list
        else:
            raise Exception('tune_param_list elements must be contained in self.p_names[3:]')

    def set_tune_params(self, p_list):
        """ set tune parameters to values in vector p_list """
        p_list_ptr = 0
        for p_name in self.tune_param_list:
            if p_name == 'overlap radius':
                for i in range(len(self.radius[2:][::2])):  # start with 2nd level and then take every other layer
                    self.radius[2 + 2 * i] = p_list[p_list_ptr]
                    p_list_ptr += 1
            # todo replace for loops below with list comprehensions
            elif p_name == 'level twist':
                for i in range(len(self.level_twist)):
                    self.level_twist[i] = p_list[p_list_ptr]
                    p_list_ptr += 1
            elif p_name == 'interface twist':
                for i in range(len(self.interface_twist)):
                    self.interface_twist[i] = p_list[p_list_ptr]
                    p_list_ptr += 1
            elif p_name == 'interface overlap':
                for i in range(len(self.interface_overlap)):
                    self.interface_overlap[i] = p_list[p_list_ptr]
                    p_list_ptr += 1
            elif p_name == 'strut force':
                for i in range(len(self.f_strut)):
                    self.f_strut[i] = p_list[p_list_ptr]
                    p_list_ptr += 1
            elif p_name == 'interlayer tendon force':
                for i in range(len(self.f_interlayer_v_tendon)):
                    self.f_interlayer_v_tendon[i] = p_list[p_list_ptr]
                    p_list_ptr += 1

    def set_single_tune_param(self, i, p):
        p_list = self.tune_param_array
        p_list[i] = p
        self.set_tune_params(p_list)

    @property
    def build_args(self):
        return self.n, self.levels, self.height, self.radius, self.level_twist, self.interface_twist, \
            self.interface_overlap, self.f_strut, self.f_interlayer_v_tendon

    @property
    def tune_param_array(self):
        """ return a vector containing all of the tower_prism tune parameters """
        p_val_list = []
        for p_name in self.tune_param_list:
            if p_name == 'overlap radius':
                p_val_list.extend(self.radius[2:][::2])
            elif p_name == 'level twist':
                p_val_list.extend(self.level_twist)
            elif p_name == 'interface twist':
                p_val_list.extend(self.interface_twist)
            elif p_name == 'interface overlap':
                p_val_list.extend(self.interface_overlap)
            elif p_name == 'strut force':
                p_val_list.extend(self.f_strut)
            elif p_name == 'interlayer tendon force':
                p_val_list.extend(self.f_interlayer_v_tendon)
        return np.array(p_val_list)

    @property
    def tune_param_labels(self):
        """ return a list containing labels for each of the tower_prism tune parameters """
        p_label_list = []
        for p_name in self.tune_param_list:
            if p_name == 'overlap radius':
                p_label_list.extend([p_name + ' ' + str(i) for i in range(self.levels - 1)])
            elif p_name == 'level twist':
                p_label_list.extend([p_name + ' ' + str(i) for i in range(self.levels)])
            elif p_name == 'interface twist':
                p_label_list.extend([p_name + ' ' + str(i) for i in range(self.levels - 1)])
            elif p_name == 'interface overlap':
                p_label_list.extend([p_name + ' ' + str(i) for i in range(self.levels - 1)])
            elif p_name == 'strut force':
                p_label_list.extend([p_name + ' ' + str(i) for i in range(self.levels)])
            elif p_name == 'interlayer tendon force':
                p_label_list.extend([p_name + ' ' + str(i) for i in range(math.ceil((self.levels - 1) / 2))])
        return np.array(p_label_list)

    @property
    def print_tune(self):
        # for p_index in
        print('overlap radius', [p for p in self.radius[2:][::2]])
        print('level twist', [math.degrees(tw) for tw in self.level_twist])
        print('interface twist', [math.degrees(itw) for itw in self.interface_twist])
        print('interface overlap', [io for io in self.interface_overlap])
        print('strut force', [fs for fs in self.f_strut])
        print('interlayer vertical tendon force', [ft for ft in self.f_interlayer_v_tendon])


def stabilize_tower_grad_descent(tower_params, learn_rate=0.001, learn_rate_damping=0.95, max_steps=100,
                                 max_cost=0.001, min_difference=1e-5, epsilon=0.0001, verbose=False, debug=False):
    """ Adjust the initial prism_tower parameters to achieve a stable structure using home-grown gradient descent"""
    # todo should we try tuning layer overlap, radius offset and all the twists first, and then forces?
    # todo provide message if solve fails for a vertex (give vertex and step?)
    tune_p_count = len(tower_params.tune_param_array)
    cost_history = []
    f_tendon_diff_history = []
    param_history = []
    param_history.extend([tower_params])
    cost, f_tendon_diff, f_tendon, f_virtual, success = Tensegrity(*prism_tower(*param_history[0].build_args),
                                                                   name='prism').cost()
    cost_history.extend([cost])
    f_tendon_diff_history.extend([f_tendon_diff])
    termination_code = 'none'
    step = 0
    error_difference = 2 * min_difference
    if success:
        while step < max_steps and abs(error_difference) > min_difference and abs(cost_history[-1]) > max_cost:
            step += 1
            if debug:
                print('step ', step)
            # todo consider extracting the gradient from the previous two steps
            gradient = tower_gradient(param_history[step - 1], cost_history[step - 1], epsilon=epsilon)
            new_params = copy.deepcopy(param_history[step - 1])
            if len(cost_history) > 2 and cost_history[-1] > cost_history[-2]:
                learn_rate = learn_rate_damping * learn_rate

            new_params.set_tune_params(param_history[step - 1].tune_param_array -
                                       gradient * learn_rate * param_history[step - 1].tune_param_array)
            param_history.extend([new_params])
            cost, f_tendon_diff, f_tendon, f_virtual, success = \
                Tensegrity(*prism_tower(*param_history[-1].build_args), name='prism').cost()
            if success:
                cost_history.extend([cost])
                f_tendon_diff_history.extend([f_tendon_diff])
                error_difference = cost_history[step - 1] - cost_history[step]
            else:
                break
        if step >= max_steps:
            termination_code = 'max steps reached'
        if abs(error_difference) <= min_difference:
            termination_code = 'error difference <= min_difference'
        if abs(cost_history[-1]) <= max_cost:
            termination_code = 'cost < max cost'
        if not success:
            termination_code = 'force matrix solution failed'
        if verbose:
            # print('grad descent finished on step', step, 'with cost history', cost_history)
            print('*** initial parameters ***')
            param_history[0].print_tune
            print('*** initial cost', cost_history[0])
            print('*** grad descent finished on step', step, 'with termination code:', termination_code)
            print('*** final cost', cost_history[-1])
            print('*** final params')
            param_history[-1].print_tune

    return param_history, cost_history, f_tendon_diff_history, step, termination_code


def tower_gradient(tune_params, cost_of_tune_params, epsilon):
    cost_of_x_plus_epsilon = np.zeros_like(tune_params.tune_param_array)
    for i, p in enumerate(tune_params.tune_param_array):
        # todo check to see if copy.copy is needed
        t_params_w_epsilon = copy.deepcopy(tune_params)
        t_params_w_epsilon.set_single_tune_param(i, p + epsilon)
        # np.sum(np.array(Tensegrity(*prism_tower(*t_params_w_epsilon.arguments).tendon_f_difference())))
        cost_of_x_plus_epsilon[i] = \
            np.sum(np.array(Tensegrity(*prism_tower(*t_params_w_epsilon.build_args), name='prism').cost()[0]))
        gradient = (cost_of_x_plus_epsilon - cost_of_tune_params) / epsilon
    return gradient


if __name__ == '__main__':
    mode = 'stabilize prism tower'
    # mode = 'hat rack 1.0'
    if mode == 'hat rack 1.0':
        """ hourglass shaped hat rack """
        strut_count = 4
        level_count = 3
        # h_to_r = level_count * [3]
        dowel_length = 915  # in mm 1/2" wood dowel 36"
        bottom_extension = 30  # mm measured from tip of dowel to inside of hub (includes hub length)
        bot_layer_top_extension = 250
        mid_layer_top_extension = 250
        top_layer_top_extension = 200
        vv_bot_layer = dowel_length - bottom_extension - bot_layer_top_extension  # vv denotes vertex to vertex len
        vv_mid_layer = dowel_length - bottom_extension - mid_layer_top_extension
        vv_top_layer = dowel_length - bottom_extension - top_layer_top_extension
        print('vv_bot_layer:', vv_bot_layer, 'vv_mid_layer:', vv_mid_layer, 'vv_top_layer:', vv_top_layer)
        reference_radius = (vv_bot_layer / 2) * (635 / 887)  # fudge factor from trial run
        r_off = 0.3 * reference_radius
        radii = [1.1 * reference_radius + r_off, reference_radius - r_off] + \
                [reference_radius - r_off, reference_radius - r_off] + \
                [reference_radius - r_off, reference_radius + r_off]
        h = [2 * reference_radius * (635 / 627),
             (2 * reference_radius) * (635 / 545.3) * (635 / 610),
             (2 * reference_radius) * (685 / 630) * (685 / 673)]

        l_twist = level_count * [math.pi * (1 / 2 - 1 / strut_count)]
        iface_twist = (level_count - 1) * [math.pi / strut_count]
        # iface_overlap = (level_count - 1) * [0.25]
        iface_overlap = (level_count - 1) * [0.35]
        # frc_strut = (level_count - 1) * [1.0, 1.1]
        frc_strut = level_count * [1.0]
        f_interlayer_tendon = math.ceil((level_count - 1) / 2) * [0.2]  # interfaces 0, 2, 4... get f_interlayer_tendon
        learn_rate = 0.006
        # learn_rate_damping = 0.98
        learn_rate_damping = 0.97
        t_params = TowerTuneParams(n=strut_count, levels=level_count, height=h, radius=radii, level_twist=l_twist,
                                   interface_twist=iface_twist, interface_overlap=iface_overlap,
                                   f_strut=frc_strut, f_interlayer_v_tendon=f_interlayer_tendon,
                                   tune_param_list=['overlap radius', 'level twist', 'interface twist',
                                                    # 'interface overlap', 'strut force',
                                                    'interface overlap',
                                                    'interlayer tendon force'])
        param_hist, cost_hist, f_tendon_diff_hist, end_step, termination_msg = \
            stabilize_tower_grad_descent(t_params, learn_rate=learn_rate, learn_rate_damping=learn_rate_damping,
                                         max_steps=200, max_cost=0.001, min_difference=1e-7, epsilon=1e-7, verbose=True)
        thing = Tensegrity(*prism_tower(*param_hist[-1].build_args, verbose=False),
                           name='prism')
        final_cost, final_f_tendon_diff, final_f_tendon, final_f_virtual, final_success = thing.cost()
        thing.print_tendon_forces(thing, final_f_tendon)
        print('>>> f_virtual ', final_f_virtual)
        thing.print_build(thing)
        # thing.plot_gradient_descent(thing, param_hist, cost_hist, f_tendon_diff_hist)
        # thing.plot()
    elif mode == 'stabilize prism tower':
        strut_count = 4
        level_count = 3
        # h_to_r = level_count * [3]
        h = level_count * [2]
        radii = level_count * [1, 1]
        if level_count > 1:
            # radii = [[1, 1], (level_count - 1) * [0.78, 1]]
            # radii = [1, 1] + (level_count - 1) * [0.8, 1]
            radii = [1, 1] + (level_count - 1) * [1, 1]
        l_twist = level_count * [math.pi * (1 / 2 - 1 / strut_count)]
        iface_twist = (level_count - 1) * [math.pi / strut_count]
        # iface_overlap = (level_count - 1) * [0.25]
        iface_overlap = (level_count - 1) * [0.35]
        # frc_strut = (level_count - 1) * [1.0, 1.1]
        frc_strut = level_count * [1.0]
        f_interlayer_tendon = math.ceil((level_count - 1) / 2) * [0.2]  # interfaces 0, 2, 4... get f_interlayer_tendon
        if strut_count == 3:
            learn_rate = 0.01
            learn_rate_damping = 0.9
        elif level_count > 2:
            learn_rate = 0.007
            # learn_rate_damping = 0.98
            learn_rate_damping = 0.9
        else:
            learn_rate = 0.005
            learn_rate_damping = 0.99
        # todo report initial and final learn rates (plot them?)
        t_params = TowerTuneParams(n=strut_count, levels=level_count, height=h, radius=radii, level_twist=l_twist,
                                   interface_twist=iface_twist, interface_overlap=iface_overlap,
                                   f_strut=frc_strut, f_interlayer_v_tendon=f_interlayer_tendon,
                                   tune_param_list=['overlap radius', 'level twist', 'interface twist',
                                                    # 'interface overlap', 'strut force',
                                                    'interface overlap',
                                                    'interlayer tendon force'])
        param_hist, cost_hist, f_tendon_diff_hist, end_step, termination_msg = \
            stabilize_tower_grad_descent(t_params, learn_rate=learn_rate, learn_rate_damping=learn_rate_damping,
                                         max_steps=200, max_cost=0.001, min_difference=1e-7, epsilon=1e-7, verbose=True)
        thing = Tensegrity(*prism_tower(*param_hist[-1].build_args, verbose=False),
                           name='prism')
        final_cost, final_f_tendon_diff, final_f_tendon, final_f_virtual, final_success = thing.cost()
        thing.print_tendon_forces(thing, final_f_tendon)
        print('>>> f_virtual ', final_f_virtual)
        thing.plot_gradient_descent(thing, param_hist, cost_hist, f_tendon_diff_hist)
        thing.plot()
