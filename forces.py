""" forces.p implements a set of tools to calculate the forces acting on a tensegrity structure"""

import numpy as np
import math
import TgyPar as tp
import matplotlib.pyplot as plt


def other_vertex_index(member, vertex_index):
    if member[0] == vertex_index:
        return member[1]
    elif member[1] == vertex_index:
        return member[0]
    else:
        raise Exception('vertex_index not found in member')


def kite():
    vertex_list = np.array([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]])
    strut_list = [[0, 2], [1, 3]]
    tendon_list = [[0, 1], [1, 2], [2, 3], [3, 0]]
    return vertex_list, strut_list, tendon_list


def prism_n3(alpha=5 * math.pi / 6):
    vertices = np.array([[math.cos(0), math.sin(0), 0],
                      [math.cos(2*math.pi/3), math.sin(2*math.pi/3), 0],
                      [math.cos(4 * math.pi / 3), math.sin(4 * math.pi / 3), 0],
                      [math.cos(alpha), math.sin(alpha), 1],
                      [math.cos(alpha + 2 * math.pi / 3), math.sin(alpha + 2 * math.pi / 3), 1],
                      [math.cos(alpha + 4 * math.pi / 3), math.sin(alpha + 4 * math.pi / 3), 1]])
    #these numbers are indicies in the vertex_list array for the vertex_list on either end of the member
    compression_members = np.array([[0, 3],
                                    [1, 4],
                                    [2, 5]])
    tension_members = np.array([[0, 1],
                                [1, 2],
                                [2, 0],
                                [3, 4],
                                [4, 5],
                                [5, 3],
                                [1, 3],
                                [2, 4],
                                [0, 5]])
    return vertices, compression_members, tension_members


def prism1(n=3, hr_ratio=1, twist=math.pi / 2 - math.pi / 3, verbose=False):
    """ Implements a single level prism and associated functions for exploring stability space. """
    h = 1
    radius = h / hr_ratio
    vertices = []
    # create vertices
    theta_step = 2 * math.pi / n
    for h, twist in zip([0, h], [0, twist]):
        vertices.extend([[radius * math.cos(i * theta_step + twist),
                          radius * math.sin(i * theta_step + twist), h] for i in range(n)])
    vertices = np.array(vertices)
    # struts
    # struts = [[v0, v1] for v0, v1 in zip(range(n), tp.rotate_list(list(range(n, 2*n)), -1))]
    struts = [[v0, v1] for v0, v1 in zip(range(n), np.roll(list(range(n, 2*n)), -1))]
    # bottom waist tendons
    tendons = [[v0, v1] for v0, v1 in zip(range(n), np.roll(list(range(n)), 1))]
    # top waist tendons
    tendons.extend([[v0, v1] for v0, v1 in zip(range(n, 2*n), np.roll(list(range(n, 2*n)), 1))])
    # vertical tendons
    tendons.extend([[v0, v1] for v0, v1 in zip(range(n), range(n, 2*n))])
    if verbose:
        print('vertices', vertices)
        print('struts', struts)
        print('tendons', tendons)
    return vertices, struts, tendons


def prism_tower(n=3, levels=2, hr_ratio=[2, 2], level_twist=[math.pi*(1/2-1/3), math.pi*(1/2-1/3)],
                interface_twist=[math.pi/3], interface_overlap=[0.1], verbose=False):
    """ Returns an numpy array of vertex coordinates, a list of struts and a list of tendons
        len(hr_ratio) and len(level_twist) must both equal levels
        len(interface_twist) must equal levels - 1
    """
    # todo add strut_length and radius as alternative to hr_ratio
    if not isinstance(hr_ratio, np.ndarray):
        hr_ratio = np.array(hr_ratio)
    # h and radius have shape (levels,) since each element belongs to a level
    h = np.array(levels * [1])
    radius = h / hr_ratio
    # Build arrays containing z, theta_start so that we can build the array of vertices
    # theta_start_array and z_array have shape (levels, 2) since each level contains [bottom, top]
    theta_start_array = [[0, level_twist[0]]]
    z_array = [[0, h[0]]]
    # z_array = [[0, h[0]]
    if levels > 1:
        # todo put for loop below into two comprehensions
        for level in range(1, levels):
            # bottom_theta = theta_start_array[level - 1][1] - interface_twist[level - 1]
            bottom_theta = theta_start_array[level - 1][1] + interface_twist[level - 1]
            top_theta = bottom_theta + level_twist[level]
            theta_start_array.extend([[bottom_theta, top_theta]])
            bottom_z = z_array[level - 1][1] - interface_overlap[level - 1]
            top_z = bottom_z + h[level]
            z_array.extend([[bottom_z, top_z]])
    # Build theta array
    theta_step = 2 * math.pi / n
    theta_array = []
    for level in range(levels):
        theta_array.extend([[[theta_start_array[level][top] + i * theta_step for i in range(n)] for top in [0, 1]]])
    # vertices = np.empty([levels, 2, n, 3])
    # for level in range(levels):
        # for layer in [0, 1]:
    # vertices.shape is (level, layer, i, xyz)
    vertices = np.array([[[[radius[level]*math.cos(theta_array[level][layer][i]),
                            radius[level]*math.sin(theta_array[level][layer][i]),
                            z_array[level][layer]] for i in range(n)] for layer in [0, 1]] for level in range(levels)])
    """ Some vertices index offsets:
    1:e vertex in layer=0 (bottom layer) of level: 2 * level * n
    1:e vertex in layer=1 (top layer) of level: 2 * level * n + n
    1: vertex in layer=1 (top layer) of top level: 2 * levels * n - n
    """
    # struts
    struts = np.empty([levels, n, 2], dtype=int)  # shape = (levels, n, layers)
    for level in range(levels):
        struts[level] = [[v0, v1] for v0, v1 in zip(range(2 * level * n, 2 * level * n + n),
                                                    range(2 * level * n + n, 2 * level * n + 2 * n))]
    vertical_tendons = np.empty([levels, n, 2])
    cap_waist_tendons = np.empty([2, n, 2], dtype=int)  # shape = (len[top, bottom], n, len[v0, v1])
    # cap_waist_tendons
    # bottom cap waist tendons
    cap_waist_tendons[0] = [[v0, v1] for v0, v1 in zip(range(n), np.roll(list(range(n)), 1))]
    # top cap waist tendons
    cap_waist_tendons[1] = [[v0, v1] for v0, v1 in
                            zip(range(2 * levels * n - n, 2 * levels * n),
                                np.roll(list(range(2 * levels * n - n, 2 * levels * n)), 1))]
    interface_waist_tendons = np.empty([levels - 1, 2 * n, 2], dtype=int)
    for interface_layer in range(levels-1):
        # interface_layer points to the level below the interface
        # we will create the interface tendons in two steps, first lower level to upper level
        # create list of lower level and list of upper level vertices
        upper_vertices = np.array([v for v in range(2 * interface_layer * n + n, 2 * interface_layer * n + 2 * n)],
                                  dtype=int)
        lower_vertices = np.array([v for v in range(2 * (interface_layer + 1) * n, 2 * (interface_layer + 1) * n + n)],
                                  dtype=int)
        interleaved_vertices = np.empty(upper_vertices.size + lower_vertices.size, dtype=int)
        # interleaved_vertices[0::2] = np.roll(upper_vertices, 1)
        interleaved_vertices[0::2] = upper_vertices
        interleaved_vertices[1::2] = lower_vertices
        print('interleaved_vertices', interleaved_vertices)
        interface_waist_tendons[interface_layer] = [[v0, v1] for v0, v1 in
                                                    zip(interleaved_vertices, np.roll(interleaved_vertices, 1))]
    print('interface_waist_tendons', interface_waist_tendons)
    #     # vertical tendons
    #     tendons.extend([[v0, v1] for v0, v1 in zip(range(n), range(n, 2 * n))])
    cap_waist_list = cap_waist_tendons.reshape(-1, 2)
    interface_waist_list = interface_waist_tendons.reshape(-1, 2)
    # tendons = np.concatenate(cap_waist_list, interface_waist_list, axis=0)
    tendons = np.concatenate((cap_waist_list, interface_waist_list), axis=0)
    if verbose:
        # print('vertices', vertices)
        # print('struts', struts)
        print('tendons', tendons)
        vertex_list = vertices.reshape(-1, 3)
        strut_list = struts.reshape(-1, 2)
        tendon_list = tendons.reshape(-1, 2)
    return vertex_list, strut_list, tendon_list



class Tensegrity:
    def __init__(self, vertices, struts, tendons, name):
        self.vertices = vertices
        self.struts = struts
        self.tendons = tendons
        self.name = name
        self.tendon_force = 1
        self.f_unit_vectors = {}  # key is vertex index: value is scalar strut force
        # self.f_unit_vectors = {}  # key is vertex index: value is list of tendon forces
        """ Build dictionaries that let us find the struts and tendons that are connected to a particular vertex"""
        self.vertex_struts = {}  # key is vertex index: value is strut index
        self.vertex_tendons = {}  # key is vertex index: values is list of tendon indices
        for i, strut in enumerate(self.struts):
            if strut[0] in self.vertex_struts or strut[1] in self.vertex_struts:
                raise Exception('vertex has more than one associated strut')
            else:
                self.vertex_struts[strut[0]] = i
                self.vertex_struts[strut[1]] = i
        for i, tendon in enumerate(self.tendons):
            if tendon[0] in self.vertex_tendons:
                deb0 = self.vertex_tendons[tendon[0]]
                self.vertex_tendons[tendon[0]].extend([i])
            else:
                self.vertex_tendons[tendon[0]] = [i]
            if tendon[1] in self.vertex_tendons:
                self.vertex_tendons[tendon[1]].extend([i])
            else:
                self.vertex_tendons[tendon[1]] = [i]
        """ We will use the convention that strut forces are positive and tendon forces are negative:
        For each vertex the sum of all forces is equal to 0 for a stable balanced tensegrity """
        """ Find the strut force vectors """
        for strut in self.struts:
            for i, vertex_index in enumerate(strut):
                if vertex_index in self.f_unit_vectors:
                    raise Exception('multiple struts are associated with vertex', vertex_index)
                # debug
                v0 = self.vertices[vertex_index]
                v1 = self.vertices[strut[(i + 1) % 2]]
                vec = v0 - v1
                vector = self.vertices[vertex_index] - self.vertices[strut[(i + 1) % 2]]
                self.f_unit_vectors[vertex_index] = [vector / np.linalg.norm(vector)]
        # print(self.f_unit_vectors)
        """ Find the tendon force vectors """
        for tendon in self.tendons:
            for i, vertex_index in enumerate(tendon):
                vector = self.vertices[tendon[(i + 1) % 2]] - self.vertices[vertex_index]
                self.f_unit_vectors[vertex_index].extend([vector/np.linalg.norm(vector)])
        # print(self.f_unit_vectors)

    def tendon_forces(self, verbose=False):
        """ Build the force matrix to solve"""
        # we will start with vertex 0
        # forces = [self.f_unit_vectors[0] + self.f_unit_vectors[0]]
        # a = np.column_stack([self.f_unit_vectors[0], [1, 0, 0]])
        a = np.column_stack([np.transpose(self.f_unit_vectors[0])])
        if self.name == 'kite':
            b = np.array([[0], [0], [1]])
            a[2] = [1, 0, 0]
        elif self.name == 'prism':
            b = np.array([[0], [0], [0], [1]])
            a = np.vstack((a, [1, 0, 0, 0]))
        """ Solve the force matrix"""
        x = np.linalg.solve(a, b)
        total_tendon_forces = 0
        for f_unit_vector, tendon_force in zip(self.f_unit_vectors[0][1:], x[1:]):
            strut_tendon_dot = np.dot(f_unit_vector, self.f_unit_vectors[0][0])
            total_tendon_forces += strut_tendon_dot * tendon_force
        if verbose:
            print('a', a)
            print('b', b)
            print('x', x)
            print('total_tendon_forces', total_tendon_forces)
            if self.name == 'prism':
                # check to see if each of the two waist tendons have the same tension
                if math.isclose(x[1][0], x[2][0], abs_tol=10**-12):
                    print('Prism is stable')
                else:
                    print('Prism is unstable')
        return x[1:]

    def plot(self, struts=True, tendons=True, lateral_f=False, vertex_f=False, axes=False, debug=False):
        """ plots tendons and struts using matplotlib """

        # fig = plt.figure(figsize=(4, 4))
        ax = plt.axes(projection='3d')
        # ax.set_xlim(-5, 5)
        # ax.set_ylim(-5, 5)
        # ax.set_zlim(0, 10)
        x_span = 2*max(vertex[0] for vertex in self.vertices) - 2*min(vertex[0] for vertex in self.vertices)
        y_span = 2*max(vertex[1] for vertex in self.vertices) - 2*min(vertex[1] for vertex in self.vertices)
        z_span = 2*max(vertex[2] for vertex in self.vertices) - 2*min(vertex[2] for vertex in self.vertices)
        span = max(x_span, y_span, z_span)
        ax.set_xlim(-span / 2, span / 2)
        ax.set_ylim(-span / 2, span / 2)
        ax.set_zlim(0, span)
        if not axes:
            ax.set_axis_off()
        if struts:
            for strut in self.struts:
                # x = [vertex[0] for vertex in strut.vertices]
                x = [self.vertices[vertex][0] for vertex in strut]
                y = [self.vertices[vertex][1] for vertex in strut]
                z = [self.vertices[vertex][2] for vertex in strut]
                # y = [vertex[1] for vertex in strut.vertices]
                # z = [vertex[2] for vertex in strut.vertices]
                ax.plot3D(x, y, z, 'red', linewidth=3)
        if tendons:
            for tendon in self.tendons:
                # x = [vertex[0] for vertex in tendon.vertices]
                # y = [vertex[1] for vertex in tendon.vertices]
                # z = [vertex[2] for vertex in tendon.vertices]
                x = [self.vertices[vertex][0] for vertex in tendon]
                y = [self.vertices[vertex][1] for vertex in tendon]
                z = [self.vertices[vertex][2] for vertex in tendon]
                ax.plot3D(x, y, z, 'grey', linewidth=1)
        if debug:
            # for vertex in self.vertices:
            strut = self.struts[0]
            vertex = strut.vertices[0]
            tendon = vertex.tendon_list[0]
            strut_vector = np.array(strut.position_vec(strut.other_vertex_index(vertex)))
            # ax.quiver(*vertex, *strut_vector, color='orange', linewidth=3)
            coords = [val for pair in zip(vertex, strut_vector) for val in pair]
            x = [vertex[0] for vertex in [vertex, strut_vector]]
            y = [vertex[1] for vertex in [vertex, strut_vector]]
            z = [vertex[2] for vertex in [vertex, strut_vector]]
            # y = [vertex[1] for vertex in tendon.vertices]
            # z = [vertex[2] for vertex in tendon.vertices]
            ax.plot3D(x, y, z, 'orange', linewidth=5)
                      # *self.strut_list[0].position_vec(self.strut_list[0].vertices[1]), 'green', linewidth=3)
            tendon_vector = np.array(tendon.position_vec(tendon.other_vertex_index(vertex)))
            # ax.quiver(*vertex, *tendon_vector, color='orange', linewidth=3)
            tendon_strut_component = ((np.dot(tendon_vector, strut_vector) / np.linalg.norm(strut_vector) ** 2)
                                      * strut_vector)
            # ax.quiver(*vertex, *tendon_strut_component, color='green', linewidth=3)
            tendon_lateral_component = tendon_vector - tendon_strut_component
            # ax.quiver(*vertex, *tendon_lateral_component, color='blue', linewidth=3)

            # for vertex in self.vertices[0:1]:
            #     for tendon in vertex.tendon_list[0:1]:
            #         radial_vec, strut_vec, tendon_vec, tendon_strut_component = vertex.tendon_radial_vec(tendon, debug=True)

        plt.show()


def plot_prism1_stability_space(n=3):
    def waist_force_difference(hr_ratio, twist):
        t_forces = Tensegrity(*prism1(n=n, hr_ratio=hr_ratio, twist=twist), name='prism').tendon_forces(verbose=False)
        return t_forces[0][0] - t_forces[1][0]
    # twist_values = np.linspace(math.pi * (1/2 - 1/n - 1/4), math.pi * (1/2 - 1/n + 1/4), 11)
    twist_values = np.linspace(math.pi * (1/2 - 1/n - 1/10), math.pi * (1/2 - 1/n + 1/10), 21)
    hr_ratio_values = np.linspace(0.5, 3, 20)
    z_list = []
    for twist_value in twist_values:
        z_list.extend([waist_force_difference(hr_ratio=hr_ratio_value, twist=twist_value)
                       for hr_ratio_value in hr_ratio_values])
    twist_degrees = [math.degrees(twist_radians) for twist_radians in twist_values]
    x, y = np.meshgrid(hr_ratio_values, twist_degrees)
    z = np.array(z_list).reshape(x.shape)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_wireframe(x, y, z)
    ax.set_xlabel('height:radius ratio')
    ax.set_ylabel('twist (degrees)')
    ax.set_zlabel('waist tendon difference')
    plt.show()


def stabilize_prism1(n=3, verbose=False):
    """ Practice function stabilizes a prism by searching for a twist (alpha) angle that stabilizes the prism """
    min_step_size = math.pi/1000
    max_step_count = 10
    err_tol = 1e-7
    initial_step_size = 10 * min_step_size
    initial_twist = math.pi * (1 / 2 - 1 / n + 1/3)
    prism = Tensegrity(*prism1(n=n, twist=initial_twist), name='prism')
    t_forces = prism.tendon_forces()
    error = t_forces[0][0] - t_forces[1][0]
    error_history = [error]
    step_count = 0
    step_direction = 1
    step = initial_step_size * step_direction
    step_history = [step]
    twist = initial_twist + step
    slope_history = []
    while abs(error) > err_tol and step_count < max_step_count:
        t_forces = Tensegrity(*prism1(n=n, hr_ratio=1, twist=twist), name='prism').tendon_forces(verbose=False)
        error = t_forces[0][0] - t_forces[1][0]
        slope = (error_history[-1] - error) / step
        error_history.append(error)
        slope_history.append(slope)
        step = error / slope
        step_history.append(step)
        twist = twist + step
        step_count += 1
    if verbose:
        print('error history', error_history)
        print('step history', step_history)
        print('slope_history', slope_history)
        print('step_count', step_count)
        if abs(error) < err_tol:
            print('*** Stabilization Successful ***')
        else:
            print('*** Stabilization Failed ***')
        print('twist', math.degrees(twist), '\nwaist tendon force error', error)
    return twist


if __name__ == '__main__':
    # mode = 'kite'
    # mode = 'prism1'  # builds single level prism using n, hr_ratio and twist
    mode = 'prism tower'
    # mode = 'stabilize prism1'
    # mode = 'prism1 stability sweep'  # sweep twist and hr ratio, plot difference in waist tendon forces
    # mode = 'prism'
    # mode = 'unstable prism'
    # mode = 'tp prism'
    # mode = 'tp prism sweep hr ratio'
    if mode == 'kite':
        thing = Tensegrity(*kite(), mode)
        thing.tendon_forces()
    elif mode == 'prism tower':
        strut_count = 3
        thing = Tensegrity(*prism_tower(verbose=True), name='prism')
        thing.plot()
    elif mode == 'stabilize prism1':
        stabilize_prism1(n=5, verbose=True)
    elif mode == 'prism1 stability sweep':
        plot_prism1_stability_space(n=4)
        # # thing = Prism1(n=3, hr_ratio=1, twist=math.pi / 2 - math.pi / 3)
        # n = 3
        # for twist in (math.pi / 2 - math.pi / n) * np.array([0.9, 1, 1.1]):
        #     for hr_ratio in [1, 2, 3]:
        #         thing = Tensegrity(*prism1(n=3, hr_ratio=hr_ratio, twist=twist), name='prism')
        #         t_forces = thing.tendon_forces(verbose=False)
        #         # print('t_forces', t_forces)
        #         print('hr_ratio', hr_ratio, 'waist force difference', t_forces[0][0] - t_forces[1][0])
        # # thing.plot()
    elif mode == 'prism1':
        # thing = Prism1(n=3, hr_ratio=1, twist=math.pi / 2 - math.pi / 3)
        # thing = Tensegrity(*prism1(n=3, hr_ratio=2, twist=math.pi / 2 - math.pi / 3), name='prism')
        strut_count = 6
        thing = Tensegrity(*prism1(n=strut_count, hr_ratio=2, twist=math.pi / 2 - math.pi / strut_count), name='prism')
        thing.plot()
        thing.tendon_forces(verbose=True)
    elif mode == 'prism':
        thing = Tensegrity(*prism_n3(), mode)
        thing.tendon_forces(verbose=True)
    elif mode == 'unstable prism':
        thing = Tensegrity(*prism_n3(alpha=5 * math.pi / 6 + 0.01), name='prism')
        # thing = Tensegrity(*prism(alpha=5*math.pi/6), name='prism')
        thing.tendon_forces(verbose=True)
    elif mode == 'tp prism':
        """ Use TgyPar code to generate tensegrity and find forces """
        tp_thing = tp.PrismTower(n=3, levels=1, radii=[5, 5], heights=[10])
        thing = Tensegrity(tp_thing.get_vertices(), tp_thing.get_struts(), tp_thing.get_tendons(), name='prism')
        thing.tendon_forces()
    elif mode == 'tp prism sweep hr ratio':
        """ Use TgyPar code to generate tensegrity and find forces 
            Sweep height to radius ratio and write csv file with tendon forces as a fraction of strut force """
        force_file = open('C:/Users/guyde/PycharmProjects/TgyCalc/forces.csv', 'w')
        force_file.write('Ratio, Waist Tendon Tension, Vertical Tendon Tension, Waist:Vertical Tendon Tension,' +
                         ' Strut Tension \n')
        for ratio in [0.1, 0.25, 0.5, 0.85, 1, 2, 3, 5, 10, 20, 50, 100]:
            tp_thing = tp.PrismTower(n=3, levels=1, radii=[5, 5], heights=[5 * ratio])
            thing = Tensegrity(tp_thing.get_vertices(), tp_thing.get_struts(), tp_thing.get_tendons(), name='prism')
            tendon_forces = thing.tendon_forces()
            print(tendon_forces)
            print(tendon_forces[0][0])
            print(tendon_forces[2][0])
            force_file.write(str(ratio) + ', ' + str(tendon_forces[0][0]) + ', ' + str(tendon_forces[2][0]) + ', ' +
                             str(tendon_forces[0][0] / tendon_forces[2][0]) + ', 1 ' + str('\n'))
        force_file.close()

