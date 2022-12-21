import math
from numpy import concatenate, shape, array, transpose, linalg, zeros, ones, full
from xalglib import xalglib


def calculate_t3_prism_tensions(alpha):
    vertices = array([[math.cos(0), math.sin(0), 0],
                      [math.cos(2*math.pi/3), math.sin(2*math.pi/3), 0],
                      [math.cos(4 * math.pi / 3), math.sin(4 * math.pi / 3), 0],
                      [math.cos(alpha), math.sin(alpha), 1],
                      [math.cos(alpha + 2 * math.pi / 3), math.sin(alpha + 2 * math.pi / 3), 1],
                      [math.cos(alpha + 4 * math.pi / 3), math.sin(alpha + 4 * math.pi / 3), 1]])
    #these numbers are indicies in the vertices array for the vertices on either end of the member
    compression_members = array([[0,3],
                                 [1,4],
                                 [2,5]])
    tension_members = array([[0,1],
                              [1,2],
                              [2,0],
                              [3,4],
                              [4,5],
                              [5,3],
                              [1,3],
                              [2,4],
                              [0,5]])
    return solve_tensegrity_tensions(compression_members, tension_members, vertices)


def calculate_kite_tensions():
    vertices = array([[0, 1], [1, 0], [0, -2], [-1, 0]])
    #these numbers are indicies in the vertices array for the vertices on either end of the member
    compression_members = array([[2, 0],
                                 [3, 1]])
    tension_members = array([[0, 1],
                             [1, 2],
                             [2, 3],
                             [3, 0]])
    return solve_tensegrity_tensions(compression_members, tension_members, vertices)


def solve_tensegrity_tensions(compression_members, tension_members, vertices, verbose=False):
    matrix = build_tensegrity_matrix(compression_members, tension_members, vertices)
    if verbose:
        print('* alglib matrix *\n', matrix)
    return find_positive_solutions(matrix)


# Returns a matrix that takes in a vector of the scalar tensions / compressions in each member, and outputs the net
# force for each vertex. For a valid tensegrity structure the vector of tensions / compressions should be >= 0
# elementwise
def build_tensegrity_matrix(compression_members, tension_members, vertices):
    compression_block = build_tensegrity_matrix_block(compression_members, vertices, 1)
    tension_block = build_tensegrity_matrix_block(tension_members, vertices, -1)
    return concatenate((compression_block, tension_block), axis=1)


# Sign = +1 for compression, -1 for tension
def build_tensegrity_matrix_block(members, vertices, sign, verbose=False):
    # TODO: this whole method seems like it could be cleaned up significantly
    result = None
    for member in members:
        force_direction = get_unit_force_vector(member, vertices)
        column = array([])
        for i, vertex in enumerate(vertices):
            vertex_force = zeros(shape(vertex))
            if i == member[0]:
                vertex_force = sign * force_direction;
            if i == member[1]:
                vertex_force = sign * -force_direction;
            column = concatenate((column, vertex_force), axis=0)
        # TODO: These two statements are especially gross.
        column_as_column_vector = transpose(array([column]))
        result = concatenate((result, column_as_column_vector),
                             axis=1) if result is not None else column_as_column_vector
    return result


def get_unit_force_vector(member, vertices):
    force_vector = vertices[member[1]] - vertices[member[0]]
    return force_vector / linalg.norm(force_vector)


    # This method finds an x s.t.:
    # matrix*x = 0 , and
    # x >= 0 elementwise, and
    # x's elements are reasonably size (~1).
def find_positive_solutions(matrix):
    # TODO: This method is pretty gross and uses
    #  some heavy machinery that doesn't seem necessary...
    state = xalglib.minlpcreate(shape(matrix)[1])
    xalglib.minlpsetcost(state, list(full(shape(matrix)[1], -1)))
    # Totally arbitrary value; this means the highest allowed ratio of max force / min force is 1:0.1
    xalglib.minlpsetbcall(state, 0.1, 1)
    # TODO: Could probably use the sparse version for large tensegrity structures if
    #  performance becomes an issue.
    range_zeros = list(zeros((shape(matrix)[0])))
    xalglib.minlpsetlc2dense(state, matrix.tolist(), range_zeros, range_zeros)  # here is where we can set error limits
    xalglib.minlpoptimize(state)
    x, rep = xalglib.minlpresults(state)
    print('alglib termination type: ', rep.terminationtype)
    if rep.terminationtype > 0:
        return x
    else:
        #TODO: Probably is a more elegent way to handle this case...
        # raise warnings("No solution found for tensegrity structure")
        print("*** No solution found for tensegrity structure ***")


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print('kite tensions')
    print(calculate_kite_tensions())
    print('valid t3 prism tensions')
    print(calculate_t3_prism_tensions(5*math.pi/6))
    print('invalid t3 prism tensions')
    print(calculate_t3_prism_tensions(0))
# main.py
#Displaying main.py.