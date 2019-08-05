from numpy import zeros, sqrt, arange, linspace, sum, power, inf, nan, array
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


PATH_1 = 1
PATH_3_1 = sqrt(10)/3
PATH_2_1 = sqrt(5)/2
PATH_1_1 = sqrt(2)

UPWARD_1 = 1
UPWARD_3_1 = 2
UPWARD_2_1 = 3
DIAGONAL = 4
RIGHT_2_1 = 5
RIGHT_3_1 = 6
RIGHT_1 = 7


def get_area(width, height, start_padding, middle_padding):
    horizontal_d = height/width
    vertical_d = width/height

    anchor_y = height/3
    anchor_x = width/3

    anchor_y_2 = 2*height/3
    anchor_x_2 = 2*width/3

    if width >= height:
        middle_padding = middle_padding + ((1-middle_padding)*(1-(height/width)))
        slope_x = (anchor_y * horizontal_d * middle_padding)
        slope_y = (anchor_x * horizontal_d * middle_padding)
    else:
        middle_padding = middle_padding + ((1-middle_padding)*(1-(width/height)))
        slope_x = (anchor_y * vertical_d * middle_padding)
        slope_y = (anchor_x * vertical_d * middle_padding)

    downslope_x = anchor_x + slope_x
    downslope_y = anchor_y - slope_y
    upslope_x = anchor_x - slope_x
    upslope_y = anchor_y + slope_y

    downslope_x_2 = anchor_x_2 + slope_x
    downslope_y_2 = anchor_y_2 - slope_y
    upslope_x_2 = anchor_x_2 - slope_x
    upslope_y_2 = anchor_y_2 + slope_y

    upper = interp1d([0, upslope_x, upslope_x_2, width-start_padding, width],
                    [0+start_padding, upslope_y, upslope_y_2, height, height], kind='linear')
    lower = interp1d([0, start_padding, downslope_x, downslope_x_2, width],
                    [0, 0, downslope_y, downslope_y_2, height-start_padding], kind='linear')

    return upper, lower


def path_step(x, y, upward, right, cost, tag, local_matrix, cumulated_matrix, backtrace_matrix, open_cells):
    max_x = local_matrix.shape[0]
    max_y = local_matrix.shape[1]

    if x + upward >= max_x or y + right >= max_y:
        return

    if local_matrix[x+upward, y+right] == nan:
        return

    s_upward = upward
    if s_upward > 1:
        s_upward = 1

    s_right = right
    if s_right > 1:
        s_right = 1

    d = cumulated_matrix[x, y] + sum(local_matrix[x+s_upward:x+upward+1, y+s_right:y+right+1]) * cost
    if d < cumulated_matrix[x+upward, y+right]:
        cumulated_matrix[x+upward, y+right] = d
        backtrace_matrix[x+upward, y+right] = tag
        if not [x+upward, y+right] in open_cells:
            open_cells.insert(0, [x+upward, y+right])


def get_path(backtrace_matrix):
    x = backtrace_matrix.shape[0]
    y = backtrace_matrix.shape[1]

    path = [[x-1, y-1]]
    while not [0, 0] in path:
        #path 4
        if backtrace_matrix[path[0][0], path[0][1]] == DIAGONAL:
            path.insert(0, [path[0][0]-1, path[0][1]-1])
        #path 1
        elif backtrace_matrix[path[0][0], path[0][1]] == UPWARD_1:
            path.insert(0, [path[0][0]-1, path[0][1]])
        #path 7
        elif backtrace_matrix[path[0][0], path[0][1]] == RIGHT_1:
            path.insert(0, [path[0][0], path[0][1]-1])
        #path 3
        elif backtrace_matrix[path[0][0], path[0][1]] == UPWARD_2_1:
            path.insert(0, [path[0][0]-1, path[0][1]-0])
            path.insert(0, [path[1][0]-2, path[1][1]-1])
        #path 5
        elif backtrace_matrix[path[0][0], path[0][1]] == RIGHT_2_1:
            path.insert(0, [path[0][0]-0, path[0][1]-1])
            path.insert(0, [path[1][0]-1, path[1][1]-2])
        #path 2
        elif backtrace_matrix[path[0][0], path[0][1]] == UPWARD_3_1:
            path.insert(0, [path[0][0]-1, path[0][1]-0])
            path.insert(0, [path[1][0]-2, path[1][1]-0])
            path.insert(0, [path[2][0]-3, path[2][1]-1])
        #path 6
        elif backtrace_matrix[path[0][0], path[0][1]] == RIGHT_3_1:
            path.insert(0, [path[0][0]-0, path[0][1]-1])
            path.insert(0, [path[1][0]-0, path[1][1]-2])
            path.insert(0, [path[2][0]-1, path[2][1]-3])
        else:
            break

    return path


def get_backtrace_matrix(local_matrix):
    x = local_matrix.shape[0]
    y = local_matrix.shape[1]

    cumulated_matrix = zeros((x,y))
    cumulated_matrix[:] = inf
    cumulated_matrix[0, 0] = local_matrix[0, 0]

    backtrace_matrix = zeros((x,y))

    open_cells = []
    open_cells.insert(0, [0,0])

    while open_cells:
        o = open_cells.pop()

        if o[0] >= x or o[1] >= y:
            continue

        if local_matrix[o[0], o[1]] == nan:
            continue

        #path 4
        path_step(o[0], o[1], 1, 1, PATH_1_1, DIAGONAL, local_matrix, cumulated_matrix, backtrace_matrix, open_cells)

        #path 1
        #path_step(o[0], o[1], 1, 0, 1, UPWARD_1, local_matrix, cumulated_matrix, backtrace_matrix, open_cells)

        #path 7
        #path_step(o[0], o[1], 0, 1, 1, RIGHT_1, local_matrix, cumulated_matrix, backtrace_matrix, open_cells)

        #path 3
        path_step(o[0], o[1], 2, 1, PATH_2_1, UPWARD_2_1, local_matrix, cumulated_matrix, backtrace_matrix, open_cells)

        #path 5
        path_step(o[0], o[1], 1, 2, PATH_2_1, RIGHT_2_1, local_matrix, cumulated_matrix, backtrace_matrix, open_cells)

        #path 2
        path_step(o[0], o[1], 3, 1, PATH_3_1, UPWARD_3_1, local_matrix, cumulated_matrix, backtrace_matrix, open_cells)

        #path 6
        path_step(o[0], o[1], 1, 3, PATH_3_1, RIGHT_3_1, local_matrix, cumulated_matrix, backtrace_matrix, open_cells)

    return backtrace_matrix


def adtw(des_array, fit_array, start_padding=5, middle_padding=0.50, plot_distance=False):
    height = des_array.shape[0]
    width = fit_array.shape[0]

    distance_matrix = zeros((height, width))
    distance_matrix[:,:] = nan
    points = linspace(0, width, width)
    upper, lower = get_area(width, height, start_padding, middle_padding)
    upper = upper(points).astype('int')
    lower = lower(points).astype('int')

    for i in range(0, width):
        for j in range(lower[i], upper[i]):
            distance_matrix[j,i] = sqrt(sum(power((des_array[j]-fit_array[i]), 2)))


    backtrace_matrix = get_backtrace_matrix(distance_matrix)
    path = get_path(backtrace_matrix)

    return array(path)
