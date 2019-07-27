from numpy import zeros, sqrt, arange, linspace, sum, power, inf, nan
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


PATH_1 = 1
PATH_3_1 = sqrt(10)
PATH_2_1 = sqrt(5)
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


def get_backtraced_path(local_distance_matrix):
    x = local_distance_matrix.shape[0]
    y = local_distance_matrix.shape[1]

    cumulated_distance_matrix = zeros((x,y))
    cumulated_distance_matrix[:] = inf
    cumulated_distance_matrix[0, 0] = local_distance_matrix[0, 0]

    backtrace_matrix = zeros((x,y))

    open_cells = []
    open_cells.insert(0, [0,0])

    while open_cells:
        o = open_cells.pop()

        if o[0] >= x or o[1] >= y:
            continue

        #path 4
        if o[0] < x-1 and o[1] < y-1:
            d = cumulated_distance_matrix[o[0], o[1]] + (local_distance_matrix[o[0]+1, o[1]+1] * PATH_1_1)
            if d < cumulated_distance_matrix[o[0]+1, o[1]+1]:
                cumulated_distance_matrix[o[0]+1, o[1]+1] = d
                backtrace_matrix[o[0]+1, o[1]+1] = DIAGONAL
                if not [o[0]+1, o[1]+1] in open_cells:
                    open_cells.insert(0, [o[0]+1, o[1]+1])

        #path 1
        if o[0] < x-1:
            d = cumulated_distance_matrix[o[0], o[1]] + local_distance_matrix[o[0]+1, o[1]]
            if d < cumulated_distance_matrix[o[0]+1, o[1]]:
                cumulated_distance_matrix[o[0]+1, o[1]] = d
                backtrace_matrix[o[0]+1, o[1]] = UPWARD_1
                if not [o[0]+1, o[1]] in open_cells:
                    open_cells.insert(0, [o[0]+1, o[1]])

        #path 7
        if o[1] < y-1:
            d = cumulated_distance_matrix[o[0], o[1]] + local_distance_matrix[o[0], o[1]+1]
            if d < cumulated_distance_matrix[o[0], o[1]+1]:
                cumulated_distance_matrix[o[0], o[1]+1] = d
                backtrace_matrix[o[0], o[1]+1] = RIGHT_1
                if not [o[0], o[1]+1] in open_cells:
                    open_cells.insert(0, [o[0], o[1]+1])

    path = [[x-1, y-1]]
    while not [0, 0] in path:
        #path 4
        if backtrace_matrix[path[0][0], path[0][1]] == DIAGONAL:
            path.insert(0, [path[0][0]-1, path[0][1]-1])
        #path 1
        if backtrace_matrix[path[0][0], path[0][1]] == UPWARD_1:
            path.insert(0, [path[0][0]-1, path[0][1]])
        #path 7
        if backtrace_matrix[path[0][0], path[0][1]] == RIGHT_1:
            path.insert(0, [path[0][0], path[0][1]-1])

    plt.subplot(1,3,1)
    plt.pcolormesh(local_distance_matrix)
    plt.subplot(1,3,2)
    plt.pcolormesh(cumulated_distance_matrix)
    plt.subplot(1,3,3)
    plt.pcolormesh(backtrace_matrix)
    plt.show()

    return path


def adtw(des_array, fit_array, start_padding=5, middle_padding=0.50):
    height = des_array.shape[1]
    width = fit_array.shape[1]

    distance_matrix = zeros((height, width))
    distance_matrix[:,:] = nan
    points = linspace(0, width, width)
    upper, lower = get_area(width, height, start_padding, middle_padding)
    upper = upper(points).astype('int')
    lower = lower(points).astype('int')

    for i in range(0, width):
        for j in range(lower[i], upper[i]):
            distance_matrix[j,i] = sqrt(sum(power((des_array.T[j]-fit_array.T[i]), 2)))


    get_backtraced_path(distance_matrix)

    return distance_matrix
