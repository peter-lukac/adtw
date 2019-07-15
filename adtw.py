from numpy import zeros, sqrt, arange, linspace, sum, power, inf, nan
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


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


def adtw(des_array, fit_array, start_padding=5, middle_padding=0.30):
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

    return distance_matrix
