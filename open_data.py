import functions as func
import numpy as np


def open_data(n, measurements, data_path):

    # open data
    sub = measurements[n]
    file_location = data_path + str(sub[0]) + "\\"
    file_name = "data_" + str(sub[1])
    file_extension = ".dat"
    file_all = file_location + file_name + file_extension

    # title - extract date of the measurements
    title = int(func.get_num(file_location))

    # headers
    f = open(file_all, 'r')
    headers_lines = f.readlines()[0]
    f.close()

    headers = headers_lines.split('\t')

    # data
    f = open(file_all, 'r')
    data_lines = f.readlines()[1:]
    f.close()

    # number of beads
    file_name = "data_" + str(sub[1])
    file_extension = ".log"
    file_all = file_location + file_name + file_extension

    # find "number of beads" column
    f = open(file_all, 'r')
    beads = f.readlines()[9]
    f.close()
    beads = int(func.get_num(beads))

    # calculate force from magnet position
    magnet = []
    time = []
    force = []
    for n, x in enumerate(data_lines):
        magnet.append(float(x.split()[headers.index('Stepper shift (mm)')]))
        time.append(float(x.split()[headers.index('Time (s)')]))
    magnet = np.array(magnet)
    for i in magnet:
        force.append(func.calc_force(i))

    # calculating the first derivative of magnet to discriminate in pull/release-curve
    dx = np.diff(time)
    dy = np.diff(magnet)
    diff_magnet = np.array(np.divide(dy, dx))

    return time, force, beads, diff_magnet, headers, data_lines, title, file_name