import os

def measurements():

    # paths
    data_path = "C:\\Users\\brouw\\Desktop\\"
    subfolder = "mono 601 - batch_test2"
    analysis_path = data_path + subfolder + "\\"
    if not os.path.exists(analysis_path):
        os.makedirs(analysis_path)

    measurements = []

    # measurements

    measurements.append(['180201', '006'])
    measurements.append(['180201', '007'])
    measurements.append(['180201', '008'])
    measurements.append(['180201', '009'])
    measurements.append(['180201', '011'])
    measurements.append(['180201', '012'])
    measurements.append(['180201', '013'])
    measurements.append(['180201', '014'])
    measurements.append(['180201', '015'])

    return measurements, data_path, subfolder, analysis_path