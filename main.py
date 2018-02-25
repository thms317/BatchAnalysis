import measurements_20180201 as meas
import batch_analysis as ba
import open_data as open

measurements, data_path, subfolder, analysis_path = meas.measurements()

for n, i in enumerate(measurements):
    print("Processing measurement " + str(n) + " of " + str(len(measurements)))
    time, force, beads, diff_magnet, headers, data_lines, title, file_name = open.open_data(n, measurements, data_path)
    ba.batch_analysis(beads, data_lines, headers, time, diff_magnet, force, title, file_name, analysis_path)