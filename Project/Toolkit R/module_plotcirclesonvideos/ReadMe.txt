process_in_one_click.bat -> this is batch file that runs the .py file in this folder
process_timeseries_to_add_circles.py -> this the python file that takes input from timeseries folder, and original_videos folder, to output new videos with circles in the output_... folder
timeseries -> this is has the csv's with the timeseries PER VIDEO, such that filename is same as video (minus the extension name), and there is an 'frame' and 'x' and 'y' and 'radius' variable which is used to draw the circle
original_videos -> this contains the videos, need to be the same name as the timeseries (minus the extension name)
output_videos_with_circles -> this has the outputted videos from the .py script with circles added
