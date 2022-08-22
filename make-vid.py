import os

prefix = 'top'

n = 0
for i in range(0,110,10):
    
    in_file = prefix + '{0:03}'.format(i) + '.png'
    out_file = prefix + '{0:04}'.format(n) + '.png'
    os.rename(in_file, out_file)
    n += 1

command = 'ffmpeg -start_number 0000 -i top%04d.png output.mp4'
os.system(command)

