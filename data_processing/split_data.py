f = open('PINTS03-06', 'r')
resources = {}
lines_read = 0
files = [open('2003.dat', 'w'),
         open('2004.dat', 'w'),
         open('2005.dat', 'w'),
         open('2006.dat', 'w'),
         ]
prev = ''
for line in f:
    year = int(line[3]) - 3
    split = line.split('\t')
    new_line = line[0:10] + '\t' + split[1] + '\t' + split[2] + '\n'
    if new_line != prev:
        files[year].write(new_line)
        prev = new_line

for out_file in files:
    out_file.close()
