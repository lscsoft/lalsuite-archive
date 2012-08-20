import sys
if len(sys.argv) is not 3:
    sys.stderr.write("Usage: down_sample.py [input file] [rate]!")
    exit(1)

in_file=open(sys.argv[1],'r')
down_sample_factor=int(sys.argv[2])

for i,line in enumerate(in_file):
    if i%down_sample_factor==0:
        sys.stdout.write(line)
