import os,sys
import subprocess

tem = open(sys.argv[1]).readlines()
for line in tem:
    subprocess.Popen(line + ' &')

