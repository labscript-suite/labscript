import sys
import time
myfile = open(sys.argv[1])
start_time = time.time()
x = myfile.read()
myfile.close()
print "time taken to read:", time.time() - start_time
