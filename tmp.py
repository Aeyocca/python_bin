#!/bin/python
#Alan E. Yocca
#just seeing if things print out to slurm script

import time


def print_things(variable = "", add_maf = True):
  if add_maf:
    print("true")
  else:
    print("false")
  print("In print things: %s" % (variable))

if __name__ == "__main__":


  for i in range(0,10):
    print_things(i)
#    time.sleep(10)
  print("Finished")
