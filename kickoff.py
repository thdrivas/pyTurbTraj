#!/usr/bin/env python
import os
import time
import subprocess
import sys

def main():
	proc = restart_process()
	while True:
		time.sleep(15)

		if proc.poll() is None:
			print "Process is ok!"
		else:
			print "Process died! Restart it" 
			proc = restart_process()


def restart_process():
	proc = subprocess.Popen([sys.executable, "./particles_script.py", r"&"], stdout=sys.stdout)
	return proc


if __name__ == "__main__":
	main()