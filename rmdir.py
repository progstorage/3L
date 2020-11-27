from sys import argv
import os


def f():
	try:
		for file in os.listdir("sorted_points"):
			os.remove("sorted_points/"+file)
		os.rmdir("sorted_points")
		for file in os.listdir("levels"):
				os.remove("levels/"+file)
		os.rmdir("levels")
		for file in os.listdir("coefs"):
				os.remove("coefs/"+file)
		os.rmdir("coefs")
		return
	except:
		return
	
 
if __name__ == '__main__':
	f()