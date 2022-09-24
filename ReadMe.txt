The black hole imaging codes are developed by us: Tsaichen Lee, Jiaming Lu and Yuchen Wang
as our C161 project in 2022 Spring.

Compiled with python 3.9

Kerr_trajectory.py calculates and graphs the trajectory of photons. Can be run directly. Takes 1-2 minutes.

Kerr_imaging.py process a picture to see what happends if we put a black hole in front of it.
You need to provide it with a picture in line #246, here it uses 'berkeley.png'. Put it in the same folder with the python code. Change resolution at #243. For grid = 256 we get a image of resolution 256*256. Which takes about 18 hours. One can change the grid numer to make it faster.
