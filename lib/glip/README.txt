
   █████████  █████       █████ ███████████ 
  ███▒▒▒▒▒███▒▒███       ▒▒███ ▒▒███▒▒▒▒▒███
 ███     ▒▒▒  ▒███        ▒███  ▒███    ▒███
▒███          ▒███        ▒███  ▒██████████ 
▒███    █████ ▒███        ▒███  ▒███▒▒▒▒▒▒  
▒▒███  ▒▒███  ▒███      █ ▒███  ▒███        
 ▒▒█████████  ███████████ █████ █████       
  ▒▒▒▒▒▒▒▒▒  ▒▒▒▒▒▒▒▒▒▒▒ ▒▒▒▒▒ ▒▒▒▒▒        

*** Graph-based Library for Image Processing ***


Instructions for use:

1) Copy "glip" files to any folder on your computer.
For example, consider the folder: /home/user/prog/lib/glip

2) Create an environment variable GLIP_DIR and place it in .bashrc:
export GLIP_DIR=/home/user/prog/lib/glip

3) Enter the "glip" folder and run:
make clean
make

4) If the following error occurs "fatal error: zlib.h: No such file or directory", then you have to install the package zlib1g-dev from zlib and repeat step 3 again.

5) Finally, go to the folder of the program that uses glip and run:
make clean
make
