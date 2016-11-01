#Rydberg Molecule Potential Programme (RyMoPP)


Ke Liao 2016


Any use of this code for publication should
cite this project or the Manual.pdf.
Any commercial use should follow the
terms of GPL3 license.

##Important Note


Anyone with a developer permission to this 
repository is welcome to create your own 
branch or make a fork of this project. 


However, the master branch is under protection
just in case any developer may unintentionally
push bad files or delete
some important files. 


After you create your own branch or fork 
your own project, you are free to modify
the code in your branch. 


**Important**


Learn how to use git first before you
push any changes to your repository on GitLab


##Compile


gcc 4.8 is required. Other versions (either
old or newer ones) may not work properly.


On linux or Mac, use the following command
 in terminal to check the version of gcc:


   `gcc -v`


If it is something like 


  `gcc version 4.8.x`


then you can proceed. If not try to reinstall
the 4.8 version.


Go to ./Code, run the following script


  `./compile.sh`


make sure that this script has the x permision
on your machine. If not, use


  `chmod +x compile.sh`


to assign then running permission to it.


##Run


After a successful compile, an excutable file


  `RyMoPP`


will be generated. Just run it with


  `./RyMoPP`


If you are on a server and want to keep the
programme running after you log off, try this


  `nohup ./RyMoPP`
