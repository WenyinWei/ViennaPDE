# Executable
Now the executable is put at 
/PDE/build/src/pde_test,
which is a monotone Godunov scheme to solve the Burger's equation with i.c. as 1/3+2/3 sin(x).

# Compile
It is compiled with a NVIDIA CUDA driver and powered by ViennaCL openCL library, so I am afraid that on other machines users may need to reconfigure and rebuild for his/her system.
## CPP
One may be able to change the parameters he want to change at /PDE/src/pde_test.cpp

## ViennaPDE
If Prof. Yu or Prof. Du wants to check the concrete scheme steps, they are put at the hpp files which compose the ViennaPDE library. ViennaPDE is written by Wenyin during his study of PDE at Tsinghua University based on the ViennaCL library.

## Analyze
The following is the python script to analyze the result.

%matplotlib inline
import matplotlib.pyplot as plt
import csv

x=[]
y=[]
n=0
dx=0.1

with open('30.csv', 'r') as csvfile:
    [plots]= csv.reader(csvfile, delimiter=',')
    for pointdata in plots:
        n = n + 1
        x.append(n*dx)
        y.append(float(pointdata))


plt.plot(x,y, marker='o')

plt.title('Data from the CSV File: Burger Equation ')

plt.xlabel('Number of People')
plt.ylabel('Expenses')

plt.show()