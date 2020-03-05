# Gaussian Elimination
A program designed to compute a system of linear equations using Gaussian's naive solution or scaled-partial pivoting.

How to use:
  Have a file with the following format
  1) the first line will contain the amount of variables in the system (n)
  2) for (n) lines, write the coefficients of each respective variable
  3) Lastly, on the final line, have the constants
  
  Example of conversion from system to file:
  
  System:
  
      1x1 + 2x2 + 3x3 = 1
      4x1 + 5x2 + 6x3 = 1
      7x1 + 8x2 + 9x3 = 1
  File:
  
     3
     1 2 3
     4 5 6
     7 8 9
     1 1 1
    
Check out the 'exampleFile.lin' to get a clear look at the file.

How to run:
    Open terminal/cmd and place the program within the same folder as your system of equations file. 
    Compile the java code using:
    
    javac gaussian.java
    
    
   Run the code by using either of the following:
    
      java gaussian sys1.lin
      java gaussian --spp sys1.lin
      
    The first version will run the naive gaussian elimination.
    The second version will run the scaled-partial pivoting gaussian elimination.
 
Output:
  The result of running this program will be the creation of a new file within the directory. This file will
  have the same basename as the input file, except the extension will be '.sol'
  
Example of output file created:
    
    Input file: sys1.lin
    Output file: sys1.sol

# What I Learned
  • Keeping track of many for loops to be able to manipulate the system of equations
  
  • Improved my ability to search and find bugs
  
  • Improved familiarity with Gaussian Elimination
