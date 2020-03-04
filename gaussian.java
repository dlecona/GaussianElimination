import java.util.*;
import java.io.*;

public class gaussian{
	
	//creating variables for globally used variables
	private static int numOfVars;
	private static Vector<Float> sol;
	private static float[][] coeff;
	private static Vector<Float> constants;
	private static String file;
	private static ArrayList<String> name;
	
	public static void main(String[] args) throws FileNotFoundException{
    	//Asking user for file
		
    	if(args.length == 1) {
		//creating a file object to be able to read file data
    		file = args[0];
	        File data = new File(args[0]);
	        Scanner sc = new Scanner(data);
	        //checking to see if the file exists
	        if(data.exists()) {
	        	//creating a scanner class to read the file
	            String num = sc.next();
				numOfVars = Integer.parseInt(num);
	            float[][] coeff = new float[numOfVars][numOfVars];
	            
	            //grabbing the data from the file.
	            for(int i = 0; i < numOfVars; i++)
	            	for(int j = 0; j < numOfVars; j++)
	            		coeff[i][j] = sc.nextFloat();

	            //creating a vector to hold the constants
	            Vector<Float> constants = new Vector<>(numOfVars);

	            //adding values to the vector
	            while(sc.hasNextFloat()){
					constants.add(sc.nextFloat());
				}
	            sc.close();
	            
	            //calling the naive Gaussian function
	            naiveGaussian(coeff, constants);
	        }
	        else {
	        	System.out.println("File does not exist, make sure the file is named \"sys1.lin\".");
	        }
    	}
    	else {
    		file = args[1];
	        File data = new File(file);
	        Scanner sc = new Scanner(data);
    		if(data.exists()) {
	        	//creating a scanner class to read the file
	            String num = sc.next();
				numOfVars = Integer.parseInt(num);
	            float[][] coeff = new float[numOfVars][numOfVars];
	            
	            //grabbing the data from the file.
	            for(int i = 0; i < numOfVars; i++)
	            	for(int j = 0; j < numOfVars; j++)
	            		coeff[i][j] = sc.nextFloat();

	            //creating a vector to hold the constants
	            Vector<Float> constants = new Vector<>(numOfVars);

	            //adding values to the vector
	            while(sc.hasNextFloat()){
					constants.add(sc.nextFloat());
				}
	            sc.close();
	            
	            //calling the naive Gaussian function
	            SPPGaussian(coeff, constants);
    		}
    	}
              try {
				file = file.substring(0, file.lastIndexOf('.'));
				FileWriter myWriter = new FileWriter(file + ".sol");
				for(int i = 0; i < numOfVars; i++)
					myWriter.write(Float.toString(sol.get(i)) + " ");
	            myWriter.close();
	            System.out.println("Solution file created successfully");
			} catch (IOException e) {e.printStackTrace();}
	}

	//Naive Guassian elimination set up
    public static void naiveGaussian(float[][] coeff, Vector<Float> constants){
    	sol = new Vector<>();
    	sol.setSize(numOfVars);
    	fwdElimination(coeff, constants);
    	backSubst(coeff, constants, sol);
    }
    
    //Forward elimination 
	private static void fwdElimination(float[][] coeff, Vector<Float> constants) {
		//triple nested for loop to get the pivot, row, and certain coefficient
		for(int i = 0; i < numOfVars; i++) {
			for(int j = i+1; j < numOfVars; j++) {
				//This grabs the multiple needed to make the next rows 0
				float mult = coeff[j][i] / coeff[i][i];
				
				//THE PLANTED BUG!!! ITS SUPPOSED TO BE STARTING AT 0, NOT AT 1
				for(int k = i; k < numOfVars; k++) {
					coeff[j][k] = coeff[j][k] - (mult * coeff[i][k]);
				}
				//updating the constants to be calculated with the rest of the row.
				float result = constants.get(j) - (mult * constants.get(i));
				constants.set(j, result);
			}
		}
	}
    
	private static void backSubst(float[][] coeff, Vector<Float> constants, Vector<Float> sol) {
		//setting size of solution vector
		sol.setSize(numOfVars);
		
		//loading values into the solution
		for(int i = 0; i < numOfVars; i++)
			sol.set(i, constants.get(i) / coeff[i][i]);
		
		//back substituting the variables to be able to come to a final solution
		for(int i = numOfVars-1; i >= 0; i--) {
			float sum = constants.get(i);
			for(int j = i+1; j < numOfVars; j++) {
				sum = sum - (coeff[i][j] * sol.get(j));
			}
			
			//final solution
			sol.set(i, sum/coeff[i][i]);
		}
	}
	
	//Scaled Partial Pivoting Guassian Method!
	private static void SPPGaussian(float[][] coeff, Vector<Float> constants) {
		sol = new Vector<>();
    	sol.setSize(numOfVars);
		Vector<Integer> index = new Vector<>();
		index.setSize(numOfVars);
		sol.setSize(numOfVars);
		//initializing index vector
		for(int i = 0; i < numOfVars; i++)
			index.set(i, i);
		
		//calling the methods to execute scaled partial pivoting guassian method
		SPPFwdElimination(coeff, constants, index);
		SPPBackSubstitution(coeff, constants, sol, index);
	}

	private static void SPPFwdElimination(float[][] coeff, Vector<Float> constants, Vector<Integer> index) {
		Vector<Float> scaling = new Vector<>();
		scaling.setSize(numOfVars);
		
		float max = 0;
		for(int i = 0; i < numOfVars; i++) {
			max = 0;
			//finding coeffient with the greatest absolute value
			for(int j = 0; j < numOfVars; j++)
				max = Math.max(max, Math.abs(coeff[i][j]));
			scaling.set(i, max);
		}
		
		//finding the max of the Xn coefficients.
		for(int i = 0; i < numOfVars-1; i++) {
			float rmax = 0;
			int maxIndex = i;
			for(int j = i; j < numOfVars; j++) {
				// ratio of coefficient to scaling factor
				float ratio = Math.abs(coeff[index.get(j)][i] / scaling.get(index.get(j)));
				if(ratio > rmax) {
					rmax = ratio;
					maxIndex = j;
				}
			}

			//Swapping the max value for Xn to be the top most equation.
			int temp = index.get(maxIndex);
			index.set(maxIndex, i);
			index.set(i, temp);
			
			for(int j = i + 1; j < numOfVars; j++) {
				float mult = coeff[index.get(j)][i] / coeff[index.get(i)][i];
				
				for(int k = i + 1; k < numOfVars; k++) {
					coeff[index.get(j)][k] = coeff[index.get(j)][k] - (mult * coeff[index.get(i)][k]);
				}
				float result = constants.get(index.get(j)) - (mult* constants.get(index.get(i)));
				constants.set(index.get(j), result); 
			}
		}
	}
	
	private static void SPPBackSubstitution(float[][] coeff, Vector<Float> constants, Vector<Float> sol,
											Vector<Integer> index) {

		float result = constants.get(numOfVars-1) / coeff[index.get(numOfVars-1)][numOfVars-1];
		sol.set(numOfVars-1, result);
		
		for(int i = numOfVars-1; i >= 0; i--) {
			float sum = constants.get(index.get(i));
			for(int j = i + 1; j < numOfVars; j++) {
				sum = sum - coeff[index.get(i)][j] * sol.get(j);
			}
			sol.set(i, (sum/ coeff[index.get(i)][i]));
		}
	}
}
