/* /////////////////////////////////////////////////////////////////////
//
// CLASS       	: BackPropagationTester
// AUTHOR(S)   	: Houle, Daniel B
// DESCRIPTION 	: Reads formated information from the Network.printNetworkInfo()
//							method. It will take the weight & threshold information
//							& update the them accordingly. 
// DEPENDENCIES   : Network must be in the same directory.
//							The information file must be in the same directory.
//
///////////////////////////////////////////////////////////////////// */
import java.io.Console;
import java.util.Arrays;
import java.io.IOException;
import java.util.*;
import java.util.ArrayList.*;

public class BackPropagationTester{

public static void main(String[] args) throws IOException{
	
	Console c = System.console();
  	if (c == null) {
		System.err.println("No console.");
      System.exit(1);
	}
	Stats.initializeVariables();
	String inputNum = c.readLine("Enter # of inputs: ");
	int in = (new Integer(inputNum.trim())).intValue();
	Stats.setInputNumber(in);
	String outputNum = c.readLine("Enter # of outputs: ");
	int out = (new Integer(outputNum.trim())).intValue();
	Stats.setOutputNumber(out);
	String numLayers = c.readLine("Enter # of layers: ");
	int num = (new Integer(numLayers.trim())).intValue();
	Stats.setNumberOfHiddenLayers(num - 1);
	if(num > 1){
		String hiddenLayers = c.readLine("Enter # of Neurons per layer, seperated by a comma(,): ");
		String[] stringLayers = hiddenLayers.split(",");
		int[] intLayers = new int[stringLayers.length];
		for(int i = 0; i < stringLayers.length; i++){ intLayers[i] = (new Integer(stringLayers[i].trim())).intValue();}
		Stats.setNumPerLayer(intLayers);
	}
	String activation = c.readLine("Enter activation function: ");
	Stats.setActivationFunction(activation.trim());
	String weights = c.readLine("Enter weights, seperated by a colon(:): ");
	String[] indWeights = weights.split(":");
	double[] eachWeights = new double[indWeights.length];
	//go through and convert every weight into a double
	for(int i = 0; i < indWeights.length; i++){ eachWeights[i] = (new Double(indWeights[i].trim())).doubleValue();}
	String thresholds = c.readLine("Enter thresholds, seperated by a colon(:): ");
	String[] indThetas = thresholds.split(":");
	double[] eachTheta = new double[indThetas.length];
	for(int i = 0; i < indThetas.length; i++){ eachTheta[i] = (new Double(indThetas[i].trim())).doubleValue();}
	Network ann = new Network();
	ann.updateNetworkWeights(eachWeights);
	ann.updateNetworkThresholds(eachTheta);
	System.out.println("Exit by entering \"Exit Now\"");
	do{
		String pattern = c.readLine("Enter a pattern, seperated by a comma(,): ");
		if(pattern.equalsIgnoreCase("exit now")){ break;}
		String[] inputPattern = pattern.split(",");
		if(inputPattern.length != in){ System.out.println("Incorrect number of elements in the pattern. Try again"); continue;}
		double[] doublePattern = new double[in];
		for(int i = 0; i < in; i++){ doublePattern[i] = (new Double(inputPattern[i].trim())).doubleValue();}
		double[] output = ann.think(doublePattern);
		if(output.length != out){ System.out.println("Something went wrong."); break;}
		System.out.print("\nOutput(s): ");
		for(int i = 0; i < out; i++){ System.out.print(output[i] + "	");}
		System.out.println("");
	}while(true);
	System.exit(0);
}


//end class BackPropagationTester
}