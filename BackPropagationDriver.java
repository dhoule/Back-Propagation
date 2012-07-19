/* /////////////////////////////////////////////////////////////////////
//
// CLASS       	: BackPropagationDriver
// AUTHOR(S)   	: Houle, Daniel B
// DESCRIPTION 	: Reads the information files and runs the ANN for the 
//							number of runs or the ANN has reached the allowed
//							error rate. The Network is built depending on the 
//							name/value pairs; seperated by colons; in the 
//							information; "info" for short; file. The patterns
//							the ANN is to learn & the corresponding output(s) are
//							to be in the data file. The symbols of each pattern
//							& it's corresponding output(s) are to be seperated
//							by commas(,). White space is not required but is
//							allowed in order for the data to be more human 
//							readable.
// DEPENDENCIES   : Network & Stats must be in the same directory.
//							The information file must be in the same directory.
//							The data file must be in the same directory.
//
///////////////////////////////////////////////////////////////////// */
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.Files;
import java.util.*;
import java.util.ArrayList.*;

public class BackPropagationDriver{

/**
* Receives the information required to create and train the ANN.
* @param args Command line arguments. 1) Information file 2) Data file
**/
public static void main(String[] args){

	if(args.length != 2){
		System.err.println("\n\nusage: java BackPropagationDriver <info_file> <data_file>\n\n");
		System.exit(1);
	}

	Path fileOne = Paths.get(args[0]);
	Path fileTwo = Paths.get(args[1]);
	Path infoPath = fileOne.toAbsolutePath();
	Path dataPath = fileTwo.toAbsolutePath();
	if(!Files.isRegularFile(infoPath) || !Files.isReadable(infoPath)){
		System.err.println("\n\nFile " + infoPath + " is not readable or doesn't exist.\n\n");
		System.exit(1);
	}
	if(!Files.isRegularFile(dataPath) || !Files.isReadable(dataPath)){
		System.err.println("\n\nFile " + dataPath + " is not readable or doesn't exist.\n\n");
		System.exit(1);
	}
	//The UTF-8 charset is specified by RFC 2279; the transformation format upon which it is based is 
	//specified in Amendment 2 of ISO 10646-1 and is also described in the Unicode Standard.
	Charset charset = Charset.forName("UTF-8");
	
	//try to read the information file
	BufferedReader infoInput = null;
	try{
		infoInput = Files.newBufferedReader(infoPath, charset);		
		String line = null;
		String[] nameValue = null;
		Stats.initializeVariables();
		while ((line = infoInput.readLine()) != null) {
			nameValue = line.split(":");
			if(nameValue.length != 2){
				System.err.println("\n\nThe line \"" + line + "\", in the info_file, is not in the correct format.");
				System.err.println("Correct format is:");
				System.err.println("name: value\n\n");
				System.exit(1);
			}
			setupStatsInfo(nameValue[0].trim(), nameValue[1].trim());
		}
	}
	catch(IOException e){ System.err.println("\n\nIOException: " + e.toString() + ".\n\n");}
	
	//done with file. Time to close it.
	try{ infoInput.close();}
	catch(IOException e){ System.err.println("\n\nIOException: " + e.toString() + ".\n\n");}
	
	
	ArrayList<double[]> patterns = new ArrayList<double[]>();
	//try to read the data file
	BufferedReader dataInput = null;
	try{
		dataInput = Files.newBufferedReader(dataPath, charset);		
		String line = null;
		int totalElementsCount = Stats.getInputNumber() + Stats.getOutputNumber();
		String[] patternOutput = null;
		while ((line = dataInput.readLine()) != null) {
			//pattern must be seperated by commas
			if(line.contains(",")){
				patternOutput = line.split(",");
				//need to make sure the pattern and it's output(s) == the input and output #
				if(patternOutput.length != totalElementsCount){
					try{ dataInput.close();}
					catch(IOException e){ System.err.println("\n\nIOException: " + e.toString() + ".\n\n");}
					System.err.println("Pattern elements does not equal the # of inputs and outputs given.");
					System.exit(1);
				}
				//convert the elements into doubles
				double[] elements = new double[totalElementsCount];
				for(int i = 0; i < totalElementsCount; i++){ elements[i] = (new Double(patternOutput[i])).doubleValue();}
				//add this pattern to patterns
				patterns.add(elements);
			}
			else{
				try{ dataInput.close();}
				catch(IOException e){ System.err.println("\n\nIOException: " + e.toString() + ".\n\n");}
				System.err.println("Pattern of data file is incorrect. Pattern must be seperated by commas.");
				System.exit(1);
			}
		}
	}
	catch(IOException e){ System.err.println("\n\nIOException: " + e.toString() + "\n\n");}
	
	//done with file. Time to close it.
	try{ dataInput.close();}
	catch(IOException e){ System.err.println("\n\nIOException: " + e.toString() + ".\n\n");}
	
	//The information has been collected, formated, and setup.
	receiving(patterns);
	
	System.exit(0);
}

/**
* Sets up the information in the static Stats class that is used to build and run the ANN.
* @param name The name property
* @param value The value property
**/
private static void setupStatsInfo(String name, String value){
	if(name.equalsIgnoreCase("Input")){
		Stats.setInputNumber((new Integer(value)).intValue());
	}
	else if(name.equalsIgnoreCase("Output")){
		Stats.setOutputNumber((new Integer(value)).intValue());
	}
	else if(name.equalsIgnoreCase("Hidden Layers") || name.equalsIgnoreCase("HiddenLayers")){
		Stats.setNumberOfHiddenLayers((new Integer(value)).intValue());
	}
	else if(name.equalsIgnoreCase("Activation Function") || name.equalsIgnoreCase("ActivationFunction")){
		Stats.setActivationFunction(value);
	}
	else if(name.equalsIgnoreCase("Alpha")){
		Stats.setLearningRate((new Double(value)).doubleValue());
	}
	else if(name.equalsIgnoreCase("Error Rate") || name.equalsIgnoreCase("ErrorRate")){
		Stats.setAllowedErrorRate((new Double(value)).doubleValue());
	}
	else if(name.equalsIgnoreCase("Neurons Per Layer") || name.equalsIgnoreCase("Neurons PerLayer") || 
			  name.equalsIgnoreCase("NeuronsPer Layer") || name.equalsIgnoreCase("NeuronsPerLayer")){
		int[] layers = null;
		if(value.contains(",")){
			String[] num = value.split(",");
			layers = new int[num.length];
			for(int i = 0; i < num.length; i++){ layers[i] = (new Integer(num[i].trim())).intValue();}
		}
		else{
			layers = new int[1];
			layers[0] = new Integer(value.trim()).intValue();
		}
		Stats.setNumPerLayer(layers);
	}
	else{
		System.err.println("\n\nThe name property " + name + " is unknown.\n\n");
		System.exit(2);
	}
}

private static void receiving(ArrayList<double[]> patterns){
	
	//create a new ANN
	Network ann = new Network();
	//get the allowed error rate
	double errorRate = Stats.getAllowedErrorRate();
	//get the # of inputs
	int numInput = Stats.getInputNumber();
	//get the # of outputs
	int numOutput = Stats.getOutputNumber();
	//determine the # of patterns
	int numPatterns = patterns.size();
	//this will keep track of the # of times the patterns are run through
	int epochCount = 0;
	//need to go through at least ONE run through the pattern
	do{
		epochCount++;
	}while(run(patterns, numPatterns, ann, errorRate, numInput, numOutput));
	
	System.out.println("\n\nIt took " + epochCount + " runs.\n");
	ann.printNetworkInfo();
	System.out.println("\n\n");
}

private static boolean run(ArrayList<double[]> patterns, int numPatterns, 
								  Network ann, double errorRate, int numInput, int numOutput){
								  
	//double[] runError = new double[numPatterns];
	boolean errorFound = false;
	//go through the patterns
	for(int i = 0; i < numPatterns; i++){
		//need to get the inputs & outs from the pattern
		double[] inputs = new double[numInput];
		double[] desiredOutput = new double[numOutput];
		double[] currentPattern = patterns.get(i);
		//seperate the inputs from the pattern
		for(int q = 0; q < numInput; q++){ inputs[q] = currentPattern[q];}
		//seperate the outputs from the pattern
		for(int q = 0; q < numOutput; q++){ desiredOutput[q] = currentPattern[q + numInput];}
		
		//double currentError = 0.0;
		double[] actualOutput = ann.think(inputs);
		double[] individualErrors = new double[numOutput];
		//need to compare the actualOutput against the desiredOutput & use the sum of the squares
		for(int q = 0; q < numOutput; q++){ 
			individualErrors[q] = desiredOutput[q] - actualOutput[q];
			if(Math.abs(individualErrors[q]) > errorRate){ errorFound = true;}
			//currentError += individualErrors[q] * individualErrors[q];
		}
		//runError[i] = currentError;
		if(errorFound){ ann.backPropagateError(individualErrors, actualOutput);}
	}
	
	return errorFound;
}
//end class BackPropagationDriver
}