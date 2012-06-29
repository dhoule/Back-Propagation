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
	try(BufferedReader infoInput = Files.newBufferedReader(infoPath, charset)){		
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
	catch(IOException e){
		System.err.println("\n\nIOException: " + e.toString() + ".\n\n");
	}
	//try to read the data file
	try(BufferedReader dataInput = Files.newBufferedReader(dataPath, charset)){		
		String line = null;
		while ((line = dataInput.readLine()) != null) {
			System.out.println(line);
		}
	}
	catch(IOException e){
		System.err.println("\n\nIOException: " + e.toString() + "\n\n");
	}
}

/**
* Sets up the information in the static Stats class that is used to build and run the ANN.
* @param property An array of 2 elements that simulates a name/value pairing. It must match
*						certain criteria or the program will exit with an error code of 2.
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
	else if(name.equalsIgnoreCase("Neurons Per Hidden Layer") || name.equalsIgnoreCase("Neurons Per Hidden Layer") || 
			  name.equalsIgnoreCase("Neurons Per Hidden Layer") || name.equalsIgnoreCase("Neurons Per Hidden Layer") ||
			  name.equalsIgnoreCase("Neurons Per Hidden Layer") || name.equalsIgnoreCase("NeuronsPerHiddenLayer")){
		int[] layers = null;
		if(value.contains(",")){
			String[] num = value.split(",");
			for(int i = 0; i < num.length; i++){ layers[i] = (new Integer(num[i])).intValue();}
		}
		else{
			layers = new int[1];
			layers[0] = new Integer(value).intValue();
		}
		Stats.setNumPerHiddenLayer(layers);
	}
	else{
		System.err.println("\n\nThe name property " + name + " is unknown.\n\n");
		System.exit(2);
	}
}

//end class BackPropagationDriver
}