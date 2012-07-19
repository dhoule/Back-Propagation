/* /////////////////////////////////////////////////////////////////////
//
// CLASS       	: Network
// AUTHOR(S)   	: Houle, Daniel B
// DESCRIPTION 	: Creates the ANN. A Network object is composed of 1/+ 
//							NeuronLayer objects
// DEPENDENCIES   : NeuronLayer & Stats must be in the same directory
//
///////////////////////////////////////////////////////////////////// */
import java.util.*;
import java.util.ArrayList.*;


public class Network{
/** store the NeuronLayers of the Network **/
private ArrayList<NeuronLayer> layers;
/** number of inputs into this Network **/
private int inputs;
/** number of outputs from this Network **/
private int outputs;
/** number of hidden layers in this Network **/
private int hidden_layers;
/** array of the number of Neurons in the hidden layers from top to down **/
private int[] neuronsInHiddenLayers;
/** the last outputs of this Network **/
private double[] output;
/** the total number of weights in this Network **/
private int num_weights;
/** activation function that is to be used for this object **/
private String current_activation_function;
/** the learning rate for every Neuron in this Network **/
private double alpha;

/**
* Constructor sets up the layers from the information it gets from Stats.
**/
public Network(){
	//get the number of inputs for this Network
	inputs = Stats.getInputNumber();
	//get the number of outputs for this Network
	outputs = Stats.getOutputNumber();
	//get the number of hidden layers for this Network
	hidden_layers = Stats.getNumberOfHiddenLayers();
	//get the size of each hidden layer in this Network
	neuronsInHiddenLayers = Stats.getNumPerLayer();
	//get the activation function that will be used for this Network
	current_activation_function = Stats.getActivationFunction();
	//get the learning rate for this Network
	alpha = Stats.getLearningRate();
	//initialize the ArrayList
	layers = new ArrayList<NeuronLayer>();
	
	//if there are actually hidden layers in the network
	//create a fully connected ANN
	if(hidden_layers > 0){
		//input layer doesn't really matter. the 1st hidden layer, or the "last" layer is what actually does the math.
		//the # of inputs into the network, will be the number of inputs into the neurons of the 1st "thinking" layer			
		//used to setup the inputs of a layer based on the number of neurons in the prior layer
		int newInput = 0;			
		//add the hidden layer(s) plus the output Layer. 'hidden_layers' has 1 added to it just incase there are no hidden layers
		for(int i = 0; i < (hidden_layers + 1); i++){
			//the # of inputs of the current layer, is based upon the # of neurons in the last layer
			if(i == 0){ newInput = inputs;}
			else{ newInput = neuronsInHiddenLayers[i - 1];}
			//need to create the new NeuronLayer objects and add them to layers
			if(i < hidden_layers){ layers.add(new NeuronLayer(neuronsInHiddenLayers[i], newInput, current_activation_function, alpha));}
			else{ layers.add(new NeuronLayer(outputs, newInput, current_activation_function, alpha));}
		}			
	}
	else{
		//there are no hidden layers.
		layers.add(new NeuronLayer(outputs, inputs, current_activation_function, alpha));
	}
	
	//accumulator
	num_weights = 0;
	//traverse each layer
	for(int i = 0; i < (hidden_layers + 1); i++){
		num_weights += layers.get(i).getNumWeights();
	}	
//end constructor
}

/**
* Returns the number of layers in this Network
* @return int
**/
public int getLayerNum(){ return (hidden_layers + 1);}

/**
* Returns a single Layer of this Network
* @return NeuronLayer
**/
public NeuronLayer getIndividualLayer(int i){ return layers.get(i);}

/**
* Returns the weights of the Network in order by layer & left to right
* @return double[]
**/
public double[] getNetworkWeights(){		
	double[] r_weights = new double[0];			
	for(int i = 0; i < (hidden_layers + 1); i++){ r_weights = extendArray(r_weights, layers.get(i).getLayerWeights());}		
	return r_weights;
//end getWeights()
}

/**
* Returns the thresholds of the Network in order by layer & left to right
* @return double[]
**/
public double[] getNetworkThresholds(){
	double[] r_thresholds = new double[0];			
	for(int i = 0; i < (hidden_layers + 1); i++){
		r_thresholds = extendArray(r_thresholds, layers.get(i).getLayerThresholds());
	}		
	return r_thresholds;
}

/**
* Used to extend a primitve array
* @return double[]
**/
private double[] extendArray(double[] original, double[] dirty){
	int size = original.length + dirty.length;
	int dirtyCt = 0;
	double[] eArray = new double[size];		
	for(int i = 0; i < size; i++){
		if(i < original.length){ eArray[i] = original[i];}
		else{
			eArray[i] = dirty[dirtyCt];
			dirtyCt++;
		}
	}		
	return eArray;
}

/**
* Updates the weights of this Network by layer & left to right
* @param new_weights The new weights for the Neurons of this Network
**/
public void updateNetworkWeights(double[] new_weights){
	//used to determine the number of weights in each NeuronLayer
	int layerWeights = 0;
	//this is used to itterate throuh new_weights
	int current = 0;		
	//traverse the layers & find out how many weights/Neuron there are
	for(int i = 0; i < (hidden_layers + 1); i++){
		layerWeights = layers.get(i).getNumWeights();
		double[] iWeights = new double[layerWeights];
		//get the weights needed for the current layer
		for(int q = 0; q < layerWeights; q++){
			iWeights[q] = new_weights[current];
			current++;
		}			
		layers.get(i).setLayerWeights(iWeights);
	}		
//end setWeights()
}
	
/**
* Updates the thresholds of this Network by layer & left to right
* @param new_thetas The new thresholds for the Neurons of this Network
**/
public void updateNetworkThresholds(double[] new_thetas){
	//used to determine the # of Perceptons in each layer
	int layerNum = 0;
	//this is used to itterate throuh new_thetas
	int current = 0;
	for(int i = 0; i < (hidden_layers + 1); i++){
		layerNum = layers.get(i).getNumLayer();
		double[] iTheta = new double[layerNum];
		//get the thresholds needed for the current layer
		for(int q = 0; q < layerNum; q++){
			iTheta[q] = new_thetas[current];
			current++;
		}
		layers.get(i).setLayerThresholds(iTheta);
	}
}

/**
* Returns the number of weights in this Network
* @return int
**/
public int getTotalWeightNum(){ return num_weights;}

/**
* Updates the Network with new inputs
* @param temp_inputs new inputs for this Network
* @return double[]
**/
public double[] think(double[] temp_inputs){	
	if(temp_inputs.length != inputs){
		System.err.println("New inputs doesn't equal original number");
		System.exit(1);
	}
	//this is used to capture the output of each layer
	output = new double[outputs];	
			
	//go through the layers of this Network
	for(int i = 0; i < (hidden_layers + 1); i++){
		output = layers.get(i).think(temp_inputs);
		//reassign temp_inputs to be used again			
		temp_inputs = new double[output.length];	
		for(int q = 0; q < temp_inputs.length; q++){ temp_inputs[q] = output[q];}		
	}				
	//returns the output of the last layer
	return output;
//end update()
}

/**
* Returns an array containing the number of neurons in each layer. 
* Each element contains that layer's number of Neurons.
* @return int[]
**/
public int[] getNumNeuronsPerLayer(){
	int[] ret = new int[hidden_layers + 1];
	//traverse each layer
	for(int i = 0; i < (hidden_layers + 1); i++){ ret[i] = layers.get(i).getNumLayer();}
	return ret;
}

/**
* Return the number of inputs of the current Network.
* @return int
**/
public int getNetworkInputs(){ return inputs;}
	
/**
* Takes the error and moves it back through the Network
**/
public void backPropagateError(double[] error, double[] output){
	
	int[] layers = this.getNumNeuronsPerLayer();
	int numInputs = this.getNetworkInputs();
	int numLayer = layers.length;
	int error_num = error.length;
	double[] error_prop, current_errorProp, current_weights;
	//ct is a general counter
	int ct = 0, currentNum = -1, previousNum = -1;
	NeuronLayer current = new NeuronLayer(), previous = new NeuronLayer();
	//calculate the error gradiant for the output layer
	if(current_activation_function.equalsIgnoreCase("sigmoid")){
		//delta_x = Y_x * (1 - Y_x) * e_x
		//e_x = error of Neuron X
		//Y_x = output og Neuron X
		//X = Neuron in Output Layer
		//output[i] * (1.0 - output[i]) is the derivative of the sigmoid function
		for(int i = 0; i < error_num; i++){ error[i] *= output[i] * (1.0 - output[i]);}
	}	
	//need to update the error gradient of the last layer
	this.layers.get(numLayer - 1).setLayerPropagationError(error);
	//need to get the last inputs from the output layer
	double[] lastLayerInputs = this.layers.get(numLayer - 1).getLayerInput();		
	//need to get the number of weights in the output layer
	int outputWeightNum = this.layers.get(numLayer - 1).getNumWeights();
	//used to hold the delta_weights of the output layer
	double[] deltaWeights = new double[outputWeightNum];
	//calculate the delta weights for the output layer
	//go through each Neuron in the last layer
	for(int i = 0; i < layers[layers.length - 1]; i++){
		//go through each input of the last layer
		for(int j = 0; j < lastLayerInputs.length; j++){
			deltaWeights[ct] = alpha * lastLayerInputs[j] * error[i];
			ct++;
		}
	}
	//need to update the delta weights of the last layer
	this.layers.get(numLayer - 1).setLayerDeltaWeights(deltaWeights);
	//calculate the delta thetas for the output layer
	double[] deltaThetas = new double[layers[layers.length - 1]];
	//go through each Neuron in the last layer
	for(int i = 0; i < layers[layers.length - 1]; i++){ deltaThetas[i] = alpha * (-1.0) * error[i];}
	//need to update the delta thetas of the last layer
	this.layers.get(numLayer - 1).setLayerDeltaThetas(deltaThetas);
	
	
	//go through the other layers and do the same things as above
	for(int i = numLayer - 2; i >= 0; i--){
		//need to get the "current" layer and "previous" layer			
		current = this.getIndividualLayer(i);
		currentNum = this.layers.get(i).getNumLayer();
		previous = this.getIndividualLayer(i + 1);
		previousNum = this.layers.get(i + 1).getNumLayer();
		//the output of the current layer is the input of the previous layer
		double[] currentOuput = previous.getLayerInput();
		double[] currentInput = this.layers.get(i).getLayerInput();
		deltaWeights = new double[this.layers.get(i).getNumWeights()];
		double[] previousErrorGradient = previous.getLayerPropagationError();
		double[] previousWeights = previous.getLayerWeights();
		//need to calculate the error gradients of the current layer
		error = new double[currentNum];
		//reset the general counter
		ct = 0;
		//calculate the error gradiant for the current layer
		if(current_activation_function.equalsIgnoreCase("sigmoid")){
			//go through the current layer
			for(int q = 0; q < currentNum; q++){
				double derivative = currentOuput[q] * (1.0 - currentOuput[q]);
				//go through the previous layer
				for(int w = 0; w < previousNum; w++){
					error[q] += derivative * previousErrorGradient[w] * previousWeights[ct];
					ct++;
				}
			}
		}
		//need to update the error gradient of the current layer
		this.layers.get(i).setLayerPropagationError(error);
		//reset the general counter
		ct = 0; 	
		//calculate the delta weights for the current layer
		//go through each Neuron in the current layer
		for(int q = 0; q < currentNum; q++){
			//go through each input of the previous layer
			for(int j = 0; j < currentInput.length; j++){					
				deltaWeights[ct] = alpha * currentInput[j] * error[q];
				//this is to check for -0.0
				if(Math.abs(deltaWeights[ct]) == 0.0){ deltaWeights[ct] = 0.0;}
				ct++;
			}
		}
		//need to update the delta weights of the last layer
		this.layers.get(i).setLayerDeltaWeights(deltaWeights);
		
		//calculate the delta thetas for the current layer
		deltaThetas = new double[currentNum];
		//go through each Neuron in the current layer
		for(int q = 0; q < currentNum; q++){ deltaThetas[q] = alpha * (-1.0) * error[q];}
		//need to update the delta thetas of the current layer
		this.layers.get(i).setLayerDeltaThetas(deltaThetas);			
	}
	propagateWeightsAndThetas();
}
	
/**
* Goes through each layer and forces them to update their weights
**/
private void propagateWeightsAndThetas(){
	//traverse the layers & forces them to update their weights & thresholds
	for(int i = 0; i < (hidden_layers + 1); i++){ layers.get(i).propagateWeightsAndThetas();}
}

/**
* Prints the weights of this Network to standard output.
**/
public void printWeights(){
	System.out.println("Weights:");
	double[] temp_weights = this.getNetworkWeights();
	for(int i = 0; i < temp_weights.length; i++){
		System.out.print(i + ": " + temp_weights[i] + "\t");
	}
	System.out.println("\n");
}
	
/**
* Prints the thresholds of this Network to standard output.
**/
public void printThresholds(){
	System.out.println("Thresholds:");
	double[] temp_thresholds = this.getNetworkThresholds();
	for(int i = 0; i < temp_thresholds.length; i++){
		System.out.print(i + ": " + temp_thresholds[i] + "\t");
	}
	System.out.println("\n");
}

/**
* Prints the # of inputs, # of outputs, individual weights and thresholds,
* the learning rate, # of layers, and the # of Neurons in each layer for this Network.
**/
public void printNetworkInfo(){
	System.out.println("Network information");
	System.out.println("Activation Function: " + current_activation_function);
	System.out.println("# of inputs: " + inputs);
	System.out.println("# of outputs: " + outputs);
	this.printLayerCountInfo();
	this.printWeights();
	this.printThresholds();
	System.out.println("___________________");
}
	
/**
* Prints the LATEST output of this Network.
**/
public void printOutput(int ct){
	System.out.println("*** Network Outputs for " + ct + " ***");
	for(int i = 0; i < outputs; i++){
		System.out.print(i + ": " + output[i]);
	}
	System.out.println("");
}

/**
* Prints the # of Neurons in every Layer in this Network.
**/
public void printLayerCountInfo(){
	int layerNum = layers.size();
	System.out.println("Layer information:");
	for(int i = 0; i < layerNum; i++){System.out.println("Layer " + i + ": " + layers.get(i).getNumLayer());}
}

/* /////////////////////////////////////////////////////////////////////
//
// CLASS       	: NeuronLayer
// AUTHOR(S)   	: Houle, Daniel B
// DESCRIPTION 	: Holds mulptiple Neurons for the Network.
// DEPENDENCIES   : Neuron must be in the same directory
//
///////////////////////////////////////////////////////////////////// */
public class NeuronLayer{

/** Holds the individual Neurons for this layer **/
private ArrayList<Neuron> neurons;
/** The # of Neurons in this layer **/
private int num_neurons;
/** The # of weights into this layer **/
private int num_weights;

/** This is here to make the compiler happy. **/
public NeuronLayer(){}

/** 
* Constructor must be told the number of Neurons that are in this layer, 
* how many inputs each Neuron takes, and the name of the activation 
* function that each Neuron will use. 
* @param num # of Neurons in this Layer
* @param inputs # of inputs into this Layer. Every Neuron recieves the same input.
* @param func The activation function for each Neuron in this Layer. 
* @param alpha The learning rate for each Neuron in this Layer.
**/
public NeuronLayer(int num, int inputs, String func, double alpha){
	neurons = new ArrayList<Neuron>();
	//set the # of Neurons in this Layer
	num_neurons = num;
	//create the Neurons for this Layer
	for(int i = 0; i < num; i++){ neurons.add(new Neuron(inputs, (new Random()).nextLong(), func, alpha));}
	//set the # of weights for this Layer
	num_weights = num * inputs;
}

/**
* Returns the # of Neurons in this Layer.
* @return	The # of Neurons in this Layer.
**/
public int getNumLayer(){ return num_neurons;}

/**
* This is pretty much here only for standards purposes.
* @param nothing 
**/
public void setNumLayer(int nothing){ num_neurons = nothing;}

/**
* Returns the weights of this layer for each Neuron from left to right.
* @return	This Layer's weights.
**/
public double[] getLayerWeights(){
	double[] weights = new double[num_weights];
	int ct = 0;
	for(int i = 0; i < num_neurons; i++){
		double[] current = neurons.get(i).getWeights();
		for(int q = 0; q < current.length; q++){
			weights[ct] = current[q];
			ct++;
		}
	}
	return weights;
}

/**
* Updates the weights of this layer for each Neuron from left to right.
* @param weights The new weights for this Layer
**/
public void setLayerWeights(double[] weights){
	int weightsPerNeuron = num_weights / num_neurons;
	int ct = 0;
	//go through each Neuron
	for(int i = 0; i < num_neurons; i++){
		double[] individualWeights = new double[weightsPerNeuron];
		//gather the individual weights
		for(int q = 0; q < weightsPerNeuron; q++){ individualWeights[q] = weights[ct]; ct++;}
		neurons.get(i).setWeights(individualWeights);
	}
}

/**
* Updates the propagation error for every Neuron in thist Layer.
* @param error The propagation errors for this Layer's Neurons.
**/
public void setLayerPropagationError(double[] error){ for(int i = 0; i < num_neurons; i++){ neurons.get(i).setPropagateError(error[i]);}} 

/**
* Returns an array consisting of the propagation errors for every Neuron in this Layer.
* @return	An array of the propagation erros for this Layer.
**/
public double[] getLayerPropagationError(){
	double[] errors = new double[num_neurons];
	for(int i = 0; i < num_neurons; i++){ errors[i] = neurons.get(i).getPropagateError();}
	return errors;
}

/**
* Updates the thresholds of each Neuron in this Layer.
* @param thetas The new thresholds for each Neuron.
**/
public void setLayerThresholds(double[] thetas){ for(int i = 0; i < num_neurons; i++){ neurons.get(i).setTheta(thetas[i]);}}

/**
* Returns the thresholds for each Neuron in this layer.
* @return	An array of the thresholds for this Layer
**/
public double[] getLayerThresholds(){
	double[] thetas = new double[num_neurons];
	for(int i = 0; i < num_neurons; i++){ thetas[i] = neurons.get(i).getTheta();}
	return thetas;
}

/**
* Returns the # of weights in this Layer.
* @return	# of weights.
**/
public int getNumWeights(){ return num_weights;}

/**
* Returns this Layer as a primitive array.
* @return	A primitive array of the Neurons in this Layer.
**/
public Neuron[] getLayer(){ 
	Neuron[] primitive = new Neuron[num_neurons];
	primitive = neurons.toArray(primitive);
	return primitive; 
}

/**
* Receives the inputs & returns the outputs of this Layer
* @param in The inputs for this Layer
* @return	The output of every Neuron in this Layer from left to right
**/
public double[] think(double[] in){
	//Each Neuron receives multiple inputs but will only give ONE output
	double[] outputs = new double[num_neurons];
	//Each Neuron in this layer will receive the same inputs
	for(int i = 0; i < num_neurons; i++){ outputs[i] = neurons.get(i).think(in);}
	return outputs;
}

/**
* Returns this Layer's last inputs. Each Neuron in a Layer receives the same input, so only have to return the last input
* of the Neuron at position 0.
* @return	The last inputs of this Layer.
**/
public double[] getLayerInput(){ return neurons.get(0).getInput();}

/**
* Updates the delta weights of this Layer for Back Propagation purposes.
* @param delta The offset for each weight in this Layer.
**/
public void setLayerDeltaWeights(double[] delta){
	int weightsPerNeuron = num_weights / num_neurons;
	int ct = 0;
	//go through each Neuron
	for(int i = 0; i < num_neurons; i++){
		double[] deltaWeights = new double[weightsPerNeuron];
		//gather the individual weights
		for(int q = 0; q < weightsPerNeuron; q++){ deltaWeights[q] = delta[ct]; ct++;}
		neurons.get(i).setDeltaWeights(deltaWeights);
	}
}

/**
* Returns an array of the delta weights for this Layer for Back Propagation.
* @return	An array of the delta weights for this Layer.
**/
public double[] getLayerDeltaWeights(){
	double[] delta = new double[num_weights];
	int ct = 0;
	for(int i = 0; i < num_neurons; i++){
		double[] current = neurons.get(i).getDeltaWeights();
		for(int q = 0; q < current.length; q++){
			delta[ct] = current[q];
			ct++;
		}
	}
	return delta;
}

/**
* Updates the delta thresholds of this Layer for Back Propagation purposes.
* @param delta The offset for each threshold in this Layer.
**/
public void setLayerDeltaThetas(double[] delta){
	if(delta.length != num_neurons){
		System.err.println("\n\nThe # of delta thresholds sent to this Layer is incorrect.\n\n");
		System.exit(1);
	}
	//go through each Neuron
	for(int i = 0; i < num_neurons; i++){ neurons.get(i).setDeltaTheta(delta[i]);}
}

/**
* Returns an array of the delta thresholds for this Layer for Back Propagation.
* @return	An array of the delta thresholds for this Layer.
**/
public double[] getLayerDeltaThetas(){
	double[] delta = new double[num_neurons];
	for(int i = 0; i < num_neurons; i++){ delta[i] = neurons.get(i).getDeltaTheta();}
	return delta;
}

/**
* Updates the weights and thresholds for every neuron in this Layer.
**/
public void propagateWeightsAndThetas(){
	for(int i = 0; i < num_neurons; i++){ neurons.get(i).propagateWeightsAndThetas();}
}


/* /////////////////////////////////////////////////////////////////////
//
// CLASS       	: Neuron
// AUTHOR(S)   	: Houle, Daniel B
// DESCRIPTION 	: Creates an individual Neuron for the larger Network.
// DEPENDENCIES   : Mthods must be called correctly.
//
///////////////////////////////////////////////////////////////////// */
private class Neuron{

/** allowed activation functions: Step, Sign, Sigmoid, Linear, and Hyperbolic Tangent **/
private String[] activationFunction = {"step", "sign", "sigmoid", "linear", "hyperbolic tangent"};
/** activation function that is to be used for this Neuron **/
private String currentActivationFunction;
/** number of inputs into this Neuron **/
private int num_inputs;
/** the last inputs for the Neuron **/
private double[] input;
/** the last output of this Neuron **/
private double output;
/** the learning rate for this Neuron **/
private double alpha;
/** the threshold for this Neuron **/
private double theta;
/** the weights for each input for this Neuron **/
private double[] weight;
/** the propagation error for this Neuron **/
private double propagate_error;
/** the delta weights that will be used for Back Proagation **/
private double[] deltaWeights;
/** the delta threshold that will be used for Back Proagation **/
private double deltaTheta;

/** 
* Constructor must be told how many inputs it will be recieving, a seed to randomize the initial 
* weights & threshold, the name of the activation function it will be using, and the learning rate. 
* @param num # of inputs into this Neuron.
* @param seed Used to seed the Random object used for randomizing the weights and threshold of this Neuron.
* @param func The activation function for this Neuron.
* @param learningRate The learning rate for this Neuron.
**/
public Neuron(int num, long seed, String func, double learningRate){
	boolean found = false;
	//need to check if the activation function is valid
	for(int i = 0; i < activationFunction.length; i++){
		if(func.equalsIgnoreCase(activationFunction[i])){ 
			currentActivationFunction = func;
			found = true;
			break;
		}
	}
	if(!found){
		System.err.println("\n\nActivation function not recognized.\n\n");
		System.exit(1);
	}
	alpha = learningRate;
	//set the size of the inputs
	num_inputs = num;
	//set the new threshold
	Random offset = new Random(seed);
	theta = offset.nextDouble();
	//set the new weights
	weight = new double[num];
	deltaWeights = new double[num];
	//uniformly distributed double value between -0.5 and 0.5
	for(int i = 0; i < num; i++){ weight[i] = offset.nextDouble() - 0.5;}		
}

//Standard Setters and Getters

/**
* Updates the threshold for this Neuron.
* @param threshold New threshold 
**/
public void setTheta(double threshold){ theta = threshold;}

/**
* Returns the threshold for this Neuron.
* @return 	The current threshold for this Neuron 
**/
public double getTheta(){ return theta;}

/**
* Updates the weights of this Neuron.
* $param newWeights New weights
**/
public void setWeights(double[] newWeights){ weight = newWeights;}

/**
* Returns the current weights of this Neuron.
* @return	The current weights.
**/
public double[] getWeights(){ return weight;}

/**
* Returns the # of weights for this Neuron.
* @return	# of weights
**/
public int getNumWeights(){ return num_inputs;}

/**
* Returns the # of inputs for this Neuron.
* @return	# of inputs
**/
public int getNumInputs(){ return num_inputs;}

/**
* Returns the last output of this Neuron.
* @return	Last computed output
**/
public double getOutput(){ return output;}

/**
* Returns the last input into this Neuron.
* @return	An array of the last inputs this Neuron recieved.
**/
public double[] getInput(){ return input;}

//End standard Setters and Getters

/**
* Receives the inputs & returns the output of this Neuron.
* Sums the multiplication of the inputs and weights, then determines which activation function to use.
* @param in The new inputs for this Neuron.
* @return	The single output for this Neuron.
**/
public double think(double[] in){
	//update the input private member variable
	input = in;
	
	double sum = 0.0;
	//sum the multiplication of the inputs and weights
	for(int i = 0; i < num_inputs; i++){ sum += (weight[i] * in[i]);}
	
	if(currentActivationFunction.equalsIgnoreCase("linear")){ sum = lineFunc(sum - theta);}
	else if(currentActivationFunction.equalsIgnoreCase("step")){ sum = stepFunc(sum - theta);}
	else if(currentActivationFunction.equalsIgnoreCase("sign")){ sum = signFunc(sum - theta);}
	else if(currentActivationFunction.equalsIgnoreCase("sigmoid")){ sum = sigmFunc(sum - theta);}
	else if(currentActivationFunction.equalsIgnoreCase("hyperbolic tangent")){ sum = tanhFunc(sum - theta);}
	
	return sum;
}

//Activation Functions

/**
* Linear Function Y = X.
* @param sum 
* @return	The value of the Linear Function.
**/
private double lineFunc(double sum){ return sum;}

/**
* Step Function Y = X >= 0 ? 1.0 : 0.0.
* @param sum 
* @return	The value of the Step Function.
**/
private double stepFunc(double sum){ return sum >= 0.0 ? 1.0 : 0.0;}

/**
* Sign Function Y = X >= 0 ? 1.0 : -1.0
* @param sum
* @return	The value of the Sign Function.
**/
private double signFunc(double sum){ return sum >= 0.0 ? 1.0 : -1.0;}

/**
* Sigmoid Function Y = 1.0 / (1.0 + e^(-X)).
* @param sum
* @return	The value of the Sigmoid Function.
**/
private double sigmFunc(double sum){ return (1.0 / (1.0 + Math.exp(-1.0 * sum)));}

/**
* Hyperbolic Tangent Function Y = (e^(2X) - 1.0)/( e^(2X) + 1.0)
* @param sum
* @return	The value of the Hyperbolic Tangent Function.
**/
private double tanhFunc(double sum){ return ((Math.exp(2.0 * sum) - 1.0) / (Math.exp(2.0 * sum) + 1.0));}

//End Activation Functions

//Back Propagation Elements

/**
* Updates the delta weights that are used for Back Propagation.
* @param delta The offset for the new weights.
**/
public void setDeltaWeights(double[] delta){ deltaWeights = delta;}

/**
* This is pretty much here only for standards purposes.
* It is never called, as I; caine2003; don't give a damn about these values.
* @return	The last delta weights of this Neuron for Back Propagation.
**/
public double[] getDeltaWeights(){ return deltaWeights;}

/**
* Updates the delta threshold that is used for Back Propagation.
* @param delta The offset for the new threshold.
**/
public void setDeltaTheta(double delta){ deltaTheta = delta;}

/**
* This is pretty much here only for standards purposes.
* It is never called, as I; caine2003; don't give a damn about the value.
* @return	The last delta threshold of this Neuron for Back Propagation.
**/
public double getDeltaTheta(){ return deltaTheta;}

/**
* Updates the propagation error for this Neuron.
* @param error New progagation error for this Neuron
**/
public void setPropagateError(double error){ propagate_error = error;}

/**
* Returns the current propagation error for this Neuron.
* @return	The current propagation error for this Neuron
**/
public double getPropagateError(){ return propagate_error;}

/**
* Updates the weights and threshold of this Neuron by using the Delta Weights and Delta Threshold for this Neuron.
**/
public void propagateWeightsAndThetas(){
	for(int i = 0; i < num_inputs; i++){ weight[i] += deltaWeights[i];}
	theta += deltaTheta;
}

//end class Neuron
}
//end class NeuronLayer
}
//end Network class
}