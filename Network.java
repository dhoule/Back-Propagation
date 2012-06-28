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
	/**activation function that is to be used for this object **/
	private String current_activation_function;
	
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
		neuronsInHiddenLayers = Stats.getNumPerHiddenLayer();
		//get the activation function that will be used for this Network
		current_activation_function = Stats.getActivationFunction();
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
				if(i < hidden_layers){ layers.add(new NeuronLayer(neuronsInHiddenLayers[i], newInput, current_activation_function));}
				else{ layers.add(new NeuronLayer(outputs, newInput, current_activation_function));}
			}			
		}
		else{
			//there are no hidden layers.
			layers.add(new NeuronLayer(outputs, inputs, current_activation_function));
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
		for(int i = 0; i < (hidden_layers + 1); i++){
			r_weights = extendArray(r_weights, layers.get(i).getLayerWeights());
		}		
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
		//get the learning rate for this network
		double alpha = Stats.getLearningRate();
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
		this.layers.get(numLayer - 1).setDeltaThetas(deltaThetas);
		
		
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
			this.layers.get(i).updateLayerPropagationError(error);
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
			this.layers.get(i).setDeltaThetas(deltaThetas);			
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
		System.out.println("# of inputs " + inputs);
		System.out.println("# of outputs " + outputs);
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
//end Network class
}