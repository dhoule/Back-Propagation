

import java.util.*;
import java.util.ArrayList.*;
import java.lang.*;

public class Neuron{

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

/** Constructor must be told how many inputs it will be recieving, a seed to randomize the initial 
	 weights & threshold, and the name of the activation function it will be using. 
 @param num # of inputs into this Neuron.
 @param seed Used to seed the Random object used for randomizing the weights and threshold of this Neuron.
 @param func The activation function for this Neuron.
**/
public Neuron(int num, int seed, String func){
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
	//set the size of the inputs
	num_inputs = num;
	//set the new threshold
	Random offset = new Random((long) seed);
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
* @param double threshold New threshold 
**/
public void setTheta(double threshold){ theta = threshold;}

/**
* Returns the threshold for this Neuron.
* @return 	The current threshold for this Neuron 
**/
public double getTheta(){ return theta;}

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
* Updates the weights of this Neuron.
* $param newWeights New weights
**/
public void setWeights(double[] newWeights){ 
	if(newWeights.length != num_inputs){
		System.err.println("\n\nThe number of new weights is not correct for this Neuron.\n\n");
		System.exit(1);
	}
	weight = newWeights;
}

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

//End standard Setters and Getters

/**
* Recieves the inputs & returns the output of this Neuron.
* Sums the multiplication of the inputs and weights, then determines which activation function to use.
* @param in The new inputs for this Neuron.
* @return	The single output for this Neuron.
**/
public double think(double[] in){
	if(in.length != num_inputs){
		System.err.println("\n\nThe number of inputs is not correct for this Neuron.\n\n");
		System.exit(1);
	}
	double sum = 0.0;
	//sum the multiplication of the inputs and weights
	for(int i = 0; i < num_inputs; i++){ sum += (weight[i] * in[i]);}
	
	if(currentActivationFunction.equalsIgnoreCase("linear")){ sum = lineFunc(sum - theta);}
	else if(currentActivationFunction.equalsIgnoreCase("step")){ sum = stepFunc(sum - theta);}
	else if(currentActivationFunction.equalsIgnoreCase("sign")){ sum = signFunc(sum - theta);}
	else if(currentActivationFunction.equalsIgnoreCase("sigmoid")){ sum = sigmFunc(sum - theta);}
	
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
private double sigmFunc(double sum){return (1.0 / (1.0 + Math.exp(-1.0 * sum)));}

//End Activation Functions

//Back Propagation Elements

/**
* Updates the delta weights that are used for Back Propagation.
* @param delta The offset for the new weights.
**/
public void setDeltaWeights(double[] delta){ 
	if(delta.length != num_inputs){
		System.err.println("\n\nThe number of delta weights is not correct for this Neuron.\n\n");
		System.exit(1);
	}	
	deltaWeights = delta;
}

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
* Updates the weights and threshold of this Neuron by using the Delta Weights and Delta Threshold for this Neuron.
**/
public void propagateWeightsAndThetas(){
	for(int i = 0; i < num_inputs; i++){ weight[i] += deltaWeights[i];}
	theta += deltaTheta;
}
//end class Neuron
}