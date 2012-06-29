/* /////////////////////////////////////////////////////////////////////
//
// CLASS       	: NeuronLayer
// AUTHOR(S)   	: Houle, Daniel B
// DESCRIPTION 	: Holds mulptiple Neurons for the Network.
// DEPENDENCIES   : Neuron must be in the same directory
//
///////////////////////////////////////////////////////////////////// */

import java.util.*;
import java.util.ArrayList.*;
import java.lang.*;

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
	//set the # of Neurons in this Layer
	num_neurons = num;
	//create the Neurons for this Layer
	for(int i = 0; i < num; i++){ neurons.add(new Neuron(inputs, (new Random()).nextInt(), func, alpha));}
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
	if(weights.length != num_weights){
		System.err.println("\n\nThe # of weights sent to this Layer is incorrect.\n\n");
		System.exit(1);
	}
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
* Returns a single Neuron.
* @param i The index of the Neuron requested.
* @return	The specified Neuron
**/
public Neuron getSingleNeuron(int i){ 
	if(i > num_neurons){
		System.err.println("\n\nThe index is out of range for this Layer.\n\n");
		System.exit(1);
	}
	return neurons.get(i);
}

/**
* Updates the propagation error for every Neuron in thist Layer.
* @param error The propagation errors for this Layer's Neurons.
**/
public void setLayerPropagationError(double[] error){
	if(error.length != num_neurons){
		System.err.println("\n\nThe # of propagation errors sent to this Layer is incorrect.\n\n");
		System.exit(1);
	}
	for(int i = 0; i < num_neurons; i++){ neurons.get(i).setPropagateError(error[i]);}
} 

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
public void setLayerThresholds(double[] thetas){
	if(thetas.length != num_neurons){
		System.err.println("\n\nThe # of thresholds sent to this Layer is incorrect.\n\n");
		System.exit(1);
	}
	for(int i = 0; i < num_neurons; i++){ neurons.get(i).setTheta(thetas[i]);}
}

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
	if(in.length != (num_weights/ num_neurons)){
		System.err.println("\n\nThe # of inputs sent to this Layer is incorrect.\n\n");
		System.exit(1);
	}
	//Each Neuron receives multiple inputs but will only give ONE output
	double[] outputs = new double[num_neurons];
	//Each Neuron in this layer will receive the same inputs
	for(int i = 0; i < num_neurons; i++){ outputs[i] = neurons.get(i).think(in);}
	return outputs;
}

/**
* Returns this Layer's last inputs.
* @return	The last inputs of this Layer.
**/
public double[] getLayerInput(){
	double[] input = new double[num_weights];
	int ct = 0;
	int numInput = num_weights/ num_neurons;
	//go through each Neuron and get it's last input
	for(int i = 0; i < num_neurons; i++){
		double[] temp = neurons.get(i).getInput();
		//go through each Neurons inputs & add them to input
		for(int q = 0; q < numInput; q++){
			input[ct] = temp[q];
			ct++;
		}
	}
	return input;
}

/**
* Updates the delta weights of this Layer for Back Propagation purposes.
* @param delta The offset for each weight in this Layer.
**/
public void setLayerDeltaWeights(double[] delta){
	if(delta.length != num_weights){
		System.err.println("\n\nThe # of delta weights sent to this Layer is incorrect.\n\n");
		System.exit(1);
	}
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
//end class NeuronLayer
}