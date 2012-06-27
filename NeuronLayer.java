
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

/** Constructor must be told the number of Neurons that are in this layer, 
    how many inputs each Neuron takes, and the name of the activation 
	 function that each Neuron will use. 
 @param num # of Neurons in this Layer
 @param inputs # of inputs into this Layer. Every Neuron recieves the same input.
 @param func The activation function for each Neuron in this Layer. 
**/
public NeuronLayer(int num, int inputs, String func){
	//set the # of Neurons in this Layer
	num_neurons = num;
	//create the Neurons for this Layer
	for(int i = 0; i < num; i++){ neurons.add(new Neuron(inputs, (new Random()).nextInt(), func));}
	//set the # of weights for this Layer
	num_weights = num * inputs;
}

/**
* Returns the # of Neurons in this Layer.
* @ return	The # of Neurons in this Layer.
**/
public int getNumLayer(){ return num_neurons;}

/**
* This is pretty much here only for standards purposes.
* @param nothing 
**/
public void setNumLayer(int nothing){ num_neurons = nothing;}
//end class NeuronLayer
}