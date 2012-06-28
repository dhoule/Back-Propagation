/* /////////////////////////////////////////////////////////////////////
//
// CLASS       	: Stats
// AUTHOR(S)   	: Houle, Daniel B
// DESCRIPTION 	: Holds static information so it may be used through
//							out the project.
// DEPENDENCIES   : Methods must be called correctly.
//
///////////////////////////////////////////////////////////////////// */

public class Stats{

/** # of inputs into the Network **/
private static int inputs;
/** # of outputs from the Network & the # of Neurons in the last Layer **/
private static int outputs;
/** # of hidden layers in the Network **/
private static int hidden_layers;
/** Activation function to be used by every Neuron in the Network **/
private static String activation_function;
/** The learning rate for every Neuron in the Network. 0.0 < alpha < 1.0 **/
private static double alpha;
/** The acceptable error rate for the Network. Usually 0.001. **/
private static double allowed_error_rate;
/** An array of the # of Neurons in each hidden layer from top-down **/
private static int[] num_in_hidden_layers;


/**
* Updates the learning rate that will be used by every Neuron in the Network.
* @param lr The new learning rate.
**/
public static void setLearningRate(double lr){ alpha = lr;}

/**
* Returns the learning rate that the Network uses.
* @return	The learning rate.
**/
public static double getLearningRate(){ return alpha;}

/**
* Updates the # of inputs the Network will have.
* @param num The new # of inputs
**/
public static void setInputNumber(int num){ inputs = num;}

/**
* Returns the # of inputs the Network has.
* @return	The # of inputs the Network has.
**/
public static int getInputNumber(){ return inputs;}

/**
* Updates the # of outputs the Network will have. This is also the # of Neurons in the last layer.
* @param num The # of outputs from the Network & Neurons in the last Layer. 
**/
public static void setOutputNumber(int num){ outputs = num;}

/**
* Returns the # of outputs from the Network and the # of Neurons in the output Layer.
* @return	The # of outputs from the Network and the # of Neurons in the output Layer.
**/
public static int getOutputNumber(){ return outputs;}

/**
* Updates the Activation Function that every Neuron will be using.
* @param func The new Activation Function the Network will use.
**/
public static void setActivationFunction(String func){ activation_function = func;}

/**
* Returns the Activation Function that the Network will use.
* @return	The Activation Function that will be used.
**/
public static String getActivationFunction(){ return activation_function;}

/**
* Updates the allowed error rate to confirm when the Network has reached a satisfactory output.
* @param error The new error rate.
**/
public static void setAllowedErrorRate(double error){ allowed_error_rate = error;}

/**
* Returns the allowed error rate that the Network is to use.
* @return	The allowed error rate.
**/
public static double getAllowedErrorRate(){ return allowed_error_rate;}

/**
* Updates the # of hidden layers that will be in the Network.
* @param num The # of hidden layers that will be in the Network.
**/
public static void setNumberOfHiddenLayers(int num){ hidden_layers = num;}

/**
* Returns the # of hidden layers that the Network will have.
* @return	The # of hidden layers in the Network.
**/
public static int getNumberOfHiddenLayers(){ return hidden_layers;}

/**
* Updates the # of Neurons per hidden layer.
* @param layers An array of the # of Neurons per hidden layer.
**/
public static void setNumPerHiddenLayer(int[] layers){ num_in_hidden_layers = layers;}

/**
* Returns the # of Neurons per hidden layer.
* @return	An array of the # of Neurons per hidden layer.
**/
public static int[] getNumPerHiddenLayer(){ return num_in_hidden_layers;}
//end class Stats
}