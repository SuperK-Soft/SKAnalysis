# Data Model
*************************

conceptually the DataModel class is a transient data storage class thats accessable from every Tool within a ToolChain. It is the primary location for storing and passing data between Tools.

The DataModel class can therefore be adapted and expanded to include any data objects or classes that you require in your applicaion

In it's initial blank state it comes with just two variables:

 [ Store vars; ]: An ascii storage class which maps string keys to any variable type (see Store class in ToolFramework). This class houses by defualt all the variables passed to a ToolChain in the ToolChain config as well as all the command line argumnets at execution. It also has some flag vaiables which can be read about int he manual for things like stopping execution when looping.
 
 and
 
 [ Logging *Log; ]: A lass for sending logg messages similar to c++ cout with a redirectable stream (see Logging class in ToolFramework/ToolDAQ)



A TTree map with getter and setter functions is provided and can be uncommented if required for those who wish to organise root TTrees (note: requires root dependancy)

 
