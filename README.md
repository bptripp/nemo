nemo
====

A neural simulator for Matlab that is based on the Neural Engineering Framework (NEF).


Relationship with Other NEF Simulators 
--------------------------------------

There are other NEF-based simulators, including: 

* NESim (Eliasmith & Anderson's original Matlab code)
* Nengo (There are two things that have been called Nengo, Java/Python code that was used to build Spaun, and newer Python code that was designed to address the frustrations of building Spaun)

Nemo is a port of the Java version of Nengo to Matlab. It is meant to be simple and to allow use of Matlab toolboxes with NEF models. 

It is SLOW. This is partly because Matlab OO is slow, and partly because the design favours certain kinds of modifiability over speed. It is most useful for tinkering with NEF-related theory in small models, e.g. finding decoders using the optimization toolbox. 

