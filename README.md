# 3DOF_Hand_Vibration_Absorber

This repository contains the formulas, equations and simulation codes for a 3 DOF hand model vibration absorber.

To model the human hand, we have taken it to be a double pendulum. we have fixed the shoulder joint and have taken all the motions in thewrist and forearm to be planar. This will represent the Flexion-Extension motion of the wrist and elbow. On the forearm of the patient, a vibration absorber is placed to absorbe the vibration of the hand.

Absorber is a SMA (Shape Memory Alloy) spring connected to the forearm and from ots other side, a mass is connected to it. The SMA spring is in the PE (Pseudo Elastic) phase and because if the SMA's histersetic behaviour in PE phase, it can be modeled as a spring-damper system.

------------------------------------------------------------------
The folder schemes:
1- Model: Contains the formulas, Proofs, Equations and mathematica formuations.
2- Simulation: Contains the MATLAB coeds for brinson's model, its application and a SIMULINK simulation of the entire model using the said MATLAB model.
3- Optimization: The codes for optimizing the model parameters.

------------------------------------------------------------------
The generated chatrts by this simulation:
1-Displacement fo the hand and forearm with respect to frequency and amplitude of the vibration.
2-Amplitude of vibration with respect to absorber mass.
3-Acceleration vs time and frequency
4-Comparison between acceleration/amplitude in damped and undamped system
5-POWER SPECTRAL DENSITY OF THE JOINT(?)
6-The charts should be in angle units as well (obviously).
