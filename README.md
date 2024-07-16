# Modeling-and-Control-of-a-Multi--Prismatic-Robotic-System
This project was developed as a part of the MAE547 Modeling and Control of Robots Course. 

In this project, we developed an n-prismatic links robotic system with indirect force control under compliance and impedance control.


**Part 1: Joint Variables Plots and Torque equations**

With DH parameters as inputs, the code assists the user in generating transformation matrices, position vectors, and the jacobians associated with each of the links and motors and evaluates the equations of dynamics of motion.
The final output of the first part of the code generates the torque equations along with the plots of joint positions, velocities and accelerations with respect to time.


**Part 2: Compilance Control**

For performing compliance control, the user can input proportional and derivative gain values along with the desired end effector position. On applying PD control with gravity compensation, we can perform compliance control and get the plots of desired vs actual end effector positions as well as the contact forces.


**Part 3: Impedance Control**

For performing impedance control, the user can give input for proportional and derivative gains along with any other mechanical forces being applied to the joints. When applying impedance control, the user can plot desired and actual end effector positions as well as contact forces.


**Graphical User Interface (GUI):**

Grid layout, push buttons, table and edit numerical elements are used to generate the graphical user interface with functions and call back functions integrated into the script. The user can provide all the necessary input parameters using the easy-to-understand interface and access scripts to obtain the necessary plots.
