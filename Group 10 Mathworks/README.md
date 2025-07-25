Group 10's Optimization Battery Profile - Mathworks Project

## Project Description

This project models and simulates the **AC and DC fast-charging behavior** of the BMW i7's battery using MATLAB.
    It includes RC (resistor-capacitor) circuit-based approximations to evaluate voltage and energy profiles during the charging process.
    Users can define their own battery specifications to observe how different parameters affect the charging time and behavior.

## Methods & Results

AC Simulation (AC_Charger_BMW.m): Models Level-2 charging (~11 kW) using an RC circuit approximation to simulate voltage rise over time.
The script estimates how long it takes to reach 80% and 100% charge, and displays a voltage vs. time plot.

DC Simulation (DC_Charger_BMW.m): Simulates high-speed fast-charging behavior with nonlinear voltage growth, optimized for rapid 80% and full charge milestones.
The script highlights key voltage targets and calculates charging durations.

User Input Module (User_Battery_Input.m): Allows users to enter custom battery parameters (capacity, resistance, initial voltage), 
which are then used to run simulations with personalized conditions.

Outputs: Each simulation generates voltage-time, power-time. current-time, battery temperature-time, resistance-time, and rate of change-time graphs, marking 80% and 100% charge levels. 
These can be saved to the media/ folder for reporting or presentation.

## ️ How to Use

1. **Launch MATLAB.**
2. **Navigate** to the `Group 10 Mathworks` folder.
3. Run the desired script (choose one at a time to run):
   - `AC_Charger_BMW.m` — for standard AC charging simulations.
   - `DC_Charger_BMW.m` — for fast DC charging evaluation.
   - `User_Battery_Input.m` — to define custom input values like resistance, capacity, or target voltage.

##  Student Contributions

| Student Name   | Contribution                                                |
|----------------|-------------------------------------------------------------|
| Elan Boganim   | Developed DC and User simulation simulation scripts, worked on the majority of the slide show |
| Jordan Kapner  | Contributed to the slide show and worked on debugging   |
| Denise Estana  | Developed the AC simulation script added detail to the slide show |
| Jonathan Herrera-Aguayo | Worked on debugging code and helped with the AC simulation script |
