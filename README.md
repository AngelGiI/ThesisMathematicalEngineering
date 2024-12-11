**Optimizing Bicycle Sharing Systems in Urban Networks**

This repository contains the code and models developed for my degree thesis, which focused on optimizing bicycle sharing systems in urban networks. The thesis addressed challenges in rebalancing bicycle availability across stations, ensuring service reliability, and minimizing operational costs.

### Overview:
- **Main Objectives:**
  - Develop a mathematical programming model to optimize the redistribution of bicycles in a city.
  - Incorporate dynamic and spatial aspects of the problem using graph-based modeling.
  - Balance operational costs with user satisfaction by minimizing shortages and surpluses at docking stations.

- **Key Components:**
  1. **Mathematical Models:**  
     Designed and implemented optimization models based on flow networks using Gurobi.
  2. **Simulation Algorithms:**  
     Built heuristic and metaheuristic approaches, including my own proposal: AVET (Algoritmo Voraz Espacio-Temporal i.e. "Greedy Spatio-Temporal Algorithm"), to handle large-scale instances and compare with optimal solutions.

- **Projects Highlight:**
  - **Dynamic Redistribution Model:**  
    Modeled the temporal and spatial flow of bicycles in the network to ensure station balance over time, accounting for fluctuating user demands and vehicle constraints.
  - **Case Study with BiciMAD:**  
    Applied the models to Madrid's bike-sharing system (BiciMAD), analyzing real-world data for station demand and travel patterns. The study included scenarios with randomized instances for validation.

This repository demonstrates the integration of mathematical optimization, simulation, and urban mobility analytics to address a practical, impactful problem.
