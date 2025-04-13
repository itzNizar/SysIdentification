🔧 System Identification of a Dynamic System Using MATLAB
📌 Project Summary
This project focuses on identifying the dynamic behavior of an unknown system using experimental input-output data. Two ARX (Auto-Regressive with eXogenous input) models of different orders were implemented and evaluated using both graphical methods and the Least Squares (LS) algorithm.

🎯 Objectives
Load and analyze time-series data from a physical system.

Perform graphical estimation of system dynamics (gain & time constant).

Implement the Least Squares algorithm manually in MATLAB.

Compare ARX(1,1,1) and ARX(2,2,1) model structures.

Validate the models and analyze the goodness of fit.

🛠️ Methods Used
Graphical Estimation of ARX(1,1,1) parameters.

Manual LS Implementation for parametric identification.

Model Validation using 50% of the data for testing.

System Identification Toolbox (arx, compare) for model validation.

📊 Results
Model	Fit (%) on Validation	Cost on Estimation Data
ARX(1,1,1)	95.26%	19.38
ARX(2,2,1)	94.37%	14.58
The ARX(1,1,1) model offered better generalization on unseen data.

The ARX(2,2,1) model better captured finer dynamics (lower cost), but with a small risk of overfitting.

🧠 What I Learned
Practical use of system identification techniques from raw data.

Hands-on experience in parameter estimation and model validation.

Trade-offs between model complexity and generalization.

How to cross-check manual implementations with built-in MATLAB tools.

🧪 Tools & Technologies
MATLAB

System Identification Toolbox

Data visualization & time-series analysis
