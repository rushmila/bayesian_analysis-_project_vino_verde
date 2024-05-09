Assessment 3 

Bayesian analysis project: Vinho Verde

For this assessment, you will perform a logistic regression on a true dataset. 

Wine dataset
This datasets is related to red variants of the Portuguese "Vinho Verde" wine. The dataset is described in the publication by Cortez, P., Cerdeira, A., Almeida, F., Matos, T., & Reis, J. (2009). Modeling wine preferences by data mining from physicochemical properties.

Dataset available here.

The input variables (based on physicochemical tests) are:

- fixed acidity

- volatile acidity

- citric acid

- residual sugar

- chlorides

- free sulfur dioxide

- total sulfur dioxide

- density

- pH

- sulphates

- alcohol

The output variable (based on sensory data) is quality (a score between 0 and 10). 




1. Read the dataset into R
Check if there are missing values (NA) and, in case there are, remove them.


2. We want to implement a logistic regression, therefore we want a response variable which assume values either 0 or 1. Suppose we consider "good" a wine with quality above 6.5 (included).


3. Run a frequentist analysis on the logistic model, using the glm() function. What are the significant coefficients?


4. Estimate the probabilities of having a "success": fix each covariate at its mean level, and compute the probabilities for a wine to score "good" varying 
total.sulfur.dioxide, and plot the results.


5. Perform a Bayesian analysis of the logistic model for the dataset, i.e. approximate the posterior distributions of the regression coefficients, following these steps: 

- Write an R function for the log posterior distribution.

- Fix the number of simulation at 10^4.

- Choose 4 different initialisations for the coefficients.

- For each initialisation, run a Metropolis–Hastings algorithm.

- Plot the chains for each coefficients (the 4 chains on the same plot) and comment.


6. Approximate the posterior predictive distribution of an unobserved variable characterised by 

- fixed acidity: 7.5

- volatile acidity: 0.6

- citric acid: 0.0

- residual sugar: 1.70

- chlorides: 0.085

- free sulfur dioxide: 5

- total sulfur dioxide: 45

- density: 0.9965

- pH: 3.40

- sulphates: 0.63

- alcohol: 12 

Plot the approximate posterior predictive distribution.


Use the metrop() function available in the  mcmc package to perform the same analysis on the posterior distribution you have approximated for Question 6. 
Choose again 10^4 simulations and compare the results with the results obtained with your code. 
