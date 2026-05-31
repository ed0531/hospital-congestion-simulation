\# Hospital Congestion Simulation



This project implements a Monte Carlo simulation of emergency department congestion using a stochastic arrival process with Hawkes-like self-excitation dynamics.



Patient arrivals follow a time-varying intensity with daily seasonality and short-term clustering effects. Incoming patients are assigned to triage categories (white, green, yellow, and red) according to empirical probabilities calibrated from real emergency department data.



The model incorporates congestion-dependent discharge rates: as the number of patients increases, service efficiency decreases, generating potential overcrowding effects. Repeated simulations are used to estimate the expected occupancy level, triage composition, and the probability of exceeding critical capacity thresholds throughout the day.



The project is intended as an exploratory framework for studying emergency department dynamics, patient flow, and overcrowding risk under stochastic demand.



