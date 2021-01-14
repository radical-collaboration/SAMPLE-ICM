## Title: SAMPLE: A framework to facilitate ensemble sampling (needs work)

### Introduction
- Ensembles require large amounts of computational power to process
- Computational resource usage can be minimized through intelligent sampling with minimal impact to the accuracy of results
- We developed the SAMPLE framework to facilitate the sampling of ensembles.
- SAMPLE leverages the RADICAL-toolkit as a backend for efficient task distribution and scheduling


### Methods
- Global design (sample ensembles, run, evaluate, terminate/back to step 1)
- Server side
     - communicates with radical to launch tasks
     - determines whether sampling results have converges or need to proceed
- Client side
    - evaluation of ensemble results
    - communication to server of the results

- Experiments: ice sheet modelling
    - sampling algorithms used

### Results
- Results of experiments
    - how many resources were saved
    - how much time was saved
    - accuracy of results

### Discussion
- different strategies to evaluate ensemble results should be explored


### Conclusion
