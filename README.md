# FishPoliceAnalysis

This repository holds the anonymized data from 3 FishPolice Experiments together with the matlab programs to analyze that data. 

## loadData
This reads in the data and stores it in tables e1 e2 e2pro e3 e3pro and a table containing all of the others e

## runDemo
This allows one to run various analysis tools on an of the tables. In the near future we will document all of the analysis tools here.

Here is an example of showing a visualization of the log files of all novice players in the 2nd experiment
```
runDemo(e2,'gameplay5')
```


  ![Novice Game Players Log](/images/day2noviceusers.png "Day 2 Novice Users")
  
