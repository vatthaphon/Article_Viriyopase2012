# Feasibility of long-range zero-lag synchronization

<p align="justify">Investigate the synchronization of a network of three neuronal groups, in which the two outer groups communicate through the inner group.</p>

<p align="justify">This repository contains the source code accompanying the paper titled "When long-range zero-lag synchronization is feasible in cortical networks", Frontiers in Computational Neuroscience. 
The goal of this study is to theoretically and numerically investigate when zero-lag synchroniztion between the two outer neuronal groups is possible.</p>

## Table of Contents  
[1. Introduction](#Introduction)  
[2. Installation](#Installation)  
[3. Usage](#Usage)  
[4. File Structure](#FileStructure)  
[5. License](#License)  
          
## Introduction<a name="Introduction"/>
<p align="justify">This repository provides the implementation of the methods and algorithms described in the paper "When long-range zero-lag synchronization is feasible in cortical networks". 
The main goal of this study is to theoretically and numerically investigate at what parameters of the network, e.g. delays and synaptic strengths, the zero-lag synchronization between the two outer neuronal groups is possible. 
We generally find that the long-range zero-lag synchronization between the two outer groups appears when the connections between the inner neuronal group and the two outer neuronal groups are effectively similar.</p>

## Installation<a name="Installation"/>
- Matlab 2010.

## Usage<a name="Usage"/>
- Run the Matlab files in each directory to generate dynamics of the network.

## File Structure<a name="FileStructure"/>
```plaintext
├── sim_3MS_EventBased/                # Generate dynamics using the event-based approach.
├── sim_3MS_TimeBased/                 # Generate dynamics using the time-based approach.  
├── sim_3MS_TimeBased_with_Bi_Poo_LW/  # Generate dynamics using the time-based approach with STDP.
├── draw_all_line_separation.m         # Draw the line in Figure A1 of the paper.   
├── Viriyopase2012.pdf                 # The article.   
├── README.md  
└── LICENSE
```

## License<a name="License"/>
This project is licensed under the MIT License - see the LICENSE file for details.

