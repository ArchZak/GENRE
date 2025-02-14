# GENRE (Gene Expression Network Reading with Explanations)

## Overview
GENRE is a tool that integrates LLM feedback into prokaryotic genome annotation. This repo is an example of its proof of concept before I work to fully deploy this as a web application that anyone who is either a professional, or simply interested in prokaryotic genome annotation can use.

## Current Implementation
The proof of concept currently demonstrates:

- Integration of multiple bioinformatics tools onto a single platform:
  - Prokka (start and stop site prediction + overlap)
  - GeneMarkS2 (start and stop site prediction + overlap and coding potential analysis)
  - BLASTP (functional analysis)
  - HHSearch (functional analysis)
- Parallel execution of analysis tools using ThreadPoolExecutor
- Basic comparison of start and stop between Prokka and GeneMarkS2 predictions
- Functional analysis comparison between HHSearch and BLASTp

## Proof of Concept Features
- Automated submission and retrieval of results from Galaxy tools, BLASTP over the internet, and locally running GenMarkS2
- Concurrent processing of multiple analysis tools
- Basic parsing and comparison of gene predictions
- Returning LLM feedback regarding the top three functional predictions onto the terminal with a terminal user interface

## Work in Progress
Currently working on:
- Implementing user control structures and enchancing error handling
- Developing a more sophisticated user interface in the form of a web app

## Future Plans
Aside from the web development. I plan to incorporate many quality of life features:
- A cloud system using MongoDB to store previous jobs if the FASTA sequence is recognized 
- I would like to reintegrate the AUGUSTUS tool to see if I can incorporate eukaryotic annotations
- I would like to locally host BLASTp for faster searchers
- If resources ever permit me to, I would like to toy with the idea of using my own model. I think it might be a good idea to use a transformer, and then train it on data from the  model I'm working on frokm OpenAI
- Add jobs or listeners to periodically check galaxy and see how the storage is doing

## Setup

source .venv/bin/activate
requirements.txt coming soon

## Note for Readers
This project demonstrates:
- Integration of complex bioinformatics tools
- Parallel processing implementation
- API interaction with Galaxy
- Error handling in bioinformatics contexts
- Modern Python practices (async operations, type hints)

The current implementation serves as a foundation for a more comprehensive gene analysis platform, showcasing both technical capabilities and bioinformatics domain knowledge.