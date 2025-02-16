# GENRE (Gene Expression Network Reading with Explanations)

## Overview
GENRE is a tool that integrates LLM feedback into prokaryotic genome annotation. This repo is an example of its proof of concept before I work to fully deploy this as a web application that anyone who is either a professional, or simply interested in prokaryotic genome annotation can use.

Currently, researchers all around the world are annotating genomes, and this can be a time-consuming process. During my second semester of research, I was manually annotating the Jollymon phage genome. This proved to be a menial process, and I wondered if this could be automated. I then thought about my internship in LLM feedback for OSCEs, and then GENRE was born! Now, as I'm getting more experienced in genetics and computer science conventions, I am trying to replicate new base-lines using code and LLM feedback.

This is currently the first iteration. I wanted to see if this was possible, since I hadn't seen anyone else try this yet. As an aspiring computational biologist, I believe I have a duty to create tools that can support my fellow experimentalists. GENRE will be no exception to that duty. I want to try and push this project to be more than just a resume piece, but something that can serve as a building block for something greater. 

## Current Implementation (Proof of Concept)
The current implementation is a temporary proof of concept that demonstrates:

- Integration of multiple bioinformatics tools onto a single platform:
  - Prokka (start and stop site prediction + overlap)
  - GeneMarkS2 (start and stop site prediction + overlap and coding potential analysis)
  - BLASTP (functional analysis)
  - HHSearch (functional analysis)
- Parallel execution of analysis tools using ThreadPoolExecutor
- Basic comparison of start, stop, and overlap between Prokka and GeneMarkS2 predictions
- Functional analysis comparison between HHSearch and BLASTp
- Temporary file-based storage system (to be replaced with cloud database)
- Basic command-line interface (to be replaced with web UI)

## Planned Architecture Improvements
- Moving away from file-based storage to a proper database system:
  - Implementation of MongoDB/some cloud storage for job storage and retrieval
  - Caching of previous results for identical FASTA sequences
  - Structured storage of analysis results and LLM feedback
- Replacing temporary command-line interface with a full web application
- Implementing proper user authentication and job management
- Moving from local file handling to cloud storage for scalability

## Proof of Concept Features
- Automated submission and retrieval of results from Galaxy tools, BLASTP over the internet, and locally running GenMarkS2
- Concurrent processing of multiple analysis tools
- Basic parsing and comparison of gene predictions
- Returning LLM feedback regarding the top three functional predictions onto the terminal with a terminal user interface

## Work in Progress
Currently working on:
- Implementing user control structures and enchancing error handling
- Developing a more sophisticated user interface in the form of a web app
- Using my own model to generate feedback rather than OpenAI's API

## Future Plans
Aside from the web development. I plan to incorporate many quality of life features:
- A cloud system using MongoDB to store previous jobs if the FASTA sequence is recognized 
- I would like to reintegrate the AUGUSTUS tool to see if I can incorporate eukaryotic annotations
- I would like to locally host BLASTp for faster searchers
- If resources ever permit me to, I would like to toy with the idea of using my own model. I think it might be a good idea to use a transformer, and then train it on data from the  model I'm working on frokm OpenAI
- Add jobs or listeners to periodically check galaxy and see how the storage is doing

## Tool Setup

source .venv/bin/activate
pip install -r requirements.txt

In the same folder as `app`, you should create a `blastp_downloads` and `galaxy_downloads` folder. 

### Environment 
Create a `.env` file in the root directory with the following variables:

GALAXY_API_KEY=your_galaxy
HISTORY_ID = your_history_id

The set up is somewhat unconventional after this.

You are going to need to locally download GeneMarkS2. Make sure you put your key in the same folder as the one you plan to run code in.
You will need to create a usegalaxy account, and then create an api key for your own history and galaxy use. The storage is very generous with 250 gigabytes of free storage.

After you do both, integrate a model of your own choice. To simply test out the code, you could use Openai's API. This is a proof of concept. 

## Note for Readers
This project currently demonstrates:
- Integration of various bioinformatics tools across different platforms being used asynchronously
- Parallel processing implementation
- API interaction with Galaxy
- Error handling in bioinformatics contexts
- Modern Python practices (async operations, type hints)
- I have included a random fasta file to test with

The current implementation serves as a foundation for a more comprehensive gene analysis platform. Many current implementations (file handling, CLI interface, etc.) are temporary and will be replaced with more robust solutions in the production version.