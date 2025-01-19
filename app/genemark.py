import subprocess
import os
import argparse

def fix_file(fasta_file):
    with open(fasta_file, "r") as file:
        lines = file.readlines() 
        sequence = [] 
        header_exists = False
        
        for line in lines: 
            line = line.strip()
            if line.startswith(">"): #needs > line or else genemark returns error
                header_exists = True
                header = line
            else:
                sequence.append(line)

    sequence = "".join(sequence) #basically create one massive line
    
    new_file = f"fixed_{os.path.basename(fasta_file)}"
    with open(new_file, "w") as output_file:
        if not header_exists:
            output_file.write(">sequence\n") #need to add a sequence line or else genemark wont take it
        else:
            output_file.write(f"{header}\n")

        for i in range(0, len(sequence), 60):
            output_file.write(sequence[i:i+60] + "\n") #proper fasta file practices of 60 chrs per line
    
    return new_file #return str name, not physical file

def genemark(fasta_file: str, prokaryotic: bool = True):
    if not os.path.exists(fasta_file): #checks to see if fasta_file exists, should return page error later
        raise FileNotFoundError(f"given fasta file not found: {fasta_file}")

    fixed_fasta_file = fix_file(fasta_file) #file needs to be fixed or else it all breaks, also needs to be str
    
    command = ["perl", "../gms2_macos/gms2.pl"] #run subprocess to access tool
    command.extend([
        "--seq", fixed_fasta_file,
        "--genome-type", "bacteria" if prokaryotic else "auto",
        "--ps-output"
    ])
    
    try:
        subprocess.run(
            command,
            text=True,
            capture_output=True,
            check=True
        )
        result_file = "gms2.lst"
        if os.path.exists(result_file):
            print(f"results saved to {result_file}")
            
    except subprocess.CalledProcessError as e: #needs this specific syntax for some reason, dont touch excepts
        print(f"Error running GeneMarkS: {e}")
        if e.stderr: 
            print("\nError output:") 
            print(e.stderr)
        if e.stdout:
            print("\nStandard output:")
            print(e.stdout)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

def main(): #just an arg parser to retreive file, will need to automate this when creatin web tool
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_file")
    args = parser.parse_args()
    genemark(args.fasta_file)

if __name__ == "__main__": 
    main()