from importlib.resources import files
import os

def main():
    #checking sequins data from the resources
    gtf_file = files("data").joinpath("rnasequin_annotation_2.4.gtf")
    count = 0
    with open(gtf_file, 'r') as mygtf:
        for line in mygtf:
            print(line)
            count +=1
            if count >= 10:
                break
    print("**finish_checking_the_files**")
