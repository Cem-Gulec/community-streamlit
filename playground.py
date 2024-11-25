import subprocess

# List of strings in Python
#input_strings = ["T:CALM1_HSPC:PDE1B", "T:CALM2_HSPC:PDE1B", "HSPC:CALM2_HSPC:PDE1B", "HSPC:CALM1_HSPC:PDE1B", "HSPC:CD47_HSPC:SIRPA", "B:CALM2_HSPC:PDE1B", "B:CALM1_HSPC:PDE1B", "Mono:CALM1_HSPC:PDE1B", "Mono:CALM3_HSPC:PDE1B", "Mono:NRG1_Gran:NETO2", "Mono:NRG1_Gran:MS4A4A", "Mono:NRG1_Gran:HLA-DPB1", "Mono:TNFSF14_Gran:LTBR", "Mono:TNFSF14_Gran:TNFRSF14", "Mono:TNFSF14_Mono:TNFRSF14", "Mono:NRG1_Mono:SIGLEC7", "Mono:HP_Gran:CD163", "Mono:HP_Gran:ASGR2", "Mono:HP_Gran:ASGR1"]
input_strings = []

# Join the list into a single string with spaces (or any separator you like)
string_arg = " ".join(input_strings)

# Run the R script with subprocess, passing the strings as arguments
try:
    subprocess.run(
        ["Rscript", "heatmap.R", string_arg],
        text=True,
        check=True
    )
except subprocess.CalledProcessError as e:
    print("An error occurred:", e)
    print("R script error output:", e.stderr)
