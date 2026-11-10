# Generate macros/HSCPSelector.h and macros/HSCPSelector.C from the templates
# macros/HSCPTemplate.h and macros/HSCPTemplate.C.
#
# The templates are copied line by line. After each marker line, the code generated
# from the preselections listed in cfg/HSCPpreSelections.csv is inserted:
#   ADD-HSCP-SELECTION     (HSCPTemplate.h)  declarations of the PassHSCPpresel_<label> methods
#   ADD-SELECTION-METHODS  (HSCPTemplate.C)  their implementations
#   FILL-SELECTION-VECTOR  (HSCPTemplate.C)  registration of each method and of its label
#
# Must be run from the python/ directory (the paths below are relative):
#   cd python && python3 CreateSelector.py

from Functions import *
import csv


#######################################
# Take care of the HSCP preselection
#######################################
# Read config files
#######################################

# list of (label, expression) pairs; the first line of the file is a header
with open("../cfg/HSCPpreSelections.csv", newline='') as f:
    rows = [row for row in csv.reader(f, skipinitialspace=True) if row]
HSCPsel = [(row[0], row[1]) for row in rows[1:]]
print('Selections are:')
for label, instruction in HSCPsel:
    print('  '+label)

# header files: template and new one
h_temp = open("../macros/HSCPTemplate.h",'r')
h_out = open("../macros/HSCPSelector.h",'w')

for line in h_temp:
    #copy the line of the code
    h_out.write(line)
    #add code to declare data members corresponding to the HSCP preselections
    if line.find("ADD-HSCP-SELECTION")>0:
        content = Code_HSCPSelectorProto(HSCPsel)
        h_out.write(content)

h_out.close()
h_temp.close()

#  C files: template and new one
c_temp = open("../macros/HSCPTemplate.C",'r')
c_out = open("../macros/HSCPSelector.C",'w')

for line in c_temp:
    #copy the line of the code
    c_out.write(line)
    #add code for the HSCP preselection methods
    if line.find("ADD-SELECTION-METHODS")>0:
        content = Code_HSCPSelectorImpl(HSCPsel)
        c_out.write(content)
    #add code to add labels & members for HSCP preselection into vectors
    if line.find("FILL-SELECTION-VECTOR")>0:
        content = Code_HSCPSelectorAddLabelsAndPointers2Vector(HSCPsel)
        c_out.write(content)

c_out.close()
c_temp.close()
