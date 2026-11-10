cd python
python3 CreateSelector.py
cd ../macros
time root -l -q -b macro.cc

rm -f ../macros/HSCPSelector_C*