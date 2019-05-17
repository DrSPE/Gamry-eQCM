# Gamry-eQCM
Scripts to import and resort Data from the Gamry file format .DTA into a better analysible format

A collection of scripts to convert the output from the Gamry eQCM setup (Gamry’s eQCM 10M + Potentiostat e.g. Interface 1010, see https://www.gamry.com/qcm/quartz-crystal-microbalance/eqcm-quartz-crystal-microbalance-2/) into a file that can be plotted easily 
without having to use Garmins plotting software first.

Contains a script to sort files into their own folder.
The QCM import script splits the Garmin file format which contains the frequency data and the echem data (CV only) in different blocks and unequal spacing and merges the same time point into the same row while interpolating missing points to make sure they are equidistant.


Used for the analysis performed in my Master thesis as well as the paper for which the results were used as well: 
[1] S. P. Emge, "Einblicke in das kapazitive Ladeverhalten von porösen Kohlenstoffelektroden in verschiedenen ionischen Flüssigkeiten durch EQCM Messungen" (Master thesis), 2016, Philipps-Universität Marburg, Marburg, Germany.

[2] N. Jäckel, S. P. Emge, B. Krüner, B. Roling, V. Presser, "Quantitative Information about Electrosorption of Ionic Liquids in Carbon Nanopores from Electrochemical Dilatometry and Quartz Crystal Microbalance Measurements" J. Phys. Chem. C 2017, 121, 19120.
DOI: 10.1021/acs.jpcc.7b06915
