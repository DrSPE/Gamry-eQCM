# Gamry-eQCM
Scripts to import and resort Data from the Gamry file format .DTA into a better analysible format

A collection of scripts to convert the output from the Gamry eQCM setup (Garmin QCM + Potentiostat, see https://www.gamry.com/qcm/quartz-crystal-microbalance/eqcm-quartz-crystal-microbalance-2/) into a file that can be plotted easily 
without having to use Garmins plotting software first.

Contains a script to sort files into their own folder.
The QCM import script splits the Garmin file format which contains the frequency data and the echem data (CV only) in different blocks and unequal spacing
and merges the same data point into the same row while interpolating missing points to make sure they are equidistant.
